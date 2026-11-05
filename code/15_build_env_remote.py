#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""15_build_env_remote.py — 在新 GPU 实例上搭建 B 线计算环境（可重复执行）。

背景
----
旧 AutoDL 实例已报废，环境需重建。MatPool 实例（hz-4）实测：
  - Ubuntu 22.04 / 32 核 / 240 GB RAM / 800 GB 可用盘 / 4×RTX A4000 16GB / CUDA 12.2
  - **GitHub 与 ftp.gromacs.org 不可达**（GROMACS 源码编译路线走不通）
  - 清华 pip/conda 镜像可达（pypi 约 14 MB/s）

因此采用 **AmberTools 单栈方案**：
  tleap(建拓扑) → pmemd.cuda(GPU MD) → cpptraj(去水) → MMPBSA.py(PB)
即不经过 ParmEd→GROMACS 转换，省掉编译 GROMACS 这一步，也少一处出错点。
若远端确实能拿到 GROMACS（GPU 版），再切回 GROMACS 路线。

本脚本做的事
------------
1. conda 换清华源（提速）
2. 创建 env `gamd`：python 3.10 + ambertools + openbabel + parmed + mdtraj + numpy/pandas
3. 校验关键可执行文件是否存在（tleap / antechamber / parmchk2 / cpptraj / MMPBSA.py / pmemd.cuda）
4. 缺少 mdtraj 时用 pip(清华源) 补装
5. 结果写本地 results/15_env_build.json

用法
----
/opt/anaconda3/bin/python rebuild_20261025/code/15_build_env_remote.py \
    --host hz-4.matpool.com --port 26255 --user root --password <pwd> [--poll-min 2]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

try:
    import paramiko
except ImportError:
    sys.exit("缺少 paramiko：/opt/anaconda3/bin/python -m pip install paramiko")

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "rebuild_20261025" / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

CONDA = "/root/miniconda3/bin/conda"
REMOTE_ENV_SH = "/root/autodl_env_build.sh"
REMOTE_LOG = "/root/env_build.log"

BUILD_SH = r"""#!/bin/bash
set -x
export PATH=/root/miniconda3/bin:$PATH
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ --system 2>/dev/null
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ --system 2>/dev/null
conda config --set show_channel_urls yes --system 2>/dev/null
if conda env list | grep -q '^gamd'; then
  echo "[skip] env gamd 已存在"
else
  conda create -y -n gamd -c conda-forge python=3.10 ambertools openbabel parmed mdtraj numpy pandas matplotlib 2>&1
fi
source /root/miniconda3/etc/profile.d/conda.sh
conda activate gamd
python -c "import mdtraj" 2>/dev/null || pip install -q -i https://pypi.tuna.tsinghua.edu.cn/simple mdtraj
echo "=== VERSION ==="
which tleap antechamber parmchk2 cpptraj MMPBSA.py pmemd pmemd.cuda obabel python 2>/dev/null
python --version
echo "=== BUILD_DONE ==="
"""


def poll(c, log: str, phase: str, every_s: int = 30) -> None:
    seen = 0
    while True:
        _, o, _ = c.exec_command(f"tail -n +{seen + 1} {log} 2>/dev/null | tail -40", timeout=60)
        chunk = o.read().decode(errors="replace")
        if chunk:
            print(chunk, end="")
            seen += chunk.count("\n")
        _, o2, _ = c.exec_command(f"grep -c 'BUILD_DONE' {log} 2>/dev/null", timeout=30)
        if o2.read().decode().strip() not in ("", "0"):
            print(f"\n[{phase}] 环境构建完成")
            return
        time.sleep(every_s)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", required=True)
    ap.add_argument("--port", type=int, default=22)
    ap.add_argument("--user", default="root")
    ap.add_argument("--password", default=os.environ.get("SSH_PASS", ""))
    ap.add_argument("--poll-min", type=float, default=1.0)
    args = ap.parse_args()

    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(args.host, port=args.port, username=args.user, password=args.password,
              timeout=30, look_for_keys=False, allow_agent=False)

    def run(cmd, timeout=300):
        _, o, e = c.exec_command(cmd, timeout=timeout)
        return o.read().decode(errors="replace").strip(), e.read().decode(errors="replace").strip()

    sftp = c.open_sftp()
    with sftp.open(REMOTE_ENV_SH, "w") as fh:
        fh.write(BUILD_SH)
    sftp.close()
    run(f"chmod +x {REMOTE_ENV_SH}")

    print("启动环境构建（conda create gamd: ambertools + openbabel + parmed + mdtraj）…")
    run(f"cd /root && setsid bash -c '{REMOTE_ENV_SH} > {REMOTE_LOG} 2>&1' < /dev/null & echo LAUNCHED", timeout=60)
    poll(c, REMOTE_LOG, "env", max(15, int(args.poll_min * 60)))

    checks = ["tleap", "antechamber", "parmchk2", "cpptraj", "MMPBSA.py", "pmemd.cuda", "obabel"]
    found = {}
    for tool in checks:
        out, _ = run(f"source /root/miniconda3/etc/profile.d/conda.sh; conda activate gamd; which {tool} 2>/dev/null || echo MISSING")
        found[tool] = out.strip().splitlines()[-1] if out.strip() else "MISSING"
    ver, _ = run("source /root/miniconda3/etc/profile.d/conda.sh; conda activate gamd; python --version")
    found["python"] = ver.strip()

    print("\n=== 校验结果 ===")
    for k, v in found.items():
        flag = "OK  " if v != "MISSING" else "MISS"
        print(f"  [{flag}] {k:<12} {v}")

    has_pmemd_cuda = found.get("pmemd.cuda", "MISSING") != "MISSING"
    summary = {
        "host": f"{args.host}:{args.port}",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "tools": found,
        "amber_gpu_md_available": has_pmemd_cuda,
        "plan": ("AmberTools 单栈：tleap → pmemd.cuda → cpptraj → MMPBSA.py（不编译 GROMACS）"
                 if has_pmemd_cuda else
                 "pmemd.cuda 不可用 → 需另寻 GPU MD 引擎（另开 GROMACS GPU 版或换镜像）"),
    }
    out_json = RESULTS / "15_env_build.json"
    out_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n写入 {out_json}")
    c.close()
    return 0 if has_pmemd_cuda else 2


if __name__ == "__main__":
    raise SystemExit(main())
