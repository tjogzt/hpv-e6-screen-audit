#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""17_deploy_pilot_v3.py — 用 Amber 单栈管线部署并监控 pilot。

pilot 目的（不是出结果，是过门槛）
----------------------------------
1. 验证 fallback Zn 参数能否维持双 Zn 四配位（Zn–S ≤2.8 Å、每 Zn ≥4 S、合格帧 ≥95%）
2. 实测吞吐 ns/day（A4000 比 4090 慢约 2.5–3 倍，据此重排 Tier 2/3 日程与预算）
3. 实测 replica 间 ΔH 标准差，供后续预注册 replica 数
4. 验证 v3 管线（NPT 生产、严格解析、配体电荷溯源、replica 种子、pmemd.cuda）跑通

默认 PROD_NS=20 —— 配位崩塌通常在最初几 ns 就暴露，先用短测省机时；
门槛通过后阶段 B 再拉长到 100 ns。

用法
----
/opt/anaconda3/bin/python rebuild_20261025/code/17_deploy_pilot_v3.py \
    --host <host> --port <port> --user root --password <pwd> [--ns 20] [--replicas 3] [--poll]

--check 只做环境体检；--ns/--replicas/--gpus 控制起法。
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

try:
    import paramiko
except ImportError:
    sys.exit("缺少 paramiko：/opt/anaconda3/bin/python -m pip install paramiko")

REBUILD = Path(__file__).resolve().parents[1]
LOCAL_PIPE = REBUILD / "code" / "16_pipeline_v3_amber.py"
REMOTE_BASE = "/root/autodl-tmp/HPV_structure/candidates_run"
REMOTE_PIPE = f"{REMOTE_BASE}/16_pipeline_v3_amber.py"
CONDA_SH = "/root/miniconda3/etc/profile.d/conda.sh"


def ssh_run(c, cmd, timeout=300):
    _, o, e = c.exec_command(cmd, timeout=timeout)
    return o.read().decode(errors="replace").strip(), e.read().decode(errors="replace").strip()


def check_env(c) -> bool:
    print("=== 环境体检 ===")
    out, _ = ssh_run(c, f"nvidia-smi --query-gpu=index,name,memory.used,memory.total,utilization.gpu "
                        f"--format=csv,noheader", timeout=90)
    print(out)
    need = ["tleap", "antechamber", "parmchk2", "cpptraj", "MMPBSA.py", "pmemd.cuda", "obabel"]
    ok = True
    for tool in need:
        out, _ = ssh_run(c, f"source {CONDA_SH}; conda activate gamd; which {tool} 2>/dev/null || echo MISSING",
                         timeout=120)
        miss = out.strip().endswith("MISSING") or not out.strip()
        ok &= not miss
        print(f"  [{'OK  ' if not miss else 'MISS'}] {tool:<12} {out.strip().splitlines()[-1] if out.strip() else ''}")
    out, _ = ssh_run(c, f"source {CONDA_SH}; conda activate gamd; "
                        f"python -c 'import mdtraj,numpy,parmed; print(mdtraj.__version__)' 2>&1 | tail -1",
                     timeout=120)
    print(f"  [{'OK  ' if 'Error' not in out and 'Traceback' not in out else 'MISS'}] mdtraj      {out}")
    out, _ = ssh_run(c, f"ls {REMOTE_BASE} 2>/dev/null; df -h / | tail -1", timeout=90)
    print(f"=== 数据与磁盘 ===\n{out}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", required=True)
    ap.add_argument("--port", type=int, default=22)
    ap.add_argument("--user", default="root")
    ap.add_argument("--password", default=os.environ.get("SSH_PASS", ""))
    ap.add_argument("--ligand", default="CHEMBL46917_pose1")
    ap.add_argument("--replicas", type=int, default=3)
    ap.add_argument("--ns", type=float, default=20.0)
    ap.add_argument("--gpus", default="0,1,2")
    ap.add_argument("--check", action="store_true", help="只体检，不启动")
    ap.add_argument("--poll", action="store_true", help="启动后轮询监控")
    ap.add_argument("--poll-sec", type=int, default=120)
    ap.add_argument("--max-min", type=float, default=180.0)
    args = ap.parse_args()
    gpus = [int(x) for x in args.gpus.split(",") if x.strip()]

    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(args.host, port=args.port, username=args.user, password=args.password,
              timeout=30, look_for_keys=False, allow_agent=False)

    ok = check_env(c)
    if args.check:
        c.close()
        return 0 if ok else 2
    if not ok:
        print("[warn] 环境有缺失，仍继续（失败会写进 SUMMARY_v3.csv）")

    sftp = c.open_sftp()
    sftp.put(str(LOCAL_PIPE), REMOTE_PIPE)
    sftp.close()
    print(f"\n已上传 v3 管线 → {REMOTE_PIPE}")

    for r in range(args.replicas):
        gpu = gpus[r % len(gpus)]
        cmd = (f"cd {REMOTE_BASE} && setsid bash -c '"
               f"GPU_ID={gpu} LIGDIR={REMOTE_BASE}/ligands_vina LIG_ONLY={args.ligand} "
               f"REPLICAS=1 PROD_NS={args.ns} "
               f"python {REMOTE_PIPE} > {REMOTE_BASE}/pilotv3_rep{r+1}.log 2>&1' "
               f"< /dev/null > /dev/null 2>&1 & sleep 1; echo LAUNCHED_REP{r+1}_GPU{gpu}")
        out, err = ssh_run(c, cmd, timeout=60)
        print(f"  {out or err[:150]}")

    if args.poll:
        print(f"\n=== 监控（每 {args.poll_sec}s，最长 {args.max_min} 分钟）===")
        t_end = time.time() + args.max_min * 60
        try:
            while time.time() < t_end:
                time.sleep(args.poll_sec)
                out, _ = ssh_run(c, f"tail -n 8 {REMOTE_BASE}/run_v3.log 2>/dev/null; "
                                    f"echo '--- SUMMARY ---'; "
                                    f"cut -c1-150 {REMOTE_BASE}/SUMMARY_v3.csv 2>/dev/null", timeout=90)
                print(f"\n[{time.strftime('%H:%M:%S')}]\n{out}")
        except KeyboardInterrupt:
            print("\n监控退出（远端任务继续）")
    c.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
