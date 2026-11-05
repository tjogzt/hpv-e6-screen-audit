#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""14_deploy_pilot.py — 部署并监控 B 线 pilot（双 Zn 验收门槛）。

pilot 的目的（不是出结果，是过门槛）
------------------------------------
1. 验证 fallback Zn 参数（非键 Li/Merz）能否在 100 ns 内维持双 Zn 四配位
   —— 合格判据：Zn–S ≤ 2.8 Å、每个 Zn ≥4 个 S、两个 Zn 同时合格的帧 ≥95%
2. 实测吞吐量（ns/day），用于重排 Tier 2/3 的日程与预算
3. 实测 replica 间 ΔH 标准差 σ_rep 与自相关时间，用于预注册后续 replica 数
4. 验证 v2 管线（NPT 生产、严格解析器、配体电荷溯源、replica 种子）能跑通

用法
----
/opt/anaconda3/bin/python rebuild_20261025/code/14_deploy_pilot.py \
    --host <host> --port <port> --user root --password <pwd> \
    [--ligand CHEMBL46917] [--replicas 3] [--ns 100] [--gpus 0,1,2] [--poll 30]

先跑 --precheck 只做环境体检（不启动生产），确认 gmx/tleap/antechamber/MMPBSA.py/mdtraj 齐备。
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from pathlib import Path

try:
    import paramiko
except ImportError:
    sys.exit("缺少 paramiko：/opt/anaconda3/bin/python -m pip install paramiko")

ROOT = Path(__file__).resolve().parents[2]
REBUILD = Path(__file__).resolve().parents[1]
LOCAL_PIPE = REBUILD / "code" / "12_pipeline_v2.py"
LOG_CSV = REBUILD / "results" / "pipeline_run_log.csv"
REMOTE_BASE = "/root/autodl-tmp/HPV_structure/candidates_run"
REMOTE_PIPE = f"{REMOTE_BASE}/12_pipeline_v2.py"


def ssh_run(c, cmd, timeout=300):
    _, o, e = c.exec_command(cmd, timeout=timeout)
    return o.read().decode(errors="replace").strip(), e.read().decode(errors="replace").strip()


def precheck(c) -> bool:
    print("=== 环境体检 ===")
    checks = {
        "gmx": "source /root/autodl-tmp/miniforge3/etc/profile.d/conda.sh; conda activate mmpbsa; gmx --version 2>/dev/null | head -3",
        "tleap": "source /root/autodl-tmp/miniforge3/etc/profile.d/conda.sh; conda activate gamd; which tleap antechamber parmchk2 cpptraj MMPBSA.py",
        "obabel": "which obabel",
        "mdtraj": "source /root/autodl-tmp/miniforge3/etc/profile.d/conda.sh; conda activate gamd; python -c 'import mdtraj, numpy; print(mdtraj.version.version)'",
        "gpu": "nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader",
    }
    all_ok = True
    for key, cmd in checks.items():
        out, err = ssh_run(c, cmd, timeout=180)
        ok = bool(out.strip())
        all_ok &= ok
        print(f"  [{'OK' if ok else 'MISSING'}] {key}: {out.splitlines()[0] if out else (err[:120] or '未找到')}")
        if out:
            for line in out.splitlines()[1:]:
                print(f"         {line}")
    if not all_ok:
        print("\n缺失组件的补救顺序：")
        print("  1) mdtraj → conda activate gamd; pip install mdtraj（Zn QC 依赖，缺则 QC 记为 SKIPPED）")
        print("  2) gmx/tleap → 需先建环境（miniforge + gamd/mmpbsa env），见 docs/gpu_requirements.md 第 4 节")
    return all_ok


def launch(c, ligand: str, replicas: int, ns: float, gpus: list[int]) -> None:
    for r in range(replicas):
        one = (f"cd {REMOTE_BASE} && setsid bash -c '"
               f"GPU_ID={gpus[r % len(gpus)]} LIGDIR={REMOTE_BASE}/ligands_vina "
               f"LIG_ONLY={ligand} REPLICAS=1 PROD_NS={ns} SAVE_PS=20 "
               f"python {REMOTE_PIPE} > {REMOTE_BASE}/pilot_rep{r+1}.log 2>&1' "
               f"< /dev/null & echo LAUNCHED_REP{r+1}_GPU{gpus[r % len(gpus)]}")
        out, err = ssh_run(c, one, timeout=60)
        print(f"  {out or err[:150]}")
        time.sleep(1)


def poll(c, every_min: int, max_hours: float) -> None:
    print(f"\n=== 监控（每 {every_min} 分钟，最长 {max_hours} 小时；Ctrl-C 退出监控不会停远端任务）===")
    t_end = time.time() + max_hours * 3600
    try:
        while time.time() < t_end:
            out, _ = ssh_run(c, f"tail -n 6 {REMOTE_BASE}/run_v2.log 2>/dev/null; "
                                f"echo '--- SUMMARY ---'; cat {REMOTE_BASE}/SUMMARY_v2.csv 2>/dev/null | cut -c1-160",
                             timeout=60)
            print(f"\n[{time.strftime('%H:%M:%S')}]\n{out}")
            if os.path.exists("/tmp/stop_pilot_poll"):
                break
            time.sleep(every_min * 60)
    except KeyboardInterrupt:
        print("\n监控已退出（远端任务继续运行）")


def main() -> int:
    ap = argparse.ArgumentParser(description="部署/监控 B 线 pilot")
    ap.add_argument("--host", required=True)
    ap.add_argument("--port", type=int, default=22)
    ap.add_argument("--user", default="root")
    ap.add_argument("--password", default=os.environ.get("SSH_PASS", ""))
    ap.add_argument("--ligand", default="CHEMBL46917")
    ap.add_argument("--replicas", type=int, default=3)
    ap.add_argument("--ns", type=float, default=100.0)
    ap.add_argument("--gpus", default="0,1,2")
    ap.add_argument("--precheck", action="store_true", help="只做环境体检，不启动生产")
    ap.add_argument("--poll", type=int, default=30, help="监控间隔（分钟）")
    ap.add_argument("--max-hours", type=float, default=24.0)
    args = ap.parse_args()
    gpus = [int(x) for x in args.gpus.split(",") if x.strip() != ""]

    if not args.password:
        sys.exit("未提供密码：--password 或环境变量 SSH_PASS")

    t0 = time.time()
    try:
        c = paramiko.SSHClient()
        c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        c.connect(args.host, port=args.port, username=args.user, password=args.password,
                  timeout=30, look_for_keys=False, allow_agent=False)
    except Exception as exc:  # noqa: BLE001
        print(f"[error] SSH 连接失败：{exc}")
        return 1

    try:
        ok = precheck(c)
        if args.precheck:
            return 0 if ok else 2
        if not ok:
            print("[warn] 环境体检有缺失，仍继续（失败会写入 SUMMARY_v2.csv）")

        sftp = c.open_sftp()
        sftp.put(str(LOCAL_PIPE), REMOTE_PIPE)
        sftp.close()
        print(f"\n已上传管线 v2 → {REMOTE_PIPE}")

        # mdtraj 缺失时自动补装（Zn QC 依赖）
        out, err = ssh_run(c, "source /root/autodl-tmp/miniforge3/etc/profile.d/conda.sh; "
                              "conda activate gamd; python -c 'import mdtraj' 2>/dev/null || "
                              "pip install -q mdtraj", timeout=600)
        print(f"mdtraj 检查/安装：{out or '已处理'}")

        print(f"\n启动 pilot：{args.ligand} × {args.replicas} replica × {args.ns} ns，GPU {gpus}")
        launch(c, args.ligand, args.replicas, args.ns, gpus)
        poll(c, args.poll, args.max_hours)
    finally:
        c.close()

    LOG_CSV.parent.mkdir(parents=True, exist_ok=True)
    new = not LOG_CSV.exists()
    with open(LOG_CSV, "a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["timestamp", "script", "host", "uploaded", "verified", "failed", "seconds"])
        w.writerow([time.strftime("%Y-%m-%d %H:%M:%S"), "14_deploy_pilot.py",
                    args.host, 1, 1 if not args.precheck else 0, 0, round(time.time() - t0, 1)])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
