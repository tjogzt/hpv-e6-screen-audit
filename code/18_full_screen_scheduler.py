#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
18_full_screen_scheduler.py — 14 配体 × 3 replica 全量筛选调度器（4 GPU 槽位）

设计
----
- 每个 job = (配体, replica)，4 个 GPU 槽位并发，每槽位一次只跑一个 job
- 跳过已有完成标记 SUMMARY_v3_{lig}_rep{r}.csv 的 job（pilot 结果以此标记免重跑）
- 每个 job 以显式 env 启动 12_pipeline_v3.py（HOME/PATH/LD_LIBRARY_PATH，见 P8）
- flock 防重复启动；全部完成后合并 SUMMARY_v3_all.csv
- AUTOSHUTDOWN=1 时全完成后自动关机（sync + shutdown + kill PID1 兜底）

环境变量
--------
BASE        默认 /root/autodl-tmp/HPV_structure/candidates_run
PROD_NS     默认 20
N_GPU       默认 4
NTOMP       默认 8（4 槽 × 8 = 32 核）
REPLICAS    默认 3
AUTOSHUTDOWN 默认 0；=1 全部完成后关机
"""
from __future__ import annotations

import csv
import datetime
import glob
import os
import subprocess
import sys
import time

BASE = os.environ.get("BASE", "/root/autodl-tmp/HPV_structure/candidates_run")
LIGDIR = os.environ.get("LIGDIR", os.path.join(BASE, "ligands_vina"))
PROD_NS = os.environ.get("PROD_NS", "20")
N_GPU = int(os.environ.get("N_GPU", "4"))
NTOMP = os.environ.get("NTOMP", "8")
REPLICAS = int(os.environ.get("REPLICAS", "3"))
AUTOSHUTDOWN = os.environ.get("AUTOSHUTDOWN", "0") == "1"
PIPELINE = os.path.join(BASE, "12_pipeline_v3.py")
LOG = os.path.join(BASE, "scheduler.log")
PY = "/root/miniconda3/envs/gamd/bin/python"


def log(msg: str) -> None:
    line = f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def marker(lig: str, rep: int) -> str:
    return os.path.join(BASE, f"SUMMARY_v3_{lig}_rep{rep}.csv")


def job_env(lig: str, rep: int, gpu: int) -> dict:
    env = os.environ.copy()
    env.update({
        "HOME": "/root",
        "PATH": "/root/miniconda3/envs/gamd/bin:/usr/local/cuda/bin:"
                "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
        "LD_LIBRARY_PATH": "/usr/local/cuda/lib64:" + env.get("LD_LIBRARY_PATH", ""),
        "BASE": BASE, "LIGDIR": LIGDIR,
        "CONDA_SH": "/root/miniconda3/etc/profile.d/conda.sh",
        "ENV_AMBER": "gamd", "ENV_GMX": "gamd",
        "LIG_ONLY": lig, "REPLICA_ID": str(rep), "REPLICAS": "1",
        "PROD_NS": PROD_NS, "GPU_ID": str(gpu), "NTOMP": NTOMP,
    })
    return env


def main() -> int:
    # flock 防重复启动
    import fcntl
    lock_fd = open(os.path.join(BASE, ".scheduler_lock"), "w")
    try:
        fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print("已有调度器在跑，退出", flush=True)
        return 1

    ligs = sorted(os.path.splitext(os.path.basename(p))[0]
                  for p in glob.glob(os.path.join(LIGDIR, "*.pdb")))
    # LIG_ORDER：逗号分隔的优先配体名，排到队列最前（用于先跑对照/回归实验）
    prio = [s.strip() for s in os.environ.get("LIG_ORDER", "").split(",") if s.strip()]
    prio = [p for p in prio if p in ligs]
    rest = [l for l in ligs if l not in prio]
    ligs = prio + rest
    if prio:
        log(f"LIG_ORDER 优先：{prio}")
    jobs = [(lig, rep) for lig in ligs for rep in range(1, REPLICAS + 1)]
    todo = [(l, r) for (l, r) in jobs if not os.path.exists(marker(l, r))]
    log(f"配体 {len(ligs)} 个，job 共 {len(jobs)} 个，待跑 {len(todo)} 个 "
        f"（N_GPU={N_GPU}, PROD_NS={PROD_NS}, NTOMP={NTOMP}）")

    slots: dict[int, tuple[subprocess.Popen, str, int]] = {}  # gpu -> (proc, lig, rep)
    results: dict[tuple[str, int], int] = {}
    while todo or slots:
        for gpu in range(N_GPU):
            if gpu not in slots and todo:
                lig, rep = todo.pop(0)
                jlog = open(os.path.join(BASE, f"job_{lig}_rep{rep}.log"), "w")
                p = subprocess.Popen([PY, PIPELINE], stdout=jlog, stderr=subprocess.STDOUT,
                                     cwd=BASE, env=job_env(lig, rep, gpu))
                slots[gpu] = (p, lig, rep)
                log(f"启动 {lig} rep{rep} -> GPU{gpu} (pid {p.pid})")
        for gpu in list(slots):
            p, lig, rep = slots[gpu]
            rc = p.poll()
            if rc is None:
                continue
            results[(lig, rep)] = rc
            ok = os.path.exists(marker(lig, rep))
            log(f"完成 {lig} rep{rep} (GPU{gpu}, rc={rc}, 标记{'有' if ok else '缺失!'})")
            del slots[gpu]
        time.sleep(60)

    # 合并汇总
    all_csv = os.path.join(BASE, "SUMMARY_v3_all.csv")
    rows = []
    header = None
    for f in sorted(glob.glob(os.path.join(BASE, "SUMMARY_v3_*_rep*.csv"))):
        with open(f, newline="") as fh:
            rd = csv.DictReader(fh)
            header = header or rd.fieldnames
            rows.extend(list(rd))
    if header:
        with open(all_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=header)
            w.writeheader()
            w.writerows(rows)
    log(f"全部完成：{len(results)} job，合并 {len(rows)} 行 -> {all_csv}")

    if AUTOSHUTDOWN:
        log("AUTOSHUTDOWN=1：sync + shutdown")
        subprocess.run(["sync"])
        subprocess.run(["bash", "-c",
                        "shutdown -h now 2>/dev/null; sleep 10; kill -TERM 1; sleep 5; kill -9 1"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
