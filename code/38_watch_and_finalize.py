#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""38_watch_and_finalize.py — 监视 batch 2，跑满 48 run 后自动收尾

流程：每 INTERVAL 秒 SSH 查远端 SUMMARY_v3_*.csv 数量 + 在跑进程数；
当 48 个全部就绪且无进程在跑时，依次执行：
    1) 36_pull_full_screen.py      拉取 + 对照校准汇总
    2) 37_manuscript_38_block.py   生成 §3.8 文本块
    3) Rscript 34_fig2_zn_pose.R + 35_fig3_startpose_deltaH.R  重画 Fig2/3
    4) 写 results/batch2_finalized.done 哨兵文件

日志：results/batch2_watch.log（追加）。
用法：
    /opt/homebrew/bin/python3.10 code/38_watch_and_finalize.py \
        --interval 900 --max-rounds 40
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
RES = REBUILD / "results"
LOG = RES / "batch2_watch.log"
DONE = RES / "batch2_finalized.done"

PY = "/opt/homebrew/bin/python3.10"
N_EXPECTED = 48

import os
HOST = os.environ.get("MATPOOL_HOST", "hz-4.matpool.com")
PORT = int(os.environ.get("MATPOOL_PORT", "26255"))
USER = os.environ.get("MATPOOL_USER", "root")
PASS = os.environ.get("MATPOOL_PASS", "wmKG=3YRsZaHmij0")
BASE = "/root/autodl-tmp/HPV_structure/candidates_run"


def log(msg: str) -> None:
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def remote_status() -> tuple[int, int] | None:
    """返回 (summary 数, 在跑进程数)；SSH 失败返回 None。"""
    import paramiko
    cli = paramiko.SSHClient()
    cli.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        cli.connect(HOST, port=PORT, username=USER, password=PASS, timeout=25)
        cmd = (
            f"ls {BASE}/SUMMARY_v3_*.csv 2>/dev/null | wc -l; "
            "ps -eo cmd | grep -cE '[1]2_pipeline_v3|[g]mx mdrun'"
        )
        _, so, _ = cli.exec_command(cmd, timeout=40)
        out = so.read().decode().split()
        return int(out[0]), int(out[1])
    except Exception as e:
        log(f"SSH 失败：{e}")
        return None
    finally:
        try:
            cli.close()
        except Exception:
            pass


def run(cmd: list[str], tag: str) -> bool:
    rc = subprocess.run(cmd, cwd=REBUILD, capture_output=True, text=True)
    tail = "\n".join((rc.stdout + rc.stderr).splitlines()[-8:])
    log(f"{tag} rc={rc.returncode}\n{tail}")
    return rc.returncode == 0


def finalize() -> bool:
    ok = True
    ok &= run([PY, "code/36_pull_full_screen.py"], "36 拉取+汇总")
    ok &= run([PY, "code/37_manuscript_38_block.py"], "37 §3.8 文本块")
    ok &= run(["Rscript", "code/34_fig2_zn_pose.R"], "34 Figure 2")
    ok &= run(["Rscript", "code/35_fig3_startpose_deltaH.R"], "35 Figure 3")
    if ok:
        with open(DONE, "w") as fh:
            fh.write(f"finalized at {datetime.now():%Y-%m-%d %H:%M:%S}\n")
        log("收尾完成 -> results/batch2_finalized.done")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--interval", type=int, default=900, help="秒，默认 900")
    ap.add_argument("--max-rounds", type=int, default=40, help="最大轮数，默认 40")
    args = ap.parse_args()

    if DONE.exists():
        log("哨兵已存在，直接退出")
        return 0
    log(f"开始监视：interval={args.interval}s max_rounds={args.max_rounds}")
    for i in range(1, args.max_rounds + 1):
        st = remote_status()
        if st is None:
            log(f"round {i}: 查询失败，下轮重试")
        else:
            n_sum, n_proc = st
            log(f"round {i}: summary={n_sum}/{N_EXPECTED} 在跑={n_proc}")
            if n_sum >= N_EXPECTED and n_proc == 0:
                log("48/48 就绪，开始收尾")
                return 0 if finalize() else 2
        if i < args.max_rounds:
            time.sleep(args.interval)
    log("达到最大轮数仍未跑满，退出（不收尾）")
    return 1


if __name__ == "__main__":
    sys.exit(main())
