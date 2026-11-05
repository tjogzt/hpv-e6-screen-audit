#!/usr/bin/env python3
"""46_watch_gamd.py — GaMD 收尾守望：完成 → 分析 → 拉取 → 重加权 → 图 5 → §3.9 块。

每 15 min 查一次远端。9 个作业（runner_gpu0..3 的 log 里出现 ALLDONE_GPU*）
全部结束后，自动执行：
    42 deploy / run   远端轨迹分析（trjconv -pbc mol + worker）
    42 pull           拉回 frames.csv / resocc.csv / gamd-reweight.dat / gamd-para.dat
    43                重加权 + 汇总量
    44                Figure 5 (PDF)
    45                §3.9 文本块
然后写哨兵 results/gamd_finalized.done，日志 results/gamd_watch.log。

用法：python3 code/46_watch_gamd.py
"""
from __future__ import annotations

import subprocess
import os
import sys
import time
from pathlib import Path

import paramiko

HOST, PORT, USER, PW = "hz-4.matpool.com", 26255, "root", "wmKG=3YRsZaHmij0"
# GAMD_TAG：gamd（第一版 σ0=6）/ gamd_v2（标定后可重加权版）
TAG = os.environ.get("GAMD_TAG", "gamd").strip() or "gamd"
GBASE = f"/root/autodl-tmp/HPV_structure/candidates_run/{TAG}"
RSUFFIX = "_v2" if TAG != "gamd" else ""     # launch2 用 *_v2.sh / *_v2.log
NGPU = int(os.environ.get("GAMD_NGPU", "4"))  # 本轮实际占用的 GPU 数

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
LOG = REBUILD / "results" / f"gamd_watch_{TAG}.log"
DONE = REBUILD / "results" / f"gamd_finalized_{TAG}.done"
CODE = HERE
ENV = dict(os.environ, GAMD_TAG=TAG)

SYSTEMS = ["CHEMBL848_pose1", "CHEMBL1129_pose1", "LxxLL_peptide"]
REPLICAS = [1, 2, 3]


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with LOG.open("a") as fh:
        fh.write(line + "\n")


def cli() -> paramiko.SSHClient:
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, port=PORT, username=USER, password=PW, timeout=30)
    return c


def sh(c, cmd, timeout=120):
    _, o, e = c.exec_command(cmd, timeout=timeout)
    return o.read().decode(errors="replace").strip()


def run_local(args, desc):
    log(f"== {desc} ==")
    r = subprocess.run(args, cwd=REBUILD, capture_output=True, text=True,
                       timeout=7200, env=ENV)
    tail = (r.stdout or "")[-1500:]
    log(f"  rc={r.returncode}\n{tail}")
    if r.returncode != 0:
        log("  STDERR:\n" + (r.stderr or "")[-1500:])
    return r.returncode


def main() -> int:
    if DONE.exists():
        log("哨兵已存在，退出。")
        return 0
    while True:
        try:
            c = cli()
            alldone = sh(c, f"cat {GBASE}/runner_gpu*{RSUFFIX}.log 2>/dev/null | "
                            f"grep -c ALLDONE_GPU")
            try:
                n_done = int(alldone.splitlines()[-1])
            except ValueError:
                n_done = 0
            # 各作业步数
            steps = sh(c, f"for d in {GBASE}/*_rep*; do grep -oE '^ *[0-9]+ +[0-9.]+$' "
                          f"$d/gamd.log 2>/dev/null | tail -1 | awk '{{print $1}}'; done | "
                          f"paste -sd+ - | bc 2>/dev/null", 120)
            log(f"GPU 队列完成 {n_done}/{NGPU}；累计步数 {steps.strip()}")
            c.close()
            if n_done >= NGPU:
                break
        except Exception as ex:
            log(f"查询失败：{ex}；300 s 后重试")
        time.sleep(900)

    log("全部 GaMD 完成，开始分析收尾")
    # 用本解释器的绝对路径：nohup/非交互 shell 下 `python3` 可能解析到没有
    # paramiko 的 3.13（本机 zsh 里 python3 是 alias 到 3.10），会让子进程 ImportError。
    PY = sys.executable
    for args, desc in (
        ([PY, str(CODE / "42_gamd_analyse.py"), "deploy"], "上传分析 worker"),
        ([PY, str(CODE / "42_gamd_analyse.py"), "run"], "远端轨迹分析（后台启动）"),
    ):
        if run_local(args, desc):
            log(f"!! {desc} 失败")
            return 1

    # 等远端分析跑完（watch analyse.log 的 DONE 数）
    while True:
        try:
            c = cli()
            n = sh(c, f"grep -l DONE {GBASE}/*_rep*/analyse.log 2>/dev/null | wc -l", 120)
            runner = sh(c, f"grep -c ANALYSE_ALLDONE {GBASE}/runner_analyse.log 2>/dev/null", 60)
            log(f"分析完成 {n.splitlines()[-1]}/9；runner={runner.splitlines()[-1]}")
            c.close()
            if int(runner.splitlines()[-1] or 0) >= 1:
                break
        except Exception as ex:
            log(f"分析查询失败：{ex}")
        time.sleep(300)

    for args, desc in (
        ([PY, str(CODE / "42_gamd_analyse.py"), "pull"], "拉取分析结果"),
        ([PY, str(CODE / "43_gamd_reweight.py")], "重加权 + 汇总"),
        (["Rscript", str(CODE / "44_fig5_gamd.R"), str(REBUILD / "results" / TAG),
          str(REBUILD / "results" / f"Fig5_gamd_contact_pmf_{TAG}.pdf")], "Figure 5"),
        ([PY, str(CODE / "45_manuscript_39_block.py")], "§3.9 文本块"),
    ):
        if run_local(args, desc):
            log(f"!! {desc} 失败")
            return 1

    DONE.write_text(f"done {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    log("收尾完成。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
