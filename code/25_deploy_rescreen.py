#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
25_deploy_rescreen.py — 修复 Zn 模型后的全量重跑部署（run2）

与 run1 的差别（provenance 必须可区分）：
  - ZN_MODEL 由环境变量指定（默认 zaff_cym = bonded ZAFF + CYM），受体换成 rec_cym.pdb，
    两个 Zn 改写为同链(C) 残基 144/145 并与其 4 个配位 Cys 成键。
  - ZN_S_CUTOFF 已修正为 0.28 nm（run1 的 2.8 被当作 Å 实际是 28 Å，门从未生效）。
  - run1 的产物（SUMMARY_v3_*_rep*.csv + CHEMBL*_pose1/ 工作目录）整体归档到
    archive_run1_nonbond_cys/，绝不与 run2 混写。

用法（本机执行，实例须已开机）：
    ZN_MODEL=zaff_cym python3 25_deploy_rescreen.py
    ZN_MODEL=zaff_cym AUTOSHUTDOWN=1 python3 25_deploy_rescreen.py   # 跑完自动关机
"""
from __future__ import annotations

import os
import sys
import time

import paramiko

HOST = os.environ.get("MATPOOL_HOST", "hz-4.matpool.com")
PORT = int(os.environ.get("MATPOOL_PORT", "26255"))
PASS = os.environ.get("MATPOOL_PASS", "wmKG=3YRsZaHmij0")
BASE = "/root/autodl-tmp/HPV_structure/candidates_run"
HERE = os.path.dirname(os.path.abspath(__file__))
PY = "/root/miniconda3/envs/gamd/bin/python"

ZN_MODEL = os.environ.get("ZN_MODEL", "zaff_cym")
PROD_NS = os.environ.get("PROD_NS", "20")
N_GPU = os.environ.get("N_GPU", "4")
NTOMP = os.environ.get("NTOMP", "8")
REPLICAS = os.environ.get("REPLICAS", "3")
AUTOSHUTDOWN = os.environ.get("AUTOSHUTDOWN", "0")


def conn():
    for i in range(5):
        try:
            c = paramiko.SSHClient()
            c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            c.connect(HOST, port=PORT, username="root", password=PASS,
                      timeout=25, banner_timeout=30)
            return c
        except Exception as ex:  # noqa: BLE001
            print(f"connect retry {i}: {type(ex).__name__} {str(ex)[:60]}")
            time.sleep(5)
    sys.exit("SSH 连接失败：检查实例是否已开机、端口是否变化")


def main() -> int:
    c = conn()

    def run(cmd: str, t: int = 120) -> str:
        _, o, e = c.exec_command(cmd, timeout=t)
        return o.read().decode(errors="replace") + e.read().decode(errors="replace")

    # 0) 前置检查：Zn 验证结果必须存在，且选定模型 PASS
    sel = f"{BASE}/zn_model_test/{ZN_MODEL}/zn_result.json"
    out = run(f"cat {sel} 2>/dev/null | head -20")
    if '"verdict"' not in out:
        sys.exit(f"[中止] 未找到 Zn 模型验证结果 {sel}；请先跑 24_zn_model_test.py 并确认 PASS")
    print("== Zn 模型验证结果 ==")
    print(out)
    if '"verdict": "FAIL"' in out:
        sys.exit(f"[中止] ZN_MODEL={ZN_MODEL} 验证未通过（FAIL），不得用于全量重跑")

    # 1) 上传最新管线
    sftp = c.open_sftp()
    for fn in ("12_pipeline_v3.py", "18_full_screen_scheduler.py", "21_aggregate_screen.py"):
        sftp.put(os.path.join(HERE, fn), os.path.join(BASE, fn))
        print("上传", fn)
    sftp.close()
    print(run(f"{PY} -m py_compile {BASE}/12_pipeline_v3.py "
              f"{BASE}/18_full_screen_scheduler.py && echo COMPILE_OK"))

    # 2) 归档 run1（只移工作目录与 SUMMARY，templates/ligands_vina 保留复用）
    print(run(f"cd {BASE} && mkdir -p archive_run1_nonbond_cys && "
              f"mv SUMMARY_v3_*.csv archive_run1_nonbond_cys/ 2>/dev/null; "
              f"mv CHEMBL*_pose1 archive_run1_nonbond_cys/ 2>/dev/null; "
              f"mv SUMMARY_v3.csv archive_run1_nonbond_cys/ 2>/dev/null; "
              f"ls archive_run1_nonbond_cys | wc -l"))

    # 3) 启动调度器（显式 env，含 ZN_MODEL）
    env = (f"HOME=/root PATH=/root/miniconda3/envs/gamd/bin:/usr/local/cuda/bin:"
           "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin "
           "LD_LIBRARY_PATH=/usr/local/cuda/lib64 "
           f"BASE={BASE} LIGDIR={BASE}/ligands_vina "
           "CONDA_SH=/root/miniconda3/etc/profile.d/conda.sh "
           "ENV_AMBER=gamd ENV_GMX=gamd "
           f"PROD_NS={PROD_NS} N_GPU={N_GPU} NTOMP={NTOMP} REPLICAS={REPLICAS} "
           f"ZN_MODEL={ZN_MODEL} AUTOSHUTDOWN={AUTOSHUTDOWN}")
    launch = (f"cd {BASE} && setsid nohup env {env} {PY} {BASE}/18_full_screen_scheduler.py "
              f"> {BASE}/scheduler_run2.log 2>&1 < /dev/null & sleep 15; "
              f"ps -ef | grep -c '[1]8_full_screen'")
    print(run(launch, t=90))

    print("\n== 已启动 ==")
    print(f"  ZN_MODEL={ZN_MODEL}  PROD_NS={PROD_NS}  REPLICAS={REPLICAS}  "
          f"N_GPU={N_GPU}  NTOMP={NTOMP}")
    print(f"  日志：{BASE}/scheduler_run2.log  /  {BASE}/scheduler.log")
    print(f"  旧批次已归档：{BASE}/archive_run1_nonbond_cys/")
    c.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
