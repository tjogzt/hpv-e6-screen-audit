#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
19_deploy_full_screen.py — 全量筛选一键部署（MatPool 实例开机后在本机执行）

步骤：上传 12_pipeline_v3.py + 18_full_screen_scheduler.py → 启动调度器（4 卡并发）。
注意：MatPool 重启后端口/host 可能变化 → 用环境变量覆盖：
    MATPOOL_HOST / MATPOOL_PORT / MATPOOL_PASS
用法：
    python3 19_deploy_full_screen.py
"""
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


def main() -> None:
    c = conn()

    def run(cmd: str, t: int = 60) -> str:
        _, o, e = c.exec_command(cmd, timeout=t)
        return o.read().decode(errors="replace")

    # 1) 上传 v3 管线 + 调度器
    sftp = c.open_sftp()
    for fn in ("12_pipeline_v3.py", "18_full_screen_scheduler.py"):
        sftp.put(os.path.join(HERE, fn), os.path.join(BASE, fn))
        print("上传", fn)
    sftp.close()
    print("语法检查:", run(f"{PY} -m py_compile {BASE}/12_pipeline_v3.py "
                            f"{BASE}/18_full_screen_scheduler.py && echo OK"))

    # 2) 启动调度器（flock 防重复）
    env = (f"HOME=/root PATH={PY.rsplit('/',1)[0]}:/usr/local/cuda/bin:/usr/local/sbin:"
           "/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin "
           "LD_LIBRARY_PATH=/usr/local/cuda/lib64 "
           "PROD_NS=20 N_GPU=4 NTOMP=8 REPLICAS=3 AUTOSHUTDOWN=1")
    print(run(f"cd {BASE} && setsid nohup env {env} {PY} 18_full_screen_scheduler.py "
              f"> scheduler_nohup.log 2>&1 < /dev/null & echo launched"))
    time.sleep(20)
    print("调度器进程:", run("pgrep -af 18_full_screen_scheduler | grep -v pgrep | head -2"))
    print(run("sleep 5; head -5 /root/autodl-tmp/HPV_structure/candidates_run/scheduler.log "
              "2>/dev/null"))
    c.close()
    print("\n部署完成。监控：ssh 后 tail -f candidates_run/scheduler.log")


if __name__ == "__main__":
    main()
