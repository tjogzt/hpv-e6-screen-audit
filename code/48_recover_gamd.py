#!/usr/bin/env python3
"""48_recover_gamd.py — 云端实例恢复后的一键收尾（GaMD v2 / Figure 5 / §3.9）。

背景：MatPool 实例停机或重启后 SSH 端口会轮换。原先所有脚本把
`hz-4.matpool.com:26255` 写死，端口一变全线失联。现在统一走环境变量：
    HPV_SSH_HOST / HPV_SSH_PORT / HPV_SSH_USER / HPV_SSH_PW

用法：
    HPV_SSH_PORT=<新端口> python3 code/48_recover_gamd.py status
    HPV_SSH_PORT=<新端口> python3 code/48_recover_gamd.py all

status：只读，检查 9 条轨迹 + 分析产物的完整度。
all    ：缺分析就远端补跑并等待 → 拉回 → 43 重加权 → 44 Figure 5 → 45 §3.9。
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import paramiko

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
GBASE = "/root/autodl-tmp/HPV_structure/candidates_run/gamd_v2"
JOBS = [f"{s}_rep{r}" for s in
        ("CHEMBL848_pose1", "CHEMBL1129_pose1", "LxxLL_peptide") for r in (1, 2, 3)]
ENV = dict(os.environ, GAMD_TAG="gamd_v2")
PY = sys.executable


def cli():
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(os.environ.get("HPV_SSH_HOST", "hz-4.matpool.com"),
              port=int(os.environ.get("HPV_SSH_PORT", "26255")),
              username=os.environ.get("HPV_SSH_USER", "root"),
              password=os.environ.get("HPV_SSH_PW", "wmKG=3YRsZaHmij0"),
              timeout=30)
    return c


def sh(c, cmd, timeout=180):
    _, o, e = c.exec_command(cmd, timeout=timeout)
    return o.read().decode(errors="replace").strip()


def status(c):
    out = sh(c, f"""cd {GBASE} 2>/dev/null || {{ echo 'NO_DIR'; exit 0; }}
for j in {' '.join(JOBS)}; do
  d={GBASE}/$j
  n=$(grep -oE '^ *[0-9]+ +[0-9.]+$' $d/gamd.log 2>/dev/null|tail -1|awk '{{print $1}}')
  fin=$(grep -c 'Finished mdrun' $d/gamd.log 2>/dev/null)
  f=$( [ -s $d/frames.csv ] && tail -n +2 $d/frames.csv | wc -l || echo 0 )
  rc=$( [ -s $d/rescontact.csv ] && echo Y || echo N )
  rw=$( [ -s $d/gamd-reweight.dat ] && echo Y || echo N )
  echo "$j step=${{n:-0}} ns=$(awk -v x=${{n:-0}} 'BEGIN{{printf "%.0f", x*0.002/1000}}') fin=$fin frames=$f bitmap=$rc reweight=$rw"
done""")
    print("# GaMD v2 状态")
    print(out)
    return out


def run_local(args, desc):
    print(f"== {desc} ==", flush=True)
    r = subprocess.run(args, cwd=REBUILD, env=ENV, capture_output=True, text=True,
                       timeout=7200)
    print((r.stdout or "")[-1200:])
    if r.returncode != 0:
        print("STDERR:", (r.stderr or "")[-1200:])
    return r.returncode


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "status"
    try:
        c = cli()
    except Exception as ex:
        print(f"❌ 连不上云端：{ex}")
        print("   去 MatPool 控制台启动实例并核对最新 host:port，然后：")
        print("   HPV_SSH_PORT=<新端口> python3 code/48_recover_gamd.py status")
        return 1
    try:
        s = status(c)
    finally:
        c.close()
    if mode != "all":
        return 0
    if "NO_DIR" in s:
        print("❌ 远端 gamd_v2 目录不存在——数据盘可能已释放，需重建。")
        return 1

    need_analysis = any(("frames=0" in ln) or ("bitmap=N" in ln)
                        for ln in s.splitlines() if "_rep" in ln)
    if need_analysis:
        if run_local([PY, str(HERE / "42_gamd_analyse.py"), "run"], "远端补跑分析"):
            return 1
        # 等 9 条都产出逐帧接触位图（不依赖 runner 标记，避免旧标记误判）
        c = cli()
        try:
            for k in range(120):                 # 上限 2 h
                time.sleep(60)
                s = status(c)
                ok = all(("frames=0" not in ln) and ("bitmap=Y" in ln)
                         for ln in s.splitlines() if "_rep" in ln)
                print(f"  [{k + 1} min] 9 条分析完成 = {ok}", flush=True)
                if ok:
                    break
            else:
                print("❌ 分析超时")
                return 1
        finally:
            c.close()

    for args, desc in (
        ([PY, str(HERE / "42_gamd_analyse.py"), "pull"], "拉回分析产物"),
        ([PY, str(HERE / "43_gamd_reweight.py")], "重加权 + 汇总"),
        (["Rscript", str(HERE / "44_fig5_gamd.R"),
          str(REBUILD / "results" / "gamd_v2"),
          str(REBUILD / "results" / "Fig5_gamd_contact_pmf_gamd_v2.pdf")], "Figure 5"),
        ([PY, str(HERE / "45_manuscript_39_block.py")], "§3.9 文本块"),
    ):
        if run_local(args, desc):
            print(f"❌ {desc} 失败")
            return 1
    print("✅ 收尾完成")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
