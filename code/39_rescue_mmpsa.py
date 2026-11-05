#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""39_rescue_mmpsa.py — 给「过阳性对照下限但被 500 帧硬门刷掉」的副本补算 ΔH_PB

为什么做这件事
--------------
batch 2 里 CHEMBL848 三个副本的 pose 保留率分别是 0.363 / 0.381 / 0.303，
全部高于阳性对照（LxxLL 肽）的下限 0.215，却因为没到 500/1001 帧而**从未跑
MM-PBSA**，因此表里 ΔH 是 NA。而阳性对照自己也有 1/3 副本低于 500 帧。
这个硬门的事后暴露说明判据需要校准，而那些副本的能量**可以**补算。

为什么不用重跑 MD
-----------------
管线只是 ``n_kept < MIN_KEEP_FRAMES`` 时跳过了第 7 步（cpptraj 抽帧）和
第 8 步（MMPBSA.py）；这两步的输入都完整留在副本目录里：
    com_nobox.nc / com_nobox.prmtop / rec_nobox.prmtop / lig_nobox.prmtop
所以本脚本在远端只补这两步，并**直接复用 12_pipeline_v3.py 里的常量
(LIG_RMSD_CUTOFF, MMPBSA_IN) 和严格解析函数 parse_dg_bind_strict**，
保证补算值与正式批次数值口径完全一致。

用法
----
    python code/39_rescue_mmpsa.py deploy     # 选目标 + 上传 worker + 启动
    python code/39_rescue_mmpsa.py status     # 查进度
    python code/39_rescue_mmpsa.py collect    # 拉取结果 -> results/39_rescue_deltaH.csv
    python code/39_rescue_mmpsa.py wait       # 等待并自动 collect（后台用）
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
RES = REBUILD / "results"
F_FULL = RES / "36_screen_full.csv"
OUT_CSV = RES / "39_rescue_deltaH.csv"

HOST = os.environ.get("MATPOOL_HOST", "hz-4.matpool.com")
PORT = int(os.environ.get("MATPOOL_PORT", "26255"))
USER = os.environ.get("MATPOOL_USER", "root")
PASS = os.environ.get("MATPOOL_PASS", "wmKG=3YRsZaHmij0")
BASE = "/root/autodl-tmp/HPV_structure/candidates_run"
N_PAR = 8        # 并行度（32 核机器，单个 MMPBSA.PB 单线程）

POS_CTRL = "LxxLL_peptide"

WORKER = '''#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""rescue_mmpsa_worker.py — 单个副本补算 ΔH_PB（第 7、8 步）"""
import importlib.util, json, os, subprocess, sys

BASE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("p3", os.path.join(BASE, "12_pipeline_v3.py"))
p3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p3)


def main() -> int:
    lig, rep = sys.argv[1].split(":")
    d = os.path.join(BASE, lig, "rep%s" % rep)
    out = {"ligand": lig, "replica": int(rep), "dir": d, "ok": False}
    try:
        # ---- 第 7 步：帧过滤（判据与主管线同一常量）----
        keep, n_all = [], 0
        with open(os.path.join(d, "lig_rmsd.dat"), errors="replace") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.split()
                n_all += 1
                if float(cols[1]) < p3.LIG_RMSD_CUTOFF:
                    keep.append(int(cols[0]))
        out["n_frames"], out["kept_frames"] = n_all, len(keep)
        if not keep:
            raise RuntimeError("no kept frame at RMSD<%s" % p3.LIG_RMSD_CUTOFF)
        fin = os.path.join(d, "filt_rescue.in")
        with open(fin, "w") as fh:
            fh.write("parm %s\\ntrajin %s\\ntrajout %s onlyframes %s\\n" % (
                os.path.join(d, "com_nobox.prmtop"),
                os.path.join(d, "com_nobox.nc"),
                os.path.join(d, "com_nobox_filt_rescue.nc"),
                ",".join(map(str, keep))))
        rc = subprocess.run("cpptraj -i %s > %s 2>&1" % (fin, os.path.join(d, "cpptraj_filt_rescue.log")),
                            shell=True).returncode
        if rc != 0 or not os.path.exists(os.path.join(d, "com_nobox_filt_rescue.nc")):
            raise RuntimeError("cpptraj 抽帧失败 rc=%s" % rc)

        # ---- 第 8 步：MMPBSA.py PB（cwd=d，避免临时文件互踩 P7）----
        minp = os.path.join(d, "mmpbsa_rescue.in")
        with open(minp, "w") as fh:
            fh.write(p3.MMPBSA_IN)
        cmd = ("MMPBSA.py -O -i %s -o %s -cp %s -rp %s -lp %s -y %s > %s 2>&1" % (
            minp, os.path.join(d, "FINAL_RESULTS_RESCUE.dat"),
            os.path.join(d, "com_nobox.prmtop"),
            os.path.join(d, "rec_nobox.prmtop"),
            os.path.join(d, "lig_nobox.prmtop"),
            os.path.join(d, "com_nobox_filt_rescue.nc"),
            os.path.join(d, "mmpbsa_rescue.log")))
        rc = subprocess.run(cmd, shell=True, cwd=d).returncode
        fres = os.path.join(d, "FINAL_RESULTS_RESCUE.dat")
        if rc != 0 or not os.path.exists(fres):
            raise RuntimeError("MMPBSA 失败 rc=%s" % rc)

        # ---- 严格解析（与主管线同一函数）----
        txt = open(fres, errors="replace").read()
        dH, gas, solv, comps = p3.parse_dg_bind_strict(txt)
        if dH is None:
            raise RuntimeError("严格解析失败：缺 Differences 段")
        out.update(ok=True, delta_H_PB_kcal_mol=dH, dG_gas=gas, dG_solv=solv)
    except Exception as e:
        out["error"] = str(e)[:300]
    outdir = os.path.join(BASE, "rescue_out")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "%s_rep%s.json" % (out["ligand"], out["replica"])), "w") as fh:
        json.dump(out, fh, ensure_ascii=False)
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
'''


def targets() -> list[tuple[str, str]]:
    """读取本地汇总，挑出「Ks ≥ 阳性对照下限 且 ΔH 缺失」的副本。"""
    rows = list(csv.DictReader(open(F_FULL, newline="")))
    floor = min(float(r["kept_frac"]) for r in rows
                if r["CHEMBL"] == POS_CTRL and r["kept_frac"])
    out = []
    for r in rows:
        if not r["kept_frac"]:
            continue
        kf = float(r["kept_frac"])
        dh = (r.get("delta_H_PB_kcal_mol") or "").strip()
        if kf >= floor and dh in ("", "NA"):
            out.append((r["CHEMBL"], str(r["replica"]), kf))
    out.sort(key=lambda x: -x[2])
    return [(a, b) for a, b, _ in out], floor


def cli():
    import paramiko
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, port=PORT, username=USER, password=PASS, timeout=25)
    return c


def run(c, cmd: str, timeout: int = 120):
    _, so, se = c.exec_command(cmd, timeout=timeout)
    return so.read().decode(), se.read().decode()


def run_bg(c, cmd: str) -> None:
    """启动完全脱离终端的后台任务。

    坑：``nohup ... &`` 后若在**同一 channel** 里读 stdout，paramiko 会一直阻塞
    到任务结束（此处表现为 PipeTimeout）。正确做法是 setsid + 重定向 + < /dev/null，
    并且**不读**这个 channel 的输出——想看日志就另开一条命令 tail。
    """
    c.exec_command(f"cd {BASE} && setsid nohup {cmd} < /dev/null > /dev/null 2>&1 &",
                   timeout=10)


def deploy() -> int:
    tgt, floor = targets()
    print(f"阳性对照下限 = {floor:.3f}；待补算副本 {len(tgt)} 个：")
    for lig, rep in tgt:
        print(f"  {lig} rep{rep}")
    c = cli()
    sftp = c.open_sftp()
    with sftp.open(f"{BASE}/rescue_mmpsa_worker.py", "w") as fh:
        fh.write(WORKER)
    with sftp.open(f"{BASE}/rescue_targets.txt", "w") as fh:
        fh.write("\n".join(f"{a}:{b}" for a, b in tgt) + "\n")
    runner = (
        "#!/bin/bash\n"
        "source /root/miniconda3/etc/profile.d/conda.sh\n"
        "conda activate gamd\n"
        f"cd {BASE}\n"
        "mkdir -p rescue_out\n"
        f"cat rescue_targets.txt | xargs -P {N_PAR} -I{{}} python rescue_mmpsa_worker.py {{}}\n"
        "echo RESCUE_DONE\n"
    )
    with sftp.open(f"{BASE}/run_rescue.sh", "w") as fh:
        fh.write(runner)
    sftp.close()
    out, _ = run(c, f"pgrep -fc 'rescue_mmpsa_worker' || true")
    if int((out.strip() or "0").splitlines()[0]) > 0:
        print(f"[deploy] 已有 {out.strip()} 个 worker 在跑，跳过重复启动")
        c.close()
        return 0
    run_bg(c, f"bash {BASE}/run_rescue.sh > {BASE}/rescue.log 2>&1")
    time.sleep(5)
    n = run(c, f"pgrep -fc 'rescue_mmpsa_worker' || true")[0].strip()
    print(f"[deploy] 启动完成，worker 进程数 = {n}")
    c.close()
    return 0


def status() -> int:
    c = cli()
    out, _ = run(c, f"ls {BASE}/rescue_out/*.json 2>/dev/null | wc -l; "
                    f"cat {BASE}/rescue_targets.txt | wc -l; "
                    f"tail -2 {BASE}/rescue.log; "
                    f"grep -l '\"ok\": false' {BASE}/rescue_out/*.json 2>/dev/null | head")
    print(out)
    c.close()
    return 0


def collect() -> int:
    """拉取远端 rescue_out JSON -> results/39_rescue_deltaH.csv"""
    c = cli()
    sftp = c.open_sftp()
    try:
        names = [n for n in sftp.listdir(f"{BASE}/rescue_out") if n.endswith(".json")]
    except Exception as e:
        print(f"[FAIL] 读取远端结果失败：{e}")
        c.close()
        return 1
    local = RES / "rescue_out"
    local.mkdir(exist_ok=True)
    for n in names:
        sftp.get(f"{BASE}/rescue_out/{n}", str(local / n))
    sftp.close()
    c.close()
    rows = []
    for n in sorted(names):
        rows.append(json.load(open(local / n, encoding="utf-8")))
    rows.sort(key=lambda r: (r["ligand"], r["replica"]))
    cols = ["ligand", "replica", "kept_frames", "n_frames", "kept_frac",
            "delta_H_PB_kcal_mol", "dG_gas", "dG_solv", "status", "error"]
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            nf, kf = r.get("n_frames"), r.get("kept_frames")
            r["kept_frac"] = round(kf / nf, 4) if nf else "NA"
            r["status"] = "OK" if r.get("ok") else "FAILED"
            w.writerow(r)
    ok = sum(1 for r in rows if r.get("ok"))
    print(f"[collect] {ok}/{len(rows)} 个副本补算成功 -> {OUT_CSV}")
    for r in rows:
        if r.get("ok"):
            print(f"  {r['ligand']:20s} rep{r['replica']}  "
                  f"kept {r['kept_frames']}/{r['n_frames']} ({r['kept_frac']:.3f})  "
                  f"ΔH={r['delta_H_PB_kcal_mol']}")
        else:
            print(f"  {r['ligand']:20s} rep{r['replica']}  FAILED: {r.get('error')}")
    return 0


def wait_and_collect(interval: int = 120, max_rounds: int = 90) -> int:
    c = cli()
    out, _ = run(c, f"cat {BASE}/rescue_targets.txt | wc -l")
    n_target = int(out.strip().splitlines()[0])
    c.close()
    print(f"[{datetime.now():%H:%M:%S}] 目标 {n_target} 个，开始轮询")
    for i in range(1, max_rounds + 1):
        try:
            c = cli()
            out, _ = run(c, f"ls {BASE}/rescue_out/*.json 2>/dev/null | wc -l; "
                            f"tail -1 {BASE}/rescue.log")
            lines = out.strip().splitlines()
            done = int(lines[0]) if lines and lines[0].isdigit() else 0
            tail = lines[1] if len(lines) > 1 else ""
            c.close()
            print(f"[{datetime.now():%H:%M:%S}] round {i}: {done}/{n_target} {tail[:60]}")
            if done >= n_target:
                return collect()
        except Exception as e:
            print(f"[WARN] 轮询异常：{e}")
        time.sleep(interval)
    print("达到最大轮数，未完成")
    return 1


def main() -> int:
    cmd = sys.argv[1] if len(sys.argv) > 1 else "deploy"
    return {"deploy": deploy, "status": status, "collect": collect,
            "wait": wait_and_collect}[cmd]()


if __name__ == "__main__":
    raise SystemExit(main())
