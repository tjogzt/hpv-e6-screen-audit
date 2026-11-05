#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
11_seed_remote.py — B 线重建第一步：把本地资产重传到全新 AutoDL 实例。

用途
----
原 AutoDL 实例已不可恢复。本脚本把重建所需的全部本地资产上传至新实例的
/root/autodl-tmp/HPV_structure/candidates_run，并逐文件校验（size + md5）。

上传内容
--------
1) 模板   : rec_clean.pdb, zn1.pdb, zn2.pdb, ZN.mol2, zinc_nonbond_fallback.frcmod
2) 对照   : peptide_LxxLL.pdb (阳性), DMSO.mol2 (阴性)
3) 配体   : docking_results_full/4GIZ/CHEMBL*.pdb|.pdbqt|.smi
4) 管线   : .workbuddy/candidates_pipeline.py (已验证的 440 行管线)
5) 试点   : rebuild_20261025/data/pilot/CHEMBL46917_pose1.pdb

用法
----
/opt/anaconda3/bin/python rebuild_20261025/code/11_seed_remote.py \
    --host <host> --port <port> --user root --password <pwd>
# 只列清单不连接：加 --dry-run
# 密码也可由环境变量 SSH_PASS 提供（避免写进命令行历史）

设计约束（遵循研究规范）
------------------------
- 路径一律 pathlib / os.path.join，不用字符串拼斜杠
- 每个网络调用 try/except 包住，失败给降级指南
- 校验：本地 md5 与远端 md5 逐文件比对，打印覆盖率日志，不做静默丢弃
- 运行结束追加一行到 results/pipeline_run_log.csv
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import sys
import time
from pathlib import Path

try:
    import paramiko
except ImportError:
    sys.exit("缺少 paramiko：请先 /opt/anaconda3/bin/python -m pip install paramiko")

ROOT = Path(__file__).resolve().parents[2]           # HPV_structure/
REBUILD = ROOT / "rebuild_20261025"
LOG_CSV = REBUILD / "results" / "pipeline_run_log.csv"
REMOTE_BASE = "/root/autodl-tmp/HPV_structure/candidates_run"


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def collect_files() -> list[tuple[Path, str]]:
    """返回 [(本地路径, 远端相对子目录)]，按脚本中列出的顺序收集。"""
    items: list[tuple[Path, str]] = []

    # 1) 受体模板 + Zn 参数
    for name in ("rec_clean.pdb", "zn1.pdb", "zn2.pdb", "ZN.mol2"):
        p = ROOT / ".workbuddy" / name
        items.append((p, "templates"))
    items.append((REBUILD / "assets" / "zinc_nonbond_fallback.frcmod", "templates"))

    # 2) 对照体系
    for name in ("peptide_LxxLL.pdb", "DMSO.mol2"):
        items.append((ROOT / ".workbuddy" / name, "templates"))

    # 3) 配体（原始输入，仅供审计对照）
    lig_dir = ROOT / "docking_results_full" / "4GIZ"
    if lig_dir.is_dir():
        for ext in (".pdb", ".pdbqt", ".smi"):
            for p in sorted(lig_dir.glob(f"CHEMBL*{ext}")):
                items.append((p, "ligands"))
    else:
        print(f"[warn] 配体目录不存在：{lig_dir}")

    # 3b) Vina MODEL 1 姿态（13 号脚本生成，MD 起始结构必须用这些，不能用上面未对接的 pdb）
    poses_dir = REBUILD / "data" / "poses"
    if poses_dir.is_dir():
        for p in sorted(poses_dir.glob("*_pose1.pdb")):
            items.append((p, "ligands_vina"))
    else:
        print(f"[warn] 姿态目录不存在：{poses_dir}（先跑 13_prepare_vina_poses.py）")

    # 4) 管线脚本
    items.append((ROOT / ".workbuddy" / "candidates_pipeline.py", ""))

    # 5) 试点姿态
    items.append((REBUILD / "data" / "pilot" / "CHEMBL46917_pose1.pdb", "pilot"))

    return items


def main() -> int:
    ap = argparse.ArgumentParser(description="重建 B 线：上传本地资产到新 AutoDL 实例")
    ap.add_argument("--host", required=True)
    ap.add_argument("--port", type=int, default=22)
    ap.add_argument("--user", default="root")
    ap.add_argument("--password", default=os.environ.get("SSH_PASS", ""))
    ap.add_argument("--dry-run", action="store_true", help="只列出待上传清单，不连接")
    args = ap.parse_args()

    items = collect_files()
    missing = [(str(p), sub) for p, sub in items if not p.exists()]
    present = [(p, sub) for p, sub in items if p.exists()]

    print(f"[清单] 待上传 {len(present)} 个文件，缺失 {len(missing)} 个")
    for name, sub in missing:
        print(f"  [missing] {name}  (->{sub or REMOTE_BASE})")
    if missing:
        print("  降级指南：缺失文件需先在本地补齐，否则远端 tleap 会因缺模板/参数而失败。")

    if args.dry_run:
        for p, sub in present:
            print(f"  {p.stat().st_size:>10d}  {md5(p)}  {p}  -> {sub or '.'}")
        return 0

    if not args.password:
        sys.exit("未提供密码：用 --password 或环境变量 SSH_PASS")

    t0 = time.time()
    uploaded, verified, failed = 0, 0, []

    try:
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect(args.host, port=args.port, username=args.user,
                       password=args.password, timeout=30,
                       look_for_keys=False, allow_agent=False)
    except Exception as exc:  # noqa: BLE001
        print(f"[error] SSH 连接失败：{exc}")
        print("  降级指南：确认实例已开机、端口为最新值（AutoDL 控制台 → 容器实例 → SSH 登录）")
        return 1

    try:
        sftp = client.open_sftp()

        def mkdirs(path: str) -> None:
            parts, cur = path.strip("/").split("/"), ""
            for part in parts:
                cur = f"{cur}/{part}"
                try:
                    sftp.stat(cur)
                except IOError:
                    sftp.mkdir(cur)

        mkdirs(REMOTE_BASE)
        for sub in ("templates", "ligands", "ligands_vina", "pilot"):
            mkdirs(os.path.join(REMOTE_BASE, sub))

        for local, sub in present:
            remote_dir = os.path.join(REMOTE_BASE, sub) if sub else REMOTE_BASE
            remote = os.path.join(remote_dir, local.name)
            try:
                sftp.put(str(local), remote)
                uploaded += 1
            except Exception as exc:  # noqa: BLE001
                failed.append((str(local), f"upload: {exc}"))
                continue
            # 远端 md5 校验
            try:
                _, out, _ = client.exec_command(f"md5sum {remote} | awk '{{print $1}}'", timeout=30)
                rmd5 = out.read().decode().strip()
                if rmd5 == md5(local):
                    verified += 1
                else:
                    failed.append((str(local), f"md5 mismatch local={md5(local)} remote={rmd5}"))
            except Exception as exc:  # noqa: BLE001
                failed.append((str(local), f"verify: {exc}"))

        sftp.close()

        # 远端落盘确认
        _, out, err = client.exec_command(
            f"cd {REMOTE_BASE} && echo '--- templates ---' && ls -1 templates && "
            f"echo '--- ligands ---' && ls -1 ligands | wc -l && "
            f"echo '--- pipeline ---' && ls -l candidates_pipeline.py", timeout=30)
        print(out.read().decode(errors="replace").strip())
    finally:
        client.close()

    print(f"\n[结果] 上传 {uploaded}，md5 校验通过 {verified}，失败 {len(failed)}")
    for name, why in failed:
        print(f"  [fail] {name}: {why}")

    # 运行日志
    LOG_CSV.parent.mkdir(parents=True, exist_ok=True)
    new_file = not LOG_CSV.exists()
    with open(LOG_CSV, "a", newline="") as fh:
        w = csv.writer(fh)
        if new_file:
            w.writerow(["timestamp", "script", "host", "uploaded", "verified", "failed", "seconds"])
        w.writerow([time.strftime("%Y-%m-%d %H:%M:%S"), "11_seed_remote.py",
                    args.host, uploaded, verified, len(failed), round(time.time() - t0, 1)])
    print(f"[日志] 已追加 {LOG_CSV}")

    return 0 if not failed else 2


if __name__ == "__main__":
    raise SystemExit(main())
