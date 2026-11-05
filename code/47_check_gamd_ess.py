#!/usr/bin/env python3
"""47_check_gamd_ess.py — GaMD 重加权**有效性验收**（写进 Methods 的那条判据）。

背景（本项目踩过两次的坑）：
  σ₀ 越小 boost 越弱，重加权越可信，但增强也越弱。只报 ⟨ΔV⟩ 或 β·σ_ΔV
  是不够的——**决定成败的是重加权的有效样本量**：

      w_i ∝ exp(β·ΔV_i)，  ESS = (Σw)² / Σw²  （Kish 有效样本量）
      若 ΔV 近似正态且 s = β·σ_ΔV，则  ESS/N ≈ exp(−s²)

  实测对照：
      σ₀=6.0  → βs≈11  → ESS/N 天文数字级地趋 0   → PMF 完全不可信
      σ₀=0.5  → βs≈2.3 → ESS/N = 0.06–0.21%       → 90% 权重落在 3–15 帧，不可用
      σ₀=0.25 → 目标 βs≤1.2 → ESS/N ≥ 20%

验收阈值（本项目的口径，写进 Methods）：
  ESS/N ≥ 15% 且 max(β·σ_ΔV) ≤ 1.5  → PASS
  否则请降低 σ₀ 后重跑

用法：
  python3 code/47_check_gamd_ess.py            # 查远端 gamd_v2
  python3 code/47_check_gamd_ess.py <远端目录名> [最少样本数]
"""
from __future__ import annotations

import os
import sys

import numpy as np
import paramiko

# 云端实例停机/重启后 SSH 端口会轮换：用环境变量覆盖，别再改源码
HOST = os.environ.get("HPV_SSH_HOST", "hz-4.matpool.com")
PORT = int(os.environ.get("HPV_SSH_PORT", "26255"))
USER = os.environ.get("HPV_SSH_USER", "root")
PW = os.environ.get("HPV_SSH_PW", "wmKG=3YRsZaHmij0")
REMOTE_BASE = "/root/autodl-tmp/HPV_structure/candidates_run"
KT = 0.001987204259 * 300.0
BETA = 1.0 / KT

PASS_ESS_FRAC = 0.15      # ESS/N 下限
PASS_BETA_SD = 1.5        # β·σ_ΔV 上限


def analyse_reweight(path: str):
    """返回 (n, mean, sd, beta_sd, max, ess, ess_frac, k90)。"""
    v = []
    try:
        with open(path) as fh:
            for ln in fh:
                if ln.startswith("#") or not ln.strip():
                    continue
                p = ln.split()
                if len(p) < 9:
                    continue
                try:
                    v.append(float(p[7]) + float(p[8]))
                except ValueError:
                    continue
    except FileNotFoundError:
        return None
    if not v:
        return None
    v = np.asarray(v)
    w = np.exp(BETA * (v - v.max()))
    w /= w.sum()
    ess = 1.0 / (w ** 2).sum()
    o = np.sort(w)[::-1]
    k90 = int(np.searchsorted(np.cumsum(o), 0.9)) + 1
    return (len(v), float(v.mean()), float(v.std(ddof=1)),
            float(BETA * v.std(ddof=1)), float(v.max()),
            float(ess), float(ess / len(v)), k90)


def main() -> int:
    tag = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("GAMD_TAG", "gamd_v2")
    min_n = int(sys.argv[2]) if len(sys.argv) > 2 else 0

    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, port=PORT, username=USER, password=PW, timeout=30)
    sftp = c.open_sftp()
    base = f"{REMOTE_BASE}/{tag}"
    try:
        jobs = sorted(d for d in sftp.listdir(base) if "_rep" in d)
    except IOError:
        print(f"远端目录不存在：{base}")
        return 1

    rows, verdicts = [], []
    for j in jobs:
        rp = f"{base}/{j}/gamd-reweight.dat"
        tmp = f"/tmp/_ess_{j}.dat"
        try:
            sftp.get(rp, tmp)
        except IOError:
            rows.append((j, "prep", None))
            continue
        r = analyse_reweight(tmp)
        os.remove(tmp)
        if r is None or r[0] < min_n:
            rows.append((j, f"n={0 if r is None else r[0]}（样本不足）", None))
            continue
        n, mean, sd, bsd, mx, ess, frac, k90 = r
        rows.append((j, None, (n, mean, sd, bsd, mx, ess, frac, k90)))
        ok = frac >= PASS_ESS_FRAC and bsd <= PASS_BETA_SD
        verdicts.append((j, ok))

    print(f"# GaMD 重加权验收 —— {tag}   (KT={KT:.4f} kcal/mol, beta={BETA:.4f})")
    print(f"# 阈值：ESS/N >= {PASS_ESS_FRAC:.0%} 且 beta*sigma <= {PASS_BETA_SD}")
    print(f"{'job':<24}{'n':>7}{'<dV>':>8}{'sd':>7}{'beta*sd':>9}"
          f"{'max':>8}{'ESS':>7}{'ESS/N':>9}{'前90%wt帧数':>12}  判定")
    for j, note, r in rows:
        if r is None:
            print(f"{j:<24}{note}")
            continue
        n, mean, sd, bsd, mx, ess, frac, k90 = r
        ok = frac >= PASS_ESS_FRAC and bsd <= PASS_BETA_SD
        print(f"{j:<24}{n:>7}{mean:>8.2f}{sd:>7.2f}{bsd:>9.2f}{mx:>8.2f}"
              f"{ess:>7.0f}{frac:>8.1%}{k90:>12}  {'PASS' if ok else 'FAIL'}")
    if verdicts:
        bad = [j for j, ok in verdicts if not ok]
        print("-" * 90)
        if bad:
            print(f"❌ {len(bad)}/{len(verdicts)} 不通过：{', '.join(bad)}")
            print("   → 降低 sigma0 后重跑（经验上 sigma0 减半 ≈ beta*sigma 减半）。")
            return 2
        print(f"✅ 全部 {len(verdicts)} 个副本通过，可安全重加权出 PMF。")
    sftp.close()
    c.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
