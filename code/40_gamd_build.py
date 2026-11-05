#!/usr/bin/env python3
"""构建 GaMD 版 GROMACS（math-diff/gromacs-gamd, GROMACS 2025.4 + GaMD v2.x）。

目的：稿件 §3.9 需要在「修好的受体」上对 CHEMBL848（主）/ CHEMBL1129（对照）
跑 GaMD，产出接触占据率 + boost-校正 PMF（Figure 5）。

远端实例：MatPool hz-4.matpool.com:26255（4×RTX A4000）。
原版 /opt/gmx2025.4 是**不含 GaMD** 的 GROMACS 2025.4，必须另建一份。

已验证的基线（本机 /root/build_gmx.sh + /root/gromacs-2025.4/build/CMakeCache.txt）：
    cmake .. -DCMAKE_INSTALL_PREFIX=... -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA \
             -DGMX_MPI=OFF -DGMX_OPENMP=ON -DREGRESSIONTEST_DOWNLOAD=OFF \
             -DCMAKE_BUILD_TYPE=Release
本脚本沿用同一套参数，只换源码与安装前缀，避免环境差异导致构建失败。

坑（本机实测）：
  1. GitHub `git clone` 会报 "RPC failed; curl 16 Error in the HTTP2 framing layer"；
     改用 codeload 直链 + `--http1.1` 的 tarball。
  2. `archive/refs/heads/master.tar.gz` 会 404 —— 该仓库默认分支是 **main**。
  3. paramiko：后台任务用 `setsid nohup ... < /dev/null > log 2>&1 &`
     且**绝不读**该 channel 的 stdout（读了会 PipeTimeout）。

用法：
    python3 code/40_gamd_build.py download    # 后台下载源码包
    python3 code/40_gamd_build.py status      # 查看各阶段状态
    python3 code/40_gamd_build.py extract     # 解压
    python3 code/40_gamd_build.py configure   # cmake 配置（后台）
    python3 code/40_gamd_build.py build       # make -j（后台）
    python3 code/40_gamd_build.py install     # make install
    python3 code/40_gamd_build.py verify      # 校验 gmx 是否为 GaMD 版
    python3 code/40_gamd_build.py all         # 顺序跑完（每步后台 + 轮询）
"""
import sys
import time

import paramiko

HOST, PORT, USER, PW = "hz-4.matpool.com", 26255, "root", "wmKG=3YRsZaHmij0"

SRC_TAR = "/root/gamd_src_main.tar.gz"
SRC_DIR = "/root/gromacs-gamd"
BUILD_DIR = "/root/gamd-build"
PREFIX = "/usr/local/gromacs-gamd"

URL = "https://codeload.github.com/math-diff/gromacs-gamd/tar.gz/refs/heads/main"

# 注意：nvcc 不在 paramiko exec_command 的默认 PATH 里（.bashrc 不被加载），
# 必须显式给出，否则 cmake 报 "No CMAKE_CUDA_COMPILER could be found"。
# 这两个值直接取自本机已成功构建的原版 CMakeCache.txt。
CMAKE_CMD = f"""cmake -S . -B {BUILD_DIR} \\
  -DCMAKE_INSTALL_PREFIX={PREFIX} \\
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc \\
  -DCMAKE_CUDA_ARCHITECTURES=86 \\
  -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA \\
  -DGMX_MPI=OFF -DGMX_OPENMP=ON -DREGRESSIONTEST_DOWNLOAD=OFF \\
  -DCMAKE_BUILD_TYPE=Release"""


def cli() -> paramiko.SSHClient:
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, port=PORT, username=USER, password=PW, timeout=30)
    return c


def sh(c, cmd, timeout=120):
    """前台执行并读回输出。"""
    _, o, e = c.exec_command(cmd, timeout=timeout)
    out = o.read().decode(errors="replace").strip()
    err = e.read().decode(errors="replace").strip()
    return out, err


def run_bg(c, cmd):
    """完全脱离终端的后台任务；不读该 channel。"""
    c.exec_command(f"setsid nohup {cmd} < /dev/null > /dev/null 2>&1 &", timeout=10)


def p(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- 阶段

def do_download(c):
    sh(c, "pkill -f 'curl.*gamd_src' ; rm -f /root/gamd_src.tar.gz /root/gamd_src_main.tar.gz")
    run_bg(c, f"bash -c 'curl -sSL --http1.1 --retry 5 --retry-delay 3 "
               f"-o {SRC_TAR} {URL} > /root/dl_main.log 2>&1'")
    p("下载已启动，目标：", SRC_TAR)
    for k in range(60):
        time.sleep(10)
        out, _ = sh(c, f"ls -la {SRC_TAR} 2>/dev/null | awk '{{print $5}}'; "
                       f"pgrep -f 'curl.*gamd_src_main' | wc -l", 30)
        lines = out.split()
        size = int(lines[0]) if lines and lines[0].isdigit() else 0
        alive = int(lines[1]) if len(lines) > 1 and lines[1].isdigit() else 0
        p(f"  [{(k + 1) * 10}s] size={size/1048576:.1f} MiB  curl_alive={alive}")
        if size > 0 and alive == 0:
            out, _ = sh(c, f"tar tzf {SRC_TAR} > /dev/null 2>&1 && echo TAR_OK || echo TAR_BAD", 180)
            p("  完整性：" + out)
            return 0 if out == "TAR_OK" else 1
    p("  超时未完成，稍后用 status 再查")
    return 1


def do_extract(c):
    """解压前**必须**先校验 tarball 完整。

    坑（本机已踩）：下载是后台进行的，若在 curl 还在写时就解压，
    `tar xzf` 不报错、只给出一棵残缺的目录树，随后 cmake 会以
    "include could not find requested file src/testutils/TestMacros.cmake"
    这种**看起来像源码有问题**的方式失败，极难定位。
    因此这里先 `tar tzf`（校验）再解压，并显式检查几个关键子目录。
    """
    out, _ = sh(c, f"tar tzf {SRC_TAR} > /dev/null 2>&1 && echo TAR_OK || echo TAR_BAD", 600)
    p("tarball 校验：" + out)
    if out != "TAR_OK":
        return 1
    out, err = sh(c, f"rm -rf {SRC_DIR} && mkdir -p {SRC_DIR} && "
                     f"tar xzf {SRC_TAR} -C {SRC_DIR} --strip-components=1 && "
                     f"ls {SRC_DIR}/CMakeLists.txt {SRC_DIR}/src/testutils/TestMacros.cmake "
                     f"{SRC_DIR}/src/gromacs/applied_forces/gamd/gamdcontroller.cpp "
                     f"2>&1 && echo '---' && ls {SRC_DIR} | head -8", 600)
    p(out)
    return 0 if "No such file" not in out else 1


def do_configure(c):
    pre = (f"rm -rf {BUILD_DIR}; "
           "source /root/miniconda3/etc/profile.d/conda.sh; "
           "export PATH=/root/miniconda3/bin:/usr/local/cuda/bin:$PATH; "
           f"command -v cmake; cmake --version | head -1; cd {SRC_DIR}")
    run_bg(c, f"bash -c \"{pre} && {CMAKE_CMD} > /root/gamd_cmake.log 2>&1; "
              f"echo rc=$? >> /root/gamd_cmake.log\"")
    p("cmake 配置已启动")
    for k in range(60):
        time.sleep(15)
        out, _ = sh(c, "tail -3 /root/gamd_cmake.log 2>/dev/null; "
                       f"ls {BUILD_DIR}/Makefile 2>/dev/null | wc -l; "
                       "pgrep -f 'cmake' | wc -l", 60)
        p(f"  [{(k + 1) * 15}s] {out.splitlines()[-3:] if out else ''}")
        if "rc=" in out:
            tail, _ = sh(c, "tail -12 /root/gamd_cmake.log", 60)
            p(tail)
            return 0 if "rc=0" in out else 1
    return 1


def do_build(c):
    pre = ("source /root/miniconda3/etc/profile.d/conda.sh; "
           "export PATH=/root/miniconda3/bin:$PATH; "
           f"cd {BUILD_DIR}")
    run_bg(c, f"bash -c \"{pre} && make -j 28 > /root/gamd_make.log 2>&1; "
              f"echo rc=$? >> /root/gamd_make.log\"")
    p("make -j 28 已启动（GROMACS 全量编译，预计 20–40 min）")
    for k in range(160):
        time.sleep(30)
        out, _ = sh(c, "grep -c 'Building CXX' /root/gamd_make.log 2>/dev/null; "
                       "tail -1 /root/gamd_make.log 2>/dev/null | cut -c1-110; "
                       "pgrep -x cc1plus | wc -l", 60)
        p(f"  [{(k + 1) * 30}s] {out.strip()}")
        if "rc=" in out:
            tail, _ = sh(c, "tail -20 /root/gamd_make.log", 60)
            p(tail)
            return 0 if "rc=0" in out else 1
    return 1


def do_install(c):
    out, err = sh(c, f"cd {BUILD_DIR} && make install > /root/gamd_install.log 2>&1; "
                     f"echo rc=$?; tail -5 /root/gamd_install.log; "
                     f"ls {PREFIX}/bin/GMXRC 2>&1", 900)
    p(out)
    return 0 if "rc=0" in out else 1


def do_verify(c):
    out, _ = sh(c, f"source {PREFIX}/bin/GMXRC > /dev/null 2>&1; "
                   f"command -v gmx; gmx --version 2>&1 | head -25", 120)
    p(out)
    ok_bin = PREFIX in out or "/usr/local/gromacs-gamd" in out
    # GaMD 版会把 gamd 相关参数写进 mdp 帮助
    out2, _ = sh(c, f"source {PREFIX}/bin/GMXRC > /dev/null 2>&1; "
                    f"gmx grompp -h 2>&1 | grep -i gamd | head -5; "
                    f"echo '--- mdp options ---'; "
                    f"grep -rl 'gamd' {SRC_DIR}/src/ 2>/dev/null | head -5", 180)
    p(out2)
    p("VERDICT:", "GaMD 构建成功" if ok_bin else "未找到新二进制")
    return 0 if ok_bin else 1


def do_status(c):
    checks = [
        f"ls -la {SRC_TAR} 2>&1 | awk '{{print $5, $9}}'",
        f"ls -d {SRC_DIR} 2>&1; ls {SRC_DIR}/CMakeLists.txt 2>&1",
        f"ls {BUILD_DIR}/Makefile 2>&1",
        f"ls {PREFIX}/bin/GMXRC 2>&1",
        "tail -2 /root/gamd_cmake.log 2>/dev/null",
        "tail -2 /root/gamd_make.log 2>/dev/null",
        "pgrep -af 'cmake|cc1plus|make' | head -5",
        "nvidia-smi --query-gpu=index,utilization.gpu --format=csv,noheader",
    ]
    for cmd in checks:
        out, err = sh(c, cmd, 60)
        p("$", cmd)
        p(out if out else "(empty)")
        p("-" * 40)
    return 0


def do_all(c):
    for name, fn in (("download", do_download), ("extract", do_extract),
                     ("configure", do_configure), ("build", do_build),
                     ("install", do_install), ("verify", do_verify)):
        p(f"\n===== {name} =====")
        if fn(c) != 0:
            p(f"阶段 {name} 未通过，停止。")
            return 1
    return 0


def main() -> int:
    cmd = sys.argv[1] if len(sys.argv) > 1 else "status"
    fn = {"download": do_download, "extract": do_extract, "configure": do_configure,
          "build": do_build, "install": do_install, "verify": do_verify,
          "status": do_status, "all": do_all}.get(cmd)
    if fn is None:
        print(__doc__)
        return 2
    c = cli()
    try:
        return fn(c)
    finally:
        c.close()


if __name__ == "__main__":
    sys.exit(main())
