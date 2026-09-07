#!/usr/bin/env python3
"""
扫描不同 ka 值，运行 EFIE RCS 求解器并与 Mie 级数对比。
球半径固定 a=1.0 m，通过改变波长 lambda = 2*pi/ka 来改变 ka。
"""
import subprocess
import re
import os
import sys
import numpy as np

# 扫描的 ka 值列表（避开已知的球腔体谐振点 ka≈2.744, 4.493）
KA_VALUES = [
    0.5, 0.7, 1.0, 1.2, 1.5, 1.8, 2.0, 2.2, 2.4,
    2.6, 2.8, 3.0, 3.1416, 3.3, 3.5, 3.8, 4.0, 4.2,
    4.5, 4.8, 5.0, 5.5, 6.0
]

A = 1.0
MAIN_RCS = "MAIN_RCS.F90"
SUMMARY_FILE = "mie_mom_summary.txt"

def clear_summary():
    if os.path.exists(SUMMARY_FILE):
        os.remove(SUMMARY_FILE)

def set_lambda(lambda_val):
    """修改 MAIN_RCS.F90 中的 LAMBDA 参数值"""
    with open(MAIN_RCS, 'r', encoding='utf-8') as f:
        content = f.read()
    content_new = re.sub(
        r"REAL, PARAMETER :: LAMBDA = [0-9.]+",
        f"REAL, PARAMETER :: LAMBDA = {lambda_val:.6f}",
        content
    )
    with open(MAIN_RCS, 'w', encoding='utf-8') as f:
        f.write(content_new)

def run_case(ka):
    lambda_val = 2.0 * np.pi * A / ka
    print(f"\n=== ka = {ka:.4f}, lambda = {lambda_val:.4f} m ===")
    set_lambda(lambda_val)

    # 确保可执行文件未被锁定
    for exe in ['rcs_solver.exe']:
        if os.path.exists(exe):
            try:
                os.remove(exe)
            except OSError:
                pass

    # 编译
    result = subprocess.run(
        ["mingw32-make", "rcs_solver"],
        capture_output=True, text=True, encoding='utf-8', errors='replace'
    )
    if result.returncode != 0:
        print(f"  编译失败: {result.stderr[-500:]}")
        return None

    # 运行 RCS
    result = subprocess.run(
        ["./rcs_solver.exe", "网格文件/point.txt", "网格文件/triangel_point.txt"],
        capture_output=True, text=True, encoding='utf-8', errors='replace'
    )
    if result.returncode != 0:
        print(f"  RCS 运行失败: {result.stderr[-500:]}")
        return None

    # 运行对比（使用当前 Python 解释器，避免环境不一致）
    result = subprocess.run(
        [sys.executable, "mie_compare.py"],
        capture_output=True, text=True, encoding='utf-8', errors='replace'
    )
    if result.returncode != 0:
        print(f"  对比脚本失败: {result.stderr[-500:]}")
        return None

    # 从摘要文件读取最后一行
    with open(SUMMARY_FILE, 'r', encoding='utf-8') as f:
        lines = f.read().strip().split('\n')
    last = lines[-1]
    parts = last.split()
    # parts: k wavelength rmse max_err rmse2 max_err2
    return {
        'ka': float(parts[0]) * A,
        'lambda': float(parts[1]),
        'rmse': float(parts[2]),
        'max_err': float(parts[3]),
        'rmse2': float(parts[4]),
        'max_err2': float(parts[5])
    }

def main():
    clear_summary()
    results = []
    for ka in KA_VALUES:
        res = run_case(ka)
        if res is not None:
            results.append(res)
            print(f"  RMSE(|cosθ|>0.1) = {res['rmse2']:.3f} dB, max_err = {res['max_err2']:.3f} dB")

    print("\n" + "=" * 70)
    print("EFIE 扫频结果摘要")
    print("=" * 70)
    print(f"{'ka':>8} {'lambda':>10} {'RMSE_all':>10} {'max_all':>10} {'RMSE_cos':>10} {'max_cos':>10}")
    for r in results:
        print(f"{r['ka']:8.4f} {r['lambda']:10.4f} {r['rmse']:10.3f} {r['max_err']:10.3f} {r['rmse2']:10.3f} {r['max_err2']:10.3f}")
    print("=" * 70)

    # 保存摘要表格
    with open('ka_scan_summary.txt', 'w', encoding='utf-8') as f:
        f.write("# ka lambda RMSE_all[dB] max_all[dB] RMSE(|cos|>0.1)[dB] max(|cos|>0.1)[dB]\n")
        for r in results:
            f.write(f"{r['ka']:.6f} {r['lambda']:.6f} {r['rmse']:.6f} {r['max_err']:.6f} {r['rmse2']:.6f} {r['max_err2']:.6f}\n")
    print("摘要已保存到 ka_scan_summary.txt")

if __name__ == '__main__':
    main()
