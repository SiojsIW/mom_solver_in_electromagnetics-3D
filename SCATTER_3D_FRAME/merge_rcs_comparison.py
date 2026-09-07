#!/usr/bin/env python3
"""
合并 QBX、奇异提取法、Mie 级数三条 RCS 曲线到单一文件，方便绘图对比。

入射条件：+z 方向，x 极化，PEC 球半径 a=1.0 m
波长 lambda = pi m，波数 k = 2.0 rad/m，ka = 2.0
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from scipy.special import spherical_jn, spherical_yn


def _angular_functions(n, theta_rad):
    """计算 Mie 角函数 pi_n(theta) 和 tau_n(theta)，在 theta=0, pi 处稳定。"""
    n = np.asarray(n, dtype=int)
    theta = float(theta_rad)
    mu = np.cos(theta)
    n_max = int(n.max())

    pi = np.zeros(n_max + 1)
    tau = np.zeros(n_max + 1)

    pi[0] = 0.0
    pi[1] = 1.0
    for nn in range(2, n_max + 1):
        pi[nn] = ((2.0 * nn - 1.0) * mu * pi[nn - 1] - nn * pi[nn - 2]) / (nn - 1.0)

    for nn in range(1, n_max + 1):
        tau[nn] = nn * mu * pi[nn] - (nn + 1.0) * pi[nn - 1]

    return pi[n], tau[n]


def mie_rcs_sphere(a, wavelength, theta_deg, n_max=120):
    """计算 PEC 金属球在 phi=0 平面内的 theta 极化双站 RCS [m^2]。"""
    k = 2.0 * np.pi / wavelength
    x = k * a
    theta_deg = np.asarray(theta_deg, dtype=float)
    n_theta = theta_deg.size

    n = np.arange(1, n_max + 1)
    jn = spherical_jn(n, x)
    yn = spherical_yn(n, x)
    hn = jn - 1j * yn

    jn_d = spherical_jn(n, x, derivative=True)
    yn_d = spherical_yn(n, x, derivative=True)
    hn_d = jn_d - 1j * yn_d

    a_n = jn / hn
    b_n = (jn + x * jn_d) / (hn + x * hn_d)
    coeff = (2.0 * n + 1.0) / (n * (n + 1.0))

    rcs_theta = np.zeros(n_theta)

    for i in range(n_theta):
        theta = np.deg2rad(theta_deg[i])
        pi_n, tau_n = _angular_functions(n, theta)
        S_parallel = np.sum(coeff * (a_n * pi_n + b_n * tau_n))
        rcs_theta[i] = (wavelength ** 2 / np.pi) * np.abs(S_parallel) ** 2

    return rcs_theta


def load_rcs(filename):
    """读取 rcs_results.txt 格式，返回 theta[deg], rcs_theta[dBsm], rcs_phi[dBsm], rcs_total[dBsm], k。"""
    k = None
    with open(filename, 'r', encoding='utf-8') as f:
        lines = []
        for ln in f:
            s = ln.strip()
            if s.startswith('#'):
                if '波数 k' in s or 'k (rad/m)' in s:
                    parts = s.split(':')
                    if len(parts) >= 2:
                        try:
                            k = float(parts[1].strip().split()[0])
                        except ValueError:
                            pass
            elif s:
                lines.append(s)
    data = np.loadtxt(lines)
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3], k


def dbsm(x):
    return 10.0 * np.log10(np.maximum(x, 1e-20))


def main():
    a = 1.0

    theta_qbx, rcs_theta_qbx, _, _, k = load_rcs('rcs_results_qbx_pure.txt')
    theta_ext, rcs_theta_ext, _, _, _ = load_rcs('rcs_results_extraction.txt')

    if k is None:
        k = 2.0  # 默认 fallback
    wavelength = 2.0 * np.pi * a / k

    if not np.allclose(theta_qbx, theta_ext, atol=1e-6):
        raise ValueError("QBX 与提取法角度采样不一致")

    theta = theta_qbx

    rcs_theta_mie = mie_rcs_sphere(a, wavelength, theta, n_max=120)
    rcs_theta_mie_db = dbsm(rcs_theta_mie)

    out = np.column_stack([
        theta,
        rcs_theta_qbx,
        rcs_theta_ext,
        rcs_theta_mie_db
    ])

    header = (
        "# 金属球双站 RCS 对比\n"
        f"# 条件：PEC 球 a=1.0 m, lambda={wavelength:.6f} m, k={k:.4f} rad/m, ka={k*a:.4f}\n"
        "# 入射：+z 方向，x 极化，观察面 phi=0 (xz-plane)\n"
        "# theta[deg]  RCS_QBX[dBsm]  RCS_Extraction[dBsm]  RCS_Mie[dBsm]\n"
    )

    with open('rcs_comparison_qbx_extract_mie.txt', 'w', encoding='utf-8') as f:
        f.write(header)
        np.savetxt(f, out, fmt='%10.4f %16.6f %16.6f %16.6f')

    # 计算 QBX vs 提取法 误差
    diff_qbx_ext = rcs_theta_qbx - rcs_theta_ext
    rmse_qbx_ext = np.sqrt(np.mean(diff_qbx_ext ** 2))
    max_err_qbx_ext = np.max(np.abs(diff_qbx_ext))

    # 计算 QBX vs Mie 误差
    diff_qbx_mie = rcs_theta_qbx - rcs_theta_mie_db
    rmse_qbx_mie = np.sqrt(np.mean(diff_qbx_mie ** 2))
    max_err_qbx_mie = np.max(np.abs(diff_qbx_mie))

    # 计算 提取法 vs Mie 误差
    diff_ext_mie = rcs_theta_ext - rcs_theta_mie_db
    rmse_ext_mie = np.sqrt(np.mean(diff_ext_mie ** 2))
    max_err_ext_mie = np.max(np.abs(diff_ext_mie))

    print("=" * 60)
    print("三条 RCS 曲线已合并到 rcs_comparison_qbx_extract_mie.txt")
    print("=" * 60)
    print(f"角度点数: {len(theta)}")
    print(f"QBX vs 奇异提取法  RMSE = {rmse_qbx_ext:.4f} dB, 最大误差 = {max_err_qbx_ext:.4f} dB")
    print(f"QBX vs Mie          RMSE = {rmse_qbx_mie:.4f} dB, 最大误差 = {max_err_qbx_mie:.4f} dB")
    print(f"提取法 vs Mie       RMSE = {rmse_ext_mie:.4f} dB, 最大误差 = {max_err_ext_mie:.4f} dB")
    print("=" * 60)


if __name__ == '__main__':
    main()
