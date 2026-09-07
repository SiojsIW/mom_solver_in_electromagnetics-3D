#!/usr/bin/env python3
"""
金属球双站 RCS 的 Mie 级数参考解与 MoM 结果对比

入射条件：+z 方向，x 极化，PEC 球半径 a=1.0 m，波长 lambda=2.0 m (ka=pi)
观察面：phi=0 (xz 平面)，theta 从 0° 到 180°
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from scipy.special import spherical_jn, spherical_yn, assoc_legendre_p


def _angular_functions(n, theta_rad):
    """
    计算 Mie 角函数 pi_n(theta) 和 tau_n(theta)。
    使用递推公式，在 theta=0, pi 处也稳定。
    """
    n = np.asarray(n, dtype=int)
    theta = float(theta_rad)
    mu = np.cos(theta)
    n_max = int(n.max())

    pi = np.zeros(n_max + 1)
    tau = np.zeros(n_max + 1)

    # 递推初值
    pi[0] = 0.0
    pi[1] = 1.0

    # 递推 pi_n = (2n-1)/(n-1) * mu * pi_{n-1} - n/(n-1) * pi_{n-2}
    for nn in range(2, n_max + 1):
        pi[nn] = ((2.0 * nn - 1.0) * mu * pi[nn - 1] - nn * pi[nn - 2]) / (nn - 1.0)

    # tau_n = n * mu * pi_n - (n+1) * pi_{n-1}
    tau[0] = 0.0
    for nn in range(1, n_max + 1):
        tau[nn] = nn * mu * pi[nn] - (nn + 1.0) * pi[nn - 1]

    return pi[n], tau[n]


def mie_rcs_sphere(a, wavelength, theta_deg, n_max=100):
    """
    计算 PEC 金属球在 phi=0 平面内的双站 RCS（theta 与 phi 极化分量）。

    参数
    ----
    a : float
        球半径 [m]
    wavelength : float
        波长 [m]
    theta_deg : array_like
        散射角 theta [度]
    n_max : int
        Mie 级数截断阶数

    返回
    ----
    rcs_theta : ndarray
        theta 极化 RCS [m^2]
    rcs_phi : ndarray
        phi 极化 RCS [m^2]
    """
    k = 2.0 * np.pi / wavelength
    x = k * a
    theta_deg = np.asarray(theta_deg, dtype=float)
    n_theta = theta_deg.size

    # 球形 Bessel / Hankel (第二类)
    n = np.arange(1, n_max + 1)
    jn = spherical_jn(n, x)
    yn = spherical_yn(n, x)
    hn = jn - 1j * yn

    # 导数
    jn_d = spherical_jn(n, x, derivative=True)
    yn_d = spherical_yn(n, x, derivative=True)
    hn_d = jn_d - 1j * yn_d

    # Mie 系数 (PEC 球)
    # 注意：此处用正号以满足光学定理 (sigma_ext = sigma_sca)。
    # RCS 只看 |S|^2，符号不影响 dBsm 对比，但影响复散射场相位。
    a_n = jn / hn
    b_n = (jn + x * jn_d) / (hn + x * hn_d)

    coeff = (2.0 * n + 1.0) / (n * (n + 1.0))

    rcs_theta = np.zeros(n_theta)
    rcs_phi = np.zeros(n_theta)

    for i in range(n_theta):
        theta = np.deg2rad(theta_deg[i])
        pi_n, tau_n = _angular_functions(n, theta)

        # 几何定义（与 MoM 一致）：
        #   入射 k = +z, E = +x, 观察面 phi = 0 (xz 平面)
        #   入射电场位于 scattering plane 内，因此应取 parallel 分量，
        #   即 NIST 约定中的 S_2。
        # 在本代码的角函数约定下，parallel 分量为：
        #   S_parallel = sum coeff * (a_n * pi_n + b_n * tau_n)
        S_parallel = np.sum(coeff * (a_n * pi_n + b_n * tau_n))

        # RCS = (lambda^2 / pi) * |S_parallel|^2
        rcs_theta[i] = (wavelength ** 2 / np.pi) * np.abs(S_parallel) ** 2

    return rcs_theta, rcs_phi


def load_mom_rcs(filename):
    """读取 rcs_results.txt，返回 theta[deg], rcs_theta[dBsm], rcs_phi[dBsm], k"""
    # 文件头含中文，显式指定 UTF-8
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
    return data[:, 0], data[:, 1], data[:, 2], k


def dbsm(x):
    """m^2 -> dBsm"""
    return 10.0 * np.log10(np.maximum(x, 1e-20))


def main():
    a = 1.0

    theta_mom, rcs_theta_mom_db, rcs_phi_mom_db, k = load_mom_rcs('rcs_results.txt')
    if k is None:
        k = 2.0 * np.pi / 2.0  # 默认 ka=pi
    wavelength = 2.0 * np.pi * a / k

    # Mie 参考解
    rcs_theta_mie, rcs_phi_mie = mie_rcs_sphere(a, wavelength, theta_mom, n_max=120)
    rcs_theta_mie_db = dbsm(rcs_theta_mie)
    rcs_phi_mie_db = dbsm(rcs_phi_mie)

    # 只比较 theta 极化（phi 极化理论为 0，MoM 中应很小）
    valid = rcs_theta_mom_db > -150.0
    diff = rcs_theta_mom_db[valid] - rcs_theta_mie_db[valid]
    rmse = np.sqrt(np.mean(diff ** 2))
    max_err = np.max(np.abs(diff))
    mean_err = np.mean(diff)

    # 排除 theta=90° 附近（cosθ 接近 0，数值敏感区）
    theta_rad = np.deg2rad(theta_mom)
    valid2 = valid & (np.abs(np.cos(theta_rad)) > 0.1)
    diff2 = rcs_theta_mom_db[valid2] - rcs_theta_mie_db[valid2]
    rmse2 = np.sqrt(np.mean(diff2 ** 2))
    max_err2 = np.max(np.abs(diff2))

    print("=" * 60)
    print("MoM vs Mie 级数 RCS 对比")
    print("=" * 60)
    print(f"球半径 a = {a} m, 波长 lambda = {wavelength:.4f} m, ka = {k*a:.4f}")
    print(f"角度点数: {len(theta_mom)}")
    print(f"theta 极化 RMSE (全角度):   {rmse:.3f} dB")
    print(f"theta 极化最大误差 (全角度): {max_err:.3f} dB")
    print(f"theta 极化平均偏差 (全角度): {mean_err:.3f} dB")
    print(f"theta 极化 RMSE (|cosθ|>0.1):   {rmse2:.3f} dB")
    print(f"theta 极化最大误差 (|cosθ|>0.1): {max_err2:.3f} dB")
    print("-" * 60)
    print(f"{'theta':>8} {'MoM_theta':>14} {'Mie_theta':>14} {'err':>10}")
    idx = np.where(valid)[0]
    for i in idx[::10]:
        print(f"{theta_mom[i]:8.1f} {rcs_theta_mom_db[i]:14.3f} {rcs_theta_mie_db[i]:14.3f} {rcs_theta_mom_db[i]-rcs_theta_mie_db[i]:10.3f}")
    print("=" * 60)

    # 保存对比数据
    out = np.column_stack([
        theta_mom,
        rcs_theta_mom_db, rcs_theta_mie_db,
        rcs_phi_mom_db, rcs_phi_mie_db
    ])
    header = "theta[deg]  MoM_theta[dBsm]  Mie_theta[dBsm]  MoM_phi[dBsm]  Mie_phi[dBsm]"
    np.savetxt('mie_mom_comparison.txt', out, fmt='%.4f', header=header)
    print("详细对比已保存到 mie_mom_comparison.txt")

    # 保存摘要（供扫描脚本收集）
    summary = f"{k:.6f} {wavelength:.6f} {rmse:.6f} {max_err:.6f} {rmse2:.6f} {max_err2:.6f}\n"
    with open('mie_mom_summary.txt', 'a', encoding='utf-8') as f:
        f.write(summary)

    # 尝试绘图
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10, 5))
        plt.plot(theta_mom, rcs_theta_mom_db, 'b-', label='MoM theta')
        plt.plot(theta_mom, rcs_theta_mie_db, 'r--', label='Mie theta')
        plt.plot(theta_mom, rcs_phi_mom_db, 'g-', alpha=0.5, label='MoM phi')
        plt.xlabel('Theta [deg]')
        plt.ylabel('RCS [dBsm]')
        plt.title(f'MoM vs Mie, ka={2*np.pi*a/wavelength:.2f}, RMSE={rmse:.2f} dB')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig('mie_mom_comparison.png', dpi=150)
        print("对比图已保存到 mie_mom_comparison.png")
    except ImportError:
        print("未安装 matplotlib，跳过绘图")


if __name__ == '__main__':
    main()
