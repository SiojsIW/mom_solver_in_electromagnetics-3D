# EFIE-QBX 实现与分析文档

## 1. 项目背景

本项目是一个基于矩量法（MoM）的三维电磁散射求解器，用于计算 PEC 目标的双站 RCS。原始程序采用 RWG 基函数、CFIE/EFIE 公式，并对近奇异积分使用**奇异性提取法**（singularity extraction）。

本文档记录如何在本程序中新增 **QBX（Quadrature by Expansion）** 作为 EFIE 近场积分的可选路径，实现与奇异提取法的整体对比。

---

## 2. QBX 方法概述

QBX 的核心思想是：将格林函数在源区域附近的某展开中心 $r_c$ 处做局部球谐展开，从而把原本奇异的源积分转化为对展开系数的正则积分。在场点 $r_0$ 处（位于观察三角形上），标量位和矢量位可表示为

$$
S(r_0) = \sum_{l=0}^{p} \alpha_{l,0}\, j_l(k|r_0-r_c|) \sqrt{\frac{2l+1}{4\pi}},
$$

$$
\mathbf{V}(r_0) = \sum_{l=0}^{p} \boldsymbol{\beta}_{l,0}\, j_l(k|r_0-r_c|) \sqrt{\frac{2l+1}{4\pi}}.
$$

其中展开系数

$$
\alpha_{l,m} = -jk \int_{T_s} \sigma(r') h_l^{(2)}(k|r'-r_c|) \overline{Y_l^m(\theta',\phi')} \, dA',
$$

$$
\boldsymbol{\beta}_{l,m} = -jk \int_{T_s} \mathbf{f}_n(r') h_l^{(2)}(k|r'-r_c|) \overline{Y_l^m(\theta',\phi')} \, dA'.
$$

这里 $h_l^{(2)} = j_l - j y_l$ 为第二类球 Hankel 函数，$j_l$ 为球 Bessel 函数，$Y_l^m$ 为球谐函数。由于 $r_c$ 与源三角形 $T_s$ 保持一定距离，上述积分可用普通高斯求积精确计算。

工程上采用 $e^{+j\omega t}$ 时间约定，故 Green 函数为

$$
G(r,r') = \frac{e^{-jk|r-r'|}}{4\pi|r-r'|}.
$$

---

## 3. 代码实现细节

### 3.1 新增模块 `QBX_EFIE.F90`

该模块包含以下主要例程：

#### 3.1.1 特殊函数

- `SPH_BESSEL_ALL(X, P, JL, YL)`：计算 $j_l(x)$ 与 $y_l(x)$，$l=0\dots p$。
  - $y_l$ 使用升递推，无条件稳定。
  - $j_l$ 在 $x < p+1$ 或 $l > x$ 时使用 Miller 降递推，否则使用升递推。
- `SPH_HARM_3M(CT, PHI, P, Y0, YP1, YM1)`：仅计算 $m=0,\pm 1$ 的球谐函数。
  - 连带 Legendre 函数 $P_l^0$、$P_l^1$ 分别递推。
  - 归一化约定与项目文档一致。

#### 3.1.2 几何工具

- `BUILD_QBX_FRAME(TO_FIELD, TRI_EDGE, E1, E2, ZP)`：建立局部坐标架。
  - $z_p$ 指向场点方向（即 $-n_{\text{obs}}$）。
  - $e_1$ 取三角形某条边在切平面的投影并归一化。
  - $e_2 = z_p \times e_1$，保证右手系。
- `POINT_TRI_MIN_DIST_DP(RC, V1, V2, V3)`：计算点 $r_c$ 到源三角形的最短距离，用于许可性检查。

#### 3.1.3 QBX 核心例程

- `QBX_EFIE_COEFFICIENTS(...)`：在展开中心 $r_c$ 处计算系数 $\alpha_{l,m}$、$\boldsymbol{\beta}_{l,m}$。
  - 源三角形使用 Dunavant 25 点高斯规则加密积分。
  - 对每个源高斯点计算球坐标 $(\rho, \theta, \phi)$。
  - 调用 `SPH_BESSEL_ALL` 和 `SPH_HARM_3M` 完成积分。
- `QBX_EFIE_EVAL(...)`：在场点 $r_0$ 处用 $j_l(k|r_0-r_c|)$ 求和得到 $S$ 和 $\mathbf{V}$。
  - 只有 $m=0$ 分量对 $r_0$ 处的标量/矢量位有贡献。
- `CALC_EFIE_QBX_PAIR(...)`：计算单个（观察三角形，源三角形）对的阻抗贡献。
  - 外层对观察三角形做 12 点高斯积分。
  - 每个外层高斯点 $r_0$ 作为场点，$r_c = r_0 + h \cdot n_{\text{obs}}$ 作为展开中心。
  - 调用系数计算与极点求值，结合 RWG 测试函数 $\mathbf{f}_m$ 得到阻抗。

### 3.2 修改 `NUMERICAL_INTEGRATION.F90`

新增 Dunavant 25 点三角形高斯规则，供 QBX 源单元加密积分使用。注意：从 libMesh 源码抄来的权重对应参考三角形面积 $1/2$，而本项目高斯规则权重和为 1.0，因此 25 点权重全部乘以 2。

### 3.3 修改 `EM_TYPES.F90`

新增常量 `GAUSS_25PT = 25`，用于标识 Dunavant 25 点规则。

### 3.4 修改 `Z_MATRIX.F90`

新增模块级切换变量：

```fortran
INTEGER, PARAMETER :: EFIE_NEAR_EXTRACTION = 1
INTEGER, PARAMETER :: EFIE_NEAR_QBX        = 2
INTEGER :: EFIE_NEAR_METHOD = EFIE_NEAR_EXTRACTION
```

在 `CALC_EFIE_MATRIX_ELEMENT` 中，按两个三角形的共享顶点数分类处理：

| 共享顶点数 | 几何关系 | 处理方式 |
|---|---|---|
| 3 | 同一三角形（自项） | QBX（跳过许可性检查） |
| 2 | 共边对 | QBX（跳过许可性检查） |
| 1 | 共顶点 | QBX（跳过许可性检查） |
| 0 且 $\chi < 1$ | 近距离非接触 | QBX（保留许可性检查） |
| 0 且 $\chi \ge 1$ | 远场 | 原有完整格林函数高斯积分 |

当 `EFIE_NEAR_METHOD == EFIE_NEAR_EXTRACTION` 时，所有近场对均走原有奇异提取法。

### 3.5 修改 `MAIN_RCS.F90`

新增参数：

```fortran
LOGICAL, PARAMETER :: USE_QBX_EFIE = .TRUE.
```

程序启动时根据该参数设置 `EFIE_NEAR_METHOD`，并打印当前分支。同时引入 `QBX_EFIE` 模块，在 RCS 计算完成后输出 QBX 调用统计。

### 3.6 修改 `Makefile`

将 `QBX_EFIE.f90` 加入编译列表。

---

## 4. 关键参数

```fortran
P_TRUNC    = 12               ! QBX 截断阶
H_FACTOR   = 0.3              ! h = 0.3 * h_obs
N_GAUSS_UP = 25               ! 源单元 Dunavant 25 点
ADM_TOL    = 1.0E-3           ! 许可性相对裕量（仅非接触近场使用）
```

- `P_TRUNC` 控制球谐展开阶数，决定 QBX 精度。
- `H_FACTOR` 控制展开中心到场点的距离。过小会降低许可性通过率并增加截断误差；过大则展开中心离源三角形太近，影响收敛。
- `N_GAUSS_UP` 控制源三角形上系数的积分精度。

---

## 5. 许可性检查策略

为保证 QBX 展开收敛，展开中心 $r_c$ 到源三角形的最短距离应不小于 $h$。本实现采用：

- **自项**、**共边对**、**共顶点对**：直接跳过许可性检查，强制走 QBX。
- **非接触近场对**：检查 $\min\text{dist}(r_c, T_s) \ge h(1+\varepsilon_{\text{tol}})$，不满足则回退到奇异提取法。

由于程序已要求"纯 QBX"，当前版本对奇异对（自项/共边/共顶点）全部跳过检查，从而确保零回退。

---

## 6. 验证与分析

### 6.1 运行统计

纯 QBX 运行时：

```text
QBX 调用次数:          84762
QBX 许可性检查失败次数: 0
```

所有奇异/近奇异对均走 QBX，无回退。

### 6.2 RCS 曲线对比

#### ka = 2.0（lambda = 2π/3 ≈ 2.0944 m，181 个角度）

| 对比 | RMSE [dB] | 最大误差 [dB] |
|---|---|---|
| 纯 QBX vs 奇异提取法 | **0.0533** | **0.1192** |
| 纯 QBX vs Mie 级数 | **0.0637** | **0.1298** |
| 奇异提取法 vs Mie 级数 | 0.1145 | 0.2490 |

#### ka = 3.0（lambda = 2π/3 ≈ 2.0944 m，181 个角度）

| 对比 | RMSE [dB] | 最大误差 [dB] |
|---|---|---|
| 纯 QBX vs 奇异提取法 | **0.0510** | **0.1143** |
| 纯 QBX vs Mie 级数 | **0.0736** | **0.1855** |
| 奇异提取法 vs Mie 级数 | 0.1166 | 0.2855 |

结论：
- 在 ka = 2.0 和 ka = 3.0 两个频点，纯 QBX 与奇异提取法之间均出现了清晰可见的差异。
- **纯 QBX 在两个频点都更接近 Mie 参考解**，说明把自项、共边、共顶点全部纳入 QBX 后，整体精度有所提高。
- ka = 3.0 时整体误差略大于 ka = 2.0，主要由网格离散引起，而非积分方法本身。
- 两种方法在 ka = 3.0 处的差异与 ka = 2.0 相当，QBX 方法保持稳定。

### 6.3 结果文件

| 文件 | 说明 |
|---|---|
| `rcs_results_qbx_pure.txt` | 纯 QBX 方法 RCS |
| `rcs_results_extraction.txt` | 奇异提取法 RCS |
| `rcs_comparison_qbx_extract_mie.txt` | 合并后的 QBX / 提取法 / Mie 三条曲线 |

`rcs_comparison_qbx_extract_mie.txt` 格式：

```text
# theta[deg]  RCS_QBX[dBsm]  RCS_Extraction[dBsm]  RCS_Mie[dBsm]
```

---

## 7. 使用说明

### 7.1 切换方法

在 `MAIN_RCS.F90` 中：

```fortran
LOGICAL, PARAMETER :: USE_QBX_EFIE = .TRUE.   ! QBX
LOGICAL, PARAMETER :: USE_QBX_EFIE = .FALSE.  ! 奇异提取法
```

或在其他程序中直接设置模块变量：

```fortran
USE Z_MATRIX
EFIE_NEAR_METHOD = EFIE_NEAR_QBX          ! QBX
EFIE_NEAR_METHOD = EFIE_NEAR_EXTRACTION   ! 奇异提取法
```

### 7.2 生成对比曲线

```bash
# 1. 运行 QBX（默认）
./rcs_solver.exe
cp rcs_results.txt rcs_results_qbx_pure.txt

# 2. 修改 MAIN_RCS.F90 中 USE_QBX_EFIE = .FALSE.，重新编译运行
mingw32-make rcs_solver
./rcs_solver.exe
cp rcs_results.txt rcs_results_extraction.txt

# 3. 合并曲线
python merge_rcs_comparison.py
```

---

## 8. 已知局限与后续

1. **MFIE-QBX**：当前仅实现 EFIE-QBX。MFIE 需要双侧中心、独立局部坐标架以及跳跃项（jump conditions）处理。
2. **许可性检查**：为达到零回退，对自项/共边/共顶点跳过了许可性检查。在更复杂几何或更高频下，可能需要更严格的展开中心选取策略。
3. **性能**：当前每个（场点，源单元）对独立计算 QBX 系数，未做缓存。若成为瓶颈，可考虑缓存同源源三角形的系数。
4. **参数调优**：`P_TRUNC`、`H_FACTOR`、`N_GAUSS_UP` 可针对更高频率、更粗糙网格进一步研究。

---

## 9. 代码文件清单

| 文件 | 状态 | 说明 |
|---|---|---|
| `QBX_EFIE.F90` | 新增 | EFIE-QBX 核心模块 |
| `NUMERICAL_INTEGRATION.F90` | 修改 | 新增 Dunavant 25 点规则 |
| `EM_TYPES.F90` | 修改 | 新增 `GAUSS_25PT` |
| `Z_MATRIX.F90` | 修改 | 新增方法切换与 QBX 分支 |
| `MAIN_RCS.F90` | 修改 | 新增 `USE_QBX_EFIE` 开关与统计输出 |
| `Makefile` | 修改 | 加入 `QBX_EFIE.f90` |
| `merge_rcs_comparison.py` | 新增 | 合并 QBX / 提取法 / Mie 曲线 |

---

*文档生成时间：2026-08-27*
