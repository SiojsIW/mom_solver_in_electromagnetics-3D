# SCATTER_3D_FRAME 项目交接文档

## 1. 项目是什么

这是一个**三维电磁散射矩量法（MoM）求解器**，用于计算理想导体（PEC）目标在平面波照射下的双站雷达散射截面（RCS）。

- **方法**：RWG（Rao-Wilton-Glisson）基函数 + 三角形表面网格
- **核心方程**：EFIE（电场积分方程），可选 CFIE（组合场积分方程，含 MFIE 部分）
- **验证手段**：金属球 Mie 级数解析解
- **语言**：Fortran 90（求解器）+ Python（后处理 / Mie 对比）

## 2. 关键文件说明

| 文件 | 作用 |
|---|---|
| `MAIN_RCS.F90` | 主程序：读网格 → 建几何 → 组装 EFIE 矩阵 → 求解电流 → 计算双站 RCS |
| `EM_TYPES.F90` | 自定义数据类型（节点、三角形、边、RWG 基函数等） |
| `MESH_GEOMETRY.F90` | 网格读取、法向校正、边提取、几何量计算 |
| `RWG_BASIS_BUILD.F90` | 从公共边构建 RWG 基函数 |
| `Z_MATRIX.F90` | 组装 EFIE/MFIE 阻抗矩阵，含奇异性提取；近场方法由 `EFIE_NEAR_METHOD` 选择（提取/QBX/DIRECTFN） |
| `DIRECTFN_BRIDGE.F90` | 第三种近场方法（NEAR_METHOD=3）：经 ISO_C_BINDING 调 C++ 库 DIRECTFN 做 4D 全数值奇异积分 |
| `DIRECTFN_WRAPPER.cpp` | DIRECTFN 的 extern "C" 包装（`dfn_tri_iss`），链 `DIRECTFN-master/build/libdirectfn.a` |
| `GREEN_FUNCTIONS.F90` | 格林函数及其梯度、平滑余项 |
| `SINGULAR_INTEFRAL.F90` | 1/R 标量/线性势解析公式、MFIE 静态主值矩 |
| `RHS.F90` | 右端向量（入射平面波与 RWG 测试） |
| `FAR_FIELD.F90` | 远区散射场 / 双站 RCS 计算 |
| `NUMERICAL_INTEGRATION.F90` | 三角形高斯积分规则（1/3/4/7/12 点） |
| `mie_compare.py` | **Mie 级数参考解 + MoM 结果对比** |
| `test_singular.F90` | 奇异积分单元测试 |
| `Makefile` | 构建脚本（使用 `mingw32-make`） |
| `.claude/WORKLOG.md` | 详细工作日志 |

## 3. 如何构建和运行

```bash
# 构建
mingw32-make clean && mingw32-make

# 运行金属球 RCS（默认读取根目录 point.txt + triangle_point.txt）
./rcs_solver.exe

# 或使用子目录中的网格文件
./rcs_solver.exe "网格文件/point.txt" "网格文件/triangel_point.txt"

# 与 Mie 级数对比
python mie_compare.py
```

输出文件：
- `rcs_results.txt`：MoM 计算的双站 RCS
- `mie_mom_comparison.txt`：MoM vs Mie 逐角度对比
- `mie_mom_summary.txt`：每次对比的摘要（追加）

## 4. 关键配置参数

在 `MAIN_RCS.F90` 顶部：

```fortran
REAL, PARAMETER :: LAMBDA = 4.188790   ! 波长 (m)，ka = 2π/LAMBDA * a
REAL, PARAMETER :: ETA0 = 120.0 * PI   ! 自由空间波阻抗
REAL, PARAMETER :: E0 = 1.0            ! 入射电场幅度
REAL, PARAMETER :: ALPHA_CFIE = 1.0    ! 1.0 = 纯 EFIE，0.5 = CFIE
INTEGER, PARAMETER :: N_ANGLES = 181   ! 0°~180°，步长 1°
```

入射条件（固定）：
- 入射方向 `k_i = +z`
- 极化方向 `E_i = +x`
- 观察面 `φ = 0`（xz 平面）

球半径 `a = 1.0 m`，因此 `ka = 2π / LAMBDA`。

## 5. Mie 级数对比：极化分量选择（重点！）

这是本项目踩过的最大坑，**必须注意**：

### 几何与极化对应关系

- 入射 `E_i = +x`，沿 `+z` 传播
- 观察面 `φ = 0` 即 xz 平面
- 入射电场 **位于散射面（scattering plane）内**
- 因此 MoM 算的是 **parallel 分量**，对应 NIST 约定中的 **S₂**

### `mie_compare.py` 中的实现

```python
# Mie 系数（PEC 球），正号使光学定理 sigma_ext = sigma_sca 成立
a_n = jn / hn
b_n = (jn + x * jn_d) / (hn + x * hn_d)

# 对当前 MoM 几何，取 parallel 分量（S_2）
S_parallel = np.sum(coeff * (a_n * pi_n + b_n * tau_n))

# RCS
rcs_theta = (wavelength ** 2 / np.pi) * np.abs(S_parallel) ** 2
```

**注意**：早期版本错误地使用了 `a_n*tau_n + b_n*pi_n`，这实际上算的是 perpendicular 分量（S₁），导致中间角度出现 2–6 dB 的“伪误差”。

### 验证结果

修正极化后，纯 EFIE 与 Mie 级数对比：

| ka | RMSE | 最大误差 |
|---|---|---|
| 1.5 | 0.097 dB | 0.170 dB |
| 3.0 | 0.117 dB | 0.285 dB |

## 6. 近场奇异积分：三种可选方法

`MAIN_RCS.F90` 顶部 `NEAR_METHOD` 选择近场（奇异/近奇异）阻抗积分的处理方法：

```fortran
INTEGER, PARAMETER :: NEAR_METHOD = 3   ! 1=奇异提取法  2=QBX  3=DIRECTFN
```

| 方法 | 原理 | ka=3 时 RCS vs Mie |
|---|---|---|
| 1 奇异提取法 | 把 Green 函数解析奇异性提取后用闭式 1/R 势公式 | RMSE 0.117 dB，最大 0.286 dB |
| 2 QBX | 展开求积法（级数展开 + 边界项） | RMSE 0.074 dB，最大 0.186 dB |
| 3 DIRECTFN | 4 重 Gauss-Legendre 全数值积分（第三方 C++ 库） | RMSE 0.115 dB，最大 0.282 dB |

三种方法精度同量级，详细逐角度数据见 `rcs_comparison_three_methods_mie.txt`
（及单方法原始结果 `rcs_results_{extraction,qbx_pure,directfn}.txt`）。

### DIRECTFN 集成（方法 3）

- 库：`DIRECTFN-master/`（C++，已编译出 `build/libdirectfn.a`）
- 接口：`DIRECTFN_WRAPPER.cpp` 提供 `extern "C" dfn_tri_iss(adj, pts, k0, n_gauss, iss_rwg, iss_const)`
  - `adj`：3=ST 共三角形（3 点）/ 2=EA 共边（4 点）/ 1=VA 共顶点（5 点）
  - contour 顺序（桥接 `DIRECTFN_BRIDGE.F90` 已按此构造）：
    VA = `(公共点, 观察其余两点, 源其余两点)`，EA = `(公共边e1, 公共边e2, 观察第三顶点, 源第三顶点)`
- 调用开关：`Z_MATRIX.F90` 中 `EFIE_NEAR_METHOD == EFIE_NEAR_DIRECTFN` 时走 `CALC_EFIE_DIRECTFN_PAIR`
- 高斯阶数：`DIRECTFN_N_GAUSS_ST / _EA / _VA`（按邻接类型分别设置，默认 8/8/6）。
  单对代价 ∝ N⁴ 且 ST:EA:VA ≈ 24:9:2。ka=3 标定结果：N=12、全 8、(8,8,6) 三档的
  RCS vs Mie 完全相同（RMSE 0.115 dB），故默认取 8/8/6
- 性能（ka=3 全量求解墙钟）：原始逐 RWG 对调用 45 min → tri-pair 缓存+对象复用后
  N=12 时 6.4 min → N=(8,8,6) 时约 2.1 min（总加速 ~21×）。
  关键优化（均已验证结果逐位/逐角度不变）：
  1. **tri-pair 缓存**（~9×）：同一 (观察,源) 三角形对只积分一次，结果按哈希表
     （键=(观察,源) 三角形号）缓存，3×3 个 RWG 对共享（`DIRECTFN_RESET_CACHE` /
     `DIRECTFN_PRINT_STATS`，MAIN_RCS 在填充前/后调用）。83754 次 RWG 对调用
     去重为 9306 次唯一三角形对（716 ST + 2148 EA + 6442 VA）
  2. **wrapper 对象复用**：6 个算法对象（ST/EA/VA × RWG/Constant）静态持有，
     每对象独立跟踪 (k0,N)，变化时才重设高斯节点（`DIRECTFN_WRAPPER.cpp`）
- 精度仲裁（k=3 逐元）：DFN 的 Z 误差平均 2.1e-07 / 最大 2.3e-06，优于提取法
  （3.3e-04）与 QBX（7.5e-03）

**已踩过的坑（重要）**：
1. **Iss 9 分量索引转置**：DIRECTFN 的 `Iss` 线性布局是 **(源槽位, 观察槽位) 行优先**
   （见 `directfn_kernel_tri.cpp` 的 `pmf_rwg_val_` 顺序），Fortran 侧取
   (观察=IO, 源=JO) 分量必须读 `ISS_RWG(JO + 3*(IO-1))`。读错的故障指纹：
   对角项全对（~1e-11）、非对角项全是转置值（差 30%~600%）、ST 因对称不受影响。
2. **K_D 作用域**：缓存命中路径不经过 wrapper 调用，凡命中路径也要用的标量
   （如 K_D）必须在缓存查询前赋值，否则读到未初始化局部变量（NaN 电流）。
3. 调试：`DIRECTFN_DEBUG > 0` 打印 `@@DFNRAW` 原始积分行供外部参考程序解析。

## 7. 程序当前状态

- **EFIE 实现已验证正确**：在 ka=1.5 和 ka=3.0 均与 Mie 级数高度吻合。
- **MFIE 实现存在历史问题**：此前尝试 CFIE 时发现 MFIE 有系统性偏差，目前 `ALPHA_CFIE = 1.0`（纯 EFIE）。
- **网格**：360 节点 / 716 三角形 / 1074 RWG，节点在单位球面上，质量良好。
- **远场积分**：已启用 1→4 三角形细分，每个子三角 12 点高斯。

## 8. 常见坑点

1. **网格文件名**：三角形文件是 `triangel_point.txt`（注意拼写，不是 `triangle`）。
2. **Mie 极化**：见第 5 节，是最容易复现的误判来源。
3. **DIRECTFN 分量索引**：`Iss` 9 分量内存布局为 (源槽, 观察槽)，Fortran 取 (观察 IO, 源 JO) 须读
   `ISS_RWG(JO + 3*(IO-1))`，见第 6 节。读错的指纹是对角全对、非对角全错。
4. **内谐振**：纯 EFIE 在闭合 PEC 体的某些 ka 处会出现内谐振（对金属球，TM 约在 ka≈2.744，TE 约在 ka≈4.493），但 ka=1.5 和 ka=3.0 的验证误差极小，说明当前关心的频点不受此影响。
5. **构建命令**：Windows 下用 `mingw32-make`，不是 `make`。

## 9. 快速验证流程

```bash
# 1. 确认 LAMBDA 对应想要验证的 ka
# 2. 构建并运行
mingw32-make clean && mingw32-make
./rcs_solver.exe "网格文件/point.txt" "网格文件/triangel_point.txt"
# 3. 对比
python mie_compare.py
```

如果 `MoM_theta` 与 `Mie_theta` 的 RMSE 在 0.1–0.3 dB 量级，说明 EFIE 结果可信。
