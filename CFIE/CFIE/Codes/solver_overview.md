# CFIE/RWG 矩量法求解器概览

## 1. 程序定位

本代码是一个基于 **RWG（Rao–Wilton–Glisson）边元** 的 **Method of Moments（MoM）** 电磁散射求解器，使用 **Fortran 90** 编写。它支持：

- 理想导体（PEC）三角面片目标的散射分析
- 电场积分方程（EFIE）、磁场积分方程（MFIE）及其线性组合 **CFIE**
- 平面波入射激励
- 双站 RCS（雷达散射截面）远场计算

## 2. 文件清单与职责

| 文件 | 类型 | 主要功能 |
|------|------|----------|
| `MAIN.f90` | 主程序 | 流程编排：输入 → 读网格 → 预处理 → 填充矩阵 → 激励 → 求解 → 散射场 |
| `Global_module.f90` | 模块 | 公共数组、物理常数、Dunavant 积分点、远场参数 |
| `Input_setting.f90` | 子程序 | 读取频率、入射波方向、RCS 采样设置、迭代求解器参数 |
| `Load_triangle_tetrahedron_mesh.F90` | 子程序 | 读取 Nastran（`.nas`）格式网格 |
| `preprocessing_triangle.f90` | 子程序 | 计算三角形面积、重心，构造三角形-RWG 关联表，判断是否为闭合并启用 CFIE |
| `form_RWG_Triangle.f90` | 子程序 | 基于堆排序识别公共边，生成 RWG 边索引 `side_index` |
| `MoM_SS.f90` | 子程序 | 填充阻抗矩阵（EFIE/MFIE/CFIE），主双重循环遍历三角形对 |
| `MoM_MVP.f90` | 子程序 | 矩阵-向量乘积，供 BiCGstab 调用 |
| `Comp_eltsI.f90` | 子程序 | 计算标量/矢量势积分 `I1~I4`，含奇异提取与解析奇异积分 `Spara_cal` |
| `Comp_eltsU.f90` | 子程序 | 计算 Green 函数梯度相关积分 `U1~U5`，含 `gradG_S` 与解析奇异积分 `Spara_cal_2` |
| `excitation_filling.f90` | 子程序 | 根据入射平面波填充右端项（激励向量） |
| `Bicgstab.f90` | 子程序 | BiCGstab 迭代法求解线性方程组 |
| `scattering_field.f90` | 子程序 | 远场 `Eθ`、`Eφ` 计算并输出 RCS（dB） |

## 3. 控制方程

### 3.1 组合场积分方程（CFIE）

当目标为闭合 PEC 体时（`3×Num_triangle == 2×Num_side`），自动启用 CFIE：

```
Z_CFIE = α · Z_EFIE + η0 · (1 − α) · Z_MFIE
```

其中 `α = 0.5`（默认值），`η0 = 120π Ω` 为自由空间波阻抗。

### 3.2 阻抗矩阵元素

对每个场三角形 `m`、源三角形 `n` 以及它们上的 RWG 基函数对 `(j, i)`，电部分为：

```
Amn_E = lj · li · [ coef1 · ( (rmj·rni) I1 − rni·I2 ) + ctmp ]
ctmp  = coef1 · (I4 − rmj·I3) − coef2 · I1
coef1 = jωμ0 / (16π)
coef2 = j / (4πωε0)
```

磁部分在非自三角形时：

```
Amn_H = coef3 · lj · li · [ ( (rmj·rni) n − (rni·n) rmj ) · U1
                          + (rni·n) U3 − rni·U5 + ctmp1 ]
coef3 = 1 / (16π)
ctmp1 = − (rmj × n) · U2 − n · U4
```

当 `m == n` 时，MFIE 主值采用 Duffy/面积重心近似解析式：

```
Amn_H = (lj·li) / (8·Am) · [ A1 + rmj·rni − centroid·(rmj+rni) ]
A1    = (1/6) Σ ri·rj   （三角形顶点坐标二次项和）
```

## 4. 程序流程

```
main
 ├── Input_setting(1)              # 读取用户输入
 ├── Load_triangle_tetrahedron_mesh # 读取 .nas 网格
 ├── preprocessing_triangle        # 生成 RWG、面积、重心、判断 CFIE
 ├── Input_setting(2)              # 计算角频率、波长、波数、近场阈值
 ├── MoM_SS                        # 填充阻抗矩阵
 ├── excitation_filling            # 右端项
 ├── allocate Inn
 ├── bicgstab                      # 迭代求解电流系数
 └── scattering_field              # 计算并输出 RCS
```

## 5. 关键数据结构

### 5.1 公共数组（`Common_array`）

- `point_cor(3, Num_point)`：节点坐标
- `triangle_point(3, Num_triangle)`：三角形节点编号
- `side_index(4, Num_side)`：RWG 边信息 `[node1, node2, tri_plus, tri_minus]`
- `triangle_area(Num_triangle)`、`triangle_centroid(3, Num_triangle)`
- `triangle_rwgs(Num_triangle)`：每个三角形关联的 RWG 列表
- `SparseValue(Num_side, Num_side)`：稠密阻抗矩阵（物理意义上为稠密）
- `Inn(Num_side)`：未知电流系数
- `excitation(Num_side)`：右端激励向量

### 5.2 公共常数（`Common_constants`）

```fortran
PI, MU0, EPS0, ETA0
iu = (0,1)
CFSIE = .TRUE. / .FALSE.   ! 是否启用 CFIE
alpha = 0.5                ! CFIE 混合系数
threshold_R = 0.1 * lambda ! 远/近场阈值
ratio_limit1 = 1e-6 * lambda  ! 判定 w0≈0 的阈值
ratio_limit2 = 3e-2           ! 防止 atan(有限/0) 的阈值
```

### 5.3 积分规则（`Dunavant_integration_formulism`）

采用 4 点 Dunavant 面积分规则：

```fortran
SNum_samp = 4
weightS = (/ -0.5625, 1.5625/3.0, 1.5625/3.0, 1.5625/3.0 /)
(epsilS, etaS, xiS)  # 面积坐标
```

## 6. 数值方法要点

### 6.1 RWG 基函数

RWG 函数定义在共享一条公共边的两个相邻三角形上，表达式为：

```
f_n(r) = ± (ln / 2An) · (r − r_n^±)
```

其中 `ln` 为公共边长度，`An` 为三角形面积，`r_n^±` 为对顶点。

### 6.2 矩阵填充

`MoM_SS.f90` 使用双重循环遍历所有三角形对。对每个三角形对：

1. 计算两三角形重心距离 `dist`
2. 根据 `dist >= threshold_R` 选择远场/近场分支
3. 调用 `Comp_eltsI` 计算势积分
4. 若启用 CFIE 且 `m != n`，调用 `Comp_eltsU` 计算梯度积分
5. 组装 EFIE/MFIE 矩阵元并累加到 `SparseValue(j,i)`

### 6.3 迭代求解

使用 `BiCGstab` 求解：

```
Z · Inn = excitation
```

矩阵-向量乘积由 `MVP_MoM.f90` 实现，对每个未知量做直接向量点积（本代码未做快速算法加速）。

### 6.4 远场计算

`scattering_field.f90` 对每个远场方向 `(θ, φ)` 计算辐射积分：

```
E_far(θ,φ) ∝ −jωμ0 / (4π) · Σ  In · ln · (ρc · ê) · exp(j k · rc)
```

其中 `ρc = rc − rn`，`ê` 为 `θ̂` 或 `φ̂`。

## 7. 编译构建

项目使用 CMake，主程序文件需包含 `PROGRAM` 关键字。`CMakeLists.txt` 自动：

1. 扫描所有 `.f90`
2. 将不含 `PROGRAM` 的文件编译为静态库 `common_modules`
3. 每个含 `PROGRAM` 的文件生成独立可执行目标

```bash
cmake -B cmake-build-debug
cmake --build cmake-build-debug
```

## 8. 输入/输出

### 8.1 输入

- `Meshes/<name>.nas`：Nastran 格式三角面片网格
- `Setting.txt`（可选）：批量输入参数

### 8.2 输出

- `Output/FixedPhi<i>_RCS_db.txt`
- `Output/FixedTheta<i>_RCS_db.txt`

每行格式：`角度(度)  Eθ_RCS_dB  Eφ_RCS_dB  总RCS_dB`

## 9. 奇异处理说明（详见 `singularity_treatment.md`）

奇异/近奇异积分是 MoM 核心难点，本代码主要策略为：

1. **距离阈值分离**：`threshold_R = 0.1λ` 区分远场（纯数值积分）与近场/奇异（解析提取）。
2. **Green 函数奇异性提取**：将 `exp(−ikR)/R` 拆分为 `1/R + F1(R)`。
3. **解析面积分**：利用三角形局部坐标和边参数，解析计算 `1/R` 及 `R` 的各阶矩积分。
4. **梯度主值提取**：对 `∇G` 做 `k³I30 − (I31 + k²/2 · I32)/An` 分解，其中 `I31`、`I32` 由 `Spara_cal_2` 解析给出。
5. **数值鲁棒性处理**：对 `w0≈0`、`R0≈0`、`Rm+lm≈0` 等极限情况做分支保护。

具体公式、源码摘录与调用关系见 **`singularity_treatment.md`**。
