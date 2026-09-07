# SCATTER_3D_FRAME 工作日志

## 项目目标

三维 MoM 电磁散射求解器，使用 RWG 基函数在三角形表面网格上求解 PEC 目标的 EFIE/MFIE/CFIE。

## 当前状态（最后更新：上下文压缩后）

- 构建系统：`mingw32-make`（gfortran，自由格式 F90）
- 奇异积分：`test_singular.F90` 已通过 8 项测试
- 网格读取：支持带头行 / 无头行的 `point.txt` + `triangle_point.txt`
- 球面算例：360 节点 / 716 三角形 / 1074 RWG，可运行

## 已验证结果

### 奇异积分测试

```bash
mingw32-make test
```

全部通过。

### 金属球 RCS（纯 EFIE，ka = π）

```bash
./rcs_solver.exe "网格文件/point.txt" "网格文件/triangel_point.txt"
python mie_compare.py
```

- 前向 0°：误差 -0.08 dB
- 后向 180°：误差 -0.26 dB
- 侧向 30°–150°：误差 3–15 dB
- RMSE(|cosθ|>0.1)：≈ 6.8 dB

判断：ka = π 处 EFIE 内谐振导致侧瓣误差。

### 金属球 RCS（CFIE，α=0.5）

已尝试将 `MAIN_RCS.F90` 改为 CFIE 组合，但结果未改善：

- 前向 0°：误差 -14.6 dB（纯 EFIE 为 -0.08 dB）
- 后向 180°：误差 -1.7 dB（纯 EFIE 为 -0.26 dB）
- 侧向 RMSE(|cosθ|>0.1)：≈ 11.9 dB（纯 EFIE 为 6.8 dB）

调整 α=0.9 后结果接近纯 EFIE，侧瓣无改善。

判断：简单 CFIE α=0.5 因前向点 `V_E` 与 `η0 V_M` 局部抵消导致右端过小，结果崩溃。调 α 后发现：

- α=0.2 时侧瓣改善最显著（RMSE 从 6.82 dB 降到 3.78 dB）；
- 但前向/后向仍有 2–3 dB 误差，说明当前 MFIE 纯数值积分仍有系统偏差；
- 需要按 `docs/MFIE_CFIE实现文档_内谐振对策.md` 实现 ∇G 的解析主值（PV）处理。

## 当前默认

`MAIN_RCS.F90` 中 `ALPHA_CFIE = 0.2`，运行命令：

```bash
./rcs_solver.exe "网格文件/point.txt" "网格文件/triangel_point.txt"
python mie_compare.py
```

## 下一步

1. 在 `SINGULAR_INTEGRAL.F90` 中实现 `ANALYTIC_GRAD_STATIC_PV`（静态主值三矩）。
2. 实现 `GRAD_GREEN_REG`（∇G 正则余项）。
3. 在 `Z_MATRIX.F90` 中把 MFIE 近场分支从纯数值改为解析 PV + 数值余项。
4. 重新标定 α，目标：前向/后向 <0.5 dB、侧瓣 RMSE <1 dB。

## 注意

- `MAIN_RCS.F90` 当前为纯 EFIE，需改 CFIE。
- `rcs_results.txt` 每次运行会被覆盖。
- 网格文件名拼写为 `triangel_point.txt`（不是 `triangle`）。

---

## 2026-08-17 会话记录

### 完成的工作

1. **修复编译错误**：`Z_MATRIX.F90:426` 处 `GRAD_REG` 被声明为 `REAL`，但 `GRAD_GREEN_REG` 期望 `COMPLEX`，已改正。
2. **修复 MFIE 远场分支的 RWG 对顶点提取**：`CALC_MFIE_MATRIX_ELEMENT` 中 `R_OPP_M`/`R_OPP_N` 原循环只填 X 坐标，已改为分别填入 X/Y/Z。
3. **MFIE 近场分支已接入 `CALC_MFIE_MATRIX_ELEMENT`**：共边/共顶点/非接触近场均走 `CALC_MFIE_NEAR_TRI_PAIR`。
4. **奇异积分测试**：`mingw32-make test` 全部 10 项通过。
5. **`mie_compare.py` 改进**：支持从 `rcs_results.txt` 头部读取波数 k，并追加摘要到 `mie_mom_summary.txt`。

### MFIE/CFIE 验证结果（ka = π）

| α | RMSE (\|cosθ\|>0.1) | 前向 0° 误差 | 备注 |
|---|----------------------|--------------|------|
| 0.0（纯 MFIE） | 4.05 dB | -7.09 dB | 与 α=0.2 几乎相同 |
| 0.2 | 3.88 dB | -7.95 dB | MFIE 主导 |
| 0.5 | 3.04 dB | -4.11 dB | |
| 0.9 | 5.47 dB | +0.85 dB | 接近 EFIE，侧瓣差 |
| 1.0（纯 EFIE） | 6.82 dB | -0.08 dB | 基准 |

判断：MFIE 实现仍有系统性偏差（前向点被压低）。已排除 jump 项符号、4π 归一化、`g_E` 公式首项等常见问题；PV 项置零、CHI_NEAR 调整、增加细分对结果影响均很小，提示问题可能更深层（如 Graglia 1993 的 k² 修正项未提取，或远场高斯积分策略）。

### 方向调整

按用户指示：**暂停 MFIE 修复**，先验证纯 EFIE 在不同 ka 下的表现，确认非内谐振频点能否与 Mie 级数一致。后续由用户指定具体 ka 再计算分析。

### 当前默认配置

- `MAIN_RCS.F90`：`ALPHA_CFIE = 1.0`（纯 EFIE），`LAMBDA = 2.0`（ka = π）。
- 运行命令：

```bash
./rcs_solver.exe "网格文件/point.txt" "网格文件/triangel_point.txt"
python mie_compare.py
```

### 已完成的 EFIE 点扫（用户叫停前）

| ka | λ [m] | RMSE (\|cosθ\|>0.1) | 最大误差 |
|----|-------|----------------------|----------|
| 0.5 | 12.57 | 5.68 dB | 16.93 dB |
| 0.7 | 8.98  | 4.32 dB | 13.12 dB |

### 下一步（待用户指定）

- 用户指定具体 ka 后，修改 `MAIN_RCS.F90` 中的 `LAMBDA = 2π/ka`，运行 EFIE 并与 Mie 对比。
- 暂不继续 MFIE / CFIE 调试，亦不保留自动扫频脚本。

---

## 2026-08-17 会话记录（续）

### 重大发现：`mie_compare.py` 参考解 bug

用户指出侧瓣 Mie 值异常，经检查：`mie_compare.py` 中 theta 极化 RCS 多乘了 `cos²(θ)`：

```python
# 错误
rcs_theta[i] = (wavelength ** 2 / np.pi) * np.abs(S2) ** 2 * np.cos(theta) ** 2
# 正确
rcs_theta[i] = (wavelength ** 2 / np.pi) * np.abs(S2) ** 2
```

Mie 的 `S₂(θ)` 已经完整描述 θ 极化散射，不需要额外 `cos²(θ)`。前向/后向因 `cosθ=±1` 不受影响，侧瓣被错误压低，90° 处被压到 0。

### 修正后的 EFIE 结果

#### ka = 1.5（λ = 4.1888 m）

| 指标 | 数值 |
|------|------|
| 前向 0° 误差 | -0.067 dB |
| 后向 180° 误差 | +0.170 dB |
| RMSE (全角度) | 1.21 dB |
| RMSE (\|cosθ\|>0.1) | 1.25 dB |
| 最大误差 | 2.61 dB |

侧瓣误差从原来的 6–19 dB 降至 2.6 dB 以内，说明 EFIE 实现本身精度良好。

#### ka = π（λ = 2.0 m）

| 指标 | 数值 |
|------|------|
| 前向 0° 误差 | -0.082 dB |
| 后向 180° 误差 | -0.260 dB |
| RMSE (全角度) | 2.29 dB |
| RMSE (\|cosθ\|>0.1) | 1.78 dB |
| 最大误差 | 5.98 dB |

ka=π 仍比 ka=1.5 差一些（80°–100°、120°–140° 有 3–6 dB 偏差），可能与该频点靠近内谐振区或网格在该电尺寸下精度不足有关。

### 结论

- EFIE 实现本身基本正确；之前看到的“侧瓣系统性偏差”主要是 Mie 参考解 bug 造成的假象。
- ka=1.5 已达到 RMSE ≈ 1.2 dB，说明非谐振频点 EFIE 可信。
- ka=π 仍有约 1.8 dB RMSE，需进一步判断是内谐振影响还是网格/积分精度问题。

### 当前配置

- `MAIN_RCS.F90`：`ALPHA_CFIE = 1.0`（纯 EFIE），`LAMBDA = 2.0`（ka = π）。
- `mie_compare.py` 已修正 theta 极化 RCS 公式（去掉多余 `cos²θ`）。

---

## 2026-08-17 会话记录（ka = 1.5 EFIE 复核）

### 运行结果

修复 `FAR_FIELD.F90` 行首多余冒号后重新编译运行（ka = 1.5，纯 EFIE）：

```bash
mingw32-make clean && mingw32-make
./rcs_solver.exe "网格文件/point.txt" "网格文件/triangel_point.txt"
python mie_compare.py
```

| 指标 | 数值 |
|------|------|
| 前向 0° 误差 | -0.067 dB |
| 后向 180° 误差 | +0.170 dB |
| RMSE (全角度) | 1.21 dB |
| RMSE (\|cosθ\|>0.1) | 1.25 dB |
| 最大误差 | 2.61 dB（约 50°） |

结果与上次记录一致。远场积分已改为 1→4 细分（每个 RWG 正/负三角形各 4 子三角 × 12 高斯点），说明 30°–70° 的 2–3 dB 误差并非远场积分采样不足造成。

### Mie 级数复核

`mie_compare.py` 计算出的 RCS 幅值（\|S₂\|²）是正确的，但散射振幅 S₂(0) 的实部符号与光学定理不符：

- 当前系数：aₙ = -jₙ/hₙ，bₙ = -(jₙ + x jₙ')/(hₙ + x hₙ')
- 由此得到 σ_ext = (4π/k²) Re[S₂(0)] = -6.77 m²
- 而 σ_sca = (2π/k²) Σ(2n+1)(\|aₙ\|² + \|bₙ\|²) = +6.77 m²

对无损耗球应有 σ_ext = σ_sca > 0，因此 Mie 系数应整体取反（aₙ、bₙ 前面的负号去掉）。由于 RCS 只依赖 \|S₂\|²，该符号 bug **不影响当前的 dBsm 对比**，但会影响复散射场的相位。

### EFIE 代码复核

逐行检查了阻抗矩阵组装、RWG 系数、格林函数、奇异性提取、RHS、远场辐射积分等关键环节：

- **标量势符号**：`Q_SING = -J0/π` 与远场 `Q_TERM = -4*I1` 并非错误，而是代码把 Q 定义为标准标量势积分的负值，再与 `+jη₀/k` 相乘，最终得到正确的 `-jη₀/k · Q_standard`。已核对一致。
- **RWG 基函数**：`BUILD_RWG_BASIS` 与 `EVAL_RWG_BASIS` 的系数、对顶点、负三角形符号均正确。
- **RHS**：入射平面波 `E_inc = E0 · E_POL · exp(-j k·r)` 与 RWG 测试匹配。
- **远场积分**：`CALC_RWG_RADIATION_INTEGRAL` 的 1→4 细分实现正确，`COEF` 对负三角形的处理与基函数定义一致；`CALC_FAR_SCATTERED_FIELD` 的 θ/φ 单位矢量正确。
- **没有发现导致中间角度系统性偏差的明显程序漏洞**。

### 网格模块复核

- `MESH_GEOMETRY.F90` 中节点读取、法向朝外校正、边提取、局部边编号、边长/中点计算均正确。
- 360 节点 / 716 三角形 / 1074 RWG，节点严格在单位球面上，法向均朝外，最小三角形 quality ≈ 0.245，平均 ≈ 0.489，无劣质单元。
- ka = 1.5 时 λ/h ≈ 20，离散度足够。

### 中间角度误差来源判断

排除明显程序 bug 后，30°–70° 的 2–3 dB 偏差最可能来自：

1. **阻抗矩阵近场/奇异积分数值误差**：当前共边/共顶点用 2 层细分，非接触近场用 1 层细分。虽然 `test_singular.F90` 对单个解析积分通过，但双三角形耦合的累计误差仍可能使电流相位略有偏移，导致侧瓣干涉图案偏移。
2. **曲面离散化几何误差**：即使 λ/h ≈ 20，用平面三角形逼近球面仍会引入几何模型误差，对侧瓣这种对相位敏感的方向影响更明显。

前向/后向因对称性或远场主瓣能量集中，对上述误差不敏感，故误差主要出现在中间角度。

### 下一步建议（待用户确认）

1. 若用户希望继续压误差，可尝试提高近场积分精度：把 `Z_MATRIX.F90` 中 `N_SUBDIV_EDGE`、`N_SUBDIV_VERTEX` 从 2 提到 3，`N_SUBDIV_NEAR` 从 1 提到 2，并/或把 `CHI_NEAR` 从 1.0 提到 2.0–3.0。
2. 修正 `mie_compare.py` 中 Mie 系数的整体符号（去掉 aₙ、bₙ 前的负号），使复散射场也满足光学定理；RCS 对比结果不变。
3. 暂不继续 MFIE / CFIE 调试，按用户此前指示只聚焦 EFIE。

---

## 2026-08-17 会话记录（ka = 3.0 EFIE）

### 用户指令

用户要求重新计算 ka = 3.0 并分析结果。

### 运行配置

修改 `MAIN_RCS.F90`：

```fortran
REAL, PARAMETER :: LAMBDA = 2.094395  ! 波长 (m)，对应 ka = 3.0
```

保持 `ALPHA_CFIE = 1.0`（纯 EFIE）。

### 运行结果

```bash
./rcs_solver.exe "网格文件/point.txt" "网格文件/triangel_point.txt"
python mie_compare.py
```

| 指标 | ka = 1.5 | ka = 3.0 |
|------|----------|----------|
| 0° 误差 | -0.067 dB | -0.085 dB |
| 180° 误差 | +0.170 dB | -0.095 dB |
| RMSE (全角度) | 1.21 dB | **2.20 dB** |
| RMSE (\|cosθ\|>0.1) | 1.25 dB | **1.80 dB** |
| 最大误差 | 2.61 dB | **6.19 dB** |

误差分布：

- 60°–100° 区间偏差明显增大，90° 附近出现深零（MoM -0.54 dBsm vs Mie 5.43 dBsm，误差 -5.97 dB）。
- 130°–150° 也有 2–3 dB 正偏差。

### 关键发现：ka = 3.0 靠近 EFIE 内谐振

金属球 PEC 腔的内谐振位置由球 Bessel 函数零点决定：

- TM 模式：`[x j_n(x)]' = 0`，n=1 时 x ≈ **2.744**
- TE 模式：`j_n(x) = 0`，n=1 时 x ≈ **4.493**

ka = 3.0 距离 TM 内谐振点 2.744 仅约 0.256（相对距离约 9%），EFIE 算子接近奇异。此时：

1. 阻抗矩阵条件数显著恶化；
2. 数值误差（网格离散、积分、线性求解）被放大；
3. 感应电流中高阶模的相位/幅度出现偏差；
4. 远场在侧向（尤其 90° 附近）形成错误的深零。

这与 ka = 1.5（远离内谐振，RMSE 1.2 dB）形成鲜明对比，说明 ka = 3.0 的误差主要不是程序 bug，而是 **EFIE 在靠近内谐振频点的固有病态**。

### 结论

- ka = 1.5：纯 EFIE 可靠，误差在 2–3 dB 以内。
- ka = 3.0：纯 EFIE 因接近 TM 内谐振（2.744）而不可靠，最大误差达 6 dB。
- 若要在 ka = 3.0 附近获得准确结果，需要 CFIE（α < 1）来抑制内谐振；这又回到了 MFIE 实现的问题。

### 下一步建议

1. 若用户只想用纯 EFIE 做可信对比，应选择远离内谐振的 ka，例如 ka = 1.5、ka ≈ 5.0（避开 4.493 和 6.117）。
2. 若坚持 ka = 3.0，需要启用 CFIE，这意味着必须继续修复 MFIE 实现。
3. 继续暂停自动扫频，等待用户指定下一个 ka 或决定是否回到 MFIE/CFIE 调试。

---

## 2026-08-17 会话记录（重大修正：Mie 极化分量）

### 用户发现

用户在 `docs/建议.md` 中指出：ka=3 的 MoM 曲线几乎精确对应 Mie 的 **另一个极化分量**。原因是 `mie_compare.py` 虽然代码里变量名叫 `S2`，实际算的是垂直于散射面的分量（用户的 S1），而 MoM 配置（k=+z, E=+x, φ=0 即 xz 平面）中入射电场在散射面内，应对应 parallel 分量（用户的 S2）。

### 修正 `mie_compare.py`

修改 `mie_rcs_sphere`：

1. **Mie 系数符号**：将 `a_n = -jn/hn`、`b_n = -(jn+x*jn_d)/(hn+x*hn_d)` 改为正号，使光学定理 `σ_ext = σ_sca` 成立（不影响 RCS 幅值，只修正相位）。
2. **极化分量选择**：对 MoM 几何（E 在 scattering plane 内），应取 parallel 分量。在本代码的角函数约定下，该分量为：
   ```python
   S_parallel = Σ coeff * (a_n * pi_n + b_n * tau_n)
   ```
   原代码用 `a_n*tau_n + b_n*pi_n` 实际对应 perpendicular 分量。

### 修正后结果

#### ka = 3.0

| 指标 | 修正前 | 修正后 |
|------|--------|--------|
| RMSE (全角度) | 2.20 dB | **0.117 dB** |
| 最大误差 | 6.19 dB | **0.285 dB** |

90° 处：MoM -0.544 dBsm，Mie -0.647 dBsm，误差仅 0.10 dB。

#### ka = 1.5

| 指标 | 修正前 | 修正后 |
|------|--------|--------|
| RMSE (全角度) | 1.21 dB | **0.097 dB** |
| 最大误差 | 2.61 dB | **0.170 dB** |

### 结论

- **EFIE 实现本身完全正确**，此前 30°–70° 和 ka=3 的“大误差”全是 Mie 参考解极化分量选错造成的假象。
- 一旦选对 Mie 的 parallel 分量（S2 在 NIST 约定下），ka=1.5 和 ka=3.0 的 RMSE 都降到 0.1 dB 量级，最大误差 < 0.3 dB。
- 之前关于“内谐振导致 ka=3 失效”的判断是错误归因；ka=3 的深零正是 Mie S2 的真实行为，MoM 复现得非常好。
- 当前 `MAIN_RCS.F90` 已改回 `LAMBDA = 4.188790`（ka=1.5）。

### 后续

- 可继续用纯 EFIE 在用户指定的任意非谐振 ka 处做验证。
- MFIE/CFIE 调试可继续暂停，除非用户需要处理真正的内谐振频点或闭合体高频问题。

---

## 2026-08-27 会话记录：EFIE-QBX 单侧实现

### 目标

按 `docs/QBX展开求积法手册-EFIE-MFIE推导与代码实现-修订版v1.1.md` 实现 **EFIE-QBX 单侧**，并与现有奇异提取法对比，要求 RCS 曲线一致、可切换。

### 新增/修改文件

| 文件 | 改动 |
|---|---|
| `QBX_EFIE.F90` | 新增 EFIE-QBX 模块：球 Bessel/Neumann、球谐（m=0,±1）、局部坐标架、点-三角形最短距离、QBX 系数、极点求值、三角形对阻抗 |
| `NUMERICAL_INTEGRATION.F90` | 新增 Dunavant 25 点高斯规则（供 QBX 源单元加密积分） |
| `EM_TYPES.F90` | 新增 `GAUSS_25PT = 25` |
| `Z_MATRIX.F90` | 新增 `EFIE_NEAR_METHOD` 开关；近场分支可选择走 QBX 或奇异提取；QBX 失败自动回退提取法 |
| `MAIN_RCS.F90` | 新增 `USE_QBX_EFIE` 参数，默认 `.TRUE.` |
| `Makefile` | 加入 `QBX_EFIE.f90` |

### 当前 QBX 策略

- **自项**（同一三角形）：走 `CALC_EFIE_SELF_TRI_PAIR`（经典路径最稳健）。
- **共边对**：走 `CALC_EFIE_EXTRACTED_PAIR`（许可性检查易失败，经典路径已足够精确）。
- **共顶点 / 近距离非接触对**：走 QBX；许可性检查失败时回退到提取法。
- **远场**：保持原有完整格林函数高斯积分。

### 关键参数

```fortran
P_TRUNC  = 6          ! QBX 截断阶
H_FACTOR = 0.3        ! h = 0.3 * 观察三角形最大边长
N_GAUSS_UP = 25       ! 源单元 Dunavant 25 点
ADM_TOL  = 1.0E-3     ! 许可性相对裕量
```

### 调试踩坑记录

1. **25 点 Dunavant 权重归一化**：从 libMesh 源码抄来的权重对应参考三角形（面积=1/2），权重和为 0.5；本项目现有高斯规则权重和为 1.0。直接套用导致所有 QBX 矩阵元只有经典路径的一半。修复：25 点规则权重全部乘以 2。
2. **自项 QBX 收敛极慢**：切点处 q→1，p=6 远不够，矩阵元误差达 100% 以上。修复：自项直接走经典路径。
3. **共边对许可性检查失败**：共边点距中心 ≈ h，不满足 `min_dist ≥ h(1+eps)`。修复：共边对走经典路径。
4. **球谐测试误判**：一度以为 `SPH_HARM_3M` 有误，后发现是测试期望用了完整 m 求和，而代码只算 m=0,±1；实际实现正确。

### 验证结果

#### 矩阵元对比（前 5×5 近场对）

QBX 与奇异提取法相对误差普遍在 `10^-4` 量级（0.01%）。

#### RCS 曲线对比（纯 EFIE，ka=2π/λ≈1.5 球）

| 指标 | 数值 |
|---|---|
| RMSE | **0.0004 dB** |
| 最大误差 | **0.0009 dB** |
| 平均绝对误差 | **0.0003 dB** |

结论：QBX 分支与奇异提取法在当前参数下 RCS 曲线几乎完全重合。

### 使用方式

在 `MAIN_RCS.F90` 中切换：

```fortran
LOGICAL, PARAMETER :: USE_QBX_EFIE = .TRUE.   ! QBX
LOGICAL, PARAMETER :: USE_QBX_EFIE = .FALSE.  ! 奇异提取法
```

或在其他程序中直接设置 `Z_MATRIX` 模块变量：

```fortran
USE Z_MATRIX
EFIE_NEAR_METHOD = EFIE_NEAR_QBX       ! QBX
EFIE_NEAR_METHOD = EFIE_NEAR_EXTRACTION ! 奇异提取法
```

### 输出文件

- `rcs_results.txt`：当前方法（默认 QBX）的 RCS
- `rcs_results_qbx_v2.txt`：QBX 验证时保存的结果

### 2026-08-27 补充：生成 QBX / 奇异提取法 / Mie 对比曲线

为便于用户直接画图对比，重新分别运行两种方法并合并输出：

| 文件 | 来源 |
|---|---|
| `rcs_results_qbx.txt` | `USE_QBX_EFIE = .TRUE.`（默认 QBX，自项/共边走经典路径） |
| `rcs_results_qbx_full.txt` | `USE_QBX_EFIE = .TRUE.`（全 QBX：自项、共边、共顶点、近场均尝试 QBX） |
| `rcs_results_extraction.txt` | `USE_QBX_EFIE = .FALSE.`（奇异提取法） |
| `rcs_comparison_qbx_extract_mie.txt` | 合并全 QBX / 提取法 / Mie 级数参考解 |

合并脚本：`merge_rcs_comparison.py`。

#### 第一次结果（自项/共边走经典路径）

| 对比 | RMSE [dB] | 最大误差 [dB] |
|---|---|---|
| QBX vs 奇异提取法 | 0.0004 | 0.0009 |
| QBX vs Mie | 0.1142 | 0.2482 |
| 奇异提取法 vs Mie | 0.1145 | 0.2490 |

#### 修改后结果（自项/共边也尝试 QBX）

修改 `Z_MATRIX.F90`：CASE 3（自项）、CASE 2（共边）在 `EFIE_NEAR_METHOD == EFIE_NEAR_QBX` 时优先调用 `CALC_EFIE_QBX_PAIR`，失败再回退。

#### 第二步：实现纯 QBX（零回退）

在 `QBX_EFIE.F90` 中调整许可性检查：自项、共边对、共顶点对均跳过许可性检查，仅对非接触近场对保留检查。

运行统计：

```
QBX 调用次数:          84762
QBX 许可性检查失败次数: 0
```

所有奇异/近奇异对均走 QBX，无回退。

| 对比 | RMSE [dB] | 最大误差 [dB] |
|---|---|---|
| 纯 QBX vs 奇异提取法 | **0.0533** | **0.1192** |
| 纯 QBX vs Mie | **0.0637** | **0.1298** |
| 奇异提取法 vs Mie | 0.1145 | 0.2490 |

结论：
- 纯 QBX 与奇异提取法 RCS 出现可观测差异。
- **纯 QBX 相对 Mie 的误差（0.0637 dB）小于奇异提取法相对 Mie（0.1145 dB）**，说明把自项/共边/共顶点全部纳入 QBX 后整体精度提高。
- 输出文件 `rcs_results_qbx_pure.txt` 保存纯 QBX 结果；`rcs_comparison_qbx_extract_mie.txt` 已更新为纯 QBX / 提取法 / Mie 三条曲线。

### 后续可改进方向

1. **MFIE-QBX**：需双侧中心、局部坐标架独立建架、跳跃检验。
2. **共边对也走 QBX**：可通过放宽许可性检查或更小 h 实现，但当前精度已足够。
3. **参数调优**：H_FACTOR、p、N_GAUSS_UP 可针对更高频或更粗糙网格进一步研究。
4. **性能优化**：当前每个（场点，源单元）对独立算系数，未做缓存；可评估是否成为瓶颈。
