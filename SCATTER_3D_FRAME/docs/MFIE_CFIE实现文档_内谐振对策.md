# MFIE / CFIE 实现文档（内谐振问题对策）

> **定位**：本文档与《EFIE 奇异与近奇异积分处理文档（修订版 v2）》配套使用。两者共享：§2 距离驱动分类、§4 解析 $1/R$ 势、§7 高斯积分方案、§8 阈值参数。本文档新增：$\nabla G$ 的奇异性拆分、主值（PV）解析积分、MFIE 矩阵组装、jump 项、CFIE 组合。
>
> **背景**：求解器在 $ka \approx 3.14$ 附近出现异常，判定为封闭体 EFIE 内谐振。结论先行：**单独实现 MFIE 不能解决问题**（MFIE 在相同频率同样奇异），必须实现 CFIE 组合。MFIE 是 CFIE 的必需组件，故本文档以 MFIE 实现为主体，CFIE 组装为收尾。

---

## 目录

1. [物理背景：内谐振与 CFIE 原理](#1-物理背景内谐振与-cfie-原理)
2. [MFIE 方程与 RWG 离散化](#2-mfie-方程与-rwg-离散化)
3. [∇G 的 k 幂次拆分](#3-g-的-k-幂次拆分)
4. [近场分支：解析主值积分](#4-近场分支解析主值积分)
5. [远场分支](#5-远场分支)
6. [MFIE 矩阵元素组装](#6-mfie-矩阵元素组装)
7. [CFIE 组合与参数](#7-cfie-组合与参数)
8. [实现步骤（算法）](#8-实现步骤算法)
9. [代码模块对应关系](#9-代码模块对应关系)
10. [验证标准与测试清单](#10-验证标准与测试清单)

---

## 1. 物理背景：内谐振与 CFIE 原理

### 1.1 机理（逻辑链条）

1. EFIE 只约束散射体表面的**切向电场**：$\hat{\mathbf{n}} \times (\mathbf{E}^i + \mathbf{E}^s) = 0$；
2. 当频率等于"与散射体同形的 PEC 空腔"的谐振频率时，存在非零表面电流分布 $\mathbf{J}_{\text{cav}}$（腔体模），它在腔内产生谐振场、在腔外产生**零场**；
3. 于是 $\hat{\mathbf{n}} \times \mathbf{E}^s(\mathbf{J}_{\text{cav}})|_S = 0$，即 EFIE 算子有**非平凡零空间**——矩阵在该频率病态，解可任意叠加 $\mathbf{J}_{\text{cav}}$，表现为 RCS 曲线的非物理尖峰/凹陷；
4. 该谐振频率由几何唯一确定，**加密网格不会消除**，只随网格收敛到更尖锐。

### 1.2 诊断清单（确认 $ka \approx 3.14$ 是内谐振而非代码错误）

| # | 检查 | 内谐振的表现 |
|---|------|-------------|
| 1 | 矩阵条件数–频率曲线 | 谐振点处出现尖峰，网格加密后位置不变、峰更窄 |
| 2 | RCS–频率曲线 | 非物理凹陷/尖峰，位置不随网格移动 |
| 3 | 解的稳定性 | 对右端项微小扰动极度敏感 |
| 4 | 换 MFIE 算同一频点 | **同样失效**（两算子共享同一腔体模零空间）——这是确认内谐振的决定性判据 |

> **备注（易错点）**：若散射体是半径 $a$ 的球，第一内谐振应在 $ka \approx 2.744$（腔体 TM$_{011}$ 模），其次 $ka \approx 4.493$（TE$_{011}$）；$ka = \pi$ 并**不是**球的腔体谐振。若在球上于 $ka \approx 3.14$ 看到异常，先按上表排查是否为实现错误（如相邻对精度、$G_{\text{smooth}}$ 对角点），不要直接归因内谐振。非球几何则以其内腔尺寸估算谐振点。

### 1.3 为什么必须是 CFIE

- MFIE 算子 $\frac{1}{2}\mathcal{I} - \mathcal{K}$ 在**相同**的内谐振频率奇异（零空间对应同一批腔体模），单独换用 MFIE 只是换了一种失效方式，且 MFIE 的 $1/R^2$ 核对积分精度更敏感；
- CFIE 组合 $\alpha\,\text{EFIE} + (1-\alpha)\,\eta_0\,\text{MFIE}$（$\alpha \in (0,1)$ 实数）对应的内部问题是**带阻抗边界条件的腔体**，其本征频率含非零虚部，对任何实频率算子均非奇异——内谐振被彻底消除；
- MFIE 仅适用于**封闭体**（开放结构无 $\hat{\mathbf{n}} \times \mathbf{H}$ 方程可写）。

---

## 2. MFIE 方程与 RWG 离散化

### 2.1 方程

时间约定 $e^{j\omega t}$，$G = e^{-jkR}/(4\pi R)$。MFIE（从散射体外侧趋近表面）：

$$\frac{1}{2}\mathbf{J}(\mathbf{r}) - \hat{\mathbf{n}} \times \mathrm{PV}\!\int_S \nabla G(\mathbf{r},\mathbf{r}') \times \mathbf{J}(\mathbf{r}')\, dS' = \hat{\mathbf{n}} \times \mathbf{H}^i(\mathbf{r})$$

其中 $1/2$ 为 jump 项（$\nabla G$ 积分的法向跃变），PV 为柯西主值，

$$\nabla G = -\frac{(1 + jkR)\, e^{-jkR}}{4\pi R^3}\,(\mathbf{r} - \mathbf{r}')$$

> **符号警告**：jump 项与 PV 项的符号依赖时间约定与 $\nabla G$/$\nabla' G$ 约定。上式对应 $e^{j\omega t}$ + $G = e^{-jkR}/(4\pi R)$ + 外法向 $\hat{\mathbf{n}}$。实现后以 §10 测试 2、6 校验符号。

### 2.2 核展开

RWG 展开 $\mathbf{J} = \sum_n I_n \mathbf{f}_n$，用 $\mathbf{f}_m$ 伽辽金测试。利用 BAC–CAB 恒等式（$\mathbf{f}_n$ 在源面片内，$\hat{\mathbf{n}} = \hat{\mathbf{n}}_f$ 为**场**三角形法向）：

$$\hat{\mathbf{n}}_f \times (\nabla G \times \mathbf{f}_n) = \nabla G\,(\hat{\mathbf{n}}_f \cdot \mathbf{f}_n) - \mathbf{f}_n\,(\hat{\mathbf{n}}_f \cdot \nabla G)$$

矩阵元素：

$$Z^M_{mn} = \underbrace{\frac{1}{2}\langle \mathbf{f}_m, \mathbf{f}_n \rangle}_{\text{jump 项}} - C_M C_N \iint \Big[ (\tilde{\mathbf{f}}_m \cdot \nabla G)(\hat{\mathbf{n}}_f \cdot \tilde{\mathbf{f}}_n) - (\tilde{\mathbf{f}}_m \cdot \tilde{\mathbf{f}}_n)(\hat{\mathbf{n}}_f \cdot \nabla G) \Big] dS' dS$$

其中 $\tilde{\mathbf{f}}_m = \mathbf{r} - \mathbf{r}_{\text{opp}}^M$，$\tilde{\mathbf{f}}_n = \mathbf{r}' - \mathbf{r}_{\text{opp}}^N$（$M$ 源、$N$ 场，与 EFIE 文档约定一致）。

右端项：$V^M_m = \langle \mathbf{f}_m,\ \hat{\mathbf{n}} \times \mathbf{H}^i \rangle$；平面波 $\mathbf{H}^i = (\hat{\mathbf{k}} \times \mathbf{E}^i)/\eta_0$。

### 2.3 Jump 项（解析）

$\frac{1}{2}\langle \mathbf{f}_m, \mathbf{f}_n \rangle$ 仅当两 RWG 支撑重叠（同一 RWG，或共用同一三角形的 RWG 对）时非零，用重心坐标积分解析计算：

$$\int_T \lambda_i \lambda_j\, dS = \frac{A}{12}(1 + \delta_{ij}), \qquad \mathbf{r} - \mathbf{r}_{\text{opp}} = \sum_i \lambda_i (\mathbf{v}_i - \mathbf{r}_{\text{opp}})$$

---

## 3. ∇G 的 k 幂次拆分

利用 $(1 + jx)e^{-jx} = 1 + \frac{x^2}{2} - \frac{jx^3}{3} - \frac{x^4}{8} + \frac{jx^5}{30} + \frac{x^6}{144} + \cdots$（$x = kR$），得

$$\nabla G = \frac{1}{4\pi}\Big[ \underbrace{-\frac{\mathbf{r}-\mathbf{r}'}{R^3}}_{\text{静态主值项}} \;\underbrace{-\frac{k^2}{2}\frac{\mathbf{r}-\mathbf{r}'}{R}}_{k^2\text{ 修正项}} \;+ \underbrace{k^3(\mathbf{r}-\mathbf{r}')\Big[\frac{kR}{8}\Big(1-\frac{(kR)^2}{18}\Big) + \frac{j}{3}\Big(1-\frac{(kR)^2}{10}\Big)\Big]}_{\text{正则余项 } \mathbf{F}_{\text{reg}}} + O(k^6 R^3) \Big]$$

| 部分 | 奇异强度 | 处理 |
|------|---------|------|
| 静态项 $-(\mathbf{r}-\mathbf{r}')/R^3$ | $1/R^2$ 强奇异（需 PV） | **解析**（§4.2） |
| $k^2$ 修正项 $-\frac{k^2}{2}(\mathbf{r}-\mathbf{r}')/R$ | $1/R$ 弱奇异 | **解析**（§4.3，复用 EFIE 文档 §4 势函数 + 二阶矩） |
| 正则余项 $\mathbf{F}_{\text{reg}}$ | 有界（$O(k^3 R)$） | 数值高斯（§4.4） |

> **归一化警告（承接 EFIE 文档 §3）**：老师版 `gradG` 不含 $1/4\pi$，且其装配 `gradG = k³I_30 − (I_31 + k²/2·I_32)/An` 中面积归一化在外层约定不同。移植时逐项核对 $4\pi$ 与面积因子，这是两套代码对接最高频的错误源。

---
## 4. 近场分支：解析主值积分

### 4.1 内层矩定义

对每个场高斯点 $\mathbf{r} \in T_f$，需要源三角形 $T_s$ 上的三个内层矩（$\nabla G$ 按其拆分逐项计算）：

| 矩 | 定义 | 类型 |
|------|------|------|
| $\mathbf{g}_A(\mathbf{r})$ | $\int_{T_s} \nabla G\, dS'$ | 矢量 |
| $g_B(\mathbf{r})$ | $\hat{\mathbf{n}}_f \cdot \mathbf{g}_A$ | 标量（由 $\mathbf{g}_A$ 收缩） |
| $\mathbf{g}_C(\mathbf{r})$ | $\int_{T_s} \mathbf{r}'\, [\hat{\mathbf{n}}_f \cdot \nabla G]\, dS'$ | 矢量 |
| $\mathbf{g}_E(\mathbf{r})$ | $\int_{T_s} (\hat{\mathbf{n}}_f \cdot \mathbf{r}')\, \nabla G\, dS'$ | 矢量 |

### 4.2 静态主值项（解析）

全部在源三角形局部系中计算，**复用 EFIE 文档 §4 的输出**：$I_0$、$\mathbf{I}_{\text{shifted}}$、$\hat{\mathbf{m}}_e^{\text{out}}$、$d$、$\mathbf{p}$、$|\Omega|$（§4.3），并定义 $\sigma(d) = \mathrm{sgn}(d)$（**规定 $\sigma(0) = 0$**），$\hat{\mathbf{n}}_{f,\parallel} = \hat{\mathbf{n}}_f - (\hat{\mathbf{n}}_f \cdot \hat{\mathbf{n}}_s)\hat{\mathbf{n}}_s$。

**所需边线积分原语**（均为闭式）：

| 原语 | 定义 | 闭式 |
|------|------|------|
| $g_{2,e}$ | $\int_e \frac{dl}{R}$ | $\ln\dfrac{R_e^+ + s_e^+}{R_e^- + s_e^-}$（EFIE 文档 §4 已有） |
| $\mathbf{L}_e$ | $\int_e \frac{\mathbf{r}'}{R} dl$ | $\mathbf{v}_e^- g_{2,e} + \hat{\mathbf{t}}_e\big[(R_e^+ - R_e^-) + s_{\text{start},e}\, g_{2,e}\big]$ |

**三个静态主值矩**（推导恒等式见本节末注）：

$$\mathbf{g}_A^{\text{st}} = -\sum_e \hat{\mathbf{m}}_e^{\text{out}}\, g_{2,e} \;-\; \hat{\mathbf{n}}_s\, \sigma(d)\, |\Omega|$$

$$\mathbf{g}_C^{\text{st}} = (\hat{\mathbf{n}}_f\cdot\hat{\mathbf{n}}_s)\big[ -\sigma(d)|\Omega|\, \mathbf{p} + d \sum_e \hat{\mathbf{m}}_e^{\text{out}} g_{2,e} \big] - \sum_e (\hat{\mathbf{n}}_{f,\parallel}\cdot\hat{\mathbf{m}}_e^{\text{out}})\, \mathbf{L}_e + \hat{\mathbf{n}}_{f,\parallel}\, I_0$$

$$\mathbf{g}_E^{\text{st}} = (\hat{\mathbf{n}}_f\cdot\hat{\mathbf{n}}_s)(\hat{\mathbf{n}}_s\cdot\mathbf{p})\, \mathbf{g}_A^{\text{st}} - \sum_e (\hat{\mathbf{n}}_{f,\parallel}\cdot\hat{\mathbf{m}}_e^{\text{out}})\big(\mathbf{L}_e - \mathbf{p}\, g_{2,e}\big) + \hat{\mathbf{n}}_{f,\parallel}\, I_0 + d\,\big(\hat{\mathbf{n}}_{f,\parallel}\cdot\textstyle\sum_e \hat{\mathbf{m}}_e^{\text{out}} g_{2,e}\big)\,\hat{\mathbf{n}}_s$$

**退化情形（重要，覆盖自作用与共面相邻对）**：$d = 0$ 且 $\hat{\mathbf{n}}_{f,\parallel} = \mathbf{0}$ 时，

$$\mathbf{g}_C^{\text{st}} = \mathbf{0}, \qquad g_B^{\text{st}} = \hat{\mathbf{n}}_f \cdot \mathbf{g}_A^{\text{st}} = 0$$

即 MFIE 自作用与共面相邻对的 PV 贡献恒为零，只剩 jump 项——这是平面面片 MFIE 的经典结论，可作为实现的强校验（§10 测试 2）。此时代码可直接跳过 PV 计算（优化）。

> **推导恒等式**（供核查）：$\nabla_\mathbf{r}(1/R) = (\mathbf{r}' - \mathbf{r})/R^3 = (\mathbf{u}' - d\hat{\mathbf{n}}_s)/R^3$（$\mathbf{u}' = \mathbf{r}' - \mathbf{p}$ 面内）；$\int \mathbf{u}'/R^3\, dS' = -\oint (1/R)\hat{\mathbf{m}}^{\text{out}} dl$；$-d\int 1/R^3 dS' = -\sigma(d)|\Omega|$；$\int \mathbf{r}'_i\, \partial'_j(1/R) dS' = \oint (1/R) r'_i m_j dl - \delta^{\parallel}_{ij} I_0$。所有含 $1/R^3$ 的量都通过以上恒等式化为正则边线积分或立体角，PV 意义自动满足。

### 4.3 $k^2$ 修正项（解析）

核为 $-\frac{k^2}{2}(\mathbf{r}-\mathbf{r}')/R$（弱奇异）。除 EFIE 文档 §4 的矢量势外，需要**一个新原语：二阶矩** $\mathbf{U}_2 = \int_{T_s} \mathbf{u}' \otimes \mathbf{u}'/R\, dS'$（面内对称 $2\times2$ 张量，3 个独立分量）：

$$\boxed{\;(\mathbf{U}_2)_{ij} = \sum_e m^{\text{out}}_{e,i}\, \mathbf{H}_{e,j} - \delta^{\parallel}_{ij}\, J_R\;}$$

其中新增闭式原语：

| 原语 | 定义 | 闭式 |
|------|------|------|
| $G_{3,e}$ | $\int_e R\, dl$ | $\frac12\big[s_e^+ R_e^+ - s_e^- R_e^- + \rho_e^2\, g_{2,e}\big]$（即 EFIE 文档 $\mathbf{I}_{\text{shifted}}$ 的边项） |
| $H_{3,e}$ | $\int_e s' R\, dl$ | $\frac{1}{3}\big[(R_e^+)^3 - (R_e^-)^3\big] + s_{\text{start},e}\, G_{3,e}$ |
| $\mathbf{H}_e$ | $\int_e \mathbf{u}' R\, dl$ | $(\mathbf{v}_e^- - \mathbf{p})\, G_{3,e} + \hat{\mathbf{t}}_e H_{3,e}$ |
| $J_R$ | $\int_{T_s} R\, dS'$ | $\frac13\Big[\sum_e d_e^{\text{out}}\, G_{3,e} + d^2 I_0\Big]$，$d_e^{\text{out}} = (\mathbf{p}-\mathbf{v}_e^-)\cdot\hat{\mathbf{m}}_e^{\text{out}}$ |

$k^2$ 修正矩（$\mathbf{c} = \hat{\mathbf{n}}_{f,\parallel}$）：

$$\mathbf{g}_A^{k^2} = -\frac{k^2}{2}\big( d\, I_0\, \hat{\mathbf{n}}_s - \mathbf{I}_{\text{shifted}} \big)$$

$$\mathbf{g}_C^{k^2} = -\frac{k^2}{2}\Big[ (\hat{\mathbf{n}}_f\cdot\mathbf{r})\, \mathbf{I}_{\text{vec}} - (\hat{\mathbf{n}}_f\cdot\hat{\mathbf{n}}_s)(\hat{\mathbf{n}}_s\cdot\mathbf{p})\, \mathbf{I}_{\text{vec}} - \mathbf{p}\,(\mathbf{c}\cdot\mathbf{I}_{\text{shifted}}) - \mathbf{U}_2\,\mathbf{c} \Big]$$

$$\mathbf{g}_E^{k^2} = -\frac{k^2}{2}\Big[ (\hat{\mathbf{n}}_f\cdot\hat{\mathbf{n}}_s)(\hat{\mathbf{n}}_s\cdot\mathbf{p})\big(d\,I_0\,\hat{\mathbf{n}}_s - \mathbf{I}_{\text{shifted}}\big) + d\,\hat{\mathbf{n}}_s\,(\mathbf{c}\cdot\mathbf{I}_{\text{shifted}}) - \mathbf{U}_2\,\mathbf{c} \Big]$$

> 注：自作用时 $d=0$、$\mathbf{c}=\mathbf{0}$，$\mathbf{g}_C^{k^2} = \mathbf{0}$、$\hat{\mathbf{n}}_s\cdot\mathbf{g}_A^{k^2} = 0$，与 §4.2 退化结论一致（$k^2$ 项同样不贡献自作用 PV）。

### 4.4 正则余项（数值）

对源三角形做高斯积分（12 点），核有界：

$$\mathbf{g}_{\{\cdot\}}^{\text{reg}} = A_s \sum_q w_q\; \{\cdot\}\text{-加权}\; \mathbf{F}_{\text{reg}}(\mathbf{r}, \mathbf{r}'_q), \qquad \mathbf{g}_C^{\text{reg}} = A_s \sum_q w_q\, \mathbf{r}'_q\,(\hat{\mathbf{n}}_f\cdot\mathbf{F}_{\text{reg}}), \quad \mathbf{g}_E^{\text{reg}} = A_s \sum_q w_q\, (\hat{\mathbf{n}}_f\cdot\mathbf{r}'_q)\,\mathbf{F}_{\text{reg}}$$

$\mathbf{F}_{\text{reg}}$ 按 §3 公式逐项求值（无 $0/0$ 点，无需特殊分支）。

### 4.5 汇总与归一化

$$\mathbf{g}_{\{\cdot\}} = \frac{1}{4\pi}\big( \mathbf{g}_{\{\cdot\}}^{\text{st}} + \mathbf{g}_{\{\cdot\}}^{k^2} + \mathbf{g}_{\{\cdot\}}^{\text{reg}} \big)$$

§4.2/§4.3 的解析公式基于 $1/R$ 势（不含 $4\pi$），统一在最后乘 $1/(4\pi)$；正则余项公式已含 $k^3$ 尺度，同样除以 $4\pi$。

---

## 5. 远场分支

完整 $\nabla G$ 直接求值，双重高斯积分：

$$\nabla G = -\frac{(1+jkR)e^{-jkR}}{4\pi R^3}(\mathbf{r}-\mathbf{r}')$$

内层矩 $\mathbf{g}_A, \mathbf{g}_C, \mathbf{g}_E$ 按 §4.1 定义直接数值计算，外层矩与组装（§6）与近场分支共用同一代码。

---

## 6. MFIE 矩阵元素组装

### 6.1 外层矩（场三角形 $T_f$，12 点高斯）

| 矩 | 定义 | 类型 |
|------|------|------|
| $\mathbf{P}_1$ | $A_f \sum_i w_i\, \mathbf{g}_A(\mathbf{r}_i)$ | 矢量 |
| $P_2$ | $A_f \sum_i w_i\, \mathbf{r}_i\cdot\mathbf{g}_A$ | 标量 |
| $\mathbf{P}_3$ | $A_f \sum_i w_i\, \mathbf{g}_E$ | 矢量 |
| $P_4$ | $A_f \sum_i w_i\, \mathbf{r}_i\cdot\mathbf{g}_E$ | 标量 |
| $\mathbf{P}_5$ | $A_f \sum_i w_i\, \mathbf{g}_C$ | 矢量 |
| $P_6$ | $A_f \sum_i w_i\, \mathbf{r}_i\cdot\mathbf{g}_C$ | 标量 |
| $P_7$ | $A_f \sum_i w_i\, g_B$ | 标量 |
| $\mathbf{P}_8$ | $A_f \sum_i w_i\, \mathbf{r}_i\, g_B$ | 矢量 |

### 6.2 组装公式

$$A_{\text{项}} = P_4 - (\hat{\mathbf{n}}_f\cdot\mathbf{r}_{\text{opp}}^N)\, P_2 - \mathbf{r}_{\text{opp}}^M\cdot\mathbf{P}_3 + (\hat{\mathbf{n}}_f\cdot\mathbf{r}_{\text{opp}}^N)(\mathbf{r}_{\text{opp}}^M\cdot\mathbf{P}_1)$$

$$B_{\text{项}} = P_6 - \mathbf{r}_{\text{opp}}^N\cdot\mathbf{P}_8 - \mathbf{r}_{\text{opp}}^M\cdot\mathbf{P}_5 + (\mathbf{r}_{\text{opp}}^M\cdot\mathbf{r}_{\text{opp}}^N)\, P_7$$

$$Z^M_{\text{pair}} = \frac{1}{2}\langle \mathbf{f}_m, \mathbf{f}_n \rangle_{T_f}\; - \; C_M C_N\,(A_{\text{项}} - B_{\text{项}})$$

jump 项按 §2.3 解析，仅在两 RWG 支撑重叠的三角形对上非零；PV 项对全部三角形对按分支计算。

> **配对规则**（与 EFIE 文档 §5.2 同源）：展开 $\tilde{\mathbf{f}}_m \cdot \tilde{\mathbf{f}}_n$ 时，源对顶点 $\mathbf{r}_{\text{opp}}^M$ 配场坐标矩，场对顶点 $\mathbf{r}_{\text{opp}}^N$ 配源坐标矩；上式已按此规则展开，实现后不得调换。

---

## 7. CFIE 组合与参数

$$Z^{\text{CFIE}} = \alpha\, Z^E + (1-\alpha)\,\eta_0\, Z^M, \qquad \mathbf{V}^{\text{CFIE}} = \alpha\, \mathbf{V}^E + (1-\alpha)\,\eta_0\, \mathbf{V}^M$$

| 参数 | 建议值 | 说明 |
|------|--------|------|
| $\alpha$ | 0.2（默认） | 常用范围 0.2–0.5；$\alpha$ 越小 MFIE 权重越大，精度略降但迭代收敛略好 |
| $\eta_0$ 因子 | 必须保留 | 量纲匹配：$Z^E$ 为欧姆量级，$Z^M$ 无量纲，$\eta_0\mathbf{H}$ 与 $\mathbf{E}$ 同量纲——**漏乘 $\eta_0$ 是 CFIE 实现最高频错误** |

MFIE 分支/阈值与 EFIE 共用（$d_{\text{th}} = 0.1\lambda$，同一对三角形两算子走同一分支，与老师版 `dist>=threshold_R` 同时控制 `Comp_eltsI`/`Comp_eltsU` 的做法一致）。

---
## 8. 实现步骤（算法）

### 8.1 静态主值内层 `ANALYTIC_GRAD_STATIC_PV`

输入：场点 $\mathbf{r}$、场法向 $\hat{\mathbf{n}}_f$、源三角形顶点。输出：$\mathbf{g}_A^{\text{st}}, \mathbf{g}_C^{\text{st}}, \mathbf{g}_E^{\text{st}}$。

1. 复用标量/矢量势例程的几何量：$\hat{\mathbf{n}}_s$、$d$、$\mathbf{p}$、$I_0$、$|\Omega|$、各边 $\hat{\mathbf{t}}_e, \hat{\mathbf{m}}_e^{\text{out}}, s_e^\pm, R_e^\pm$；
2. 每边算原语 $g_{2,e}$、$\mathbf{L}_e$（§4.2 表）；
3. 退化检测：若 $|d| < \varepsilon_{\log}$ 且 $\|\hat{\mathbf{n}}_{f,\parallel}\| < 10^{-12}$，直接返回 $\mathbf{g}_A^{\text{st}} = -\sum \hat{\mathbf{m}}_e^{\text{out}} g_{2,e}$、$\mathbf{g}_C^{\text{st}} = \mathbf{0}$、$\mathbf{g}_E^{\text{st}} = (\hat{\mathbf{n}}_s\cdot\mathbf{p})\,\mathbf{g}_A^{\text{st}}$；
4. 否则按 §4.2 三式组装（$\sigma(d)$ 取符号函数，$\sigma(0)=0$）。

### 8.2 $k^2$ 修正内层 `ANALYTIC_GRAD_K2`

1. 复用 $\mathbf{I}_{\text{shifted}}, \mathbf{I}_{\text{vec}}, I_0$；
2. 每边算 $G_{3,e}$、$H_{3,e}$、$\mathbf{H}_e$；累加得 $J_R$ 与 $\mathbf{U}_2$（§4.3 表与公式，注意 $\delta^{\parallel}_{ij}$ 只作用于面内分量）；
3. 按 §4.3 三式组装。

### 8.3 近场分支 `CALC_MFIE_NEAR_TRI_PAIR(TRI_M, TRI_N, RWG_M, RWG_N, K, ETA0, ZM_PAIR)`

1. 提取 $C_M, \mathbf{r}_{\text{opp}}^M, C_N, \mathbf{r}_{\text{opp}}^N$ 与场法向 $\hat{\mathbf{n}}_f$；
2. 外层 12 点：每点 $\mathbf{r}_i$ 调 §8.1、§8.2 与 §4.4 正则余项数值积分，合成 $\mathbf{g}_A, \mathbf{g}_C, \mathbf{g}_E$（乘 $1/4\pi$），$g_B \leftarrow \hat{\mathbf{n}}_f\cdot\mathbf{g}_A$；
3. 累加外层矩 $\mathbf{P}_1 \sim \mathbf{P}_8$（乘 $A_f$ 与权重）；
4. 按 §6.2 组装 PV 部分；若两 RWG 共用本三角形，加 jump 项（§2.3）。

### 8.4 远场分支与入口

1. 远场：双重高斯 + 完整 $\nabla G$，内层矩数值计算，§6 组装共用；
2. 入口 `CALC_MFIE_MATRIX_ELEMENT`：与 EFIE 共用同一分支判据（§2 分类，$N_{\text{share}} \geq 1$ 或 $d_c < d_{\text{th}}$）；
3. CFIE 开关：`IF (CFIE) Z = α·Z_E + (1-α)·η0·Z_M`（右端项同），建议编译期/输入文件双开关 `USE_MFIE`、`USE_CFIE(α)`。

### 8.5 分期实施建议

| 阶段 | 内容 | 出口判据 |
|------|------|---------|
| 1 | jump 项 + 远场 $\nabla G$ + CFIE 组装（近场暂用高阶数值凑） | 非谐振频点 CFIE 与 EFIE 一致（§10 测试 7） |
| 2 | §8.1/§8.2 解析内层替换数值 | §10 测试 1–5 全过 |
| 3 | 球扫频验收 | §10 测试 6（谐振点条件数与 RCS） |

---

## 9. 代码模块对应关系

| 模块文件 | 函数 | 功能 | 备注 |
|---------|------|------|------|
| `GREEN_FUNCTIONS.F90` | `GRAD_GREEN_FUNC` | 完整 $\nabla G$（§5） | 原版拼写 `GARD_` 建议改名 |
| | `GRAD_GREEN_REG` | 正则余项 $\mathbf{F}_{\text{reg}}$（§3） | 新增 |
| `SINGULAR_INTEGRAL.F90` | `ANALYTIC_GRAD_STATIC_PV` | 静态主值三矩（§4.2） | 新增 |
| | `ANALYTIC_GRAD_K2` | $k^2$ 修正三矩（§4.3） | 新增；内含 $J_R$、$\mathbf{U}_2$ |
| `Z_MATRIX.F90` | `RWG_OVERLAP_INTEGRAL` | jump 项 $\langle \mathbf{f}_m,\mathbf{f}_n\rangle$（§2.3） | 新增 |
| | `CALC_MFIE_NEAR_TRI_PAIR` | MFIE 近场对（§8.3） | 新增 |
| | `CALC_MFIE_MATRIX_ELEMENT` | MFIE 统一入口 | 新增 |
| | CFIE 组装 | $Z^C, V^C$（§7） | 新增 |
| `test_singular.F90` | MFIE 测试组 | §10 | 新增 |

---

## 10. 验证标准与测试清单

**单元测试：**

| # | 测试 | 方法 | 通过标准 |
|---|------|------|---------|
| 1 | $\mathbf{g}_A^{\text{st}}$（一般面外点，含 $\hat{\mathbf{n}}_{f,\parallel}\neq0$ 情形） | vs 加密数值 PV 参考 | 各分量相对误差 $<0.1\%$ |
| 2 | 自作用/共面退化：$\mathbf{g}_C \equiv \mathbf{0}$、$g_B \equiv 0$、$\hat{\mathbf{n}}_s\cdot\mathbf{F}_{\text{reg}}\equiv0$ | 直接检验 | $<10^{-12}$（相对量级） |
| 3 | $J_R = \int R\,dS'$、$(\mathbf{U}_2)_{ij}$ | vs 高阶数值 | $<0.1\%$；$\mathbf{U}_2$ 对称性到机器精度 |
| 4 | 远点渐近：$\mathbf{g}_A \to \nabla(A/R)$ | 解析极限 | $<0.1\%$（$|\mathbf{r}|\gg h$） |
| 5 | 阈值连续性（MFIE 版）：同一对横跨 $d_{\text{th}}$ 强制两分支各算一次 | 相互对比 | 相对差 $<0.1\%$ |

**物理验收（直接针对 $ka\approx3.14$ 问题）：**

| # | 测试 | 通过标准 |
|---|------|---------|
| 6 | 封闭体扫频（覆盖问题频点）：EFIE / MFIE / CFIE 三曲线 | EFIE 与 MFIE 条件数在同一频率尖峰（确认内谐振）；CFIE 条件数全程平稳 |
| 7 | 非谐振频点：CFIE vs EFIE 的 RCS | 偏差 $<0.5\%$ |
| 8 | 谐振频点：CFIE 的 RCS vs 参考解（球用 Mie 级数） | 偏差 $<1\%$ |
| 9 | $\alpha$ 无关性：$\alpha = 0.2$ 与 $0.5$ 全频段结果对比 | 偏差 $<0.5\%$ |
| 10 | EFIE 原有 12 项测试回归 | 全过（MFIE 不得破坏 EFIE 路径） |

> **测试 6 的读法**：若 CFIE 在谐振点仍异常，先查 §7 的 $\eta_0$ 因子与 jump 项符号（两处最高频错误），再查 $\mathbf{g}_C/\mathbf{g}_E$ 的 $4\pi$ 归一化。

---

## 参考文献

1. Rao, S. M., Wilton, D. R., and Glisson, A. W. "Electromagnetic scattering by surfaces of arbitrary shape." *IEEE Trans. Antennas Propag.*, 1982.
2. Graglia, R. D. "On the numerical integration of the linear shape functions times the 3-D Green's function or its gradient on a plane triangle." *IEEE Trans. Antennas Propag.*, 1993.（$\nabla G$ 主值处理的原始出处）
3. Wilton, D. R., et al. "Potential integrals for uniform and linear source distributions on polygonal and polyhedral domains." *IEEE Trans. Antennas Propag.*, 1984.
4. 老师版求解器文档（`Comp_eltsU` / `Spara_cal_2`：$\nabla G$ 的 $k$ 幂次拆分与 $w_0$ 分支结构来源）。
