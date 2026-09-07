# EFIE 奇异与近奇异积分处理文档（修订版 v2）

> **本版修订要点**
> 1. 分类准则由"拓扑分类（同一/共边/其他）"改为"**距离驱动 + 共顶点兜底**"，自作用、共边、共顶点、近距离三角形对统一走解析提取分支；
> 2. 修正原版 §5.2 中 $\mathbf{r}_{\text{opp}}$ 配对规则的自相矛盾（见 §5.2）；
> 3. 统一 $Q$ 算子符号约定，全文唯一（见 §1）；
> 4. `G_smooth` 在 $R=0$ 处显式赋值 $-jk/(4\pi)$，禁止返回 0（见 §3）；
> 5. 更正原版 Bug #9 的影响归因（自作用场点恒在面内，立体角项对其无影响，见 §12）。

---

## 目录

1. [背景：EFIE 矩量法公式与符号约定](#1-背景efie-矩量法公式与符号约定)
2. [奇异性来源与统一处理策略](#2-奇异性来源与统一处理策略)
3. [格林函数拆分与正则余项求值](#3-格林函数拆分与正则余项求值)
4. [解析 1/R 势：边界积分公式](#4-解析-1r-势边界积分公式)
5. [近奇异/奇异分支：统一阻抗计算](#5-近奇异奇异分支统一阻抗计算)
6. [远区数值积分](#6-远区数值积分)
7. [数值积分方案](#7-数值积分方案)
8. [阈值参数汇总](#8-阈值参数汇总)
9. [代码模块对应关系](#9-代码模块对应关系)
10. [实现逻辑（算法步骤）](#10-实现逻辑算法步骤)
11. [验证标准与测试清单](#11-验证标准与测试清单)
12. [与原版的差异汇总](#12-与原版的差异汇总)

---

## 1. 背景：EFIE 矩量法公式与符号约定

对 PEC 散射问题，EFIE 经伽辽金测试（时间约定 $e^{j\omega t}$）后，阻抗矩阵元素为：

$$Z_{mn} = j\omega\mu \underbrace{\iint \mathbf{f}_m(\mathbf{r}) \cdot \mathbf{f}_n(\mathbf{r}') G(\mathbf{r},\mathbf{r}') dS' dS}_{P\text{ 项（矢量势）}} + \frac{1}{j\omega\varepsilon} \underbrace{\iint (\nabla \cdot \mathbf{f}_m(\mathbf{r})) (\nabla' \cdot \mathbf{f}_n(\mathbf{r}')) G(\mathbf{r},\mathbf{r}') dS' dS}_{\hat{Q}\text{ 项（标量势）}}$$

自由空间格林函数（**含 $1/4\pi$ 归一化**）：

$$G(\mathbf{r},\mathbf{r}') = \frac{e^{-jkR}}{4\pi R}, \quad R = |\mathbf{r} - \mathbf{r}'|, \quad k = \frac{2\pi}{\lambda}$$

波数与波阻抗换算：

$$j\omega\mu = jk\eta_0, \quad \frac{1}{j\omega\varepsilon} = -\frac{j\eta_0}{k}, \quad \eta_0 = 120\pi\ \Omega$$

> **符号约定（全文唯一，修订点）**：代码与本文档中 $Q$ 算子定义为**带负号**的散度积分
>
> $$Q := -\hat{Q} = -\iint (\nabla \cdot \mathbf{f}_m)(\nabla' \cdot \mathbf{f}_n) G\, dS' dS$$
>
> 从而组装式中两项系数形式统一为 $+j$：
>
> $$Z_{mn} = jk\eta_0\, P + \frac{j\eta_0}{k}\, Q$$
>
> 后续所有 $Q_{\text{sing}}$、$Q_{\text{smooth}}$ 均按此带负号定义给出，不再中途变号。物理检验：自作用 $Q$ 的 $1/R$ 部分为负实数，乘以 $j\eta_0/k$ 得负虚部，对应电容性贡献，方向正确。

RWG 基函数在三角形 $T$ 上：

$$\mathbf{f}(\mathbf{r}) = C \cdot (\mathbf{r} - \mathbf{r}_{\text{opp}}), \quad C = \pm \frac{L}{2A}$$

正三角形取 $C = +L/2A$，负三角形取 $-L/2A$，$\mathbf{r}_{\text{opp}}$ 为公共边对顶点，散度 $\nabla \cdot \mathbf{f} = 2C$ 为常数。

**全文指标约定（修订点）**：$M$ 恒为**源** RWG（坐标 $\mathbf{r}'$），$N$ 恒为**场** RWG（坐标 $\mathbf{r}$）。原版 §5.2 曾在一节之内混用两种约定，本版统一。

---

## 2. 奇异性来源与统一处理策略

$G \propto 1/R$ 的奇异性强弱取决于**源场点间最小距离相对于三角形尺寸**，而不是两三角形是否共享顶点。因此分类准则改为距离驱动：

**判据（按优先级）：**

1. 两三角形共用顶点数 $N_{\text{share}} \geq 1$（含同一三角形 $N_{\text{share}}=3$、共边 $=2$、共顶点 $=1$）→ **近场分支**；
2. 否则，重心距 $d_c = \|\mathbf{r}_c^{M} - \mathbf{r}_c^{N}\| < d_{\text{th}}$（默认 $d_{\text{th}} = 0.1\lambda$）→ **近场分支**；
3. 其余 → **远场分支**。

| 分支 | 覆盖情形 | 处理方法 | 高斯阶数 |
|------|---------|---------|---------|
| 近场（奇异/近奇异） | 自作用、共边、共顶点、近距离对 | 内层 $1/R$ 解析（§4）+ 外层高斯 + $G_{\text{smooth}}$ 双重数值（§5） | 外层 12 点，$G_{\text{smooth}}$ 12×12 |
| 远场（非奇异） | 其余 | 完整 $G$ 双重高斯（§6） | 调用者指定 |

> **与原版差异**：原版按拓扑分类，自作用用解析提取、共边对用 12 点乘积高斯蛮力、共顶点对未分类、近距离但不接触的对落入普通高斯。乘积高斯对弱奇异核只有代数级慢收敛，且近距离对的精度完全不可控。本版将这些情形统一纳入近场分支，自作用仅作为 $w_0 = 0$ 的特例存在，不再有独立代码路径。
>
> **判据改进备注**：重心距对大三角形可能漏判（重心相距远但边很近）。如需更稳，可将判据 2 替换为两三角形最小间距 $d_{\min}$，或归一化判据 $d_c / \max(h_M, h_N) < 0.6$（$h$ 为三角形特征尺寸）。初版实现用重心距即可，保留接口以便替换。

---

## 3. 格林函数拆分与正则余项求值

奇异性提取：

$$G(\mathbf{r},\mathbf{r}') = \underbrace{\frac{1}{4\pi R}}_{G_{\text{sing}}\text{（解析处理，§4）}} + \underbrace{\frac{e^{-jkR} - 1}{4\pi R}}_{G_{\text{smooth}}\text{（数值积分）}}$$

$G_{\text{smooth}}$ 在 $R \to 0$ 处有限：$G_{\text{smooth}}(0) = -jk/(4\pi)$。

**求值方式（两种，实现时任选其一，建议方式 B）：**

**方式 A：直接求值（任意 $kR$ 精确）**

- $R \geq \varepsilon_R$：$G_{\text{smooth}} = (e^{-jkR}-1)/(4\pi R)$；
- $R < \varepsilon_R$：**必须显式赋值** $G_{\text{smooth}} = -jk/(4\pi)$（$\varepsilon_R$ 建议取 $10^{-12} \cdot h$，$h$ 为三角形特征尺寸）。

> **修订点**：原版代码对 $R < 10^{-10}$ 返回 $(0,0)$，静默丢弃对角虚部 $-jk/(4\pi)$ 的贡献。自作用 $G_{\text{smooth}}$ 双重积分的 $12\times12$ 点对中有 12 对 $R \equiv 0$（同一组高斯点自 pairing），必须走显式赋值路径。注意方式 A 在 $kR \ll 1$ 时存在浮点相消（$e^{-jkR}-1$ 损失有效数字），$kR \gtrsim 10^{-4}$ 时可忽略，更小则需方式 B。

**方式 B：泰勒多项式求值（近场分支内推荐，适用 $kR \lesssim 2$）**

$$G_{\text{smooth}}(R) = \frac{k}{4\pi} \sum_{i=0}^{7} a_i\, x^i, \quad x = kR$$

| $i$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|-----|---|---|---|---|---|---|---|---|
| $a_i$ | $-j$ | $-1/2$ | $j/6$ | $1/24$ | $-j/120$ | $-1/720$ | $j/5040$ | $1/40320$ |

系数即 $a_i = (-j)^{i+1}/(i+1)!$。$R = 0$ 时自动得 $-jk/(4\pi)$，无分支、无相消。用 Horner 法求值。

> **借用注意**：老师版代码中 `Green(R) = e^{-jkR}/R` **不含** $1/4\pi$，其 `F1` 也不含。本文档 $G$ 含 $1/4\pi$，借用该多项式时必须整体除以 $4\pi$（上式已含）。单位约定不一致是两套代码对接时最易引入的隐性错误。

**8 阶截断的适用范围**：近场分支内 $R \leq d_{\text{th}} + h_M + h_N$，典型网格下 $kR \lesssim 2$，截断误差 $\sim (kR)^8/9! \lesssim 10^{-3}$ 相对量级。若将来调大 $d_{\text{th}}$，需重新核算或增加项数。

---
## 4. 解析 1/R 势：边界积分公式

对平面三角形 $T$（顶点 $\mathbf{v}_1,\mathbf{v}_2,\mathbf{v}_3$，单位法向 $\hat{\mathbf{n}}$），用边界积分把面积分化为三条边的求和。本节公式对**任意场点**（面内/面外/边上）成立，同时服务于自作用（$w_0=0$）、相邻对与近奇异对（$w_0 \neq 0$）。

### 4.1 标量势 $I_0$

$$I_0(\mathbf{r}) = \int_T \frac{1}{|\mathbf{r} - \mathbf{r}'|} dS' = -\sum_{i=1}^{3} h_i \ln\left(\frac{R_i^+ + s_i^+}{R_i^- + s_i^-}\right) - |d| \cdot \Omega$$

| 符号 | 定义 | 含义 |
|------|------|------|
| $\hat{\mathbf{t}}_i$ | $(\mathbf{v}_i^+ - \mathbf{v}_i^-)/L_i$ | 边单位切向量 |
| $\hat{\mathbf{m}}_i$ | $\hat{\mathbf{n}} \times \hat{\mathbf{t}}_i$（按 $\hat{\mathbf{m}}_i \cdot (\mathbf{v}_{\text{opp}} - \mathbf{v}_i^-) > 0$ 定向） | 边内法向 |
| $h_i$ | $(\mathbf{r} - \mathbf{v}_i^-) \cdot \hat{\mathbf{m}}_i$ | 场点到边所在直线的有符号面内距离 |
| $s_i^-$ / $s_i^+$ | $(\mathbf{r} - \mathbf{v}_i^{\mp}) \cdot \hat{\mathbf{t}}_i$ | 场点沿边切向到起/终点的有符号投影 |
| $R_i^-$ / $R_i^+$ | $\|\mathbf{r} - \mathbf{v}_i^{\mp}\|$ | 场点到起/终点的 3D 距离 |
| $d$ | $(\mathbf{r} - \mathbf{v}_1) \cdot \hat{\mathbf{n}}$ | 场点到平面的有符号距离 |
| $\Omega$ | 见 §4.3 | 三角形对场点张的立体角 |

公式前的负号保证 $\mathbf{r}$ 在三角形内部时 $I_0 > 0$；面内点 $d = 0$ 时立体角项自动归零。

### 4.2 矢量势与线性势

矢量势：

$$\mathbf{I}_{\text{vec}}(\mathbf{r}) = \int_T \frac{\mathbf{r}'}{R} dS' = \mathbf{p}\, I_0 + \mathbf{I}_{\text{shifted}}, \quad \mathbf{p} = \mathbf{r} - d\,\hat{\mathbf{n}}$$

平移矢量势（边界积分）：

$$\mathbf{I}_{\text{shifted}} = \frac12 \sum_{e=1}^{3} \hat{\mathbf{m}}_e^{\text{out}} \left[ \rho_e^2 \ln\left(\frac{R_{\text{start}} + s_{\text{start}}}{R_{\text{end}} + s_{\text{end}}}\right) + s_{\text{start}} R_{\text{start}} - s_{\text{end}} R_{\text{end}} \right]$$

- $\hat{\mathbf{m}}_e^{\text{out}} = -\hat{\mathbf{m}}_e$（外法向）；
- $\rho_e^2 = d^2 + d_e^2$，$d_e = (\mathbf{p} - \mathbf{v}_{\text{start}})\cdot\hat{\mathbf{m}}_e^{\text{out}}$；
- $s_{\text{start}}, s_{\text{end}}$ 为投影点 $\mathbf{p}$ 沿边切向到起/终点的有符号投影，$R_{\text{start}}, R_{\text{end}}$ 为场点到起/终点的 3D 距离；
- **注意对数比的方向与标量势相反**（start/end 而非 end/start），且 $sR$ 项同号配套，二者必须同时使用，不可只换其一。

线性势由矢量势解出，恒等式 $c_1 + c_2 + c_3 = I_0$ 由构造保证：

$$c_1 = \frac{[(\mathbf{I}_{\text{vec}} - \mathbf{v}_3 I_0) \times (\mathbf{v}_2 - \mathbf{v}_3)] \cdot \hat{\mathbf{n}}}{2A}, \quad c_2 = \frac{[(\mathbf{v}_1 - \mathbf{v}_3) \times (\mathbf{I}_{\text{vec}} - \mathbf{v}_3 I_0)] \cdot \hat{\mathbf{n}}}{2A}, \quad c_3 = I_0 - c_1 - c_2$$

> 不采用"每条边向两端点分配"的旧公式（不满足 $\Sigma c_i = I_0$，曾导致自作用 $P$ 项系统性错误）。矢量势法经测试与直接法一致。

### 4.3 立体角 $\Omega$

两种方式等价，任选其一：

**方式 A：Van Oosterom & Strackee 整体公式**

$$\Omega = 2\,\mathrm{atan2}\left( |\mathbf{r}_1 \cdot (\mathbf{r}_2 \times \mathbf{r}_3)|,\ r_1 r_2 r_3 + r_1(\mathbf{r}_2\cdot\mathbf{r}_3) + r_2(\mathbf{r}_3\cdot\mathbf{r}_1) + r_3(\mathbf{r}_1\cdot\mathbf{r}_2) \right)$$

其中 $\mathbf{r}_i = \mathbf{v}_i - \mathbf{r}$。面内点分子为 0，$\Omega$ 自动归零。

**方式 B：逐边 atan 求和（老师版实现）**

$$\Omega = \sum_{e=1}^{3} \beta_e, \quad \beta_e = \arctan\frac{P_{0,e}\, l_e^+}{R_{0,e}^2 + |d| R_e^+} - \arctan\frac{P_{0,e}\, l_e^-}{R_{0,e}^2 + |d| R_e^-}$$

其中 $P_{0,e}$ 为投影点到边的有符号面内距离，$R_{0,e}^2 = P_{0,e}^2 + d^2$。当 $R_{0,e}^2 \approx 0$（投影点落在边线上且 $d\approx0$）时置 $\beta_e = 0$，防止 $\arctan(0/0)$。

### 4.4 数值稳定化

利用 $R^2 - s^2 = \rho^2$，即 $R + s = \rho^2/(R - s)$，处理 $R + s \to 0$（场点靠近顶点）时的对数不稳定：

| 情形 | 标量势（$\rho^2 = d^2 + h_i^2$） | 矢量势（$\rho^2 = d^2 + d_e^2$） |
|------|------|------|
| 正常 | $h_i[\ln(R_e{+}s_e) - \ln(R_s{+}s_s)]$ | $\ln(R_s{+}s_s) - \ln(R_e{+}s_e)$ |
| 靠近起点（$R_s{+}s_s \to 0$） | $h_i[\ln(R_e{+}s_e) + \ln(R_s{-}s_s) - \ln\rho^2]$ | $\ln\rho^2 - \ln(R_s{-}s_s) - \ln(R_e{+}s_e)$ |
| 靠近终点（$R_e{+}s_e \to 0$） | $h_i[\ln\rho^2 - \ln(R_e{-}s_e) - \ln(R_s{+}s_s)]$ | $\ln(R_s{+}s_s) + \ln(R_e{-}s_e) - \ln\rho^2$ |

保护规则：$|h_i| < \varepsilon_{\log}$ 时该边标量势贡献取 0（$0 \cdot \ln\infty$ 型极限为 0）；$\rho^2$ 与 $R - s$ 取下限保护。$\varepsilon_{\log}$ 建议取 $10^{-12} \cdot h$。

---

## 5. 近奇异/奇异分支：统一阻抗计算

对进入近场分支的三角形对（场三角形 $T_N$ 承载场 RWG $N$，源三角形 $T_M$ 承载源 RWG $M$；自作用时 $T_M = T_N$），按"内层解析 + 外层数值 + 正则余项双重数值"计算。

### 5.1 矩的定义

**解析 1/R 部分（内层解析 × 外层 12 点高斯）：**

对每个场三角形高斯点 $\mathbf{r}_i \in T_N$，调用 §4 公式得 $I_0(\mathbf{r}_i) = \int_{T_M} 1/R\, dS'$ 与 $\mathbf{I}_{\text{vec}}(\mathbf{r}_i) = \int_{T_M} \mathbf{r}'/R\, dS'$（由 $c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + c_3\mathbf{v}_3$ 重建），然后：

| 矩 | 定义 | 加权坐标 |
|------|------|---------|
| $J_0$ | $A_N \sum_i w_i\, I_0(\mathbf{r}_i)$ | 无 |
| $\mathbf{J}_{\text{rio}}$ | $A_N \sum_i w_i\, \mathbf{r}_i\, I_0(\mathbf{r}_i)$ | **场坐标** |
| $\mathbf{J}_{\text{vec}}$ | $A_N \sum_i w_i\, \mathbf{I}_{\text{vec}}(\mathbf{r}_i)$ | **源坐标** |
| $J_{\text{rdot}}$ | $A_N \sum_i w_i\, \mathbf{r}_i \cdot \mathbf{I}_{\text{vec}}(\mathbf{r}_i)$ | 混合 |

**正则余项部分（双重 12×12 高斯，调用 $G_{\text{smooth}}$ 积分例程）：**

| 矩 | 定义 | 加权坐标 |
|------|------|---------|
| $I_1$ | $\iint G_{\text{smooth}}\, dS' dS$ | 无 |
| $\mathbf{I}_2$ | $\iint \mathbf{r}\, G_{\text{smooth}}\, dS' dS$ | **场坐标** |
| $\mathbf{I}_3$ | $\iint \mathbf{r}'\, G_{\text{smooth}}\, dS' dS$ | **源坐标** |
| $I_4$ | $\iint (\mathbf{r} \cdot \mathbf{r}')\, G_{\text{smooth}}\, dS' dS$ | 混合 |

### 5.2 配对规则（核心修正）

由 $\mathbf{f}_M(\mathbf{r}')\cdot\mathbf{f}_N(\mathbf{r}) = C_M C_N\, (\mathbf{r}' - \mathbf{r}_{\text{opp}}^M)\cdot(\mathbf{r} - \mathbf{r}_{\text{opp}}^N)$ 展开：

$$\mathbf{r}\cdot\mathbf{r}' \;-\; \mathbf{r}_{\text{opp}}^M \cdot \mathbf{r} \;-\; \mathbf{r}_{\text{opp}}^N \cdot \mathbf{r}' \;+\; \mathbf{r}_{\text{opp}}^M \cdot \mathbf{r}_{\text{opp}}^N$$

> **规则**：**源对顶点 $\mathbf{r}_{\text{opp}}^M$ 配场坐标矩（$\mathbf{J}_{\text{rio}}$ / $\mathbf{I}_2$）；场对顶点 $\mathbf{r}_{\text{opp}}^N$ 配源坐标矩（$\mathbf{J}_{\text{vec}}$ / $\mathbf{I}_3$）。**
>
> 原版 §5.2 中 $P_{\text{sing}}$ 的配对与此相反（源自该行展开式把 $\mathbf{f}_M$ 写在了场坐标上），仅因自作用对称性 $\mathbf{J}_{\text{rio}} \equiv \mathbf{J}_{\text{vec}}$ 才数值无害。本版统一后，同一规则对自作用、相邻、近奇异对全部成立——对后两者 $\mathbf{J}_{\text{rio}} \neq \mathbf{J}_{\text{vec}}$，配对写错就是真错误。

### 5.3 组装公式

$$P_{\text{sing}} = \frac{C_M C_N}{4\pi}\Big[ J_{\text{rdot}} - \mathbf{r}_{\text{opp}}^M\cdot\mathbf{J}_{\text{rio}} - \mathbf{r}_{\text{opp}}^N\cdot\mathbf{J}_{\text{vec}} + (\mathbf{r}_{\text{opp}}^M\cdot\mathbf{r}_{\text{opp}}^N)\, J_0 \Big]$$

$$Q_{\text{sing}} = -\frac{C_M C_N}{\pi}\, J_0$$

$$P_{\text{smooth}} = C_M C_N \Big[ (\mathbf{r}_{\text{opp}}^M\cdot\mathbf{r}_{\text{opp}}^N)\, I_1 - \mathbf{r}_{\text{opp}}^M\cdot\mathbf{I}_2 - \mathbf{r}_{\text{opp}}^N\cdot\mathbf{I}_3 + I_4 \Big]$$

$$Q_{\text{smooth}} = -4 C_M C_N\, I_1$$

$$Z_{\text{pair}} = jk\eta_0\,(P_{\text{sing}} + P_{\text{smooth}}) + \frac{j\eta_0}{k}\,(Q_{\text{sing}} + Q_{\text{smooth}})$$

### 5.4 自作用特例说明

$T_M = T_N$ 时无需特殊代码路径，以下性质自动成立，可作为校验：

- 所有场点在面内，$d = 0$（到舍入精度），立体角项为零；
- 由 $\mathbf{r} \leftrightarrow \mathbf{r}'$ 对称性：$\mathbf{J}_{\text{rio}} \equiv \mathbf{J}_{\text{vec}}$，$\mathbf{I}_2 \equiv \mathbf{I}_3$（解析恒等，数值上相差外层积分误差，应 $< 0.1\%$）；
- $G_{\text{smooth}}$ 双重积分中 12 对重合点 $R \equiv 0$，走 §3 的显式赋值/泰勒路径。

---
## 6. 远区数值积分

远场分支对完整 $G$ 做双重高斯积分，逐 RWG 三角形对累加：

$$Z_{mn} = \sum_{i,j=1}^{2} \left[ jk\eta_0\, P_{ij} + \frac{j\eta_0}{k}\, Q_{ij} \right]$$

四个格林积分（$G$ 完整形式，定义同 §5.1 下表）：

$$I_1 = \iint G\, dS' dS, \quad \mathbf{I}_2 = \iint \mathbf{r}\, G\, dS' dS, \quad \mathbf{I}_3 = \iint \mathbf{r}'\, G\, dS' dS, \quad I_4 = \iint (\mathbf{r}\cdot\mathbf{r}')\, G\, dS' dS$$

$$P_{ij} = C_M C_N \Big[ (\mathbf{r}_{\text{opp}}^M\cdot\mathbf{r}_{\text{opp}}^N)\, I_1 - \mathbf{r}_{\text{opp}}^M\cdot\mathbf{I}_2 - \mathbf{r}_{\text{opp}}^N\cdot\mathbf{I}_3 + I_4 \Big]$$

$$Q_{ij} = -4 C_M C_N\, I_1$$

配对规则与 §5.2 完全一致（源对顶点配场坐标矩，场对顶点配源坐标矩），全程序只此一种约定。

> **建议**：远场与近场分支的 $P/Q$ 组装共用同一段代码，仅替换"四矩来源"（完整 $G$ 数值 vs 解析+正则拆分），避免两套组装公式日后漂移不一致。

---

## 7. 数值积分方案

Dunavant 三角形对称高斯规则：

| 名称 | 点数 | 多项式精度 | 用途 |
|------|------|-----------|------|
| `GAUSS_3PT` | 3 | 2 | 简单测试 |
| `GAUSS_7PT` | 7 | 5 | 远场可选 |
| **`GAUSS_12PT`** | **12** | **6** | **近场分支默认（外层及 $G_{\text{smooth}}$ 双重积分）** |

双重积分：

$$\iint_{T_M \times T_N} f(\mathbf{r}, \mathbf{r}') dS' dS \approx A_M A_N \sum_{i,j} w_i w_j f(\mathbf{r}_i, \mathbf{r}'_j)$$

---

## 8. 阈值参数汇总

| 参数 | 建议值 | 作用 |
|------|--------|------|
| $d_{\text{th}}$ | $0.1\lambda$ | 远/近场分支切换（重心距判据） |
| $\varepsilon_R$ | $10^{-12}\cdot h$ | $G_{\text{smooth}}$ 小 $R$ 显式赋值阈值（方式 A） |
| $\varepsilon_{\log}$ | $10^{-12}\cdot h$ | §4.4 对数稳定化与 $|h_i|$ 判零阈值 |
| $\varepsilon_\beta$ | $(3\times10^{-2}\, L_e)^2$ | 方式 B 立体角中 $R_0^2$ 判零、跳过 $\beta_e$ |

---

## 9. 代码模块对应关系

| 模块文件 | 主要函数 | 功能 | 与原版差异 |
|---------|---------|------|-----------|
| `SINGULAR_INTEGRAL.F90`（建议改名，原版拼写 `INTEFRAL` 有误） | `ANALYTIC_SCALAR_POT_1_OVER_R` | 标量势边界积分（§4.1） | 不变 |
| | `ANALYTIC_LINEAR_POT_1_OVER_R` | 矢量势/线性势边界积分（§4.2） | 不变 |
| `GREEN_FUNCTIONS.F90` | `GREEN_FUNC_SMOOTH` | $G_{\text{smooth}}$，$R \to 0$ 显式赋值或泰勒求值（§3） | **修正 $R=0$ 返回值** |
| | `GREEN_FUNC` | 完整 $G$ | 不变 |
| | `CALC_GREEN_INTEGALS` | 完整 $G$ 四积分（§6） | 不变 |
| | `CALC_GREEN_SMOOTH_INTEGALS` | $G_{\text{smooth}}$ 四积分（§5.1） | 不变 |
| `Z_MATRIX.F90` | `COUNT_SHARED_VERTICES` | 共用顶点计数（替代原 `ARE_TRIANGLES_ADJACENT`） | **新增** |
| | `CALC_EFIE_NEAR_TRI_PAIR` | 近场分支统一阻抗（§5），泛化自原 `CALC_EFIE_SELF_TRI_PAIR` | **泛化** |
| | `CALC_EFIE_MATRIX_ELEMENT` | 统一入口，两分支（近场/远场） | **删除相邻分支** |
| `NUMERICAL_INTEGRATION.F90` | `INIT_GAUSS_TRI` / `GET_TRI_GLOBAL_GAUSS_POINTS` | 高斯数据 | 不变；`INIT_GAUSS_TRI(12)` 提到循环外只调一次 |
| `test_singular.F90` | `TEST_SINGULAR` | 验证程序（§11） | **扩充测试清单** |

---

## 10. 实现逻辑（算法步骤）

### 10.1 标量势 `ANALYTIC_SCALAR_POT_1_OVER_R`

输入：场点 $\mathbf{r}$，顶点 $\mathbf{v}_1,\mathbf{v}_2,\mathbf{v}_3$。输出：$I_0(\mathbf{r})$。

1. $\hat{\mathbf{n}} \leftarrow (\mathbf{v}_2-\mathbf{v}_1)\times(\mathbf{v}_3-\mathbf{v}_1)$ 归一化；$d \leftarrow (\mathbf{r}-\mathbf{v}_1)\cdot\hat{\mathbf{n}}$；$\text{RESULT} \leftarrow 0$；
2. 对每条边：算 $\hat{\mathbf{t}}$、$\hat{\mathbf{m}}$（按指向内部定向）、$h_i$、$s^{\pm}$、$R^{\pm}$；
3. 按 §4.4 表格分类计算该边贡献（$|h_i|<\varepsilon_{\log}$ 取 0；正常/靠近起点/靠近终点三分支）；累加；
4. 按 §4.3 算 $\Omega$（方式 A 或 B）；
5. $I_0 \leftarrow -\text{RESULT} - |d|\,\Omega$。

### 10.2 矢量势/线性势 `ANALYTIC_LINEAR_POT_1_OVER_R`

1. 调标量势得 $I_0$；$\mathbf{p} \leftarrow \mathbf{r} - d\hat{\mathbf{n}}$；
2. 对每条边：$\hat{\mathbf{m}}^{\text{out}}$、$s_{\text{start/end}}$、$R_{\text{start/end}}$、$d_e$、$\rho^2$，按 §4.4 稳定化算对数项，累加 $\frac12 \hat{\mathbf{m}}^{\text{out}}[\rho^2 \ln + s_{\text{start}}R_{\text{start}} - s_{\text{end}}R_{\text{end}}]$ 得 $\mathbf{I}_{\text{shifted}}$；
3. $\mathbf{I}_{\text{vec}} \leftarrow \mathbf{p}\, I_0 + \mathbf{I}_{\text{shifted}}$；
4. 按 §4.2 解出 $c_1, c_2, c_3$ 返回。

### 10.3 近场分支 `CALC_EFIE_NEAR_TRI_PAIR(MESH, TRI_M, TRI_N, RWG_M, RWG_N, K, ETA0, Z_PAIR)`

1. 提取 $C_M$、$\mathbf{r}_{\text{opp}}^M$（在 $T_M$ 上）与 $C_N$、$\mathbf{r}_{\text{opp}}^N$（在 $T_N$ 上）；
2. 外层 12 点高斯：对每个 $\mathbf{r}_i \in T_N$ 调 §10.1/§10.2（源三角形为 $T_M$），累加得 $J_0, \mathbf{J}_{\text{rio}}, \mathbf{J}_{\text{vec}}, J_{\text{rdot}}$（乘 $A_N$ 与权重）；
3. 按 §5.3 组装 $P_{\text{sing}}, Q_{\text{sing}}$（注意 §5.2 配对规则）；
4. 调 `CALC_GREEN_SMOOTH_INTEGALS`($T_N$, $T_M$, 12pt) 得 $I_1, \mathbf{I}_2, \mathbf{I}_3, I_4$；
5. 组装 $P_{\text{smooth}}, Q_{\text{smooth}}$，合成 `Z_PAIR`。

### 10.4 统一入口 `CALC_EFIE_MATRIX_ELEMENT`

```
对 RWG_M 的每个三角形 i，RWG_N 的每个三角形 j：
    N_share ← COUNT_SHARED_VERTICES(TRI_M(i), TRI_N(j))
    若 N_share ≥ 1 或 重心距 < d_th：
        Z_MN += CALC_EFIE_NEAR_TRI_PAIR(...)        ! 近场分支
    否则：
        普通高斯 + CALC_GREEN_INTEGALS，按 §6 组装   ! 远场分支
```

---

## 11. 验证标准与测试清单

**解析公式单元测试（沿用并重新解释原版）：**

| # | 测试 | 参考/判据 | 通过标准 |
|---|------|----------|---------|
| 1 | $I_0$ 面内/面外点 | 超高阶参考（Duffy 变换 + 高阶高斯，或 $\geq 500\times500$ 加密）；附远点渐近校验 $I_0 \to A/R$ | 相对误差 $< 0.1\%$ |
| 2 | $G_{\text{smooth}}$ 恒等式 $= G - 1/(4\pi R)$；$G_{\text{smooth}}(0) = -jk/(4\pi)$ | 解析 | $< 10^{-8}$；$R=0$ 处精确 |
| 3 | $\Sigma c_i = I_0$（含面外点） | 恒等式 | 机器精度 |
| 4 | $\mathbf{I}_{\text{vec}}$ | 同测试 1 参考 | $< 0.1\%$ |

> **参考解误差提醒**：原版以 $200\times200$ 乘积高斯作"真值"，其对奇异核自身误差即约 $0.1\%\sim0.5\%$，故原报告的 0.34% 偏差主要反映的是参考解误差而非解析公式误差（解析公式是精确的）。参考解阶数应远高于被测对象，或用 Duffy 变换消除奇异性后再加密。

**分支正确性测试（新增，本次修订的核心验证）：**

| # | 测试 | 方法 | 通过标准 |
|---|------|------|---------|
| 5 | 自作用对称性 $\mathbf{J}_{\text{rio}} = \mathbf{J}_{\text{vec}}$、$\mathbf{I}_2 = \mathbf{I}_3$ | 同一自作用对两种矩对比 | 相对差 $< 0.1\%$（外层积分误差量级） |
| 6 | 共边相邻对 $I_1$ | 与 Duffy 变换 + 高阶参考对比 | $< 0.1\%$ |
| 7 | 近奇异对（$d_{\min} = 0.01\lambda$ 平行三角形） | 与超高阶自适应数值对比 | $< 0.1\%$ |
| 8 | **阈值连续性**：构造若干对，使其 $d_c$ 横跨 $d_{\text{th}}$，同一对分别强制走两分支计算 | 两分支结果对比 | 相对差 $< 0.1\%$ |
| 9 | 共顶点对 | 与高阶参考对比 | $< 0.5\%$ |

> 测试 8 的意义：近场分支（解析+正则拆分）与远场分支（完整 $G$ 数值）在数学上恒等，阈值处结果必须连续。若不连续，说明至少一边有错——它能同时检验两个分支而无需外部真值。

**物理回归测试：**

| # | 测试 | 通过标准 |
|---|------|---------|
| 10 | 金属球 RCS 与 Mie 级数对比 | 谐振区外偏差 $< 1\%$ |
| 11 | 阻抗矩阵对角元虚部为负（电容性标量势主导） | 全对角满足 |
| 12 | 原版 11 项测试全部回归 | 全过 |

---

## 12. 与原版的差异汇总

| # | 位置 | 原版问题 | 本版处理 |
|---|------|---------|---------|
| 1 | §2 分类 | 拓扑分类：共边对蛮力高斯、共顶点对未分类、近距离对无处理 | 距离驱动 + 共顶点兜底，统一近场分支 |
| 2 | §5.2 配对 | $\mathbf{r}_{\text{opp}}$ 配对在 $P_{\text{sing}}$ 与 $P_{\text{smooth}}$/远场公式间相反，靠自作用对称性掩盖 | 统一为"源对顶点配场坐标矩，场对顶点配源坐标矩" |
| 3 | $Q$ 符号 | §1 用 $+Q_{\text{raw}}$ 配 $-j\eta_0/k$，§5 用 $-Q_{\text{raw}}$ 配 $+j\eta_0/k$，中途变号未声明 | §1 一次性定义 $Q := -\hat{Q}$，全文不变号 |
| 4 | §3 $G_{\text{smooth}}$ | $R < 10^{-10}$ 返回 0，丢弃对角虚部 | 显式赋值 $-jk/(4\pi)$，或泰勒求值无分支 |
| 5 | Bug #9 归因 | 原版称缺立体角修正"导致自作用 $P_{\text{sing}}$ 出错"——自作用场点恒在面内，该项恒零，不可能影响自阻抗；实际只影响面外点测试 | 已更正；且本版近场分支使 $w_0 \neq 0$ 能力真正投入使用 |
| 6 | 参考解解释 | "解析 vs 200×200 数值差 0.34%"被记为解析公式误差 | 更正为参考解自身误差主导，升级参考解方案（§11 注） |
| 7 | 模块 | 相邻分支冗余、自作用子程序仅限同一 RWG 逻辑 | 删除相邻分支，泛化为近场对子程序；`INIT_GAUSS_TRI(12)` 移出循环 |
| 8 | 拼写 | `SINGULAR_INTEFRAL`、`GARD_GREEN_FUNC` | 建议改名 `SINGULAR_INTEGRAL`、`GRAD_GREEN_FUNC`（非强制） |

---

## 参考文献

1. Wilton, D. R., et al. "Potential integrals for uniform and linear source distributions on polygonal and polyhedral domains." *IEEE Trans. Antennas Propag.*, 1984.
2. Rao, S. M., Wilton, D. R., and Glisson, A. W. "Electromagnetic scattering by surfaces of arbitrary shape." *IEEE Trans. Antennas Propag.*, 1982.
3. Graglia, R. D. "On the numerical integration of the linear shape functions times the 3-D Green's function or its gradient on a plane triangle." *IEEE Trans. Antennas Propag.*, 1993.
4. Dunavant, D. A. "High degree efficient symmetrical Gaussian quadrature rules for the triangle." *Int. J. Numer. Methods Eng.*, 1985.
5. Van Oosterom, A., and Strackee, J. "The solid angle of a plane triangle." *IEEE Trans. Biomed. Eng.*, 1983.
6. 老师版求解器奇异/近奇异处理文档（`Comp_eltsI` / `Comp_eltsU` / `Spara_cal` 架构，距离阈值 $0.1\lambda$ 方案来源）。
