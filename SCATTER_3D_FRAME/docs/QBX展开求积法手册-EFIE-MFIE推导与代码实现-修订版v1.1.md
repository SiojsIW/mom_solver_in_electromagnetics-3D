# QBX（展开求积法）手册：EFIE/MFIE 推导与代码实现（修订版 v1.1）

> 面向 RWG-Galerkin 矩量法（MoM）的 QBX 奇异性/近奇异性处理完整参考。
> 内容：数学推导（不跳步）→ 代码编写步骤 → 验证清单 → 特殊函数附录。
> 所有核心公式均经过数值实验验证（机器精度）。

> **修订版说明（v1.1，2026-08）**：本版依据独立数值审查（22 项实测，含 4 组随机配置的加法定理、极点取值/梯度、MFIE 旋度结构、双侧跳跃关系、特殊函数闭式等）修订，改动如下：
> 1. **【重要】§4.1 / §5.2 / §6.4 / §6.7：内中心坐标架。** 原版 §6.7 对内、外中心复用同一个局部架，导致内中心的场点实际落在局部架"南极"却套用北极公式，MFIE 内侧求值直接出错（实测大圆盘均匀电流：$L^-$ 算成 $+0.055$，正确值 $-0.5$；跳跃检验当场抓获）。修订为：**每个中心各自建架**，$\hat z'$ 始终取"中心指向场点"方向（内中心 $\hat z'=+\hat n$，$\hat e_2'$ 随之反向）。独立建架后 $L^\pm=\pm0.5$、跳跃 $=1.0$ 精确到机器精度。
> 2. **§8 Q7：符号与标度更正。** $h_\ell^{(2)}$ 小宗量渐近的虚部符号笔误（$-i\to+i$，与附录 B.4 一致）；乘积标度更正为 $j_\ell(kh)h_\ell^{(2)}(k\rho')\sim i\,(h/\rho')^\ell/[k\rho'(2\ell+1)]$；并补充实测确认的动态范围机制。
> 3. **§1.5：术语矛盾消除。** $\nabla G\sim1/R^2$ 由"超奇异（主值不可积）"更正为"强奇异（Cauchy 主值存在）"，与 §1.6/§2.3 统一；真正 $1/R^3$ 超奇异的是不散度转移的整体算子。
> 4. **§3.4：因果更正。** 大 $h$ 改善的是系数积分光滑性；收敛比 $q=h/\rho'_{\min}$ 随 $h$ 单调增大，场级数收敛反而变慢——不再写成"$h$ 大→收敛快"。
> 5. **§2.4：开销量级更正。** 单场点极点方案的系数开销由 $O(p^2)$ 更正为 $O(p)$（$m\in\{0,\pm1\}$ 共 $3(p+1)$ 个系数；$O(p^2)$ 是 panel-based 全 $m$ 方案）。
> 6. **小瑕疵 4 处**：§3.5 最短距离说法更正（自项最近点在内部投影点）；§3.5/§6.5/§6.7 许可性裕量方向统一为 $(1+\epsilon)$；参考文献 4 期刊更正（*SIAM J. Sci. Comput.* 40(3):A1225–A1249，原误作 *J. Comput. Phys.*）；§3.3 换序严格性加注（自项切点仅条件收敛，严格性靠 EGK 2013）。
> 7. **§3.2：补充收敛率定量估计**（自项场级数项 3D $\sim1/(k\ell^{3/2})$ 绝对收敛、2D $\sim\sin\ell/\ell$ 条件收敛，实测确认）。

---

## 第 0 章 符号、约定与总览

### 0.1 时间约定（先把坑堵住）

本手册**全程采用工程约定** $e^{j\omega t}$（与 FEKO、绝大多数 MoM 代码、你的求解器一致）：

| 量 | 工程约定（本手册） | 物理约定（Klöckner 等数学文献） |
|---|---|---|
| 时间因子 | $e^{j\omega t}$ | $e^{-i\omega t}$ |
| 外向格林函数 | $G=\dfrac{e^{-jkR}}{4\pi R}$ | $G=\dfrac{e^{ikR}}{4\pi R}$ |
| 球 Hankel | $h_\ell^{(2)}$ | $h_\ell^{(1)}$ |
| 2D Hankel | $H_\ell^{(2)}$ | $H_\ell^{(1)}$ |
| 加法定理系数 | $-jk$ | $+ik$ |

两套约定互为复共轭（实 $k$ 下）。**读数学圈论文时，把所有 $h^{(1)}$ 换成 $h^{(2)}$、$ik$ 换成 $-jk$，公式就平移过来了。** 附录 E 有完整对照表。

### 0.2 符号表

| 符号 | 含义 |
|---|---|
| $S$ | PEC 散射体表面，$\hat n$ 为指向外部区域的单位法向 |
| $k=\omega\sqrt{\mu\varepsilon}$，$\eta=\sqrt{\mu/\varepsilon}$ | 波数、波阻抗 |
| $G(\mathbf r,\mathbf r')=e^{-jkR}/(4\pi R)$，$R=\lvert\mathbf r-\mathbf r'\rvert$ | 3D 自由空间格林函数 |
| $\mathbf r$ | 场点（任意空间点） |
| $\mathbf r_0$ | **面上的**场点（积分方程的观测点、Galerkin 测试点） |
| $\mathbf r'$ | 源点（积分虚拟变量，在源三角形上跑） |
| $\mathbf r_c=\mathbf r_0+h\hat n$ | QBX 展开中心（离面距离 $h$） |
| $\boldsymbol\rho=\mathbf r-\mathbf r_c$，$\boldsymbol\rho'=\mathbf r'-\mathbf r_c$ | 相对中心的向量 |
| $p$ | QBX 截断阶（$\ell=0,\dots,p$） |
| $j_\ell,y_\ell,h_\ell^{(2)}=j_\ell-i y_\ell$ | 球 Bessel / 球 Neumann / 第二类球 Hankel |
| $Y_\ell^m$ | 归一化球谐函数，$\int Y_\ell^m\bar Y_{\ell'}^{m'}\,d\Omega=\delta_{\ell\ell'}\delta_{mm'}$，含 Condon–Shortley 相位 |
| $\mathbf f_n$ | 第 $n$ 个 RWG 基函数 |
| $T_n^\pm$ | RWG 基函数的正/负三角形 |

### 0.3 QBX 一句话

把"在奇异核旁边硬积分"换成"在离面中心处做局部展开：系数靠**光滑积分**算，场点值靠**收敛级数**算"。奇异性被"推远"了——核函数 $h_\ell^{(2)}(k\rho')$ 在 $\rho'\ge h>0$ 处是无限光滑的普通函数。

### 0.4 全局流程图（先在脑子里建这张地图）

```
Maxwell 方程组
   │  时谐 + 均匀无源区
   ▼
矢量 Helmholtz 方程 ──点源──▶ 格林函数 G ──叠加──▶ 位函数 A, φ
   │                                                    │
   │                              边界条件（PEC 表面）    │
   ▼                                                    ▼
EFIE（电场积分方程）  /  MFIE（磁场积分方程）  ──组合──▶ CFIE
   │
   │  RWG 离散 + Galerkin 测试
   ▼
矩阵元 = 外层测试积分 × 内层源积分  ←── 自项/近项奇异在这里
   │
   │  内层源积分 ← QBX（本手册核心）
   ▼
线性方程组 ZI = V
```

---

## 第 1 章 从 Maxwell 方程到 EFIE / MFIE

这一章把积分方程的来龙去脉补全。推导不跳步，但目标是"够用到能看懂 QBX 作用在哪个量上"，更细的泛函分析从略。

### 1.1 时谐场与散射问题设定

时谐约定 $e^{j\omega t}$ 下，Maxwell 旋度方程为

$$\nabla\times\mathbf E=-j\omega\mu\mathbf H,\qquad \nabla\times\mathbf H=j\omega\varepsilon\mathbf E+\mathbf J.$$

散射问题：入射场 $(\mathbf E^i,\mathbf H^i)$ 照射 PEC 体 $S$，在 $S$ 上感应出面电流 $\mathbf J$，感应电流在均匀无界空间（$\mu,\varepsilon$）中辐射散射场 $(\mathbf E^s,\mathbf H^s)$。外部区域**无源**（$\mathbf J$ 只活在面 $S$ 上），总场 $=$ 入射场 $+$ 散射场。

### 1.2 位函数与 Helmholtz 方程

外部无源区 $\nabla\cdot\mathbf B=0$，故可引入磁矢位 $\mathbf B^s=\nabla\times\mathbf A$。把 $\mathbf E^s=-j\omega\mathbf A-\nabla\phi$ 代入 Ampère 定律并取 Lorenz 规范 $\nabla\cdot\mathbf A=-j\omega\mu\varepsilon\,\phi$，得到两个**非齐次**标量 Helmholtz 方程：

$$\nabla^2\mathbf A+k^2\mathbf A=-\mu\mathbf J,\qquad \nabla^2\phi+k^2\phi=-\rho_s/\varepsilon,$$

其中 $\rho_s$ 为面电荷密度，由连续性方程（对 $e^{j\omega t}$）：

$$\nabla_{\!s}\cdot\mathbf J=-j\omega\rho_s\quad\Longrightarrow\quad \rho_s=\frac{j}{\omega}\nabla_{\!s}\cdot\mathbf J.$$

### 1.3 格林函数：点源响应

定义标量格林函数为点源响应：

$$\nabla^2 G+k^2G=-\delta(\mathbf r-\mathbf r'),\qquad G=\frac{e^{-jkR}}{4\pi R}.$$

**验证（不跳步）：** 对 $R>0$，$G=g(R)$ 只是 $R$ 的函数，球坐标径向拉氏算子给出

$$\frac{1}{R^2}\frac{d}{dR}\!\Big(R^2\frac{dg}{dR}\Big)+k^2g=0
\;\Longrightarrow\; \frac{d^2(Rg)}{dR^2}+k^2(Rg)=0,$$

故 $Rg=Ae^{-jkR}+Be^{+jkR}$。$e^{-jkR}\cdot e^{j\omega t}=e^{j(\omega t-kR)}$ 是**外向行波**（等相面 $R=\omega t/k$ 随时间外扩），$e^{+jkR}$ 是内向波，由辐射条件（能量只能向外流）弃之，得 $G=Ae^{-jkR}/R$。定 $A$：以 $\mathbf r'$ 为心作小球 $B_\epsilon$，对 $\nabla^2G+k^2G=-\delta$ 积分，

$$\int_{B_\epsilon}\nabla^2G\,dV=\oint_{\partial B_\epsilon}\frac{\partial G}{\partial R}\,dS
=Ae^{-jk\epsilon}\Big(-\frac{jk}{\epsilon}-\frac1{\epsilon^2}\Big)\cdot4\pi\epsilon^2
\xrightarrow{\epsilon\to0}-4\pi A,$$

而 $\int k^2G\,dV\to0$（$G\sim1/R$ 可积），$-\int\delta\,dV=-1$，故 $A=1/(4\pi)$。$\blacksquare$

### 1.4 叠加原理 → 位函数的积分表示

非齐次 Helmholtz 方程是线性的，面电流 $\mathbf J$ 是连续分布的点源 $\mathbf J(\mathbf r')\,dS'$，故

$$\mathbf A(\mathbf r)=\mu\int_S G(\mathbf r,\mathbf r')\mathbf J(\mathbf r')\,dS',\qquad
\phi(\mathbf r)=\frac1\varepsilon\int_S G(\mathbf r,\mathbf r')\rho_s(\mathbf r')\,dS'.$$

### 1.5 EFIE

散射电场

$$\mathbf E^s=-j\omega\mathbf A-\nabla\phi
=-j\omega\mu\int_S G\,\mathbf J\,dS'-\frac{j}{\omega\varepsilon}\nabla\int_S G\,(\nabla'_{\!s}\cdot\mathbf J)\,dS'.$$

用 $\omega\mu=k\eta$、$1/(\omega\varepsilon)=\eta/k$ 化简：

$$\boxed{\;\mathbf E^s(\mathbf r)=-jk\eta\int_S\Big[\,G\,\mathbf J+\frac1{k^2}\big(\nabla G\big)\big(\nabla'_{\!s}\cdot\mathbf J\big)\Big]dS'\;}$$

PEC 边界条件：面上切向总电场为零，$\hat n\times(\mathbf E^i+\mathbf E^s)=\mathbf 0$。于是

$$\text{EFIE:}\qquad jk\eta\,\hat n\times\int_S\Big[G\,\mathbf J+\frac1{k^2}\nabla G\,(\nabla'_{\!s}\cdot\mathbf J)\Big]dS'=\hat n\times\mathbf E^i,\quad \mathbf r_0\in S.$$

- 第一项核 $G\sim1/R$：**弱奇异**（面积元 $R\,dR$ 可积）；
- 第二项核 $\nabla G\sim1/R^2$：**强奇异**（Cauchy 主值意义下**存在**——对 RWG 的分片常数散度，领头奇项的角向均值为零，与 §1.6、§2.3 的定性一致；只是数值上不如弱奇异好处理），标准做法仍是 Galerkin 测试时用分部积分把一个 $\nabla$ 转移到测试函数上（见 §2.2），剩下 $1/R$ 弱奇异。注意：真正 $1/R^3$ 超奇异、连主值也不存在的，是不做散度转移时的整体算子 $\nabla\!\displaystyle\int G\,(\nabla'_{\!s}\cdot\mathbf J)\,dS'$（见 §2.3 表格）。【v1.1 修订：原版此处误作"超奇异（主值不可积）"，与 §1.6/§2.3 自相矛盾】

### 1.6 MFIE 与跳跃关系

散射磁场

$$\mathbf H^s=\frac1\mu\nabla\times\mathbf A=\nabla\times\int_S G\,\mathbf J\,dS'
=\int_S \nabla G\times\mathbf J\,dS',$$

最后一步用 $\nabla\times(G\mathbf J)=\nabla G\times\mathbf J+G(\nabla\times\mathbf J)$，而 $\mathbf J$ 是 $\mathbf r'$ 的函数、$\nabla$ 对 $\mathbf r$ 作用，故 $\nabla\times\mathbf J=\mathbf 0$。

定义

$$K(\mathbf J)(\mathbf r)=\int_S\nabla G(\mathbf r,\mathbf r')\times\mathbf J(\mathbf r')\,dS'.$$

核 $\nabla G\sim1/R^2$：**强奇异**，只以 Cauchy 主值（PV）存在。电流边界条件 $\mathbf J=\hat n\times(\mathbf H^i+\mathbf H^s)\big|_{S^+}$ 中，$\hat n\times K$ 从两侧逼近面时有跳跃：

$$\hat n\times K\Big|_{S^\pm}=\pm\frac12\mathbf J+\hat n\times\mathrm{PV}\,K(\mathbf J).$$

**用无限大平面验证符号（不跳步）：** 面上均匀电流 $\mathbf J=J_0\hat x$，$\hat n=\hat z$。由对称性，上半空间 $\mathbf H^s=\tfrac12\mathbf J\times\hat n=-\tfrac{J_0}2\hat y$（下半空间反号）。平面的主值积分为零（左右对称抵消），故 $\hat n\times\mathbf H^s\big|_+=\hat z\times(-\tfrac{J_0}2\hat y)=+\tfrac12J_0\hat x=+\tfrac12\mathbf J$。$\blacksquare$

代入边界条件（取 $+$ 侧，即外部）：

$$\text{MFIE:}\qquad \frac12\mathbf J-\hat n\times\mathrm{PV}\int_S\nabla G\times\mathbf J\,dS'=\hat n\times\mathbf H^i,\quad \mathbf r_0\in S.$$

### 1.7 CFIE

EFIE 与 MFIE 在 PEC 内谐振频率处各自有非唯一解（内谐振问题），线性组合可消除：

$$\text{CFIE}=\alpha\cdot\text{EFIE}+(1-\alpha)\,\eta\cdot\text{MFIE},\qquad 0<\alpha\le1\ (\text{典型 }0.2\sim0.5).$$

$$\alpha\, jk\eta\,\hat n\times\!\!\int\!\Big[G\mathbf J+\tfrac1{k^2}\nabla G(\nabla'\!\cdot\mathbf J)\Big]dS'
+(1-\alpha)\eta\Big[\tfrac12\mathbf J-\hat n\times\mathrm{PV}\!\!\int\!\nabla G\times\mathbf J\,dS'\Big]
=\alpha\,\hat n\times\mathbf E^i+(1-\alpha)\eta\,\hat n\times\mathbf H^i.$$

**对 QBX 的含义：** CFIE 需要的全部内层积分只有三类——$\int G\mathbf f\,dS'$（弱奇异）、$\int G\,(\nabla'\cdot\mathbf f)\,dS'$（弱奇异）、$\int\nabla G\times\mathbf f\,dS'$（强奇异、主值）。QBX 用**同一套展开机制**统一处理这三类，这正是它相对"奇异提取 + Duffy 各管一段"的结构性优势。

---

## 第 2 章 RWG-Galerkin 离散与奇异积分问题

### 2.1 RWG 基函数

第 $n$ 条内边共享三角形对 $T_n^+,T_n^-$（面积 $A_n^\pm$，边长 $\ell_n$），RWG 基函数

$$\mathbf f_n(\mathbf r)=\begin{cases}\dfrac{\ell_n}{2A_n^+}\big(\mathbf r-\mathbf v_n^+\big),&\mathbf r\in T_n^+,\\[6pt] \dfrac{\ell_n}{2A_n^-}\big(\mathbf v_n^--\mathbf r\big),&\mathbf r\in T_n^-,\\[6pt] \mathbf 0,&\text{其他}.\end{cases}$$

关键性质：**面散度在每个三角形上是常数**

$$\nabla_{\!s}\cdot\mathbf f_n=\begin{cases}+\ell_n/A_n^+,&T_n^+,\\ -\ell_n/A_n^-,&T_n^-.\end{cases}$$

（由 $\nabla_{\!s}\cdot(\mathbf r-\mathbf v)=2$ 立得。）电流展开 $\mathbf J\approx\sum_n I_n\mathbf f_n$。

### 2.2 Galerkin 测试与矩阵元

用 $\mathbf f_m$ 测试 EFIE（$\mathbf f_m$ 本身切向，故直接点乘即可），标量位那一项用面积分形式的分部积分：

$$\int_{T_m}\mathbf f_m\cdot\nabla\psi\,dS=-\int_{T_m}(\nabla_{\!s}\cdot\mathbf f_m)\psi\,dS+\oint_{\partial T_m}\psi\,\mathbf f_m\cdot\hat\nu\,dl,$$

边界项在基函数的两个三角形间逐边抵消（$\mathbf f_m$ 的法向分量跨边连续、外法向反向），故 EFIE 矩阵元为

$$Z_{mn}^E=\underbrace{jk\eta\int_{T_m}\mathbf f_m\cdot\overbrace{\int_{T_n}G\,\mathbf f_n\,dS'}^{\text{矢量内层}}dS}_{\text{矢量位部分}}
-\underbrace{\frac{j\eta}{k}\int_{T_m}(\nabla_{\!s}\cdot\mathbf f_m)\overbrace{\int_{T_n}G\,(\nabla'_{\!s}\cdot\mathbf f_n)\,dS'}^{\text{标量内层}}dS}_{\text{标量位部分}}.$$

**要点：散度转移之后，EFIE-RWG 里只剩 $1/R$ 弱奇异内层积分。** 超奇异被"分摊"成了两个弱奇异积分的乘积结构。QBX 对 EFIE 只需要求**内层积分的值**，不需要梯度。

MFIE 矩阵元（用 $\mathbf f_m\cdot(\hat n\times\mathbf v)=\mathbf v\cdot(\mathbf f_m\times\hat n)$ 改写测试项）：

$$Z_{mn}^M=\frac12\int_{T_m}\mathbf f_m\cdot\mathbf f_n\,dS-\int_{T_m}\big(\mathbf f_m\times\hat n\big)\cdot\mathrm{PV}\underbrace{\int_{T_n}\nabla G\times\mathbf f_n\,dS'}_{=K(\mathbf f_n)}\,dS.$$

第一项是普通质量矩阵元（无奇异，常规高斯积分）。第二项的内层是 $1/R^2$ 强奇异主值积分——**这是 MFIE 的硬骨头，也是 QBX 梯度公式（§5）的用武之地**。

### 2.3 奇异性分类汇总

| 内层积分 | 核的渐近 | 奇异级别 | 出现位置 | 经典处理 |
|---|---|---|---|---|
| $\int G\mathbf f\,dS'$ | $1/R$ | 弱奇异（可积） | EFIE 矢量位 | 奇异提取、Duffy |
| $\int G\,(\nabla'\cdot\mathbf f)\,dS'$ | $1/R$ | 弱奇异 | EFIE 标量位 | 同上 |
| $\int\nabla G\times\mathbf f\,dS'$ | $1/R^2$ | 强奇异（主值） | MFIE | 主值扣除 + 解析项 |
| $\nabla\int G\,(\nabla'\cdot\mathbf f)\,dS'$ | $1/R^3$ 结构 | 超奇异 | EFIE（不散度转移时）、Nyström | 一般回避 |

**近奇异**（场点离源单元很近但不重合）比自项更难缠：不奇异但剧烈变化，固定阶高斯积分精度崩塌，且没有解析结构可利用。QBX 对自项、近项用**同一套代码路径**处理——这是第二个结构性优势。

### 2.4 现有工程方法（对比基线，写论文 related work 用）

1. **奇异项提取（singularity extraction/subtraction）**：把 $G$ 拆成 $1/(4\pi R)+[G-1/(4\pi R)]$，静态 $1/R$ 部分用解析公式（Ylä-Oijala & Taskinen, IEEE TAP 2003 是 CFIE 矩阵元的标准参考文献），余项数值积分。便宜、成熟；但每个核、每个算子都要单独推导解析式，阶数固定在低阶。
2. **Duffy 变换**：把三角形映射到消去 $1/R$ 雅可比的参数域，配合高斯点。通用；但对近奇异需要自适应细分，阶数仍受限。
3. **DIRECTFN 类（Polimeridis 等）**：解析内层 + 数值外层，高精度但逐对推导。
4. **QBX（本手册）**：一种机制通吃所有核（$\nabla$ 多少次都行——对局部展开逐项微分即可），阶数 $p$ 可调，自项近项统一；代价是实现复杂度与每对近场交互的系数开销——单场点极点方案只需 $m\in\{0,\pm1\}$ 共 $3(p+1)$ 个系数，即 $O(p)$（panel-based 全 $m$ 方案才是 $O(p^2)$，见 §4.5 实现要点与 §6.8）。【v1.1 修订：原版误写 $O(p^2)$】

---

## 第 3 章 QBX 原理

### 3.1 核心恒等式：3D Helmholtz 加法定理

设 $\mathbf r_c$ 为展开中心，$\boldsymbol\rho=\mathbf r-\mathbf r_c$，$\boldsymbol\rho'=\mathbf r'-\mathbf r_c$。当 $|\boldsymbol\rho|<|\boldsymbol\rho'|$ 时（完整推导见附录 D）：

$$\boxed{\;\frac{e^{-jk|\mathbf r-\mathbf r'|}}{4\pi|\mathbf r-\mathbf r'|}
=-jk\sum_{\ell=0}^{\infty}\sum_{m=-\ell}^{\ell}
j_\ell(k\rho)\,Y_\ell^m(\hat\rho)\;\cdot\;
h_\ell^{(2)}(k\rho')\,\bar Y_\ell^m(\hat\rho')\;}$$

**结构解读（这是整个 QBX 的命门）：** 右边把"场点函数" $j_\ell(k\rho)Y_\ell^m(\hat\rho)$ 与"源点函数" $h_\ell^{(2)}(k\rho')\bar Y_\ell^m(\hat\rho')$ **完全分离**。场点侧的 $j_\ell Y_\ell^m$ 处处光滑（$j_\ell(k\rho)\sim\rho^\ell$，附录 B）；源点侧的 $h_\ell^{(2)}$ 只在 $\rho'=0$ 奇异，只要源点离中心足够远就是普通光滑函数。

### 3.2 几何：切球与收敛半径

面上场点 $\mathbf r_0$，法向 $\hat n$，取中心 $\mathbf r_c=\mathbf r_0+h\hat n$（$h>0$ 为参数）。以 $\mathbf r_c$ 为心、$h$ 为半径作球 $B_h$：

- $B_h$ 与面 $S$ 在 $\mathbf r_0$ 处**相切**（中心在法线上，半径恰好到面）；
- **场点在球内/球面上**：$\mathbf r_0$ 本身在球面上（$|\boldsymbol\rho|=h$）；
- **源点在球外/球面上**：对局部平坦的源三角形，$|\boldsymbol\rho'|^2=h^2+d^2\ge h^2$（$d$ 为源点到 $\mathbf r_0$ 的面内距离），等号仅在切点成立。

级数对每一对 $(\mathbf r,\mathbf r')$ 要求 $|\boldsymbol\rho|<|\boldsymbol\rho'|$。最坏情形是切点处 $|\boldsymbol\rho|=|\boldsymbol\rho'|=h$，比值趋于 1——收敛变慢但**仍然收敛**（Epstein–Greengard–Klöckner, SIAM J. Numer. Anal. 2013 对光滑面上的 Laplace 核给出了切点收敛的严格证明；Helmholtz 核的奇异性结构相同，工程实践一致支持）。**收敛半径 $=\mathrm{dist}(\mathbf r_c,S)=h$。**

> **收敛速度的定量感受（v1.1 补充，实测确认）：** 自项场级数的项渐近 $\sim 1/(k\,\ell^{3/2})$（3D，**绝对收敛**）；2D 情形 $\sim\sin\ell/\ell$（**条件收敛**，对数慢）。两者都是"收敛"，没有过度承诺——只是切点附近想压到很低误差需要相应提高 $p$。

> **误区澄清 1：** "球和源点相切"≠"球里只能装一个源点"。球不是用来"装源点"的——恰恰相反，规则是**场点在球内、源点在球外**。切点 $\mathbf r_0$ 是球与面的唯一公共点（局部），其余所有源点都在球外，满足 $|\boldsymbol\rho'|\ge h$。

### 3.3 两个"代入"的区别（最容易绕晕的地方）

对内层积分 $u(\mathbf r)=\displaystyle\int_T G(\mathbf r,\mathbf r')\sigma(\mathbf r')\,dS'$ 应用加法定理：

$$u(\mathbf r)=\sum_{\ell=0}^\infty\sum_{m=-\ell}^\ell
\underbrace{\Big[-jk\int_T h_\ell^{(2)}(k\rho')\,\bar Y_\ell^m(\hat\rho')\,\sigma(\mathbf r')\,dS'\Big]}_{\text{系数 }\alpha_{\ell m}\ \leftarrow\ \textbf{第一次代入：源点 }}
\cdot\;\underbrace{j_\ell(k\rho)Y_\ell^m(\hat\rho)}_{\textbf{第二次代入：场点}}.$$

- **第一次代入（源点 $\mathbf r'\to\boldsymbol\rho'$）：发生在积分号内**，$h_\ell^{(2)}(k\rho')\bar Y_\ell^m$ 对 $\rho'\ge h>0$ 是光滑普通函数，**没有任何收敛性问题**，系数积分用（加密）高斯积分就能算准。积分与求和换序的合法性由级数在 $|\boldsymbol\rho'|\ge h+\epsilon$ 上的一致收敛保证。【v1.1 加注：严格说，一致收敛论证对 $\min\rho'>h$ 的**邻居单元**成立；**自项**的积分域含切点 $\rho'=h$，级数在那里只有条件收敛，换序与取值的严格性依赖 EGK 2013 对切点的分析（§3.2 已引）——工程上由第 7 章测试 2–4 兜底。】
- **第二次代入（场点 $\mathbf r\to\boldsymbol\rho$）：发生在级数上**，这才是"收敛半径"约束作用的地方：$|\boldsymbol\rho|$ 必须小于（最坏等于）所有 $|\boldsymbol\rho'|$。

> **误区澄清 2：** "源三角形上有的高斯点到中心距离大于收敛半径，是不是违规？"——不违规。$|\boldsymbol\rho'|$ **越大越好**（$h_\ell^{(2)}$ 越大越平滑、比值 $|\boldsymbol\rho|/|\boldsymbol\rho'|$ 越小收敛越快）。要担心的是 $|\boldsymbol\rho'|<h$ 的源点（折叠的邻居三角形、薄结构对面），那才破坏 $|\boldsymbol\rho|<|\boldsymbol\rho'|$。

### 3.4 截断误差与参数 $p$、$h$

截断到 $\ell\le p$ 的尾部受比值 $q=\max|\boldsymbol\rho|/\min|\boldsymbol\rho'|$ 控制。远离切点时 $q<1$，误差 $\sim O(q^{p+1}/(1-q))$，几何收敛；切点附近 $q\to1$，退化为代数收敛。工程经验：

| 参数 | 典型取值 | 说明 |
|---|---|---|
| $h$ | $(0.5\sim1.5)\times$ 单元尺寸 $h_T$ | $h$ 大→源点离中心更远→**系数积分更光滑**、动态范围更小（源侧 $h_\ell^{(2)}(k\rho')$ 变化更缓，加密高斯点数可少）；但收敛比 $q=\lvert\boldsymbol\rho\rvert/\lvert\boldsymbol\rho'\rvert=h/\sqrt{h^2+d^2}$（$d$ 为源点面内距离）随 $h$ **单调增大**（$h\to\infty$ 时 $q\to1$），**场级数收敛变慢、所需 $p$ 略增**——两者权衡，故取 $0.5\sim1.5$ 倍单元尺寸；且球大容易罩住别的几何（薄结构）【v1.1 修订：原版"$h$ 大→收敛快"因果倒置】 |
| $p$ | $4\sim10$ | 2–4 位有效数字的常用区间；$p$ 每加 1 成本约线性增 |
| 源积分加密 | 16–37 点 Dunavant 规则 | 系数积分要在单元上分辨 $h_\ell^{(2)}(k\rho')\bar Y_\ell^m$ 的角向振荡（$\sim\ell$ 个波长），故需比常规装配更高阶 |

**精度天花板提醒：** QBX 自身可以做到很高阶，但 RWG（分片线性电流、分片常数电荷）+ Galerkin 的离散误差是 $O(h_T)\sim O(h_T^2)$。把 $p$ 堆到 15 不会让你的 RCS 更准——积分精度只需不拖离散误差的后腿。论文里不要宣称"谱精度"，要宣称"统一、稳健、参数可控"。

### 3.5 许可性检查（admissibility）

对每个（场点，源单元）对，在算系数**之前**检查：

$$\min_{\mathbf r'\in T}\;|\mathbf r'-\mathbf r_c|\ \overset{?}{\ge}\ h\,(1+\epsilon),\qquad \epsilon\sim10^{-3}\text{ 的相对裕量}.$$

实际实现：严格做法是求 $\mathbf r_c$ 到源三角形的**最短距离**——最近点可能在三角形**内部**（自项情形的最近点恰是内部投影点 $\mathbf r_0$，$\min\rho'=h$），不一定在顶点【v1.1 修订：原版"凸组合下最小值出现在顶点附近"不成立】。简便起见可检查三顶点与加密高斯点到 $\mathbf r_c$ 的距离作保守筛查，但对自项/近共面单元应显式计算点-三角形最短距离。

违规情形与对策：

1. **折叠邻居**（二面角内凹，邻居三角形折进球内）→ 该对退回 Duffy/奇异提取，或对该场点缩小 $h$；
2. **薄结构**（对面壁距离 $<2h$）→ 同上；
3. **非流形边、T 型接头** → 网格预处理时标记，该对走经典路径。

**QBX 不是全网格无脑替换，而是"近场交互的可插拔高精度引擎"：通过检查的近场对走 QBX，违规对走经典路径，远场对走普通低阶高斯。** 这个混合策略也是你论文里"稳健性"故事的支点。

---

## 第 4 章 EFIE 的 QBX 推导

### 4.1 局部坐标架与"极点"约定

对每个面上场点 $\mathbf r_0$（即测试三角形上的一个高斯点）：

1. 中心 $\mathbf r_c=\mathbf r_0+h\hat n$（外中心；MFIE 另需内中心 $\mathbf r_0-h\hat n$，见 §5.2）；
2. 建局部正交架 $\{\hat e_1',\hat e_2',\hat z'\}$：**取 $\hat z'$ 沿"中心指向场点"方向**，即

$$\hat z'=\frac{\mathbf r_0-\mathbf r_c}{h}.$$

对外中心 $\hat z'=-\hat n$；对内中心 $\hat z'=+\hat n$——**每个中心各自建架**【v1.1 修订强调】。$\hat e_1'$ 取切平面内任一单位向量（如测试三角形第一条边方向归一化），$\hat e_2'=\hat z'\times\hat e_1'$（内中心的 $\hat e_2'$ 随之反向）。

于是场点在中心坐标系中的位置是 $\boldsymbol\rho_0=\mathbf r_0-\mathbf r_c=h\hat z'$——**恰好在半径 $h$ 球面的北极**（$\theta'=0$）。下面所有求值公式都利用"北极"这一特殊位置，这是单场点 QBX 实现的关键简化来源。

> **【v1.1 警示】** "北极"是全部求值/梯度公式（§4.2、§4.3、§4.5、§5.1）的**隐含前提**：场点必须位于自己局部架的 $\theta'=0$ 处。若对内中心复用外中心的架（$\hat z'=-\hat n$ 不变），则 $\mathbf r_0-\mathbf r_c^-=+h\hat n=-h\hat z'$，场点实际在**南极**（$\theta'=\pi$）却套用北极公式，结果直接错误（实测大圆盘均匀电流：内侧极限 $L^-$ 算出 $+0.055$，正确值 $-0.5$）。§5.2 的跳跃检验能零成本抓获此类错误——这也是它被列为首选调试手段的原因。

### 4.2 标量内层（EFIE 标量位部分）

目标：对源三角形 $T$ 上的密度 $\sigma(\mathbf r')=\nabla'_{\!s}\cdot\mathbf f_n$（**常数**，见 §2.1），算

$$S(\mathbf r_0)=\int_T G(\mathbf r_0,\mathbf r')\,\sigma(\mathbf r')\,dS'.$$

**第一步（写系数）。** 加法定理代入并交换积分求和：

$$S(\mathbf r)=\sum_{\ell=0}^\infty\sum_{m=-\ell}^\ell \alpha_{\ell m}\,j_\ell(k\rho)Y_\ell^m(\hat\rho),\qquad
\alpha_{\ell m}=-jk\int_T h_\ell^{(2)}(k\rho')\,\bar Y_\ell^m(\hat\rho')\,\sigma(\mathbf r')\,dS'.$$

**第二步（数值算系数）。** 系数积分光滑，用源三角形上的加密高斯规则（节点 $\mathbf r'_q$、权 $w_q$、雅可比已含在 $w_q$）：

$$\alpha_{\ell m}\approx-jk\,\sigma\sum_q w_q\,h_\ell^{(2)}(k\rho'_q)\,\bar Y_\ell^m(\hat\rho'_q),
\qquad \boldsymbol\rho'_q=\mathbf r'_q-\mathbf r_c,$$

其中 $(\rho'_q,\theta'_q,\varphi'_q)$ 是 $\boldsymbol\rho'_q$ 在**局部架**里的球坐标：

$$\cos\theta'_q=\frac{\boldsymbol\rho'_q\cdot\hat z'}{\rho'_q},\qquad
\varphi'_q=\mathrm{atan2}\big(\boldsymbol\rho'_q\cdot\hat e_2',\ \boldsymbol\rho'_q\cdot\hat e_1'\big).$$

**第三步（极点求值）。** 在北极 $\theta'=0$ 处，$Y_\ell^m(0,\varphi')=0\ (m\neq0)$，且

$$Y_\ell^0(0,\varphi')=\sqrt{\frac{2\ell+1}{4\pi}}$$

（因为 $Y_\ell^m\propto P_\ell^m(\cos\theta')e^{im\varphi'}$，而 $P_\ell^m(1)=0$ 对 $m\ge1$——$P_\ell^m(x)=(-1)^m(1-x^2)^{m/2}\frac{d^mP_\ell}{dx^m}$ 含因子 $(1-x^2)^{m/2}$；$P_\ell(1)=1$。）

$$\boxed{\;S(\mathbf r_0)=\sum_{\ell=0}^{p}\alpha_{\ell 0}\;j_\ell(kh)\sqrt{\frac{2\ell+1}{4\pi}}\;}$$

**只有 $m=0$ 的系数参与求值。**（已数值验证至机器精度。）

### 4.3 矢量内层（EFIE 矢量位部分）

目标：$\mathbf V(\mathbf r_0)=\displaystyle\int_T G(\mathbf r_0,\mathbf r')\,\mathbf f_n(\mathbf r')\,dS'$。

逐分量照搬 §4.2，系数变成**矢量系数**：

$$\boldsymbol\beta_{\ell m}=-jk\int_T h_\ell^{(2)}(k\rho')\,\bar Y_\ell^m(\hat\rho')\,\mathbf f_n(\mathbf r')\,dS'
\approx-jk\sum_q w_q\,h_\ell^{(2)}(k\rho'_q)\bar Y_\ell^m(\hat\rho'_q)\,\mathbf f_n(\mathbf r'_q),$$

$$\mathbf V(\mathbf r_0)=\sum_{\ell=0}^p\boldsymbol\beta_{\ell 0}\,j_\ell(kh)\sqrt{\frac{2\ell+1}{4\pi}}.$$

### 4.4 组装进 EFIE 矩阵元

对近场对 $(T_m,T_n)$，外层测试积分仍用常规高斯规则（节点 $\mathbf r_{q'}$、权 $w_{q'}$），但内层值由 QBX 供给：

$$Z_{mn}^{E}\approx jk\eta\sum_{q'}w_{q'}\,\mathbf f_m(\mathbf r_{q'})\cdot\mathbf V_n(\mathbf r_{q'})
-\frac{j\eta}{k}\sum_{q'}w_{q'}\,(\nabla_{\!s}\cdot\mathbf f_m)\,S_n(\mathbf r_{q'}),$$

其中每个 $\mathbf r_{q'}$ 都有自己的中心 $\mathbf r_c(\mathbf r_{q'})$ 和自己的一组系数（§7.6 装配主循环）。

> **组织结构铁律（"一场点一中心，一源单元一组系数"）：** 系数 $\alpha_{\ell m},\boldsymbol\beta_{\ell m}$ 依赖于（中心，源单元）这一对组合。装配矩阵元 $Z_{mn}$ 时，只能用"当前场点的中心 $+$ 当前源单元 $T_n$"算出的系数；**把别的单元的系数混进来就是错误**（不是浪费，是交叉污染）。只有在求解后的总势求值（matvec、近场计算）中，才可以把同一中心下多个源单元的系数累加成一套再求值——那时合并是合法的，因为目标本来就是 $\int_S G\sigma\,dS'$ 全体源的总和。

### 4.5 需要梯度的时候（可选路径）

如果走"不散度转移"的路线（Nyström、或未来做高阶元、或对标 Peng–Lee 的直接超奇异求值），需要 $\nabla S$。局部展开逐项可微（$j_\ell Y_\ell^m$ 光滑，微分与求和可交换）：

$$\nabla S(\mathbf r_0)=\sum_{\ell,m}\alpha_{\ell m}\,\nabla\big[j_\ell(k\rho)Y_\ell^m(\hat\rho)\big]\Big|_{\text{极点}}.$$

极点处球坐标梯度 $\nabla u=\hat z'\dfrac{\partial u}{\partial\rho}+\dfrac1\rho\Big(\hat\theta'\dfrac{\partial u}{\partial\theta'}+\dfrac{\hat\varphi'}{\sin\theta'}\dfrac{\partial u}{\partial\varphi'}\Big)$。逐项分析（完整推导不跳步）：

**(a) 径向分量。** $\dfrac{\partial}{\partial\rho}$ 只动 $j_\ell(k\rho)$，极点处仍只有 $m=0$：

$$\frac{\partial S}{\partial\rho}\Big|_{\text{极点}}=\sum_{\ell=0}^p \alpha_{\ell0}\,k\,j_\ell'(kh)\sqrt{\frac{2\ell+1}{4\pi}},
\qquad j_\ell'(x)=j_{\ell-1}(x)-\frac{\ell+1}{x}j_\ell(x).$$

**(b) 切向分量：只有 $m=\pm1$ 存活。** 对 $m=1$：在 $\theta'\to0$ 附近，$P_\ell^1(\cos\theta')=-\sin\theta'\,P_\ell'(\cos\theta')\approx-\theta'\,P_\ell'(1)$。由 Legendre 方程 $(1-x^2)P_\ell''-2xP_\ell'+\ell(\ell+1)P_\ell=0$ 在 $x\to1$ 取极限：$-2P_\ell'(1)+\ell(\ell+1)P_\ell(1)=0$，即 $P_\ell'(1)=\ell(\ell+1)/2$，故 $P_\ell^1(\cos\theta')\approx-\theta'\ell(\ell+1)/2$。于是（$N_{\ell m}$ 为 $Y$ 的归一化常数）

$$\Big(\hat\theta'\partial_{\theta'}+\tfrac{\hat\varphi'}{\sin\theta'}\partial_{\varphi'}\Big)\Big[N_{\ell1}P_\ell^1(\cos\theta')e^{i\varphi'}\Big]\Big|_{\theta'\to0}
=-\frac{\ell(\ell+1)}2\,N_{\ell1}\,e^{i\varphi'}\big(\hat\theta'+i\hat\varphi'\big).$$

极点处 $\varphi'$ 本来看似无定义，但组合 $e^{i\varphi'}(\hat\theta'+i\hat\varphi')=\hat e_1'+i\hat e_2'$ 是**常矢量**（把 $\hat\theta'=\cos\varphi'\hat e_1'+\sin\varphi'\hat e_2'$、$\hat\varphi'=-\sin\varphi'\hat e_1'+\cos\varphi'\hat e_2'$ 代入即得），奇异性完全消掉。$m=-1$ 同理（用 $P_\ell^{-1}=-P_\ell^1/[\ell(\ell+1)]$，得 $e^{-i\varphi'}(\hat\theta'-i\hat\varphi')=\hat e_1'-i\hat e_2'$）。注意到 $N_{\ell,1}\cdot\frac{\ell(\ell+1)}2=N_{\ell,-1}\cdot\frac12=c_\ell$，其中

$$c_\ell=\frac12\sqrt{\frac{(2\ell+1)\,\ell(\ell+1)}{4\pi}},$$

最终

$$\boxed{\;\nabla_t S\Big|_{\text{极点}}=\frac1h\sum_{\ell=1}^p j_\ell(kh)\,c_\ell\,
\Big[\big(\alpha_{\ell,-1}-\alpha_{\ell,1}\big)\hat e_1'
-i\big(\alpha_{\ell,1}+\alpha_{\ell,-1}\big)\hat e_2'\Big]\;}$$

$$\boxed{\;\nabla S(\mathbf r_0)=\nabla_t S+\hat z'\frac{\partial S}{\partial\rho}
=\nabla_t S-\hat n\,\frac{\partial S}{\partial\rho}\;}$$

（$\hat z'=-\hat n$。$\ell=0$ 项切向贡献自动为零，因为 $c_0=0$。已数值验证至机器精度。）

**实现要点：单场点 QBX 只需 $m\in\{0,\pm1\}$ 三类系数**——每个 $\ell$ 算 3 个系数而非 $2\ell+1$ 个，系数总数 $3(p+1)$ 而非 $(p+1)^2$。这是把 QBX 塞进工程代码时最重要的成本节约。（若做"一中心覆盖多个场点"的 panel-based QBX，才需要全部 $m$。）

---

## 第 5 章 MFIE 的 QBX 推导

### 5.1 旋度落到局部展开上

MFIE 内层 $K(\mathbf f_n)(\mathbf r)=\displaystyle\int_T\nabla G\times\mathbf f_n\,dS'=\nabla\times\int_T G\,\mathbf f_n\,dS'=\nabla\times\mathbf V(\mathbf r)$。

代入 §4.3 的矢量局部展开，逐项取旋度。系数 $\boldsymbol\beta_{\ell m}$ 是**常矢量**，用恒等式 $\nabla\times(\psi\mathbf b)=\nabla\psi\times\mathbf b$（$\mathbf b$ 常矢量）：

$$K(\mathbf f_n)(\mathbf r)=\sum_{\ell,m}\nabla\big[j_\ell(k\rho)Y_\ell^m(\hat\rho)\big]\times\boldsymbol\beta_{\ell m}.$$

**梯度结构和 EFIE 标量位完全相同**——§4.5 的极点梯度公式（m=0 法向 + m=±1 切向）原样适用，只是把标量系数 $\alpha_{\ell m}$ 换成矢量系数 $\boldsymbol\beta_{\ell m}$，最后再叉乘：

$$\frac{\partial}{\partial\rho}:\quad \mathbf g_r=\sum_{\ell=0}^p k\,j_\ell'(kh)\sqrt{\tfrac{2\ell+1}{4\pi}}\;\hat z'\times\boldsymbol\beta_{\ell0}$$

$$\text{切向}:\quad \mathbf g_t=\frac1h\sum_{\ell=1}^p j_\ell(kh)\,c_\ell\,
\Big[\hat e_1'\times\big(\boldsymbol\beta_{\ell,-1}-\boldsymbol\beta_{\ell,1}\big)
-i\,\hat e_2'\times\big(\boldsymbol\beta_{\ell,1}+\boldsymbol\beta_{\ell,-1}\big)\Big]$$

$$K(\mathbf f_n)(\mathbf r_0)=\mathbf g_t+\mathbf g_r.$$

**结论：MFIE 与 EFIE 矢量位共享同一套系数 $\boldsymbol\beta_{\ell m}$。** 一个源单元对一个场点，系数算一次，EFIE 矢量位（取值）和 MFIE（取值 + 叉乘结构）同时使用——CFIE 的成本大头被摊薄了一半，这是 CFIE-QBX 的第三个性结构优势。

### 5.2 强奇异的处理：双侧中心与主值

$K$ 的核 $\sim1/R^2$，主值积分的值依赖趋近方向。QBX 用**两个中心**复现这个结构：

$$\mathbf r_c^+=\mathbf r_0+h\hat n\ (\text{外侧}),\qquad \mathbf r_c^-=\mathbf r_0-h\hat n\ (\text{内侧}).$$

**【v1.1 修订：两个中心，两套局部架】** §4.5/§5.1 的极点公式是**北极版**（场点在 $\hat z'$ 方向、$\theta'=0$）。对外中心 $\mathbf r_0-\mathbf r_c^+=-h\hat n=h\hat z'$，前提成立；但对内中心 $\mathbf r_0-\mathbf r_c^-=+h\hat n$，场点落在"外中心架"的**南极**。因此每个中心必须各自建架（§4.1 的约定 $\hat z'=(\mathbf r_0-\mathbf r_c)/h$ 逐中心执行）：

$$\text{外中心：}\ \hat z'=\frac{\mathbf r_0-\mathbf r_c^+}{h}=-\hat n;\qquad
\text{内中心：}\ \hat z'=\frac{\mathbf r_0-\mathbf r_c^-}{h}=+\hat n,\quad \hat e_2'=\hat z'\times\hat e_1'\ \text{随之反向}.$$

若坚持两个中心共用一个架，则内中心必须改用**南极版公式**（级数值乘 $(-1)^\ell$、径向梯度方向取 $-\hat z'$、切向组合符号翻转）——能推，但独立建架是干净得多的修法，出错面更小。实测对比（大圆盘均匀电流，$h=0.5$，$p=60$）：共架时 $L^-$ 算成 $+0.055$（正确值 $-0.5$），跳跃检验直接崩溃；独立建架后 $L^+=+0.500000$、$L^-=-0.500000$、跳跃 $=1.000000$，精确到机器精度。

各算各的系数、各求各的值，得到两个单侧极限 $L^+,L^-$（$L^\pm=\hat n\times K$ 的对应侧极限）。然后：

$$\mathrm{PV}\,[\hat n\times K]=\frac{L^++L^-}{2},\qquad
\text{跳跃检验：}\ L^+-L^-\overset{?}{=}\mathbf J(\mathbf r_0)\ \text{（被测试的密度本身）}.$$

跳跃关系 $\hat n\times K|_{S^\pm}=\pm\frac12\mathbf J+\mathrm{PV}$（§1.6）保证两侧之差恰好是密度——**这是零成本的实现自检**：跳跃检验不通过，说明系数、球谐、符号、局部架至少有一处写错了。调试 MFIE-QBX 时先跑这个。

### 5.3 组装进 MFIE / CFIE 矩阵元

$$Z_{mn}^{M}\approx\frac12\int_{T_m}\mathbf f_m\cdot\mathbf f_n\,dS
-\sum_{q'}w_{q'}\,\big(\mathbf f_m(\mathbf r_{q'})\times\hat n\big)\cdot K_n^{\mathrm{PV}}(\mathbf r_{q'}),$$

其中 $K_n^{\mathrm{PV}}(\mathbf r_{q'})=\tfrac12\big[K_n^{+}(\mathbf r_{q'})+K_n^{-}(\mathbf r_{q'})\big]$。CFIE 矩阵元 $Z_{mn}=\alpha Z_{mn}^E+(1-\alpha)\eta Z_{mn}^M$。

**成本说明：** MFIE 需要双侧中心，系数开销是 EFIE 单侧的两倍；但 $\frac12\mathbf J$ 项无奇异，走普通路径。对外侧问题，若采用"只取外中心 + 解析补 $\frac12\mathbf J$"的混合方案（$K|_{S^+}=\mathrm{PV}+\frac12\mathbf J$ 的单侧形式），可以只用单侧中心——但那样就放弃了跳跃自检，调试期不建议。

---

## 第 6 章 代码编写步骤（Fortran 工程视角）

### 6.1 模块划分

```
qbx_params.f90      ! 参数: p, h_factor, 加密规则, 许可性裕量
qbx_special.f90     ! j_l, y_l, h_l^(2), P_l^m, Y_l^m (递推, 见 §6.3)
qbx_frame.f90       ! 局部坐标架构建
qbx_coef.f90        ! 系数计算 (单中心 × 单源单元)
qbx_eval.f90        ! 极点求值 / 极点梯度 / 叉乘(MFIE)
qbx_assembly.f90    ! 装配主循环改造 + 许可性检查 + 经典路径回退
qbx_test.f90        ! 跳跃检验、与奇异提取对拍、单元测试
```

### 6.2 参数块

```fortran
module qbx_params
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  integer,  parameter :: dp       = real64
  complex(dp), parameter :: jj    = (0.0_dp, 1.0_dp)   ! 虚数单位 j
  real(dp), parameter :: PI       = 3.14159265358979323846_dp

  integer,  parameter :: P_TRUNC  = 6        ! QBX 截断阶 (4~10)
  real(dp), parameter :: H_FACTOR = 1.0_dp   ! h = H_FACTOR * 测试单元尺寸
  integer,  parameter :: N_GAUSS_UP = 25     ! 源单元加密高斯点数 (Dunavant 25 点)
  real(dp), parameter :: ADM_TOL  = 1.0e-3_dp! 许可性相对裕量
  real(dp), parameter :: R_NEAR   = 3.0_dp   ! 近场判定: 距离 < R_NEAR*单元尺寸 走 QBX
end module qbx_params
```

### 6.3 特殊函数例程（稳定递推）

**球 Bessel/Neumann**（附录 B 的递推关系）：

```fortran
! 返回 j_l(x), y_l(x), l=0..p。x>0。
subroutine sph_bessel_all(x, p, jl, yl)
  ! j_l: x < p+1 时用降递推(Miller 算法)保稳定; x 较大时升递推即可
  ! y_l: 升递推无条件稳定
  !   j_0 = sin(x)/x,  j_1 = sin(x)/x^2 - cos(x)/x
  !   y_0 = -cos(x)/x, y_1 = -cos(x)/x^2 - sin(x)/x
  !   递推: f_{l+1} = (2l+1)/x * f_l - f_{l-1}
  ! h_l^(2) = j_l - i*y_l  (由调用方组合)
end subroutine
```

**连带 Legendre 与球谐**（附录 C）：只需 $m\in\{0,\pm1\}$ 三列。

```fortran
! 计算 Y_l^m(ct, phi), m in {0,+1,-1}, l=0..p; 同时给出共轭(系数用)
subroutine sph_harm_3m(ct, phi, p, Y0, Yp1, Ym1)
  ! P_l(ct):        升递推 (l+1)P_{l+1} = (2l+1)*ct*P_l - l*P_{l-1}
  ! P_l^1(ct):      P_l^1(x) = -(1-x^2)^{1/2} * P_l'(x),
  !                 或用稳定递推; 注意 Condon-Shortley 相位 (-1)^m 已含
  ! P_l^{-1}:       P_l^{-m} = (-1)^m (l-m)!/(l+m)! P_l^m
  ! Y_l^m = sqrt((2l+1)(l-m)!/(4*pi*(l+m)!)) * P_l^m(ct) * exp(i*m*phi)
  ! 阶乘比用连乘积算, 不要用 gamma(阶乘大会溢出)
end subroutine
```

**单元测试先行**：$j_\ell,y_\ell$ 对已知闭式值（附录 B）；$Y_\ell^m$ 对归一性 $\sum_m|Y_\ell^m|^2=(2\ell+1)/(4\pi)$（球谐加法定理的特例）。这一步不过，后面全白搭。

### 6.4 局部坐标架

```fortran
! 输入: "中心指向场点"的单位向量 to_field = (r0-rc)/h；输出: 局部架 e1, e2, zp (= to_field)
!   【v1.1 修订】输入语义从 n_hat 改为 to_field，使内外中心都能正确建架:
!     外中心: call build_frame(-n_hat, ...) ；内中心: call build_frame(+n_hat, ...)
!   两个中心必须各自调用、不可复用同一个架（否则内中心场点落在南极却套北极公式, §5.2）
subroutine build_frame(to_field, e1, e2, zp)
  real(dp), intent(in)  :: to_field(3)
  real(dp), intent(out) :: e1(3), e2(3), zp(3)
  zp = to_field
  ! 任取不与 zp 平行的向量叉乘出 e1 (如测试三角形的第一条边方向)
  e1 = ... ; e1 = e1 - dot_product(e1,zp)*zp ; e1 = e1/norm2(e1)
  e2 = cross(zp, e1)   ! 右手系
end subroutine
```

### 6.5 系数计算例程（一个中心 × 一个源单元）

```fortran
! 输入: 中心 rc, 局部架(e1,e2,zp), 波数 k, 源三角形 verts(3,3),
!       源密度数据 (sigma 常数; fn 在加密高斯点上的矢量值)
! 输出: alpha(0:p, -1:1) 标量系数, beta(0:p, -1:1, 3) 矢量系数
subroutine qbx_coefficients(rc, e1, e2, zp, k, verts, sigma, fn_q, alpha, beta)
  call get_up_gauss(verts, N_GAUSS_UP, rq, wq)   ! 加密高斯节点/权(含雅可比)
  alpha = (0.0_dp, 0.0_dp) ; beta = (0.0_dp, 0.0_dp)
  do q = 1, N_GAUSS_UP
     d   = rq(:,q) - rc
     rho = norm2(d)
     ct  = dot_product(d, zp)/rho
     phi = atan2( dot_product(d,e2), dot_product(d,e1) )
     call sph_bessel_all(k*rho, P_TRUNC, jl, yl)
     h2  = jl - jj*yl                                    ! h_l^(2)
     call sph_harm_3m(ct, phi, P_TRUNC, Y0, Yp1, Ym1)
     ! m = 0
     alpha(:,0) = alpha(:,0) - jj*k*wq(q)*sigma * h2(:)*conjg(Y0(:))
     beta(:,0,:) = beta(:,0,:) - jj*k*wq(q)*h2(:)*conjg(Y0(:)) * spread(fn_q(:,q),1,P_TRUNC+1)
     ! m = +1, -1 同理
     ...
  end do
end subroutine
```

**许可性检查放在循环前**：`if (min_dist(verts, rc) < (1+ADM_TOL)*h) → 置回退标志，该对走经典路径`。【v1.1 修订：裕量方向与 §3.5 的 $\min\rho'\ge h(1+\epsilon)$ 统一——原版写成 $(1-\text{ADM\_TOL})$，相当于允许 $\rho'$ 略小于 $h$，方向相反】

### 6.6 极点求值/梯度例程

```fortran
! 求值 (EFIE 用): 只需 m=0
subroutine qbx_value(alpha, kh, p, u)
  call sph_bessel_all(kh, p, jl, yl)
  u = sum( alpha(:,0) * jl(:) * sqrt((2*l+1)/(4*pi)) )   ! l=0..p 逐项
end subroutine

! 梯度 (MFIE / 直接超奇异路径用): m=0 法向 + m=±1 切向
subroutine qbx_gradient(alpha, kh, h, e1, e2, zp, grad)
  call sph_bessel_all(kh, p, jl, yl)
  jlp(0) = cos(kh)/kh - sin(kh)/kh**2          ! j_0' = -j_1, 或直接用递推
  do l = 1, p ; jlp(l) = jl(l-1) - (l+1)/kh*jl(l) ; end do
  dudrho = sum( alpha(:,0) * k*jlp(:) * sqrt((2*l+1)/(4*pi)) )
  grad_t = 0
  do l = 1, p
     cl = 0.5_dp*sqrt((2*l+1)*l*(l+1)/(4*pi))
     grad_t = grad_t + jl(l)*cl*( (alpha(l,-1)-alpha(l,1))*e1 &
             - jj*(alpha(l,1)+alpha(l,-1))*e2 )
  end do
  grad = grad_t/h + zp*dudrho
end subroutine
```

MFIE 叉乘：把上面 `alpha` 换成矢量 `beta`（逐分量算梯度矢量），再按 §5.1 做 $\hat z'\times\boldsymbol\beta_{\ell0}$、$\hat e_1'\times(\cdots)$ 的叉乘组合。也可以更直接：先对每个 $\ell,m$ 算出标量梯度因子 $g_{\ell m}$（$\nabla[j_\ell Y_\ell^m]$ 在极点的值，3 分量复数），然后 $\mathbf K=\sum g_{\ell m}\times\boldsymbol\beta_{\ell m}$——**预计算 $g_{\ell m}$ 表，MFIE 就只是查表叉乘**。

### 6.7 装配主循环改造

```fortran
do m = 1, N_edges                     ! 测试基函数
  do n = 1, N_edges                   ! 源基函数
    dist = triangle_pair_distance(Tm, Tn)
    if (dist > R_NEAR*h_T) then
       Z(m,n) = standard_gauss_assembly(m,n)          ! 远场: 原路径不动
    else
       Z_E = 0 ; Z_M = 0
       do qp = 1, n_gauss_test                        ! 外层测试高斯点 r0
          r0 = test_node(Tm, qp) ; n_hat = test_normal(Tm, qp)
          ok = .true.
          ! ---- 外中心 ----
          rc = r0 + H_FACTOR*h_T*n_hat
          call build_frame(-n_hat, e1, e2, zp)   ! 【v1.1】外中心自己的架: zp=(r0-rc)/h=-n_hat
          if (min_dist(Tn_verts, rc) < (1+ADM_TOL)*H_FACTOR*h_T) ok = .false.  ! 【v1.1】裕量方向同 §3.5
          if (ok) then
             call qbx_coefficients(rc, ..., Tn, sigma_n, fn_q, alpha, beta)
             call qbx_value(alpha, k*H_FACTOR*h_T, P_TRUNC, S_val)
             call qbx_value_vec(beta, ..., V_val)               ! EFIE 矢量位
             if (need_MFIE) call qbx_K(beta, ..., K_plus)       ! 查表叉乘
          end if
          ! ---- MFIE 需要内中心 ----
          if (ok .and. need_MFIE) then
             rc = r0 - H_FACTOR*h_T*n_hat
             call build_frame(+n_hat, e1, e2, zp)   ! 【v1.1 关键修复】内中心必须重建架: zp=(r0-rc)/h=+n_hat
             if (min_dist(Tn_verts, rc) < (1+ADM_TOL)*H_FACTOR*h_T) ok = .false.
             if (ok) then
                (系数 + K_minus，全部用内中心自己的这套架)
             end if
             K_PV = 0.5_dp*(K_plus + K_minus)
             ! 跳跃自检(调试模式): jump = n_hat×(K_plus - K_minus) ≈ f_n(r0)
             !   ↑ 坐标架用错（内外中心共架）时此检验当场崩溃，是它定位了 v1.1 的问题 1
          end if
          if (.not. ok) then
             (S_val, V_val, K_PV) = classic_path(r0, Tn)       ! Duffy/奇异提取
          end if
          Z_E = Z_E + w(qp)*( jj*k*eta*dot(fm, V_val) - jj*eta/k*divfm*S_val )
          Z_M = Z_M + w(qp)*( -dot(cross(fm,n_hat), K_PV) )
       end do
       Z(m,n) = alpha_cfie*Z_E + (1-alpha_cfie)*eta*(Z_M + 0.5*mass(m,n))
    end if
  end do
end do
```

### 6.8 性能注意

1. **系数复用层级**：同一（场点，源单元）的 $\boldsymbol\beta$ 被 EFIE 矢量位与 MFIE 共享（§5.1）；同一源单元对**不同**场点的系数不能复用（中心不同）。若内存/时间吃紧，可把同一测试三角形的所有外层高斯点共用一个中心（panel-based QBX），$h$ 取大到罩住整个测试三角形仍满足许可性——这时需要全部 $m$ 的系数（$(p+1)^2$ 个），是一次"系数数量 vs 中心数量"的权衡。
2. **远场不要动**：QBX 只替换近场对（通常占矩阵元的 5–15%），装配总开销增幅可控；它**不加速迭代求解**，加速的是填充阶段的精度/稳健性。
3. **宽频带接口（第二篇论文的种子）**：$h_\ell^{(2)}(x)$ 有有限闭式 $h_\ell^{(2)}(x)=\dfrac{e^{-ix}}{x}\times(\ell\text{ 次 }1/x\text{ 多项式})$（附录 B 给出 $\ell=0,1,2$ 显式）。系数积分里频率依赖全部锁在 $e^{-jk\rho'}$ 与 $k$ 的幂次中；对中心频点 $k_0$ 预存"矩" $\int_T \sigma(\mathbf r')\,\rho'^{-q}\,(\cdots)\,dS'$，新频点用 $e^{-jk\rho'}=e^{-jk_0\rho'}\sum_s\frac{(-j\Delta k\rho')^s}{s!}$ 泰勒平移——QBX 系数变成**频率无关矩表**的线性组合，近场填充可跨频率复用。近场带宽 $|\Delta k|\rho'_{\max}\sim$ 几个弧度，折算 $\Delta f/f_0$ 相当宽，近场不是瓶颈（瓶颈在远场/矩阵-向量乘）。

---

## 第 7 章 验证清单（按顺序做，每步不过就停）

| # | 测试 | 通过标准 |
|---|---|---|
| 1 | $j_\ell,y_\ell,Y_\ell^m$ 单元测试 | 闭式值/归一性到 $10^{-14}$ |
| 2 | **加法定理单点测试**：单点源 $\sigma\delta(\mathbf r'-\mathbf r_p)$（球外），QBX 级数值 vs 解析 $G$ | 相对误差 $<10^{-12}$（$p=30$） |
| 3 | **极点梯度单点测试**：同上，比 $\nabla G$ 解析值 | $<10^{-12}$ |
| 4 | **跳跃检验**：平板/球面上 RWG 密度，$L^+-L^-$ vs $\mathbf f_n(\mathbf r_0)$ | $<10^{-3}$（工程精度） |
| 5 | **与奇异提取对拍**：你现有 EFIE 代码的近场矩阵元 vs QBX 矩阵元 | 相对误差 $<10^{-4}\sim10^{-6}$（随 $p$ 改善） |
| 6 | 许可性检查触发测试：人为构造折叠网格/薄板 | 违规对正确回退经典路径 |
| 7 | **基准算例**：PEC 球 RCS vs Mie 级数（CFIE，含内谐振频点） | RCS 误差 $<1\%$；谐振频点无伪解 |
| 8 | $p$ 收敛曲线、$h$ 敏感度曲线 | 画出论文核心图：$p=4\sim10$ 平台区、$h/h_T=0.5\sim1.5$ 稳定区 |

测试 2、3 就是本手册写作前所做的数值验证（单点源、$k=1.7$、$h=1$、$p=40$），取值与梯度均达机器精度——公式本身是对的，你的实现若不过，问题在代码不在公式。

---

## 第 8 章 常见误区 FAQ（历史踩坑汇总）

**Q1：球和源三角形相切，那球里不是只能装一个源点吗？**
A：球不是用来装源点的。规则是场点在球内、源点在球外。切点是球与面的唯一公共点，其余源点都在球外——这正是加法定理要求的配置。

**Q2：源三角形上有的高斯点到中心距离超过收敛半径 $h$，违规吗？**
A：不违规，$|\boldsymbol\rho'|$ 越大越好（收敛越快）。要担心的是 $|\boldsymbol\rho'|<h$（折叠邻居、薄结构对面），许可性检查管的就是这个。

**Q3：系数到底用"场点所在三角形"的源点算，还是"附近所有三角形"的？**
A：看场景。**装配矩阵元 $Z_{mn}$**：中心属于当前场点，系数按源单元逐个算——$Z_{mn}$ 只用 $T_n$ 的源点算一套系数，不能混入别的单元（混入即错误）。**总势求值**（matvec、后处理近场）：同一中心下所有近场源单元的系数可以累加成一套再求值，合法且省时。

**Q4：代入场点时为什么只能代"这一个点"？换个场点为什么不对？**
A：级数 $u(\mathbf r)=\sum\alpha_{\ell m}j_\ell(k\rho)Y_\ell^m$ 只在球 $|\boldsymbol\rho|<\min|\boldsymbol\rho'|$ 内收敛。换的场点若落在球外，级数发散或收敛到错误值。所以一场点一中心（或 panel-based：一中心覆盖一个小面片，§6.8）。

**Q5：QBX 能加快求解速度吗？**
A：不能。QBX 替换的是**矩阵填充**阶段的近场积分；迭代求解阶段（matvec）若配合快速多极子等，QBX 系数结构可以被复用（Wala–Klöckner 2019 做的就是这个），但那属于快速算法论文的地盘。你的卖点是**统一性**（一种机制通吃所有算子所有奇异级别）、**稳健性**（许可性检查 + 回退）、**参数可控**（$p,h$ 规则），不是速度。

**Q6：有了奇异提取，为什么还要 QBX？**
A：奇异提取对每个核逐个推导解析式，阶数固定、近奇异仍要配自适应；QBX 一种机制处理所有核（任意阶微分都行）、自项近项同一路径、$p$ 可调。它是"下一代统一框架"的候选——也是 Peng–Lee（2016，IEDG 框架）和你（RWG-Galerkin 框架）先后看中它的原因。

**Q7：低频（$k\to0$）会出问题吗？**
A：$j_\ell(kh)\sim(kh)^\ell/(2\ell+1)!!$ 快速变小，$h_\ell^{(2)}(k\rho')\sim +i(2\ell-1)!!/(k\rho')^{\ell+1}$ 快速变大【v1.1 修订：原版虚部误作 $-i$。由 $h^{(2)}=j-iy$、$y_\ell\sim-(2\ell-1)!!/x^{\ell+1}$ 得 $-iy_\ell\sim+i(2\ell-1)!!/x^{\ell+1}$，与附录 B.4 一致，实测 $\ell=0..3$ 比值 $1.000000$】，乘积良性：

$$j_\ell(kh)\,h_\ell^{(2)}(k\rho')\;\sim\;\frac{i\,(h/\rho')^\ell}{k\rho'\,(2\ell+1)}\qquad(k\to0),$$

【v1.1 修订：原版乘积标度"$\rho'^{\,-\ell-1}\cdot(h/\rho')^\ell$"把 $\rho'$ 的幂次多算了一阶；上式由 $(2\ell+1)!!=(2\ell+1)(2\ell-1)!!$ 化简即得】但中间量动态范围大，$p\le10$、双精度下无忧。**实测确认机制（v1.1 补充）**：$h=0.5$ 时 $\ell=60$ 的系数中间量可达 $\sim10^{120}$，但对应场权重 $j_{60}(kh)\sim10^{-123}$，乘积良性；$\ell$ 稍大后系数中的双精度对消噪声（实测 $\ell=20$ 时 $m=\pm1$ 系数相对噪声 $\sim7\%$）完全被场权重压死，不影响结果。$k\to0$ 时 QBX 退化为 Laplace 核的 QBX，与静态奇异提取可比——测试 5 可以在 $k$ 很小时做，更灵敏。

---

## 参考文献

1. A. Klöckner, A. Barnett, L. Greengard, M. O'Neil, "Quadrature by expansion: A new method for the evaluation of layer potentials," *J. Comput. Phys.* 252:332–349, 2013.（QBX 原始论文，2D Helmholtz）
2. C. L. Epstein, L. Greengard, A. Klöckner, "On the convergence of local expansions of layer potentials," *SIAM J. Numer. Anal.* 51(5):2660–2679, 2013.（切点收敛的严格证明）
3. A. H. Barnett, "Evaluation of layer potentials close to the boundary for Laplace and Helmholtz problems on analytic planar domains," *SIAM J. Sci. Comput.* 36(2), 2014.
4. L. af Klinteberg, A.-K. Tornberg, "Adaptive Quadrature by Expansion for Layer Potential Evaluation in Two Dimensions," *SIAM J. Sci. Comput.* 40(3):A1225–A1249, 2018. DOI: 10.1137/17M1121615.（加密策略、panel-based QBX）【v1.1 修订：期刊名更正，原版误作 *J. Comput. Phys.*】
5. M. Wala, A. Klöckner, "A fast algorithm for Quadrature by Expansion, I: Conceptually," *J. Comput. Phys.* 388:655–689, 2019.（QBX-FMM，**速度类宣称的禁区**）
6. S. Peng, J.-F. Lee, "Direct Evaluation of Hyper-Singularity in Surface Integral Equation Using Quadrature by Expansion," IEEE AP-S/URSI, 2016；及 S. Peng 博士论文（Ohio State, 2019）.（**电磁圈最近先行工作**：并矢格林球谐展开 + QBX，IEDG 框架，必须引用并差异化）
7. P. Ylä-Oijala, M. Taskinen, "Calculation of CFIE impedance matrix elements with RWG and n̂×RWG functions," *IEEE Trans. Antennas Propag.* 51(8), 2003.（奇异提取基线）
8. S. M. Rao, D. R. Wilton, A. W. Glisson, "Electromagnetic scattering by surfaces of arbitrary shape," *IEEE Trans. Antennas Propag.* 30(3), 1982.（RWG）
9. F. Vico et al., "A New Quadrature Method for Singular Integrals of Boundary Integral Equations in Electromagnetism," IEEE AP-S/URSI, 2020.（Ewald 分裂 + 渐近展开，非 QBX，对比方法）
10. NIST Digital Library of Mathematical Functions (DLMF), Ch. 10, 18.（特殊函数公式权威来源）

---

## 附录 A 二维问题：Bessel 函数、Hankel 函数与 Graf 加法定理

### A.1 从 Helmholtz 方程到 Bessel 方程（分离变量，不跳步）

2D Helmholtz 方程 $\nabla^2\psi+k^2\psi=0$。极坐标 $(\rho,\varphi)$ 下拉氏算子

$$\nabla^2\psi=\frac{\partial^2\psi}{\partial\rho^2}+\frac1\rho\frac{\partial\psi}{\partial\rho}+\frac1{\rho^2}\frac{\partial^2\psi}{\partial\varphi^2}.$$

设 $\psi=R(\rho)\Phi(\varphi)$，代入并乘 $\rho^2/(R\Phi)$：

$$\frac{\rho^2R''+\rho R'}{R}+k^2\rho^2=-\frac{\Phi''}{\Phi}.$$

左边只含 $\rho$、右边只含 $\varphi$，两边必等于同一常数。记为 $\ell^2$：

- **角向**：$\Phi''=-\ell^2\Phi\Rightarrow\Phi=e^{i\ell\varphi}$。物理要求 $\psi$ 单值（$\varphi$ 绕一圈回原值）$\Rightarrow$ $\ell$ 必须是**整数**。
- **径向**：$\rho^2R''+\rho R'+(k^2\rho^2-\ell^2)R=0$。令 $x=k\rho$，得 **Bessel 方程**

$$x^2R''+xR'+(x^2-\ell^2)R=0.$$

### A.2 两个独立解与 Hankel 函数

Bessel 方程是二阶常微分方程，有两个线性独立解：

- **第一类 Bessel 函数 $J_\ell(x)$**：在原点正则，$J_\ell(x)\sim\dfrac{x^\ell}{2^\ell\ell!}\ (x\to0)$。代表"驻波"，可在包含原点的区域使用。
- **第二类 Bessel 函数（Neumann）$Y_\ell(x)$**：在原点奇异，$Y_0\sim\frac2\pi\ln x$，$Y_\ell\sim-\dfrac{2^\ell(\ell-1)!}{\pi x^\ell}$。

类比 $\sin,\cos\to e^{\pm ix}$，定义 **Hankel 函数**

$$H_\ell^{(1)}=J_\ell+iY_\ell\sim\sqrt{\frac{2}{\pi x}}\,e^{+i(x-\ell\pi/2-\pi/4)},\qquad
H_\ell^{(2)}=J_\ell-iY_\ell\sim\sqrt{\frac{2}{\pi x}}\,e^{-i(x-\ell\pi/2-\pi/4)}\quad(x\to\infty).$$

**为什么散射问题用 $H^{(2)}$（工程约定）？** 远场 $H_\ell^{(2)}(k\rho)e^{j\omega t}\sim e^{j(\omega t-k\rho)}$ 是外向柱面波，满足辐射条件；$H^{(1)}$ 是内向波。数学圈（$e^{-i\omega t}$）则反过来用 $H^{(1)}$——附录 E。

### A.3 2D 格林函数与 Graf 加法定理

$$\nabla^2G+k^2G=-\delta\ \Rightarrow\ G(\mathbf r,\mathbf r')=\frac{1}{4j}H_0^{(2)}(kR)=\frac{-i}{4}H_0^{(2)}(kR),\quad R=|\mathbf r-\mathbf r'|.$$

（归一化验证思路同 §1.3：小圆上积分散度项。）

**Graf 加法定理**（2D 版 QBX 的核心恒等式）：$|\boldsymbol\rho|<|\boldsymbol\rho'|$ 时

$$\boxed{\;H_0^{(2)}(k|\mathbf r-\mathbf r'|)=\sum_{\ell=-\infty}^{\infty}J_\ell(k\rho)\,H_\ell^{(2)}(k\rho')\,e^{i\ell(\varphi-\varphi')}\;}$$

结构与 3D 版完全平行：场点侧 $J_\ell e^{i\ell\varphi}$ 处处正则（可验证 $J_\ell(k\rho)e^{i\ell\varphi}$ 逐项满足 Helmholtz 方程——这正是"局部展开在球内是良解"的本质），源点侧 $H_\ell^{(2)}$ 光滑（$\rho'\ge h$），求和积分换序得系数公式

$$\alpha_\ell=\frac{1}{4j}\int_T H_\ell^{(2)}(k\rho')\,e^{-i\ell\varphi'}\,\sigma(\mathbf r')\,dS',
\qquad u(\mathbf r)=\sum_\ell\alpha_\ell J_\ell(k\rho)e^{i\ell\varphi}.$$

2D 代码是 3D 的"热身版"：建议先写 2D QBX 跑通，再上 3D。

---

## 附录 B 球 Bessel / 球 Hankel 函数（3D 径向函数）

### B.1 来源

3D Helmholtz 球坐标分离变量（附录 C 管角向），径向方程

$$\frac{d}{d\rho}\Big(\rho^2\frac{dR}{d\rho}\Big)+\big[k^2\rho^2-\ell(\ell+1)\big]R=0
\ \xrightarrow{x=k\rho}\ 
x^2R''+2xR'+\big[x^2-\ell(\ell+1)\big]R=0,$$

即**球 Bessel 方程**。其解与半阶 Bessel 函数的关系：$j_\ell(x)=\sqrt{\dfrac{\pi}{2x}}J_{\ell+1/2}(x)$，$y_\ell(x)=\sqrt{\dfrac{\pi}{2x}}Y_{\ell+1/2}(x)$。

### B.2 显式（Rayleigh 公式）

$$j_\ell(x)=(-x)^\ell\Big(\frac1x\frac{d}{dx}\Big)^\ell\frac{\sin x}{x},\qquad
y_\ell(x)=-(-x)^\ell\Big(\frac1x\frac{d}{dx}\Big)^\ell\frac{\cos x}{x}.$$

低阶显式（实现单元测试用）：

$$j_0=\frac{\sin x}{x},\quad j_1=\frac{\sin x}{x^2}-\frac{\cos x}{x},\quad j_2=\Big(\frac3{x^3}-\frac1x\Big)\sin x-\frac{3\cos x}{x^2},$$

$$y_0=-\frac{\cos x}{x},\quad y_1=-\frac{\cos x}{x^2}-\frac{\sin x}{x},\quad y_2=-\Big(\frac3{x^3}-\frac1x\Big)\cos x-\frac{3\sin x}{x^2}.$$

### B.3 球 Hankel 与有限闭式（3D 独有，宽频带矩表的根基）

$$h_\ell^{(2)}=j_\ell-i y_\ell:\qquad
h_0^{(2)}=\frac{i\,e^{-ix}}{x},\quad
h_1^{(2)}=e^{-ix}\Big(\frac{i}{x^2}-\frac1x\Big),\quad
h_2^{(2)}=e^{-ix}\Big(\frac{3i}{x^3}-\frac3{x^2}-\frac{i}{x}\Big).$$

（已数值验证。）**一般结构：$h_\ell^{(2)}(x)=\dfrac{e^{-ix}}{x}\times\big[\ell\text{ 阶的 }1/x\text{ 多项式}\big]$——有限项、精确的**，这是 3D 球 Hankel 与 2D Hankel 的本质区别（2D 只有渐近级数），也是 §6.8 频率无关矩表可行（3D-only）的数学根基。

注意 $h_0^{(2)}(kR)=i\,e^{-jkR}/(kR)$，于是

$$G(\mathbf r,\mathbf r')=\frac{e^{-jkR}}{4\pi R}=\frac{k}{4\pi i}\,h_0^{(2)}(kR)=-\frac{jk}{4\pi}h_0^{(2)}(kR),$$

即格林函数本身就是 $\ell=0$ 的球 Hankel——加法定理就是"把 $\ell=0$ 的偏心奇点展开成以 $\mathbf r_c$ 为心的全部 $\ell$ 的叠加"。

### B.4 小宗量行为（奇异性来源）

$$j_\ell(x)\sim\frac{x^\ell}{(2\ell+1)!!},\qquad
y_\ell(x)\sim-\frac{(2\ell-1)!!}{x^{\ell+1}},\qquad
h_\ell^{(2)}(x)\sim\frac{i\,(2\ell-1)!!}{x^{\ell+1}}\quad(x\to0).$$

- 场点侧 $j_\ell(k\rho)\sim\rho^\ell$：**正则**，级数每项在中心光滑（多极子场）；
- 源点侧 $h_\ell^{(2)}(k\rho')\sim\rho'^{-\ell-1}$：**只在 $\rho'=0$ 奇异**，离面中心保证 $\rho'\ge h>0$，奇异被移出积分域。这就是 QBX "把奇异性推远"的机制本身。

### B.5 递推与 Wronskian（数值实现用）

$$f_{\ell+1}=\frac{2\ell+1}{x}f_\ell-f_{\ell-1}\quad(f=j,y,h^{(1,2)}\text{ 均满足}),$$

$$j_\ell'(x)=j_{\ell-1}(x)-\frac{\ell+1}{x}j_\ell(x),\qquad
W\{j_\ell,y_\ell\}=j_\ell y_\ell'-j_\ell'y_\ell=\frac1{x^2}.$$

稳定性：$y_\ell$ 升递推无条件稳定；$j_\ell$ 在 $\ell>x$ 时升递推不稳定（真解按阶乘衰减被舍入误差淹没），用**降递推（Miller 算法）**：从远大于 $x$ 的 $\ell_{\max}$ 处令 $j=0,1$ 往下递推，再用 $j_0=\sin x/x$ 归一化。

---

## 附录 C Legendre 方程、Legendre 多项式与球谐函数

### C.1 角向分离 → Legendre 方程

3D Helmholtz 球坐标分离 $\psi=R(\rho)\Theta(\theta)\Phi(\varphi)$，角向给出

$$\frac{1}{\sin\theta}\frac{d}{d\theta}\Big(\sin\theta\frac{d\Theta}{d\theta}\Big)+\Big[\ell(\ell+1)-\frac{m^2}{\sin^2\theta}\Big]\Theta=0,\qquad \Phi=e^{im\varphi}.$$

$m=0$、令 $x=\cos\theta$ 即 **Legendre 方程**

$$(1-x^2)P''-2xP'+\ell(\ell+1)P=0.$$

### C.2 为什么 $\ell$ 必须是整数（量子化）

Legendre 方程对任意参数 $\ell$ 都有级数解，但一般解在 $x=\pm1$（**南北极点**）对数发散。物理场在极点必须有限 $\Rightarrow$ 级数必须中断为多项式 $\Rightarrow$ $\ell=0,1,2,\dots$，解为 **Legendre 多项式**

$$P_\ell(x)=\frac1{2^\ell\ell!}\frac{d^\ell}{dx^\ell}(x^2-1)^\ell\quad(\text{Rodrigues 公式}),\qquad P_0=1,\ P_1=x,\ P_2=\tfrac12(3x^2-1).$$

正交性 $\displaystyle\int_{-1}^1P_\ell P_{\ell'}\,dx=\frac{2}{2\ell+1}\delta_{\ell\ell'}$；端点值 $P_\ell(1)=1$，$P_\ell'(1)=\ell(\ell+1)/2$（§4.5 用过）。

### C.3 连带 Legendre 与球谐

$m\neq0$ 时解为连带 Legendre 函数（含 Condon–Shortley 相位）：

$$P_\ell^m(x)=(-1)^m(1-x^2)^{m/2}\frac{d^mP_\ell}{dx^m},\qquad
P_\ell^{-m}=(-1)^m\frac{(\ell-m)!}{(\ell+m)!}P_\ell^m,\quad m=0,\dots,\ell.$$

**球谐函数**（本手册归一化约定）

$$Y_\ell^m(\theta,\varphi)=\sqrt{\frac{(2\ell+1)(\ell-m)!}{4\pi(\ell+m)!}}\,P_\ell^m(\cos\theta)\,e^{im\varphi},\qquad
\int Y_\ell^m\bar Y_{\ell'}^{m'}\,dΩ=\delta_{\ell\ell'}\delta_{mm'}.$$

极点值（QBX 求值公式的来源）：$P_\ell^m(1)=0\ (m\ge1)$，故 $Y_\ell^m(0,\varphi)=\delta_{m0}\sqrt{(2\ell+1)/(4\pi)}$。

### C.4 球谐加法定理（不要和 Helmholtz 加法定理混淆）

设 $\gamma$ 为方向 $\hat a,\hat b$ 的夹角，则

$$P_\ell(\cos\gamma)=\frac{4\pi}{2\ell+1}\sum_{m=-\ell}^\ell Y_\ell^m(\hat a)\,\bar Y_\ell^m(\hat b).$$

这是纯几何恒等式（把 $P_\ell$ 这个"轴对称核"拆成球谐分量），附录 D 用它把 Legendre 级数转成 $\sum_{\ell m}$ 形式。特例（$\hat a=\hat b$）：$\sum_m|Y_\ell^m|^2=(2\ell+1)/(4\pi)$——§6.3 单元测试用的就是它。

### C.5 $\ell$ 的物理意义：多极子

$\ell$ 是角动量/多极阶数：$\ell=0$ 单极（点电荷），$\ell=1$ 偶极，$\ell=2$ 四极……QBX 级数 $u=\sum\alpha_{\ell m}j_\ell Y_\ell^m$ 的物理图像：**源分布对中心产生的多极展开**，$\alpha_{\ell m}$ 是各阶多极矩，$j_\ell(k\rho)Y_\ell^m$ 是该阶多极在球内的正则场。截断 $p$ = 保留到 $2^p$ 极子。

---

## 附录 D 3D Helmholtz 加法定理：完整推导（不跳步）

**目标**：证明 $|\boldsymbol\rho|<|\boldsymbol\rho'|$ 时

$$\frac{e^{-jk|\mathbf r-\mathbf r'|}}{4\pi|\mathbf r-\mathbf r'|}
=-jk\sum_{\ell=0}^\infty\sum_{m=-\ell}^\ell j_\ell(k\rho)Y_\ell^m(\hat\rho)\,h_\ell^{(2)}(k\rho')\bar Y_\ell^m(\hat\rho').$$

**第 1 步：角向展开。** 固定 $\mathbf r'$，把 $G$ 看作 $\mathbf r$ 的函数。它只依赖 $r=|\mathbf r|$ 与 $\mathbf r,\mathbf r'$ 夹角 $\gamma$，按 Legendre 多项式展开：

$$G=\sum_{\ell=0}^\infty g_\ell(r,r')\,P_\ell(\cos\gamma).$$

**第 2 步：径向 ODE。** $\mathbf r\neq\mathbf r'$ 时 $(\nabla^2+k^2)G=0$。对轴对称情形径向方程为（附录 B.1）

$$\frac1{r^2}\frac{d}{dr}\Big(r^2\frac{dg_\ell}{dr}\Big)+\Big(k^2-\frac{\ell(\ell+1)}{r^2}\Big)g_\ell=0
\ \Rightarrow\ g_\ell=A_\ell j_\ell(kr)+B_\ell y_\ell(kr).$$

**第 3 步：定边界条件。**
- $r<r'$ 区域包含原点：$y_\ell$ 在原点奇异，场在原点必须有限 $\Rightarrow B_\ell=0$，$g_\ell=A_\ell j_\ell(kr)$；
- $r>r'$ 区域延伸到无穷远：辐射条件要求外向波 $\Rightarrow g_\ell=C_\ell h_\ell^{(2)}(kr)$。
合并写法（自动满足两侧）：$g_\ell=D_\ell\,j_\ell(kr_<)\,h_\ell^{(2)}(kr_>)$，$r_<=\min(r,r')$，$r_>=\max(r,r')$。

**第 4 步：跃变条件定 $D_\ell$。** $\delta(\mathbf r-\mathbf r')$ 的角向展开（完备性）

$$\delta(\mathbf r-\mathbf r')=\frac{\delta(r-r')}{r^2}\sum_{\ell}\frac{2\ell+1}{4\pi}P_\ell(\cos\gamma).$$

把径向 ODE 写成自伴形式并在 $[r'-\epsilon,r'+\epsilon]$ 上积分：$\dfrac{d}{dr}\Big(r^2\dfrac{dg_\ell}{dr}\Big)+\big(k^2r^2-\ell(\ell+1)\big)g_\ell=-\dfrac{2\ell+1}{4\pi}\delta(r-r')$，得导数跃变

$$r'^2\Big[g_\ell'(r'^+)-g_\ell'(r'^-)\Big]=-\frac{2\ell+1}{4\pi}\quad(g_\ell\text{ 本身连续}).$$

代入第 3 步形式：

$$g_\ell'(r'^+)-g_\ell'(r'^-)=D_\ell k\Big[j_\ell(kr')h_\ell^{(2)\prime}(kr')-j_\ell'(kr')h_\ell^{(2)}(kr')\Big]
=D_\ell k\,W\{j_\ell,h_\ell^{(2)}\}.$$

由 $h^{(2)}=j-iy$ 与 $W\{j,y\}=1/x^2$（附录 B.5）：$W\{j,h^{(2)}\}=-i/x^2=-i/(kr')^2$。于是

$$D_\ell k\cdot\frac{-i}{k^2r'^2}=-\frac{2\ell+1}{4\pi r'^2}
\ \Rightarrow\ 
D_\ell=\frac{(2\ell+1)k}{4\pi i}=-\frac{ik(2\ell+1)}{4\pi}.$$

**第 5 步：合成并转球谐。**

$$G=-\frac{ik}{4\pi}\sum_\ell(2\ell+1)\,j_\ell(kr_<)h_\ell^{(2)}(kr_>)\,P_\ell(\cos\gamma).$$

用球谐加法定理（附录 C.4）$P_\ell(\cos\gamma)=\dfrac{4\pi}{2\ell+1}\sum_mY_\ell^m(\hat a)\bar Y_\ell^m(\hat b)$，取 $\hat a=\hat\rho,\hat b=\hat\rho'$、$r_<=\rho,r_>=\rho'$（条件 $|\boldsymbol\rho|<|\boldsymbol\rho'|$），即得目标公式。$\blacksquare$

**一致性检验**：令 $\boldsymbol\rho\to0$（场点就在中心），$j_\ell(0)=\delta_{\ell0}$，$Y_0^0=1/\sqrt{4\pi}$，右式退化为 $-jk\,h_0^{(2)}(k\rho')/(4\pi)=e^{-jk\rho'}/(4\pi\rho')$——正是 $\mathbf r=\mathbf r_c$ 处的格林函数值。✓（本手册写作时另以随机点做了数值验证，机器精度。）

---

## 附录 E 时间约定对照表

| 对象 | 物理约定（$e^{-i\omega t}$；Klöckner、Epstein、Wala 等数学文献） | 工程约定（$e^{j\omega t}$；本手册、你的代码、FEKO） |
|---|---|---|
| 外向波因子 | $e^{+ikR}$ | $e^{-jkR}$ |
| 3D 格林函数 | $e^{ikR}/4\pi R$ | $e^{-jkR}/4\pi R$ |
| 2D 格林函数 | $\frac{i}{4}H_0^{(1)}(kR)$ | $\frac{1}{4j}H_0^{(2)}(kR)$ |
| 径向外向函数 | $h_\ell^{(1)}=j_\ell+iy_\ell$ / $H_\ell^{(1)}$ | $h_\ell^{(2)}=j_\ell-iy_\ell$ / $H_\ell^{(2)}$ |
| 加法定理前因子 | $+ik$ | $-jk$ |
| $h_0$ 显式 | $h_0^{(1)}=-i\,e^{ix}/x$ | $h_0^{(2)}=i\,e^{-ix}/x$ |
| 转换规则 | — | 全文取复共轭（$j\to-i$，实数不变） |

**两条军规**：① 同一篇文档/同一份代码内只用一种约定；② 从数学文献抄公式时先整体共轭再抄，抄完用附录 D 的一致性检验（$\rho\to0$ 退化）验算。

---

*文档结束（修订版 v1.1，2026-08）。核心公式（加法定理取值、极点梯度 m=0/±1 结构、$h_1^{(2)},h_2^{(2)}$ 闭式）均经数值验证至机器精度；本版关键修订——内/外中心各自独立建架——经跳跃检验实测确认（$L^\pm=\pm0.500000$、跳跃 $=1.000000$，机器精度）。实现时按第 7 章清单逐级测试。*
