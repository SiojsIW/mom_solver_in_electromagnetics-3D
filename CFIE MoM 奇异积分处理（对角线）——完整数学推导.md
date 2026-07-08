以下是重新整理的版本，**所有公式均独立成行，不嵌在表格里**，复制粘贴时不会截断。

---

# CFIE MoM 奇异积分处理（对角线）——完整数学推导

---

## 总览：七阶段计算依赖链

```
阶段1: 格林函数分解（纯定义，无依赖）
    ↓
阶段2: 解析标量势 I₀(r)（独立几何计算）
    ↓
阶段3: 解析线性势 Iⱼ(r)（独立几何计算，用 I₀ 验证）
    ↓
阶段4: EFIE 单三角形自阻抗（依赖 1, 2, 3）
    ↓
阶段5: EFIE 相邻三角形阻抗（依赖 1，数值双重积分）
    ↓
阶段6: MFIE 对角线元素（依赖 1，半解析+数值）
    ↓
阶段7: CFIE 组合（依赖 4, 5, 6）
```

---

## 阶段1：格林函数分解

**动机**：$G(R)$ 在 $R \to 0$ 时以 $1/R$ 强奇异，无法直接用高斯积分。将奇异部分解析提取，剩余部分光滑数值积分。

**分解公式**：

$$G(R)=\frac{e^{-jkR}}{4\pi R}=\underbrace{\frac{1}{4\pi R}}_{\text{强奇异，解析积分}}-\underbrace{\frac{k^{2}R}{8\pi}}_{\text{一阶导跳跃，通常数值处理}}+\underbrace{G_{\text{smooth2}}(R)}_{\text{二阶光滑，数值积分}}$$

**光滑核定义**：

$$G_{\text{smooth2}}(R)=\frac{e^{-jkR}-1+\dfrac{k^{2}R^{2}}{2}}{4\pi R}$$

**极限值**（供 $R=0$ 时直接返回）：

$$\lim_{R\to 0}G_{\text{smooth2}}(R)=-\frac{jk}{4\pi}$$

**Taylor 验证**：将 $e^{-jkR}=1-jkR-\frac{k^2R^2}{2}+\frac{jk^3R^3}{6}+O(R^4)$ 代入分子，得：

$$\text{分子}=\frac{jk^3R^3}{6}+O(R^4)\quad\Rightarrow\quad G_{\text{smooth2}}(R)=-\frac{jk}{4\pi}+\frac{jk^3R^2}{24\pi}+O(R^3)$$

分子除以 $R$ 后剩 $R^2$ 项，故二阶导数连续，可用普通高斯积分。

**输出**：供后续阶段调用 `GREEN_FUNC_SMOOTH2(R, K)`。

---

## 阶段2：解析标量势 $I_0(\mathbf{r})$

**动机**：计算 $\displaystyle\int_T \frac{1}{|\mathbf{r}-\mathbf{r}'|}dS'$，这是 $1/R$ 奇异在单三角形上的解析原函数，供后续双重积分的外层使用。

**输入**：场点 $\mathbf{r}$，三角形顶点 $\mathbf{v}_1,\mathbf{v}_2,\mathbf{v}_3$。

### 2.1 几何预处理

$$\hat{\mathbf{n}}=\frac{(\mathbf{v}_2-\mathbf{v}_1)\times(\mathbf{v}_3-\mathbf{v}_1)}{|(\mathbf{v}_2-\mathbf{v}_1)\times(\mathbf{v}_3-\mathbf{v}_1)|}$$

$$A=\frac{1}{2}|(\mathbf{v}_2-\mathbf{v}_1)\times(\mathbf{v}_3-\mathbf{v}_1)|$$

### 2.2 对第 $i$ 条边（从 $\mathbf{v}_i$ 到 $\mathbf{v}_{i+1}$，下标模3）

$$\hat{\mathbf{t}}_i=\frac{\mathbf{v}_{i+1}-\mathbf{v}_i}{L_i},\quad L_i=|\mathbf{v}_{i+1}-\mathbf{v}_i|$$

$$\hat{\mathbf{m}}_i=\hat{\mathbf{n}}\times\hat{\mathbf{t}}_i\quad\text{（边内法向，指向三角形内部）}$$

**方向修正**：取对顶点 $\mathbf{v}_{\text{opp}}$（不是该边端点的那个顶点），计算：

$$\text{DOT}=(\mathbf{v}_{\text{opp}}-\mathbf{v}_i)\cdot\hat{\mathbf{m}}_i$$

若 $\text{DOT}<0$，则 $\hat{\mathbf{m}}_i\leftarrow -\hat{\mathbf{m}}_i$。

### 2.3 有符号几何量

$$h_i=(\mathbf{r}-\mathbf{v}_i)\cdot\hat{\mathbf{m}}_i\quad\text{（场点到边所在直线的有符号距离）}$$

$$s_i^{\text{start}}=(\mathbf{r}-\mathbf{v}_i)\cdot\hat{\mathbf{t}}_i,\quad s_i^{\text{end}}=s_i^{\text{start}}-L_i$$

$$R_i^{\text{start}}=|\mathbf{r}-\mathbf{v}_i|,\quad R_i^{\text{end}}=|\mathbf{r}-\mathbf{v}_{i+1}|$$

### 2.4 单条边贡献

若 $R_i^{\text{start}}+s_i^{\text{start}}>10^{-12}$：

$$I_0^{(i)}(\mathbf{r})=h_i\ln\!\left(\frac{R_i^{\text{end}}+s_i^{\text{end}}}{R_i^{\text{start}}+s_i^{\text{start}}}\right)$$

否则（数值稳定形式）：

$$I_0^{(i)}(\mathbf{r})=h_i\ln\!\left(\frac{R_i^{\text{start}}-s_i^{\text{start}}}{R_i^{\text{end}}-s_i^{\text{end}}}\right)$$

### 2.5 总标量势

$$I_0(\mathbf{r})=\sum_{i=1}^{3}I_0^{(i)}(\mathbf{r})$$

**输出**：$I_0(\mathbf{r})$ 供阶段4的 Q 算子解析部分使用。

---

## 阶段3：解析线性势 $I_j(\mathbf{r})$

**动机**：RWG 基函数是线性的，需要 $\displaystyle\int_T \frac{\xi_j(\mathbf{r}')}{|\mathbf{r}-\mathbf{r}'|}dS'$，其中 $\xi_j$ 是面积坐标（在顶点 $j$ 处为1，其余为0）。

**输入**：同阶段2。

### 3.1 面积坐标与顶点关系

$$\mathbf{r}'=\sum_{j=1}^{3}\xi_j(\mathbf{r}')\mathbf{v}_j,\qquad \sum_{j=1}^{3}\xi_j=1$$

### 3.2 线性势定义

$$I_j(\mathbf{r})=\int_T\frac{\xi_j(\mathbf{r}')}{|\mathbf{r}-\mathbf{r}'|}dS',\quad j=1,2,3$$

**实现逻辑**：复用阶段2的边循环框架，对每条边按 Graglia 1993 Eq.(25)–(27) 为三个顶点分别累加贡献。每条边的贡献形式为：

$$\text{CONTRIB}_{i,j}= \frac{h_i^2}{2}\xi_j^{\text{(edge)}}\ln\!\left(\frac{R_i^{\text{end}}+s_i^{\text{end}}}{R_i^{\text{start}}+s_i^{\text{start}}}\right)+\frac{1}{2}\xi_j^{\text{(edge)}}\Big(s_i^{\text{end}}R_i^{\text{end}}-s_i^{\text{start}}R_i^{\text{start}}\Big)+\text{(顶点位置线性项)}$$

其中 $\xi_j^{\text{(edge)}}$ 在边起点处为 $\delta_{j,\text{start}}$，在终点处为 $\delta_{j,\text{end}}$。

### 3.3 验证恒等式（最关键）

$$I_1(\mathbf{r})+I_2(\mathbf{r})+I_3(\mathbf{r})\equiv I_0(\mathbf{r})\qquad\text{（因 }\xi_1+\xi_2+\xi_3=1\text{）}$$

**对称性验证**：等边三角形 + 重心场点 $\Rightarrow I_1=I_2=I_3=I_0/3$。

若上述两条不满足，说明公式实现有误。

**输出**：$I_1,I_2,I_3$ 供阶段4构造向量势 $\mathbf{I}_{\text{vec}}$。

---

## 阶段4：EFIE 单三角形自阻抗 $(T,T)$

**动机**：计算一个 RWG 基函数在单个三角形上"自己对自己"的阻抗贡献。这是整个对角线中最奇异的项，必须分解为 $1/R$ 解析部分 + $G_{\text{smooth2}}$ 数值部分。

**输入**：三角形 $T$（顶点 $\mathbf{v}_1,\mathbf{v}_2,\mathbf{v}_3$，面积 $A$），RWG 参数，波数 $k$，本征阻抗 $\eta_0$。

### 4.1 提取 RWG 局部参数

$$\mathbf{f}(\mathbf{r})=C(\mathbf{r}-\mathbf{r}_{\text{opp}}),\qquad \nabla\cdot\mathbf{f}=2C$$

其中：
- 若 $T=T^+$（正三角形）：$C=\dfrac{l}{2A}$，$\mathbf{r}_{\text{opp}}=\mathbf{r}_{\text{opp}}^+$
- 若 $T=T^-$（负三角形）：$C=-\dfrac{l}{2A}$，$\mathbf{r}_{\text{opp}}=\mathbf{r}_{\text{opp}}^-$

> **注意符号**：负三角形的 $C$ 本身带负号，这是 RWG 连续性要求的。

### 4.2 预计算常数

$$\text{JW\_MU}=jk\eta_0$$

$$\text{J\_WE}=\frac{j\eta_0}{k}$$

### 4.3 解析部分：外层高斯积分 + 内层解析积分

**核心思想**：对 $1/R$ 项，内层（对 $\mathbf{r}'$）用解析公式 $I_0,I_j$，外层（对 $\mathbf{r}$）用高斯积分。

#### 向量势构造

$$\mathbf{I}_{\text{vec}}(\mathbf{r})=\int_T\frac{\mathbf{r}'}{|\mathbf{r}-\mathbf{r}'|}dS'=\mathbf{v}_1I_1(\mathbf{r})+\mathbf{v}_2I_2(\mathbf{r})+\mathbf{v}_3I_3(\mathbf{r})$$

#### 外层高斯累积（对 $1/R$ 项）

取 $N_p$ 个高斯点 $\mathbf{r}_p$（建议7点），权重 $w_p$，统一乘以面积 $A$：

$$J_0=A\sum_{p=1}^{N_p}w_p\,I_0(\mathbf{r}_p)=\int_T\int_T\frac{1}{|\mathbf{r}-\mathbf{r}'|}dS'dS$$

$$\mathbf{J}_{\mathbf{r}I_0}=A\sum_{p=1}^{N_p}w_p\,\mathbf{r}_p\,I_0(\mathbf{r}_p)=\int_T\int_T\frac{\mathbf{r}}{|\mathbf{r}-\mathbf{r}'|}dS'dS$$

$$\mathbf{J}_{\text{vec}}=A\sum_{p=1}^{N_p}w_p\,\mathbf{I}_{\text{vec}}(\mathbf{r}_p)=\int_T\int_T\frac{\mathbf{r}'}{|\mathbf{r}-\mathbf{r}'|}dS'dS$$

$$J_{\mathbf{r}\cdot\text{vec}}=A\sum_{p=1}^{N_p}w_p\,\big[\mathbf{r}_p\cdot\mathbf{I}_{\text{vec}}(\mathbf{r}_p)\big]=\int_T\int_T\frac{\mathbf{r}\cdot\mathbf{r}'}{|\mathbf{r}-\mathbf{r}'|}dS'dS$$

### 4.4 Q 算子 $1/R$ 解析部分

Q 算子（散度-散度项）的原始形式：

$$Q_{\text{raw}}=\int_T\int_T (\nabla\cdot\mathbf{f})(\nabla'\cdot\mathbf{f})\frac{1}{4\pi|\mathbf{r}-\mathbf{r}'|}dS'dS=4C^2\int_T\int_T\frac{1}{4\pi|\mathbf{r}-\mathbf{r}'|}dS'dS$$

文档把负号和 $4C^2$ 一并封装：

$$Q_{\text{sing}}=-4C^2\cdot\frac{1}{4\pi}\cdot J_0=-\frac{C^2}{\pi}J_0$$

> **逻辑**：Q 算子在 EFIE 中整体带负号（$-\frac{j\eta_0}{k}Q_{\text{raw}}$），文档把 $-4C^2$ 提前放入 $Q_{\text{sing}}$，后续总公式用加号。

### 4.5 P 算子 $1/R$ 解析部分

P 算子（向量-向量项）的核展开：

$$(\mathbf{r}-\mathbf{r}_{\text{opp}})\cdot(\mathbf{r}'-\mathbf{r}_{\text{opp}})=\mathbf{r}\cdot\mathbf{r}'-\mathbf{r}\cdot\mathbf{r}_{\text{opp}}-\mathbf{r}'\cdot\mathbf{r}_{\text{opp}}+|\mathbf{r}_{\text{opp}}|^2$$

对应四项积分，分别记为 TERM1 ~ TERM4：

**TERM1**（对应 $\mathbf{r}\cdot\mathbf{r}'$ 项）：

$$\text{TERM1}=\int_T\int_T\frac{\mathbf{r}\cdot\mathbf{r}'}{|\mathbf{r}-\mathbf{r}'|}dS'dS=J_{\mathbf{r}\cdot\text{vec}}$$

**TERM2**（对应 $-\mathbf{r}\cdot\mathbf{r}_{\text{opp}}$ 项）：

$$\text{TERM2}=-\mathbf{r}_{\text{opp}}\cdot\int_T\int_T\frac{\mathbf{r}}{|\mathbf{r}-\mathbf{r}'|}dS'dS=-\mathbf{r}_{\text{opp}}\cdot\mathbf{J}_{\mathbf{r}I_0}$$

**TERM3**（对应 $-\mathbf{r}'\cdot\mathbf{r}_{\text{opp}}$ 项）：

$$\text{TERM3}=-\mathbf{r}_{\text{opp}}\cdot\int_T\int_T\frac{\mathbf{r}'}{|\mathbf{r}-\mathbf{r}'|}dS'dS=-\mathbf{r}_{\text{opp}}\cdot\mathbf{J}_{\text{vec}}$$

**TERM4**（对应 $|\mathbf{r}_{\text{opp}}|^2$ 项）：

$$\text{TERM4}=|\mathbf{r}_{\text{opp}}|^2\int_T\int_T\frac{1}{|\mathbf{r}-\mathbf{r}'|}dS'dS=|\mathbf{r}_{\text{opp}}|^2 J_0$$

**P 算子解析部分**：

$$P_{\text{sing}}=C^2\big(\text{TERM1}+\text{TERM2}+\text{TERM3}+\text{TERM4}\big)$$

> **注意**：这里 $P_{\text{sing}}$ 是否缺少 $1/4\pi$ 因子需要与现有代码核对。若你之前的非对角线 $P$ 包含 $1/4\pi$，则此处四项也应乘以 $1/4\pi$。文档在此处写法可能省略了该因子，建议以格林函数 $G$ 的完整定义为准。

### 4.6 $G_{\text{smooth2}}$ 数值部分（双重高斯积分）

对光滑核 $G_{\text{smooth2}}$，$(T,T)$ 上直接做双重高斯积分（如 $7\times 7$ 点）：

$$I_1^{\text{sm}}=\int_T\int_T G_{\text{smooth2}}(|\mathbf{r}-\mathbf{r}'|)\,dS'dS$$

$$\mathbf{I}_2^{\text{sm}}=\int_T\int_T \mathbf{r}\,G_{\text{smooth2}}\,dS'dS$$

$$\mathbf{I}_3^{\text{sm}}=\int_T\int_T \mathbf{r}'\,G_{\text{smooth2}}\,dS'dS$$

$$I_4^{\text{sm}}=\int_T\int_T (\mathbf{r}\cdot\mathbf{r}')\,G_{\text{smooth2}}\,dS'dS$$

然后按非对角线公式结构组合：

$$P_{\text{smooth2}}=C^2\Big[|\mathbf{r}_{\text{opp}}|^2 I_1^{\text{sm}}-\mathbf{r}_{\text{opp}}\cdot\mathbf{I}_2^{\text{sm}}-\mathbf{r}_{\text{opp}}\cdot\mathbf{I}_3^{\text{sm}}+I_4^{\text{sm}}\Big]$$

$$Q_{\text{smooth2}}=-4C^2\cdot I_1^{\text{sm}}$$

### 4.7 总单三角形阻抗

$$Z_{(T,T)}^E=C^2\Big[\text{JW\_MU}\cdot(P_{\text{sing}}+P_{\text{smooth2}})+\text{J\_WE}\cdot(Q_{\text{sing}}+Q_{\text{smooth2}})\Big]$$

展开即：

$$Z_{(T,T)}^E=jk\eta_0 C^2(P_{\text{sing}}+P_{\text{smooth2}})+\frac{j\eta_0}{k}(Q_{\text{sing}}+Q_{\text{smooth2}})$$

**输出**：$Z_{(T,T)}^E$ 供阶段7的 EFIE 对角线组装。

---

## 阶段5：EFIE 相邻三角形阻抗 $(T^+,T^-)$ 与 $(T^-,T^+)$

**动机**：共享一条边的两个三角形之间是弱奇异（$1/R$ 可积），不需要解析提取，直接用完整格林函数数值积分即可。

**输入**：RWG 基函数（含 $T^+,T^-$），$G_{\text{full}}=\dfrac{e^{-jkR}}{4\pi R}$。

### 5.1 提取系数

$$C_{\text{pos}}=\frac{l}{2A^+}$$

$$C_{\text{neg}}=-\frac{l}{2A^-}$$

### 5.2 双重高斯积分（完整格林函数）

对 $(T^+,T^-)$ 用 $7\times 7$ 高斯点计算：

$$Z_{(T^+,T^-)}^E=\text{按非对角线公式，核用 }G_{\text{full}}$$

对 $(T^-,T^+)$ 同理得 $Z_{(T^-,T^+)}^E$。

### 5.3 相邻对总和

$$Z_{\text{adj}}^E=Z_{(T^+,T^-)}^E+Z_{(T^-,T^+)}^E$$

**输出**：$Z_{\text{adj}}^E$ 供阶段7。

---

## 阶段6：MFIE 对角线元素 $Z_{mm}^M$

**动机**：MFIE 的矩阵元包含 K 算子（主值积分）和 $1/2$ 内积项。对角线需要区分"同一三角形对"与"相邻三角形对"。

**标准形式**（供对照）：

$$\hat{\mathbf{n}}\times\mathbf{H}^{\text{scat}}=\mathbf{K}(\mathbf{J})-\frac{1}{2}\mathbf{J}$$

测试后：

$$Z_{mn}^M=\langle\mathbf{f}_m,\mathbf{K}(\mathbf{f}_n)\rangle-\frac{1}{2}\langle\mathbf{f}_m,\mathbf{f}_n\rangle$$

文档采用等价形式（符号已调整）：

$$Z_{mm}^M=\underbrace{\text{HALF\_TERM}}_{1/2\text{ 内积项}}+\underbrace{K_{\text{TERM}}}_{\text{K算子}}$$

### 6.1 同一三角形对 $(T^+,T^+)$ 与 $(T^-,T^-)$

**K 算子主值**：**严格置 0**。

> 逻辑：平面三角形上，K 算子的柯西主值在自片积分为零（$\nabla G\times\mathbf{f}$ 的奇异性在对称积分下抵消）。

**$1/2$ 内积项**：

$$\int_{T^+}|\mathbf{f}^+|^2dS=C_{\text{pos}}^2\int_{T^+}|\mathbf{r}-\mathbf{r}_{\text{opp}}^{\text{pos}}|^2dS=C_{\text{pos}}^2\,A^+\sum_p w_p|\mathbf{r}_p-\mathbf{r}_{\text{opp}}^{\text{pos}}|^2$$

$$\int_{T^-}|\mathbf{f}^-|^2dS=C_{\text{neg}}^2\int_{T^-}|\mathbf{r}-\mathbf{r}_{\text{opp}}^{\text{neg}}|^2dS=C_{\text{neg}}^2\,A^-\sum_p w_p|\mathbf{r}_p-\mathbf{r}_{\text{opp}}^{\text{neg}}|^2$$

$$\text{HALF\_TERM}=\frac{1}{2}\left(\int_{T^+}|\mathbf{f}^+|^2dS+\int_{T^-}|\mathbf{f}^-|^2dS\right)$$

### 6.2 相邻三角形对 $(T^+,T^-)$ 与 $(T^-,T^+)$

**$1/2$ 内积项**：**严格为 0**。

> 逻辑：$\mathbf{f}^+$ 与 $\mathbf{f}^-$ 的支撑域交集仅为公共边（一维测度为零），面积分为零。

**K 算子**：复用现有 `CALC_MFIE_MATRIX_ELEMENT` 的 K 算子逻辑，对 $(T^+,T^-)$ 和 $(T^-,T^+)$ 分别做高斯积分：

$$\mathbf{K}(\mathbf{f})=\int_T \nabla G(|\mathbf{r}-\mathbf{r}'|)\times\mathbf{f}(\mathbf{r}')\,dS'$$

然后测试积分：

$$K_{(T^+,T^-)}=\int_{T^+}\mathbf{f}^+(\mathbf{r})\cdot\Big[\int_{T^-}\nabla G\times\mathbf{f}^-(\mathbf{r}')\,dS'\Big]dS$$

$$K_{(T^-,T^+)}=\int_{T^-}\mathbf{f}^-(\mathbf{r})\cdot\Big[\int_{T^+}\nabla G\times\mathbf{f}^+(\mathbf{r}')\,dS'\Big]dS$$

$$K_{\text{TERM}}=K_{(T^+,T^-)}+K_{(T^-,T^+)}$$

### 6.3 总 MFIE 对角线

$$Z_{mm}^M=\text{HALF\_TERM}+K_{\text{TERM}}$$

**验证标准**：$Z_{mm}^M$ 应为**正实数**（$1/2$ 项占主导），模值约 $10^{-2}\sim 10^{-1}$（正四面体尺度）。

**输出**：$Z_{mm}^M$ 供阶段7。

---

## 阶段7：CFIE 组合与主程序集成

**动机**：CFIE 消除内谐振，对角线元素直接组合 EFIE 与 MFIE。

### 7.1 EFIE 对角线完整组装

对每个 RWG 基函数 $m$：

$$Z_{mm}^E=Z_{(T_m^+,T_m^+)}^E+Z_{(T_m^-,T_m^-)}^E+Z_{\text{adj},m}^E$$

其中第一项和第二项来自阶段4，第三项来自阶段5。

### 7.2 MFIE 对角线

直接取阶段6结果 $Z_{mm}^M$。

### 7.3 CFIE 组合

$$Z_{mm}^{\text{CFIE}}=\alpha\,Z_{mm}^E+(1-\alpha)\,\eta_0\,Z_{mm}^M$$

通常取 $\alpha=0.5$。

---

## 关键易错点与概念陷阱

| 陷阱                               | 说明                                                         |
| :--------------------------------- | :----------------------------------------------------------- |
| $P_{\text{sing}}$ 的 $1/4\pi$ 因子 | 阶段4的 $P_{\text{sing}}$ 公式若严格对应格林函数 $G=\frac{1}{4\pi R}+\dots$，则四项累积后应显式乘以 $\frac{1}{4\pi}$。请核对你代码中非对角线 $P$ 的定义，确保对角线与非对角线的 $1/4\pi$ 处理一致。 |
| $C$ 的符号                         | 负三角形的 $C$ 为负值（$C=-l/2A^-$），但 $C^2$ 为正。在 P 算子的双线性展开中，$(\mathbf{r}-\mathbf{r}_{\text{opp}})$ 的向量方向不要额外取反，符号已由 $C$ 携带。 |
| $Q_{\text{sing}}$ 的负号           | 文档把 $-4C^2$ 封装进 $Q_{\text{sing}}$，因此总公式中是 `+ J_WE * (Q_SING + Q_SMOOTH2)`，而非减号。若你习惯标准形式 $-\frac{j\eta_0}{k}Q_{\text{raw}}$，请确认符号封装方式。 |
| 解析线性势的恒等式                 | $I_1+I_2+I_3\equiv I_0$ 是**必要条件**。若数值不满足，100% 是 Graglia 公式实现错误，不要怀疑几何数据。 |
| MFIE 的 K 算子自片                 | $(T^+,T^+)$ 和 $(T^-,T^-)$ 的 K 算子主值必须**严格置零**，不要试图数值计算（会发散且主值为零）。 |
| $G_{\text{smooth2}}$ 的 $R=0$ 极限 | 当高斯点重合（同一三角形的双重积分中可能出现），必须返回 $-jk/(4\pi)$，不能留空或设零。 |

---

## 执行顺序与检查清单

| 阶段 | 任务                               | 验证点                                                       |
| :--- | :--------------------------------- | :----------------------------------------------------------- |
| 1    | GREEN_FUNC_SMOOTH2                 | $R=0$ 时虚部 $=-k/(4\pi)$                                    |
| 2    | ANALYTIC_SCALAR_POT                | 等边+重心 $>0$，与50点高斯误差 $<10^{-10}$                   |
| 3    | ANALYTIC_LINEAR_POT                | $I_1+I_2+I_3=I_0$；等边+重心 $I_1=I_2=I_3$                   |
| 4    | CALC_MFIE_DIAGONAL_ELEMENT         | 结果为正实数，模值 $10^{-2}\sim10^{-1}$                      |
| 5    | CALC_EFIE_SELF_TRI_PAIR（仅Q算子） | $Q_{\text{sing}}$ 为负实数，$|Q_{\text{sing}}|>|Q_{\text{smooth2}}|$ |
| 6    | CALC_EFIE_SELF_TRI_PAIR（完整P+Q） | $Z_{\text{pair}}$ 复数，实部为负，模值 $\gg$ 非对角          |
| 7    | CALC_EFIE_ADJACENT_TRI_PAIR        | 结果复数，模值 $<$ 同一三角形                                |
| 8    | 主程序集成 + CFIE                  | 对角线填充后矩阵条件数 $< 10^6$                              |

---

**建议**：按阶段 $1\to 2\to 3\to 4$ 的顺序逐个验证，每通过一个再进入下一个。把阶段2和阶段3的验证结果确认后，再进入阶段4的 MFIE 对角线。