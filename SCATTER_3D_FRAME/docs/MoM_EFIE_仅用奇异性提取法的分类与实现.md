# MoM–EFIE 仅使用奇异性提取法时的分类与实现建议

## 1. 适用前提

本文假设只采用以下格林函数奇异性提取：

\[
G(R)
=
\frac{e^{-jkR}}{4\pi R}
=
\underbrace{\frac{1}{4\pi R}}_{G_{\mathrm{stat}}}
+
\underbrace{\frac{e^{-jkR}-1}{4\pi R}}_{G_{\mathrm{rem}}}.
\]

其中：

\[
G_{\mathrm{stat}}=\frac{1}{4\pi R},
\qquad
G_{\mathrm{rem}}=\frac{e^{-jkR}-1}{4\pi R}.
\]

\(G_{\mathrm{rem}}\) 在 \(R\to0\) 时有限：

\[
G_{\mathrm{rem}}(0)=-\frac{jk}{4\pi}.
\]

本方案不使用 Taylor–Duffy 等其他专用奇异积分公式，而是统一采用：

1. 静态 \(1/R\) 部分：对源三角形做解析内积分；
2. 剩余部分：用数值双重积分；
3. 对外层三角形根据真奇异位置或最近区域进行自适应细分。

---

## 2. 统一的静态内积分

对固定场点 \(\mathbf r\)，解析计算源三角形 \(T_s\) 上的零阶和一阶矩：

\[
I_0(\mathbf r)
=
\int_{T_s}
\frac{1}{|\mathbf r-\mathbf r'|}
\,dS',
\]

\[
\mathbf I_{\mathrm{vec}}(\mathbf r)
=
\int_{T_s}
\frac{\mathbf r'}{|\mathbf r-\mathbf r'|}
\,dS'.
\]

对于 RWG 矢量项，还可以直接形成：

\[
\int_{T_s}
\frac{\mathbf r'-\mathbf r_{\mathrm{opp},s}}
{|\mathbf r-\mathbf r'|}
\,dS'
=
\mathbf I_{\mathrm{vec}}(\mathbf r)
-
\mathbf r_{\mathrm{opp},s}I_0(\mathbf r).
\]

这样可将原来的四维静态积分转换成场三角形上的二维外积分。

---

## 3. 仅用奇异性提取时仍需五类分类

虽然各近场类型都使用同一种奇异性提取思想，但外层积分的非光滑位置不同，仍应分成五类：

\[
\boxed{
\text{同一三角形、共边、共顶点、非接触近场、普通远场}
}
\]

分类不是为了更换奇异性提取公式，而是为了决定：

- 是否使用奇异性提取；
- 外层三角形向哪里细分；
- 外层积分阶数和容差；
- 是否需要进行源场互换复核。

---

## 4. 第一类：同一三角形

### 判据

\[
N_{\mathrm{common}}=3.
\]

### 处理

对每个场三角形外层积分点 \(\mathbf r\)：

1. 解析计算源三角形的 \(I_0(\mathbf r)\)；
2. 解析计算 \(\mathbf I_{\mathrm{vec}}(\mathbf r)\)；
3. 形成静态 \(P\)、\(Q\) 项；
4. 对 \(G_{\mathrm{rem}}\) 使用双重数值积分。

即：

\[
Z_{\mathrm{self}}
=
Z_{\mathrm{stat,self}}
+
Z_{\mathrm{rem,self}}.
\]

### 外层积分要求

解析内积分消除了原始 \(1/R\) 点奇异，但所得势函数在三角形边界附近仍不是普通多项式。建议：

- 外层使用高阶三角形高斯积分；
- 比较 7 点、12 点和更高阶结果；
- 必要时对场三角形自适应细分。

不能只凭 7 点和12点结果接近就断定已收敛。

---

## 5. 第二类：共边三角形

### 判据

\[
N_{\mathrm{common}}=2.
\]

### 处理

仍采用相同分解：

\[
G=G_{\mathrm{stat}}+G_{\mathrm{rem}}.
\]

静态部分：

- 对每个场点 \(\mathbf r\in T_o\)，解析计算源三角形 \(T_s\) 的 \(I_0(\mathbf r)\) 和 \(\mathbf I_{\mathrm{vec}}(\mathbf r)\)；
- 再对场三角形做外层数值积分。

剩余部分：

- 对 \(G_{\mathrm{rem}}\) 做双重数值积分。

### 关键问题

奇异性提取虽然去掉了内层 \(1/R\) 奇异，但外层被积函数在共享边附近仍可能具有较强的梯度或非光滑性。因此不能简单地把整个场三角形固定使用 \(12\) 点积分。

### 推荐外层细分

沿共享边对场三角形和/或源三角形进行定向细分，使子三角形逐渐远离公共边。例如：

1. 找出公共边；
2. 将三角形按公共边中点或若干分点划分；
3. 对靠近公共边的子三角形继续递归细分；
4. 用相邻两级结果差控制误差。

---

## 6. 第三类：共顶点三角形

### 判据

\[
N_{\mathrm{common}}=1.
\]

### 处理

与共边类型相同，仍使用：

\[
G=G_{\mathrm{stat}}+G_{\mathrm{rem}}.
\]

静态内积分继续使用解析三角形势；剩余部分使用数值积分。

### 推荐外层细分

递归细分应集中在公共顶点附近：

1. 找出公共顶点；
2. 从公共顶点向对边划分子三角形；
3. 对包含公共顶点的子三角形继续细分；
4. 其余子三角形按普通高阶积分处理。

不能把共顶点三角形归入普通非奇异分支。

---

## 7. 第四类：无公共顶点但距离很近

### 判据

\[
N_{\mathrm{common}}=0,
\qquad
\chi=
\frac{d_{\min}}{h_{\mathrm{pair}}}
<
\chi_{\mathrm{near}}.
\]

其中：

\[
d_{\min}
=
\operatorname{dist}(T_o,T_s),
\qquad
h_{\mathrm{pair}}
=
\max(h_o,h_s).
\]

### 处理

仍使用同一种奇异性提取：

\[
G=G_{\mathrm{stat}}+G_{\mathrm{rem}}.
\]

静态部分采用解析内积分，可以显著减小近距离 \(1/R\) 峰值对双重高斯采样的要求。

### 外层自适应策略

1. 计算三角形对最近点；
2. 确定最近点位于顶点、边还是三角形内部；
3. 将场三角形向最近点区域细分；
4. 必要时也细分源三角形；
5. 递归计算低阶和高阶结果；
6. 达到误差容限后停止。

### 注意

奇异性提取会改善近场积分，但并不意味着固定12点一定足够。解析内积分之后，外层势函数仍可能在最近区域快速变化。

---

## 8. 第五类：普通远场

### 判据

\[
N_{\mathrm{common}}=0,
\qquad
\chi\geq\chi_{\mathrm{near}}.
\]

### 推荐处理

从效率角度，可以直接对完整格林函数使用普通双重高斯积分：

\[
G(R)=\frac{e^{-jkR}}{4\pi R}.
\]

若代码希望所有类别完全统一，也可以继续使用奇异性提取，但远场中这样会增加计算量：

\[
I_{\mathrm{full}}
=
I_{\mathrm{stat}}
+
I_{\mathrm{rem}}.
\]

两种写法理论上等价。建议：

- 近场和真奇异类别：强制奇异性提取；
- 普通远场：直接完整格林函数积分；
- 调试阶段可同时计算两种结果，检查二者一致性。

远场高斯阶数还应考虑：

\[
kh_{\mathrm{pair}}.
\]

---

## 9. 推荐的统一分类表

| 类型 | 判据 | 是否奇异性提取 | 静态部分 | 剩余部分 | 外层处理 |
|---|---|---:|---|---|---|
| 同一三角形 | \(N_{\mathrm{common}}=3\) | 必须 | 解析 \(I_0,\mathbf I_{\mathrm{vec}}\) | \(G_{\mathrm{rem}}\) 双重积分 | 高阶或自适应 |
| 共边 | \(N_{\mathrm{common}}=2\) | 必须 | 解析内积分 | \(G_{\mathrm{rem}}\) 双重积分 | 向公共边定向细分 |
| 共顶点 | \(N_{\mathrm{common}}=1\) | 必须 | 解析内积分 | \(G_{\mathrm{rem}}\) 双重积分 | 向公共顶点定向细分 |
| 非接触近场 | \(N_{\mathrm{common}}=0,\chi<\chi_{\mathrm{near}}\) | 推荐强制 | 解析内积分 | \(G_{\mathrm{rem}}\) 高阶/自适应 | 向最近区域细分 |
| 普通远场 | \(N_{\mathrm{common}}=0,\chi\geq\chi_{\mathrm{near}}\) | 可选 | 解析或不分解 | 普通高斯 | 根据 \(\chi,kh\) 选阶 |

---

## 10. 推荐程序流程

```fortran
N_COMMON = NUMBER_OF_COMMON_VERTICES(TRI_OBS, TRI_SRC)

IF (N_COMMON == 3) THEN

    MODE = EXTRACT_SELF

ELSE IF (N_COMMON == 2) THEN

    MODE = EXTRACT_COMMON_EDGE

ELSE IF (N_COMMON == 1) THEN

    MODE = EXTRACT_COMMON_VERTEX

ELSE

    D_MIN = TRIANGLE_TRIANGLE_DISTANCE(TRI_OBS, TRI_SRC)

    H_OBS = MAX_EDGE_LENGTH(TRI_OBS)
    H_SRC = MAX_EDGE_LENGTH(TRI_SRC)
    H_PAIR = MAX(H_OBS, H_SRC)

    CHI = D_MIN / H_PAIR

    IF (CHI < CHI_NEAR) THEN
        MODE = EXTRACT_NEAR
    ELSE
        MODE = REGULAR_FAR
    END IF

END IF
```

积分入口：

```fortran
SELECT CASE (MODE)

CASE (EXTRACT_SELF)

    CALL STATIC_ANALYTIC_INNER_INTEGRAL(...)
    CALL SMOOTH_REMAINDER_INTEGRAL(...)
    CALL OUTER_ADAPTIVE_SELF(...)

CASE (EXTRACT_COMMON_EDGE)

    CALL STATIC_ANALYTIC_INNER_INTEGRAL(...)
    CALL SMOOTH_REMAINDER_INTEGRAL(...)
    CALL OUTER_ADAPTIVE_TOWARD_EDGE(...)

CASE (EXTRACT_COMMON_VERTEX)

    CALL STATIC_ANALYTIC_INNER_INTEGRAL(...)
    CALL SMOOTH_REMAINDER_INTEGRAL(...)
    CALL OUTER_ADAPTIVE_TOWARD_VERTEX(...)

CASE (EXTRACT_NEAR)

    CALL STATIC_ANALYTIC_INNER_INTEGRAL(...)
    CALL SMOOTH_REMAINDER_INTEGRAL(...)
    CALL OUTER_ADAPTIVE_TOWARD_CLOSEST_REGION(...)

CASE (REGULAR_FAR)

    QUAD_ORDER = SELECT_ORDER_FROM_CHI_AND_KH(...)
    CALL REGULAR_FULL_GREEN_INTEGRAL(...)

END SELECT
```

---

## 11. 静态部分的统一形式

对 RWG 基函数

\[
\mathbf f_o(\mathbf r)
=
C_o(\mathbf r-\mathbf a_o),
\qquad
\mathbf f_s(\mathbf r')
=
C_s(\mathbf r'-\mathbf a_s),
\]

静态矢量项为：

\[
P_{\mathrm{stat}}
=
\frac{C_oC_s}{4\pi}
\int_{T_o}
(\mathbf r-\mathbf a_o)
\cdot
\left[
\mathbf I_{\mathrm{vec}}(\mathbf r)
-
\mathbf a_s I_0(\mathbf r)
\right]
dS.
\]

静态散度项为：

\[
Q_{\mathrm{raw,stat}}
=
\frac{4C_oC_s}{4\pi}
\int_{T_o}
I_0(\mathbf r)\,dS.
\]

之后再按照程序统一的 \(Q\) 符号约定乘频率系数。必须避免在 \(Q\) 定义和频率系数中重复吸收负号。

---

## 12. 外层自适应误差控制

对任意需要提取的三角形对，建议比较父三角形结果和一次细分后的结果：

\[
Z_{\mathrm{coarse}},
\qquad
Z_{\mathrm{refined}}.
\]

定义误差指标：

\[
\epsilon
=
\frac{
|Z_{\mathrm{refined}}-Z_{\mathrm{coarse}}|
}{
\max(|Z_{\mathrm{refined}}|,Z_{\mathrm{scale}})
}.
\]

若

\[
\epsilon>\epsilon_{\mathrm{quad}},
\]

则继续细分最靠近奇异集合或最近区域的子三角形。

建议同时对 \(P\) 项和 \(Q\) 项分别检查，而不是只看最终相加后的 \(Z\)，因为两项可能发生抵消。

---

## 13. 推荐阈值标定方式

先选取一组代表性三角形对：

- 自三角形；
- 共边；
- 共顶点；
- \(\chi=0.05,0.1,0.25,0.5,1,2\) 的非接触三角形；
- 不同 \(kh\) 的普通远场三角形。

然后比较：

- 固定7点；
- 固定12点；
- 更高阶；
- 一层细分；
- 多层自适应。

根据目标相对误差确定：

\[
\chi_{\mathrm{near}}
\quad\text{和}\quad
\epsilon_{\mathrm{quad}}.
\]

不要直接把某个经验阈值视为所有网格和频率下都适用。

---

## 14. 实施优先级

若希望在当前代码基础上逐步修改，建议顺序为：

1. 增加公共顶点数量判断；
2. 增加真实三角形最小距离和 \(\chi\)；
3. 将解析 \(I_0,\mathbf I_{\mathrm{vec}}\) 内积分扩展到所有需要提取的三角形对；
4. 共边外层向公共边细分；
5. 共顶点外层向公共顶点细分；
6. 非接触近场向最近区域细分；
7. 增加父级/细分级误差估计；
8. 最后根据 \(\chi\) 和 \(kh\) 优化远场高斯阶数。

---

## 15. 重要限制

仅做格林函数奇异性提取，并不等于所有奇异和近奇异积分已经自动准确。

奇异性提取主要解决：

\[
\frac{1}{R}
\]

在内层积分中的发散或尖峰问题。

它不能自动保证：

- 共边外层积分收敛；
- 共顶点外层积分收敛；
- 极近非接触三角形的固定高斯规则准确；
- \(G_{\mathrm{rem}}\) 在高频大单元上的积分阶数足够。

因此，“只用奇异性提取法”的完整实现仍需要：

\[
\boxed{
\text{拓扑分类}
+
\text{距离分类}
+
\text{外层定向自适应}
+
\text{误差控制}
}
\]

---

## 16. 最终推荐

若当前程序只实现这一种奇异性提取方法，可采用以下实用策略：

\[
\boxed{
\begin{aligned}
&N_{\mathrm{common}}=3:
&&\text{提取 + 自项外层高阶/自适应};\\
&N_{\mathrm{common}}=2:
&&\text{提取 + 向公共边细分};\\
&N_{\mathrm{common}}=1:
&&\text{提取 + 向公共顶点细分};\\
&N_{\mathrm{common}}=0,\ \chi<\chi_{\mathrm{near}}:
&&\text{提取 + 向最近区域细分};\\
&N_{\mathrm{common}}=0,\ \chi\geq\chi_{\mathrm{near}}:
&&\text{普通完整格林函数高斯积分}.
\end{aligned}
}
\]

这样不需要引入新的奇异积分理论框架，但仍能把同一种奇异性提取方法正确应用到不同拓扑和距离类型。
