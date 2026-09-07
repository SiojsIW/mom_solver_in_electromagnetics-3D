# 奇异与近奇异积分处理方法

本求解器的 MoM 阻抗矩阵填充中，最核心、最易数值失稳的部分是三角形面片之间的 **Green 函数及其梯度** 的面积分。当场三角形 `m` 与源三角形 `n` 距离很近甚至重合时，被积函数出现 `1/R` 或 `1/R³` 奇异性，必须采用解析提取 + 数值正则化的混合策略。

---

## 1. 整体策略

程序在 `Input_setting(2)` 中设置距离阈值：

```fortran
threshold_R = 0.1 * lambda
```

在 `MoM_SS.f90` 中，对每一对三角形先算重心距：

```fortran
dist = re_Norm(rcm - rcn)
```

然后以 `dist >= threshold_R` 作为逻辑开关，分别调用奇异分支和非奇异分支：

```fortran
call Comp_eltsI(I1ss, I2ss, I3ss, I4ss, dist>=threshold_R, ...)
if(CFSIE .and. m/=n) then
    call Comp_eltsU(U1ss, U2ss, U3ss, U4ss, U5ss, dist>=threshold_R, ...)
endif
```

`dist >= threshold_R` 为 `.TRUE.` 时走**远场纯数值积分**；为 `.FALSE.` 时走**近奇异/奇异解析提取**。

---

## 2. 标量 Green 函数积分：`Comp_eltsI.f90`

### 2.1 要计算的积分

`Comp_eltsI` 输出四个量，对应 EFIE 中的标量势和矢量势积分：

```
I1 = ∫_Tm ∫_Tn G(R) dS' dS
I2 = ∫_Tm ∫_Tn r · G(R) dS' dS
I3 = ∫_Tm ∫_Tn r' · G(R) dS' dS
I4 = ∫_Tm ∫_Tn r · r' · G(R) dS' dS
R  = |r − r'|
G(R) = exp(−jkR) / R
```

### 2.2 远场分支（非奇异）

```fortran
if(nonsingularity) then
    do q=1,Nq
        wq = weightS(q)
        rq = ris(:,q)
        R  = re_Norm(rp-rq)
        ctmp = wq * Green(R)
        inner   = inner   + ctmp
        innerv  = innerv  + ctmp * rq
    enddo
endif
```

其中 `Green(R)` 直接计算自由空间 Green 函数：

```fortran
pure complex function Green(R)
    Green = exp(cmplx(0.0, -k_wnum*R)) / R
endfunction
```

### 2.3 近场/奇异分支（奇异提取）

核心思想：把 Green 函数拆成 `1/R` 奇异部分 + 正则余项 `F1(R)`：

```
exp(−jkR)/R = 1/R + F1(R)
```

```fortran
else
    !--- 数值计算正则余项 F1(R)
    do q=1,Nq
        wq = weightS(q)
        rq = ris(:,q)
        R  = re_Norm(rp-rq)
        F1_temp = wq * F1(R)
        inner   = inner   + F1_temp
        innerv  = innerv  + F1_temp * rq
    enddo

    !--- 解析计算 1/R 的奇异面积分
    call Spara_cal(rp, r1s, r2s, r3s, avi, t3, t4, w0, ns)
    t3 = t3 / avi
    t4 = t4 / avi
    innerv = innerv + (rp - ns*w0) * t3

    inner  = inner  + t3
    innerv = innerv + t4
endif
```

这里 `avi = An` 是源三角形面积，`rp` 是场点，`r1s/r2s/r3s` 是源三角形顶点。

### 2.4 正则余项 `F1(R)` 的泰勒展开

```fortran
pure complex function F1(R)
    real,intent(in):: R
    real x
    complex,parameter:: GreenExpan_coeff(0:7) = &
        (/(0,-1.0), (-0.5,0), (0,0.1666667), (4.1666667e-2,0), &
          (0.0,-8.3333333e-3), (-1.3888889e-3,0), (0,1.9841270e-4), &
          (2.4801587e-5,0)/)
    x = k_wnum * R
    F1 = GreenExpan_coeff(7)
    do i=6,0,-1
        F1 = F1 * x
        F1 = F1 + GreenExpan_coeff(i)
    enddo
    F1 = k_wnum * F1
endfunction
```

该 7 次多项式对应 `exp(−jx)/x − 1/x` 在 `x=0` 附近的泰勒展开，保证 `R→0` 时无奇异性。

### 2.5 解析奇异积分：`Spara_cal`

`Spara_cal` 计算源三角形上 `1/R` 的解析面积分，返回：

- `t3 = ∫_Tn 1/R dS'`
- `t4 = ∫_Tn r' / R dS'` 的相关分量
- `w0`：场点到源三角形平面的有向距离
- `ns`：源三角形单位法向

其实现基于**三角形局部坐标参数化 + 边积分**，把面积分化为对三条边的求和。

#### 局部坐标系建立

```fortran
call unit_normal(P1, P2, P3, ns)       ! 源三角形单位法向

abs_l(1) = re_Norm(P3-P2)              ! 边 1 长度
abs_l(2) = re_Norm(P1-P3)              ! 边 2 长度
abs_l(3) = re_Norm(P2-P1)              ! 边 3 长度

l1 = (P3-P2)/abs_l(1)                  ! 边方向
l2 = (P1-P3)/abs_l(2)
l3 = (P2-P1)/abs_l(3)

call Re_cross_product(l1, ns, m1)      ! 边外法向
call Re_cross_product(l2, ns, m2)
call Re_cross_product(l3, ns, m3)

u = l3
w = ns
call Re_cross_product(ns, u, v)        ! 建立 (u,v,w) 局部正交坐标系
```

#### 投影参数

```fortran
r_P1 = r - P1
P3_P1 = P3 - P1

u0 = dot_product(u, r_P1)              ! 场点投影的 u 坐标
v0 = dot_product(v, r_P1)              ! 场点投影的 v 坐标
w0 = dot_product(w, r_P1)              ! 场点到平面的有向距离
u3 = dot_product(u, P3_P1)             ! 第三个顶点 u 坐标
v3 = 2.0*areai / abs_l(3)              ! 第三个顶点 v 坐标
```

#### 每条边的参数

对三条边分别计算：

```fortran
lm(1) = -((abs_l(3)-u0)*(abs_l(3)-u3) + v0*v3) / abs_l(1)
lm(2) = -(u3*(u3-u0) + v3*(v3-v0)) / abs_l(2)
lm(3) = -u0

lp(1) = ((u3-u0)*(u3-abs_l(3)) + v3*(v3-v0)) / abs_l(1)
lp(2) = (u0*u3 + v0*v3) / abs_l(2)
lp(3) = abs_l(3) - u0

P0(1) = (v0*(u3-abs_l(3)) + v3*(abs_l(3)-u0)) / abs_l(1)
P0(2) = (u0*v3 - v0*u3) / abs_l(2)
P0(3) = v0
```

`lm`、`lp` 分别是场点投影到各边所在直线时，沿边方向的相对距离参数。

#### 解析公式

对每条边 `m`：

```fortran
R0 = P0(m)**2 + w0**2
Rm = sqrt(lm(m)**2 + P0(m)**2 + w0**2)
Rp = sqrt(lp(m)**2 + P0(m)**2 + w0**2)

!--- 对数项 g2
if(lm(m)>0 .and. lp(m)>0) then
    g2 = log((Rp+lp(m)) / (Rm+lm(m)))
    g3 = (lp(m)*Rp - lm(m)*Rm) + R0*g2
elseif(lm(m)<=0 .and. lp(m)>=0) then
    if(Rm+lm(m)==0) then
        singular = .TRUE.
    else
        g2 = log((Rp+lp(m)) / (Rm+lm(m)))
        g3 = (lp(m)*Rp - lm(m)*Rm) + R0*g2
    endif
elseif(lm(m)<0 .and. lp(m)<0) then
    g2 = log((Rm-lm(m)) / (Rp-lp(m)))
    g3 = (lp(m)*Rp - lm(m)*Rm) + R0*g2
endif

!--- 反正切项 beta（当 R0 很小时置 0 防止 NaN）
if(R0 <= (0.0e-2*abs_l(m))**2) then
    beta = 0.0
else
    beta = atan(P0(m)*lp(m) / (R0 + abs(w0)*Rp)) &
         - atan(P0(m)*lm(m) / (R0 + abs(w0)*Rm))
endif

!--- 合成该边贡献
if(.not.singular) then
    t1(m) = P0(m)*g2 - abs(w0)*beta
    t2(m) = 0.5*g3
else
    t1(m) = -abs(w0)*beta
    t2(m) = 0.5 * (lp(m)*Rp - lm(m)*Rm)
endif
```

最后三条边求和：

```fortran
t3 = sum(t1)
t4 = m1*t2(1) + m2*t2(2) + m3*t2(3)
```

这里 `t3` 就是 `∫_Tn 1/R dS'`，`t4` 用于后续恢复 `∫_Tn r'/R dS'`。

---

## 3. Green 函数梯度积分：`Comp_eltsU.f90`

### 3.1 要计算的积分

`Comp_eltsU` 输出对应 MFIE 中 `∇G` 相关的五个积分：

```
U1 = ∫ ∇G dS'
U2 = ∫ r × ∇G dS'
U3 = ∫ r · ∇G dS'
U4 = ∫ r × (r × ∇G) dS'
U5 = ∫ r · (n · ∇G) dS'
```

具体实现中通过 `gradG_S` 计算每个场点 `rp` 处源三角形上的 `∫_Tn ∇G dS'`。

### 3.2 远场分支

```fortran
if(FarNear) then
    gradG = (0.0,0.0)
    do q=1,SNum_samp
        rq = ris(:,q)
        R  = rp - rq
        abs_R = re_Norm(R)
        temp_gradient_G = -cmplx(1.0, k_wnum*abs_R) &
                          * exp(cmplx(0.0, -k_wnum*abs_R)) / (abs_R**3) * R
        wq = weightS(q)
        gradG = gradG + wq * temp_gradient_G
    enddo
endif
```

### 3.3 近场/奇异分支

```fortran
else
    I_30 = (0.0,0.0)

    !--- 普通正则项（数值积分）
    do q=1,SNum_samp
        rq = ris(:,q)
        R  = rp - rq
        abs_R = re_Norm(R)
        temp_gradient_G = R * cmplx(k_wnum*abs_R*(1.0 - k_wnum**2*abs_R**2/18.0)/8.0, &
                                    (1.0 - 0.1*k_wnum**2*abs_R**2)/3.0)
        wq = weightS(q)
        I_30 = I_30 + wq * temp_gradient_G
    enddo

    !--- 奇异项：主值（PV）解析计算
    call Spara_cal_2(rp, P1, P2, P3, An, t1, t3, t4, t5, w0, ns)

    I_31 = 0.0
    if(abs(w0) <= ratio_limit1) then       ! w0 近似为 0
        do i=1,3
            I_31 = I_31 + t5(:,i)
        enddo
    else
        do i=1,3
            I_31 = I_31 + t4(:,i) + t5(:,i)
        enddo
    endif

    I_32 = 0.0
    if(abs(w0) <= ratio_limit1) then       ! w0 近似为 0
        do i=1,3
            I_32 = I_32 - t3(:,i)
        enddo
    else
        do i=1,3
            I_32 = I_32 + t1(i)*w0*ns - t3(:,i)
        enddo
    endif

    gradG = k_wnum**3 * I_30 - (I_31 + k_wnum**2/2.0 * I_32) / An
endif
```

这里将 `∇G` 的积分按 `k` 的幂次分离：

- `I_30`：正则高阶项，数值积分
- `I_31`、`I_32`：解析计算的奇异/主值项

### 3.4 解析梯度主值：`Spara_cal_2`

`Spara_cal_2` 与 `Spara_cal` 几何参数化完全相同，但返回更多张量：

```fortran
subroutine Spara_cal_2(r, P1, P2, P3, An, t1, t3, t4, t5, w0, ns)
    real,intent(out):: t1(3), t3(3,3), t4(3,3), t5(3,3)
```

- `t1(3)`：与 `atan` 项相关的标量
- `t3(3,3)`：含 `R0²` 的向量项
- `t4(3,3)`：法向 `ns` 分量（仅 `w0≠0`）
- `t5(3,3)`：边外法向 `m_i` 分量

#### 主要区别：`w0≈0` 分支

```fortran
if(abs(w0) <= ratio_limit1) then     ! w0 为 0
    do i=1,3
        if(lm(i)>0 .and. lp(i)>0) then
            g2(i) = log((Rp(i)+lp(i)) / (Rm(i)+lm(i)))
            t3(:,i) = 0.5*m_u(:,i)*(R0(i)**2*g2(i) + Rp(i)*lp(i) - Rm(i)*lm(i))
            t5(:,i) = m_u(:,i)*g2(i)
        elseif(lm(i)<=0 .and. lp(i)>=0) then
            if(Rm(i)+lm(i)==0) then
                t3(:,i) = 0.5*m_u(:,i)*(Rp(i)*lp(i) - Rm(i)*lm(i))
                t5(:,i) = 0.0*m_u(:,i)
            else
                g2(i) = log((Rp(i)+lp(i)) / (Rm(i)+lm(i)))
                t3(:,i) = 0.5*m_u(:,i)*(R0(i)**2*g2(i) + Rp(i)*lp(i) - Rm(i)*lm(i))
                t5(:,i) = m_u(:,i)*g2(i)
            endif
        elseif(lm(i)<0 .and. lp(i)<0) then
            g2(i) = log((Rm(i)-lm(i)) / (Rp(i)-lp(i)))
            t3(:,i) = 0.5*m_u(:,i)*(R0(i)**2*g2(i) + Rp(i)*lp(i) - Rm(i)*lm(i))
            t5(:,i) = m_u(:,i)*g2(i)
        endif
    enddo
```

当 `w0=0` 时，含 `atan` 的 `t4` 项与 `t1` 项消失，避免 `0/0` 型奇异性。

#### `w0≠0` 分支

```fortran
else                                   ! w0 不为 0
    do i=1,3
        if(lm(i)>0 .and. lp(i)>0) then
            g2(i) = log((Rp(i)+lp(i)) / (Rm(i)+lm(i)))

            if(R0(i) <= ratio_limit2*abs_l(i)) then
                beta(i) = 0.0
            else
                beta(i) = atan(P0(i)*lp(i) / (R0(i)**2 + abs(w0)*Rp(i))) &
                        - atan(P0(i)*lm(i) / (R0(i)**2 + abs(w0)*Rm(i)))
            endif

            t1(i)    = P0(i)*g2(i) - abs(w0)*beta(i)
            t4(:,i)  = ns * sgn(w0) * beta(i)
            t3(:,i)  = 0.5*m_u(:,i)*(R0(i)**2*g2(i) + Rp(i)*lp(i) - Rm(i)*lm(i))
            t5(:,i)  = m_u(:,i)*g2(i)
        ...
```

`ratio_limit2 = 3.0e-2` 用于防止 `atan(有限/0)` 产生 NaN。

#### 分支汇总

| `lm`、`lp` 关系 | 对数处理 | 说明 |
|-----------------|----------|------|
| `lm>0, lp>0`    | `log((Rp+lp)/(Rm+lm))` | 标准情形 |
| `lm≤0, lp≥0`    | 若 `Rm+lm==0` 跳过对数 | 投影点落在边线段上，数值零保护 |
| `lm<0, lp<0`    | `log((Rm-lm)/(Rp-lp))` | 转换到正区间 |

---

## 4. 阈值参数

在 `Global_module.f90` / `Input_setting.f90` 中：

```fortran
real, parameter:: ratio_limit2 = 3.0e-2
real threshold_R, ratio_limit1

! 在 Input_setting(2) 中：
threshold_R = 0.1 * lambda
ratio_limit1 = 1.0e-6 * lambda   ! 判断 w0 是否近似为 0
```

这些阈值直接影响：

- `threshold_R`：多远算“近场”，触发奇异分支
- `ratio_limit1`：场点是否被认为落在源三角形平面内
- `ratio_limit2`：投影点是否被认为落在某条边上，决定是否跳过 `atan` 项

---

## 5. 调用关系图

```
MoM_SS
 ├── Comp_eltsI
 │    ├── Green(R)          (远场)
 │    ├── F1(R)             (近场正则余项)
 │    └── Spara_cal         (近场 1/R 解析积分)
 │         └── unit_normal
 │         └── Re_cross_product
 │
 └── Comp_eltsU             (仅 CFIE 且 m≠n)
      ├── gradG_S
      │    ├── Spara_cal_2  (近场梯度主值)
      │    │    ├── unit_normal
      │    │    ├── Re_cross_product
      │    │    └── sgn
      │    └── F1 (局部副本)
      └── ReCom_cross_product
```

---

## 6. 关键代码摘录

### 6.1 Green 函数拆分

```fortran
pure complex function Green(R)
    use Common_constants
    implicit none
    real,intent(in):: R
    Green = exp(cmplx(0.0, -k_wnum*R)) / R
endfunction
```

### 6.2 正则余项 `F1`

```fortran
pure complex function F1(R)
    use Common_constants
    implicit none
    real,intent(in):: R
    integer i
    real x
    complex,parameter:: GreenExpan_coeff(0:7) = &
        (/(0,-1.0), (-0.5,0), (0,0.1666667), (4.1666667e-2,0), &
          (0.0,-8.3333333e-3), (-1.3888889e-3,0), (0,1.9841270e-4), &
          (2.4801587e-5,0)/)
    x = k_wnum * R
    F1 = GreenExpan_coeff(7)
    do i=6,0,-1
        F1 = F1 * x
        F1 = F1 + GreenExpan_coeff(i)
    enddo
    F1 = k_wnum * F1
endfunction
```

### 6.3 `Comp_eltsI` 奇异分支

```fortran
do q=1,Nq
    wq = weightS(q)
    rq = ris(:,q)
    R  = re_Norm(rp-rq)
    F1_temp = wq * F1(R)
    inner   = inner   + F1_temp
    innerv  = innerv  + F1_temp * rq
enddo

call Spara_cal(rp, r1s, r2s, r3s, avi, t3, t4, w0, ns)
t3 = t3 / avi
t4 = t4 / avi
innerv = innerv + (rp - ns*w0) * t3

inner  = inner  + t3
innerv = innerv + t4
```

### 6.4 `Comp_eltsU` 奇异分支

```fortran
I_30 = (0.0,0.0)
do q=1,SNum_samp
    rq = ris(:,q)
    R  = rp - rq
    abs_R = re_Norm(R)
    temp_gradient_G = R * cmplx(k_wnum*abs_R*(1.0-k_wnum**2*abs_R**2/18.0)/8.0, &
                                (1.0-0.1*k_wnum**2*abs_R**2)/3.0)
    wq = weightS(q)
    I_30 = I_30 + wq * temp_gradient_G
enddo

call Spara_cal_2(rp, P1, P2, P3, An, t1, t3, t4, t5, w0, ns)

I_31 = 0.0
if(abs(w0) <= ratio_limit1) then
    do i=1,3
        I_31 = I_31 + t5(:,i)
    enddo
else
    do i=1,3
        I_31 = I_31 + t4(:,i) + t5(:,i)
    enddo
endif

I_32 = 0.0
if(abs(w0) <= ratio_limit1) then
    do i=1,3
        I_32 = I_32 - t3(:,i)
    enddo
else
    do i=1,3
        I_32 = I_32 + t1(i)*w0*ns - t3(:,i)
    enddo
endif

gradG = k_wnum**3 * I_30 - (I_31 + k_wnum**2/2.0 * I_32) / An
```

### 6.5 `Spara_cal` 边积分公式

```fortran
do m=1,3
    R0 = P0(m)**2 + w0**2
    Rm = sqrt(lm(m)**2 + P0(m)**2 + w0**2)
    Rp = sqrt(lp(m)**2 + P0(m)**2 + w0**2)

    singular = .FALSE.
    if(lm(m)>0 .and. lp(m)>0) then
        g2 = log((Rp+lp(m))/(Rm+lm(m)))
        g3 = (lp(m)*Rp-lm(m)*Rm) + R0*g2
    elseif(lm(m)<=0 .and. lp(m)>=0) then
        if(Rm+lm(m)==0) then
            singular = .TRUE.
        else
            g2 = log((Rp+lp(m))/(Rm+lm(m)))
            g3 = (lp(m)*Rp-lm(m)*Rm) + R0*g2
        endif
    elseif(lm(m)<0 .and. lp(m)<0) then
        g2 = log((Rm-lm(m))/(Rp-lp(m)))
        g3 = (lp(m)*Rp-lm(m)*Rm) + R0*g2
    endif

    if(R0 <= (0.0e-2*abs_l(m))**2) then
        beta = 0.0
    else
        beta = atan(P0(m)*lp(m)/(R0+abs(w0)*Rp)) &
             - atan(P0(m)*lm(m)/(R0+abs(w0)*Rm))
    endif

    if(.not.singular) then
        t1(m) = P0(m)*g2 - abs(w0)*beta
        t2(m) = 0.5*g3
    else
        t1(m) = -abs(w0)*beta
        t2(m) = 0.5*((lp(m)*Rp-lm(m)*Rm))
    endif
enddo

t3 = sum(t1)
t4 = m1*t2(1) + m2*t2(2) + m3*t2(3)
```

---

## 7. 小结

本求解器的奇异/近奇异处理可概括为：

1. **Green 函数拆分**：`exp(−jkR)/R = 1/R + F1(R)`，其中 `F1(R)` 用 7 阶泰勒多项式数值计算。
2. **解析面积分**：`Spara_cal` 计算 `1/R` 的源三角形解析积分；`Spara_cal_2` 计算梯度主值所需的各阶矩。
3. **局部坐标几何参数化**：将任意三角形映射到 `(u,v,w)` 坐标系，把面积分化为三条边的解析求和。
4. **极限保护**：对 `w0≈0`、`R0≈0`、`Rm+lm≈0` 等退化情形设置阈值分支，避免 `log(0)`、`atan(∞/0)` 等数值问题。
5. **距离阈值控制**：`0.1λ` 作为远/近场切换点，保证远场效率与近场精度。
