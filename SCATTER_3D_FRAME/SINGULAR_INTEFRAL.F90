MODULE SINGULAR_INTEGRAL
    USE EM_TYPES
    IMPLICIT NONE
CONTAINS

    ! 计算两个三维向量的叉乘
    FUNCTION CROSS_PRODUCT(A, B) RESULT(C)
        REAL, INTENT(IN) :: A(3), B(3)
        REAL :: C(3)
        C(1) = A(2)*B(3) - A(3)*B(2)
        C(2) = A(3)*B(1) - A(1)*B(3)
        C(3) = A(1)*B(2) - A(2)*B(1)
    END FUNCTION CROSS_PRODUCT

    ! 对固定场点r，计算源三角形T上的标量势 ∫ 1/|r-r'| dS' 的边界积分
    ! 使用 Wilton 公式: I_S = -Σ H_I * log((R_e+s_e)/(R_s+s_s)) - |d|·Ω
    SUBROUTINE ANALYTIC_SCALAR_POT_1_OVER_R(R_PT, V1, V2, V3, RESULT)
        REAL, INTENT(IN) :: R_PT(3) ! 场点坐标
        REAL, INTENT(IN) :: V1(3), V2(3), V3(3) ! 三角形顶点坐标
        REAL, INTENT(OUT) :: RESULT

        REAL :: NORMAL(3) ! 三角形法向量
        INTEGER :: I ! 循环变量
        REAL :: V_START(3), V_END(3) ! 边的起点和终点
        REAL :: EDGE_VEC(3) ! 边向量
        REAL :: L_I ! 边长度
        REAL :: T_HAT(3) ! 边的单位切向量
        REAL :: M_HAT(3) ! 边的单位法向量
        REAL :: H_I ! 场点到边的垂直距离(有符号距离)
        REAL :: S_START, S_END ! 场点到边的起点和终点的投影长度(有符号)
        REAL :: R_START, R_END ! 场点到边的起点和终点的距离
        REAL :: CONTRIB     ! 当前边的积分贡献
        REAL :: V_OPP(3) ! 当前边的对边顶点坐标
        REAL :: DOT_CHECK ! 验证内法向是否指向内部的点积
        REAL :: H_PT        ! 场点到三角形平面的有符号距离
        REAL :: R1(3), R2(3), R3(3)   ! 场点到三个顶点的向量
        REAL :: R1M, R2M, R3M       ! 场点到三个顶点的距离
        REAL :: SOLID_ANGLE         ! 三角形对场点张的立体角

        NORMAL(1) = (V2(2) - V1(2)) * (V3(3) - V1(3)) - (V2(3) - V1(3)) * (V3(2) - V1(2))
        NORMAL(2) = (V2(3) - V1(3)) * (V3(1) - V1(1)) - (V2(1) - V1(1)) * (V3(3) - V1(3))
        NORMAL(3) = (V2(1) - V1(1)) * (V3(2) - V1(2)) - (V2(2) - V1(2)) * (V3(1) - V1(1))
        NORMAL = NORMAL / SQRT(SUM(NORMAL ** 2)) ! 把法向量归一化

        ! 场点到平面的有符号距离（用于立体角修正）
        H_PT = DOT_PRODUCT(R_PT - V1, NORMAL)

        RESULT = 0.0

        DO I = 1, 3
            SELECT CASE (I)
                CASE (1)
                    V_START = V1
                    V_END = V2
                    V_OPP = V3
                CASE (2)
                    V_START = V2
                    V_END = V3
                    V_OPP = V1
                CASE (3)
                    V_START = V3
                    V_END = V1
                    V_OPP = V2
            END SELECT

            EDGE_VEC = V_END - V_START
            L_I = SQRT(SUM(EDGE_VEC ** 2))
            T_HAT = EDGE_VEC / L_I ! 计算边的单位切向量

            ! 计算边的内法向 M_HAT = NORMAL x T_HAT
            M_HAT(1) = NORMAL(2) * T_HAT(3) - NORMAL(3) * T_HAT(2)
            M_HAT(2) = NORMAL(3) * T_HAT(1) - NORMAL(1) * T_HAT(3)
            M_HAT(3) = NORMAL(1) * T_HAT(2) - NORMAL(2) * T_HAT(1)

            ! 验证内法向是否指向内部
            DOT_CHECK = DOT_PRODUCT(M_HAT, V_OPP - V_START)
            IF (DOT_CHECK < 0.0) THEN
                M_HAT = -M_HAT
            END IF

            ! 计算场点到边的垂直距离 H_I
            H_I = DOT_PRODUCT(R_PT - V_START, M_HAT)
            S_START = DOT_PRODUCT(R_PT - V_START, T_HAT)
            S_END = DOT_PRODUCT(R_PT - V_END, T_HAT)
            R_START = SQRT(SUM((R_PT - V_START) ** 2))
            R_END = SQRT(SUM((R_PT - V_END) ** 2))

            ! 计算当前边的积分贡献
            IF (ABS(H_I) < 1.0E-12) THEN
                ! 场点在边的延长线上，贡献为0
                CONTRIB = 0.0
            ELSE IF (R_START + S_START > 1.0E-12 .AND. R_END + S_END > 1.0E-12) THEN
                ! 正常情况：使用标准对数公式
                CONTRIB = H_I * LOG((R_END + S_END) / (R_START + S_START))
            ELSE
                ! 近奇异性情况：场点靠近某顶点，使用恒等式 R²-S² = ρ² = h²+H_I² 进行稳定计算
                ! log(R±S) = log(ρ²) - log(R∓S)
                ! ρ² = H_PT² + H_I²（3D垂直距离平方，与线性势一致）
                IF (R_START + S_START <= 1.0E-12) THEN
                    ! 靠近起点：R_START+S_START ≈ ρ²/(R_START-S_START)
                    CONTRIB = H_I * (LOG(R_END + S_END) &
                        + LOG(R_START - S_START) - LOG(MAX(H_PT**2 + H_I**2, 1.0E-12)))
                ELSE
                    ! 靠近终点：R_END+S_END ≈ ρ²/(R_END-S_END)
                    CONTRIB = H_I * (LOG(MAX(H_PT**2 + H_I**2, 1.0E-12)) &
                        - LOG(R_END - S_END) - LOG(R_START + S_START))
                END IF
            END IF

            RESULT = RESULT + CONTRIB
        END DO

        ! 计算立体角 Ω（Van Oosterom & Strackee 公式）
        ! Ω = 2 * atan2( |r1·(r2×r3)|, r1r2r3 + r1(r2·r3) + r2(r3·r1) + r3(r1·r2) )
        R1 = R_PT - V1;  R2 = R_PT - V2;  R3 = R_PT - V3
        R1M = SQRT(SUM(R1**2)); R2M = SQRT(SUM(R2**2)); R3M = SQRT(SUM(R3**2))
        SOLID_ANGLE = 2.0 * ATAN2( ABS(DOT_PRODUCT(R1, CROSS_PRODUCT(R2, R3))),  &
            R1M*R2M*R3M + R1M*DOT_PRODUCT(R2,R3) + R2M*DOT_PRODUCT(R3,R1) + R3M*DOT_PRODUCT(R1,R2) )

        ! Wilton公式符号约定：积分值需取负号使 ∫ 1/R dS' > 0
        ! 完整公式: I₀ = -Σ edge_contrib - |H_PT| * Ω
        RESULT = -RESULT - ABS(H_PT) * SOLID_ANGLE

    END SUBROUTINE ANALYTIC_SCALAR_POT_1_OVER_R

    ! 计算三个线性势 I_i = ∫_T λ_i / |r-r'| dS'，其中 λ_i 为重心坐标
    ! 方法：先通过边界积分求矢量势 ∫_T r'/R dS'，再解线性方程组得到 I_i。
    ! 此实现保证 Σ I_i = I_0（标量势），消除了原实现中的 1/2 因子错误。
    SUBROUTINE ANALYTIC_LINEAR_POT_1_OVER_R(R_PT, V1, V2, V3, RESULT_3)
        REAL, INTENT(IN) :: R_PT(3) ! 场点坐标
        REAL, INTENT(IN) :: V1(3), V2(3), V3(3) ! 三角形顶点坐标
        REAL, INTENT(OUT) :: RESULT_3(3)

        REAL :: I_S            ! 标量势 ∫ 1/R dS'
        REAL :: NORMAL(3)      ! 单位法向量
        REAL :: AREA           ! 三角形面积
        REAL :: H              ! 场点到三角形平面的有符号距离
        REAL :: P_PROJ(3)      ! 场点在三角形平面上的投影
        REAL :: I_SHIFTED(3)   ! ∫_T (r' - P_PROJ)/R dS'
        REAL :: I_VEC(3)       ! 矢量势 ∫_T r'/R dS'
        INTEGER :: I
        REAL :: V_START(3), V_END(3), V_OPP(3)
        REAL :: EDGE_VEC(3), T_HAT(3), M_HAT(3), M_OUT(3)
        REAL :: L_I            ! 边长
        REAL :: S_START, S_END ! 投影长度
        REAL :: R_START, R_END ! 到场点距离
        REAL :: D_I            ! 投影点到边的平面内距离（带符号）
        REAL :: RHO2           ! h^2 + d^2
        REAL :: LOG_TERM, TERM
        REAL :: H2_D2
        REAL, PARAMETER :: EPS = 1.0E-12

        ! 法向量和面积
        NORMAL = CROSS_PRODUCT(V2 - V1, V3 - V1)
        AREA = 0.5 * SQRT(SUM(NORMAL**2))
        NORMAL = NORMAL / SQRT(SUM(NORMAL**2))

        ! 标量势
        CALL ANALYTIC_SCALAR_POT_1_OVER_R(R_PT, V1, V2, V3, I_S)

        ! 场点在三角形平面上的投影
        H = DOT_PRODUCT(R_PT - V1, NORMAL)
        P_PROJ = R_PT - H * NORMAL

        ! 边界积分计算 I_SHIFTED = ∫_T (r' - P_PROJ)/R dS'
        I_SHIFTED = 0.0
        DO I = 1, 3
            SELECT CASE (I)
                CASE (1)
                    V_START = V1; V_END = V2; V_OPP = V3
                CASE (2)
                    V_START = V2; V_END = V3; V_OPP = V1
                CASE (3)
                    V_START = V3; V_END = V1; V_OPP = V2
            END SELECT

            EDGE_VEC = V_END - V_START
            L_I = SQRT(SUM(EDGE_VEC**2))
            T_HAT = EDGE_VEC / L_I
            M_HAT = CROSS_PRODUCT(NORMAL, T_HAT)

            ! M_HAT = n x t。验证并调整使其指向三角形内部（与标量势约定一致）
            IF (DOT_PRODUCT(M_HAT, V_OPP - V_START) < 0.0) THEN
                M_HAT = -M_HAT
            END IF
            ! 边界积分公式需要边的外法向
            M_OUT = -M_HAT

            S_START = DOT_PRODUCT(P_PROJ - V_START, T_HAT)
            S_END   = DOT_PRODUCT(P_PROJ - V_END,   T_HAT)
            R_START = SQRT(SUM((R_PT - V_START)**2))
            R_END   = SQRT(SUM((R_PT - V_END)**2))
            D_I     = DOT_PRODUCT(P_PROJ - V_START, M_OUT)
            RHO2    = H**2 + D_I**2

            ! 稳定计算对数项
            IF (R_START + S_START > EPS .AND. R_END + S_END > EPS) THEN
                LOG_TERM = LOG((R_START + S_START)/(R_END + S_END))
            ELSE IF (R_START + S_START <= EPS) THEN
                H2_D2 = MAX(RHO2, EPS)
                LOG_TERM = LOG(H2_D2) - LOG(MAX(R_START - S_START, EPS)) - LOG(R_END + S_END)
            ELSE
                H2_D2 = MAX(RHO2, EPS)
                LOG_TERM = LOG(R_START + S_START) + LOG(MAX(R_END - S_END, EPS)) - LOG(H2_D2)
            END IF

            TERM = RHO2 * LOG_TERM + S_START * R_START - S_END * R_END
            I_SHIFTED = I_SHIFTED + 0.5 * M_OUT * TERM
        END DO

        ! 矢量势
        I_VEC = P_PROJ * I_S + I_SHIFTED

        ! 解 I_VEC = V1*I1 + V2*I2 + V3*I3，且 I1+I2+I3 = I_S
        RESULT_3(1) = DOT_PRODUCT(CROSS_PRODUCT(I_VEC - V3*I_S, V2 - V3), NORMAL) / (2.0 * AREA)
        RESULT_3(2) = DOT_PRODUCT(CROSS_PRODUCT(V1 - V3, I_VEC - V3*I_S), NORMAL) / (2.0 * AREA)
        RESULT_3(3) = I_S - RESULT_3(1) - RESULT_3(2)

    END SUBROUTINE ANALYTIC_LINEAR_POT_1_OVER_R

    ! =========================================================================
    ! MFIE 静态主值（PV）矩：g_A^st, g_C^st, g_E^st
    ! 参考 docs/MFIE_CFIE实现文档_内谐振对策.md §4.2
    ! 输入：R_PT 场点，N_F 场三角形单位法向，V1/V2/V3 源三角形顶点
    ! 输出：G_A_ST, G_C_ST, G_E_ST（均不含 1/(4π) 因子）
    ! =========================================================================
    SUBROUTINE ANALYTIC_GRAD_STATIC_PV(R_PT, N_F, V1, V2, V3, &
                                        G_A_ST, G_C_ST, G_E_ST)
        REAL, INTENT(IN)  :: R_PT(3), N_F(3), V1(3), V2(3), V3(3)
        REAL, INTENT(OUT) :: G_A_ST(3), G_C_ST(3), G_E_ST(3)

        REAL :: N_S(3), P_PROJ(3), R1(3), R2(3), R3(3)
        REAL :: D, H_PT, OMEGA, I0
        REAL :: N_F_PAR(3), N_F_DOT_N_S, N_S_DOT_P, SIGMA_D
        REAL :: SUM_M_G2(3), SUM_L(3), SUM_NFP_M_G2, SUM_NFP_L(3)
        REAL :: TERM_C(3)
        REAL :: N_F_DOT_P
        INTEGER :: I
        REAL :: V_START(3), V_END(3), V_OPP(3), EDGE_VEC(3)
        REAL :: T_HAT(3), M_HAT(3), M_OUT(3)
        REAL :: L_E, S_START, S_END, R_START, R_END, G2, H_I
        REAL :: L_VEC(3)
        REAL, PARAMETER :: EPS = 1.0E-12

        ! 源三角形法向与面积
        N_S = CROSS_PRODUCT(V2 - V1, V3 - V1)
        N_S = N_S / SQRT(SUM(N_S**2))

        ! 场点到源平面的有符号距离
        D = DOT_PRODUCT(R_PT - V1, N_S)
        P_PROJ = R_PT - D * N_S
        H_PT = D

        ! 标量势 I0 = ∫ 1/R dS'
        CALL ANALYTIC_SCALAR_POT_1_OVER_R(R_PT, V1, V2, V3, I0)

        ! 立体角 |Ω|（Van Oosterom & Strackee）
        R1 = R_PT - V1; R2 = R_PT - V2; R3 = R_PT - V3
        OMEGA = 2.0 * ATAN2(ABS(DOT_PRODUCT(R1, CROSS_PRODUCT(R2, R3))), &
                            SQRT(SUM(R1**2))*SQRT(SUM(R2**2))*SQRT(SUM(R3**2)) + &
                            SQRT(SUM(R1**2))*DOT_PRODUCT(R2, R3) + &
                            SQRT(SUM(R2**2))*DOT_PRODUCT(R3, R1) + &
                            SQRT(SUM(R3**2))*DOT_PRODUCT(R1, R2))

        ! 场法向的源面内分量
        N_F_DOT_N_S = DOT_PRODUCT(N_F, N_S)
        N_F_PAR = N_F - N_F_DOT_N_S * N_S
        N_S_DOT_P = DOT_PRODUCT(N_S, P_PROJ)

        ! σ(d)：符号函数，d≈0 时取 0
        IF (ABS(D) < EPS) THEN
            SIGMA_D = 0.0
        ELSE
            SIGMA_D = SIGN(1.0, D)
        END IF

        ! 边线积分累加
        SUM_M_G2 = 0.0
        SUM_L = 0.0
        SUM_NFP_M_G2 = 0.0
        SUM_NFP_L = 0.0

        DO I = 1, 3
            SELECT CASE (I)
                CASE (1); V_START = V1; V_END = V2; V_OPP = V3
                CASE (2); V_START = V2; V_END = V3; V_OPP = V1
                CASE (3); V_START = V3; V_END = V1; V_OPP = V2
            END SELECT

            EDGE_VEC = V_END - V_START
            L_E = SQRT(SUM(EDGE_VEC**2))
            T_HAT = EDGE_VEC / L_E

            ! 内法向
            M_HAT = CROSS_PRODUCT(N_S, T_HAT)
            IF (DOT_PRODUCT(M_HAT, V_OPP - V_START) < 0.0) M_HAT = -M_HAT
            ! 外法向
            M_OUT = -M_HAT

            ! 投影量
            S_START = DOT_PRODUCT(P_PROJ - V_START, T_HAT)
            S_END   = DOT_PRODUCT(P_PROJ - V_END,   T_HAT)
            R_START = SQRT(SUM((R_PT - V_START)**2))
            R_END   = SQRT(SUM((R_PT - V_END)**2))
            H_I     = DOT_PRODUCT(R_PT - V_START, M_HAT)

            ! g_{2,e} = ∫_e dl/R = ln((R_start+s_start)/(R_end+s_end))，带稳定化
            IF (R_START + S_START > EPS .AND. R_END + S_END > EPS) THEN
                G2 = LOG((R_START + S_START) / (R_END + S_END))
            ELSE IF (R_START + S_START <= EPS) THEN
                G2 = LOG(MAX(H_PT**2 + H_I**2, EPS)) - LOG(MAX(R_START - S_START, EPS)) - &
                     LOG(R_END + S_END)
            ELSE
                G2 = LOG(R_START + S_START) + LOG(MAX(R_END - S_END, EPS)) - &
                     LOG(MAX(H_PT**2 + H_I**2, EPS))
            END IF

            ! L_e = v_start g2 + t_hat [(R_end - R_start) + s_start g2]
            L_VEC = V_START * G2 + T_HAT * ((R_END - R_START) + S_START * G2)

            SUM_M_G2 = SUM_M_G2 + M_OUT * G2
            SUM_L    = SUM_L    + L_VEC
            SUM_NFP_M_G2 = SUM_NFP_M_G2 + DOT_PRODUCT(N_F_PAR, M_OUT) * G2
            SUM_NFP_L    = SUM_NFP_L    + DOT_PRODUCT(N_F_PAR, M_OUT) * L_VEC
        END DO

        ! 退化检测：d≈0 且 n_f_parallel≈0（自作用、共面相邻对）
        IF (ABS(D) < EPS .AND. SQRT(SUM(N_F_PAR**2)) < 1.0E-10) THEN
            G_A_ST = -SUM_M_G2
            G_C_ST = 0.0
            G_E_ST = N_F_DOT_N_S * N_S_DOT_P * G_A_ST
            RETURN
        END IF

        ! g_A^st = -Σ m̂_out g2 - n_s σ(d)|Ω|
        G_A_ST = -SUM_M_G2 - N_S * SIGMA_D * OMEGA

        ! g_C^st = (n_f·n_s)[ -σ(d)|Ω| p + d Σ m̂_out g2 ] - Σ(n_f_par·m̂_out)L_e + n_f_par I0
        TERM_C = -SIGMA_D * OMEGA * P_PROJ + D * SUM_M_G2
        G_C_ST = N_F_DOT_N_S * TERM_C - SUM_NFP_L + N_F_PAR * I0

        ! g_E^st = I0 n_f^∥ - Σ(n_f^∥·m)L + (n_f·p) g_A^st + r (n_f^∥·Σm g2)
        ! 参考文档 §4.2 并修正其首项为 (n_f·p) g_A（非 (n_f·n_s)(n_s·p) g_A）
        N_F_DOT_P = DOT_PRODUCT(N_F, P_PROJ)
        G_E_ST = N_F_PAR * I0 - SUM_NFP_L + N_F_DOT_P * G_A_ST + &
                 R_PT * DOT_PRODUCT(N_F_PAR, SUM_M_G2)

    END SUBROUTINE ANALYTIC_GRAD_STATIC_PV


END MODULE SINGULAR_INTEGRAL
