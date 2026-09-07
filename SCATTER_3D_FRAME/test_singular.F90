! 奇异积分公式验证程序
! 对应 docs/奇异积分处理文档.md 中的测试项目
PROGRAM TEST_SINGULAR
    USE EM_TYPES
    USE NUMERICAL_INTEGRATION
    USE GREEN_FUNCTIONS
    USE SINGULAR_INTEGRAL
    IMPLICIT NONE

    REAL :: V1(3), V2(3), V3(3)
    REAL :: R_PT(3)
    REAL :: I0_ANA, I0_REF
    REAL :: C_ANA(3)
    REAL :: IVEC_ANA(3), IVEC_REF(3)
    REAL :: ERR, MAX_ERR
    LOGICAL :: ALL_PASS

    ALL_PASS = .TRUE.

    ! 测试三角形：边长为 1 的正三角形，位于 z=0 平面
    V1 = [0.0, 0.0, 0.0]
    V2 = [1.0, 0.0, 0.0]
    V3 = [0.5, SQRT(3.0)/2.0, 0.0]

    PRINT *, '============================================================'
    PRINT *, '  奇异积分公式验证'
    PRINT *, '============================================================'

    ! ----------------------------------------------------------------
    ! Test 1: 标量势 I0 解析 vs 高密度数值（场点在三角形上方）
    ! ----------------------------------------------------------------
    R_PT = [0.5, SQRT(3.0)/6.0, 0.5]   ! 重心正上方 0.5
    CALL ANALYTIC_SCALAR_POT_1_OVER_R(R_PT, V1, V2, V3, I0_ANA)
    CALL NUMERIC_SCALAR_POT_1_OVER_R(R_PT, V1, V2, V3, I0_REF)
    ERR = ABS(I0_ANA - I0_REF) / MAX(ABS(I0_REF), 1.0E-12)
    PRINT *, 'Test 1: I0 解析 vs 数值（面外点）'
    PRINT *, '  解析值 = ', I0_ANA, ' 数值值 = ', I0_REF
    PRINT *, '  相对误差 = ', ERR
    IF (ERR < 0.01) THEN
        PRINT *, '  PASSED'
    ELSE
        PRINT *, '  FAILED'
        ALL_PASS = .FALSE.
    END IF

    ! ----------------------------------------------------------------
    ! Test 2: 标量势 I0 解析 vs 高密度数值（场点在三角形平面内但在外部）
    ! 该点不在源三角形上，因此 1/R 数值参考积分可靠，同时检验面内解析公式。
    ! ----------------------------------------------------------------
    R_PT = [1.5, 0.5, 0.0]
    CALL ANALYTIC_SCALAR_POT_1_OVER_R(R_PT, V1, V2, V3, I0_ANA)
    CALL NUMERIC_SCALAR_POT_1_OVER_R(R_PT, V1, V2, V3, I0_REF)
    ERR = ABS(I0_ANA - I0_REF) / MAX(ABS(I0_REF), 1.0E-12)
    PRINT *, 'Test 2: I0 解析 vs 数值（面外点，平面外）'
    PRINT *, '  解析值 = ', I0_ANA, ' 数值值 = ', I0_REF
    PRINT *, '  相对误差 = ', ERR
    IF (ERR < 0.01) THEN
        PRINT *, '  PASSED'
    ELSE
        PRINT *, '  FAILED'
        ALL_PASS = .FALSE.
    END IF

    ! ----------------------------------------------------------------
    ! Test 3: 线性势归一化恒等式 Σ c_i = I0
    ! ----------------------------------------------------------------
    R_PT = [0.7, 0.2, 0.3]
    CALL ANALYTIC_SCALAR_POT_1_OVER_R(R_PT, V1, V2, V3, I0_ANA)
    CALL ANALYTIC_LINEAR_POT_1_OVER_R(R_PT, V1, V2, V3, C_ANA)
    ERR = ABS(SUM(C_ANA) - I0_ANA) / MAX(ABS(I0_ANA), 1.0E-12)
    PRINT *, 'Test 3: Σ c_i = I0 恒等式'
    PRINT *, '  I0 = ', I0_ANA, ' Σc_i = ', SUM(C_ANA)
    PRINT *, '  相对误差 = ', ERR
    IF (ERR < 1.0E-10) THEN
        PRINT *, '  PASSED'
    ELSE
        PRINT *, '  FAILED'
        ALL_PASS = .FALSE.
    END IF

    ! ----------------------------------------------------------------
    ! Test 4: 矢量势 I_vec 解析 vs 数值
    ! ----------------------------------------------------------------
    IVEC_ANA = C_ANA(1)*V1 + C_ANA(2)*V2 + C_ANA(3)*V3
    CALL NUMERIC_VECTOR_POT_1_OVER_R(R_PT, V1, V2, V3, IVEC_REF)
    MAX_ERR = MAXVAL(ABS(IVEC_ANA - IVEC_REF)) / MAX(MAXVAL(ABS(IVEC_REF)), 1.0E-12)
    PRINT *, 'Test 4: I_vec 解析 vs 数值'
    PRINT *, '  解析值 = ', IVEC_ANA
    PRINT *, '  数值值 = ', IVEC_REF
    PRINT *, '  最大相对误差 = ', MAX_ERR
    IF (MAX_ERR < 0.01) THEN
        PRINT *, '  PASSED'
    ELSE
        PRINT *, '  FAILED'
        ALL_PASS = .FALSE.
    END IF

    ! ----------------------------------------------------------------
    ! Test 5: 质心点对称性 c2 = c3（正三角形）
    ! ----------------------------------------------------------------
    R_PT = [0.5, SQRT(3.0)/6.0, 0.0]   ! 重心
    CALL ANALYTIC_LINEAR_POT_1_OVER_R(R_PT, V1, V2, V3, C_ANA)
    ERR = ABS(C_ANA(2) - C_ANA(3))
    PRINT *, 'Test 5: 质心点 c2 = c3 对称性'
    PRINT *, '  c1 = ', C_ANA(1), ' c2 = ', C_ANA(2), ' c3 = ', C_ANA(3)
    PRINT *, '  |c2 - c3| = ', ERR
    IF (ERR < 1.0E-6) THEN
        PRINT *, '  PASSED'
    ELSE
        PRINT *, '  FAILED'
        ALL_PASS = .FALSE.
    END IF

    ! ----------------------------------------------------------------
    ! Test 6: G_smooth 恒等式 G_smooth = G - 1/(4πR)
    ! ----------------------------------------------------------------
    CALL TEST_GREEN_SMOOTH_IDENTITY(ALL_PASS)

    ! ----------------------------------------------------------------
    ! Test 7: G_smooth(0) 极限 = -jk/(4π)
    ! ----------------------------------------------------------------
    CALL TEST_GREEN_SMOOTH_LIMIT(ALL_PASS)

    ! ----------------------------------------------------------------
    ! Test 8: J0 收敛性（3/7/12 点外层高斯）
    ! ----------------------------------------------------------------
    CALL TEST_J0_CONVERGENCE(V1, V2, V3, ALL_PASS)

    ! ----------------------------------------------------------------
    ! Test 9: MFIE 静态主值退化情形（g_C^st = 0, g_E^st = 0）
    ! ----------------------------------------------------------------
    CALL TEST_MFIE_GRAD_STATIC_DEGENERATE(V1, V2, V3, ALL_PASS)

    ! ----------------------------------------------------------------
    ! Test 10: MFIE 静态主值 + 正则余项 ≡ 完整 ∇G（面外点）
    ! ----------------------------------------------------------------
    CALL TEST_MFIE_GRAD_CONSISTENCY(V1, V2, V3, ALL_PASS)

    ! ----------------------------------------------------------------
    ! 汇总
    ! ----------------------------------------------------------------
    PRINT *, '============================================================'
    IF (ALL_PASS) THEN
        PRINT *, '  全部测试通过'
    ELSE
        PRINT *, '  存在失败的测试'
        STOP 1
    END IF
    PRINT *, '============================================================'

CONTAINS

    ! 将三角形 1→4 等面积细分（连接三边中点）
    SUBROUTINE SUBDIVIDE_TRI_4(VA, VB, VC, SUB_V)
        REAL, INTENT(IN)  :: VA(3), VB(3), VC(3)
        REAL, INTENT(OUT) :: SUB_V(3, 3, 4)
        REAL :: MAB(3), MBC(3), MCA(3)
        INTEGER :: K
        DO K = 1, 3
            MAB(K) = (VA(K) + VB(K)) / 2.0
            MBC(K) = (VB(K) + VC(K)) / 2.0
            MCA(K) = (VC(K) + VA(K)) / 2.0
        END DO
        SUB_V(:, 1, 1) = VA;  SUB_V(:, 2, 1) = MAB; SUB_V(:, 3, 1) = MCA
        SUB_V(:, 1, 2) = VB;  SUB_V(:, 2, 2) = MBC; SUB_V(:, 3, 2) = MAB
        SUB_V(:, 1, 3) = VC;  SUB_V(:, 2, 3) = MCA; SUB_V(:, 3, 3) = MBC
        SUB_V(:, 1, 4) = MAB; SUB_V(:, 2, 4) = MBC; SUB_V(:, 3, 4) = MCA
    END SUBROUTINE SUBDIVIDE_TRI_4

    ! 递归地用 12 点高斯求 1/R 的标量势（场点不在源三角形上）
    RECURSIVE SUBROUTINE NUMERIC_SCALAR_POT_RECURSIVE(VA, VB, VC, R_PT, &
                                                       LEVEL, MAX_LEVEL, GDATA, ACC)
        REAL, INTENT(IN)    :: VA(3), VB(3), VC(3), R_PT(3)
        INTEGER, INTENT(IN) :: LEVEL, MAX_LEVEL
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        REAL, INTENT(INOUT) :: ACC
        REAL :: SUB_V(3, 3, 4), AREA_SUB, CROSS(3), R_VEC(3)
        INTEGER :: I, K
        IF (LEVEL >= MAX_LEVEL) THEN
            CROSS(1) = (VB(2)-VA(2))*(VC(3)-VA(3)) - (VB(3)-VA(3))*(VC(2)-VA(2))
            CROSS(2) = (VB(3)-VA(3))*(VC(1)-VA(1)) - (VB(1)-VA(1))*(VC(3)-VA(3))
            CROSS(3) = (VB(1)-VA(1))*(VC(2)-VA(2)) - (VB(2)-VA(2))*(VC(1)-VA(1))
            AREA_SUB = 0.5 * SQRT(SUM(CROSS**2))
            DO I = 1, GDATA%N_POINTS
                DO K = 1, 3
                    R_VEC(K) = VA(K)*GDATA%UVW(1,I) + VB(K)*GDATA%UVW(2,I) + VC(K)*GDATA%UVW(3,I)
                END DO
                ACC = ACC + GDATA%WEIGHTS(I) * AREA_SUB / SQRT(SUM((R_PT - R_VEC)**2))
            END DO
        ELSE
            CALL SUBDIVIDE_TRI_4(VA, VB, VC, SUB_V)
            DO I = 1, 4
                CALL NUMERIC_SCALAR_POT_RECURSIVE(SUB_V(:,1,I), SUB_V(:,2,I), SUB_V(:,3,I), &
                                                   R_PT, LEVEL+1, MAX_LEVEL, GDATA, ACC)
            END DO
        END IF
    END SUBROUTINE NUMERIC_SCALAR_POT_RECURSIVE

    SUBROUTINE NUMERIC_SCALAR_POT_1_OVER_R(R_PT, VA, VB, VC, VAL)
        REAL, INTENT(IN)  :: R_PT(3), VA(3), VB(3), VC(3)
        REAL, INTENT(OUT) :: VAL
        TYPE(GAUSS_TRI_DATA) :: GDATA
        VAL = 0.0
        CALL INIT_GAUSS_TRI(GAUSS_12PT, GDATA)
        CALL NUMERIC_SCALAR_POT_RECURSIVE(VA, VB, VC, R_PT, 0, 2, GDATA, VAL)
    END SUBROUTINE NUMERIC_SCALAR_POT_1_OVER_R

    ! 递归地用 12 点高斯求矢量势 ∫ r'/R dS'
    RECURSIVE SUBROUTINE NUMERIC_VECTOR_POT_RECURSIVE(VA, VB, VC, R_PT, &
                                                       LEVEL, MAX_LEVEL, GDATA, ACC)
        REAL, INTENT(IN)    :: VA(3), VB(3), VC(3), R_PT(3)
        INTEGER, INTENT(IN) :: LEVEL, MAX_LEVEL
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        REAL, INTENT(INOUT) :: ACC(3)
        REAL :: SUB_V(3, 3, 4), AREA_SUB, CROSS(3), R_VEC(3)
        INTEGER :: I, K
        IF (LEVEL >= MAX_LEVEL) THEN
            CROSS(1) = (VB(2)-VA(2))*(VC(3)-VA(3)) - (VB(3)-VA(3))*(VC(2)-VA(2))
            CROSS(2) = (VB(3)-VA(3))*(VC(1)-VA(1)) - (VB(1)-VA(1))*(VC(3)-VA(3))
            CROSS(3) = (VB(1)-VA(1))*(VC(2)-VA(2)) - (VB(2)-VA(2))*(VC(1)-VA(1))
            AREA_SUB = 0.5 * SQRT(SUM(CROSS**2))
            DO I = 1, GDATA%N_POINTS
                DO K = 1, 3
                    R_VEC(K) = VA(K)*GDATA%UVW(1,I) + VB(K)*GDATA%UVW(2,I) + VC(K)*GDATA%UVW(3,I)
                END DO
                ACC = ACC + GDATA%WEIGHTS(I) * AREA_SUB * R_VEC / SQRT(SUM((R_PT - R_VEC)**2))
            END DO
        ELSE
            CALL SUBDIVIDE_TRI_4(VA, VB, VC, SUB_V)
            DO I = 1, 4
                CALL NUMERIC_VECTOR_POT_RECURSIVE(SUB_V(:,1,I), SUB_V(:,2,I), SUB_V(:,3,I), &
                                                   R_PT, LEVEL+1, MAX_LEVEL, GDATA, ACC)
            END DO
        END IF
    END SUBROUTINE NUMERIC_VECTOR_POT_RECURSIVE

    SUBROUTINE NUMERIC_VECTOR_POT_1_OVER_R(R_PT, VA, VB, VC, VAL)
        REAL, INTENT(IN)  :: R_PT(3), VA(3), VB(3), VC(3)
        REAL, INTENT(OUT) :: VAL(3)
        TYPE(GAUSS_TRI_DATA) :: GDATA
        VAL = 0.0
        CALL INIT_GAUSS_TRI(GAUSS_12PT, GDATA)
        CALL NUMERIC_VECTOR_POT_RECURSIVE(VA, VB, VC, R_PT, 0, 2, GDATA, VAL)
    END SUBROUTINE NUMERIC_VECTOR_POT_1_OVER_R

    SUBROUTINE TEST_GREEN_SMOOTH_IDENTITY(ALL_PASS)
        LOGICAL, INTENT(INOUT) :: ALL_PASS
        REAL :: R1(3), R2(3)
        COMPLEX :: G, G_SMOOTH, G_DIRECT
        REAL :: K, ERR
        INTEGER :: I
        K = 2.0 * PI   ! 波数任意取
        R1 = [0.1, 0.2, 0.3]
        R2 = [0.4, 0.1, 0.5]
        PRINT *, 'Test 6: G_smooth = G - 1/(4πR) 恒等式'
        DO I = 1, 5
            R2 = R2 + [0.01*I, 0.0, 0.0]
            G = GREEN_FUNC(R1, R2, K)
            G_SMOOTH = GREEN_FUNC_SMOOTH(R1, R2, K)
            G_DIRECT = G - 1.0 / (4.0 * PI * SQRT(SUM((R1 - R2)**2)))
            ERR = ABS(G_SMOOTH - G_DIRECT)
            PRINT *, '  点 ', I, ' 误差 = ', ERR
            IF (ERR > 1.0E-7) THEN
                PRINT *, '  FAILED'
                ALL_PASS = .FALSE.
                RETURN
            END IF
        END DO
        PRINT *, '  PASSED'
    END SUBROUTINE TEST_GREEN_SMOOTH_IDENTITY

    SUBROUTINE TEST_GREEN_SMOOTH_LIMIT(ALL_PASS)
        LOGICAL, INTENT(INOUT) :: ALL_PASS
        REAL :: R1(3), R2(3), K
        COMPLEX :: G_SMOOTH, LIMIT
        REAL :: ERR
        R1 = [0.0, 0.0, 0.0]
        R2 = R1
        K = 2.0 * PI
        G_SMOOTH = GREEN_FUNC_SMOOTH(R1, R2, K)
        LIMIT = -(0.0, 1.0) * K / (4.0 * PI)
        ERR = ABS(G_SMOOTH - LIMIT)
        PRINT *, 'Test 7: G_smooth(0) 极限 = -jk/(4π)'
        PRINT *, '  G_smooth(0) = ', G_SMOOTH, ' 极限 = ', LIMIT
        PRINT *, '  误差 = ', ERR
        IF (ERR < 1.0E-8) THEN
            PRINT *, '  PASSED'
        ELSE
            PRINT *, '  FAILED'
            ALL_PASS = .FALSE.
        END IF
    END SUBROUTINE TEST_GREEN_SMOOTH_LIMIT

    SUBROUTINE TEST_J0_CONVERGENCE(V1, V2, V3, ALL_PASS)
        REAL, INTENT(IN) :: V1(3), V2(3), V3(3)
        LOGICAL, INTENT(INOUT) :: ALL_PASS
        REAL :: J0_3, J0_7, J0_12
        REAL :: ERR_7, ERR_12

        CALL COMPUTE_J0(V1, V2, V3, GAUSS_3PT,  J0_3)
        CALL COMPUTE_J0(V1, V2, V3, GAUSS_7PT,  J0_7)
        CALL COMPUTE_J0(V1, V2, V3, GAUSS_12PT, J0_12)

        ERR_7  = ABS(J0_7  - J0_12) / MAX(ABS(J0_12), 1.0E-12)
        ERR_12 = ABS(J0_3  - J0_12) / MAX(ABS(J0_12), 1.0E-12)

        PRINT *, 'Test 8: J0 = ∬_T 1/R dS dS'' 收敛性'
        PRINT *, '  J0_3pt  = ', J0_3
        PRINT *, '  J0_7pt  = ', J0_7
        PRINT *, '  J0_12pt = ', J0_12
        PRINT *, '  |J0_7 - J0_12| / |J0_12| = ', ERR_7
        PRINT *, '  |J0_3 - J0_12| / |J0_12| = ', ERR_12
        IF (ERR_7 < 0.01 .AND. ERR_12 < 0.05) THEN
            PRINT *, '  PASSED'
        ELSE
            PRINT *, '  FAILED'
            ALL_PASS = .FALSE.
        END IF
    END SUBROUTINE TEST_J0_CONVERGENCE

    SUBROUTINE COMPUTE_J0(V1, V2, V3, ORDER, J0)
        REAL, INTENT(IN)  :: V1(3), V2(3), V3(3)
        INTEGER, INTENT(IN) :: ORDER
        REAL, INTENT(OUT) :: J0
        TYPE(GAUSS_TRI_DATA) :: GDATA
        REAL :: R_VEC(3), I0_VAL
        REAL :: AREA, CROSS(3)
        INTEGER :: I, K

        CROSS(1) = (V2(2)-V1(2))*(V3(3)-V1(3)) - (V2(3)-V1(3))*(V3(2)-V1(2))
        CROSS(2) = (V2(3)-V1(3))*(V3(1)-V1(1)) - (V2(1)-V1(1))*(V3(3)-V1(3))
        CROSS(3) = (V2(1)-V1(1))*(V3(2)-V1(2)) - (V2(2)-V1(2))*(V3(1)-V1(1))
        AREA = 0.5 * SQRT(SUM(CROSS**2))

        CALL INIT_GAUSS_TRI(ORDER, GDATA)
        J0 = 0.0
        DO I = 1, GDATA%N_POINTS
            DO K = 1, 3
                R_VEC(K) = V1(K)*GDATA%UVW(1,I) + V2(K)*GDATA%UVW(2,I) + V3(K)*GDATA%UVW(3,I)
            END DO
            CALL ANALYTIC_SCALAR_POT_1_OVER_R(R_VEC, V1, V2, V3, I0_VAL)
            J0 = J0 + GDATA%WEIGHTS(I) * I0_VAL * AREA
        END DO
    END SUBROUTINE COMPUTE_J0

    SUBROUTINE TEST_MFIE_GRAD_STATIC_DEGENERATE(V1, V2, V3, ALL_PASS)
        REAL, INTENT(IN) :: V1(3), V2(3), V3(3)
        LOGICAL, INTENT(INOUT) :: ALL_PASS
        REAL :: N_F(3), R_PT(3)
        REAL :: G_A_ST(3), G_C_ST(3), G_E_ST(3)
        REAL :: ERR_C, ERR_E, ERR_B
        REAL, PARAMETER :: TOL = 1.0E-6

        ! 源三角形法向
        N_F = CROSS_PRODUCT(V2 - V1, V3 - V1)
        N_F = N_F / SQRT(SUM(N_F**2))

        ! 退化 1：场点在三角形平面内（重心），n_f = n_s
        R_PT = [0.5, SQRT(3.0)/6.0, 0.0]
        CALL ANALYTIC_GRAD_STATIC_PV(R_PT, N_F, V1, V2, V3, G_A_ST, G_C_ST, G_E_ST)
        ERR_C = SQRT(SUM(G_C_ST**2))
        ERR_E = SQRT(SUM(G_E_ST**2))
        ERR_B = ABS(DOT_PRODUCT(N_F, G_A_ST))

        PRINT *, 'Test 9: MFIE 静态主值退化情形'
        PRINT *, '  场点重心，n_f = n_s'
        PRINT *, '  |g_C^st| = ', ERR_C
        PRINT *, '  |g_E^st| = ', ERR_E
        PRINT *, '  |n_f·g_A^st| = ', ERR_B
        IF (ERR_C < TOL .AND. ERR_E < TOL .AND. ERR_B < TOL) THEN
            PRINT *, '  PASSED'
        ELSE
            PRINT *, '  FAILED'
            ALL_PASS = .FALSE.
        END IF
    END SUBROUTINE TEST_MFIE_GRAD_STATIC_DEGENERATE

    SUBROUTINE TEST_MFIE_GRAD_CONSISTENCY(V1, V2, V3, ALL_PASS)
        REAL, INTENT(IN) :: V1(3), V2(3), V3(3)
        LOGICAL, INTENT(INOUT) :: ALL_PASS
        REAL :: N_F(3), R_PT(3), AREA, CROSS(3)
        REAL :: ERR_A, ERR_C, ERR_E
        INTEGER :: I
        TYPE(GAUSS_TRI_DATA) :: GDATA
        REAL, PARAMETER :: K_TEST = 2.0 * PI
        REAL, PARAMETER :: TOL = 0.05

        CROSS = CROSS_PRODUCT(V2 - V1, V3 - V1)
        AREA = 0.5 * SQRT(SUM(CROSS**2))
        R_PT = [0.5, SQRT(3.0)/6.0, 0.5]
        CALL INIT_GAUSS_TRI(GAUSS_12PT, GDATA)

        PRINT *, 'Test 10: MFIE 静态主值+余项 ≡ 完整 ∇G'

        ! 配置 1：n_f = n_s（g_C 退化）
        N_F = CROSS / SQRT(SUM(CROSS**2))
        CALL COMPARE_MFIE_GRAD_MOMENTS(R_PT, N_F, V1, V2, V3, K_TEST, AREA, GDATA, &
                                       ERR_A, ERR_C, ERR_E)
        PRINT *, '  10a n_f=n_s: g_A err=', ERR_A, ' g_C err=', ERR_C, ' g_E err=', ERR_E
        IF (ERR_A < TOL .AND. ERR_C < TOL .AND. ERR_E < TOL) THEN
            PRINT *, '    PASSED'
        ELSE
            PRINT *, '    FAILED'
            ALL_PASS = .FALSE.
        END IF

        ! 配置 2：n_f 倾斜（非退化，全矩）
        N_F = [0.5, 0.0, SQRT(3.0)/2.0]
        N_F = N_F / SQRT(SUM(N_F**2))
        CALL COMPARE_MFIE_GRAD_MOMENTS(R_PT, N_F, V1, V2, V3, K_TEST, AREA, GDATA, &
                                       ERR_A, ERR_C, ERR_E)
        PRINT *, '  10b n_f tilted: g_A err=', ERR_A, ' g_C err=', ERR_C, ' g_E err=', ERR_E
        IF (ERR_A < TOL .AND. ERR_C < TOL .AND. ERR_E < TOL) THEN
            PRINT *, '    PASSED'
        ELSE
            PRINT *, '    FAILED'
            ALL_PASS = .FALSE.
        END IF
    END SUBROUTINE TEST_MFIE_GRAD_CONSISTENCY

    SUBROUTINE COMPARE_MFIE_GRAD_MOMENTS(R_PT, N_F, V1, V2, V3, K, AREA, GDATA, &
                                          ERR_A, ERR_C, ERR_E)
        REAL, INTENT(IN) :: R_PT(3), N_F(3), V1(3), V2(3), V3(3), K, AREA
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        REAL, INTENT(OUT) :: ERR_A, ERR_C, ERR_E
        REAL :: G_A_ST(3), G_C_ST(3), G_E_ST(3)
        COMPLEX :: G_A_ANA(3), G_A_REG(3), GRAD_REG(3), GRAD_G(3)
        COMPLEX :: G_C_ANA(3), G_C_REG(3)
        COMPLEX :: G_E_ANA(3), G_E_REG(3)
        COMPLEX :: G_A_NUM(3), G_C_NUM(3), G_E_NUM(3)
        REAL :: G_A_STATIC_NUM(3), G_C_STATIC_NUM(3), G_E_STATIC_NUM(3)
        REAL :: R_VEC(3), DIST, GRAD_ST(3)
        REAL :: ERR_A_ST, ERR_C_ST, ERR_E_ST
        INTEGER :: I, KIDX
        REAL :: R_J(3)

        CALL ANALYTIC_GRAD_STATIC_PV(R_PT, N_F, V1, V2, V3, G_A_ST, G_C_ST, G_E_ST)

        G_A_REG = (0.0, 0.0); G_C_REG = (0.0, 0.0); G_E_REG = (0.0, 0.0)
        G_A_NUM = (0.0, 0.0); G_C_NUM = (0.0, 0.0); G_E_NUM = (0.0, 0.0)
        G_A_STATIC_NUM = 0.0; G_C_STATIC_NUM = 0.0; G_E_STATIC_NUM = 0.0
        DO I = 1, GDATA%N_POINTS
            DO KIDX = 1, 3
                R_J(KIDX) = V1(KIDX)*GDATA%UVW(1,I) + V2(KIDX)*GDATA%UVW(2,I) + &
                            V3(KIDX)*GDATA%UVW(3,I)
            END DO
            CALL GRAD_GREEN_REG(R_PT, R_J, K, GRAD_REG)
            CALL GARD_GREEN_FUNC(R_PT, R_J, K, GRAD_G)

            G_A_REG = G_A_REG + GDATA%WEIGHTS(I) * GRAD_REG
            G_C_REG = G_C_REG + GDATA%WEIGHTS(I) * R_J * DOT_PRODUCT(N_F, GRAD_REG)
            G_E_REG = G_E_REG + GDATA%WEIGHTS(I) * DOT_PRODUCT(N_F, R_J) * GRAD_REG

            G_A_NUM = G_A_NUM + GDATA%WEIGHTS(I) * GRAD_G
            G_C_NUM = G_C_NUM + GDATA%WEIGHTS(I) * R_J * DOT_PRODUCT(N_F, GRAD_G)
            G_E_NUM = G_E_NUM + GDATA%WEIGHTS(I) * DOT_PRODUCT(N_F, R_J) * GRAD_G

            ! 静态部分直接参考：-(r-r')/R^3（不含 4π）
            R_VEC = R_PT - R_J
            DIST = SQRT(SUM(R_VEC**2))
            GRAD_ST = -R_VEC / DIST**3
            G_A_STATIC_NUM = G_A_STATIC_NUM + GDATA%WEIGHTS(I) * GRAD_ST
            G_C_STATIC_NUM = G_C_STATIC_NUM + GDATA%WEIGHTS(I) * R_J * DOT_PRODUCT(N_F, GRAD_ST)
            G_E_STATIC_NUM = G_E_STATIC_NUM + GDATA%WEIGHTS(I) * DOT_PRODUCT(N_F, R_J) * GRAD_ST
        END DO

        G_A_ANA = CMPLX(G_A_ST)/(4.0*PI) + G_A_REG * AREA
        G_C_ANA = CMPLX(G_C_ST)/(4.0*PI) + G_C_REG * AREA
        G_E_ANA = CMPLX(G_E_ST)/(4.0*PI) + G_E_REG * AREA

        G_A_NUM = G_A_NUM * AREA
        G_C_NUM = G_C_NUM * AREA
        G_E_NUM = G_E_NUM * AREA

        ERR_A = MAXVAL(ABS(G_A_ANA - G_A_NUM)) / MAX(MAXVAL(ABS(G_A_NUM)), 1.0E-12)
        ERR_C = MAXVAL(ABS(G_C_ANA - G_C_NUM)) / MAX(MAXVAL(ABS(G_C_NUM)), 1.0E-12)
        ERR_E = MAXVAL(ABS(G_E_ANA - G_E_NUM)) / MAX(MAXVAL(ABS(G_E_NUM)), 1.0E-12)

        ERR_A_ST = MAXVAL(ABS(G_A_ST - G_A_STATIC_NUM*AREA)) / MAX(MAXVAL(ABS(G_A_STATIC_NUM*AREA)), 1.0E-12)
        ERR_C_ST = MAXVAL(ABS(G_C_ST - G_C_STATIC_NUM*AREA)) / MAX(MAXVAL(ABS(G_C_STATIC_NUM*AREA)), 1.0E-12)
        ERR_E_ST = MAXVAL(ABS(G_E_ST - G_E_STATIC_NUM*AREA)) / MAX(MAXVAL(ABS(G_E_STATIC_NUM*AREA)), 1.0E-12)

        PRINT *, '    静态误差: g_A=', ERR_A_ST, ' g_C=', ERR_C_ST, ' g_E=', ERR_E_ST
        PRINT *, '    g_A_st =', G_A_ST
        PRINT *, '    g_A_num=', G_A_STATIC_NUM*AREA
        PRINT *, '    g_C_st =', G_C_ST
        PRINT *, '    g_C_num=', G_C_STATIC_NUM*AREA
        PRINT *, '    g_E_st =', G_E_ST
        PRINT *, '    g_E_num=', G_E_STATIC_NUM*AREA
    END SUBROUTINE COMPARE_MFIE_GRAD_MOMENTS

END PROGRAM TEST_SINGULAR
