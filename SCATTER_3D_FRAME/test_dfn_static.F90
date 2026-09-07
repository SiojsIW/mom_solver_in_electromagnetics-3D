! ============================================================================
! test_dfn_static.F90
!
! DIRECTFN 核归一化静态验证（k0 = 0）：
! 对正三角形用已验证的解析 1/R 势（SINGULAR_INTEGRAL 模块）构造
! ∬(r-v_i)·(r'-v_j)/R dS dS' 的参考值，检验 DIRECTFN 返回的
! iss_rwg(i,j) 是否严格等于 (l_i*l_j/(4A^2)) * P_ref。
! 同时检验 iss_const 是否等于 J0/(4π)。
! ============================================================================

PROGRAM TEST_DFN_STATIC
    USE ISO_C_BINDING
    USE EM_TYPES
    USE NUMERICAL_INTEGRATION
    USE SINGULAR_INTEGRAL
    IMPLICIT NONE

    INTERFACE
        FUNCTION C_DFN_TRI_ISS(ADJ, PTS, K0, NGAUSS, ISS_RWG, ISS_CONST) &
                BIND(C, NAME='dfn_tri_iss') RESULT(IRC)
            IMPORT :: C_INT, C_DOUBLE, C_DOUBLE_COMPLEX
            INTEGER(C_INT), VALUE :: ADJ
            REAL(C_DOUBLE), INTENT(IN) :: PTS(3, 7)
            REAL(C_DOUBLE), VALUE :: K0
            INTEGER(C_INT), VALUE :: NGAUSS
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: ISS_RWG(9)
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: ISS_CONST(1)
            INTEGER(C_INT) :: IRC
        END FUNCTION C_DFN_TRI_ISS
    END INTERFACE

    REAL(C_DOUBLE) :: PTS(3, 7)
    COMPLEX(C_DOUBLE_COMPLEX) :: ISS_RWG(9), ISS_CONST(1)
    INTEGER(C_INT) :: RC

    REAL :: V(3, 3)          ! 三角形顶点
    REAL :: AREA, EL(3)      ! 面积、对边边长 l_i
    REAL :: J0               ! ∬ 1/R dS dS'
    REAL :: S_LIN(3)         ! S_i = ∬∫ λ_i(r)/R dS dS'（二次矩的中间量）
    REAL :: AI(3)            ! AI_a = ∬ r_a/R dS dS'
    REAL :: L_MAT(3, 3)      ! L_ij = ∬ λ_i(r) λ_j(r')/R dS dS'
    REAL :: T_AB             ! ∬ r·r'/R = Σ_a T_aa
    REAL :: P_REF, RWG_REF, COEF_IJ
    REAL :: ERR, MAX_ERR
    INTEGER :: I, J, A, II, JJ
    TYPE(GAUSS_TRI_DATA) :: GDATA
    REAL :: R_PT(3), C_ANA(3)
    INTEGER :: K

    ! 正三角形：边长 1，位于 z=0 平面
    V(:, 1) = [0.0, 0.0, 0.0]
    V(:, 2) = [1.0, 0.0, 0.0]
    V(:, 3) = [0.5, SQRT(3.0)/2.0, 0.0]
    AREA = SQRT(3.0)/4.0
    EL = 1.0   ! 正三角形三边等长

    CALL INIT_GAUSS_TRI(GAUSS_12PT, GDATA)

    ! ---- 参考矩：J0 ----
    J0 = 0.0
    DO I = 1, GDATA%N_POINTS
        R_PT = 0.0
        DO A = 1, 3
            R_PT = R_PT + V(:, A) * GDATA%UVW(A, I)
        END DO
        CALL ANALYTIC_SCALAR_POT_1_OVER_R(R_PT, V(:,1), V(:,2), V(:,3), C_ANA(1))
        J0 = J0 + GDATA%WEIGHTS(I) * C_ANA(1) * AREA
    END DO

    ! ---- 参考矩：S_i = ∫ c_i(r) dS（c_i 为解析线性势） ----
    S_LIN = 0.0
    L_MAT = 0.0
    DO I = 1, GDATA%N_POINTS
        R_PT = 0.0
        DO A = 1, 3
            R_PT = R_PT + V(:, A) * GDATA%UVW(A, I)
        END DO
        CALL ANALYTIC_LINEAR_POT_1_OVER_R(R_PT, V(:,1), V(:,2), V(:,3), C_ANA)
        S_LIN = S_LIN + GDATA%WEIGHTS(I) * C_ANA * AREA
        DO J = 1, 3
            L_MAT(:, J) = L_MAT(:, J) + GDATA%WEIGHTS(I) * GDATA%UVW(:, I) * C_ANA(J) * AREA
        END DO
    END DO

    ! AI_a = ∬ r_a/R = Σ_i v_i,a * S_i
    AI = 0.0
    DO A = 1, 3
        AI(A) = DOT_PRODUCT(V(A, :), S_LIN)
    END DO

    ! ---- 调用 DIRECTFN（k0 = 0） ----
    PTS = 0.0_C_DOUBLE
    DO I = 1, 3
        PTS(:, I) = REAL(V(:, I), C_DOUBLE)
    END DO
    RC = C_DFN_TRI_ISS(3_C_INT, PTS, 0.0_C_DOUBLE, 12_C_INT, ISS_RWG, ISS_CONST)
    IF (RC /= 0_C_INT) THEN
        PRINT *, "DIRECTFN 调用失败 RC=", RC
        STOP 1
    END IF

    PRINT *, "============================================================"
    PRINT *, " DIRECTFN 核归一化静态验证（k0=0, 正三角形）"
    PRINT *, "============================================================"

    ! ---- 检验 1：iss_const = J0/(4π) ----
    ! 注：ANALYTIC_SCALAR/LINEAR_POT 在场点位于三角形内部（本项目自项路径）时
    ! 存在 ~2.5e-3 的系统误差（已与 DIRECTFN/Python 独立参考对比确认），
    ! 因此此处容差取 1e-2，主要用于捕捉归一化错误（如 4π、12 倍因子）。
    ERR = ABS(REAL(ISS_CONST(1), KIND(AREA)) - J0/(4.0*PI)) / (J0/(4.0*PI))
    PRINT '(A)', " [检验 1] iss_const vs J0/(4π)"
    PRINT '(A, E20.12)', "   DIRECTFN = ", REAL(ISS_CONST(1))
    PRINT '(A, E20.12)', "   参考 J0/(4π) = ", J0/(4.0*PI)
    PRINT '(A, ES12.3)', "   相对误差 = ", ERR
    IF (ERR < 1.0E-2) THEN
        PRINT '(A)', "   PASSED"
    ELSE
        PRINT '(A)', "   FAILED"
    END IF

    ! ---- 检验 2：iss_rwg(i,j) = (l_i l_j/4A²)·P_ref(i,j) ----
    PRINT '(A)', " [检验 2] iss_rwg(i,j) vs (l_i*l_j/4A^2)*P_ref"
    MAX_ERR = 0.0
    DO II = 1, 3
        DO JJ = 1, 3
            ! P_ref = ∬(r-v_i)·(r'-v_j)/R dS dS'
            ! = Σ_a T_aa - v_j·AI - v_i·AI + v_i·v_j·J0
            ! 其中 Σ_a T_aa = Σ_ij v_i·v_j L_ij
            T_AB = 0.0
            DO I = 1, 3
                DO J = 1, 3
                    T_AB = T_AB + DOT_PRODUCT(V(:, I), V(:, J)) * L_MAT(I, J)
                END DO
            END DO
            P_REF = T_AB &
                  - DOT_PRODUCT(V(:, JJ), AI) &
                  - DOT_PRODUCT(V(:, II), AI) &
                  + DOT_PRODUCT(V(:, II), V(:, JJ)) * J0
            COEF_IJ = EL(II) * EL(JJ) / (4.0 * AREA * AREA)
            ! 解析势为不带 4π 的裸 1/R 势，格林函数 G = 1/(4πR)
            RWG_REF = COEF_IJ * P_REF / (4.0 * PI)
            ERR = ABS(REAL(ISS_RWG(II + 3*(JJ-1)), KIND(AREA)) - RWG_REF) / MAX(ABS(RWG_REF), 1.0E-30)
            MAX_ERR = MAX(MAX_ERR, ERR)
            PRINT '(A, I0, A, I0, A, E16.8, A, E16.8, A, ES10.2)', &
                "   (", II, ",", JJ, ") DIRECTFN=", REAL(ISS_RWG(II + 3*(JJ-1))), &
                "  参考=", RWG_REF, "  相对误差=", ERR
        END DO
    END DO
    IF (MAX_ERR < 1.0E-2) THEN
        PRINT '(A)', "   PASSED"
    ELSE
        PRINT '(A)', "   FAILED"
    END IF
    PRINT *, "============================================================"

END PROGRAM TEST_DFN_STATIC
