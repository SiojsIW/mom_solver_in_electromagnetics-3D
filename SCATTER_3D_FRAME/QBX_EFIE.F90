MODULE QBX_EFIE
    USE EM_TYPES
    USE NUMERICAL_INTEGRATION
    USE RWG_BASIS_BUILD
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL64
    IMPLICIT NONE

    ! ========================================================================
    ! EFIE-QBX 单侧实现（工程约定 e^{j omega t}）
    ! 参考 docs/QBX展开求积法手册-EFIE-MFIE推导与代码实现-修订版v1.1.md
    ! ========================================================================

    ! ---- 可调参数 ----
    INTEGER, PARAMETER :: P_TRUNC = 12                ! QBX 截断阶
    REAL(REAL64), PARAMETER :: H_FACTOR = 0.3_REAL64  ! h = H_FACTOR * h_obs
    INTEGER, PARAMETER :: N_GAUSS_UP = 25             ! 源单元加密高斯点数
    REAL(REAL64), PARAMETER :: ADM_TOL = 1.0E-3_REAL64! 许可性相对裕量

    ! 虚数单位（双精度）
    COMPLEX(REAL64), PARAMETER :: JJ_DP = (0.0_REAL64, 1.0_REAL64)

    ! ---- 调用统计 ----
    INTEGER :: QBX_CALL_COUNT = 0
    INTEGER :: QBX_ADM_FAIL_COUNT = 0

CONTAINS

    ! ===================================================================
    ! 局部：统计两三角形共享顶点数（避免模块循环依赖）
    ! ===================================================================
    FUNCTION COUNT_COMMON_VERTICES_QBX(MESH, TRI_A, TRI_B) RESULT(N)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_A, TRI_B
        INTEGER :: N
        INTEGER :: I, J, VA(3), VB(3)

        VA = MESH%TRIANGLES(TRI_A)%VERTEX_3D
        VB = MESH%TRIANGLES(TRI_B)%VERTEX_3D
        N = 0
        DO I = 1, 3
            DO J = 1, 3
                IF (VA(I) == VB(J)) N = N + 1
            END DO
        END DO
    END FUNCTION COUNT_COMMON_VERTICES_QBX

    ! ===================================================================
    ! 工具：三维叉乘
    ! ===================================================================
    PURE FUNCTION CROSS_DP(A, B) RESULT(C)
        REAL(REAL64), INTENT(IN) :: A(3), B(3)
        REAL(REAL64) :: C(3)
        C(1) = A(2)*B(3) - A(3)*B(2)
        C(2) = A(3)*B(1) - A(1)*B(3)
        C(3) = A(1)*B(2) - A(2)*B(1)
    END FUNCTION CROSS_DP

    ! ===================================================================
    ! 球 Bessel j_l / Neumann y_l，l = 0..P
    ! ===================================================================
    SUBROUTINE SPH_BESSEL_ALL(X, P, JL, YL)
        REAL(REAL64), INTENT(IN) :: X
        INTEGER, INTENT(IN) :: P
        REAL(REAL64), INTENT(OUT) :: JL(0:P), YL(0:P)

        INTEGER :: L, LMAX
        REAL(REAL64), ALLOCATABLE :: JTMP(:)
        REAL(REAL64) :: J0_TRUE, SCALE

        IF (P < 0) RETURN

        IF (X <= 0.0_REAL64) THEN
            JL = 0.0_REAL64
            YL = 0.0_REAL64
            IF (P >= 0) JL(0) = 1.0_REAL64
            RETURN
        END IF

        ! ---- y_l：升递推无条件稳定 ----
        YL(0) = -COS(X) / X
        IF (P >= 1) YL(1) = -COS(X)/(X*X) - SIN(X)/X
        DO L = 1, P - 1
            YL(L+1) = REAL(2*L + 1, REAL64) / X * YL(L) - YL(L-1)
        END DO

        ! ---- j_l：x 较小或 l > x 时用 Miller 降递推 ----
        IF (X < REAL(P + 1, REAL64) .AND. P > 0) THEN
            LMAX = P + 20
            ALLOCATE(JTMP(0:LMAX))
            JTMP(LMAX) = 0.0_REAL64
            JTMP(LMAX-1) = 1.0_REAL64
            DO L = LMAX - 1, 1, -1
                JTMP(L-1) = REAL(2*L + 1, REAL64) / X * JTMP(L) - JTMP(L+1)
            END DO
            J0_TRUE = SIN(X) / X
            SCALE = J0_TRUE / JTMP(0)
            JL = JTMP(0:P) * SCALE
            DEALLOCATE(JTMP)
        ELSE
            JL(0) = SIN(X) / X
            IF (P >= 1) JL(1) = SIN(X)/(X*X) - COS(X)/X
            DO L = 1, P - 1
                JL(L+1) = REAL(2*L + 1, REAL64) / X * JL(L) - JL(L-1)
            END DO
        END IF
    END SUBROUTINE SPH_BESSEL_ALL

    ! ===================================================================
    ! 球谐函数 Y_l^m，仅 m = 0, +1, -1，归一化约定见文档附录 C
    ! ===================================================================
    SUBROUTINE SPH_HARM_3M(CT, PHI, P, Y0, YP1, YM1)
        REAL(REAL64), INTENT(IN) :: CT, PHI
        INTEGER, INTENT(IN) :: P
        COMPLEX(REAL64), INTENT(OUT) :: Y0(0:P), YP1(0:P), YM1(0:P)

        INTEGER :: L
        REAL(REAL64) :: P_L(0:P), P_L1(0:P), ST, F0, F1
        COMPLEX(REAL64) :: EP, EM

        EP = CMPLX(COS(PHI), SIN(PHI), KIND=REAL64)
        EM = CONJG(EP)

        ST = SQRT(MAX(0.0_REAL64, 1.0_REAL64 - CT*CT))

        ! ---- 连带 Legendre P_l^0 ----
        P_L = 0.0_REAL64
        P_L(0) = 1.0_REAL64
        IF (P >= 1) P_L(1) = CT
        DO L = 1, P - 1
            P_L(L+1) = (REAL(2*L + 1, REAL64) * CT * P_L(L) - REAL(L, REAL64) * P_L(L-1)) / REAL(L + 1, REAL64)
        END DO

        ! ---- 连带 Legendre P_l^1 ----
        P_L1 = 0.0_REAL64
        IF (P >= 1) THEN
            P_L1(1) = -ST
            IF (P >= 2) P_L1(2) = -3.0_REAL64 * CT * ST
            DO L = 2, P - 1
                P_L1(L+1) = (REAL(2*L + 1, REAL64) * CT * P_L1(L) - REAL(L + 1, REAL64) * P_L1(L-1)) / REAL(L, REAL64)
            END DO
        END IF

        ! ---- 球谐 ----
        DO L = 0, P
            F0 = SQRT(REAL(2*L + 1, REAL64) / (4.0_REAL64 * PI))
            Y0(L) = F0 * CMPLX(P_L(L), 0.0_REAL64, KIND=REAL64)
            IF (L >= 1) THEN
                F1 = SQRT(REAL(2*L + 1, REAL64) / (4.0_REAL64 * PI * REAL(L*(L+1), REAL64)))
                YP1(L) = F1 * CMPLX(P_L1(L), 0.0_REAL64, KIND=REAL64) * EP
                YM1(L) = -F1 * CMPLX(P_L1(L), 0.0_REAL64, KIND=REAL64) * EM
            ELSE
                YP1(L) = (0.0_REAL64, 0.0_REAL64)
                YM1(L) = (0.0_REAL64, 0.0_REAL64)
            END IF
        END DO
    END SUBROUTINE SPH_HARM_3M

    ! ===================================================================
    ! 局部坐标架：zp = (r0 - rc)/h，e2 = zp x e1（右手系）
    ! ===================================================================
    SUBROUTINE BUILD_QBX_FRAME(TO_FIELD, TRI_EDGE, E1, E2, ZP)
        REAL(REAL64), INTENT(IN) :: TO_FIELD(3)
        REAL, INTENT(IN) :: TRI_EDGE(3)
        REAL(REAL64), INTENT(OUT) :: E1(3), E2(3), ZP(3)

        REAL(REAL64) :: AUX(3), NORM, DOTP

        ZP = TO_FIELD
        NORM = SQRT(SUM(ZP**2))
        ZP = ZP / NORM

        AUX = REAL(TRI_EDGE, REAL64)
        DOTP = DOT_PRODUCT(AUX, ZP)
        AUX = AUX - DOTP * ZP
        NORM = SQRT(SUM(AUX**2))

        IF (NORM < 1.0E-14_REAL64) THEN
            ! 边与 zp 平行，改用全局 x
            AUX = [1.0_REAL64, 0.0_REAL64, 0.0_REAL64]
            DOTP = DOT_PRODUCT(AUX, ZP)
            AUX = AUX - DOTP * ZP
            NORM = SQRT(SUM(AUX**2))
            IF (NORM < 1.0E-14_REAL64) THEN
                AUX = [0.0_REAL64, 1.0_REAL64, 0.0_REAL64]
                DOTP = DOT_PRODUCT(AUX, ZP)
                AUX = AUX - DOTP * ZP
                NORM = SQRT(SUM(AUX**2))
            END IF
        END IF
        E1 = AUX / NORM
        E2 = CROSS_DP(ZP, E1)
    END SUBROUTINE BUILD_QBX_FRAME

    ! ===================================================================
    ! 点到线段距离
    ! ===================================================================
    REAL(REAL64) FUNCTION POINT_SEGMENT_DIST_DP(P, A, B) RESULT(D)
        REAL(REAL64), INTENT(IN) :: P(3), A(3), B(3)
        REAL(REAL64) :: AB(3), T
        AB = B - A
        T = DOT_PRODUCT(P - A, AB) / DOT_PRODUCT(AB, AB)
        T = MAX(0.0_REAL64, MIN(1.0_REAL64, T))
        D = SQRT(SUM((A + T*AB - P)**2))
    END FUNCTION POINT_SEGMENT_DIST_DP

    ! ===================================================================
    ! 点到三角形最短距离
    ! ===================================================================
    REAL(REAL64) FUNCTION POINT_TRI_MIN_DIST_DP(RC, V1, V2, V3) RESULT(DMIN)
        REAL(REAL64), INTENT(IN) :: RC(3), V1(3), V2(3), V3(3)

        REAL(REAL64) :: E0(3), E1(3), N(3), D, S, T, DET, P(3), CLOSEST(3)
        REAL(REAL64) :: D1, D2, D3

        E0 = V2 - V1
        E1 = V3 - V1
        N = CROSS_DP(E0, E1)
        N = N / SQRT(SUM(N**2))

        ! 面上投影
        D = DOT_PRODUCT(RC - V1, N)
        P = RC - D * N

        ! 解 P = V1 + s*E0 + t*E1
        DET = DOT_PRODUCT(E0, E0)*DOT_PRODUCT(E1, E1) - DOT_PRODUCT(E0, E1)**2
        S = (DOT_PRODUCT(P - V1, E0)*DOT_PRODUCT(E1, E1) - DOT_PRODUCT(P - V1, E1)*DOT_PRODUCT(E0, E1)) / DET
        T = (DOT_PRODUCT(P - V1, E1)*DOT_PRODUCT(E0, E0) - DOT_PRODUCT(P - V1, E0)*DOT_PRODUCT(E0, E1)) / DET

        IF (S >= 0.0_REAL64 .AND. T >= 0.0_REAL64 .AND. (S + T) <= 1.0_REAL64) THEN
            CLOSEST = P
        ELSE
            D1 = POINT_SEGMENT_DIST_DP(RC, V1, V2)
            D2 = POINT_SEGMENT_DIST_DP(RC, V2, V3)
            D3 = POINT_SEGMENT_DIST_DP(RC, V3, V1)
            DMIN = MIN(D1, MIN(D2, D3))
            RETURN
        END IF

        DMIN = SQRT(SUM((RC - CLOSEST)**2))
    END FUNCTION POINT_TRI_MIN_DIST_DP

    ! ===================================================================
    ! QBX 系数：alpha_{l,m}, beta_{l,m}（m = 0, +/-1）
    ! ===================================================================
    SUBROUTINE QBX_EFIE_COEFFICIENTS(RC, E1, E2, ZP, K, A_SRC, V_SRC, &
                                     SIGMA, R_OPP_S, C_S, GAUSS_UP, &
                                     ALPHA, BETA)
        REAL(REAL64), INTENT(IN) :: RC(3), E1(3), E2(3), ZP(3)
        REAL, INTENT(IN) :: K
        REAL, INTENT(IN) :: A_SRC
        REAL, INTENT(IN) :: V_SRC(3, 3)
        REAL, INTENT(IN) :: SIGMA
        REAL, INTENT(IN) :: R_OPP_S(3)
        REAL, INTENT(IN) :: C_S
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GAUSS_UP
        COMPLEX(REAL64), INTENT(OUT) :: ALPHA(0:P_TRUNC, -1:1)
        COMPLEX(REAL64), INTENT(OUT) :: BETA(0:P_TRUNC, -1:1, 3)

        INTEGER :: IQ, L, KK
        REAL(REAL64) :: RQ(3), DR(3), RHO, CT, PHI, KR
        REAL(REAL64) :: JL(0:P_TRUNC), YL(0:P_TRUNC)
        COMPLEX(REAL64) :: H2(0:P_TRUNC)
        COMPLEX(REAL64) :: Y0(0:P_TRUNC), YP1(0:P_TRUNC), YM1(0:P_TRUNC)
        REAL(REAL64) :: WQ
        COMPLEX(REAL64) :: PHASE
        REAL :: F_N(3)

        ALPHA = (0.0_REAL64, 0.0_REAL64)
        BETA = (0.0_REAL64, 0.0_REAL64)

        DO IQ = 1, GAUSS_UP%N_POINTS
            ! 源高斯点全局坐标
            DO KK = 1, 3
                RQ(KK) = REAL(V_SRC(KK,1), REAL64) * REAL(GAUSS_UP%UVW(1, IQ), REAL64) + &
                         REAL(V_SRC(KK,2), REAL64) * REAL(GAUSS_UP%UVW(2, IQ), REAL64) + &
                         REAL(V_SRC(KK,3), REAL64) * REAL(GAUSS_UP%UVW(3, IQ), REAL64)
            END DO

            DR = RQ - RC
            RHO = SQRT(SUM(DR**2))
            CT = DOT_PRODUCT(DR, ZP) / RHO
            PHI = ATAN2(DOT_PRODUCT(DR, E2), DOT_PRODUCT(DR, E1))
            KR = REAL(K, REAL64) * RHO

            CALL SPH_BESSEL_ALL(KR, P_TRUNC, JL, YL)
            H2 = JL - JJ_DP * YL

            CALL SPH_HARM_3M(CT, PHI, P_TRUNC, Y0, YP1, YM1)

            WQ = REAL(GAUSS_UP%WEIGHTS(IQ), REAL64)
            PHASE = -JJ_DP * REAL(K, REAL64) * WQ * REAL(A_SRC, REAL64)

            ! 标量密度 sigma
            ALPHA(:, 0) = ALPHA(:, 0) + PHASE * REAL(SIGMA, REAL64) * H2(:) * CONJG(Y0(:))
            ALPHA(:, 1) = ALPHA(:, 1) + PHASE * REAL(SIGMA, REAL64) * H2(:) * CONJG(YP1(:))
            ALPHA(:, -1) = ALPHA(:, -1) + PHASE * REAL(SIGMA, REAL64) * H2(:) * CONJG(YM1(:))

            ! 矢量源 f_n = C_S * (r' - r_opp)
            DO KK = 1, 3
                F_N(KK) = C_S * (REAL(RQ(KK)) - R_OPP_S(KK))
            END DO
            DO L = 0, P_TRUNC
                BETA(L, 0, :) = BETA(L, 0, :) + PHASE * H2(L) * CONJG(Y0(L)) * REAL(F_N, REAL64)
                BETA(L, 1, :) = BETA(L, 1, :) + PHASE * H2(L) * CONJG(YP1(L)) * REAL(F_N, REAL64)
                BETA(L, -1, :) = BETA(L, -1, :) + PHASE * H2(L) * CONJG(YM1(L)) * REAL(F_N, REAL64)
            END DO
        END DO
    END SUBROUTINE QBX_EFIE_COEFFICIENTS

    ! ===================================================================
    ! QBX 极点求值：S(r0), V(r0)（仅 m=0 参与）
    ! ===================================================================
    SUBROUTINE QBX_EFIE_EVAL(ALPHA, BETA, KH, S_VAL, V_VAL)
        COMPLEX(REAL64), INTENT(IN) :: ALPHA(0:P_TRUNC, -1:1)
        COMPLEX(REAL64), INTENT(IN) :: BETA(0:P_TRUNC, -1:1, 3)
        REAL(REAL64), INTENT(IN) :: KH
        COMPLEX(REAL64), INTENT(OUT) :: S_VAL
        COMPLEX(REAL64), INTENT(OUT) :: V_VAL(3)

        INTEGER :: L
        REAL(REAL64) :: JL(0:P_TRUNC), YL(0:P_TRUNC), FACTOR

        CALL SPH_BESSEL_ALL(KH, P_TRUNC, JL, YL)

        S_VAL = (0.0_REAL64, 0.0_REAL64)
        V_VAL = (0.0_REAL64, 0.0_REAL64)

        DO L = 0, P_TRUNC
            FACTOR = JL(L) * SQRT(REAL(2*L + 1, REAL64) / (4.0_REAL64 * PI))
            S_VAL = S_VAL + ALPHA(L, 0) * FACTOR
            V_VAL = V_VAL + BETA(L, 0, :) * FACTOR
        END DO
    END SUBROUTINE QBX_EFIE_EVAL

    ! ===================================================================
    ! 提取三角形最大边长（单精度）
    ! ===================================================================
    REAL FUNCTION TRI_MAX_EDGE_LEN_QBX(V) RESULT(HMAX)
        REAL, INTENT(IN) :: V(3, 3)
        REAL :: E1, E2, E3
        E1 = SQRT(SUM((V(:,1) - V(:,2))**2))
        E2 = SQRT(SUM((V(:,2) - V(:,3))**2))
        E3 = SQRT(SUM((V(:,3) - V(:,1))**2))
        HMAX = MAX(E1, MAX(E2, E3))
    END FUNCTION TRI_MAX_EDGE_LEN_QBX

    ! ===================================================================
    ! EFIE-QBX 单个三角形对阻抗（近场分支通用）
    ! ===================================================================
    SUBROUTINE CALC_EFIE_QBX_PAIR(MESH, TRI_OBS, TRI_SRC, RWG_M, RWG_N, K, ETA0, Z_PAIR, SUCCESS)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_OBS, TRI_SRC
        TYPE(RWG_BASIS), INTENT(IN) :: RWG_M, RWG_N
        REAL, INTENT(IN) :: K, ETA0
        COMPLEX, INTENT(OUT) :: Z_PAIR
        LOGICAL, INTENT(OUT) :: SUCCESS

        REAL :: C_O, C_S, R_OPP_O(3), R_OPP_S(3)
        REAL :: SIGMA_O, SIGMA_S
        REAL :: V_OBS(3, 3), V_SRC(3, 3)
        REAL :: N_OBS(3), A_OBS, A_SRC
        REAL :: H_OBS
        REAL(REAL64) :: RC(3), E1(3), E2(3), ZP(3), R0(3), KH, H_DP
        COMPLEX(REAL64) :: ALPHA(0:P_TRUNC, -1:1), BETA(0:P_TRUNC, -1:1, 3)
        REAL(REAL64) :: F_M_DP(3)
        COMPLEX(REAL64) :: S_VAL, V_VAL(3), F_DOT_V
        COMPLEX(REAL64) :: JW_MU_DP, J_OVER_WE_DP
        COMPLEX :: Z_TERM
        REAL :: TRI_EDGE(3)
        INTEGER :: IQ, KK
        TYPE(GAUSS_TRI_DATA) :: GDATA_UP, GDATA_OUTER
        REAL :: GLOBAL_PTS_O(3, 12)
        REAL :: DMIN
        LOGICAL :: ADMISSIBLE

        SUCCESS = .TRUE.
        QBX_CALL_COUNT = QBX_CALL_COUNT + 1

        ! ---- 提取 RWG 在观察/源三角形上的系数、对顶点、散度 ----
        IF (TRI_OBS == RWG_M%POS_TRI_ID) THEN
            C_O = RWG_M%POS_COEF
            R_OPP_O(1) = MESH%NODES(RWG_M%POS_OPP_VERTEX)%X
            R_OPP_O(2) = MESH%NODES(RWG_M%POS_OPP_VERTEX)%Y
            R_OPP_O(3) = MESH%NODES(RWG_M%POS_OPP_VERTEX)%Z
        ELSE
            C_O = -RWG_M%NEG_COEF
            R_OPP_O(1) = MESH%NODES(RWG_M%NEG_OPP_VERTEX)%X
            R_OPP_O(2) = MESH%NODES(RWG_M%NEG_OPP_VERTEX)%Y
            R_OPP_O(3) = MESH%NODES(RWG_M%NEG_OPP_VERTEX)%Z
        END IF
        SIGMA_O = 2.0 * C_O

        IF (TRI_SRC == RWG_N%POS_TRI_ID) THEN
            C_S = RWG_N%POS_COEF
            R_OPP_S(1) = MESH%NODES(RWG_N%POS_OPP_VERTEX)%X
            R_OPP_S(2) = MESH%NODES(RWG_N%POS_OPP_VERTEX)%Y
            R_OPP_S(3) = MESH%NODES(RWG_N%POS_OPP_VERTEX)%Z
        ELSE
            C_S = -RWG_N%NEG_COEF
            R_OPP_S(1) = MESH%NODES(RWG_N%NEG_OPP_VERTEX)%X
            R_OPP_S(2) = MESH%NODES(RWG_N%NEG_OPP_VERTEX)%Y
            R_OPP_S(3) = MESH%NODES(RWG_N%NEG_OPP_VERTEX)%Z
        END IF
        SIGMA_S = 2.0 * C_S

        ! ---- 几何量 ----
        N_OBS = MESH%TRIANGLES(TRI_OBS)%NORMAL
        A_OBS = MESH%TRIANGLES(TRI_OBS)%AREA
        A_SRC = MESH%TRIANGLES(TRI_SRC)%AREA

        DO KK = 1, 3
            V_OBS(1, KK) = MESH%NODES(MESH%TRIANGLES(TRI_OBS)%VERTEX_3D(KK))%X
            V_OBS(2, KK) = MESH%NODES(MESH%TRIANGLES(TRI_OBS)%VERTEX_3D(KK))%Y
            V_OBS(3, KK) = MESH%NODES(MESH%TRIANGLES(TRI_OBS)%VERTEX_3D(KK))%Z
            V_SRC(1, KK) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(KK))%X
            V_SRC(2, KK) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(KK))%Y
            V_SRC(3, KK) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(KK))%Z
        END DO

        H_OBS = TRI_MAX_EDGE_LEN_QBX(V_OBS)
        TRI_EDGE = V_OBS(:, 2) - V_OBS(:, 1)

        ! ---- 初始化高斯数据 ----
        CALL INIT_GAUSS_TRI(GAUSS_25PT, GDATA_UP)
        CALL INIT_GAUSS_TRI(GAUSS_12PT, GDATA_OUTER)

        ! ---- 外层测试高斯点循环 ----
        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_OBS, GDATA_OUTER, GLOBAL_PTS_O)

        JW_MU_DP = JJ_DP * REAL(ETA0 * K, REAL64)
        J_OVER_WE_DP = JJ_DP * REAL(ETA0 / K, REAL64)

        Z_TERM = (0.0, 0.0)

        DO IQ = 1, GDATA_OUTER%N_POINTS
            R0 = REAL(GLOBAL_PTS_O(:, IQ), REAL64)

            ! 中心：r_c = r0 + h * n_hat
            H_DP = H_FACTOR * REAL(H_OBS, REAL64)
            RC = R0 + H_DP * REAL(N_OBS, REAL64)

            ! 局部架：zp = (r0 - rc)/h = -n_hat
            CALL BUILD_QBX_FRAME(-REAL(N_OBS, REAL64), TRI_EDGE, E1, E2, ZP)

            ! 许可性检查：自项、共边、共顶点对均跳过；仅非接触近场检查 min_dist >= h*(1+eps)
            ADMISSIBLE = .TRUE.
            IF (TRI_OBS /= TRI_SRC) THEN
                IF (COUNT_COMMON_VERTICES_QBX(MESH, TRI_OBS, TRI_SRC) == 0) THEN
                    DMIN = REAL(POINT_TRI_MIN_DIST_DP(RC, &
                        REAL(V_SRC(:,1), REAL64), REAL(V_SRC(:,2), REAL64), REAL(V_SRC(:,3), REAL64)))
                    IF (REAL(DMIN, REAL64) < H_DP * (1.0_REAL64 + ADM_TOL)) THEN
                        ADMISSIBLE = .FALSE.
                    END IF
                END IF
            END IF

            IF (.NOT. ADMISSIBLE) THEN
                ! 由上层回退到提取法
                SUCCESS = .FALSE.
                QBX_ADM_FAIL_COUNT = QBX_ADM_FAIL_COUNT + 1
                Z_PAIR = (0.0, 0.0)
                RETURN
            END IF

            ! 计算 QBX 系数与极点求值
            KH = REAL(K, REAL64) * H_DP
            CALL QBX_EFIE_COEFFICIENTS(RC, E1, E2, ZP, K, A_SRC, V_SRC, &
                                       SIGMA_S, R_OPP_S, C_S, GDATA_UP, ALPHA, BETA)
            CALL QBX_EFIE_EVAL(ALPHA, BETA, KH, S_VAL, V_VAL)

            ! RWG 测试函数 f_m 在当前场点
            DO KK = 1, 3
                F_M_DP(KK) = REAL(C_O, REAL64) * (R0(KK) - REAL(R_OPP_O(KK), REAL64))
            END DO

            F_DOT_V = DOT_PRODUCT(F_M_DP, V_VAL)

            Z_TERM = Z_TERM + CMPLX(GDATA_OUTER%WEIGHTS(IQ) * A_OBS) * &
                     CMPLX(JW_MU_DP * F_DOT_V - J_OVER_WE_DP * REAL(SIGMA_O, REAL64) * S_VAL)
        END DO

        Z_PAIR = Z_TERM
    END SUBROUTINE CALC_EFIE_QBX_PAIR

END MODULE QBX_EFIE