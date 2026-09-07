MODULE Z_MATRIX
    USE EM_TYPES
    USE NUMERICAL_INTEGRATION
    USE GREEN_FUNCTIONS
    USE RWG_BASIS_BUILD
    USE SINGULAR_INTEGRAL
    USE QBX_EFIE
    USE DIRECTFN_BRIDGE
    IMPLICIT NONE

    ! 阈值参数
    REAL, PARAMETER :: CHI_NEAR = 1.0      ! chi = d_min / h_pair 的近场阈值
    INTEGER, PARAMETER :: N_SUBDIV_EDGE = 2    ! 共边：细分层数
    INTEGER, PARAMETER :: N_SUBDIV_VERTEX = 2  ! 共顶点：细分层数
    INTEGER, PARAMETER :: N_SUBDIV_NEAR = 1    ! 近场非接触：细分层数

    ! ---- EFIE 近场分支切换 ----
    INTEGER, PARAMETER :: EFIE_NEAR_EXTRACTION = 1
    INTEGER, PARAMETER :: EFIE_NEAR_QBX = 2
    INTEGER, PARAMETER :: EFIE_NEAR_DIRECTFN = 3
    INTEGER :: EFIE_NEAR_METHOD = EFIE_NEAR_EXTRACTION

CONTAINS

    ! ================ 分类函数 ================

    ! 统计两个三角形共享的顶点数
    FUNCTION COUNT_COMMON_VERTICES(MESH, TRI_A, TRI_B) RESULT(N)
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
    END FUNCTION COUNT_COMMON_VERTICES

    ! 判断两个三角形是否相邻（恰好共享一条边，即共享两个顶点）
    FUNCTION ARE_TRIANGLES_ADJACENT(MESH, TRI_A, TRI_B) RESULT(ADJACENT)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_A, TRI_B
        LOGICAL :: ADJACENT

        ADJACENT = (COUNT_COMMON_VERTICES(MESH, TRI_A, TRI_B) == 2)
    END FUNCTION

    ! 三角形质心之间的最小距离（用于远/近场判断）
    FUNCTION TRI_CENTROID_DIST(MESH, TRI_A, TRI_B) RESULT(D)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_A, TRI_B
        REAL :: D

        D = SQRT(SUM((MESH%TRIANGLES(TRI_A)%CENTROID - MESH%TRIANGLES(TRI_B)%CENTROID)**2))
    END FUNCTION

    ! 三角形最大边长
    FUNCTION TRI_MAX_EDGE_LEN(MESH, TRI_ID) RESULT(H_MAX)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_ID
        REAL :: H_MAX
        INTEGER :: V1, V2, V3
        REAL :: P1(3), P2(3), P3(3), E1, E2, E3

        V1 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(1)
        V2 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(2)
        V3 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(3)

        P1(1) = MESH%NODES(V1)%X; P1(2) = MESH%NODES(V1)%Y; P1(3) = MESH%NODES(V1)%Z
        P2(1) = MESH%NODES(V2)%X; P2(2) = MESH%NODES(V2)%Y; P2(3) = MESH%NODES(V2)%Z
        P3(1) = MESH%NODES(V3)%X; P3(2) = MESH%NODES(V3)%Y; P3(3) = MESH%NODES(V3)%Z

        E1 = SQRT(SUM((P1 - P2)**2))
        E2 = SQRT(SUM((P2 - P3)**2))
        E3 = SQRT(SUM((P3 - P1)**2))
        H_MAX = MAX(E1, MAX(E2, E3))
    END FUNCTION

    ! ================ 三角形细分 ================

    ! 将三角形 1→4 等面积细分（连接三边中点）
    SUBROUTINE SUBDIVIDE_TRI_4(V1, V2, V3, SUB_V)
        REAL, INTENT(IN)  :: V1(3), V2(3), V3(3)
        REAL, INTENT(OUT) :: SUB_V(3, 3, 4)  ! [xyz, vertex, sub-tri]
        REAL :: M12(3), M23(3), M31(3)
        INTEGER :: K

        DO K = 1, 3
            M12(K) = (V1(K) + V2(K)) / 2.0
            M23(K) = (V2(K) + V3(K)) / 2.0
            M31(K) = (V3(K) + V1(K)) / 2.0
        END DO

        ! Sub 1: V1, M12, M31
        SUB_V(:, 1, 1) = V1;  SUB_V(:, 2, 1) = M12; SUB_V(:, 3, 1) = M31
        ! Sub 2: V2, M23, M12
        SUB_V(:, 1, 2) = V2;  SUB_V(:, 2, 2) = M23; SUB_V(:, 3, 2) = M12
        ! Sub 3: V3, M31, M23
        SUB_V(:, 1, 3) = V3;  SUB_V(:, 2, 3) = M31; SUB_V(:, 3, 3) = M23
        ! Sub 4: M12, M23, M31
        SUB_V(:, 1, 4) = M12; SUB_V(:, 2, 4) = M23; SUB_V(:, 3, 4) = M31
    END SUBROUTINE

    ! ================ 提取法阻抗（通用于 1-4 类） ================

    ! 对任意三角形对（T_o 观察, T_s 源）用奇异性提取法计算阻抗。
    ! N_SUBDIV_LEVELS: 外层细分层数（0=不细分）
    SUBROUTINE CALC_EFIE_EXTRACTED_PAIR(MESH, TRI_OBS, TRI_SRC, &
                                         RWG_M, RWG_N, K, ETA0, N_SUBDIV_LEVELS, Z_PAIR)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_OBS, TRI_SRC
        TYPE(RWG_BASIS), INTENT(IN) :: RWG_M, RWG_N
        REAL, INTENT(IN) :: K, ETA0
        INTEGER, INTENT(IN) :: N_SUBDIV_LEVELS
        COMPLEX, INTENT(OUT) :: Z_PAIR

        REAL :: C_O, C_S
        REAL :: R_OPP_O(3), R_OPP_S(3)
        COMPLEX :: JW_MU, J_OVER_WE
        REAL :: V_OBS(3, 3), V_SRC(3, 3)
        REAL :: SUM_I0, SUM_R_I0(3), SUM_IVEC(3), SUM_RDOT
        REAL :: J0, J_R_I0(3), J_VEC(3), J_R_DOT
        COMPLEX :: P_SING, Q_SING
        COMPLEX :: I1_SMOOTH, I4_SMOOTH
        COMPLEX :: I2_SMOOTH(3), I3_SMOOTH(3)
        COMPLEX :: P_SMOOTH, Q_SMOOTH
        TYPE(GAUSS_TRI_DATA) :: GDATA_INNER, GDATA_OUTER

        ! 提取 RWG 系数和对顶点
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

        JW_MU = (0.0, 1.0) * ETA0 * K
        J_OVER_WE = (0.0, 1.0) * ETA0 / K

        ! 观察三角形顶点
        V_OBS(1, 1) = MESH%NODES(MESH%TRIANGLES(TRI_OBS)%VERTEX_3D(1))%X
        V_OBS(2, 1) = MESH%NODES(MESH%TRIANGLES(TRI_OBS)%VERTEX_3D(1))%Y
        V_OBS(3, 1) = MESH%NODES(MESH%TRIANGLES(TRI_OBS)%VERTEX_3D(1))%Z
        V_OBS(1, 2) = MESH%NODES(MESH%TRIANGLES(TRI_OBS)%VERTEX_3D(2))%X
        V_OBS(2, 2) = MESH%NODES(MESH%TRIANGLES(TRI_OBS)%VERTEX_3D(2))%Y
        V_OBS(3, 2) = MESH%NODES(MESH%TRIANGLES(TRI_OBS)%VERTEX_3D(2))%Z
        V_OBS(1, 3) = MESH%NODES(MESH%TRIANGLES(TRI_OBS)%VERTEX_3D(3))%X
        V_OBS(2, 3) = MESH%NODES(MESH%TRIANGLES(TRI_OBS)%VERTEX_3D(3))%Y
        V_OBS(3, 3) = MESH%NODES(MESH%TRIANGLES(TRI_OBS)%VERTEX_3D(3))%Z

        V_SRC(1, 1) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(1))%X
        V_SRC(2, 1) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(1))%Y
        V_SRC(3, 1) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(1))%Z
        V_SRC(1, 2) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(2))%X
        V_SRC(2, 2) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(2))%Y
        V_SRC(3, 2) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(2))%Z
        V_SRC(1, 3) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(3))%X
        V_SRC(2, 3) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(3))%Y
        V_SRC(3, 3) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(3))%Z

        ! ---- 外层高斯积分（静态部分，含可选细分） ----
        CALL INIT_GAUSS_TRI(GAUSS_12PT, GDATA_OUTER)
        SUM_I0 = 0.0; SUM_R_I0 = 0.0; SUM_IVEC = 0.0; SUM_RDOT = 0.0

        ! 递归细分并积分
        CALL INTEGRATE_STATIC_RECURSIVE(V_OBS(:,1), V_OBS(:,2), V_OBS(:,3), &
                                         V_SRC, 0, N_SUBDIV_LEVELS, GDATA_OUTER, &
                                         SUM_I0, SUM_R_I0, SUM_IVEC, SUM_RDOT)

        J0     = SUM_I0
        J_R_I0 = SUM_R_I0
        J_VEC  = SUM_IVEC
        J_R_DOT = SUM_RDOT

        ! 静态部分（公式对自项和非自项一致）
        Q_SING = -J0 / PI
        P_SING = (J_R_DOT &
                  - DOT_PRODUCT(R_OPP_S, J_R_I0) &
                  - DOT_PRODUCT(R_OPP_O, J_VEC) &
                  + DOT_PRODUCT(R_OPP_O, R_OPP_S) * J0) / (4.0 * PI)

        ! ---- G_rem 数值部分（光滑，普通 12 点即可） ----
        CALL INIT_GAUSS_TRI(GAUSS_12PT, GDATA_INNER)
        CALL CALC_GREEN_SMOOTH_INTEGALS(MESH, TRI_OBS, TRI_SRC, GDATA_INNER, K, &
                                         I1_SMOOTH, I2_SMOOTH, I3_SMOOTH, I4_SMOOTH)
        ! 注意：CALC_GREEN_SMOOTH_INTEGALS(MESH, TRI_OBS, TRI_SRC, ...) 中
        ! I2 由 TRI_SRC（源）坐标加权，I3 由 TRI_OBS（场）坐标加权。
        ! 因此源对顶点 R_OPP_S 应与 I3（场加权）配对，场对顶点 R_OPP_O 应与 I2（源加权）配对。
        P_SMOOTH = DOT_PRODUCT(R_OPP_O, R_OPP_S) * I1_SMOOTH &
                   - DOT_PRODUCT(R_OPP_S, I3_SMOOTH) &
                   - DOT_PRODUCT(R_OPP_O, I2_SMOOTH) + I4_SMOOTH
        Q_SMOOTH = -4.0 * I1_SMOOTH

        Z_PAIR = C_O * C_S * (JW_MU * (P_SING + P_SMOOTH) + J_OVER_WE * (Q_SING + Q_SMOOTH))
    END SUBROUTINE CALC_EFIE_EXTRACTED_PAIR

    ! 递归静态积分辅助函数
    RECURSIVE SUBROUTINE INTEGRATE_STATIC_RECURSIVE(V1, V2, V3, TRI_SRC_V, &
                                                     LEVEL, MAX_LEVEL, GDATA, &
                                                     ACC_I0, ACC_RI0, ACC_IVEC, ACC_RDOT)
        REAL, INTENT(IN) :: V1(3), V2(3), V3(3)
        REAL, INTENT(IN) :: TRI_SRC_V(3, 3)
        INTEGER, INTENT(IN) :: LEVEL, MAX_LEVEL
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        REAL, INTENT(INOUT) :: ACC_I0, ACC_RI0(3), ACC_IVEC(3), ACC_RDOT

        REAL :: SUB_V(3, 3, 4), AREA_SUB, CROSS(3), R_VEC(3)
        REAL :: I0_VAL, C1, C2, C3, RESULT_3(3), IVEC(3)
        INTEGER :: I, K, P

        IF (LEVEL >= MAX_LEVEL) THEN
            ! 叶节点：计算当前子三角形面积
            CROSS(1) = (V2(2)-V1(2))*(V3(3)-V1(3)) - (V2(3)-V1(3))*(V3(2)-V1(2))
            CROSS(2) = (V2(3)-V1(3))*(V3(1)-V1(1)) - (V2(1)-V1(1))*(V3(3)-V1(3))
            CROSS(3) = (V2(1)-V1(1))*(V3(2)-V1(2)) - (V2(2)-V1(2))*(V3(1)-V1(1))
            AREA_SUB = 0.5 * SQRT(CROSS(1)**2 + CROSS(2)**2 + CROSS(3)**2)

            DO P = 1, GDATA%N_POINTS
                DO K = 1, 3
                    R_VEC(K) = V1(K) * GDATA%UVW(1, P) + &
                               V2(K) * GDATA%UVW(2, P) + &
                               V3(K) * GDATA%UVW(3, P)
                END DO

                CALL ANALYTIC_SCALAR_POT_1_OVER_R(R_VEC, TRI_SRC_V(:,1), &
                                                   TRI_SRC_V(:,2), TRI_SRC_V(:,3), I0_VAL)
                CALL ANALYTIC_LINEAR_POT_1_OVER_R(R_VEC, TRI_SRC_V(:,1), &
                                                   TRI_SRC_V(:,2), TRI_SRC_V(:,3), RESULT_3)
                C1 = RESULT_3(1); C2 = RESULT_3(2); C3 = RESULT_3(3)
                DO K = 1, 3
                    IVEC(K) = TRI_SRC_V(K,1)*C1 + TRI_SRC_V(K,2)*C2 + TRI_SRC_V(K,3)*C3
                END DO

                ACC_I0  = ACC_I0  + GDATA%WEIGHTS(P) * I0_VAL * AREA_SUB
                ACC_RI0 = ACC_RI0 + GDATA%WEIGHTS(P) * R_VEC * I0_VAL * AREA_SUB
                ACC_IVEC = ACC_IVEC + GDATA%WEIGHTS(P) * IVEC * AREA_SUB
                ACC_RDOT = ACC_RDOT + GDATA%WEIGHTS(P) * DOT_PRODUCT(R_VEC, IVEC) * AREA_SUB
            END DO
        ELSE
            ! 细分
            CALL SUBDIVIDE_TRI_4(V1, V2, V3, SUB_V)
            DO I = 1, 4
                CALL INTEGRATE_STATIC_RECURSIVE(SUB_V(:,1,I), SUB_V(:,2,I), SUB_V(:,3,I), &
                                                 TRI_SRC_V, LEVEL+1, MAX_LEVEL, GDATA, &
                                                 ACC_I0, ACC_RI0, ACC_IVEC, ACC_RDOT)
            END DO
        END IF
    END SUBROUTINE INTEGRATE_STATIC_RECURSIVE

    ! ================ 主入口：五类分派 ================

    SUBROUTINE CALC_EFIE_MATRIX_ELEMENT(MESH, RWG_M, RWG_N, GDATA, K, ETA0, Z_MN)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        TYPE(RWG_BASIS), INTENT(IN) :: RWG_M, RWG_N
        REAL, INTENT(IN) :: K, ETA0
        COMPLEX, INTENT(OUT) :: Z_MN

        INTEGER :: TRI_M(2), TRI_N(2)
        REAL :: COEF_M(2), COEF_N(2)
        REAL :: R_OPP_M(3, 2), R_OPP_N(3, 2)
        INTEGER :: I, J, N_COMMON
        REAL :: D_MIN, H_O, H_S, H_PAIR, CHI
        COMPLEX :: JW_MU, J_OVER_WE
        COMPLEX :: Z_TERM, P_TERM, Q_TERM
        COMPLEX :: I1, I4
        COMPLEX :: I2(3), I3(3)
        LOGICAL :: ADM_OK

        JW_MU = (0.0, 1.0) * ETA0 * K
        J_OVER_WE = (0.0, 1.0) * ETA0 / K
        Z_MN = (0.0, 0.0)

        TRI_M(1) = RWG_M%POS_TRI_ID; TRI_M(2) = RWG_M%NEG_TRI_ID
        TRI_N(1) = RWG_N%POS_TRI_ID; TRI_N(2) = RWG_N%NEG_TRI_ID

        COEF_M(1) = RWG_M%POS_COEF;  COEF_M(2) = -RWG_M%NEG_COEF
        COEF_N(1) = RWG_N%POS_COEF;  COEF_N(2) = -RWG_N%NEG_COEF

        R_OPP_M(1, 1) = MESH%NODES(RWG_M%POS_OPP_VERTEX)%X
        R_OPP_M(2, 1) = MESH%NODES(RWG_M%POS_OPP_VERTEX)%Y
        R_OPP_M(3, 1) = MESH%NODES(RWG_M%POS_OPP_VERTEX)%Z
        R_OPP_M(1, 2) = MESH%NODES(RWG_M%NEG_OPP_VERTEX)%X
        R_OPP_M(2, 2) = MESH%NODES(RWG_M%NEG_OPP_VERTEX)%Y
        R_OPP_M(3, 2) = MESH%NODES(RWG_M%NEG_OPP_VERTEX)%Z
        R_OPP_N(1, 1) = MESH%NODES(RWG_N%POS_OPP_VERTEX)%X
        R_OPP_N(2, 1) = MESH%NODES(RWG_N%POS_OPP_VERTEX)%Y
        R_OPP_N(3, 1) = MESH%NODES(RWG_N%POS_OPP_VERTEX)%Z
        R_OPP_N(1, 2) = MESH%NODES(RWG_N%NEG_OPP_VERTEX)%X
        R_OPP_N(2, 2) = MESH%NODES(RWG_N%NEG_OPP_VERTEX)%Y
        R_OPP_N(3, 2) = MESH%NODES(RWG_N%NEG_OPP_VERTEX)%Z

        DO I = 1, 2  ! 源 RWG 三角形
            DO J = 1, 2  ! 场 RWG 三角形
                N_COMMON = COUNT_COMMON_VERTICES(MESH, TRI_M(I), TRI_N(J))

                SELECT CASE (N_COMMON)
                CASE (3)
                    ! ---- 第一类：同一三角形 ----
                    IF (EFIE_NEAR_METHOD == EFIE_NEAR_QBX) THEN
                        CALL CALC_EFIE_QBX_PAIR(MESH, TRI_N(J), TRI_M(I), &
                                                RWG_N, RWG_M, K, ETA0, Z_TERM, ADM_OK)
                        IF (.NOT. ADM_OK) THEN
                            PRINT *, "警告：自项 QBX 失败，回退到奇异提取法"
                            CALL CALC_EFIE_SELF_TRI_PAIR(MESH, TRI_M(I), RWG_M, RWG_N, K, ETA0, Z_TERM)
                        END IF
                    ELSE IF (EFIE_NEAR_METHOD == EFIE_NEAR_DIRECTFN) THEN
                        CALL CALC_EFIE_DIRECTFN_PAIR(MESH, TRI_N(J), TRI_M(I), &
                                                     RWG_N, RWG_M, K, ETA0, Z_TERM)
                    ELSE
                        CALL CALC_EFIE_SELF_TRI_PAIR(MESH, TRI_M(I), RWG_M, RWG_N, K, ETA0, Z_TERM)
                    END IF
                    Z_MN = Z_MN + Z_TERM

                CASE (2)
                    ! ---- 第二类：共边三角形 ----
                    IF (EFIE_NEAR_METHOD == EFIE_NEAR_QBX) THEN
                        CALL CALC_EFIE_QBX_PAIR(MESH, TRI_N(J), TRI_M(I), &
                                                RWG_N, RWG_M, K, ETA0, Z_TERM, ADM_OK)
                        IF (.NOT. ADM_OK) THEN
                            PRINT *, "警告：共边对 QBX 失败，回退到奇异提取法"
                            CALL CALC_EFIE_EXTRACTED_PAIR(MESH, TRI_N(J), TRI_M(I), &
                                                           RWG_N, RWG_M, K, ETA0, N_SUBDIV_EDGE, Z_TERM)
                        END IF
                    ELSE IF (EFIE_NEAR_METHOD == EFIE_NEAR_DIRECTFN) THEN
                        CALL CALC_EFIE_DIRECTFN_PAIR(MESH, TRI_N(J), TRI_M(I), &
                                                     RWG_N, RWG_M, K, ETA0, Z_TERM)
                    ELSE
                        CALL CALC_EFIE_EXTRACTED_PAIR(MESH, TRI_N(J), TRI_M(I), &
                                                       RWG_N, RWG_M, K, ETA0, N_SUBDIV_EDGE, Z_TERM)
                    END IF
                    Z_MN = Z_MN + Z_TERM

                CASE (1)
                    ! ---- 第三类：共顶点三角形 ----
                    IF (EFIE_NEAR_METHOD == EFIE_NEAR_QBX) THEN
                        CALL CALC_EFIE_QBX_PAIR(MESH, TRI_N(J), TRI_M(I), &
                                                RWG_N, RWG_M, K, ETA0, Z_TERM, ADM_OK)
                        IF (.NOT. ADM_OK) THEN
                            CALL CALC_EFIE_EXTRACTED_PAIR(MESH, TRI_N(J), TRI_M(I), &
                                                           RWG_N, RWG_M, K, ETA0, N_SUBDIV_VERTEX, Z_TERM)
                        END IF
                    ELSE IF (EFIE_NEAR_METHOD == EFIE_NEAR_DIRECTFN) THEN
                        CALL CALC_EFIE_DIRECTFN_PAIR(MESH, TRI_N(J), TRI_M(I), &
                                                     RWG_N, RWG_M, K, ETA0, Z_TERM)
                    ELSE
                        CALL CALC_EFIE_EXTRACTED_PAIR(MESH, TRI_N(J), TRI_M(I), &
                                                       RWG_N, RWG_M, K, ETA0, N_SUBDIV_VERTEX, Z_TERM)
                    END IF
                    Z_MN = Z_MN + Z_TERM

                CASE (0)
                    ! ---- 第四/五类：按距离区分 ----
                    D_MIN = TRI_CENTROID_DIST(MESH, TRI_M(I), TRI_N(J))
                    H_O   = TRI_MAX_EDGE_LEN(MESH, TRI_M(I))
                    H_S   = TRI_MAX_EDGE_LEN(MESH, TRI_N(J))
                    H_PAIR = MAX(H_O, H_S)
                    CHI = D_MIN / H_PAIR

                    IF (CHI < CHI_NEAR) THEN
                        ! 第四类：非接触近场
                        IF (EFIE_NEAR_METHOD == EFIE_NEAR_QBX) THEN
                            CALL CALC_EFIE_QBX_PAIR(MESH, TRI_N(J), TRI_M(I), &
                                                    RWG_N, RWG_M, K, ETA0, Z_TERM, ADM_OK)
                            IF (.NOT. ADM_OK) THEN
                                CALL CALC_EFIE_EXTRACTED_PAIR(MESH, TRI_N(J), TRI_M(I), &
                                                               RWG_N, RWG_M, K, ETA0, N_SUBDIV_NEAR, Z_TERM)
                            END IF
                        ELSE
                            CALL CALC_EFIE_EXTRACTED_PAIR(MESH, TRI_N(J), TRI_M(I), &
                                                           RWG_N, RWG_M, K, ETA0, N_SUBDIV_NEAR, Z_TERM)
                        END IF
                    ELSE
                        ! 第五类：普通远场 → 完整格林函数
                        CALL CALC_GREEN_INTEGALS(MESH, TRI_M(I), TRI_N(J), GDATA, K, &
                                                  I1, I2, I3, I4)
                        P_TERM = DOT_PRODUCT(R_OPP_M(:, I), R_OPP_N(:, J)) * I1 &
                                 - DOT_PRODUCT(R_OPP_N(:, J), I3) &
                                 - DOT_PRODUCT(R_OPP_M(:, I), I2) + I4
                        Q_TERM = -4.0 * I1
                        Z_TERM = COEF_M(I) * COEF_N(J) * (JW_MU * P_TERM + J_OVER_WE * Q_TERM)
                    END IF
                    Z_MN = Z_MN + Z_TERM
                END SELECT
            END DO
        END DO
    END SUBROUTINE CALC_EFIE_MATRIX_ELEMENT

    ! ================ MFIE ================

    ! 判断 RWG_M 与 RWG_N 是否共享同一个三角形 TRI_ID
    FUNCTION RWG_SHARE_TRIANGLE(RWG_M, RWG_N, TRI_ID) RESULT(SHARE)
        TYPE(RWG_BASIS), INTENT(IN) :: RWG_M, RWG_N
        INTEGER, INTENT(IN) :: TRI_ID
        LOGICAL :: SHARE
        SHARE = (TRI_ID == RWG_M%POS_TRI_ID .AND. TRI_ID == RWG_N%POS_TRI_ID) .OR. &
                (TRI_ID == RWG_M%POS_TRI_ID .AND. TRI_ID == RWG_N%NEG_TRI_ID) .OR. &
                (TRI_ID == RWG_M%NEG_TRI_ID .AND. TRI_ID == RWG_N%POS_TRI_ID) .OR. &
                (TRI_ID == RWG_M%NEG_TRI_ID .AND. TRI_ID == RWG_N%NEG_TRI_ID)
    END FUNCTION RWG_SHARE_TRIANGLE

    ! 计算单个三角形上 RWG 基函数的 jump 项：1/2 ∫ f_m·f_n dS
    SUBROUTINE RWG_JUMP_TERM(MESH, TRI_ID, RWG_M, RWG_N, GDATA, JUMP)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_ID
        TYPE(RWG_BASIS), INTENT(IN) :: RWG_M, RWG_N
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        COMPLEX, INTENT(OUT) :: JUMP

        CHARACTER * 3 :: TRI_M_LOCAL, TRI_N_LOCAL
        REAL :: F_VAL_M(3), F_VAL_N(3), GLOBAL_PTS(3, GDATA%N_POINTS)
        REAL :: DOT_TERM
        INTEGER :: I

        IF (RWG_M%POS_TRI_ID == TRI_ID) THEN
            TRI_M_LOCAL = 'POS'
        ELSE
            TRI_M_LOCAL = 'NEG'
        END IF
        IF (RWG_N%POS_TRI_ID == TRI_ID) THEN
            TRI_N_LOCAL = 'POS'
        ELSE
            TRI_N_LOCAL = 'NEG'
        END IF

        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_ID, GDATA, GLOBAL_PTS)
        DOT_TERM = 0.0
        DO I = 1, GDATA%N_POINTS
            CALL EVAL_RWG_BASIS(MESH, RWG_M, TRI_M_LOCAL, GLOBAL_PTS(:, I), F_VAL_M)
            CALL EVAL_RWG_BASIS(MESH, RWG_N, TRI_N_LOCAL, GLOBAL_PTS(:, I), F_VAL_N)
            DOT_TERM = DOT_TERM + GDATA%WEIGHTS(I) * DOT_PRODUCT(F_VAL_M, F_VAL_N)
        END DO
        JUMP = CMPLX(0.5 * DOT_TERM * MESH%TRIANGLES(TRI_ID)%AREA)
    END SUBROUTINE RWG_JUMP_TERM

    ! 递归细分求源三角形上 ∇G 正则余项的三矩
    RECURSIVE SUBROUTINE INTEGRATE_GRAD_REG_RECURSIVE(V1, V2, V3, R_I, N_O, K, &
                                                       GDATA, LEVEL, MAX_LEVEL, &
                                                       G_A_REG, G_C_REG, G_E_REG)
        REAL, INTENT(IN) :: V1(3), V2(3), V3(3), R_I(3), N_O(3), K
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        INTEGER, INTENT(IN) :: LEVEL, MAX_LEVEL
        COMPLEX, INTENT(INOUT) :: G_A_REG(3), G_C_REG(3), G_E_REG(3)

        REAL :: SUB_V(3, 3, 4), CROSS(3), AREA_SUB, R_J(3)
        COMPLEX :: GRAD_REG(3)
        INTEGER :: I, KIDX

        IF (LEVEL >= MAX_LEVEL) THEN
            CROSS(1) = (V2(2)-V1(2))*(V3(3)-V1(3)) - (V2(3)-V1(3))*(V3(2)-V1(2))
            CROSS(2) = (V2(3)-V1(3))*(V3(1)-V1(1)) - (V2(1)-V1(1))*(V3(3)-V1(3))
            CROSS(3) = (V2(1)-V1(1))*(V3(2)-V1(2)) - (V2(2)-V1(2))*(V3(1)-V1(1))
            AREA_SUB = 0.5 * SQRT(SUM(CROSS**2))

            DO I = 1, GDATA%N_POINTS
                DO KIDX = 1, 3
                    R_J(KIDX) = V1(KIDX)*GDATA%UVW(1,I) + V2(KIDX)*GDATA%UVW(2,I) + &
                                V3(KIDX)*GDATA%UVW(3,I)
                END DO
                CALL GRAD_GREEN_REG(R_I, R_J, K, GRAD_REG)
                G_A_REG = G_A_REG + GDATA%WEIGHTS(I) * AREA_SUB * GRAD_REG
                G_C_REG = G_C_REG + GDATA%WEIGHTS(I) * AREA_SUB * R_J * DOT_PRODUCT(N_O, GRAD_REG)
                G_E_REG = G_E_REG + GDATA%WEIGHTS(I) * AREA_SUB * DOT_PRODUCT(N_O, R_J) * GRAD_REG
            END DO
        ELSE
            CALL SUBDIVIDE_TRI_4(V1, V2, V3, SUB_V)
            DO I = 1, 4
                CALL INTEGRATE_GRAD_REG_RECURSIVE(SUB_V(:,1,I), SUB_V(:,2,I), SUB_V(:,3,I), &
                                                   R_I, N_O, K, GDATA, LEVEL+1, MAX_LEVEL, &
                                                   G_A_REG, G_C_REG, G_E_REG)
            END DO
        END IF
    END SUBROUTINE INTEGRATE_GRAD_REG_RECURSIVE

    ! MFIE 近场三角形对：静态主值（解析）+ 正则余项（数值高斯）
    ! 参考 docs/MFIE_CFIE实现文档_内谐振对策.md §4-§6
    SUBROUTINE CALC_MFIE_NEAR_TRI_PAIR(MESH, TRI_OBS, TRI_SRC, RWG_M, RWG_N, &
                                        GDATA, K, Z_PAIR)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_OBS, TRI_SRC
        TYPE(RWG_BASIS), INTENT(IN) :: RWG_M, RWG_N
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        REAL, INTENT(IN) :: K
        COMPLEX, INTENT(OUT) :: Z_PAIR

        REAL :: C_O, C_S, R_OPP_O(3), R_OPP_S(3)
        REAL :: N_O(3), A_O, A_S
        REAL :: V_SRC(3, 3), GLOBAL_PTS_O(3, GDATA%N_POINTS), GLOBAL_PTS_S(3, GDATA%N_POINTS)
        REAL :: R_I(3), R_J(3)
        REAL :: G_A_ST(3), G_C_ST(3), G_E_ST(3)
        COMPLEX :: G_A_REG(3), G_C_REG(3), G_E_REG(3), GRAD_REG(3)
        COMPLEX :: G_A(3), G_C(3), G_E(3)
        COMPLEX :: G_B
        COMPLEX :: P1(3), P3(3), P5(3), P8(3)
        COMPLEX :: P2, P4, P6, P7
        REAL :: W_O
        INTEGER :: I, J, KIDX, N_COMMON_INNER, N_SUB_INNER
        COMPLEX :: A_TERM, B_TERM, PV, JUMP
        REAL :: NF_DOT_ROPP_S

        ! 提取系数与对顶点
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

        N_O = MESH%TRIANGLES(TRI_OBS)%NORMAL
        A_O = MESH%TRIANGLES(TRI_OBS)%AREA
        A_S = MESH%TRIANGLES(TRI_SRC)%AREA

        DO KIDX = 1, 3
            V_SRC(KIDX, 1) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(1))%X
            V_SRC(KIDX, 2) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(2))%X
            V_SRC(KIDX, 3) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(3))%X
            IF (KIDX == 2) THEN
                V_SRC(KIDX, 1) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(1))%Y
                V_SRC(KIDX, 2) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(2))%Y
                V_SRC(KIDX, 3) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(3))%Y
            ELSE IF (KIDX == 3) THEN
                V_SRC(KIDX, 1) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(1))%Z
                V_SRC(KIDX, 2) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(2))%Z
                V_SRC(KIDX, 3) = MESH%NODES(MESH%TRIANGLES(TRI_SRC)%VERTEX_3D(3))%Z
            END IF
        END DO

        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_OBS, GDATA, GLOBAL_PTS_O)
        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_SRC, GDATA, GLOBAL_PTS_S)

        ! 安全分支：同一三角形仅 jump；共面相邻对 PV=0
        IF (TRI_OBS == TRI_SRC) THEN
            Z_PAIR = (0.0, 0.0)
            IF (RWG_SHARE_TRIANGLE(RWG_M, RWG_N, TRI_OBS)) THEN
                CALL RWG_JUMP_TERM(MESH, TRI_OBS, RWG_M, RWG_N, GDATA, JUMP)
                Z_PAIR = Z_PAIR + JUMP
            END IF
            RETURN
        END IF
        IF (IS_COPLANAR_TRI_PAIR(MESH, TRI_OBS, TRI_SRC, 1.0E-6, 1.0E-6)) THEN
            Z_PAIR = (0.0, 0.0)
            RETURN
        END IF

        P1 = (0.0, 0.0); P2 = (0.0, 0.0); P3 = (0.0, 0.0); P4 = (0.0, 0.0)
        P5 = (0.0, 0.0); P6 = (0.0, 0.0); P7 = (0.0, 0.0); P8 = (0.0, 0.0)

        DO I = 1, GDATA%N_POINTS
            R_I = GLOBAL_PTS_O(:, I)
            W_O = GDATA%WEIGHTS(I)

            ! 解析静态主值矩（不含 1/(4π)）
            CALL ANALYTIC_GRAD_STATIC_PV(R_I, N_O, V_SRC(:,1), V_SRC(:,2), V_SRC(:,3), &
                                          G_A_ST, G_C_ST, G_E_ST)

            ! 数值正则余项矩（接触对细分以提高近奇异性精度）
            N_COMMON_INNER = COUNT_COMMON_VERTICES(MESH, TRI_OBS, TRI_SRC)
            IF (N_COMMON_INNER >= 2) THEN
                N_SUB_INNER = 2
            ELSE IF (N_COMMON_INNER == 1) THEN
                N_SUB_INNER = 1
            ELSE
                N_SUB_INNER = 0
            END IF

            G_A_REG = (0.0, 0.0); G_C_REG = (0.0, 0.0); G_E_REG = (0.0, 0.0)
            CALL INTEGRATE_GRAD_REG_RECURSIVE(V_SRC(:,1), V_SRC(:,2), V_SRC(:,3), &
                                               R_I, N_O, K, GDATA, 0, N_SUB_INNER, &
                                               G_A_REG, G_C_REG, G_E_REG)

            ! 合成完整矩（含 1/(4π)）
            G_A = CMPLX(G_A_ST) / (4.0 * PI) + G_A_REG
            G_C = CMPLX(G_C_ST) / (4.0 * PI) + G_C_REG
            G_E = CMPLX(G_E_ST) / (4.0 * PI) + G_E_REG
            G_B = DOT_PRODUCT(N_O, G_A)

            ! 累加外层矩（乘 A_O 与权重）
            P1 = P1 + W_O * A_O * G_A
            P2 = P2 + W_O * A_O * DOT_PRODUCT(R_I, CMPLX(G_A_ST)/(4.0*PI) + G_A_REG)
            P3 = P3 + W_O * A_O * G_E
            P4 = P4 + W_O * A_O * DOT_PRODUCT(R_I, G_E)
            P5 = P5 + W_O * A_O * G_C
            P6 = P6 + W_O * A_O * DOT_PRODUCT(R_I, G_C)
            P7 = P7 + W_O * A_O * CMPLX(G_B)
            P8 = P8 + W_O * A_O * R_I * CMPLX(G_B)
        END DO

        ! §6.2 组装
        NF_DOT_ROPP_S = DOT_PRODUCT(N_O, R_OPP_S)
        A_TERM = P4 - NF_DOT_ROPP_S * P2 - DOT_PRODUCT(R_OPP_O, P3) + &
                 NF_DOT_ROPP_S * DOT_PRODUCT(R_OPP_O, P1)
        B_TERM = P6 - DOT_PRODUCT(R_OPP_S, P8) - DOT_PRODUCT(R_OPP_O, P5) + &
                 DOT_PRODUCT(R_OPP_O, R_OPP_S) * P7

        PV = -CMPLX(C_O * C_S) * (A_TERM - B_TERM)

        ! jump 项：若两 RWG 共享 TRI_OBS
        Z_PAIR = PV
        IF (RWG_SHARE_TRIANGLE(RWG_M, RWG_N, TRI_OBS)) THEN
            CALL RWG_JUMP_TERM(MESH, TRI_OBS, RWG_M, RWG_N, GDATA, JUMP)
            Z_PAIR = Z_PAIR + JUMP
        END IF

    END SUBROUTINE CALC_MFIE_NEAR_TRI_PAIR

    SUBROUTINE CALC_MFIE_MATRIX_ELEMENT(MESH, RWG_M, RWG_N, GDATA, K, Z_MN)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        TYPE(RWG_BASIS), INTENT(IN) :: RWG_M, RWG_N
        REAL, INTENT(IN) :: K
        COMPLEX, INTENT(OUT) :: Z_MN
        INTEGER :: TRI_M(2), TRI_N(2)
        REAL :: COEF_M(2), COEF_N(2)
        REAL :: R_OPP_M(3, 2), R_OPP_N(3, 2)
        INTEGER :: M, N, TRI_M_ID, TRI_N_ID, N_COMMON
        REAL :: C_M, C_N, R_M(3), R_N(3), N_M(3), A_M, A_N
        REAL :: D_MIN, H_O, H_S, H_PAIR, CHI
        REAL :: GLOBAL_PTS_M(3, GDATA%N_POINTS), GLOBAL_PTS_N(3, GDATA%N_POINTS)
        REAL :: R_I(3), R_J(3), F_M(3), F_N(3), R_VEC(3), DIST
        COMPLEX :: Z_PAIR, K_ACC(3), GRAD_G(3), CROSS(3), NCK(3), DOT
        REAL :: F_VAL_M(3), F_VAL_N(3), DOT_FM_FN, HALF_TERM
        CHARACTER * 3 :: TRI_M_LOCAL, TRI_N_LOCAL
        INTEGER :: I, J_MFI

        TRI_M(1) = RWG_M%POS_TRI_ID; TRI_M(2) = RWG_M%NEG_TRI_ID
        TRI_N(1) = RWG_N%POS_TRI_ID; TRI_N(2) = RWG_N%NEG_TRI_ID
        COEF_M(1) = RWG_M%POS_COEF; COEF_M(2) = -RWG_M%NEG_COEF
        COEF_N(1) = RWG_N%POS_COEF; COEF_N(2) = -RWG_N%NEG_COEF

        R_OPP_M(1, 1) = MESH%NODES(RWG_M%POS_OPP_VERTEX)%X
        R_OPP_M(2, 1) = MESH%NODES(RWG_M%POS_OPP_VERTEX)%Y
        R_OPP_M(3, 1) = MESH%NODES(RWG_M%POS_OPP_VERTEX)%Z
        R_OPP_M(1, 2) = MESH%NODES(RWG_M%NEG_OPP_VERTEX)%X
        R_OPP_M(2, 2) = MESH%NODES(RWG_M%NEG_OPP_VERTEX)%Y
        R_OPP_M(3, 2) = MESH%NODES(RWG_M%NEG_OPP_VERTEX)%Z
        R_OPP_N(1, 1) = MESH%NODES(RWG_N%POS_OPP_VERTEX)%X
        R_OPP_N(2, 1) = MESH%NODES(RWG_N%POS_OPP_VERTEX)%Y
        R_OPP_N(3, 1) = MESH%NODES(RWG_N%POS_OPP_VERTEX)%Z
        R_OPP_N(1, 2) = MESH%NODES(RWG_N%NEG_OPP_VERTEX)%X
        R_OPP_N(2, 2) = MESH%NODES(RWG_N%NEG_OPP_VERTEX)%Y
        R_OPP_N(3, 2) = MESH%NODES(RWG_N%NEG_OPP_VERTEX)%Z

        Z_MN = (0.0, 0.0)

        DO M = 1, 2
            DO N = 1, 2
                TRI_M_ID = TRI_M(M); TRI_N_ID = TRI_N(N)
                C_M = COEF_M(M); C_N = COEF_N(N)
                R_M = R_OPP_M(:, M); R_N = R_OPP_N(:, N)
                N_M = MESH%TRIANGLES(TRI_M_ID)%NORMAL
                A_M = MESH%TRIANGLES(TRI_M_ID)%AREA
                A_N = MESH%TRIANGLES(TRI_N_ID)%AREA
                CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_M_ID, GDATA, GLOBAL_PTS_M)
                CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_N_ID, GDATA, GLOBAL_PTS_N)

                IF (TRI_M_ID == TRI_N_ID) THEN
                    HALF_TERM = 0.0
                    IF (RWG_M%POS_TRI_ID == TRI_M_ID) THEN
                        TRI_M_LOCAL = 'POS'
                    ELSE
                        TRI_M_LOCAL = 'NEG'
                    END IF
                    IF (RWG_N%POS_TRI_ID == TRI_N_ID) THEN
                        TRI_N_LOCAL = "POS"
                    ELSE
                        TRI_N_LOCAL = 'NEG'
                    END IF
                    DO I = 1, GDATA%N_POINTS
                        CALL EVAL_RWG_BASIS(MESH, RWG_M, TRI_M_LOCAL, GLOBAL_PTS_M(:, I), F_VAL_M)
                        CALL EVAL_RWG_BASIS(MESH, RWG_N, TRI_N_LOCAL, GLOBAL_PTS_N(:, I), F_VAL_N)
                        DOT_FM_FN = DOT_PRODUCT(F_VAL_M, F_VAL_N)
                        HALF_TERM = HALF_TERM + GDATA%WEIGHTS(I) * DOT_FM_FN
                    END DO
                    HALF_TERM = 0.5 * HALF_TERM * MESH%TRIANGLES(TRI_M_ID)%AREA
                    Z_MN = CMPLX(HALF_TERM) + Z_MN
                    CYCLE
                END IF

                IF (IS_COPLANAR_TRI_PAIR(MESH, TRI_M_ID, TRI_N_ID, 1.0E-6, 1.0E-6)) CYCLE

                N_COMMON = COUNT_COMMON_VERTICES(MESH, TRI_M_ID, TRI_N_ID)

                IF (N_COMMON >= 1) THEN
                    ! 共边 / 共顶点：解析 PV + 细分数值正则余项
                    CALL CALC_MFIE_NEAR_TRI_PAIR(MESH, TRI_M_ID, TRI_N_ID, &
                                                  RWG_M, RWG_N, GDATA, K, Z_PAIR)
                ELSE
                    D_MIN = TRI_CENTROID_DIST(MESH, TRI_M_ID, TRI_N_ID)
                    H_O   = TRI_MAX_EDGE_LEN(MESH, TRI_M_ID)
                    H_S   = TRI_MAX_EDGE_LEN(MESH, TRI_N_ID)
                    H_PAIR = MAX(H_O, H_S)
                    CHI = D_MIN / H_PAIR

                    IF (CHI < CHI_NEAR) THEN
                        CALL CALC_MFIE_NEAR_TRI_PAIR(MESH, TRI_M_ID, TRI_N_ID, &
                                                      RWG_M, RWG_N, GDATA, K, Z_PAIR)
                    ELSE
                        ! 远场：完整 ∇G 数值积分
                        Z_PAIR = 0.0
                        DO I = 1, GDATA%N_POINTS
                            R_I = GLOBAL_PTS_M(:, I)
                            F_M = C_M * (R_I - R_M)
                            K_ACC = (0.0, 0.0)
                            DO J_MFI = 1, GDATA%N_POINTS
                                R_J = GLOBAL_PTS_N(:, J_MFI)
                                F_N = C_N * (R_J - R_N)
                                R_VEC = R_I - R_J
                                DIST = SQRT(SUM(R_VEC ** 2))
                                IF (DIST > 1.0E-10) THEN
                                    CALL GARD_GREEN_FUNC(R_I, R_J, K, GRAD_G)
                                    CROSS(1) = F_N(2) * GRAD_G(3) - F_N(3) * GRAD_G(2)
                                    CROSS(2) = F_N(3) * GRAD_G(1) - F_N(1) * GRAD_G(3)
                                    CROSS(3) = F_N(1) * GRAD_G(2) - F_N(2) * GRAD_G(1)
                                    K_ACC = K_ACC + GDATA%WEIGHTS(J_MFI) * CROSS
                                END IF
                            END DO
                            K_ACC = K_ACC * A_N
                            NCK(1) = N_M(2) * K_ACC(3) - N_M(3) * K_ACC(2)
                            NCK(2) = N_M(3) * K_ACC(1) - N_M(1) * K_ACC(3)
                            NCK(3) = N_M(1) * K_ACC(2) - N_M(2) * K_ACC(1)
                            DOT = F_M(1) * NCK(1) + F_M(2) * NCK(2) + F_M(3) * NCK(3)
                            Z_PAIR = Z_PAIR + GDATA%WEIGHTS(I) * DOT
                        END DO
                        Z_PAIR = Z_PAIR * A_M
                    END IF
                END IF
                Z_MN = Z_MN + Z_PAIR
            END DO
        END DO
    END SUBROUTINE CALC_MFIE_MATRIX_ELEMENT

    ! 判断两三角形是否共面
    FUNCTION IS_COPLANAR_TRI_PAIR(MESH, TRI_M_ID, TRI_N_ID, TOL_NORM, TOL_DIST) RESULT(IS_COPLANAR)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_M_ID, TRI_N_ID
        REAL, INTENT(IN) :: TOL_NORM, TOL_DIST
        REAL :: DOT_N, DOT_G, R_VEC(3)
        LOGICAL :: IS_COPLANAR

        IS_COPLANAR = .TRUE.
        DOT_N = DOT_PRODUCT(MESH%TRIANGLES(TRI_M_ID)%NORMAL, MESH%TRIANGLES(TRI_N_ID)%NORMAL)
        IF (ABS(ABS(DOT_N) - 1.0) > TOL_NORM) IS_COPLANAR = .FALSE.
        R_VEC = MESH%TRIANGLES(TRI_M_ID)%CENTROID - MESH%TRIANGLES(TRI_N_ID)%CENTROID
        DOT_G = DOT_PRODUCT(MESH%TRIANGLES(TRI_M_ID)%NORMAL, R_VEC)
        IF (ABS(ABS(DOT_G)) > TOL_DIST) IS_COPLANAR = .FALSE.
    END FUNCTION

    ! ================ 自项（保持独立，复用提取逻辑） ================

    SUBROUTINE CALC_EFIE_SELF_TRI_PAIR(MESH, TRI_ID, RWG_M, RWG_N, K, ETA0, Z_PAIR)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_ID
        TYPE(RWG_BASIS), INTENT(IN) :: RWG_M, RWG_N
        REAL, INTENT(IN) :: ETA0, K
        COMPLEX, INTENT(OUT) :: Z_PAIR

        REAL :: C_M, C_N
        REAL :: R_OPP_M(3), R_OPP_N(3)
        COMPLEX :: JW_MU, J_OVER_WE
        COMPLEX :: P_SING, Q_SING, P_SMOOTH2, Q_SMOOTH2
        REAL :: I0_VAL, I1_VAL, I2_VAL, I3_VAL, IVEC_VAL(3)
        REAL :: SUM_I0, SUM_R_I0(3), SUM_IVEC(3), SUM_R_DOT_IVEC
        REAL :: J0, J_R_I0(3), J_VEC(3), J_R_DOT
        INTEGER :: IC, KC
        REAL :: R_VEC(3), GLOBAL_PTS(3, 12)
        REAL :: V1(3), V2(3), V3(3), RESULT_3(3)
        COMPLEX :: I1, I4, I2(3), I3(3)
        TYPE(GAUSS_TRI_DATA) :: GDATA12

        IF (TRI_ID == RWG_M%POS_TRI_ID) THEN
            C_M = RWG_M%POS_COEF
            R_OPP_M(1) = MESH%NODES(RWG_M%POS_OPP_VERTEX)%X
            R_OPP_M(2) = MESH%NODES(RWG_M%POS_OPP_VERTEX)%Y
            R_OPP_M(3) = MESH%NODES(RWG_M%POS_OPP_VERTEX)%Z
        ELSE
            C_M = -RWG_M%NEG_COEF
            R_OPP_M(1) = MESH%NODES(RWG_M%NEG_OPP_VERTEX)%X
            R_OPP_M(2) = MESH%NODES(RWG_M%NEG_OPP_VERTEX)%Y
            R_OPP_M(3) = MESH%NODES(RWG_M%NEG_OPP_VERTEX)%Z
        END IF

        IF (TRI_ID == RWG_N%POS_TRI_ID) THEN
            C_N = RWG_N%POS_COEF
            R_OPP_N(1) = MESH%NODES(RWG_N%POS_OPP_VERTEX)%X
            R_OPP_N(2) = MESH%NODES(RWG_N%POS_OPP_VERTEX)%Y
            R_OPP_N(3) = MESH%NODES(RWG_N%POS_OPP_VERTEX)%Z
        ELSE
            C_N = -RWG_N%NEG_COEF
            R_OPP_N(1) = MESH%NODES(RWG_N%NEG_OPP_VERTEX)%X
            R_OPP_N(2) = MESH%NODES(RWG_N%NEG_OPP_VERTEX)%Y
            R_OPP_N(3) = MESH%NODES(RWG_N%NEG_OPP_VERTEX)%Z
        END IF

        JW_MU = (0.0, 1.0) * ETA0 * K
        J_OVER_WE = (0.0, 1.0) * ETA0 / K

        SUM_I0 = 0.0; SUM_R_I0 = 0.0; SUM_IVEC = 0.0; SUM_R_DOT_IVEC = 0.0
        CALL INIT_GAUSS_TRI(GAUSS_12PT, GDATA12)
        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_ID, GDATA12, GLOBAL_PTS)

        V1(1) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(1))%X
        V1(2) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(1))%Y
        V1(3) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(1))%Z
        V2(1) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(2))%X
        V2(2) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(2))%Y
        V2(3) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(2))%Z
        V3(1) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(3))%X
        V3(2) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(3))%Y
        V3(3) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(3))%Z

        DO IC = 1, GDATA12%N_POINTS
            R_VEC = GLOBAL_PTS(:, IC)
            CALL ANALYTIC_SCALAR_POT_1_OVER_R(R_VEC, V1, V2, V3, I0_VAL)
            CALL ANALYTIC_LINEAR_POT_1_OVER_R(R_VEC, V1, V2, V3, RESULT_3)
            I1_VAL = RESULT_3(1); I2_VAL = RESULT_3(2); I3_VAL = RESULT_3(3)
            DO KC = 1, 3
                IVEC_VAL(KC) = V1(KC)*I1_VAL + V2(KC)*I2_VAL + V3(KC)*I3_VAL
            END DO
            SUM_I0 = SUM_I0 + GDATA12%WEIGHTS(IC) * I0_VAL
            SUM_R_I0 = SUM_R_I0 + GDATA12%WEIGHTS(IC) * R_VEC * I0_VAL
            SUM_IVEC = SUM_IVEC + GDATA12%WEIGHTS(IC) * IVEC_VAL
            SUM_R_DOT_IVEC = SUM_R_DOT_IVEC + GDATA12%WEIGHTS(IC) * DOT_PRODUCT(R_VEC, IVEC_VAL)
        END DO

        J0     = SUM_I0 * MESH%TRIANGLES(TRI_ID)%AREA
        J_R_I0 = SUM_R_I0 * MESH%TRIANGLES(TRI_ID)%AREA
        J_VEC  = SUM_IVEC * MESH%TRIANGLES(TRI_ID)%AREA
        J_R_DOT = SUM_R_DOT_IVEC * MESH%TRIANGLES(TRI_ID)%AREA

        Q_SING = -J0 / PI
        P_SING = (J_R_DOT - DOT_PRODUCT(R_OPP_N, J_VEC) &
                  - DOT_PRODUCT(R_OPP_M, J_R_I0) &
                  + DOT_PRODUCT(R_OPP_M, R_OPP_N) * J0) / (4.0 * PI)

        CALL CALC_GREEN_SMOOTH_INTEGALS(MESH, TRI_ID, TRI_ID, GDATA12, K, I1, I2, I3, I4)
        P_SMOOTH2 = DOT_PRODUCT(R_OPP_M, R_OPP_N) * I1 - DOT_PRODUCT(R_OPP_N, I3) &
                    - DOT_PRODUCT(R_OPP_M, I2) + I4
        Q_SMOOTH2 = -4.0 * I1

        Z_PAIR = C_M * C_N * (JW_MU * (P_SING + P_SMOOTH2) + J_OVER_WE * (Q_SING + Q_SMOOTH2))
    END SUBROUTINE CALC_EFIE_SELF_TRI_PAIR

END MODULE Z_MATRIX
