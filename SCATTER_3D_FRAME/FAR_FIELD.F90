MODULE FAR_FIELD
    USE EM_TYPES
    USE NUMERICAL_INTEGRATION
    USE RWG_BASIS_BUILD
    IMPLICIT NONE
CONTAINS

    ! 计算单个 RWG 基函数在远区的辐射积分
    ! F_n = ∫_{S_n} f_n(r') * exp(j*k * r_hat·r') dS'
    ! 输入：MESH, RWG, K, R_HAT（观察方向）, GDATA（高斯积分数据）
    ! 输出：F_VEC(3) — 复数矢量
    SUBROUTINE CALC_RWG_RADIATION_INTEGRAL(MESH, RWG, K, R_HAT, GDATA, F_VEC)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(RWG_BASIS), INTENT(IN) :: RWG
        REAL, INTENT(IN) :: K
        REAL, INTENT(IN) :: R_HAT(3)
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        COMPLEX, INTENT(OUT) :: F_VEC(3)

        INTEGER :: TRI_ID, P, N_PTS, ISUB
        REAL :: GLOBAL_PTS(3, GDATA%N_POINTS)
        REAL :: AREA, R_OPP(3), COEF, DOT_KR, R_PT(3)
        COMPLEX :: PHASE, INTEG(3)
        INTEGER :: N1, N2, N3
        REAL :: V1(3), V2(3), V3(3), SUB_V(3, 3, 4), SUB_AREA
        REAL :: SUB_PTS(3, GDATA%N_POINTS)

        F_VEC = (0.0, 0.0)
        N_PTS = GDATA%N_POINTS

        ! ---- 正三角形 ----
        TRI_ID = RWG%POS_TRI_ID
        AREA = MESH%TRIANGLES(TRI_ID)%AREA
        N1 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(1)
        N2 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(2)
        N3 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(3)
        V1 = [MESH%NODES(N1)%X, MESH%NODES(N1)%Y, MESH%NODES(N1)%Z]
        V2 = [MESH%NODES(N2)%X, MESH%NODES(N2)%Y, MESH%NODES(N2)%Z]
        V3 = [MESH%NODES(N3)%X, MESH%NODES(N3)%Y, MESH%NODES(N3)%Z]
        R_OPP(1) = MESH%NODES(RWG%POS_OPP_VERTEX)%X
        R_OPP(2) = MESH%NODES(RWG%POS_OPP_VERTEX)%Y
        R_OPP(3) = MESH%NODES(RWG%POS_OPP_VERTEX)%Z
        COEF = RWG%POS_COEF

        CALL SUBDIVIDE_TRI_4(V1, V2, V3, SUB_V)
        DO ISUB = 1, 4
            SUB_AREA = AREA / 4.0
            CALL GET_SUB_TRI_GAUSS_POINTS(SUB_V(:,1,ISUB), SUB_V(:,2,ISUB), SUB_V(:,3,ISUB), &
                                           GDATA, SUB_PTS)
            INTEG = (0.0, 0.0)
            DO P = 1, N_PTS
                R_PT = SUB_PTS(:, P)
                DOT_KR = K * DOT_PRODUCT(R_HAT, R_PT)
                PHASE = CEXP((0.0, 1.0) * DOT_KR)
                INTEG = INTEG + GDATA%WEIGHTS(P) * PHASE * (R_PT - R_OPP)
            END DO
            F_VEC = F_VEC + COEF * SUB_AREA * INTEG
        END DO

        ! ---- 负三角形 ----
        TRI_ID = RWG%NEG_TRI_ID
        AREA = MESH%TRIANGLES(TRI_ID)%AREA
        N1 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(1)
        N2 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(2)
        N3 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(3)
        V1 = [MESH%NODES(N1)%X, MESH%NODES(N1)%Y, MESH%NODES(N1)%Z]
        V2 = [MESH%NODES(N2)%X, MESH%NODES(N2)%Y, MESH%NODES(N2)%Z]
        V3 = [MESH%NODES(N3)%X, MESH%NODES(N3)%Y, MESH%NODES(N3)%Z]
        R_OPP(1) = MESH%NODES(RWG%NEG_OPP_VERTEX)%X
        R_OPP(2) = MESH%NODES(RWG%NEG_OPP_VERTEX)%Y
        R_OPP(3) = MESH%NODES(RWG%NEG_OPP_VERTEX)%Z
        COEF = -RWG%NEG_COEF  ! 负三角形系数带负号

        CALL SUBDIVIDE_TRI_4(V1, V2, V3, SUB_V)
        DO ISUB = 1, 4
            SUB_AREA = AREA / 4.0
            CALL GET_SUB_TRI_GAUSS_POINTS(SUB_V(:,1,ISUB), SUB_V(:,2,ISUB), SUB_V(:,3,ISUB), &
                                           GDATA, SUB_PTS)
            INTEG = (0.0, 0.0)
            DO P = 1, N_PTS
                R_PT = SUB_PTS(:, P)
                DOT_KR = K * DOT_PRODUCT(R_HAT, R_PT)
                PHASE = CEXP((0.0, 1.0) * DOT_KR)
                INTEG = INTEG + GDATA%WEIGHTS(P) * PHASE * (R_PT - R_OPP)
            END DO
            F_VEC = F_VEC + COEF * SUB_AREA * INTEG
        END DO
    END SUBROUTINE CALC_RWG_RADIATION_INTEGRAL

    ! 计算远区散射场（不包含 exp(-jkr)/r 因子）
    ! E_scat(r̂) = -j·k·η₀/(4π) * exp(-jkr)/r * F_scat(r̂)
    ! 返回 F_scat = (I - r̂r̂) · Σ I_n ∫ f_n exp(jk r̂·r') dS'
    SUBROUTINE CALC_FAR_SCATTERED_FIELD(MESH, RWG_CURRENTS, N_RWG, K, ETA0, &
                                        R_HAT, GDATA, E_THETA, E_PHI)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: N_RWG
        COMPLEX, INTENT(IN) :: RWG_CURRENTS(N_RWG)
        REAL, INTENT(IN) :: K, ETA0
        REAL, INTENT(IN) :: R_HAT(3)
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        COMPLEX, INTENT(OUT) :: E_THETA, E_PHI

        COMPLEX :: F_VEC(3), F_SUM(3)
        COMPLEX :: F_THETA, F_PHI
        REAL :: THETA_HAT(3), PHI_HAT(3), THETA_VAL, PHI_VAL
        INTEGER :: N

        ! 球形单位矢量
        ! r_hat = (sinθ cosφ, sinθ sinφ, cosθ)
        THETA_VAL = ACOS(R_HAT(3))  ! θ = acos(z)
        PHI_VAL = ATAN2(R_HAT(2), R_HAT(1))

        THETA_HAT(1) = COS(THETA_VAL) * COS(PHI_VAL)
        THETA_HAT(2) = COS(THETA_VAL) * SIN(PHI_VAL)
        THETA_HAT(3) = -SIN(THETA_VAL)

        PHI_HAT(1) = -SIN(PHI_VAL)
        PHI_HAT(2) = COS(PHI_VAL)
        PHI_HAT(3) = 0.0

        ! 对所有 RWG 基函数做辐射积分
        F_SUM = (0.0, 0.0)
        DO N = 1, N_RWG
            CALL CALC_RWG_RADIATION_INTEGRAL(MESH, MESH%RWG_BASES(N), K, R_HAT, GDATA, F_VEC)
            F_SUM = F_SUM + RWG_CURRENTS(N) * F_VEC
        END DO

        ! 投影到 θ̂ 和 φ̂
        F_THETA = DOT_PRODUCT(THETA_HAT, F_SUM)
        F_PHI = DOT_PRODUCT(PHI_HAT, F_SUM)

        ! E_s = -j·k·η₀/(4π) * exp(-jkr)/r * F （这里返回不含 exp(-jkr)/r 因子的部分）
        ! 实际返回的是电场的"方向图"部分
        E_THETA = (0.0, -1.0) * K * ETA0 / (4.0 * PI) * F_THETA
        E_PHI   = (0.0, -1.0) * K * ETA0 / (4.0 * PI) * F_PHI
    END SUBROUTINE CALC_FAR_SCATTERED_FIELD

    ! 计算双站 RCS：固定入射方向，扫描观察角
    ! 输出：OBS_ANGLES(N_ANGLES) 度, RCS_THETA(N_ANGLES) dBsm, RCS_PHI(N_ANGLES) dBsm
    SUBROUTINE CALC_BISTATIC_RCS(MESH, RWG_CURRENTS, N_RWG, K, ETA0, E0, &
                                  N_ANGLES, OBS_ANGLES_DEG, RCS_THETA_DBSM, RCS_PHI_DBSM, GDATA)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: N_RWG
        COMPLEX, INTENT(IN) :: RWG_CURRENTS(N_RWG)
        REAL, INTENT(IN) :: K, ETA0, E0
        INTEGER, INTENT(IN) :: N_ANGLES
        REAL, INTENT(OUT) :: OBS_ANGLES_DEG(N_ANGLES)
        REAL, INTENT(OUT) :: RCS_THETA_DBSM(N_ANGLES)
        REAL, INTENT(OUT) :: RCS_PHI_DBSM(N_ANGLES)
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA

        INTEGER :: I
        REAL :: THETA_DEG, THETA_RAD, R_HAT(3)
        COMPLEX :: E_THETA, E_PHI
        REAL :: SIG_THETA, SIG_PHI

        DO I = 1, N_ANGLES
            THETA_DEG = 0.0 + REAL(I - 1) * 180.0 / REAL(N_ANGLES - 1)
            THETA_RAD = THETA_DEG * PI / 180.0
            OBS_ANGLES_DEG(I) = THETA_DEG

            ! 在 φ=0 平面内扫描
            R_HAT(1) = SIN(THETA_RAD)
            R_HAT(2) = 0.0
            R_HAT(3) = COS(THETA_RAD)

            CALL CALC_FAR_SCATTERED_FIELD(MESH, RWG_CURRENTS, N_RWG, K, ETA0, &
                                          R_HAT, GDATA, E_THETA, E_PHI)

            ! RCS 定义：σ = 4π r² |E_s|² / |E_inc|²
            ! E_s = E_THETA * exp(-jkr)/r  （或 E_PHI * ...）
            ! 所以：σ_θ = 4π |E_THETA|² / E₀², σ_φ = 4π |E_PHI|² / E₀²
            SIG_THETA = 4.0 * PI * (ABS(E_THETA) ** 2) / (E0 ** 2)
            SIG_PHI   = 4.0 * PI * (ABS(E_PHI) ** 2) / (E0 ** 2)

            IF (SIG_THETA > 1.0E-20) THEN
                RCS_THETA_DBSM(I) = 10.0 * LOG10(SIG_THETA)
            ELSE
                RCS_THETA_DBSM(I) = -200.0
            END IF
            IF (SIG_PHI > 1.0E-20) THEN
                RCS_PHI_DBSM(I) = 10.0 * LOG10(SIG_PHI)
            ELSE
                RCS_PHI_DBSM(I) = -200.0
            END IF
        END DO
    END SUBROUTINE CALC_BISTATIC_RCS

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

    ! 由面积坐标得到子三角形上的全局高斯点
    SUBROUTINE GET_SUB_TRI_GAUSS_POINTS(S1, S2, S3, GDATA, SUB_PTS)
        REAL, INTENT(IN)  :: S1(3), S2(3), S3(3)
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        REAL, INTENT(OUT) :: SUB_PTS(3, GDATA%N_POINTS)
        INTEGER :: P, K
        REAL :: U, V, W
        DO P = 1, GDATA%N_POINTS
            U = GDATA%UVW(1, P)
            V = GDATA%UVW(2, P)
            W = GDATA%UVW(3, P)
            DO K = 1, 3
                SUB_PTS(K, P) = U * S1(K) + V * S2(K) + W * S3(K)
            END DO
        END DO
    END SUBROUTINE GET_SUB_TRI_GAUSS_POINTS

END MODULE FAR_FIELD
