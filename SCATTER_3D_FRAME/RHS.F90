MODULE RHS
    USE EM_TYPES
    USE NUMERICAL_INTEGRATION
    IMPLICIT NONE
CONTAINS

    ! 平面波场计算子程序，给定空间一点r，计算该点的Ei和Hi
    SUBROUTINE CALC_PLANE_WAVE(R_PT, K, K_HAT, E_POL, E0, ETA0, E_INC, H_INC)
        REAL, INTENT(IN) :: R_PT(3) ! 场点坐标r
        REAL, INTENT(IN) :: K
        REAL, INTENT(IN) :: K_HAT(3)
        REAL, INTENT(IN) :: E_POL(3)
        REAL, INTENT(IN) :: E0
        REAL, INTENT(IN) :: ETA0
        COMPLEX, INTENT(OUT) :: E_INC(3)
        COMPLEX, INTENT(OUT) :: H_INC(3)
        COMPLEX :: PHASE
        REAL :: CROSS(3)
        INTEGER :: I

        PHASE = CEXP((0.0, -1.0) * K * DOT_PRODUCT(K_HAT, R_PT))
        E_INC = E0 * E_POL * PHASE

        DO I = 1, 3
            CROSS(1) = K_HAT(2) * E_POL(3) - E_POL(2) * K_HAT(3)
            CROSS(2) = E_POL(1) * K_HAT(3) - K_HAT(1) * E_POL(3)
            CROSS(3) = K_HAT(1) * E_POL(2) - E_POL(1) * K_HAT(2)
        END DO
        
        H_INC = (E0 / ETA0) * CROSS * PHASE

    END SUBROUTINE

    ! 对单个RWG基函数fm，计算V_m(E)
    SUBROUTINE CALC_EFIE_RHS(MESH, RWG, GDATA, K, K_HAT, E_POL, E0, ETA0, V_M_E)
        TYPE(MESH_3D), intent(IN) :: MESH
        TYPE(RWG_BASIS), INTENT(IN) :: RWG
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        REAL, INTENT(IN) :: K
        REAL, INTENT(IN) :: ETA0
        REAL, INTENT(IN) :: K_HAT(3)
        REAL, INTENT(IN) :: E_POL(3)
        REAL, INTENT(IN) :: E0
        COMPLEX, INTENT(OUT) :: V_M_E

        INTEGER :: TRI_M(2)
        REAL :: COEF_M(2), R_OPP_M(3, 2)
        INTEGER :: M, I
        REAL :: GLOBAL_PTS_M(3, GDATA%N_POINTS)
        COMPLEX :: E_INC(3), H_INC(3)
        REAL :: R_PT(3), A_M, F_M(3)
        COMPLEX :: DOT
        COMPLEX :: V_LOCAL

        TRI_M(1) = RWG%POS_TRI_ID
        TRI_M(2) = RWG%NEG_TRI_ID
        COEF_M(1) = RWG%POS_COEF
        COEF_M(2) = -RWG%NEG_COEF
        R_OPP_M(1, 1) = MESH%NODES(RWG%POS_OPP_VERTEX)%X
        R_OPP_M(2, 1) = MESH%NODES(RWG%POS_OPP_VERTEX)%Y
        R_OPP_M(3, 1) = MESH%NODES(RWG%POS_OPP_VERTEX)%Z
        R_OPP_M(1, 2) = MESH%NODES(RWG%NEG_OPP_VERTEX)%X
        R_OPP_M(2, 2) = MESH%NODES(RWG%NEG_OPP_VERTEX)%Y
        R_OPP_M(3, 2) = MESH%NODES(RWG%NEG_OPP_VERTEX)%Z

        V_M_E = (0.0, 0.0)

        DO M = 1, 2
            A_M = MESH%TRIANGLES(TRI_M(M))%AREA
            CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_M(M), GDATA, GLOBAL_PTS_M)
            V_LOCAL = (0.0, 0.0) ! 局部累加，隔离当前三角形
            DO I = 1, GDATA%N_POINTS
                R_PT = GLOBAL_PTS_M(:, I)
                CALL CALC_PLANE_WAVE(R_PT, K, K_HAT, E_POL, E0, ETA0, E_INC, H_INC)
                F_M = COEF_M(M) * (R_PT - R_OPP_M(:, M))
                DOT = DOT_PRODUCT(F_M, E_INC)
                V_LOCAL = V_LOCAL + GDATA%WEIGHTS(I) * DOT
            END DO
            V_M_E = V_LOCAL * A_M + V_M_E
        END DO
    END SUBROUTINE

    ! 计算右端项V_m(M)
    SUBROUTINE CALC_MFIE_RHS(MESH, RWG, GDATA, K, K_HAT, E_POL, E0, ETA0, V_M_M)
        TYPE(MESH_3D), intent(IN) :: MESH
        TYPE(RWG_BASIS), INTENT(IN) :: RWG
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        REAL, INTENT(IN) :: K
        REAL, INTENT(IN) :: K_HAT(3)
        REAL, INTENT(IN) :: E_POL(3)
        REAL, INTENT(IN) :: E0
        REAL, INTENT(IN) :: ETA0
        COMPLEX, INTENT(OUT) :: V_M_M
        
        INTEGER :: TRI_M(2)
        REAL :: COEF_M(2), R_OPP_M(3, 2)
        INTEGER :: M, I
        REAL :: GLOBAL_PTS_M(3, GDATA%N_POINTS)
        COMPLEX :: E_INC(3), H_INC(3), NCH(3)
        REAL :: R_PT(3), A_M, F_M(3), N_M(3)
        COMPLEX :: DOT
        COMPLEX :: V_LOCAL

        TRI_M(1) = RWG%POS_TRI_ID
        TRI_M(2) = RWG%NEG_TRI_ID
        COEF_M(1) = RWG%POS_COEF
        COEF_M(2) = -RWG%NEG_COEF
        R_OPP_M(1, 1) = MESH%NODES(RWG%POS_OPP_VERTEX)%X
        R_OPP_M(2, 1) = MESH%NODES(RWG%POS_OPP_VERTEX)%Y
        R_OPP_M(3, 1) = MESH%NODES(RWG%POS_OPP_VERTEX)%Z
        R_OPP_M(1, 2) = MESH%NODES(RWG%NEG_OPP_VERTEX)%X
        R_OPP_M(2, 2) = MESH%NODES(RWG%NEG_OPP_VERTEX)%Y
        R_OPP_M(3, 2) = MESH%NODES(RWG%NEG_OPP_VERTEX)%Z

        V_M_M = (0.0, 0.0)

        DO M = 1, 2
            A_M = MESH%TRIANGLES(TRI_M(M))%AREA
            N_M = MESH%TRIANGLES(TRI_M(M))%NORMAL
            CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_M(M), GDATA, GLOBAL_PTS_M)
            V_LOCAL = (0.0, 0.0)
            DO I = 1, GDATA%N_POINTS
                R_PT = GLOBAL_PTS_M(:, I)
                CALL CALC_PLANE_WAVE(R_PT, K, K_HAT, E_POL, E0, ETA0, E_INC, H_INC)
                NCH(1) = N_M(2) * H_INC(3) - N_M(3) * H_INC(2)
                NCH(2) = N_M(3) * H_INC(1) - N_M(1) * H_INC(3)
                NCH(3) = N_M(1) * H_INC(2) - N_M(2) * H_INC(1)
                F_M = COEF_M(M) * (R_PT - R_OPP_M(:, M))
                DOT = DOT_PRODUCT(F_M, NCH)
                V_LOCAL = V_LOCAL + GDATA%WEIGHTS(I) * DOT
            END DO
            V_M_M = V_LOCAL * A_M + V_M_M
        END DO        
    END SUBROUTINE

END MODULE