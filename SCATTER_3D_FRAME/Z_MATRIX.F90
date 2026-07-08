MODULE Z_MATRIX
    USE EM_TYPES
    USE NUMERICAL_INTEGRATION
    USE GREEN_FUNCTIONS
    USE RWG_BASIS_BUILD
    USE SINGULAR_INTEGRAL
    IMPLICIT NONE
CONTAINS

    ! 计算两个 RWG 基函数之间的 EFIE 非对角阻抗矩阵元素。
    SUBROUTINE CALC_EFIE_MATRIX_ELEMENT(MESH, RWG_M, RWG_N, GDATA, K, ETA0, Z_MN)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        TYPE(RWG_BASIS), INTENT(IN) :: RWG_M, RWG_N ! M代表源，N代表场
        REAL, INTENT(IN) :: ETA0 ! 自由空间波阻抗
        REAL, INTENT(IN) :: K
        COMPLEX, INTENT(OUT) :: Z_MN
        INTEGER :: TRI_M(2), TRI_N(2) ! RWG基函数对应的正负三角形的ID
        REAL :: COEF_M(2), COEF_N(2) ! 正负三角形的系数
        REAL :: R_OPP_M(3, 2), R_OPP_N(3, 2) ! 正负三角形中的两个对顶点的坐标
        INTEGER :: I, J
        COMPLEX :: P_TERM, Q_TERM
        COMPLEX :: JW_MU, J_OVER_WE
        COMPLEX :: I1 ! 标量
        COMPLEX :: I4 ! 标量
        COMPLEX :: I2(3) ! 矢量
        COMPLEX :: I3(3) ! 矢量

        JW_MU = (0.0, 1.0) * ETA0 * K
        J_OVER_WE = (0.0, 1.0) * ETA0 / K
        Z_MN = (0.0, 0.0)

        TRI_M(1) = RWG_M%POS_TRI_ID
        TRI_M(2) = RWG_M%NEG_TRI_ID
        TRI_N(1) = RWG_N%POS_TRI_ID
        TRI_N(2) = RWG_N%NEG_TRI_ID

        COEF_M(1) = RWG_M%POS_COEF
        COEF_M(2) = -RWG_M%NEG_COEF
        COEF_N(1) = RWG_N%POS_COEF
        COEF_N(2) = -RWG_N%NEG_COEF

        ! 给对顶点的坐标赋值。对应ri，rj。
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

        ! 双重循环。Tm+，Tm-，Tn+，Tn-
        DO I = 1, 2 ! 源。
            DO J = 1, 2
                CALL CALC_GREEN_INTEGALS(MESH, TRI_M(I), TRI_N(J), GDATA, K, I1, I2, I3, I4)
                P_TERM = DOT_PRODUCT(R_OPP_M(:, I), R_OPP_N(:, J)) * I1 - &
                            DOT_PRODUCT(R_OPP_N(:, J), I2) - DOT_PRODUCT(R_OPP_M(:, I), I3) + I4
                Q_TERM = -4.0 * I1

                Z_MN = Z_MN + COEF_M(I) * COEF_N(J) * (JW_MU * P_TERM + J_OVER_WE * Q_TERM)
            END DO
        END DO

    END SUBROUTINE

    ! 计算两个RWG三角形之间的MFIE非对角阻抗矩阵元素
    SUBROUTINE CALC_MFIE_MATRIX_ELEMENT(MESH, RWG_M, RWG_N, GDATA, K, Z_MN)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        TYPE(RWG_BASIS), INTENT(IN) :: RWG_M, RWG_N ! m是场，n是源
        REAL, INTENT(IN) :: K
        COMPLEX, INTENT(OUT) :: Z_MN
        INTEGER :: TRI_M(2), TRI_N(2) ! RWG基函数对应的正负三角形的ID
        REAL :: COEF_M(2), COEF_N(2) ! 正负三角形的系数
        REAL :: R_OPP_M(3, 2), R_OPP_N(3, 2) ! 正负三角形中的两个对顶点的坐标
        INTEGER :: I, J, M, N

        INTEGER :: TRI_M_ID, TRI_N_ID
        REAL :: C_M, C_N
        REAL :: R_M(3), R_N(3)
        REAL :: N_M(3) ! 场三角形的单位法向量
        REAL :: A_M, A_N ! 场、源三角形的面积
        REAL :: GLOBAL_PTS_M(3, GDATA%N_POINTS), GLOBAL_PTS_N(3, GDATA%N_POINTS) ! 三角形高斯积分点的全局坐标。ri，rj
        REAL :: R_I(3), R_J(3)
        COMPLEX :: Z_PAIR  
        COMPLEX :: K_ACC(3) ! 复数向量，累加源积分
        REAL :: F_M(3), F_N(3), R_VEC(3), DIST
        COMPLEX :: GRAD_G(3), CROSS(3), NCK(3), DOT
        REAL :: F_VAL_M(3), F_VAL_N(3), DOT_FM_FN, HALF_TERM
        CHARACTER * 3 :: TRI_M_LOCAL, TRI_N_LOCAL

        TRI_M(1) = RWG_M%POS_TRI_ID
        TRI_M(2) = RWG_M%NEG_TRI_ID
        TRI_N(1) = RWG_N%POS_TRI_ID
        TRI_N(2) = RWG_N%NEG_TRI_ID

        COEF_M(1) = RWG_M%POS_COEF
        COEF_M(2) = -RWG_M%NEG_COEF
        COEF_N(1) = RWG_N%POS_COEF
        COEF_N(2) = -RWG_N%NEG_COEF
        
        ! 给对顶点的坐标赋值。对应rm，rn。
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

        DO M = 1, 2 ! 场
            DO N = 1, 2
                TRI_M_ID = TRI_M(M)
                TRI_N_ID = TRI_N(N)
                C_M = COEF_M(M)
                C_N = COEF_N(N)
                R_M = R_OPP_M(:, M)
                R_N = R_OPP_N(:, N)
                N_M = MESH%TRIANGLES(TRI_M_ID)%NORMAL
                A_M = MESH%TRIANGLES(TRI_M_ID)%AREA
                A_N = MESH%TRIANGLES(TRI_N_ID)%AREA

                ! 获取高斯点
                CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_M_ID, GDATA, GLOBAL_PTS_M)
                CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_N_ID, GDATA, GLOBAL_PTS_N)

                ! ------------同一个三角形，只算1/2项，当两个RWG基函数有共享的三角形时，此项不为0--------
                IF (TRI_M_ID == TRI_N_ID) THEN

                    HALF_TERM = 0.0

                    IF(RWG_M%POS_TRI_ID == TRI_M_ID) THEN
                        TRI_M_LOCAL = 'POS'
                    ELSE
                        TRI_M_LOCAL = 'NEG'
                    END IF

                    IF(RWG_N%POS_TRI_ID == TRI_N_ID) THEN
                        TRI_N_LOCAL = "POS"
                    ELSE
                        TRI_N_LOCAL = 'NEG'
                    END IF
                    
                    DO I = 1, GDATA%N_POINTS

                        CALL EVAL_RWG_BASIS(MESH, RWG_M, TRI_M_LOCAL, GLOBAL_PTS_M(:, I), F_VAL_M) ! 得到fm
                        CALL EVAL_RWG_BASIS(MESH, RWG_N, TRI_N_LOCAL, GLOBAL_PTS_N(:, I), F_VAL_N) ! 得到fn
                        DOT_FM_FN = DOT_PRODUCT(F_VAL_M, F_VAL_N)
                        HALF_TERM = HALF_TERM + GDATA%WEIGHTS(I) * DOT_FM_FN

                    END DO
                        
                    HALF_TERM = 0.5 * HALF_TERM * MESH%TRIANGLES(TRI_M_ID)%AREA
                    Z_MN = CMPLX(HALF_TERM) + Z_MN
                    
                    ! PRINT *, "HALF_TERM for RWG pair:", HALF_TERM

                    CYCLE
                END IF
 
                ! -------------------不同三角形，先判断是否共面----------------------------------------
                IF (IS_COPLANAR_TRI_PAIR(MESH, TRI_M_ID, TRI_N_ID, 1.0E-6, 1.0E-6)) CYCLE

                ! ---------------------------不共面，计算K算子------------------------------------------
                Z_PAIR = 0.0
                DO I = 1, GDATA%N_POINTS
                    
                    R_I = GLOBAL_PTS_M(:, I)
                    F_M = C_M * (R_I - R_M)
                    K_ACC = (0.0, 0.0)

                    DO J = 1, GDATA%N_POINTS
                        R_J = GLOBAL_PTS_N(:, J)
                        F_N = C_N * (R_J - R_N)
                        R_VEC = R_I - R_J
                        DIST = SQRT(SUM(R_VEC ** 2))
                        IF (DIST > 1.0E-10) THEN
                            CALL GARD_GREEN_FUNC(R_I, R_J, K, GRAD_G)
                            
                            ! 计算叉乘
                            CROSS(1) = F_N(2) * GRAD_G(3) - F_N(3) * GRAD_G(2)
                            CROSS(2) = F_N(3) * GRAD_G(1) - F_N(1) * GRAD_G(3)
                            CROSS(3) = F_N(1) * GRAD_G(2) - F_N(2) * GRAD_G(1)

                            K_ACC = K_ACC + GDATA%WEIGHTS(J) * CROSS
                        END IF
                    END DO
                    
                    K_ACC = K_ACC * A_N

                    ! 计算 N_CROSS_K = N_M × K_ACC（n_m x K_ACC）
                    NCK(1) = N_M(2) * K_ACC(3) - N_M(3) * K_ACC(2)
                    NCK(2) = N_M(3) * K_ACC(1) - N_M(1) * K_ACC(3)
                    NCK(3) = N_M(1) * K_ACC(2) - N_M(2) * K_ACC(1)

                    ! 计算点乘DOT = F_M · N_CROSS_K
                    DOT = F_M(1) * NCK(1) + F_M(2) * NCK(2) + F_M(3) * NCK(3)

                    Z_PAIR = Z_PAIR + GDATA%WEIGHTS(I) * DOT
                END DO
                
                Z_PAIR = Z_PAIR * A_M
                Z_MN = Z_MN + Z_PAIR
            
            END DO
        END DO
    END SUBROUTINE

    ! 判断两三角形是否共面，用于K算子相关计算。共面则K算子相关项为0
    FUNCTION IS_COPLANAR_TRI_PAIR(MESH, TRI_M_ID, TRI_N_ID, TOL_NORM, TOL_DIST) RESULT(IS_COPLANAR)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_M_ID, TRI_N_ID
        REAL, INTENT(IN) :: TOL_NORM, TOL_DIST     ! 法向平行判定阈值, 重心连线垂直判定阈值
        REAL :: DOT_N, DOT_G, R_VEC(3)
        LOGICAL :: IS_COPLANAR

        IS_COPLANAR = .TRUE.

        ! 法向平行判定
        DOT_N = DOT_PRODUCT(MESH%TRIANGLES(TRI_M_ID)%NORMAL, MESH%TRIANGLES(TRI_N_ID)%NORMAL)
        IF (ABS(ABS(DOT_N) - 1.0) > TOL_NORM) IS_COPLANAR = .FALSE.

        ! 重心连线垂直判定，验证两三角形重心连线与法向垂直
        R_VEC = MESH%TRIANGLES(TRI_M_ID)%CENTROID - MESH%TRIANGLES(TRI_N_ID)%CENTROID
        DOT_G = DOT_PRODUCT(MESH%TRIANGLES(TRI_M_ID)%NORMAL, R_VEC)
        IF (ABS(ABS(DOT_G)) > TOL_DIST) IS_COPLANAR = .FALSE.

    END FUNCTION

    ! 计算单个三角形T对自身的完整 EFIE 阻抗贡献。
    SUBROUTINE CALC_EFIE_SELF_TRI_PAIR(MESH, TRI_ID, RWG, GDATA, K, ETA0, Z_PAIR)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        INTEGER, INTENT(IN) :: TRI_ID  
        TYPE(RWG_BASIS), INTENT(IN) :: RWG
        REAL, INTENT(IN) :: ETA0 ! 自由空间波阻抗
        REAL, INTENT(IN) :: K    ! 波数
        COMPLEX, INTENT(OUT) :: Z_PAIR

        REAL :: C   ! RWG基函数系数
        REAL :: R_OPP(3) ! 对顶点坐标
        REAL :: AREA ! 三角形面积
        COMPLEX :: JW_MU, J_OVER_WE
        COMPLEX :: P_SING     ! P算子1/R解析部分
        COMPLEX :: Q_SING     ! Q算子1/R解析部分
        COMPLEX :: P_SMOOTH2  ! P算子G_smooth2数值部分
        COMPLEX :: Q_SMOOTH2  ! Q算子G_smooth2数值部分
        REAL :: I0_VAL   ! 当前场点的标量势
        REAL :: I1_VAL, I2_VAL, I3_VAL  ! 当前场点的线性势
        REAL :: IVEC_VAL(3) ! 当前场点的矢量势
        REAL :: SUM_I0   ! 外层高斯累计∑wi * I0_VAL
        REAL :: SUM_R_I0(3)   ! 外层高斯累计∑wi * ri * I0_VAL
        REAL :: SUM_IVEC(3)   ! 外层高斯累计∑wi * IVEC_VAL
        REAL :: SUM_R_DOT_IVEC   ! 外层高斯累计∑wi * ri * (IVEC_VAL · ri)
        REAL :: J0  ! SUM_I0 * AREA
        REAL :: J_R_I0(3)   ! SUM_R_I0 * AREA
        REAL :: J_VEC(3)   ! SUM_IVEC * AREA
        REAL :: J_R_DOT   ! SUM_R_DOT_IVEC * AREA
        INTEGER :: I
        REAL :: R_VEC(3) ! 当前高斯点
        REAL :: V1(3), V2(3), V3(3) ! 三角形TRI_ID的顶点坐标
        REAL :: RESULT_3(3) ! 线性势结果
        REAL :: TERM1, TERM2, TERM3, TERM4
        COMPLEX :: I1 ! 标量
        COMPLEX :: I4 ! 标量
        COMPLEX :: I2(3) ! 矢量
        COMPLEX :: I3(3) ! 矢量

        ! 提取RWG参数
        IF (TRI_ID == RWG%POS_TRI_ID) THEN
            C = RWG%POS_COEF
            R_OPP(1) = MESH%NODES(RWG%POS_OPP_VERTEX)%X
            R_OPP(2) = MESH%NODES(RWG%POS_OPP_VERTEX)%Y
            R_OPP(3) = MESH%NODES(RWG%POS_OPP_VERTEX)%Z
        ELSE IF (TRI_ID == RWG%NEG_TRI_ID) THEN
            C = -RWG%NEG_COEF
            R_OPP(1) = MESH%NODES(RWG%NEG_OPP_VERTEX)%X
            R_OPP(2) = MESH%NODES(RWG%NEG_OPP_VERTEX)%Y
            R_OPP(3) = MESH%NODES(RWG%NEG_OPP_VERTEX)%Z
        END IF

        ! 预计算常数
        JW_MU = (0.0, 1.0) * ETA0 * K
        J_OVER_WE = (0.0, 1.0) * ETA0 / K

        ! 外层高斯积分（解析部分）
        SUM_I0 = 0.0
        SUM_R_I0 = 0.0
        SUM_IVEC = 0.0
        SUM_R_DOT_IVEC = 0.0
        CALL INIT_GAUSS_TRI_DATA(GAUSS_7PT, GDATA) ! 初始化高斯积分数据
        DO I = 1, GDATA%N_POINTS
            CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_ID, GDATA, R_VEC) ! 获取当前高斯点的全局坐标
            ! 三角形TRI_ID的三个顶点
            V1(1) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(1))%X
            V1(2) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(1))%Y
            V1(3) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(1))%Z
            V2(1) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(2))%X
            V2(2) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(2))%Y
            V2(3) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(2))%Z
            V3(1) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(3))%X
            V3(2) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(3))%Y
            V3(3) = MESH%NODES(MESH%TRIANGLES(TRI_ID)%VERTEX_3D(3))%Z

            CALL ANALYTIC_SCALAR_POT_1_OVER_R(R_VEC, V1, V2, V3, I0_VAL) ! 计算标量势
            CALL ANALYTIC_LINEAR_POT_1_OVER_R(R_VEC, V1, V2, V3, RESULT_3) ! 计算线性势
            I1_VAL = RESULT_3(1)
            I2_VAL = RESULT_3(2)
            I3_VAL = RESULT_3(3)
            
            ! 计算矢量势 IVEC_VAL = I1_VAL * V1 + I2_VAL * V2 + I3_VAL * V3
            IVEC_VAL(1) = V1(1) * I1_VAL + V2(1) * I2_VAL + V3(1) * I3_VAL
            IVEC_VAL(2) = V1(2) * I1_VAL + V2(2) * I2_VAL + V3(2) * I3_VAL
            IVEC_VAL(3) = V1(3) * I1_VAL + V2(3) * I2_VAL + V3(3) * I3_VAL

            ! 累加外层高斯积分
            SUM_I0 = SUM_I0 + GDATA%WEIGHTS(I) * I0_VAL
            SUM_R_I0 = SUM_R_I0 + GDATA%WEIGHTS(I) * R_VEC * I0_VAL
            SUM_IVEC = SUM_IVEC + GDATA%WEIGHTS(I) * IVEC_VAL
            SUM_R_DOT_IVEC = SUM_R_DOT_IVEC + GDATA%WEIGHTS(I) * DOT_PRODUCT(R_VEC, IVEC_VAL)

            J0 = SUM_I0 * MESH%TRIANGLES(TRI_ID)%AREA
            J_R_I0 = SUM_R_I0 * MESH%TRIANGLES(TRI_ID)%AREA
            J_VEC = SUM_IVEC * MESH%TRIANGLES(TRI_ID)%AREA
            J_R_DOT = SUM_R_DOT_IVEC * MESH%TRIANGLES(TRI_ID)%AREA   
            
        END DO

        ! 计算解析部分的Q算子
        Q_SING = -C**2 * J0 / PI
        ! 计算解析部分的P算子
        TERM1 = J_R_DOT
        TERM2 = -DOT_PRODUCT(R_OPP, J_R_I0)
        TERM3 = -DOT_PRODUCT(R_OPP, J_VEC)
        TERM4 = DOT_PRODUCT(R_OPP, R_OPP) * J0
        P_SING = C**2 * (TERM1 + TERM2 + TERM3 + TERM4) / (4.0 * PI)

        ! G_SMOOTH数值部分，使用G_SMOOTH
        CALL CALC_GREEN_SMOOTH_INTEGALS(MESH, TRI_ID, TRI_ID, GDATA, K, I1, I2, I3, I4)
        P_SMOOTH2 = C**2 * (DOT_PRODUCT(R_OPP, R_OPP) * I1 - DOT_PRODUCT(R_OPP, I2) - DOT_PRODUCT(R_OPP, I3) + I4)
        Q_SMOOTH2 = -4.0 * C**2 * I1

        Z_PAIR = JW_MU * (P_SING + P_SMOOTH2) + J_OVER_WE * (Q_SING + Q_SMOOTH2)

    END SUBROUTINE

END MODULE