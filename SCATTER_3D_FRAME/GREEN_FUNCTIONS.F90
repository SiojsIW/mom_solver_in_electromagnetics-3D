MODULE GREEN_FUNCTIONS
    USE EM_TYPES
    USE NUMERICAL_INTEGRATION
    IMPLICIT NONE
CONTAINS

    ! 计算格林函数，输入两个坐标，计算两个坐标点的格林函数
    COMPLEX FUNCTION GREEN_FUNC(R, RP, K) RESULT(G)
        REAL, INTENT(IN):: R(3) ! 源点坐标
        REAL, INTENT(IN):: RP(3) ! 场点坐标
        REAL,INTENT(IN) :: K  ! 波数
        REAL :: DIST
        
        ! COMPLEX 类型，虚数单位是 (0.0, 1.0)。CEXP 是复数指数。PI = 4.0*ATAN(1.0)。
        DIST = SQRT(SUM((R - RP) ** 2)) ! 场点与源点距离

        IF (DIST < 1.0E-10) THEN
            G = (0.0, 0.0)
        ELSE
            G = CEXP((0.0, -1.0) * K * DIST) / (4.0 * PI * DIST)
        END IF
    END FUNCTION

    ! 计算阻抗元素Z_MN，只是计算格林函数在两个三角形上的，没有f。
    SUBROUTINE CALC_GREEN_MATRIX_ELEMENT(MESH, TM_ID, TN_ID, GDATA, K, Z_MN)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        INTEGER :: TM_ID ! 场三角形编号
        INTEGER :: TN_ID ! 源三角形编号
        INTEGER :: I, J
        REAL :: K  ! 波数
        REAL, ALLOCATABLE :: GLOBAL_PTS_M(:, :) ! 场三角形高斯积分点的全局坐标
        REAL, ALLOCATABLE :: GLOBAL_PTS_N(:, :) ! 源三角形高斯积分点的全局坐标
        COMPLEX :: G
        COMPLEX, INTENT(OUT) :: Z_MN 

        ALLOCATE(GLOBAL_PTS_M(3, GDATA%N_POINTS)) 
        ALLOCATE(GLOBAL_PTS_N(3, GDATA%N_POINTS)) 
        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TM_ID, GDATA, GLOBAL_PTS_M)  ! 得到场三角形的高斯积分点的全局坐标
        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TN_ID, GDATA, GLOBAL_PTS_N)  ! 得到源三角形的高斯积分点的全局坐标

        Z_MN = (0.0, 0.0)
        DO I = 1, GDATA%N_POINTS ! 场
            DO J = 1, GDATA%N_POINTS ! 源
                G = GREEN_FUNC(GLOBAL_PTS_M(:, I), GLOBAL_PTS_N(:, J), K) 
                Z_MN = Z_MN + GDATA%WEIGHTS(I) * GDATA%WEIGHTS(J) * G
            END DO
        END DO

        DEALLOCATE(GLOBAL_PTS_M) 
        DEALLOCATE(GLOBAL_PTS_N) 

        Z_MN = Z_MN * MESH%TRIANGLES(TM_ID)%AREA * MESH%TRIANGLES(TN_ID)%AREA

    END SUBROUTINE

    ! 计算一个三角形对上的四个积分I1,I2,I3,I4，（计算EFIE非对角矩阵元素用）
    SUBROUTINE CALC_GREEN_INTEGALS(MESH, TRI_A, TRI_B, GDATA, K, I1, I2, I3, I4)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        INTEGER, INTENT(IN) :: TRI_A ! 源三角形
        INTEGER, INTENT(IN) :: TRI_B ! 场三角形
        REAL, INTENT(IN) :: K  ! 波数
        COMPLEX, INTENT(OUT) :: I1 ! 标量
        COMPLEX, INTENT(OUT) :: I4 ! 标量
        COMPLEX, INTENT(OUT) :: I2(3) ! 矢量
        COMPLEX, INTENT(OUT) :: I3(3) ! 矢量
    
        REAL :: AREA_A, AREA_B
        REAL :: GLOBAL_PTS_A(3, GDATA%N_POINTS), GLOBAL_PTS_B(3, GDATA%N_POINTS) ! 三角形高斯积分点的全局坐标
        INTEGER :: I, J
        COMPLEX :: G
        I1 = (0.0, 0.0)
        I4 = (0.0, 0.0)
        I2 = (0.0, 0.0)   ! 标量复数自动广播到3个元素
        I3 = (0.0, 0.0)

        AREA_A = MESH%TRIANGLES(TRI_A)%AREA
        AREA_B = MESH%TRIANGLES(TRI_B)%AREA

        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_A, GDATA, GLOBAL_PTS_A)
        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_B, GDATA, GLOBAL_PTS_B)

        DO I = 1, GDATA%N_POINTS ! i代表源点
            DO J = 1, GDATA%N_POINTS
                G = GREEN_FUNC(GLOBAL_PTS_A(:, I), GLOBAL_PTS_B(:, J), K)

                I1 = I1 + GDATA%WEIGHTS(I) * GDATA%WEIGHTS(J) * G  ! 先不乘面积，循环结束后统一乘面积
                I2 = I2 + GDATA%WEIGHTS(I) * GDATA%WEIGHTS(J) * GLOBAL_PTS_B(:, J) * G
                I3 = I3 + GDATA%WEIGHTS(I) * GDATA%WEIGHTS(J) * GLOBAL_PTS_A(:, I) * G
                I4 = I4 + GDATA%WEIGHTS(I) * GDATA%WEIGHTS(J) * &
                        DOT_PRODUCT(GLOBAL_PTS_A(:, I), GLOBAL_PTS_B(:, J)) * G
                
            END DO
        END DO

        I1 = I1 * AREA_A * AREA_B
        I2 = I2 * AREA_A * AREA_B
        I3 = I3 * AREA_A * AREA_B
        I4 = I4 * AREA_A * AREA_B

    END SUBROUTINE

    ! 计算标量格林函数的梯度
    SUBROUTINE GARD_GREEN_FUNC(R, RP, K, GRAD_G)
        REAL, INTENT(IN) :: R(3) ! 场点坐标
        REAL, INTENT(IN) :: RP(3) ! 源点坐标
        REAL, INTENT(IN) :: K  ! 波数
        COMPLEX, INTENT(OUT) :: GRAD_G(3)
        COMPLEX :: FACTOR
        COMPLEX :: G
        REAL :: R_VEC(3), DIST

        R_VEC = R - RP
        DIST = SQRT(SUM(R_VEC ** 2))
        
        IF (DIST >= 1.0E-10) THEN
            G = CEXP((0.0, -1.0) * K * DIST) / (4.0 * PI * DIST)
            FACTOR = -(1.0 / DIST + (0.0, 1.0) * K) * G / DIST
            GRAD_G = FACTOR * R_VEC
        ELSE 
            GRAD_G = (0.0, 0.0)
        END IF

    END SUBROUTINE

END MODULE