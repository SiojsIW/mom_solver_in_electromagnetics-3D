MODULE NUMERICAL_INTEGRATION
    USE EM_TYPES
    IMPLICIT NONE
CONTAINS

    ! 初始化高斯数据（几点求积，面积坐标，权重）
    SUBROUTINE INIT_GAUSS_TRI(N_POINTS, GDATA)
        TYPE(GAUSS_TRI_DATA), INTENT(INOUT) :: GDATA
        INTEGER, INTENT(IN) :: N_POINTS

        IF(ALLOCATED(GDATA%UVW)) DEALLOCATE(GDATA%UVW)
        IF(ALLOCATED(GDATA%WEIGHTS)) DEALLOCATE(GDATA%WEIGHTS)

        IF(N_POINTS == GAUSS_1PT) THEN
            ALLOCATE(GDATA%UVW(3, 1))
            ALLOCATE(GDATA%WEIGHTS(1))

            GDATA%N_POINTS = GAUSS_1PT
            GDATA%UVW(1, 1) = 1.0/3.0; GDATA%UVW(2, GAUSS_1PT) = 1.0/3.0; GDATA%UVW(3, GAUSS_1PT) = 1.0/3.0
            GDATA%WEIGHTS(1) = 1.0 
        END IF

        IF(N_POINTS == GAUSS_3PT) THEN
            ALLOCATE(GDATA%UVW(3, 3))
            ALLOCATE(GDATA%WEIGHTS(3))

            GDATA%N_POINTS = GAUSS_3PT
            GDATA%UVW(1, 1) = 2.0/3.0; GDATA%UVW(2, 1) = 1.0/6.0; GDATA%UVW(3, 1) = 1.0/6.0
            GDATA%UVW(1, 2) = 1.0/6.0; GDATA%UVW(2, 2) = 2.0/3.0; GDATA%UVW(3, 2) = 1.0/6.0
            GDATA%UVW(1, 3) = 1.0/6.0; GDATA%UVW(2, 3) = 1.0/6.0; GDATA%UVW(3, 3) = 2.0/3.0
            
            GDATA%WEIGHTS(1) = 1.0/3.0; GDATA%WEIGHTS(2) = 1.0/3.0; GDATA%WEIGHTS(3) = 1.0/3.0 
        END IF

        IF(N_POINTS == GAUSS_4PT) THEN
            ALLOCATE(GDATA%UVW(3, 4))
            ALLOCATE(GDATA%WEIGHTS(4))

            GDATA%N_POINTS = GAUSS_4PT
            GDATA%UVW(1, 1) = 1.0/3.0; GDATA%UVW(2, 1) = 1.0/3.0; GDATA%UVW(3, 1) = 1.0/3.0
            GDATA%UVW(1, 2) = 0.6; GDATA%UVW(2, 2) = 0.2; GDATA%UVW(3, 2) = 0.2
            GDATA%UVW(1, 3) = 0.2; GDATA%UVW(2, 3) = 0.6; GDATA%UVW(3, 3) = 0.2
            GDATA%UVW(1, 4) = 0.2; GDATA%UVW(2, 4) = 0.2; GDATA%UVW(3, 4) = 0.6
            GDATA%WEIGHTS(1) = -27.0/48.0; GDATA%WEIGHTS(2) = 25.0/48.0; GDATA%WEIGHTS(3) = 25.0/48.0; GDATA%WEIGHTS(4) = 25.0/48.0 
        END IF       

    END SUBROUTINE

    ! 通过面积坐标，得到编号为TRI_ID的三角形的全局坐标
    SUBROUTINE GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_ID, GDATA, GLOBAL_PTS)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        INTEGER, INTENT(IN) :: TRI_ID
        REAL, INTENT(OUT) :: GLOBAL_PTS(3, GDATA%N_POINTS) ! GDATA%N_POINTS代表积分点个数
        INTEGER :: I = 0
        INTEGER :: N1, N2, N3

        N1 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(1)
        N2 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(2)
        N3 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(3)
            
        DO I = 1, GDATA%N_POINTS
            GLOBAL_PTS(1, I) = GDATA%UVW(1, I)*MESH%NODES(N1)%X + &
                                GDATA%UVW(2, I)*MESH%NODES(N2)%X + & 
                                GDATA%UVW(3, I)*MESH%NODES(N3)%X
            GLOBAL_PTS(2, I) = GDATA%UVW(1, I)*MESH%NODES(N1)%Y + &
                                GDATA%UVW(2, I)*MESH%NODES(N2)%Y + & 
                                GDATA%UVW(3, I)*MESH%NODES(N3)%Y                
            GLOBAL_PTS(3, I) = GDATA%UVW(1, I)*MESH%NODES(N1)%Z + &
                                GDATA%UVW(2, I)*MESH%NODES(N2)%Z + & 
                                GDATA%UVW(3, I)*MESH%NODES(N3)%Z    
        END DO
        
    END SUBROUTINE

    ! 检验示例函数1：对常数函数f=1在三角形TRI_ID上做高斯积分
    REAL FUNCTION INTEGRATE_CONSTANT_ON_TRI(MESH, TRI_ID, GDATA) RESULT(VAL)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        INTEGER, INTENT(IN) :: TRI_ID

        VAL = MESH%TRIANGLES(TRI_ID)%AREA * SUM(GDATA%WEIGHTS)

    END FUNCTION

    ! 检验示例函数2：对函数f(X,Y,Z)=X在三角形TRI_ID上做高斯积分
    REAL FUNCTION INTEGRATE_LINEAR_ON_TRI(MESH, TRI_ID, GDATA) RESULT(VAL)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        INTEGER, INTENT(IN) :: TRI_ID
        INTEGER :: I = 0
        REAL :: T(GDATA%N_POINTS)
        REAL :: GLOBAL_PTS(3, GDATA%N_POINTS)

        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_ID, GDATA, GLOBAL_PTS)
        DO I = 1, GDATA%N_POINTS
            T(I) = GDATA%WEIGHTS(I) * GLOBAL_PTS(1,I)
        END DO
    
        VAL = MESH%TRIANGLES(TRI_ID)%AREA * SUM(T)
    
    END FUNCTION

END MODULE