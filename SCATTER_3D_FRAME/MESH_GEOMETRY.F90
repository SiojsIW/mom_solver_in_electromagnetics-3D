MODULE MESH_GEOMETRY
    USE EM_TYPES
    IMPLICIT NONE
CONTAINS

    ! 初始化网格，设置节点数，三角形个数
    SUBROUTINE INIT_MESH_3D(MESH, N_NODE, N_TRI)
        TYPE(MESH_3D), INTENT(INOUT) :: MESH
        INTEGER, INTENT(IN) :: N_NODE, N_TRI

        MESH%NUM_NODE = N_NODE
        MESH%NUM_TRIANGLE = N_TRI

        IF (ALLOCATED(MESH%NODES)) DEALLOCATE(MESH%NODES)
        IF (ALLOCATED(MESH%TRIANGLES)) DEALLOCATE(MESH%TRIANGLES)
        IF (ALLOCATED(MESH%EDGES)) DEALLOCATE(MESH%EDGES)

        ALLOCATE(MESH%NODES(N_NODE))
        ALLOCATE(MESH%TRIANGLES(N_TRI))
        ALLOCATE(MESH%EDGES(N_TRI * 3))
    END SUBROUTINE

    ! 计算三角形外法向量与面积
    SUBROUTINE CALC_NORMAL_AREA(MESH, TRI_ID, NORMAL, AREA)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_ID
        REAL, INTENT(OUT) :: NORMAL(3), AREA
        INTEGER :: N1, N2, N3, I
        REAL :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3
        REAL :: V1(3), V2(3), CROSS(3)  ! 两条边向量

        N1 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(1)
        N2 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(2)
        N3 = MESH%TRIANGLES(TRI_ID)%VERTEX_3D(3)

        X1 = MESH%NODES(N1)%X
        Y1 = MESH%NODES(N1)%Y
        Z1 = MESH%NODES(N1)%Z
        X2 = MESH%NODES(N2)%X
        Y2 = MESH%NODES(N2)%Y
        Z2 = MESH%NODES(N2)%Z
        X3 = MESH%NODES(N3)%X
        Y3 = MESH%NODES(N3)%Y
        Z3 = MESH%NODES(N3)%Z

        V1(1) = X2 - X1
        V1(2) = Y2 - Y1
        V1(3) = Z2 - Z1
        V2(1) = X3 - X1
        V2(2) = Y3 - Y1
        V2(3) = Z3 - Z1
        ! axb=(a2b3-a3b2,a3b1-a1b3,a1b2-a2b1)
        CROSS(1) = V1(2)*V2(3) - V1(3)*V2(2)
        CROSS(2) = V1(3)*V2(1) - V1(1)*V2(3)
        CROSS(3) = V1(1)*V2(2) - V1(2)*V2(1)

        AREA = 0.5*SQRT(CROSS(1)**2 + CROSS(2)**2 + CROSS(3)**2)
        IF (AREA > 0) THEN
            DO I = 1, 3
                NORMAL(I) = CROSS(I) / SQRT(CROSS(1)**2 + CROSS(2)**2 + CROSS(3)**2)
            END DO
        ELSE
            PRINT *, "AREA IS WRONG !"
        END IF
    END SUBROUTINE

    ! 检查法向方向，调用此子程序来保证所有法向一致朝外（目前只能对于凸体）。
    ! 调用必须在UPDATE_ALL_TRIANGLES_3D之后，在EXTRACT_EDGES_FROM_MESH之前
    SUBROUTINE VERIFY_OUTWARD_NORMALS(MESH)
        TYPE(MESH_3D), INTENT(INOUT) :: MESH
        REAL :: VEC(3), CENTER(3), DOT
        INTEGER :: I
        INTEGER :: TEMP_ID
        
        ! 整个物体的几何中心。所有点坐标求和的平均值
        CENTER = 0.0

        DO I = 1, MESH%NUM_NODE
            CENTER(1) = CENTER(1) + MESH%NODES(I)%X
            CENTER(2) = CENTER(2) + MESH%NODES(I)%Y
            CENTER(3) = CENTER(3) + MESH%NODES(I)%Z
        END DO

        CENTER = CENTER / MESH%NUM_NODE

        DO I = 1, MESH%NUM_TRIANGLE
            VEC(1) = MESH%TRIANGLES(I)%CENTROID(1) - CENTER(1)
            VEC(2) = MESH%TRIANGLES(I)%CENTROID(2) - CENTER(2)
            VEC(3) = MESH%TRIANGLES(I)%CENTROID(3) - CENTER(3)
            DOT = DOT_PRODUCT(VEC, MESH%TRIANGLES(I)%NORMAL)

            IF (DOT < 0) THEN
                TEMP_ID = MESH%TRIANGLES(I)%VERTEX_3D(3)
                MESH%TRIANGLES(I)%VERTEX_3D(3) = MESH%TRIANGLES(I)%VERTEX_3D(2)
                MESH%TRIANGLES(I)%VERTEX_3D(2) = TEMP_ID
                CALL CALC_NORMAL_AREA(MESH, I, MESH%TRIANGLES(I)%NORMAL, MESH%TRIANGLES(I)%AREA)
            END IF
        END DO

    END SUBROUTINE

    ! 更新三角形的重心，面积，外法线向量
    SUBROUTINE UPDATE_ALL_TRIANGLES_3D(MESH)
        TYPE(MESH_3D), INTENT(INOUT) :: MESH
        INTEGER :: N1, N2, N3, I
        REAL :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3

        DO I = 1, MESH%NUM_TRIANGLE
            CALL CALC_NORMAL_AREA(MESH, I, MESH%TRIANGLES(I)%NORMAL, MESH%TRIANGLES(I)%AREA)
            N1 = MESH%TRIANGLES(I)%VERTEX_3D(1)
            N2 = MESH%TRIANGLES(I)%VERTEX_3D(2)
            N3 = MESH%TRIANGLES(I)%VERTEX_3D(3)

            X1 = MESH%NODES(N1)%X
            Y1 = MESH%NODES(N1)%Y
            Z1 = MESH%NODES(N1)%Z  
            X2 = MESH%NODES(N2)%X
            Y2 = MESH%NODES(N2)%Y
            Z2 = MESH%NODES(N2)%Z
            X3 = MESH%NODES(N3)%X
            Y3 = MESH%NODES(N3)%Y
            Z3 = MESH%NODES(N3)%Z

            MESH%TRIANGLES(I)%CENTROID(1) = (X1 + X2 + X3) / 3.0
            MESH%TRIANGLES(I)%CENTROID(2) = (Y1 + Y2 + Y3) / 3.0
            MESH%TRIANGLES(I)%CENTROID(3) = (Z1 + Z2 + Z3) / 3.0
        END DO
    END SUBROUTINE

    ! 遍历所有三角形，提取边、去重、建立共享关系。
    ! 展开 + 暴力匹配，不用排序，直接用双重循环找相同端点。
    SUBROUTINE EXTRACT_EDGES_FROM_MESH(MESH)
        TYPE(MESH_3D), INTENT(INOUT) :: MESH
        INTEGER :: TEMP_EDGES(3, 1000)  ! 边数组：第1-2行为端点，第3行为所属三角形编号
        INTEGER :: I, J, K, COL
        INTEGER :: V1, V2, V3
        LOGICAL :: USED(1000)  ! 标记是不是公共边,这条边有没有被收录

        COL = 0
        ! 提取边
        DO I = 1, MESH%NUM_TRIANGLE  
            V1 = MESH%TRIANGLES(I)%VERTEX_3D(1)
            V2 = MESH%TRIANGLES(I)%VERTEX_3D(2)
            V3 = MESH%TRIANGLES(I)%VERTEX_3D(3)

            COL = COL + 1
            TEMP_EDGES(1, COL) = V1
            TEMP_EDGES(2, COL) = V2
            TEMP_EDGES(3, COL) = I
            CALL NORMALIZE_EDGE(TEMP_EDGES(1, COL), TEMP_EDGES(2, COL)) ! 保证V2>V1

            COL = COL + 1
            TEMP_EDGES(1, COL) = V2
            TEMP_EDGES(2, COL) = V3
            TEMP_EDGES(3, COL) = I
            CALL NORMALIZE_EDGE(TEMP_EDGES(1, COL), TEMP_EDGES(2, COL)) ! 保证V3>V2

            COL = COL + 1
            TEMP_EDGES(1, COL) = V1
            TEMP_EDGES(2, COL) = V3
            TEMP_EDGES(3, COL) = I
            CALL NORMALIZE_EDGE(TEMP_EDGES(1, COL), TEMP_EDGES(2, COL)) ! 保证V3>V1
        END DO

        ! 根据提取的边数分配 EDGES 数组
        IF (ALLOCATED(MESH%EDGES)) DEALLOCATE(MESH%EDGES)
        ALLOCATE(MESH%EDGES(COL))

        USED(1:COL) = .FALSE.
        K = 0  ! 公共边计数器

        ! 扫描去重并建立共享关系
        DO I = 1, COL
            IF (USED(I)) CYCLE
            K = K + 1

            MESH%EDGES(K)%V1_ID = TEMP_EDGES(1, I)
            MESH%EDGES(K)%V2_ID = TEMP_EDGES(2, I)
            MESH%EDGES(K)%NUMBER_SHARED_TPI = 1 ! 共享这条边的三角形的数量
            MESH%EDGES(K)%SHARED_TRI_IDS(1) = TEMP_EDGES(3, I) ! 共享这条边的三角形的ID
            MESH%EDGES(K)%LOCAL_EDGE_IDX(1) = GET_EDGE_LOCAL_INDEX(MESH%TRIANGLES(TEMP_EDGES(3, I)), &
                                                MESH%EDGES(K)%V1_ID, MESH%EDGES(K)%V2_ID) ! 这条边在第一个三角形中的局部编号（1,2,3）

            ! 在边I的后面一个个查，查找共享同一条边的另一个三角形
            DO J = I + 1, COL
                IF (USED(J)) CYCLE
                IF (TEMP_EDGES(1, I) == TEMP_EDGES(1, J) .AND. &
                    TEMP_EDGES(2, I) == TEMP_EDGES(2, J)) THEN
                    MESH%EDGES(K)%NUMBER_SHARED_TPI = 2 
                    MESH%EDGES(K)%SHARED_TRI_IDS(2) = TEMP_EDGES(3, J)
                    MESH%EDGES(K)%LOCAL_EDGE_IDX(2) = GET_EDGE_LOCAL_INDEX(MESH%TRIANGLES(TEMP_EDGES(3, J)), &
                                                        MESH%EDGES(K)%V1_ID, MESH%EDGES(K)%V2_ID) ! 这条边在第二个三角形中的局部编号（1,2,3）

                    USED(J) = .TRUE.
                    EXIT
                END IF
            END DO

            USED(I) = .TRUE.
        END DO

        MESH%NUM_EDGE = K

    END SUBROUTINE

    SUBROUTINE NORMALIZE_EDGE(V1_TEMP, V2_TEMP)  
        INTEGER, INTENT(INOUT) :: V1_TEMP, V2_TEMP
        INTEGER :: TEMP
        IF (V1_TEMP > V2_TEMP) THEN
            TEMP = V1_TEMP
            V1_TEMP = V2_TEMP
            V2_TEMP = TEMP
        END IF
    END SUBROUTINE

    ! 局部边号约定：三角形顶点按顺序 (V1,V2,V3)，
    ! 给两个端点所对应的三角形进行边的局部编号（1,2,3）,局部边号 1 对应 (V1,V2)，2 对应 (V2,V3)，3 对应 (V3,V1)
    FUNCTION GET_EDGE_LOCAL_INDEX(TRI, V1, V2) RESULT(IDX)
        TYPE(TRIANGLE_3D), INTENT(IN) :: TRI
        INTEGER, INTENT(IN) :: V1, V2
        INTEGER :: IDX

        IF ((TRI%VERTEX_3D(1) == V1 .AND. TRI%VERTEX_3D(2) == V2) .OR. &
            (TRI%VERTEX_3D(2) == V1 .AND. TRI%VERTEX_3D(1) == V2)) THEN
            IDX = 1
        END IF
        IF ((TRI%VERTEX_3D(2) == V1 .AND. TRI%VERTEX_3D(3) == V2) .OR. &
            (TRI%VERTEX_3D(3) == V1 .AND. TRI%VERTEX_3D(2) == V2)) THEN
            IDX = 2
        END IF
        IF ((TRI%VERTEX_3D(3) == V1 .AND. TRI%VERTEX_3D(1) == V2) .OR. &
            (TRI%VERTEX_3D(1) == V1 .AND. TRI%VERTEX_3D(3) == V2)) THEN
            IDX = 3
        END IF
        IF (IDX == 0) THEN
            PRINT *, "GET_EDGE_LOCAL_INDEX: WRONG!"
        END IF
    END FUNCTION

    ! 更新边的长度与中点
    SUBROUTINE UPDATE_EDGE_GEOMETRY(MESH)
        TYPE(MESH_3D), INTENT(INOUT) :: MESH
        INTEGER :: N1, N2
        REAL :: X1 = 0.0, X2 = 0.0, Y1 = 0.0, Y2 = 0.0, Z1 = 0.0, Z2 = 0.0
        INTEGER :: I = 0

        DO I = 1, MESH%NUM_EDGE
            N1 = MESH%EDGES(I)%V1_ID
            N2 = MESH%EDGES(I)%V2_ID

            X1 = MESH%NODES(N1)%X
            Y1 = MESH%NODES(N1)%Y
            Z1 = MESH%NODES(N1)%Z
            X2 = MESH%NODES(N2)%X
            Y2 = MESH%NODES(N2)%Y
            Z2 = MESH%NODES(N2)%Z

            MESH%EDGES(I)%LENGTH = SQRT((X2 - X1) ** 2 + (Y2 - Y1) ** 2 + (Z2 - Z1) ** 2)

            MESH%EDGES(I)%MIDPOINT(1) = (X1 + X2) / 2.0
            MESH%EDGES(I)%MIDPOINT(2) = (Y1 + Y2) / 2.0
            MESH%EDGES(I)%MIDPOINT(3) = (Z1 + Z2) / 2.0
        END DO
    END SUBROUTINE

    ! 计算两个三角形重心的欧氏距离，来判断是远区/近区/重合
    REAL FUNCTION DIST_TRI_CENTROID(MESH, TRI_ID1, TRI_ID2) RESULT(D)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_ID1, TRI_ID2

        D = SQRT((MESH%TRIANGLES(TRI_ID1)%CENTROID(1) - MESH%TRIANGLES(TRI_ID2)%CENTROID(1)) ** 2 + &
                 (MESH%TRIANGLES(TRI_ID1)%CENTROID(2) - MESH%TRIANGLES(TRI_ID2)%CENTROID(2)) ** 2 + &
                 (MESH%TRIANGLES(TRI_ID1)%CENTROID(3) - MESH%TRIANGLES(TRI_ID2)%CENTROID(3)) ** 2)
    
    END FUNCTION

END MODULE