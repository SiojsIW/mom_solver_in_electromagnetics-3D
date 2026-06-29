! （1）基于PRACTICE 8，优化三个地方：
!           1.共面时K算子部分理论为0：添加IS_COPLANAR_TRI_PAIR程序，判断两三角形是否共面，如果共面则K算子相关项为0
!           2.1/2<fm,fn>在临边是不为0：在两个RWG支撑域完全不重合（四个三角形都不重合）时才为0
!           3.添加SUBROUTINE VERIFY_OUTWARD_NORMALS(MESH)。由于涉及外方向向量参与计算，应该检查是否都为外法线方向。
! （2）CFIE组合（左端项）+ 右端项
! ---------------------------------------------------------------2026/6/26---------------------

PROGRAM TETRAHENDRON_MESH
    USE EM_TYPES
    USE MESH_GEOMETRY
    USE NUMERICAL_INTEGRATION
    USE GREEN_FUNCTIONS
    USE RWG_BASIS_BUILD
    USE Z_MATRIX
    USE RHS
    IMPLICIT NONE
    
    REAL, PARAMETER :: ETA0 = 120.0 * PI
    INTEGER :: I, J, M, N ! 用于循环
    REAL, ALLOCATABLE ::GLOBAL_PTS(:, :) ! 必须用可分配大小的数组，因为目前第二维大小未知。后面再具体分配内存大小
    REAL :: VAL = 0.0 ! 高斯积分值
    REAL :: LAMBDA = 1.0 ! 波长
    REAL :: K ! 波数
    REAL :: DIST ! 两个三角形之间的欧氏距离
    REAL :: F_VAL_POS(3), F_VAL_NEG(3) ! 正负三角形的f（r）的值
    COMPLEX :: G_APPROX  ! 用格林函数近似公式计算的两个三角形Z_MN，用于验证结果对不对
    COMPLEX :: Z_MN
    COMPLEX, ALLOCATABLE :: Z_E(:, :), Z_M(:, :), Z_CFIE(:, :)
    REAL :: ALPHA
    ! 入射波平面参数定义
    REAL :: K_HAT(3) = [0.0, 0.0, 1.0]
    REAL :: E_POL(3) = [1.0, 0.0, 0.0]
    REAL :: E0 = 1.0
    COMPLEX, ALLOCATABLE :: V_M_E(:), V_M_M(:), V_CFIE(:)

    TYPE(MESH_3D) :: MESH
    TYPE(GAUSS_TRI_DATA) :: GDATA

    K = 2 * PI / LAMBDA
    
    CALL INIT_MESH_3D(MESH, 4, 4)
    CALL INIT_GAUSS_TRI(GAUSS_3PT, GDATA)

    ! 填充节点（正四面体顶点，边长为 1）
    MESH%NODES(1)%ID = 1; MESH%NODES(1)%X = 0.0; MESH%NODES(1)%Y = 0.0; MESH%NODES(1)%Z = 0.0
    MESH%NODES(2)%ID = 2; MESH%NODES(2)%X = 1.0; MESH%NODES(2)%Y = 0.0; MESH%NODES(2)%Z = 0.0
    MESH%NODES(3)%ID = 3; MESH%NODES(3)%X = 0.5; MESH%NODES(3)%Y = SQRT(3.0)/2.0; MESH%NODES(3)%Z = 0.0
    MESH%NODES(4)%ID = 4; MESH%NODES(4)%X = 0.5; MESH%NODES(4)%Y = SQRT(3.0)/6.0; MESH%NODES(4)%Z = SQRT(6.0)/3.0

    ! 填充 4 个三角形面（顶点编号按右手定则外法向朝外）
    MESH%TRIANGLES(1)%VERTEX_3D(1) = 1; MESH%TRIANGLES(1)%VERTEX_3D(2) = 3; MESH%TRIANGLES(1)%VERTEX_3D(3) = 2
    MESH%TRIANGLES(2)%VERTEX_3D(1) = 1; MESH%TRIANGLES(2)%VERTEX_3D(2) = 2; MESH%TRIANGLES(2)%VERTEX_3D(3) = 4
    MESH%TRIANGLES(3)%VERTEX_3D(1) = 2; MESH%TRIANGLES(3)%VERTEX_3D(2) = 3; MESH%TRIANGLES(3)%VERTEX_3D(3) = 4
    MESH%TRIANGLES(4)%VERTEX_3D(1) = 3; MESH%TRIANGLES(4)%VERTEX_3D(2) = 1; MESH%TRIANGLES(4)%VERTEX_3D(3) = 4

    CALL UPDATE_ALL_TRIANGLES_3D(MESH)
    CALL VERIFY_OUTWARD_NORMALS(MESH)
    CALL EXTRACT_EDGES_FROM_MESH(MESH)
    CALL UPDATE_EDGE_GEOMETRY(MESH)
    CALL BUILD_RWG_BASIS(MESH)

    PRINT *, "总边数：", MESH%NUM_EDGE
    
    DO I = 1, MESH%NUM_EDGE

        PRINT *, "第", I, "条边的第一个顶点ID: ", MESH%EDGES(I)%V1_ID, "第", I, "条边的第二个顶点ID: ", MESH%EDGES(I)%V2_ID
        PRINT *, "第", I, "条边的长度: ", MESH%EDGES(I)%LENGTH
        PRINT *, "共享第", I, "条边的三角形的个数: ", MESH%EDGES(I)%NUMBER_SHARED_TPI

    END DO

    DO I = 1, MESH%NUM_EDGE
        IF (MESH%EDGES(I)%NUMBER_SHARED_TPI == 2) THEN
            J = J + 1
        END IF
    END DO

    PRINT *, "RWG边数（被两个三角形共享的）为：", J
    
    DO I = 1, MESH%NUM_EDGE
        IF (MESH%EDGES(I)%NUMBER_SHARED_TPI == 2) THEN
            PRINT *, "第", I, "个边是公共边，他属于这两个三角形：", &
                    MESH%EDGES(I)%SHARED_TRI_IDS(1), MESH%EDGES(I)%SHARED_TRI_IDS(2), &
                    "在三角形中的局部编号分别为：", MESH%EDGES(I)%LOCAL_EDGE_IDX(1), MESH%EDGES(I)%LOCAL_EDGE_IDX(2)
        END IF
    END DO

    DO I = 1, MESH%NUM_TRIANGLE
        PRINT *, "三角形", I, "面积：", MESH%TRIANGLES(I)%AREA
        PRINT *, "法向：", MESH%TRIANGLES(I)%NORMAL(1), &
                           MESH%TRIANGLES(I)%NORMAL(2), &
                           MESH%TRIANGLES(I)%NORMAL(3)
        PRINT *, "重心：", MESH%TRIANGLES(I)%CENTROID(1), &
                           MESH%TRIANGLES(I)%CENTROID(2), &
                           MESH%TRIANGLES(I)%CENTROID(3)
    END DO

    ! 打印每个三角形的高斯积分点的全局坐标
    ALLOCATE(GLOBAL_PTS(3, GDATA%N_POINTS))
    DO I = 1, MESH%NUM_TRIANGLE
        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, I, GDATA, GLOBAL_PTS)
            PRINT *, "使用高斯", GDATA%N_POINTS, "点求积"
            
            PRINT *, "三角形", I, "的所有高斯积分点的全局X坐标为：", (GLOBAL_PTS(1, N), N = 1, GDATA%N_POINTS)
            PRINT *, "三角形", I, "的所有高斯积分点的全局Y坐标为：", (GLOBAL_PTS(2, N), N = 1, GDATA%N_POINTS)
            PRINT *, "三角形", I, "的所有高斯积分点的全局Z坐标为：", (GLOBAL_PTS(3, N), N = 1, GDATA%N_POINTS)

        PRINT *, "三角形", I, "的高斯积分点全局坐标已打印完毕"
    END DO
    DEALLOCATE(GLOBAL_PTS)

    ! 对常数函数f=1在每个三角形上做高斯积分 FUNCTION INTEGRATE_LINEAR_ON_TRI(MESH, TRI_ID, GDATA) RESULT(VAL)
    DO I = 1, MESH%NUM_TRIANGLE
        VAL = INTEGRATE_CONSTANT_ON_TRI(MESH, I, GDATA)
        PRINT *, "三角形", I, "在常数函数f=1上的", GDATA%N_POINTS, "点高斯积分值为：", VAL
    END DO
    ! 对函数f(Y,Y,Z)=X在每个三角形上做高斯积分
    DO I = 1, MESH%NUM_TRIANGLE
        VAL = INTEGRATE_LINEAR_ON_TRI(MESH, I, GDATA)
        PRINT *, "三角形", I, "在函数f(Y,Y,Z)=X上的", GDATA%N_POINTS, "点高斯积分值为：", VAL
    END DO

    ! 计算所有三角形对之间的格林函数积分
    DO I = 1, MESH%NUM_TRIANGLE
        DO J = 1, MESH%NUM_TRIANGLE
            IF (I == J) THEN
                PRINT *, "自项：暂时跳过"
                CYCLE
            END IF
            CALL CALC_GREEN_MATRIX_ELEMENT(MESH, I, J, GDATA, K, Z_MN)
            PRINT *, "阻抗元素Z_MN(", I, J, ")的模值为：", ABS(Z_MN)
        END DO
    END DO
    ! 验证

    DIST = DIST_TRI_CENTROID(MESH, 1, 2) ! 计算第一个和第二个三角形的欧氏距离
    G_APPROX = CEXP((0.0, -1.0) * K * DIST) / (4.0 * PI * DIST)

    PRINT *, "第一个三角形和第二个三角形的格林函数的模值为：", &
            ABS(G_APPROX * MESH%TRIANGLES(1)%AREA * MESH%TRIANGLES(2) % AREA)
            
    ! RWG基函数的验证：
    PRINT *, "每个RWG基函数的信息："
    DO I = 1, MESH%RWG_NUM
        PRINT *, "公共边的长度为：", MESH%RWG_BASES(I)%LENGTH
        PRINT *, "公共边的ID为：", MESH%RWG_BASES(I)%EDGE_ID
        PRINT *, "公共边", MESH%RWG_BASES(I)%EDGE_ID, "的正三角形编号为：", MESH%RWG_BASES(I)%POS_TRI_ID
        PRINT *, "公共边", MESH%RWG_BASES(I)%EDGE_ID, "的负三角形编号为：", MESH%RWG_BASES(I)%NEG_TRI_ID
        PRINT *, "公共边", MESH%RWG_BASES(I)%EDGE_ID, "的正三角形f（r）的系数为：", MESH%RWG_BASES(I)%POS_COEF
        PRINT *, "公共边", MESH%RWG_BASES(I)%EDGE_ID, "的负三角形f（r）的系数为：", MESH%RWG_BASES(I)%NEG_COEF
    END DO
    
    PRINT *, "第一个基函数的信息："
    PRINT *, "第一个基函数公共边中点坐标为：", &
                (MESH%EDGES(MESH%RWG_BASES(1)%EDGE_ID)%MIDPOINT(I), I = 1, 3)
    
    CALL EVAL_RWG_BASIS(MESH, MESH%RWG_BASES(1), "POS", &
                        MESH%TRIANGLES(MESH%RWG_BASES(1)%POS_TRI_ID)%CENTROID, F_VAL_POS)
    PRINT *, "第一个基函数正三角形重心的f（r）的值为：", F_VAL_POS

    ! 验证：公共边中点
    PRINT *, " 公共边中点验证: "
    CALL EVAL_RWG_BASIS(MESH, MESH%RWG_BASES(1), "POS", &
                        MESH%EDGES(MESH%RWG_BASES(1)%EDGE_ID)%MIDPOINT, F_VAL_POS)
    CALL EVAL_RWG_BASIS(MESH, MESH%RWG_BASES(1), "NEG", &
                        MESH%EDGES(MESH%RWG_BASES(1)%EDGE_ID)%MIDPOINT, F_VAL_NEG)

    PRINT *, "POS: ", F_VAL_POS
    PRINT *, "NEG: ", F_VAL_NEG
    PRINT *, "模: ", SQRT(SUM(F_VAL_POS**2)), SQRT(SUM(F_VAL_NEG**2))

    ! 打印EFIE阻抗元素
    DO I = 1, MESH%RWG_NUM
        DO J = 1, MESH%RWG_NUM
            IF (I == J) THEN
                PRINT *, "自项，暂时跳过"
                CYCLE
            END IF
            CALL CALC_EFIE_MATRIX_ELEMENT(MESH, MESH%RWG_BASES(I), &
                        MESH%RWG_BASES(J), GDATA, K, ETA0, Z_MN)    
            PRINT *, "EFIE Z(", I, J, ")模值: ", ABS(Z_MN)
        END DO
    END DO

    ! 打印MFIE阻抗元素
    DO M = 1, MESH%RWG_NUM
        DO N = 1, MESH%RWG_NUM
            IF (M == N) THEN
                PRINT *, "自项，暂时跳过"
                CYCLE
            END IF
            CALL CALC_MFIE_MATRIX_ELEMENT(MESH, MESH%RWG_BASES(M), &
                        MESH%RWG_BASES(N), GDATA, K, Z_MN)    
            PRINT *, "MFIE Z(", M, N, ")模值: ", ABS(Z_MN)
        END DO
    END DO
    
    ! 矩阵存储与CFIE组合
    IF(ALLOCATED(Z_E)) DEALLOCATE(Z_E)
    IF(ALLOCATED(Z_M)) DEALLOCATE(Z_M)
    IF(ALLOCATED(Z_CFIE)) DEALLOCATE(Z_CFIE)
    ALLOCATE(Z_E(MESH%RWG_NUM, MESH%RWG_NUM))
    ALLOCATE(Z_M(MESH%RWG_NUM, MESH%RWG_NUM))
    ALLOCATE(Z_CFIE(MESH%RWG_NUM, MESH%RWG_NUM))

    ! 一定要初始化！！！不然会数值爆炸
    Z_E = (0.0, 0.0)
    Z_M = (0.0, 0.0)
    Z_CFIE = (0.0, 0.0)

    DO M = 1, MESH%RWG_NUM
        DO N = 1, MESH%RWG_NUM
            IF (M == N) THEN
                PRINT *, "自项，暂时跳过"
                CYCLE
            END IF

            CALL CALC_EFIE_MATRIX_ELEMENT(MESH, MESH%RWG_BASES(M), &
                                MESH%RWG_BASES(N), GDATA, K, ETA0, Z_E(M, N))
            CALL CALC_MFIE_MATRIX_ELEMENT(MESH, MESH%RWG_BASES(M), &
                                MESH%RWG_BASES(N), GDATA, K, Z_M(M, N))
            PRINT *, "Z_E(", M, N, ")模值: ", ABS(Z_E(M, N))
            PRINT *, "Z_M(", M, N, ")模值: ", ABS(Z_M(M, N))
        END DO
    END DO

    ! CFIE组合
    ALPHA = 0.5

    DO M = 1, MESH%RWG_NUM
        DO N = 1, MESH%RWG_NUM
            IF (M == N) THEN
                PRINT *, "自项，暂时跳过"
                CYCLE
            END IF
            Z_CFIE(M, N) = ALPHA * Z_E(M, N) + (1.0 - ALPHA) * ETA0 * Z_M(M, N)
            PRINT *, "Z_CFIE(", M, N, ")模值: ", ABS(Z_CFIE(M, N))
        END DO
    END DO
    
    IF (ABS(DOT_PRODUCT(K_HAT, E_POL)) > 1.0E-6) THEN
        PRINT *, "错误：极化方向与传播方向不正交"
        STOP
    END IF

    ! 验证右端项
    IF(ALLOCATED(V_M_E)) DEALLOCATE(V_M_E)
    IF(ALLOCATED(V_M_M)) DEALLOCATE(V_M_M)
    IF(ALLOCATED(V_CFIE)) DEALLOCATE(V_CFIE)
    ALLOCATE(V_M_E(MESH%RWG_NUM))
    ALLOCATE(V_M_M(MESH%RWG_NUM))
    ALLOCATE(V_CFIE(MESH%RWG_NUM))
    DO M = 1, MESH%RWG_NUM
        CALL CALC_EFIE_RHS(MESH, MESH%RWG_BASES(M), GDATA, K, K_HAT, E_POL, E0, ETA0, V_M_E(M))
        CALL CALC_MFIE_RHS(MESH, MESH%RWG_BASES(M), GDATA, K, K_HAT, E_POL, E0, ETA0, V_M_M(M))
        V_CFIE(M) = ALPHA * V_M_E(M) + (1.0 - ALPHA) * ETA0 * V_M_M(M)
        PRINT *, "RWG", M, &
                "V_M_E模:", ABS(V_M_E(M)), &
                "V_M_M模:", ABS(V_M_M(M)), &
               "V_CFIE模:", ABS(V_CFIE(M))
    END DO
END PROGRAM