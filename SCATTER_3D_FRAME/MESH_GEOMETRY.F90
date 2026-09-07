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
        INTEGER, ALLOCATABLE :: TEMP_EDGES(:, :)  ! 边数组：第1-2行为端点，第3行为所属三角形编号
        INTEGER :: I, J, K, COL
        INTEGER :: V1, V2, V3
        LOGICAL, ALLOCATABLE :: USED(:)  ! 标记是不是公共边,这条边有没有被收录

        COL = 3 * MESH%NUM_TRIANGLE
        ALLOCATE(TEMP_EDGES(3, COL))
        ALLOCATE(USED(COL))

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

    ! 统计字符串中的空白分隔 token 数
    FUNCTION COUNT_TOKENS(LINE) RESULT(N)
        CHARACTER(LEN=*), INTENT(IN) :: LINE
        INTEGER :: N
        INTEGER :: I, LEN_USED
        LOGICAL :: IN_TOKEN
        CHARACTER :: C

        N = 0
        IN_TOKEN = .FALSE.
        LEN_USED = LEN_TRIM(LINE)
        DO I = 1, LEN_USED
            C = LINE(I:I)
            IF (C == ' ' .OR. C == '\t' .OR. C == '\r' .OR. C == '\n') THEN
                IN_TOKEN = .FALSE.
            ELSE
                IF (.NOT. IN_TOKEN) THEN
                    N = N + 1
                    IN_TOKEN = .TRUE.
                END IF
            END IF
        END DO
    END FUNCTION COUNT_TOKENS

    ! 判断一行是否只包含一个正整数（网格文件头行）
    FUNCTION IS_HEADER_COUNT(LINE, N) RESULT(IS_HEADER)
        CHARACTER(LEN=*), INTENT(IN) :: LINE
        INTEGER, INTENT(OUT) :: N
        LOGICAL :: IS_HEADER
        INTEGER :: IOS, LEN_USED
        CHARACTER(LEN=256) :: BUF

        IS_HEADER = .FALSE.
        N = 0
        BUF = ADJUSTL(LINE)
        LEN_USED = LEN_TRIM(BUF)
        IF (LEN_USED == 0) RETURN
        IF (LEN_USED > 16) RETURN   ! 头行通常很短

        ! 头行必须只有一个 token（否则会把三角形顶点行误判为头行）
        IF (COUNT_TOKENS(BUF) /= 1) RETURN

        ! 尝试只读一个整数
        READ(BUF, *, IOSTAT=IOS) N
        IF (IOS /= 0) RETURN
        IF (N <= 0) RETURN

        IS_HEADER = .TRUE.
    END FUNCTION IS_HEADER_COUNT

    ! 统计文件中非空行的数量（从当前位置到文件尾）
    FUNCTION COUNT_NONEMPTY_LINES(UNIT_NUM) RESULT(N)
        INTEGER, INTENT(IN) :: UNIT_NUM
        INTEGER :: N
        CHARACTER(LEN=512) :: LINE
        INTEGER :: IOS

        N = 0
        DO
            READ(UNIT_NUM, '(A)', IOSTAT=IOS) LINE
            IF (IOS /= 0) EXIT
            IF (LEN_TRIM(LINE) > 0) N = N + 1
        END DO
    END FUNCTION COUNT_NONEMPTY_LINES


    ! 从文件读取网格（point.txt + triangle_point.txt 格式）
    ! 支持两种格式：
    !   带头行：第一行=数量，之后每行: ID X Y Z 或 ID V1 V2 V3
    !   无头行：每行直接 X Y Z 或 V1 V2 V3，ID 按顺序 1,2,3,... 分配
    SUBROUTINE READ_MESH_FROM_FILES(MESH, NODE_FILE, TRI_FILE, IOS)
        TYPE(MESH_3D), INTENT(INOUT) :: MESH
        CHARACTER(LEN=*), INTENT(IN) :: NODE_FILE, TRI_FILE
        INTEGER, INTENT(OUT) :: IOS
        INTEGER :: N_NODES, N_TRI, I, ID, V1, V2, V3
        REAL :: X, Y, Z
        CHARACTER(LEN=512) :: LINE
        LOGICAL :: NODE_HAS_HEADER, TRI_HAS_HEADER

        ! —— 读取节点 ——
        OPEN(UNIT=10, FILE=NODE_FILE, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
        IF (IOS /= 0) THEN
            PRINT *, "错误：无法打开节点文件 ", NODE_FILE
            RETURN
        END IF

        READ(10, '(A)', IOSTAT=IOS) LINE
        IF (IOS /= 0) THEN
            PRINT *, "错误：读取节点文件第一行失败"
            CLOSE(10); RETURN
        END IF

        NODE_HAS_HEADER = IS_HEADER_COUNT(LINE, N_NODES)
        IF (NODE_HAS_HEADER) THEN
            IF (N_NODES <= 0) THEN
                PRINT *, "错误：节点数量非正"
                CLOSE(10); RETURN
            END IF
        ELSE
            ! 无头行：关闭后重新打开并统计非空行数
            CLOSE(10)
            OPEN(UNIT=10, FILE=NODE_FILE, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
            IF (IOS /= 0) THEN
                PRINT *, "错误：重新打开节点文件失败"
                RETURN
            END IF
            N_NODES = COUNT_NONEMPTY_LINES(10)
            IF (N_NODES <= 0) THEN
                PRINT *, "错误：节点文件为空"
                CLOSE(10); RETURN
            END IF
            CLOSE(10)
            ! 再次打开准备读取数据
            OPEN(UNIT=10, FILE=NODE_FILE, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
            IF (IOS /= 0) THEN
                PRINT *, "错误：重新打开节点文件失败"
                RETURN
            END IF
        END IF

        ! —— 读取三角形数量 ——
        OPEN(UNIT=11, FILE=TRI_FILE, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
        IF (IOS /= 0) THEN
            PRINT *, "错误：无法打开三角形文件 ", TRI_FILE
            CLOSE(10); RETURN
        END IF

        READ(11, '(A)', IOSTAT=IOS) LINE
        IF (IOS /= 0) THEN
            PRINT *, "错误：读取三角形文件第一行失败"
            CLOSE(10); CLOSE(11); RETURN
        END IF

        TRI_HAS_HEADER = IS_HEADER_COUNT(LINE, N_TRI)
        IF (TRI_HAS_HEADER) THEN
            IF (N_TRI <= 0) THEN
                PRINT *, "错误：三角形数量非正"
                CLOSE(10); CLOSE(11); RETURN
            END IF
        ELSE
            CLOSE(11)
            OPEN(UNIT=11, FILE=TRI_FILE, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
            IF (IOS /= 0) THEN
                PRINT *, "错误：重新打开三角形文件失败"
                CLOSE(10); RETURN
            END IF
            N_TRI = COUNT_NONEMPTY_LINES(11)
            IF (N_TRI <= 0) THEN
                PRINT *, "错误：三角形文件为空"
                CLOSE(10); CLOSE(11); RETURN
            END IF
            CLOSE(11)
            OPEN(UNIT=11, FILE=TRI_FILE, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
            IF (IOS /= 0) THEN
                PRINT *, "错误：重新打开三角形文件失败"
                CLOSE(10); RETURN
            END IF
        END IF

        ! —— 分配网格 ——
        CALL INIT_MESH_3D(MESH, N_NODES, N_TRI)

        ! —— 读取节点数据 ——
        DO I = 1, N_NODES
            IF (NODE_HAS_HEADER) THEN
                READ(10, *, IOSTAT=IOS) ID, X, Y, Z
            ELSE
                READ(10, *, IOSTAT=IOS) X, Y, Z
                ID = I
            END IF
            IF (IOS /= 0) THEN
                PRINT *, "错误：读取第", I, "个节点失败"
                CLOSE(10); CLOSE(11); RETURN
            END IF
            MESH%NODES(ID)%ID = ID
            MESH%NODES(ID)%X = X
            MESH%NODES(ID)%Y = Y
            MESH%NODES(ID)%Z = Z
        END DO
        CLOSE(10)

        ! —— 读取三角形数据 ——
        DO I = 1, N_TRI
            IF (TRI_HAS_HEADER) THEN
                READ(11, *, IOSTAT=IOS) ID, V1, V2, V3
            ELSE
                READ(11, *, IOSTAT=IOS) V1, V2, V3
                ID = I
            END IF
            IF (IOS /= 0) THEN
                PRINT *, "错误：读取第", I, "个三角形失败"
                CLOSE(11); RETURN
            END IF
            MESH%TRIANGLES(ID)%VERTEX_3D(1) = V1
            MESH%TRIANGLES(ID)%VERTEX_3D(2) = V2
            MESH%TRIANGLES(ID)%VERTEX_3D(3) = V3
        END DO
        CLOSE(11)

        PRINT *, "成功读取 ", N_NODES, " 个节点和 ", N_TRI, " 个三角形"
        IOS = 0
    END SUBROUTINE READ_MESH_FROM_FILES

END MODULE MESH_GEOMETRY
