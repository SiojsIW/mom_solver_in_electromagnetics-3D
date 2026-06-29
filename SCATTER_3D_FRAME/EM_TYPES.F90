MODULE EM_TYPES
    IMPLICIT NONE

    TYPE NODE_3D
        INTEGER :: ID = 0
        REAL :: X = 0.0, Y = 0.0, Z = 0.0
    END TYPE

    TYPE TRIANGLE_3D
        INTEGER :: VERTEX_3D(3) = 0  ! 三个节点编号
        REAL :: AREA = 0.0
        REAL :: NORMAL(3) = 0.0  ! 单位外法向向量
        REAL :: CENTROID(3) = 0.0    ! (重心/形心)
    END TYPE

    TYPE EDGE_3D
        INTEGER :: V1_ID = 0, V2_ID = 0  ! 强制 V1_ID < V2_ID
        REAL :: LENGTH = 0.0
        REAL :: MIDPOINT(3) = 0.0
        INTEGER :: NUMBER_SHARED_TPI = 0  ! 共享该边的三角形数量
        INTEGER :: SHARED_TRI_IDS(2) = 0  ! 共享该边的三角形的全局编号
        INTEGER :: LOCAL_EDGE_IDX(2) = 0  ! 该边在对应三角形中的局部编号（1,2,3）
    END TYPE

    TYPE RWG_BASIS 
        INTEGER :: EDGE_ID = 0 ! 所属公共边编号
        INTEGER :: POS_TRI_ID = 0 ! 正三角形编号
        INTEGER :: NEG_TRI_ID = 0 ! 负三角形编号
        INTEGER :: POS_OPP_VERTEX = 0 ! 正三角形对顶点编号
        INTEGER :: NEG_OPP_VERTEX = 0
        REAL :: LENGTH = 0.0  ! 公共边长度
        REAL :: POS_COEF = 0.0  ! 正三角形f（r）的系数，Ln/2*An
        REAL :: NEG_COEF = 0.0
    END TYPE

    TYPE MESH_3D
        INTEGER :: NUM_NODE = 0, NUM_TRIANGLE = 0, NUM_EDGE = 0
        INTEGER :: RWG_NUM = 0  ! 公共边数量
        TYPE(NODE_3D), ALLOCATABLE :: NODES(:)
        TYPE(TRIANGLE_3D), ALLOCATABLE :: TRIANGLES(:)
        TYPE(EDGE_3D), ALLOCATABLE :: EDGES(:)
        TYPE(RWG_BASIS), ALLOCATABLE :: RWG_BASES(:)
    END TYPE

    INTEGER, PARAMETER :: GAUSS_1PT = 1, GAUSS_3PT = 3, GAUSS_4PT = 4, GAUSS_7PT = 7
    REAL, PARAMETER :: PI = 4.0 * ATAN(1.0)

    ! 定义高斯积分数据（几点积分，面积坐标，权重）
    TYPE GAUSS_TRI_DATA
        INTEGER :: N_POINTS = 0
        REAL, ALLOCATABLE :: UVW(:, :)
        REAL, ALLOCATABLE :: WEIGHTS(:)
    END TYPE

END MODULE