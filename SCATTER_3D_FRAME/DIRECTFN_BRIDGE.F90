! ============================================================================
! DIRECTFN_BRIDGE.F90
!
! DIRECTFN（C++ 库）的 Fortran 桥接层。
! 通过 ISO_C_BINDING 调用 DIRECTFN_WRAPPER.cpp 中的 extern "C" 接口，
! 为 EFIE 阻抗矩阵提供第三种近场奇异积分处理方法（DIRECTFN 全数值法）。
!
! 物理约定：
!   DIRECTFN 返回的 iss_rwg(i,j) = ∬ f_i^obs · f_j^src * G dS dS'，
!   其中 f_i 为标准 RWG 归一化 l_i/(2A)·(r - v_i)，G = exp(-j k R)/(4πR)；
!   iss_const = ∬ G dS dS'。
!   本模块将其组合成与 Z_MATRIX 远场分支一致的阻抗贡献：
!   Z = s_o·c_o·s_s·c_s·( jωμ·P - (j/ωε)·4·I )，
!   其中 P = iss_rwg(io,jo)/(c_o·c_s)，I = iss_const，
!   c = L/(2A) 为 RWG 系数幅值，s = ±1 为所在正/负三角形的符号。
!
! 处理范围：仅自项(ST)、共边(EA)、共顶点(VA)三类奇异对；
!           非接触近场对仍由上层回退到奇异提取法。
! ============================================================================

MODULE DIRECTFN_BRIDGE
    USE ISO_C_BINDING
    USE EM_TYPES
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: CALC_EFIE_DIRECTFN_PAIR, DIRECTFN_CALL_COUNT
    PUBLIC :: DIRECTFN_N_GAUSS_ST, DIRECTFN_N_GAUSS_EA, DIRECTFN_N_GAUSS_VA
    PUBLIC :: DIRECTFN_DEBUG, DIRECTFN_RESET_CACHE, DIRECTFN_PRINT_STATS

    ! DIRECTFN 一维高斯-勒让德阶数（4 维积分每维 N 点），按邻接类型分别设置。
    ! 单对代价 ∝ N^4，且 ST:EA:VA ≈ 24:9:2，故 VA 可用更低阶数。
    ! 默认值 12 为高精度档；标定后（见 HANDOVER.md）生产可用更低值。
    INTEGER :: DIRECTFN_N_GAUSS_ST = 8
    INTEGER :: DIRECTFN_N_GAUSS_EA = 8
    INTEGER :: DIRECTFN_N_GAUSS_VA = 6
    INTEGER :: DIRECTFN_CALL_COUNT = 0
    ! 调试标志：>0 时在每次调用后打印原始积分值（@@DFNRAW 行，供 Python 参考程序解析）
    INTEGER :: DIRECTFN_DEBUG = 0

    ! ---- tri-pair 级缓存 ----
    ! 矩阵填充按 RWG 对循环，同一 (观察,源) 三角形对会被其上的 3×3 个 RWG 对
    ! 重复调用 ~9 次；缓存使 DIRECTFN 4D 积分对每个唯一三角形对只算一次（~9× 加速）。
    ! 键 = (TRI_O-1)*(NT+1) + TRI_S，直接寻址哈希表；K/ETA0 在一次填充内不变。
    TYPE :: DFN_CACHE_ENTRY
        INTEGER :: TRI_O, TRI_S
        INTEGER :: ADJ
        INTEGER :: DFN_OBS_NODES(3), DFN_SRC_NODES(3)
        COMPLEX(C_DOUBLE_COMPLEX) :: ISS_RWG(9), ISS_CONST
    END TYPE
    INTEGER, ALLOCATABLE :: DFN_HASH(:)            ! 键 → 条目号（0 = 空）
    TYPE(DFN_CACHE_ENTRY), ALLOCATABLE :: DFN_CACHE(:)
    INTEGER :: DFN_CACHE_N = 0, DFN_CACHE_CAP = 0, DFN_HASH_NT = 0
    INTEGER(C_INT64_T) :: DFN_STAT_HIT = 0, DFN_STAT_MISS = 0
    INTEGER :: DFN_STAT_ST = 0, DFN_STAT_EA = 0, DFN_STAT_VA = 0

    INTERFACE
        FUNCTION C_DFN_TRI_ISS(ADJ, PTS, K0, NGAUSS, ISS_RWG, ISS_CONST) &
                BIND(C, NAME='dfn_tri_iss') RESULT(IRC)
            IMPORT :: C_INT, C_DOUBLE, C_DOUBLE_COMPLEX
            INTEGER(C_INT), VALUE :: ADJ
            REAL(C_DOUBLE), INTENT(IN) :: PTS(3, 7)
            REAL(C_DOUBLE), VALUE :: K0
            INTEGER(C_INT), VALUE :: NGAUSS
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: ISS_RWG(9)
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: ISS_CONST(1)
            INTEGER(C_INT) :: IRC
        END FUNCTION C_DFN_TRI_ISS
    END INTERFACE

CONTAINS

    ! 矩阵填充前调用：分配缓存（NT = 三角形数）。同一三角形对只积分一次。
    SUBROUTINE DIRECTFN_RESET_CACHE(NT)
        INTEGER, INTENT(IN) :: NT
        IF (ALLOCATED(DFN_HASH)) DEALLOCATE(DFN_HASH)
        ALLOCATE(DFN_HASH(NT * (NT + 1)))
        DFN_HASH = 0
        DFN_HASH_NT = NT
        IF (ALLOCATED(DFN_CACHE)) DEALLOCATE(DFN_CACHE)
        DFN_CACHE_CAP = 1024
        ALLOCATE(DFN_CACHE(DFN_CACHE_CAP))
        DFN_CACHE_N = 0
        DFN_STAT_HIT = 0; DFN_STAT_MISS = 0
        DFN_STAT_ST = 0; DFN_STAT_EA = 0; DFN_STAT_VA = 0
    END SUBROUTINE DIRECTFN_RESET_CACHE

    ! 矩阵填充后调用：打印缓存命中与唯一三角形对统计
    SUBROUTINE DIRECTFN_PRINT_STATS()
        PRINT '(A, I10)', " DIRECTFN 唯一三角形对:  ", DFN_CACHE_N
        PRINT '(A, I10)', "   其中 ST / EA / VA 对:  ", DFN_STAT_ST
        PRINT '(A, I10)', "                          ", DFN_STAT_EA
        PRINT '(A, I10)', "                          ", DFN_STAT_VA
        PRINT '(A, I10)', " DIRECTFN 缓存命中:      ", DFN_STAT_HIT
        PRINT '(A, I10)', " DIRECTFN 实际积分次数:  ", DFN_STAT_MISS
    END SUBROUTINE DIRECTFN_PRINT_STATS

    ! 与 CALC_EFIE_QBX_PAIR 同级的（观察三角形, 源三角形）对阻抗计算
    SUBROUTINE CALC_EFIE_DIRECTFN_PAIR(MESH, TRI_OBS, TRI_SRC, RWG_OBS, RWG_SRC, K, ETA0, Z_PAIR)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: TRI_OBS, TRI_SRC
        TYPE(RWG_BASIS), INTENT(IN) :: RWG_OBS, RWG_SRC
        REAL, INTENT(IN) :: K, ETA0
        COMPLEX, INTENT(OUT) :: Z_PAIR

        INTEGER :: IDO(3), IDS(3)
        INTEGER :: E1, E2, OBS3, SRC3, OPP_O, OPP_S
        INTEGER :: DFN_OBS_NODES(3), DFN_SRC_NODES(3)
        INTEGER :: I, J, N_COMMON, IO, JO, NPTS
        INTEGER :: KEY, CIDX, ADJ_I, N_GAUSS_DFN
        LOGICAL :: CACHE_HIT
        REAL(C_DOUBLE) :: PTS(3, 7)
        REAL(C_DOUBLE) :: COEF_O, COEF_S, S_O, S_S, K_D
        COMPLEX(C_DOUBLE_COMPLEX) :: ISS_RWG(9), ISS_CONST(1)
        COMPLEX(C_DOUBLE_COMPLEX) :: P_TERM, Q_TERM, JW_MU, J_OVER_WE, JJ
        INTEGER(C_INT) :: ADJ_C, IRC

        JJ = (0.0_C_DOUBLE, 1.0_C_DOUBLE)
        DIRECTFN_CALL_COUNT = DIRECTFN_CALL_COUNT + 1
        Z_PAIR = (0.0, 0.0)

        IDO = MESH%TRIANGLES(TRI_OBS)%VERTEX_3D
        IDS = MESH%TRIANGLES(TRI_SRC)%VERTEX_3D

        ! ---- RWG 对顶点（全局节点号）、系数幅值 L/2A、正/负号 ----
        IF (TRI_OBS == RWG_OBS%POS_TRI_ID) THEN
            OPP_O  = RWG_OBS%POS_OPP_VERTEX
            COEF_O = REAL(RWG_OBS%POS_COEF, C_DOUBLE)
            S_O    = 1.0_C_DOUBLE
        ELSE
            OPP_O  = RWG_OBS%NEG_OPP_VERTEX
            COEF_O = REAL(RWG_OBS%NEG_COEF, C_DOUBLE)
            S_O    = -1.0_C_DOUBLE
        END IF
        IF (TRI_SRC == RWG_SRC%POS_TRI_ID) THEN
            OPP_S  = RWG_SRC%POS_OPP_VERTEX
            COEF_S = REAL(RWG_SRC%POS_COEF, C_DOUBLE)
            S_S    = 1.0_C_DOUBLE
        ELSE
            OPP_S  = RWG_SRC%NEG_OPP_VERTEX
            COEF_S = REAL(RWG_SRC%NEG_COEF, C_DOUBLE)
            S_S    = -1.0_C_DOUBLE
        END IF

        ! ---- 调用 DIRECTFN 的标量参数（K_D 在命中路径也要用，须在缓存查询前赋值） ----
        K_D = REAL(K, C_DOUBLE)

        ! ---- tri-pair 缓存查询：同一 (观察,源) 三角形对只积分一次 ----
        CACHE_HIT = .FALSE.
        IF (ALLOCATED(DFN_HASH)) THEN
            KEY = (TRI_OBS - 1) * (DFN_HASH_NT + 1) + TRI_SRC
            CIDX = DFN_HASH(KEY)
            IF (CIDX > 0) THEN
                CACHE_HIT = .TRUE.
                ADJ_I         = DFN_CACHE(CIDX)%ADJ
                DFN_OBS_NODES = DFN_CACHE(CIDX)%DFN_OBS_NODES
                DFN_SRC_NODES = DFN_CACHE(CIDX)%DFN_SRC_NODES
                ISS_RWG       = DFN_CACHE(CIDX)%ISS_RWG
                ISS_CONST(1)  = DFN_CACHE(CIDX)%ISS_CONST
                DFN_STAT_HIT  = DFN_STAT_HIT + 1
            END IF
        END IF

        IF (.NOT. CACHE_HIT) THEN
        ! ---- 统计公共顶点 ----
        ! E1 = 第一个公共顶点；N_COMMON==2 时 E2 = 第二个（公共边两端）；
        ! N_COMMON==1 时 E1 即唯一公共顶点
        N_COMMON = 0
        E1 = 0; E2 = 0
        DO I = 1, 3
            DO J = 1, 3
                IF (IDO(I) == IDS(J)) THEN
                    N_COMMON = N_COMMON + 1
                    IF (N_COMMON == 1) THEN
                        E1 = IDO(I)
                    ELSE IF (N_COMMON == 2) THEN
                        E2 = IDO(I)
                    END IF
                END IF
            END DO
        END DO

        ! ---- 按 DIRECTFN 约定重排 contour 顶点，并记录对顶点的局部编号 ----
        IF (N_COMMON == 3) THEN
            ! ST：同一三角形， contour = 三角形自身
            ADJ_I = 3
            DFN_OBS_NODES = IDO
            DFN_SRC_NODES = IDS
            NPTS = 3
        ELSE IF (N_COMMON == 2) THEN
            ! EA：contour = (公共边e1, 公共边e2, 观察第三顶点, 源第三顶点)
            ! DIRECTFN 内部：源 rp=(e1,e2,r3)，观察 rq=(e2,e1,r4)（set_4_pts_）
            ! Iss 分量按 (源槽, 观察槽) 行优先存放，取值时须转置索引（见下）
            ADJ_I = 2
            OBS3 = 0; SRC3 = 0
            DO I = 1, 3
                IF (IDO(I) /= E1 .AND. IDO(I) /= E2) OBS3 = IDO(I)
                IF (IDS(I) /= E1 .AND. IDS(I) /= E2) SRC3 = IDS(I)
            END DO
            DFN_OBS_NODES = [E1, E2, OBS3]
            DFN_SRC_NODES = [E2, E1, SRC3]
            NPTS = 4
        ELSE IF (N_COMMON == 1) THEN
            ! VA：contour = (公共顶点, 观察其余两点, 源其余两点)
            ! DIRECTFN 内部：源 rp=(p,r2,r3)，观察 rq=(p,r4,r5)（set_5_pts_）
            ! Iss 分量按 (源槽, 观察槽) 行优先存放，取值时须转置索引（见下）
            ADJ_I = 1
            J = 1
            DO I = 1, 3
                IF (IDO(I) /= E1) THEN
                    J = J + 1
                    DFN_OBS_NODES(J) = IDO(I)
                END IF
            END DO
            DFN_OBS_NODES(1) = E1
            J = 1
            DO I = 1, 3
                IF (IDS(I) /= E1) THEN
                    J = J + 1
                    DFN_SRC_NODES(J) = IDS(I)
                END IF
            END DO
            DFN_SRC_NODES(1) = E1
            NPTS = 5
        ELSE
            PRINT *, "DIRECTFN 警告：非奇异三角形对不应调用本例程，Z=0"
            RETURN
        END IF
        ADJ_C = INT(ADJ_I, C_INT)

        ! ---- contour 坐标（双精度） ----
        DO I = 1, 3
            PTS(:, I) = NODE_COORD_DP(MESH, DFN_OBS_NODES(I))
        END DO
        IF (N_COMMON == 2) THEN
            PTS(:, 4) = NODE_COORD_DP(MESH, DFN_SRC_NODES(3))
        ELSE IF (N_COMMON == 1) THEN
            PTS(:, 4) = NODE_COORD_DP(MESH, DFN_SRC_NODES(2))
            PTS(:, 5) = NODE_COORD_DP(MESH, DFN_SRC_NODES(3))
        END IF

        ! ---- 调用 DIRECTFN ----
        SELECT CASE (ADJ_I)
        CASE (3)
            N_GAUSS_DFN = DIRECTFN_N_GAUSS_ST
        CASE (2)
            N_GAUSS_DFN = DIRECTFN_N_GAUSS_EA
        CASE DEFAULT
            N_GAUSS_DFN = DIRECTFN_N_GAUSS_VA
        END SELECT
        IRC = C_DFN_TRI_ISS(ADJ_C, PTS, K_D, INT(N_GAUSS_DFN, C_INT), ISS_RWG, ISS_CONST)
        IF (IRC /= 0_C_INT) THEN
            PRINT *, "DIRECTFN 警告：调用失败，Z=0"
            RETURN
        END IF

        ! ---- 写入缓存 ----
        IF (ALLOCATED(DFN_HASH)) THEN
            IF (DFN_CACHE_N == DFN_CACHE_CAP) THEN
                BLOCK
                    TYPE(DFN_CACHE_ENTRY), ALLOCATABLE :: TMP(:)
                    ALLOCATE(TMP(2 * DFN_CACHE_CAP))
                    TMP(1:DFN_CACHE_CAP) = DFN_CACHE
                    CALL MOVE_ALLOC(TMP, DFN_CACHE)
                    DFN_CACHE_CAP = 2 * DFN_CACHE_CAP
                END BLOCK
            END IF
            DFN_CACHE_N = DFN_CACHE_N + 1
            DFN_CACHE(DFN_CACHE_N)%TRI_O = TRI_OBS
            DFN_CACHE(DFN_CACHE_N)%TRI_S = TRI_SRC
            DFN_CACHE(DFN_CACHE_N)%ADJ = ADJ_I
            DFN_CACHE(DFN_CACHE_N)%DFN_OBS_NODES = DFN_OBS_NODES
            DFN_CACHE(DFN_CACHE_N)%DFN_SRC_NODES = DFN_SRC_NODES
            DFN_CACHE(DFN_CACHE_N)%ISS_RWG = ISS_RWG
            DFN_CACHE(DFN_CACHE_N)%ISS_CONST = ISS_CONST(1)
            DFN_HASH(KEY) = DFN_CACHE_N
            DFN_STAT_MISS = DFN_STAT_MISS + 1
            SELECT CASE (ADJ_I)
            CASE (3); DFN_STAT_ST = DFN_STAT_ST + 1
            CASE (2); DFN_STAT_EA = DFN_STAT_EA + 1
            CASE (1); DFN_STAT_VA = DFN_STAT_VA + 1
            END SELECT
        END IF
        END IF  ! .NOT. CACHE_HIT

        ! DIRECTFN 局部编号（1..3）中的观察/源对顶点位置
        IO = 0; JO = 0
        DO I = 1, 3
            IF (DFN_OBS_NODES(I) == OPP_O) IO = I
            IF (DFN_SRC_NODES(I) == OPP_S) JO = I
        END DO
        IF (IO == 0 .OR. JO == 0) THEN
            PRINT *, "DIRECTFN 警告：RWG 对顶点不在对应三角形上，Z=0"
            RETURN
        END IF

        ! ---- 调试输出（供 Python 独立参考程序解析） ----
        IF (DIRECTFN_DEBUG > 0) THEN
            PRINT '(A, 3(1X, I0))', "@@DFNRAW1", ADJ_I, IO, JO
            PRINT '(A, 6I6)', "@@DFNNODES", DFN_OBS_NODES, DFN_SRC_NODES
            PRINT '(A, 4ES24.15)', "@@DFNRAW2", COEF_O, COEF_S, S_O, S_S
            PRINT '(A, 18ES24.15)', "@@DFNRAW3", &
                (REAL(ISS_RWG(I), C_DOUBLE), AIMAG(ISS_RWG(I)), I = 1, 9)
            PRINT '(A, 2ES24.15)', "@@DFNRAW4", &
                REAL(ISS_CONST(1), C_DOUBLE), AIMAG(ISS_CONST(1))
        END IF

        ! ---- 组合成阻抗贡献（与 Z_MATRIX 远场分支同一公式） ----
        ! iss_rwg(io,jo) = c_o·c_s·∬(r-v_io)·(r'-v_jo) G dS dS'
        !                = c_o·c_s·P_TERM
        ! 注意：DIRECTFN 的 Iss 线性布局为 (源槽, 观察槽) 行优先（见
        ! directfn_kernel_tri.cpp 的 pmf_rwg_val_ 顺序），即分量 (a,b) 位于
        ! 槽 a+3*(b-1)，a=源自由顶点槽位、b=观察自由顶点槽位。
        ! 因此取 (观察=IO, 源=JO) 分量须读 JO + 3*(IO-1)（已用独立参考验证至 1e-11）。
        P_TERM = ISS_RWG(JO + 3 * (IO - 1)) / (COEF_O * COEF_S)
        Q_TERM = -4.0_C_DOUBLE_COMPLEX * ISS_CONST(1)

        JW_MU     = JJ * REAL(ETA0, C_DOUBLE) * K_D
        J_OVER_WE = JJ * REAL(ETA0, C_DOUBLE) / K_D

        Z_PAIR = CMPLX(S_O * COEF_O * S_S * COEF_S * &
                       (JW_MU * P_TERM + J_OVER_WE * Q_TERM))
    END SUBROUTINE CALC_EFIE_DIRECTFN_PAIR

    ! 取节点坐标（双精度）
    FUNCTION NODE_COORD_DP(MESH, NID) RESULT(R)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        INTEGER, INTENT(IN) :: NID
        REAL(C_DOUBLE) :: R(3)
        R(1) = REAL(MESH%NODES(NID)%X, C_DOUBLE)
        R(2) = REAL(MESH%NODES(NID)%Y, C_DOUBLE)
        R(3) = REAL(MESH%NODES(NID)%Z, C_DOUBLE)
    END FUNCTION NODE_COORD_DP

END MODULE DIRECTFN_BRIDGE
