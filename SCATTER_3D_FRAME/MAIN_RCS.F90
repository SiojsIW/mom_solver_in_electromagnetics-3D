! 金属球 RCS 计算主程序
! 读入三角网格 → 构建 MoM 系统 → 求解电流 → 计算双站 RCS → 输出文件
PROGRAM MAIN_RCS
    USE EM_TYPES
    USE MESH_GEOMETRY
    USE NUMERICAL_INTEGRATION
    USE GREEN_FUNCTIONS
    USE RWG_BASIS_BUILD
    USE Z_MATRIX
    USE RHS
    USE FAR_FIELD
    USE QBX_EFIE
    USE DIRECTFN_BRIDGE
    IMPLICIT NONE

    ! ---- 可调参数 ----
    REAL, PARAMETER :: LAMBDA = 2.094395          ! 波长 (m), ka = 3.0
    REAL, PARAMETER :: ETA0 = 120.0 * PI     ! 自由空间波阻抗
    REAL, PARAMETER :: E0 = 1.0              ! 入射电场幅度 (V/m)
    REAL, PARAMETER :: ALPHA_CFIE = 1.0      ! 纯 EFIE：α=1.0
    INTEGER, PARAMETER :: N_ANGLES = 181     ! 观察角采样数（0°~180°，步长 1°）
    INTEGER, PARAMETER :: NEAR_METHOD = 3     ! 1=奇异提取法  2=QBX  3=DIRECTFN

    ! ---- 入射波参数 ----
    ! 平面波：E_inc = E0 * E_POL * exp(-j*k * K_HAT·r)
    ! 波从 +z 方向入射：K_HAT = (0, 0, 1)
    ! x-极化：E_POL = (1, 0, 0)
    REAL :: K_HAT(3) = [0.0, 0.0, 1.0]
    REAL :: E_POL(3) = [1.0, 0.0, 0.0]

    ! ---- 变量 ----
    TYPE(MESH_3D) :: MESH
    TYPE(GAUSS_TRI_DATA) :: GDATA, GDATA_FAR
    REAL :: K
    INTEGER :: IOS, M, N, I
    COMPLEX, ALLOCATABLE :: Z(:, :), V(:), I_CURRENT(:)
    REAL, ALLOCATABLE :: OBS_ANGLES(:), RCS_THETA(:), RCS_PHI(:)
    COMPLEX :: Z_MN_E, Z_MN_M, V_M_E, V_M_M
    ! (timing not used yet)

    ! ---- 文件路径 ----
    CHARACTER(LEN=256) :: NODE_FILE, TRI_FILE
    CHARACTER(LEN=*), PARAMETER :: OUT_FILE = "rcs_results.txt"
    INTEGER :: NARGS

    ! 命令行参数：用法 rcs_solver.exe [node_file] [tri_file]
    NARGS = COMMAND_ARGUMENT_COUNT()
    IF (NARGS >= 2) THEN
        CALL GET_COMMAND_ARGUMENT(1, NODE_FILE)
        CALL GET_COMMAND_ARGUMENT(2, TRI_FILE)
    ELSE
        NODE_FILE = "point.txt"
        TRI_FILE  = "triangel_point.txt"
    END IF

    K = 2.0 * PI / LAMBDA

    ! ==================== 1. 读入网格 ====================
    PRINT *, "=============================================="
    PRINT *, "  金属球 RCS 计算（CFIE / MoM）"
    PRINT *, "=============================================="
    PRINT '(A, F8.4, A)', " 波长 λ = ", LAMBDA, " m"
    PRINT '(A, F8.4, A)', " 波数 k = ", K, " rad/m"
    PRINT '(A, F8.4, A)', " 球半径 a ≈ 1.0 m, ka = ", K * 1.0

    CALL READ_MESH_FROM_FILES(MESH, NODE_FILE, TRI_FILE, IOS)
    IF (IOS /= 0) THEN
        PRINT *, "错误：网格读取失败，程序终止"
        STOP
    END IF

    ! ==================== 2. 构建几何 ====================
    CALL UPDATE_ALL_TRIANGLES_3D(MESH)
    CALL VERIFY_OUTWARD_NORMALS(MESH)
    CALL EXTRACT_EDGES_FROM_MESH(MESH)
    CALL UPDATE_EDGE_GEOMETRY(MESH)
    CALL BUILD_RWG_BASIS(MESH)

    PRINT '(A, I0)', " 节点数:       ", MESH%NUM_NODE
    PRINT '(A, I0)', " 三角形数:     ", MESH%NUM_TRIANGLE
    PRINT '(A, I0)', " 边数:         ", MESH%NUM_EDGE
    PRINT '(A, I0)', " RWG 基函数数: ", MESH%RWG_NUM

    ! ---- 选择 EFIE 近场积分方法 ----
    SELECT CASE (NEAR_METHOD)
    CASE (2)
        EFIE_NEAR_METHOD = EFIE_NEAR_QBX
        PRINT *, "EFIE 近场分支：QBX"
    CASE (3)
        EFIE_NEAR_METHOD = EFIE_NEAR_DIRECTFN
        PRINT *, "EFIE 近场分支：DIRECTFN"
    CASE DEFAULT
        EFIE_NEAR_METHOD = EFIE_NEAR_EXTRACTION
        PRINT *, "EFIE 近场分支：奇异提取法"
    END SELECT

    ! ==================== 3. 高斯积分数据 ====================
    CALL INIT_GAUSS_TRI(GAUSS_12PT, GDATA)
    CALL INIT_GAUSS_TRI(GAUSS_12PT, GDATA_FAR)
    PRINT '(A, I0, A, I0, A)', " 矩阵积分: ", GDATA%N_POINTS, " 点,  远场积分: ", GDATA_FAR%N_POINTS, " 点"

    ! ==================== 4. 组装阻抗矩阵 ====================
    PRINT *, "正在组装 CFIE 阻抗矩阵 ..."
    ALLOCATE(Z(MESH%RWG_NUM, MESH%RWG_NUM))
    Z = (0.0, 0.0)

    IF (NEAR_METHOD == 3) CALL DIRECTFN_RESET_CACHE(MESH%NUM_TRIANGLE)

    DO M = 1, MESH%RWG_NUM
        IF (MOD(M, 20) == 0) PRINT '(A, I0, A, I0)', "  行 ", M, " / ", MESH%RWG_NUM
        DO N = 1, MESH%RWG_NUM
            CALL CALC_EFIE_MATRIX_ELEMENT(MESH, MESH%RWG_BASES(M), &
                        MESH%RWG_BASES(N), GDATA, K, ETA0, Z_MN_E)
            CALL CALC_MFIE_MATRIX_ELEMENT(MESH, MESH%RWG_BASES(M), &
                        MESH%RWG_BASES(N), GDATA, K, Z_MN_M)
            Z(M, N) = ALPHA_CFIE * Z_MN_E + (1.0 - ALPHA_CFIE) * ETA0 * Z_MN_M
        END DO
    END DO
    PRINT *, "阻抗矩阵组装完成"

    ! ==================== 5. 组装右端向量 ====================
    PRINT *, "正在组装 CFIE 右端向量 ..."
    ALLOCATE(V(MESH%RWG_NUM))
    V = (0.0, 0.0)

    DO M = 1, MESH%RWG_NUM
        CALL CALC_EFIE_RHS(MESH, MESH%RWG_BASES(M), GDATA, K, K_HAT, E_POL, E0, ETA0, V_M_E)
        CALL CALC_MFIE_RHS(MESH, MESH%RWG_BASES(M), GDATA, K, K_HAT, E_POL, E0, ETA0, V_M_M)
        V(M) = ALPHA_CFIE * V_M_E + (1.0 - ALPHA_CFIE) * ETA0 * V_M_M
    END DO
    PRINT *, "右端向量组装完成"

    ! ==================== 6. 求解 Z * I = V ====================
    PRINT *, "正在求解线性系统 (N=", MESH%RWG_NUM, ") ..."
    ALLOCATE(I_CURRENT(MESH%RWG_NUM))
    CALL SOLVE_COMPLEX_LINEAR_SYSTEM(MESH%RWG_NUM, Z, V, I_CURRENT)
    PRINT *, "线性系统求解完成"

    ! 打印部分电流系数供验证
    PRINT *, "电流系数幅值（前 5 个）:"
    DO M = 1, MIN(5, MESH%RWG_NUM)
        PRINT '(A, I0, A, E14.6)', "  |I(", M, ")| = ", ABS(I_CURRENT(M))
    END DO

    ! ==================== 7. 计算双站 RCS ====================
    PRINT *, "正在计算双站 RCS (", N_ANGLES, " 个角度) ..."

    ALLOCATE(OBS_ANGLES(N_ANGLES), RCS_THETA(N_ANGLES), RCS_PHI(N_ANGLES))
    CALL CALC_BISTATIC_RCS(MESH, I_CURRENT, MESH%RWG_NUM, K, ETA0, E0, &
                            N_ANGLES, OBS_ANGLES, RCS_THETA, RCS_PHI, GDATA_FAR)
    PRINT *, "RCS 计算完成"

    ! 近场积分方法调用统计
    IF (NEAR_METHOD == 2) THEN
        PRINT '(A, I0)', " QBX 调用次数:          ", QBX_CALL_COUNT
        PRINT '(A, I0)', " QBX 许可性检查失败次数:", QBX_ADM_FAIL_COUNT
    ELSE IF (NEAR_METHOD == 3) THEN
        PRINT '(A, I0)', " DIRECTFN 调用次数(RWG 对): ", DIRECTFN_CALL_COUNT
        CALL DIRECTFN_PRINT_STATS()
    END IF

    ! ==================== 8. 输出结果 ====================
    OPEN(UNIT=20, FILE=OUT_FILE, STATUS='REPLACE', ACTION='WRITE')
    WRITE(20, '(A)') "# 金属球双站 RCS（CFIE MoM）"
    WRITE(20, '(A, F8.4)') "# 波长 (m): ", LAMBDA
    WRITE(20, '(A, F8.4)') "# 波数 k (rad/m): ", K
    WRITE(20, '(A, I0)') "# RWG 基函数数: ", MESH%RWG_NUM
    WRITE(20, '(A)') "# 入射方向: +z, 极化: x, 观察面: φ=0 (xz-plane)"
    WRITE(20, '(A)') "# theta(deg)  RCS_theta(dBsm)  RCS_phi(dBsm)  RCS_total(dBsm)"

    DO I = 1, N_ANGLES
        WRITE(20, '(F8.2, 3F16.6)') OBS_ANGLES(I), RCS_THETA(I), RCS_PHI(I), &
            10.0 * LOG10(10.0**(RCS_THETA(I)/10.0) + 10.0**(RCS_PHI(I)/10.0))
    END DO

    CLOSE(20)
    PRINT *, "结果已写入: ", OUT_FILE

    ! ==================== 9. 清理 ====================
    DEALLOCATE(Z, V, I_CURRENT, OBS_ANGLES, RCS_THETA, RCS_PHI)

    PRINT *, "=============================================="
    PRINT *, "  计算完成"
    PRINT *, "=============================================="

CONTAINS

    ! 复数线性方程组求解（Gauss 消元 + 部分选主元）
    SUBROUTINE SOLVE_COMPLEX_LINEAR_SYSTEM(N, A, B, X)
        INTEGER, INTENT(IN) :: N
        COMPLEX, INTENT(INOUT) :: A(N, N), B(N)
        COMPLEX, INTENT(OUT) :: X(N)

        INTEGER :: I, J, K, IP
        REAL :: MAX_VAL
        COMPLEX :: TEMP

        ! ---- LU 分解（部分选主元） ----
        DO K = 1, N - 1
            ! 选主元
            IP = K
            MAX_VAL = ABS(A(K, K))
            DO I = K + 1, N
                IF (ABS(A(I, K)) > MAX_VAL) THEN
                    MAX_VAL = ABS(A(I, K))
                    IP = I
                END IF
            END DO

            IF (MAX_VAL < 1.0E-30) THEN
                PRINT *, "错误：矩阵奇异，列 ", K, " 主元 ≈ 0"
                X = (0.0, 0.0)
                RETURN
            END IF

            ! 交换行
            IF (IP /= K) THEN
                DO J = K, N
                    TEMP = A(K, J)
                    A(K, J) = A(IP, J)
                    A(IP, J) = TEMP
                END DO
                TEMP = B(K)
                B(K) = B(IP)
                B(IP) = TEMP
            END IF

            ! 消元
            DO I = K + 1, N
                A(I, K) = A(I, K) / A(K, K)
                DO J = K + 1, N
                    A(I, J) = A(I, J) - A(I, K) * A(K, J)
                END DO
                B(I) = B(I) - A(I, K) * B(K)
            END DO
        END DO

        ! ---- 回代 ----
        DO I = N, 1, -1
            X(I) = B(I)
            DO J = I + 1, N
                X(I) = X(I) - A(I, J) * X(J)
            END DO
            X(I) = X(I) / A(I, I)
        END DO
    END SUBROUTINE SOLVE_COMPLEX_LINEAR_SYSTEM

END PROGRAM MAIN_RCS
