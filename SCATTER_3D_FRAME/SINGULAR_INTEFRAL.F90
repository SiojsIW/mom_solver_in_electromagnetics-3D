MODULE SINGULAR_INTEGRAL

    CONTAINS
    ! 对固定场点r，计算源三角形T上的标量势闭合积分
    SUBROUTINE ANALYTIC_SCALAR_POT_1_OVER_R(R_PT, V1, V2, V3, RESULT)
        REAL, INTENT(IN) :: R_PT(3) ! 场点坐标
        REAL, INTENT(IN) :: V1(3), V2(3), V3(3) ! 三角形顶点坐标
        REAL, INTENT(OUT) :: RESULT

        REAL :: NORMAL(3) ! 三角形法向量
        INTEGER :: I ! 循环变量
        REAL :: V_START(3), V_END(3) ! 边的起点和终点
        REAL :: EDGE_VEC(3) ! 边向量
        REAL :: L_I ! 边长度
        REAL :: T_HAT(3) ! 边的单位切向量
        REAL :: M_HAT(3) ! 边的单位法向量
        REAL :: H_I ! 场点到边的垂直距离(有符号距离)
        REAL :: S_START, S_END ! 场点到边的起点和终点的投影长度(有符号)
        REAL :: R_START, R_END ! 场点到边的起点和终点的距离
        REAL :: CONTRIB     ! 当前边的积分贡献
        REAL :: V_OPP(3) ! 当前边的对边顶点坐标 
        REAL :: DOT_CHECK ! 验证内法向是否指向内部的点积

        NORMAL(1) = (V2(2) - V1(2)) * (V3(3) - V1(3)) - (V2(3) - V1(3)) * (V3(2) - V1(2))
        NORMAL(2) = (V2(3) - V1(3)) * (V3(1) - V1(1)) - (V2(1) - V1(1)) * (V3(3) - V1(3))
        NORMAL(3) = (V2(1) - V1(1)) * (V3(2) - V1(2)) - (V2(2) - V1(2)) * (V3(1) - V1(1))
        NORMAL = NORMAL / SQRT(SUM(NORMAL ** 2)) ! 把法向量归一化
        RESULT = 0.0

        DO I = 1, 3
            SELECT CASE (I)
                CASE (1)
                    V_START = V1
                    V_END = V2
                    V_OPP = V3
                CASE (2)
                    V_START = V2
                    V_END = V3
                    V_OPP = V1
                CASE (3)
                    V_START = V3
                    V_END = V1
                    V_OPP = V2
            END SELECT

            EDGE_VEC = V_END - V_START
            L_I = SQRT(SUM(EDGE_VEC ** 2))
            T_HAT = EDGE_VEC / L_I ! 计算边的单位切向量

            ! 计算边的内法向 M_HAT = NORMAL x T_HAT
            M_HAT(1) = NORMAL(2) * T_HAT(3) - NORMAL(3) * T_HAT(2)
            M_HAT(2) = NORMAL(3) * T_HAT(1) - NORMAL(1) * T_HAT(3)
            M_HAT(3) = NORMAL(1) * T_HAT(2) - NORMAL(2) * T_HAT(1)

            ! 验证内法向是否指向内部
            DOT_CHECK = DOT_PRODUCT(M_HAT, V_OPP - V_START)
            IF (DOT_CHECK < 0.0) THEN
                M_HAT = -M_HAT
            END IF

            ! 计算场点到边的垂直距离 H_I
            H_I = DOT_PRODUCT(R_PT - V_START, M_HAT)
            S_START = DOT_PRODUCT(R_PT - V_START, T_HAT)
            S_END = DOT_PRODUCT(R_PT - V_END, T_HAT)
            R_START = SQRT(SUM((R_PT - V_START) ** 2))
            R_END = SQRT(SUM((R_PT - V_END) ** 2))

            ! 计算当前边的积分贡献
            IF (R_START + S_START > 1.0E-12) THEN
                CONTRIB = H_I * LOG((R_END + S_END) / (R_START + S_START))
            ELSE ! 当场点非常接近边的起点时，避免对数奇异性
                CONTRIB = H_I * (LOG(R_START - S_START) / (R_END - S_END)) ! 使用稳定替代形式
            END IF

            RESULT = RESULT + CONTRIB
        END DO

    END SUBROUTINE ANALYTIC_SCALAR_POT_1_OVER_R

    ! 计算三个线性势
    SUBROUTINE ANALYTIC_LINEAR_POT_1_OVER_R(R_PT, V1, V2, V3, RESULT_3)
        REAL, INTENT(IN) :: R_PT(3) ! 场点坐标
        REAL, INTENT(IN) :: V1(3), V2(3), V3(3) ! 三角形顶点坐标
        REAL, INTENT(OUT) :: RESULT_3(3)

        REAL :: NORMAL(3) ! 三角形法向量
        INTEGER :: I ! 循环变量
        REAL :: V_START(3), V_END(3) ! 边的起点和终点
        REAL :: EDGE_VEC(3) ! 边向量
        REAL :: L_I ! 边长度
        REAL :: T_HAT(3) ! 边的单位切向量
        REAL :: M_HAT(3) ! 边的单位法向量
        REAL :: H_I ! 场点到边的垂直距离(有符号距离)
        REAL :: S_START, S_END ! 场点到边的起点和终点的投影长度(有符号)
        REAL :: R_START, R_END ! 场点到边的起点和终点的距离
        REAL :: V_OPP(3) ! 当前边的对边顶点坐标 
        REAL :: DOT_CHECK ! 验证内法向是否指向内部的点积
        REAL :: LOG_TERM ! 对数项
        REAL :: SR_TERM ! 投影-距离项
        REAL :: C_START, C_END ! 当前边对两个端点的贡献

        NORMAL(1) = (V2(2) - V1(2)) * (V3(3) - V1(3)) - (V2(3) - V1(3)) * (V3(2) - V1(2))
        NORMAL(2) = (V2(3) - V1(3)) * (V3(1) - V1(1)) - (V2(1) - V1(1)) * (V3(3) - V1(3))
        NORMAL(3) = (V2(1) - V1(1)) * (V3(2) - V1(2)) - (V2(2) - V1(2)) * (V3(1) - V1(1))
        NORMAL = NORMAL / SQRT(SUM(NORMAL ** 2)) ! 把法向量归一化
        RESULT_3 = 0.0

        DO I = 1, 3
            SELECT CASE (I)
                CASE (1)
                    V_START = V1
                    V_END = V2
                    V_OPP = V3
                CASE (2)
                    V_START = V2
                    V_END = V3
                    V_OPP = V1
                CASE (3)
                    V_START = V3
                    V_END = V1
                    V_OPP = V2
            END SELECT

            EDGE_VEC = V_END - V_START
            L_I = SQRT(SUM(EDGE_VEC ** 2))
            T_HAT = EDGE_VEC / L_I ! 计算边的单位切向量

            ! 计算边的内法向 M_HAT = NORMAL x T_HAT
            M_HAT(1) = NORMAL(2) * T_HAT(3) - NORMAL(3) * T_HAT(2)
            M_HAT(2) = NORMAL(3) * T_HAT(1) - NORMAL(1) * T_HAT(3)
            M_HAT(3) = NORMAL(1) * T_HAT(2) - NORMAL(2) * T_HAT(1)

            ! 验证内法向是否指向内部
            DOT_CHECK = DOT_PRODUCT(M_HAT, V_OPP - V_START)
            IF (DOT_CHECK < 0.0) THEN
                M_HAT = -M_HAT
            END IF

            ! 计算场点到边的垂直距离 H_I
            H_I = DOT_PRODUCT(R_PT - V_START, M_HAT)
            S_START = DOT_PRODUCT(R_PT - V_START, T_HAT)
            S_END = DOT_PRODUCT(R_PT - V_END, T_HAT)
            R_START = SQRT(SUM((R_PT - V_START) ** 2))
            R_END = SQRT(SUM((R_PT - V_END) ** 2))

            ! 计算当前边的积分贡献
            IF (R_START + S_START > 1.0E-12) THEN
                LOG_TERM = LOG((R_END + S_END) / (R_START + S_START))
            ELSE ! 当场点非常接近边的起点时，避免对数奇异性
                LOG_TERM = LOG(R_START - S_START) - LOG(R_END - S_END) ! 使用稳定替代形式
            END IF

            ! 计算投影-距离项
            SR_TERM = S_END * R_END - S_START * R_START

            ! 计算当前边对两个端点的贡献
            C_START = 0.5 * H_I * LOG_TERM - (H_I**2 / (2 * L_I)) * LOG_TERM - SR_TERM / (2 * L_I)
            C_END = (H_I**2 / (2 * L_I)) * LOG_TERM + SR_TERM / (2 * L_I)

            ! 将贡献累加到结果中
            IDX_END = MOD(I, 3) + 1
            RESULT_3(I) = RESULT_3(I) + C_START
            RESULT_3(IDX_END) = RESULT_3(IDX_END) + C_END

        END DO

    END SUBROUTINE ANALYTIC_LINEAR_POT_1_OVER_R

    
END MODULE SINGULAR_INTEGRAL