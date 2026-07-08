MODULE GREEN_FUNCTIONS
    USE EM_TYPES
    USE NUMERICAL_INTEGRATION
    IMPLICIT NONE
CONTAINS

    ! ¼ÆËã¸ñÁÖº¯Êý£¬ÊäÈëÁ½¸ö×ø±ê£¬¼ÆËãÁ½¸ö×ø±êµãµÄ¸ñÁÖº¯Êý
    COMPLEX FUNCTION GREEN_FUNC(R, RP, K) RESULT(G)
        REAL, INTENT(IN):: R(3) ! Ô´µã×ø±ê
        REAL, INTENT(IN):: RP(3) ! ³¡µã×ø±ê
        REAL,INTENT(IN) :: K  ! ²¨Êý
        REAL :: DIST
        
        ! COMPLEX ÀàÐÍ£¬ÐéÊýµ¥Î»ÊÇ (0.0, 1.0)¡£CEXP ÊÇ¸´ÊýÖ¸Êý¡£PI = 4.0*ATAN(1.0)¡£
        DIST = SQRT(SUM((R - RP) ** 2)) ! ³¡µãÓëÔ´µã¾àÀë

        IF (DIST < 1.0E-10) THEN
            G = (0.0, 0.0)
        ELSE
            G = CEXP((0.0, -1.0) * K * DIST) / (4.0 * PI * DIST)
        END IF
    END FUNCTION

    ! ¼ÆËã×è¿¹ÔªËØZ_MN£¬Ö»ÊÇ¼ÆËã¸ñÁÖº¯ÊýÔÚÁ½¸öÈý½ÇÐÎÉÏµÄ£¬Ã»ÓÐf¡£
    SUBROUTINE CALC_GREEN_MATRIX_ELEMENT(MESH, TM_ID, TN_ID, GDATA, K, Z_MN)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        INTEGER :: TM_ID ! ³¡Èý½ÇÐÎ±àºÅ
        INTEGER :: TN_ID ! Ô´Èý½ÇÐÎ±àºÅ
        INTEGER :: I, J
        REAL :: K  ! ²¨Êý
        REAL, ALLOCATABLE :: GLOBAL_PTS_M(:, :) ! ³¡Èý½ÇÐÎ¸ßË¹»ý·ÖµãµÄÈ«¾Ö×ø±ê
        REAL, ALLOCATABLE :: GLOBAL_PTS_N(:, :) ! Ô´Èý½ÇÐÎ¸ßË¹»ý·ÖµãµÄÈ«¾Ö×ø±ê
        COMPLEX :: G
        COMPLEX, INTENT(OUT) :: Z_MN 

        ALLOCATE(GLOBAL_PTS_M(3, GDATA%N_POINTS)) 
        ALLOCATE(GLOBAL_PTS_N(3, GDATA%N_POINTS)) 
        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TM_ID, GDATA, GLOBAL_PTS_M)  ! µÃµ½³¡Èý½ÇÐÎµÄ¸ßË¹»ý·ÖµãµÄÈ«¾Ö×ø±ê
        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TN_ID, GDATA, GLOBAL_PTS_N)  ! µÃµ½Ô´Èý½ÇÐÎµÄ¸ßË¹»ý·ÖµãµÄÈ«¾Ö×ø±ê

        Z_MN = (0.0, 0.0)
        DO I = 1, GDATA%N_POINTS ! ³¡
            DO J = 1, GDATA%N_POINTS ! Ô´
                G = GREEN_FUNC(GLOBAL_PTS_M(:, I), GLOBAL_PTS_N(:, J), K) 
                Z_MN = Z_MN + GDATA%WEIGHTS(I) * GDATA%WEIGHTS(J) * G
            END DO
        END DO

        DEALLOCATE(GLOBAL_PTS_M) 
        DEALLOCATE(GLOBAL_PTS_N) 

        Z_MN = Z_MN * MESH%TRIANGLES(TM_ID)%AREA * MESH%TRIANGLES(TN_ID)%AREA

    END SUBROUTINE

    ! ¼ÆËãÒ»¸öÈý½ÇÐÎ¶ÔÉÏµÄËÄ¸ö»ý·ÖI1,I2,I3,I4£¬£¨¼ÆËãEFIE·Ç¶Ô½Ç¾ØÕóÔªËØÓÃ£©
    SUBROUTINE CALC_GREEN_INTEGALS(MESH, TRI_A, TRI_B, GDATA, K, I1, I2, I3, I4)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        INTEGER, INTENT(IN) :: TRI_A ! Ô´Èý½ÇÐÎ
        INTEGER, INTENT(IN) :: TRI_B ! ³¡Èý½ÇÐÎ
        REAL, INTENT(IN) :: K  ! ²¨Êý
        COMPLEX, INTENT(OUT) :: I1 ! ±êÁ¿
        COMPLEX, INTENT(OUT) :: I4 ! ±êÁ¿
        COMPLEX, INTENT(OUT) :: I2(3) ! Ê¸Á¿
        COMPLEX, INTENT(OUT) :: I3(3) ! Ê¸Á¿
    
        REAL :: AREA_A, AREA_B
        REAL :: GLOBAL_PTS_A(3, GDATA%N_POINTS), GLOBAL_PTS_B(3, GDATA%N_POINTS) ! Èý½ÇÐÎ¸ßË¹»ý·ÖµãµÄÈ«¾Ö×ø±ê
        INTEGER :: I, J
        COMPLEX :: G
        I1 = (0.0, 0.0)
        I4 = (0.0, 0.0)
        I2 = (0.0, 0.0)   ! ±êÁ¿¸´Êý×Ô¶¯¹ã²¥µ½3¸öÔªËØ
        I3 = (0.0, 0.0)

        AREA_A = MESH%TRIANGLES(TRI_A)%AREA
        AREA_B = MESH%TRIANGLES(TRI_B)%AREA

        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_A, GDATA, GLOBAL_PTS_A)
        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_B, GDATA, GLOBAL_PTS_B)

        DO I = 1, GDATA%N_POINTS ! i´ú±íÔ´µã
            DO J = 1, GDATA%N_POINTS
                G = GREEN_FUNC(GLOBAL_PTS_A(:, I), GLOBAL_PTS_B(:, J), K)

                I1 = I1 + GDATA%WEIGHTS(I) * GDATA%WEIGHTS(J) * G  ! ÏÈ²»³ËÃæ»ý£¬Ñ­»·½áÊøºóÍ³Ò»³ËÃæ»ý
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

    ! ¼ÆËãÒ»¸öÈý½ÇÐÎ¶ÔÉÏµÄËÄ¸ö»ý·ÖI1,I2,I3,I4£¬£¨¼ÆËãEFIE¶Ô½Ç¾ØÕóÔªËØÓÃ£¬G_SMOOTH£©
    SUBROUTINE CALC_GREEN_SMOOTH_INTEGALS(MESH, TRI_A, TRI_B, GDATA, K, I1, I2, I3, I4)
        TYPE(MESH_3D), INTENT(IN) :: MESH
        TYPE(GAUSS_TRI_DATA), INTENT(IN) :: GDATA
        INTEGER, INTENT(IN) :: TRI_A ! Ô´Èý½ÇÐÎ
        INTEGER, INTENT(IN) :: TRI_B ! ³¡Èý½ÇÐÎ
        REAL, INTENT(IN) :: K  ! ²¨Êý
        COMPLEX, INTENT(OUT) :: I1 ! ±êÁ¿
        COMPLEX, INTENT(OUT) :: I4 ! ±êÁ¿
        COMPLEX, INTENT(OUT) :: I2(3) ! Ê¸Á¿
        COMPLEX, INTENT(OUT) :: I3(3) ! Ê¸Á¿
    
        REAL :: AREA_A, AREA_B
        REAL :: GLOBAL_PTS_A(3, GDATA%N_POINTS), GLOBAL_PTS_B(3, GDATA%N_POINTS) ! Èý½ÇÐÎ¸ßË¹»ý·ÖµãµÄÈ«¾Ö×ø±ê
        INTEGER :: I, J
        COMPLEX :: G
        I1 = (0.0, 0.0)
        I4 = (0.0, 0.0)
        I2 = (0.0, 0.0)   ! ±êÁ¿¸´Êý×Ô¶¯¹ã²¥µ½3¸öÔªËØ
        I3 = (0.0, 0.0)

        AREA_A = MESH%TRIANGLES(TRI_A)%AREA
        AREA_B = MESH%TRIANGLES(TRI_B)%AREA

        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_A, GDATA, GLOBAL_PTS_A)
        CALL GET_TRI_GLOBAL_GAUSS_POINTS(MESH, TRI_B, GDATA, GLOBAL_PTS_B)

        DO I = 1, GDATA%N_POINTS ! i´ú±íÔ´µã
            DO J = 1, GDATA%N_POINTS
                G = GREEN_FUNC_SMOOTH(GLOBAL_PTS_A(:, I), GLOBAL_PTS_B(:, J), K)

                I1 = I1 + GDATA%WEIGHTS(I) * GDATA%WEIGHTS(J) * G  ! ÏÈ²»³ËÃæ»ý£¬Ñ­»·½áÊøºóÍ³Ò»³ËÃæ»ý
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

    ! ¼ÆËã±êÁ¿¸ñÁÖº¯ÊýµÄÌÝ¶È
    SUBROUTINE GARD_GREEN_FUNC(R, RP, K, GRAD_G)
        REAL, INTENT(IN) :: R(3) ! ³¡µã×ø±ê
        REAL, INTENT(IN) :: RP(3) ! Ô´µã×ø±ê
        REAL, INTENT(IN) :: K  ! ²¨Êý
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

<<<<<<< HEAD
=======
    ! ÌáÈ¡ 1/R ºóµÄÊ£ÓàºË£º(e^{-jkR} - 1) / (4¦ÐR)
    ! R=0 Ê±¼«ÏÞÎª -jk/(4¦Ð)£¬ÎÞÆæÒìÐÔ
    COMPLEX FUNCTION GREEN_FUNC_SMOOTH(R, RP, K) RESULT(G)
        REAL, INTENT(IN) :: R(3) ! ³¡µã×ø±ê
        REAL, INTENT(IN) :: RP(3) ! Ô´µã×ø±ê
        REAL, INTENT(IN) :: K  ! ²¨Êý
        REAL :: DIST
        COMPLEX :: PHASE

        DIST = SQRT(SUM((R - RP) ** 2)) ! ³¡µãÓëÔ´µã¾àÀë
        IF (DIST < 1.0E-10) THEN
            G = -(0.0, 1.0) * K / (4.0 * PI)
        ELSE
            PHASE = CEXP((0.0, -1.0) * K * DIST)
            G = (PHASE - (1.0, 0.0)) / (4.0 * PI * DIST)
        END IF
    END FUNCTION

>>>>>>> 631d457 (codeÂ£ÂºÂ¸Ã¼ÃÃ‚Ä£Â¿é£ºÂ¼Ã†Ã‹ãµ¥Â¸Ã¶ÃˆÃ½Â½Ã‡ÃÃŽTÂ¶Ã”Ã—Ã”Ã‰Ã­ÂµÃ„ÃÃªÃ•Ã» EFIE Ã—è¿¹Â¹Â±Ã×¡Â£)
END MODULE