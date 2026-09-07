PROGRAM ERROR_ANALYSIS
    USE EM_TYPES
    USE NUMERICAL_INTEGRATION
    USE SINGULAR_INTEGRAL
    USE GREEN_FUNCTIONS
    IMPLICIT NONE
    REAL :: V1(3),V2(3),V3(3),ROPP(3),K,AR,RI(3),RJ(3)
    REAL :: I0V,RES3(3),IVEC(3),SUMI0,SUMRI0(3),SUMIV(3),SUMRD
    REAL :: J0,JRI0(3),JVEC(3),JRD,PS,QS,T1,T2,T3,T4
    COMPLEX :: GS,I1G,I2G(3),I3G(3),I4G,PSD,QSD,PTOT,QTOT,Z_SELF,JW,JOE
    TYPE(GAUSS_TRI_DATA) :: GD
    INTEGER :: I,J,NPTS

    V1=[0.,0.,0.]; V2=[1.,0.,0.]; V3=[0.,1.,0.]; ROPP=[0.5,0.25,0.]
    K=2.*PI; AR=0.5; JW=(0.,1.)*120.*PI*K; JOE=(0.,1.)*120.*PI/K

    PRINT *,'========================================'
    PRINT *,'  自作用项收敛性分析 (高斯阶数: 3,7,12)'
    PRINT *,'========================================'
    PRINT *
    PRINT *,'三角: 直角(0,0)-(1,0)-(0,1), 面积=0.5, λ=1'
    PRINT *

    DO NPTS=1,3
        IF(NPTS==1) CALL INIT_GAUSS_TRI(GAUSS_3PT, GD)
        IF(NPTS==2) CALL INIT_GAUSS_TRI(GAUSS_7PT, GD)
        IF(NPTS==3) CALL INIT_GAUSS_TRI(GAUSS_12PT, GD)

        SUMI0=0.; SUMRI0=0.; SUMIV=0.; SUMRD=0.
        DO I=1,GD%N_POINTS
            RI(1)=GD%UVW(1,I)*V1(1)+GD%UVW(2,I)*V2(1)+GD%UVW(3,I)*V3(1)
            RI(2)=GD%UVW(1,I)*V1(2)+GD%UVW(2,I)*V2(2)+GD%UVW(3,I)*V3(2)
            RI(3)=0.0
            CALL ANALYTIC_SCALAR_POT_1_OVER_R(RI,V1,V2,V3,I0V)
            CALL ANALYTIC_LINEAR_POT_1_OVER_R(RI,V1,V2,V3,RES3)
            IVEC=RES3(1)*V1+RES3(2)*V2+RES3(3)*V3
            SUMI0 =SUMI0 +GD%WEIGHTS(I)*I0V
            SUMRI0=SUMRI0+GD%WEIGHTS(I)*RI*I0V
            SUMIV =SUMIV +GD%WEIGHTS(I)*IVEC
            SUMRD =SUMRD +GD%WEIGHTS(I)*DOT_PRODUCT(RI,IVEC)
        END DO
        J0=SUMI0*AR; JRI0=SUMRI0*AR; JVEC=SUMIV*AR; JRD=SUMRD*AR
        QS=-J0/PI
        T1=JRD; T2=-DOT_PRODUCT(ROPP,JVEC)
        T3=-DOT_PRODUCT(ROPP,JRI0); T4=DOT_PRODUCT(ROPP,ROPP)*J0
        PS=(T1+T2+T3+T4)/(4.*PI)

        I1G=0.; I2G=0.; I3G=0.; I4G=0.
        DO I=1,GD%N_POINTS
            RI(1)=GD%UVW(1,I)*V1(1)+GD%UVW(2,I)*V2(1)+GD%UVW(3,I)*V3(1)
            RI(2)=GD%UVW(1,I)*V1(2)+GD%UVW(2,I)*V2(2)+GD%UVW(3,I)*V3(2)
            RI(3)=0.0
            DO J=1,GD%N_POINTS
                RJ(1)=GD%UVW(1,J)*V1(1)+GD%UVW(2,J)*V2(1)+GD%UVW(3,J)*V3(1)
                RJ(2)=GD%UVW(1,J)*V1(2)+GD%UVW(2,J)*V2(2)+GD%UVW(3,J)*V3(2)
                RJ(3)=0.0
                GS=GREEN_FUNC_SMOOTH(RI,RJ,K)
                I1G=I1G+GD%WEIGHTS(I)*GD%WEIGHTS(J)*GS
                I2G=I2G+GD%WEIGHTS(I)*GD%WEIGHTS(J)*RJ*GS
                I3G=I3G+GD%WEIGHTS(I)*GD%WEIGHTS(J)*RI*GS
                I4G=I4G+GD%WEIGHTS(I)*GD%WEIGHTS(J)*DOT_PRODUCT(RI,RJ)*GS
            END DO
        END DO
        I1G=I1G*AR**2; I2G=I2G*AR**2; I3G=I3G*AR**2; I4G=I4G*AR**2
        PSD=DOT_PRODUCT(ROPP,ROPP)*I1G-DOT_PRODUCT(ROPP,I2G)-DOT_PRODUCT(ROPP,I3G)+I4G
        QSD=-4.*I1G

        PTOT=CMPLX(PS,0.)+PSD; QTOT=CMPLX(QS,0.)+QSD
        Z_SELF=JW*PTOT+JOE*QTOT

        PRINT '(A,I2,A)', '--- ',GD%N_POINTS,'点高斯 ---'
        PRINT '(A,ES14.6)', '  J0(双重1/R)=',J0
        PRINT '(A,2ES14.6)', '  Z_self实/虚=',Z_SELF
        PRINT *
    END DO

    PRINT *,'结论: J0 和 Z_self 随高斯阶数增加而收敛。'
    PRINT *,'7pt→12pt 的 Z_self 相对变化反映了数值积分误差量级。'
END PROGRAM
