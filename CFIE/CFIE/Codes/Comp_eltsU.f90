subroutine Comp_eltsU(U1,U2,U3,U4,U5,nonsingularity,AVn,nSm,rjf,ris,r1s,r2s,r3s,Np,Nq,weightF,weightS)
  use Common_constants
  implicit none
  complex,intent(out):: U1(3),U2(3),U3,U4(3),U5(3)
  logical,intent(in):: nonsingularity
  real,intent(in):: AVn
  integer,intent(in):: Np,Nq
  real,intent(in):: nSm(3),rjf(3,Np),ris(3,Nq),r1s(3),r2s(3),r3s(3)
  real,intent(in):: weightF(Np),weightS(Nq)

  
  integer p
  real wp,rp(3),rtmpv(3)
  complex gradG(3),ctmpv(3),ctmpv1(3)
  real,external:: re_Norm
    
  U1=(0.0,0.0)
  U2=(0.0,0.0)
  U3=(0.0,0.0)
  U4=(0.0,0.0)
  U5=(0.0,0.0)


	 do p=1,Np

    rp=rjf(:,p)
    wp=weightF(p)

    gradG=wp*gradG_S(AVn,rp,ris,r1s,r2s,r3s,nonsingularity) 

    
    call ReCom_cross_product(rp,gradG,ctmpv)
    
    U1=U1+gradG
    U2=U2+ctmpv
    U3=U3+dot_product(rp,gradG)
     
    call ReCom_cross_product(rp,ctmpv,ctmpv1)
    U4=U4+ctmpv1
    
    U5=U5+rp*(dot_product(nSm,gradG))
  enddo

contains

  subroutine ReCom_cross_product(a,b,c)
    implicit none
    real a(3)
    complex b(3),c(3)
    c=(/(a(2)*b(3)-a(3)*b(2)),(a(3)*b(1)-a(1)*b(3)),(a(1)*b(2)-a(2)*b(1))/)
  endsubroutine

  function gradG_S(An,rp,ris,P1,P2,P3,FarNear) result(gradG)
    use Common_constants
    use Dunavant_integration_formulism,only: SNum_samp
    implicit none
    real,intent(in):: An
    real,intent(in):: rp(3),ris(3,SNum_samp),P1(3),P2(3),P3(3)
    logical,intent(in)::  FarNear
    complex gradG(3)
  
    integer q,i
    real wq
    real I_31(3),I_32(3)
    real rq(3),R(3),abs_R,ns(3)
    real:: w0,t1(3),t3(3,3),t4(3,3),t5(3,3)
    complex I_30(3)
    complex temp_gradient_G(3)

    if(FarNear) then

      gradG=(0.0,0.0)

      do q=1,SNum_samp
        rq=ris(:,q)
	       R=rp-rq
        abs_R=re_Norm(R)
        temp_gradient_G=-cmplx(1.0,k_wnum*abs_R)*exp(cmplx(0.0,-k_wnum*abs_R))/(abs_R**3)*R
        wq=weightS(q)
        gradG=gradG+wq*temp_gradient_G
      enddo

    else

      I_30=(0.0,0.0)
      !.... ordinary term
      do q=1,SNum_samp
        rq=ris(:,q)
	       R=rp-rq
        abs_R=re_Norm(R)

        temp_gradient_G=R*cmplx(k_wnum*abs_R*(1.0-k_wnum**2*abs_R**2/18.0)/8.0,(1.0-0.1*k_wnum**2*abs_R**2)/3.0)
        wq=weightS(q)
        I_30=I_30+wq*temp_gradient_G
      enddo
      
      !.... singular term: PV
      call Spara_cal_2(rp,P1,P2,P3,An,t1,t3,t4,t5,w0,ns)  
   
      I_31=0.0
      if(abs(w0)<=ratio_limit1) then  ! w0 is zero
        do i=1,3
          I_31=I_31+t5(:,i)
        enddo
      else
        do i=1,3
          I_31=I_31+t4(:,i)+t5(:,i)
        enddo
      endif

      I_32=0.0
      if(abs(w0)<=ratio_limit1) then   ! w0 is zero
        do i=1,3
          I_32=I_32-t3(:,i)
        enddo
      else
        do i=1,3
          I_32=I_32+t1(i)*w0*ns-t3(:,i)
        enddo
      endif

      gradG=k_wnum**3*I_30-(I_31+k_wnum**2/2.0*I_32)/An
  
    endif

  endfunction

  
  
  subroutine Spara_cal_2(r,P1,P2,P3,An,t1,t3,t4,t5,w0,ns)
    use Common_constants
    implicit none
	   real,intent(in):: r(3),P1(3),P2(3),P3(3),An
    real,intent(out):: t1(3),t3(3,3),t4(3,3),t5(3,3),w0
    real,intent(out):: ns(3)
  
	   integer i
	   real:: abs_l(3)
    real:: l1(3),l2(3),l3(3),m1(3),m2(3),m3(3)
	   real:: u(3),v(3),w(3),u0,v0,u3,v3
	   real:: lm(3),lp(3),P0(3),R0(3),Rm(3),Rp(3)
	   real:: r_P1(3),P3_P1(3)
	   real:: g2(3),beta(3)
    real:: m_u(3,3)
    
    call unit_normal(P1,P2,P3,ns)                                          !unit normal of the source triangle
      
	   abs_l(1)=re_Norm(P3-P2)
	   abs_l(2)=re_Norm(P1-P3)
	   abs_l(3)=re_Norm(P2-P1)
	
   !	aver_len=(abs_l(1)+abs_l(2)+abs_l(3))/3.0


    l1 = (P3-P2)/abs_l(1)                                          !direction of the 1st edge of the source triangle
    l2 = (P1-P3)/abs_l(2)                                          !direction of the 2nd edge of the source triangle
    l3 = (P2-P1)/abs_l(3)                                          !direction of the 3rd edge of the source triangle


    call Re_cross_product(l1,ns,m1)                                        !outer unit normal of the 1st edge of the source triangle
    call Re_cross_product(l2,ns,m2)                                        !outer unit normal of the 2nd edge of the source triangle
    call Re_cross_product(l3,ns,m3)                                        !outer unit normal of the 3rd edge of the source triangle

  
    m_u(:,1)=m1
    m_u(:,2)=m2
    m_u(:,3)=m3

    u = l3                                                            !1st direction of the local coordinate 
    w = ns

    call Re_cross_product(ns,u,v)                                          !2nd direction of the local coordinate

  
    r_P1  = r-P1
    P3_P1 = P3-P1


	   u0=dot_product(u,r_P1)                                           !local coordinate of projection of the observation point (U direction)
	   v0=dot_product(v,r_P1)                                           !local coordinate of projection of the observation point (V direction)
	   w0=dot_product(w,r_P1)                                           !local coordinate of projection of the observation point (W direction)
	   
    
    u3=dot_product(u,P3_P1)                                          !local coordinate of 3rd vertice of the source triangle (U direction)
	   v3=2.0*An/abs_l(3)                                                  !local coordinate of 3rd vertice of the source triangle (V direction)

	   lm(1)=-((abs_l(3)-u0)*(abs_l(3)-u3)+v0*v3)/abs_l(1)
	   lm(2)=-(u3*(u3-u0)+v3*(v3-v0))/abs_l(2)
	   lm(3)=-u0
	   lp(1)=((u3-u0)*(u3-abs_l(3))+v3*(v3-v0))/abs_l(1)
	   lp(2)=(u0*u3+v0*v3)/abs_l(2)
	   lp(3)=abs_l(3)-u0
	   P0(1)=(v0*(u3-abs_l(3))+v3*(abs_l(3)-u0))/abs_l(1)
	   P0(2)=(u0*v3-v0*u3)/abs_l(2)
	   P0(3)=v0


 	  do i=1,3
	     R0(i)=sqrt(P0(i)**2+w0**2)
	     Rm(i)=sqrt(lm(i)**2+P0(i)**2+w0**2)
	     Rp(i)=sqrt(lp(i)**2+P0(i)**2+w0**2)
	   enddo

    if(abs(w0)<=ratio_limit1) then  ! w0 is 0
      do i=1,3
        if(lm(i)>0.and.lp(i)>0) then  !... case 1
          g2(i)=log((Rp(i)+lp(i))/(Rm(i)+lm(i)))

          t3(:,i)=0.5*m_u(:,i)*(R0(i)**2*g2(i)+Rp(i)*lp(i)-Rm(i)*lm(i))
          t5(:,i)=m_u(:,i)*g2(i)
        elseif(lm(i)<=0.and.lp(i)>=0) then  !... case 2, includes lm==0 or lp==0
          if(Rm(i)+lm(i)==0) then  !... due to round-off error, or P0 is actually zero (when w0 is actually zero)
            t3(:,i)=0.5*m_u(:,i)*(Rp(i)*lp(i)-Rm(i)*lm(i))
            t5(:,i)=0.0*m_u(:,i)  !***
       !   elseif(Rm(i)+lm(i)<0) then
        !    write (*,*) 'Error! error code: case 2 Rm(i)+lm(i)<0'
    
          else
            g2(i)=log((Rp(i)+lp(i))/(Rm(i)+lm(i)))

            t3(:,i)=0.5*m_u(:,i)*(R0(i)**2*g2(i)+Rp(i)*lp(i)-Rm(i)*lm(i))
            t5(:,i)=m_u(:,i)*g2(i)
          endif

        elseif(lm(i)<0.and.lp(i)<0) then  !... case 3, converted into case 1
          g2(i)=log((Rm(i)-lm(i))/(Rp(i)-lp(i)))

          t3(:,i)=0.5*m_u(:,i)*(R0(i)**2*g2(i)+Rp(i)*lp(i)-Rm(i)*lm(i))
          t5(:,i)=m_u(:,i)*g2(i)
        endif
            
      enddo
    else  ! w0 is not 0
      do i=1,3
        if(lm(i)>0.and.lp(i)>0) then  !... case 1
          g2(i)=log((Rp(i)+lp(i))/(Rm(i)+lm(i)))

	      !******************************
          if(R0(i)<=ratio_limit2*abs_l(i)) then  !... 2011.10.14, atan(finite/0)=nan
            beta(i)=0.0 ! actually , this is to make abs(w0)*beta(i)=0.0
          else
            beta(i)=(atan(P0(i)*lp(i)/(R0(i)**2+abs(w0)*Rp(i)))-atan(P0(i)*lm(i)/(R0(i)**2+abs(w0)*Rm(i))))
          endif
          !*******************************

	         t1(i)=P0(i)*g2(i)-abs(w0)*beta(i)

          t4(:,i)=ns*sgn(w0)*beta(i)

          t3(:,i)=0.5*m_u(:,i)*(R0(i)**2*g2(i)+Rp(i)*lp(i)-Rm(i)*lm(i))
          t5(:,i)=m_u(:,i)*g2(i)

        elseif(lm(i)<=0.and.lp(i)>=0) then  !... case 2, includes lm==0 or lp==0
          if(Rm(i)+lm(i)==0) then  !... due to round-off error, or P0 is actually zero (when w0 is actually zero)

	          !******************************
            if(R0(i)<=ratio_limit2*abs_l(i)) then  !... 2011.10.14, atan(finite/0)=nan
              beta(i)=0.0 ! actually , this is to make abs(w0)*beta(i)=0.0
            else
              beta(i)=(atan(P0(i)*lp(i)/(R0(i)**2+abs(w0)*Rp(i)))-atan(P0(i)*lm(i)/(R0(i)**2+abs(w0)*Rm(i))))
            endif
            !*******************************
            t1(i)=-abs(w0)*beta(i)

            t4(:,i)=ns*sgn(w0)*beta(i)

            t3(:,i)=0.5*m_u(:,i)*(Rp(i)*lp(i)-Rm(i)*lm(i))
            t5(:,i)=0.0*m_u(:,i)  !***

      !    elseif(Rm(i)+lm(i)<0) then
      !      write (*,*) 'Error! error code: case 2 Rm(i)+lm(i)<0'
   
          else
            g2(i)=log((Rp(i)+lp(i))/(Rm(i)+lm(i)))

	        !******************************
            if(R0(i)<=ratio_limit2*abs_l(i)) then  !... 2011.10.14, atan(finite/0)=nan
              beta(i)=0.0 ! actually , this is to make abs(w0)*beta(i)=0.0
            else
              beta(i)=(atan(P0(i)*lp(i)/(R0(i)**2+abs(w0)*Rp(i)))-atan(P0(i)*lm(i)/(R0(i)**2+abs(w0)*Rm(i))))
            endif
            !*******************************

	           t1(i)=P0(i)*g2(i)-abs(w0)*beta(i)

            t4(:,i)=ns*sgn(w0)*beta(i)

            t3(:,i)=0.5*m_u(:,i)*(R0(i)**2*g2(i)+Rp(i)*lp(i)-Rm(i)*lm(i))
            t5(:,i)=m_u(:,i)*g2(i)
          endif

        elseif(lm(i)<0.and.lp(i)<0) then  !... case 3, converted into case 1
          g2(i)=log((Rm(i)-lm(i))/(Rp(i)-lp(i)))

	      !******************************
          if(R0(i)<=ratio_limit2*abs_l(i)) then  !... 2011.10.14, atan(finite/0)=nan
            beta(i)=0.0 ! actually , this is to make abs(w0)*beta(i)=0.0
          else
            beta(i)=(atan(P0(i)*lp(i)/(R0(i)**2+abs(w0)*Rp(i)))-atan(P0(i)*lm(i)/(R0(i)**2+abs(w0)*Rm(i))))
          endif
          !*******************************

	         t1(i)=P0(i)*g2(i)-abs(w0)*beta(i)

          t4(:,i)=ns*sgn(w0)*beta(i)

          t3(:,i)=0.5*m_u(:,i)*(R0(i)**2*g2(i)+Rp(i)*lp(i)-Rm(i)*lm(i))
          t5(:,i)=m_u(:,i)*g2(i)

        endif        
      enddo
    endif
  endsubroutine
  
  pure	real function sgn(a)
	  implicit none
	  real,intent(in):: a
   if(a>0.0) then
     sgn=1.0
   elseif(a<0.0) then
     sgn=-1.0
   else
     sgn=0.0
   endif
	endfunction
 
 
  pure complex function F1(R)
	   use Common_constants
     use Dunavant_integration_formulism
	   implicit none
	   real,intent(in):: R   
	   integer i
    real x
    complex,parameter:: GreenExpan_coeff(0:7)=(/(0.0,-1.0),(-0.5,0.0),(0.0,0.166666666666667),(4.166666666666667e-002,0.0),(0.0,-8.333333333333333e-003),&
                                                (-1.388888888888889e-003,0.0),(0.0,1.984126984126984E-004),(2.480158730158730e-005,0.0)/)

    x=k_wnum*R
    F1=GreenExpan_coeff(7)
    do i=6,0,-1
      F1=F1*x
	     F1=F1+GreenExpan_coeff(i)
	   enddo
    F1=k_wnum*F1
	 endfunction
 
end