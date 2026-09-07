subroutine Comp_eltsI(I1,I2,I3,I4,nonsingularity,avi,rjf,ris,r1s,r2s,r3s,Np,Nq,weightF,weightS)
	 use Common_constants,only: k_wnum
  implicit none
  complex,intent(out):: I1,I2(3),I3(3),I4
  logical,intent(in):: nonsingularity
  real,intent(in):: avi
  integer,intent(in):: Np,Nq
  real,intent(in):: rjf(3,Np),ris(3,Nq),r1s(3),r2s(3),r3s(3)
  real,intent(in):: weightF(Np),weightS(Nq)

  
  integer p,q
	 real R,wp,wq,rp(3),rq(3),ns(3)
  real:: t3,t4(3),w0
  complex F1_temp,ctmp,ctmpv(3),inner,innerv(3)
  real,external:: re_Norm
  
  I1=(0.0,0.0)
  I2=(0.0,0.0)
  I3=(0.0,0.0)
  I4=(0.0,0.0)
    
  do p=1,Np
    rp=rjf(:,p)
    wp=weightF(p)
    inner=0
    innerv=0

    if(nonsingularity) then
      do q=1,Nq
        wq=weightS(q)
        rq=ris(:,q)
        R=re_Norm(rp-rq)
        ctmp=wq*Green(R)
        inner=inner+ctmp
        innerv=innerv+ctmp*rq
      enddo

    else
      do q=1,Nq
        wq=weightS(q)
        rq=ris(:,q)
        R=re_Norm(rp-rq)
        F1_temp=wq*F1(R)
        inner=inner+F1_temp
        innerv=innerv+F1_temp*rq
	     enddo
   
      call Spara_cal(rp,r1s,r2s,r3s,avi,t3,t4,w0,ns)
      t3=t3/avi
      t4=t4/avi
      innerv=innerv+(rp-ns*w0)*t3
    
 
      inner=inner+t3
      innerv=innerv+t4
    endif
    inner=inner*wp
    innerv=innerv*wp
    
    I1=I1+inner
    I2=I2+inner*rp
    I3=I3+innerv
    I4=I4+dot_product(rp,innerv)
  enddo

contains
  pure complex function F1(R)
	   use Common_constants
	   implicit none
	   real,intent(in):: R   
	   integer i
    real x
    complex,parameter:: GreenExpan_coeff(0:7)=(/(0,-1.0),(-0.5,0),(0,0.1666667),(4.1666667e-2,0),(0.0,-8.3333333e-3),&
                                                 (-1.3888889e-3,0),(0,1.9841270e-4),(2.4801587e-5,0)/)

    x=k_wnum*R
    F1=GreenExpan_coeff(7)
    do i=6,0,-1
      F1=F1*x
	     F1=F1+GreenExpan_coeff(i)
	   enddo
    F1=k_wnum*F1
  endfunction
 
 
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !cccc----------------Function to calculate the Green's function in free space---------------------------cccc
  pure complex function Green(R)
	   use Common_constants
	   implicit none
	   real,intent(in):: R
    Green=exp(cmplx(0.0,-k_wnum*R))/R
  endfunction
  
  
  
  subroutine Spara_cal(r,P1,P2,P3,areai,t3,t4,w0,ns)
    implicit none
    real,intent(in):: r(3),P1(3),P2(3),P3(3),areai
	   real,intent(out)::ns(3)
	   real,intent(out)::t3,t4(3),w0
  
	   integer m
    logical singular
	   real:: abs_l(3)
    real:: l1(3),l2(3),l3(3),m1(3),m2(3),m3(3)
	   real:: u(3),v(3),w(3),u0,v0,u3,v3
	   real:: lm(3),lp(3),P0(3),R0,Rm,Rp
	   real:: r_P1(3),P3_P1(3)
	   real:: g2,g3,beta,t1(3),t2(3)
  
    call unit_normal(P1,P2,P3,ns)                                          !unit normal of the source triangle
      
	   abs_l(1)=re_Norm(P3-P2)                                         !length of the 1st edge of the source triangle
	   abs_l(2)=re_Norm(P1-P3)                                         !length of the 2nd edge of the source triangle
	   abs_l(3)=re_Norm(P2-P1)                                         !length of the 3rd edge of the source triangle
	
	   l1=(P3-P2)/abs_l(1)                                          !direction of the 1st edge of the source triangle
	   l2=(P1-P3)/abs_l(2)                                          !direction of the 2nd edge of the source triangle
	   l3=(P2-P1)/abs_l(3)                                          !direction of the 3rd edge of the source triangle


    call Re_cross_product(l1,ns,m1)                                        !outer unit normal of the 1st edge of the source triangle
    call Re_cross_product(l2,ns,m2)                                        !outer unit normal of the 2nd edge of the source triangle
    call Re_cross_product(l3,ns,m3)                                        !outer unit normal of the 3rd edge of the source triangle

  
	   u=l3                                                             !1st direction of the local coordinate 
	   w=ns

    call Re_cross_product(ns,u,v)                                          !2nd direction of the local coordinate    
  
    
	   r_P1=r-P1
	   P3_P1=P3-P1

	   u0=dot_product(u,r_P1)                                           !local coordinate of projection of the observation point (U direction)
	   v0=dot_product(v,r_P1)                                           !local coordinate of projection of the observation point (V direction)
	   w0=dot_product(w,r_P1)                                           !local coordinate of projection of the observation point (W direction)
	   u3=dot_product(u,P3_P1)                                          !local coordinate of 3rd vertice of the source triangle (U direction)
	   v3=2.0*areai/abs_l(3)                                                  !local coordinate of 3rd vertice of the source triangle (V direction)

	   lm(1)=-((abs_l(3)-u0)*(abs_l(3)-u3)+v0*v3)/abs_l(1)
	   lm(2)=-(u3*(u3-u0)+v3*(v3-v0))/abs_l(2)
	   lm(3)=-u0
	   lp(1)=((u3-u0)*(u3-abs_l(3))+v3*(v3-v0))/abs_l(1)
	   lp(2)=(u0*u3+v0*v3)/abs_l(2)
	   lp(3)=abs_l(3)-u0
	   P0(1)=(v0*(u3-abs_l(3))+v3*(abs_l(3)-u0))/abs_l(1)
	   P0(2)=(u0*v3-v0*u3)/abs_l(2)
	   P0(3)=v0

	   do m=1,3
	     R0=P0(m)**2+w0**2
	     Rm=sqrt(lm(m)**2+P0(m)**2+w0**2)
	     Rp=sqrt(lp(m)**2+P0(m)**2+w0**2)

      singular=.FALSE.
      if(lm(m)>0.and.lp(m)>0) then
        g2=log((Rp+lp(m))/(Rm+lm(m)))
        g3=(lp(m)*Rp-lm(m)*Rm)+R0*g2
      elseif(lm(m)<=0.and.lp(m)>=0) then
        if(Rm+lm(m)==0) then
          singular=.TRUE.
        else
          g2=log((Rp+lp(m))/(Rm+lm(m)))
          g3=(lp(m)*Rp-lm(m)*Rm)+R0*g2
        endif
      elseif(lm(m)<0.and.lp(m)<0) then
        g2=log((Rm-lm(m))/(Rp-lp(m)))
        g3=(lp(m)*Rp-lm(m)*Rm)+R0*g2
      endif

    
      if(R0<=(0.0e-2*abs_l(m))**2) then
        beta=0.0
      else
        beta=(atan(P0(m)*lp(m)/(R0+abs(w0)*Rp))-atan(P0(m)*lm(m)/(R0+abs(w0)*Rm)))
      endif
    
    
      if(.not.singular) then
        t1(m)=P0(m)*g2-abs(w0)*beta
        t2(m)=0.5*g3
      else
        t1(m)=-abs(w0)*beta
        t2(m)=0.5*((lp(m)*Rp-lm(m)*Rm))
      endif
	   enddo
    t3=sum(t1)
    t4=m1*t2(1)+m2*t2(2)+m3*t2(3)

  endsubroutine
  
end