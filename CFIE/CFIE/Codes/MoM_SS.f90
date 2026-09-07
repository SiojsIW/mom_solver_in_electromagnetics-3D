subroutine MoM_SS
  use Common_array
  use Common_constants
  use Dunavant_integration_formulism
  implicit none
  integer n,m,node,node1,node2,ni,mj,p,q,k
  integer i,j

  real(4) rcm(3),rcn(3),dist
  real r1f(3),r2f(3),r3f(3),nSm(3),rjf(3,SNum_samp),r1s(3),r2s(3),r3s(3),ris(3,SNum_samp)
  real Am,An,lj,li,rmj(3),rni(3),rtmpv(3)
  real Jrm(3,maxRWG),Irn(3,maxRWG),Jlengthm(maxRWG),Ilengthn(maxRWG)
  real coef3
  real A1
  complex coef1,coef2
  complex Amn_E,Amn_H,Zji,I1ss,I2ss(3),I3ss(3),I4ss
  complex U1ss(3),U2ss(3),U3ss,U4ss(3),U5ss(3)
  complex ctmp,ctmp1
    
  real,external:: re_Norm
  integer*4 time1,time2,time
  
  call SYSTEM_CLOCK(time1)
  write(*,*) "Start to compute the SS part of the impedance matrix, on the way..."
  
  allocate(SparseValue(Num_side,Num_side))
  SparseValue=0
      
  coef1=cmplx(0.0,omega*MU0/(16.0*PI))
  coef2=cmplx(0.0,1.0/(4.0*PI*omega*EPS0))
  coef3=1.0/(16*PI)

 !$OMP parallel do default(none) shared(triangle_centroid,point_cor,triangle_point,threshold_R,side_index, SparseValue,triangle_area,side_length, &
 !$OMP Num_side,Num_Triangle,triangle_rwgs,CFSIE,coef1,coef2,coef3)&
 !$OMP private(n,m,rcm,rcn,Am,An,r1f,r2f,r3f,nSm,r1s,r2s,r3s,dist,I1ss,I2ss,I3ss,I4ss,ctmp,ctmp1,U1ss,U2ss,U3ss,U4ss,U5ss,&
 !$OMP  mj,ni,i,j,node1,node2,node,p,q,rmj,rni,lj,li,Amn_E,Amn_H,Zji,rjf,ris,A1,rtmpv,Jrm,Jlengthm,Ilengthn,Irn) schedule(guided)
  do m=1,Num_Triangle    
    rcm=triangle_centroid(:,m)
	   Am=triangle_area(m)
 
    r1f=point_cor(:,triangle_point(1,m))
    r2f=point_cor(:,triangle_point(2,m))
    r3f=point_cor(:,triangle_point(3,m))

	   do k=1,SNum_samp
      rjf(:,k)=epsilS(k)*r1f+etaS(k)*r2f+xiS(k)*r3f
    enddo
      
    if(CFSIE) then
      call unit_normal(r1f,r2f,r3f,nSm)
      A1=1.0/6.0*(dot_product(r1f,r1f)+dot_product(r2f,r2f)+dot_product(r3f,r3f)+dot_product(r1f,r2f)+dot_product(r1f,r3f)+dot_product(r2f,r3f))  
    endif

 
    do mj=1,triangle_rwgs(m)%N_rwg
      j=triangle_rwgs(m)%rwg(mj)
        
      if(side_index(3,j)==m) then
        !flagm=1.0
        Jlengthm(mj)=side_length(j)
      else  !if(side_index(4,j)==m) then
        !flagm=-1.0
        Jlengthm(mj)=-side_length(j)
      endif
        
   	  node1=side_index(1,j); node2=side_index(2,j)
      do p=1,3
        node=triangle_point(p,m)
        if((node/=node1).and.(node/=node2)) exit
      enddo         
      Jrm(:,mj)=point_cor(:,node)   	          
    enddo 
 
    do n=1,Num_Triangle
      rcn=triangle_centroid(:,n)
      An=triangle_area(n)

      r1s=point_cor(:,triangle_point(1,n))     !source triangle coordinates
      r2s=point_cor(:,triangle_point(2,n))
      r3s=point_cor(:,triangle_point(3,n))

      do k=1,SNum_samp
	       ris(:,k)=epsilS(k)*r1s+etaS(k)*r2s+xiS(k)*r3s
      enddo
                    
      dist=re_Norm(rcm-rcn)          
      
      call Comp_eltsI(I1ss,I2ss,I3ss,I4ss,dist>=threshold_R,An,rjf,ris,r1s,r2s,r3s,SNum_samp,SNum_samp,weightS,weightS)
       
      if(CFSIE.and.m/=n) then
        call Comp_eltsU(U1ss,U2ss,U3ss,U4ss,U5ss,dist>=threshold_R,An,nSm,rjf,ris,r1s,r2s,r3s,SNum_samp,SNum_samp,weightS,weightS)
      endif
      
      do ni=1,triangle_rwgs(n)%N_rwg
        i=triangle_rwgs(n)%rwg(ni)
		      if(side_index(3,i)==n) then
	         !flagn=1.0
          Ilengthn(ni)=side_length(i)
        else  !if(side_index(4,i)==n) then
	         !flagn=-1.0
          Ilengthn(ni)=-side_length(i)
	       endif
              
	       node1=side_index(1,i); node2=side_index(2,i)
	       do q=1,3
          node=triangle_point(q,n)
          if((node/=node1).and.(node/=node2)) exit
	       enddo
        Irn(:,ni)=point_cor(:,node)
      enddo  
     
      do mj=1,triangle_rwgs(m)%N_rwg
        j=triangle_rwgs(m)%rwg(mj)
        rmj=Jrm(:,mj)       
        lj=Jlengthm(mj) 
 
        ctmp=coef1*(I4ss-dot_product(rmj,I3ss))-coef2*I1ss
        if(CFSIE.and.n/=m) then
          
          call Re_cross_product(rmj,nSm,rtmpv)
          ctmp1=-dot_product(rtmpv,U2ss)-dot_product(nSm,U4ss)
        endif        
        do ni=1,triangle_rwgs(n)%N_rwg
          i=triangle_rwgs(n)%rwg(ni)
                    
          rni=Irn(:,ni)              
          li=Ilengthn(ni)
                    
          Amn_E=lj*li*(coef1*(dot_product(rmj,rni)*I1ss-dot_product(rni,I2ss))+ctmp)
   
          if(CFSIE) then
            if(m/=n) then
              Amn_H=coef3*lj*li*(dot_product(dot_product(rmj,rni)*nSm-dot_product(rni,nSm)*rmj,U1ss)&
                  +dot_product(rni,nSm)*U3ss-dot_product(rni,U5ss)+ctmp1)
            else
	             Amn_H=0.125*lj*li/Am*(A1+dot_product(rmj,rni)-dot_product(triangle_centroid(:,m),rmj+rni))
            endif
             
            Zji=alpha*Amn_E+ETA0*(1.0-alpha)*Amn_H
          else
            Zji=Amn_E
          endif
            
          !$OMP ATOMIC
          SparseValue(j,i)=SparseValue(j,i)+Zji     
        
	       enddo
      enddo
    enddo

	 enddo

  call SYSTEM_CLOCK(time2,time)

  write(*,"('Finished, the costing time is',F9.3,' sec',/)") (time2-time1)/float(time)

end



subroutine unit_normal(a,b,c,n)
	 implicit none
	 real,intent(in):: a(3),b(3),c(3)
  real,intent(out):: n(3)
  real rtmp1(3),rtmp2(3)
  
	 call Re_cross_product(b-a,c-a,n)

	 n=n/sqrt(sum(n*n))
end