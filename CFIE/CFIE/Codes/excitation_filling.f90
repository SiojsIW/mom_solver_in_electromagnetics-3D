subroutine excitation_filling
  use Common_array
  use Common_constants
  use Far_patterns
  implicit none
  integer j,k,node,node1,node2,node3,lj,lm,sideindexj,faceindexj
  real:: r1f(3),r2f(3),r3f(3),nf(3)
  real:: rm(3),nf_rm(3)
  real:: lengthm,aream
  complex:: V_E1,V_E2(3),V_H1,V_H2(3),I5,H_I6
  complex:: Vm_E,Vm_H,rhs
  
  write(*,*) "Start to compute the excitation vector"
 
  
  vec_k(1)=sin(theta_inc)*cos(phi_inc)
	 vec_k(2)=sin(theta_inc)*sin(phi_inc)
	 vec_k(3)=cos(theta_inc)
  vec_k=-k_wnum*vec_k

  E0(1)= E_inc_0_theta*cos(theta_inc)*cos(phi_inc)-E_inc_0_phi*sin(phi_inc)      
  E0(2)= E_inc_0_theta*cos(theta_inc)*sin(phi_inc)+E_inc_0_phi*cos(phi_inc)        
  E0(3)=-E_inc_0_theta*sin(theta_inc)

  
  H0(1)=-(sin(theta_inc)*sin(phi_inc)*E0(3)-cos(theta_inc)*E0(2))/ETA0
  H0(2)=-(cos(theta_inc)*E0(1)-sin(theta_inc)*cos(phi_inc)*E0(3))/ETA0
  H0(3)=-(sin(theta_inc)*cos(phi_inc)*E0(2)-sin(theta_inc)*sin(phi_inc)*E0(1))/ETA0

  allocate(excitation(Num_side))
  
  excitation=0
  
  if(Num_Triangle>0) then
    do j=1,Num_Triangle 
    
      r1f = point_cor(:,triangle_point(1,j))        !field triangle coordinates
      r2f = point_cor(:,triangle_point(2,j))
      r3f = point_cor(:,triangle_point(3,j))
      
      if(CFSIE) then
        call unit_normal(r1f,r2f,r3f,nf)                                          !unit normal of the source triangle
        call Vs_filling_term_by_PW_EH(j,r1f,r2f,r3f,nf,V_E1,V_E2,V_H1,V_H2)
	     else
        call Vs_filling_term_by_PW_E(j,r1f,r2f,r3f,V_E1,V_E2)
      endif
   
      do lj=1,triangle_rwgs(j)%N_rwg
        sideindexj=triangle_rwgs(j)%rwg(lj)
      
	       if(side_index(3,sideindexj)==j) then
          !flag=1.0
	         lengthm=side_length(sideindexj)
	       else  !if(side_index(4,sideindexj)==j) then
          !flag=-1.0
	         lengthm=-side_length(sideindexj)
	       endif
	      
	       node1=side_index(1,sideindexj); node2=side_index(2,sideindexj)
	       do lm=1,3
          node=triangle_point(lm,j)
	         if((node/=node1).and.(node/=node2)) exit
	       enddo
        rm=point_cor(:,node)


        I5=V_E1-dot_product(rm,V_E2)
        Vm_E=lengthm/2.0*I5
    
        
        if(CFSIE) then
          call Re_cross_product(nf,rm,nf_rm)
	         H_I6=-V_H1-dot_product(nf_rm,V_H2)
	         Vm_H=-lengthm/2.0*H_I6     
          rhs=alpha*Vm_E+ETA0*(1.0-alpha)*Vm_H   
        else
          rhs=Vm_E
        endif
        excitation(sideindexj)=excitation(sideindexj)+rhs
	     enddo
    enddo
  endif
 
 
  
  write(*,"('Finished, ',/)")
  
contains
  !_______________________________________________________________________________
  ! Subroutines to calculate terms used in the Vs filling with planewave incidence.
  ! 
  !
  !_______________________________________________________________________________


  subroutine Vs_filling_term_by_PW_E(j,Sr1f,Sr2f,Sr3f,V_E1,V_E2)
    use Dunavant_integration_formulism
    implicit none
    integer,intent(in):: j
    real,intent(in):: Sr1f(3),Sr2f(3),Sr3f(3)
    complex,intent(out):: V_E1,V_E2(3)
    
    integer jj
    real rjf(3)
    complex Einc(3),Hinc(3)

    V_E1=(0.0,0.0)
    V_E2=(0.0,0.0)

	   do jj=1,SNum_samp
      rjf = epsilS(jj)*Sr1f+etaS(jj)*Sr2f+xiS(jj)*Sr3f

      Einc=E0*exp(cmplx(0.0,-dot_product(vec_k,rjf)))
	 
      V_E1 = V_E1+weightS(jj)*dot_product(rjf,Einc)
      V_E2 = V_E2+weightS(jj)*Einc
    enddo

  endsubroutine



  subroutine Vs_filling_term_by_PW_EH(j,Sr1f,Sr2f,Sr3f,nf,V_E1,V_E2,V_H1,V_H2)
    use Dunavant_integration_formulism
    implicit none
    integer,intent(in):: j
    real,intent(in):: Sr1f(3),Sr2f(3),Sr3f(3),nf(3)
    complex,intent(out):: V_E1,V_E2(3),V_H1,V_H2(3)
    
    real rjf(3),rjf_nf(3)
    
    integer jj
    complex Einc(3),Hinc(3)


    V_E1=(0.0,0.0)
    V_E2=(0.0,0.0)
    V_H1=(0.0,0.0)
    V_H2=(0.0,0.0)


	   do jj=1,SNum_samp
      rjf = epsilS(jj)*Sr1f+etaS(jj)*Sr2f+xiS(jj)*Sr3f
	   
      Einc=E0*exp(cmplx(0.0,-dot_product(vec_k,rjf)))
      Hinc=H0*exp(cmplx(0.0,-dot_product(vec_k,rjf)))
	  
      call Re_cross_product(rjf,nf,rjf_nf)
      
      V_E1 = V_E1+weightS(jj)*dot_product(rjf,Einc)
      V_E2 = V_E2+weightS(jj)*Einc       
	     V_H1 = V_H1+weightS(jj)*dot_product(rjf_nf,Hinc)
      V_H2 = V_H2+weightS(jj)*Hinc
    enddo

  endsubroutine



end