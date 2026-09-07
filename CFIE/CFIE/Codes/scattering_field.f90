subroutine scattering_field
  use Common_constants,only: PI
  use Far_patterns
	 implicit none
  
	 integer i,j
	 real theta,phi,delta,ratio
  real P
  real E0_far(2*Num_RCS-1),E0_theta_far(2*Num_RCS-1),E0_phi_far(2*Num_RCS-1)
  complex E0_theta,E0_phi
  character*4 char
  
  write(*,*) "��ʼ����Զ��"
  
  P=4*PI*(E_inc_0_theta**2+E_inc_0_phi**2)
  
  
  delta=PI/(Num_RCS-1)
  
  if(Is_phi==1) then
    do i=1,Num_scat_phi
      write( char,'(i4)') i
      phi=phi_scats(i)*PI/180.0    

     !$OMP parallel do default(none) shared(Num_RCS,delta,phi,P,E0_theta_far,E0_phi_far,E0_far) &
     !$OMP  private(j,theta,E0_theta,E0_phi,ratio)
      do j=1,2*Num_RCS-1
	       theta=delta*(j-1)
        call E0_at_theta_phi(theta,phi,E0_theta,E0_phi)
        E0_theta_far(j)=abs(E0_theta)**2/P
        E0_phi_far(j)=abs(E0_phi)**2/P
        E0_far(j)=E0_theta_far(j)+E0_phi_far(j)

      enddo
     !$OMP end parallel do
  
       open (11,file = 'Output/FixedPhi'//trim(adjustl(char))//'_RCS_db.txt')
   
        do j=1,2*Num_RCS-1
          theta=delta*(j-1)
          write(11,*) theta*180.0/PI,10.0*log10(E0_theta_far(j)),10.0*log10(E0_phi_far(j)),10.0*log10(E0_far(j))
        enddo

      close(11)
    enddo
  endif

  
  if(Is_theta==1) then
    
    do i=1,Num_scat_theta
      write( char,'(i4)') i
      theta=theta_scats(i)*PI/180.0
   
     !$OMP parallel do default(none) shared(Num_RCS,delta,theta,P,E0_theta_far,E0_phi_far,E0_far) &
     !$OMP  private(j,phi,E0_theta,E0_phi,ratio)
      do j=1,2*Num_RCS-1
	        phi=delta*(j-1)
         call E0_at_theta_phi(theta,phi,E0_theta,E0_phi)
         E0_theta_far(j)=abs(E0_theta)**2/P
         E0_phi_far(j)=abs(E0_phi)**2/P
         E0_far(j)=E0_theta_far(j)+E0_phi_far(j)

       enddo
     !$OMP end parallel do

      
       open(11,file='Output/FixedTheta'//trim(adjustl(char))//'_RCS_db.txt')
      
        do j=1,2*Num_RCS-1
          phi=delta*(j-1)
          write(12,*) phi*180.0/PI,10.0*log10(E0_theta_far(j)),10.0*log10(E0_phi_far(j)),10.0*log10(E0_far(j))
        enddo 

      close(11)
    enddo
  endif

  
  write(*,*) "���"

contains
 

  subroutine E0_at_theta_phi(theta,phi,E0_theta,E0_phi)
    use Common_array
    use Common_constants
    implicit none
    real,intent(in):: theta,phi
    complex,intent(out):: E0_theta,E0_phi
  
    integer i,lm,li,sideindex,faceindex
    integer node1,node2,node
    real(4) E0_theta_amp,E0_phi_amp
    real(4) cartesian_to_theta(3),cartesian_to_phi(3),r_hat(3)
    real(4) Srn(3),Sc(3),rouSc(3)
    real(4) length
    complex*8 phasor,exptmp
  
    cartesian_to_theta=(/cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)/)
    cartesian_to_phi=(/-sin(phi),cos(phi),0.0/)
    r_hat=(/sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)/)

  
    E0_theta=(0.0,0.0)
    E0_phi=(0.0,0.0)
  
    do i=1,Num_Triangle
      Sc=triangle_centroid(:,i)
      exptmp=exp(cmplx(0.0,k_wnum*dot_product(Sc,r_hat)))
   

      do li=1,triangle_rwgs(i)%N_rwg
        sideindex=triangle_rwgs(i)%rwg(li)
      
        if(side_index(3,sideindex)==i) then
	         length=side_length(sideindex)
        else
	         length=-side_length(sideindex)
	       endif
      
	       node1=side_index(1,sideindex); node2=side_index(2,sideindex)
        do lm=1,3
          node=triangle_point(lm,i)
	         if((node/=node1).and.(node/=node2)) exit
        enddo
        Srn=point_cor(:,node)

        rouSc=Sc-Srn

	       E0_theta_amp=length*dot_product(rouSc,cartesian_to_theta)
	       E0_phi_amp=length*dot_product(rouSc,cartesian_to_phi) 
	       phasor=Inn(sideindex)*exptmp
      
	       E0_theta=E0_theta+E0_theta_amp*phasor
	       E0_phi=E0_phi+E0_phi_amp*phasor 
	     enddo
	   enddo

 
    E0_theta=cmplx(0.0,-omega*MU0)*E0_theta/2.0
    E0_phi=cmplx(0.0,-omega*MU0)*E0_phi/2.0
  
  endsubroutine

end