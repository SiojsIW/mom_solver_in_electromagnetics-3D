subroutine Input_setting(flag)
	 use Common_array
	 use Common_constants
	 use Dunavant_integration_formulism
  use Far_patterns
	 implicit none
  integer*1 flag
  
  integer*1 Is
  integer i,j
  integer filein 
  
  character*256 filename
  logical alive
  
  if(flag==1) then 
567 write(*,*) "Input filein, 5-Keyboard  7-text"
    read(*,*) filein

    
    if(filein/=5) open(filein,file='Setting.txt')
  !************************************
    write(*,*) 'Nastran file name input:'
    read(filein,*) nas_name

    filename='Meshes\'//trim(adjustl(nas_name))//'.nas'
    inquire(file=filename,exist=alive)
    if(alive) then
      write(*,*) '*.nas: ', trim(adjustl(nas_name)),' is to be loaded.'
      nas_name=filename
    else
      write(*,*) '*.nas: ', trim(adjustl(nas_name)),' does not exist!'
      stop
    endif
    write(*,"(/,'')")
  !************************************
     
      write(*,*) "Input the frequency(MHz)"  
      read(filein,*) freq
      freq=freq*1.0e6   

      write(*,"(/,'Input the amplitude of incident E_theta,E_phi',/)")
      write(*,"('E_theta',/)")
      read(filein,*) E_inc_0_theta
      write(*,"('E_phi',/)")
      read(filein,*) E_inc_0_phi

 
      E_inc_0_theta_pha=0
      E_inc_0_phi_pha=0
      
      
	     write(*,"(/,'Input the direction of the wave propagation',/)")
	     write(*,"(/,'Angle with respect to Z axis: theta_inc (Degree)',/)") 
	     read(filein,*) theta_inc
	     write(*,"(/,'Angle with respect to X axis: phi_inc (Degree)',/)") 
	     read(filein,*) phi_inc
      
      write(*,"(/,'The observatory planes for the Phi? 1-Yes, 0-No',/)")
      read(filein,*) Is_phi
        
        
      write(*,"(/,'The observatory planes for the Theta? 1-Yes, 0-No',/)")
      read(filein,*) Is_theta
        
        
      if(Is_phi==1) then
        write(*,"(/,'How many planes do you want to observe? (Phi)',/)")
        read(filein,*) Num_scat_phi
          
        do i=1,Num_scat_phi
         write(*,*) "The",i,"Phi (degree)"
          read(filein,*) phi_scats(i)
        enddo
      endif
        
      if(Is_theta==1) then
        write(*,"(/,'How many planes do you want to observe? (Theta)',/)")
        read(filein,*) Num_scat_theta
         
        do i=1,Num_scat_theta
          write(*,*) "The",i,"Theta (degree)"
          read(filein,*) theta_scats(i)
        enddo
      endif

        
      write(*,"(/,'RCS Pattern sampling Settings:')")
      write(*,"('enter sampling points per half slice: e.g. 181')")
      read(filein,*) Num_RCS
 
        
    write(*,"(/,'Maximum iterations:  ')") 
    read(filein,*) itMax
      
    write(*,"('Convergence within:  ')") 
    read(filein,*) tol
    write(*,*) "内存用量约为：",Num_side*(Num_side+100)*8/1024.0/1024.0,"MB"
    
    if(filein/=5) close(filein)
    
  elseif(flag==2) then
  

    
    omega=2.0*PI*freq

    lambda=3.0e8/freq
	   k_wnum=2.0*PI/lambda
  
	   threshold_R=0.1*lambda
    
    ratio_limit1=1.0e-6*lambda       !for w0 in the singularity computation，这一值的大小将会影响MFIE的精度，需要仔细设置
  endif
end