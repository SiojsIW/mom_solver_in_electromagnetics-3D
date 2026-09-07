module Common_array

  character*256 nas_name
  ! mesh node
  integer Num_point
  real,allocatable:: point_cor(:,:)

  ! RWG
  integer Num_PEC_patch
  integer Num_triangle,Num_side
  
  integer*4,allocatable:: triangle_point(:,:),side_index(:,:)
  real,allocatable:: triangle_centroid(:,:),side_length(:),triangle_area(:)
  
  complex,allocatable:: SparseValue(:,:)
  
  
  type::triangle_associated_rwg
    integer*2 N_rwg
    integer*4,allocatable::rwg(:)
  endtype
  type(triangle_associated_rwg),allocatable:: triangle_rwgs(:)
  
  complex,allocatable:: Inn(:),excitation(:)
end module





module Common_constants
  integer*2,parameter:: maxRWG=6
  
  real,parameter:: PI=4*atan(1.0),MU0=4*PI*1.0e-7,EPS0=1.0/(36*PI)*1.0e-9,ETA0=120*PI
  complex,parameter:: iu=(0.0,1.0)
  logical:: CFSIE=.TRUE.
  real,parameter:: alpha=0.5
  integer itMax
  real tol
  
  real freq,omega,lambda,k_wnum
  integer tot_freq
  real threshold_R,ratio_limit1
  real,parameter:: ratio_limit2=3.0e-2
end module



module Far_patterns
  integer*2,parameter:: maxNumscan=361
  
  integer Num_RCS
  integer Num_scat_phi,Num_scat_theta
  integer*1 Is_phi,Is_theta
  complex E0(3),H0(3)
  real vec_k(3)
  
  real theta_inc,phi_inc,E_inc_0_theta,E_inc_0_phi,E_inc_0_theta_pha,E_inc_0_phi_pha
  
  integer Num_angles_phi(maxNumscan),Num_angles_theta(maxNumscan)
  
  real phi_scat_theta_start(maxNumscan),phi_scat_theta_end(maxNumscan),theta_scat_phi_start(maxNumscan),theta_scat_phi_end(maxNumscan)
  real phi_scats(maxNumscan),theta_scats(maxNumscan)
end module



module Dunavant_integration_formulism
  integer,parameter:: SNum_samp=4
  real,parameter:: weightS(SNum_samp)=(/-0.5625,1.5625/3.0,1.5625/3.0,1.5625/3.0/)
  real,parameter:: epsilS(SNum_samp)=(/1.0/3.0,0.6,0.2,0.2/)
  real,parameter:: etaS(SNum_samp)=(/1.0/3.0,0.2,0.6,0.2/)
  real,parameter:: xiS(SNum_samp)=(/1.0/3.0,0.2,0.2,0.6/)
end module