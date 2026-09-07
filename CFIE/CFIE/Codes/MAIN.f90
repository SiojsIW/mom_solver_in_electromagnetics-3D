program main
  use Common_array
  use Common_constants
	 implicit none
  real tarry(2),et,etime

  call Input_setting(1)
  call Load_triangle_tetrahedron_mesh
  call preprocessing_triangle   
  call Input_setting(2)
  call MoM_SS
  
  call excitation_filling
  allocate(Inn(Num_side))
  call bicgstab(Num_side,Inn,excitation,tol,itMax)
  call scattering_field
  
  et=etime(tarry)
  
	 write(*,*) "Mission Completed!!! ºÄÊ±",et/60.0,"min"
  pause
  
end