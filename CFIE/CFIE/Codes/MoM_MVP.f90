subroutine MVP_MoM(X,b)
  use Common_array
  implicit none
  complex X(*)
  complex b(*)
  
  integer i,j
  complex ctmp
  
  !$OMP parallel do default(none) shared(Num_side,SparseValue,X,b) private(i,ctmp)
  do i=1,Num_side
    ctmp=0
    do j=1,Num_side
      ctmp=ctmp+SparseValue(i,j)*X(j)
    enddo
    b(i)=ctmp
  enddo  
  
end
  
  