subroutine preprocessing_triangle
  use Common_array
  use Common_constants,only: CFSIE
  implicit none
  
  integer i,j,m
  integer numtemp,Ntemp3(3),T1,T2,count   
  real Vec1(3),Vec2(3),CrossVector(3)
  real,external::re_Norm
  
  
  call form_RWG_Triangle

  allocate(side_length(Num_side))
  do i=1,Num_side
    side_length(i)=re_Norm(point_cor(:,side_index(1,i))-point_cor(:,side_index(2,i)))
  enddo  
  
  allocate(triangle_area(Num_triangle),triangle_centroid(3,Num_triangle))

  do i=1,Num_triangle
    Ntemp3=triangle_point(:,i)
        
    Vec1=point_cor(:,Ntemp3(2))-point_cor(:,Ntemp3(1))
    Vec2=point_cor(:,Ntemp3(3))-point_cor(:,Ntemp3(1))     
  
    call Re_cross_product(Vec1,Vec2,CrossVector)
      
    triangle_area(i)=re_Norm(CrossVector)/2
        
    triangle_centroid(1,i)=(point_cor(1,Ntemp3(1))+point_cor(1,Ntemp3(2))+point_cor(1,Ntemp3(3)))/3
    triangle_centroid(2,i)=(point_cor(2,Ntemp3(1))+point_cor(2,Ntemp3(2))+point_cor(2,Ntemp3(3)))/3
    triangle_centroid(3,i)=(point_cor(3,Ntemp3(1))+point_cor(3,Ntemp3(2))+point_cor(3,Ntemp3(3)))/3
  enddo

    
  allocate(triangle_rwgs(0:Num_triangle))

  do i=0,Num_triangle
    triangle_rwgs(i)%N_rwg=0
  enddo
    
    
  m=0
  do j=1,Num_side
    T1=side_index(3,j)
    T2=side_index(4,j)

    if(T1/=0) then
      triangle_rwgs(T1)%N_rwg= triangle_rwgs(T1)%N_rwg +1
      if(triangle_rwgs(T1)%N_rwg>m) m=triangle_rwgs(T1)%N_rwg
    endif
      
    if(T2/=0) then
      triangle_rwgs(T2)%N_rwg= triangle_rwgs(T2)%N_rwg +1
      if(triangle_rwgs(T2)%N_rwg>m) m=triangle_rwgs(T2)%N_rwg
    endif
  enddo
    

  do i=1,Num_triangle
    allocate(triangle_rwgs(i)%rwg(triangle_rwgs(i)%N_rwg))
    triangle_rwgs(i)%rwg=0
  enddo
      

  do j=1,Num_side
    T1=side_index(3,j)
    T2=side_index(4,j)

    do m=1,triangle_rwgs(T1)%N_rwg
      if(triangle_rwgs(T1)%rwg(m)==0) then
        triangle_rwgs(T1)%rwg(m) =j
        exit
      endif
    enddo

    do m=1,triangle_rwgs(T2)%N_rwg
      if(triangle_rwgs(T2)%rwg(m)==0) then
        triangle_rwgs(T2)%rwg(m)=j
        exit
      endif
    enddo
  enddo

  if(3*Num_triangle==2*Num_side) then
    CFSIE=.TRUE.
  else
    CFSIE=.FALSE.
  endif
  
end


subroutine Re_cross_product(a,b,c)
  implicit none
  real a(3),b(3),c(3)
  c=(/(a(2)*b(3)-a(3)*b(2)),(a(3)*b(1)-a(1)*b(3)),(a(1)*b(2)-a(2)*b(1))/)
end
 
real function re_Norm(x)
  implicit none
  real x(3)
  re_Norm=sqrt(dot_product((x),(x)))
end