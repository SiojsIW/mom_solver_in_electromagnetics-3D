subroutine Load_triangle_tetrahedron_mesh
  use Common_array
  implicit none

  integer i,j,k,Num_rows,itmp
  character char
  integer,external:: GetFileN_nas
  
  write(*,*) 'To load composite mesh.'

  open(10,file=nas_name)

    do j=1,9   
      if(j==7) then
        read(10,*) char,char,Num_triangle
      elseif(j==9) then
        read(10,*) char,char,char
      else
        read(10,*) char
      endif
    enddo

    
    itmp=GetFileN_nas(10)
    Num_point=(itmp-9)/2

    Num_rows=Num_point*2+Num_triangle+10
    
  
    
    allocate(point_cor(3,Num_point))
    if(Num_triangle>0) allocate(triangle_point(3,Num_triangle))
    
    write(*,*) 'Number of triangles: ',Num_triangle
    write(*,*) 'Number of nodes: ',Num_point
    write(*,*) 'rows in *.nas : ',Num_rows
  
    write(*,*) 'To read geometric info, Waiting...'

    do j=1,9
      read(10,*) char
    enddo
    
    if(Num_point<=9999999) then
      do j=10,Num_point*2+9  !read node data
        if(mod(j-9,2)==1) then
          i=(j-9)/2+1
          read(10,100) char,k,point_cor(1,i),point_cor(2,i),k
        else
          i=(j-9)/2
          read(10,101) char,k,point_cor(3,i)
        endif
  100   format(A5,I19,ES32.9,ES16.9,I8)
  101   format(A1,I7,ES16.9)
      enddo

    else
      do j=10,9999999*2+9  !read node data
        if(mod(j-9,2)==1) then
          i=(j-9)/2+1
          read(10,1002) char,k,point_cor(1,i),point_cor(2,i),k
        else
          i=(j-9)/2
          read(10,1012) char,k,point_cor(3,i)
        endif
 1002   format(A5,I19,ES32.9,ES16.9,I8)
 1012   format(A1,I7,ES16.9)
      enddo
    
      do j=9999999*2+10,Num_point*2+9  !read node data
        if(mod(j-9,2)==1) then
          i=(j-9)/2+1
          read(10,1001) char,k,point_cor(1,i),point_cor(2,i),k
        else
          i=(j-9)/2
          read(10,1011) char,k,point_cor(3,i)
        endif
 1001   format(A5,I19,ES32.9,ES16.9,I8)
 1011   format(A1,I8,ES16.9)
      enddo

    endif

    do j=Num_point*2+9+1,Num_point*2+Num_triangle+9  !read triangle_point data
      i=j-Num_point*2-9
      read(10,102) char,k,k,triangle_point(1,i),triangle_point(2,i),triangle_point(3,i)
102   format(A6,I10,I8,I8,I8,I8)
    enddo

  close(10)
      

  write(*,*) 'Load mesh completely.'
  write(*,"(/,'')")
 
end


integer function GetFileN_nas(iFileUnit)
	 implicit none
	 integer:: iFileUnit
	 integer:: ioS
	 character:: char

	 GetFileN_nas=0
	 rewind(iFileUnit)
	   do while(.true.)
		    read(iFileUnit,*,ioStat=ioS) char
		    if (char=='C') exit
		    GetFileN_nas=GetFileN_nas+1
	   enddo
	 rewind(iFileUnit)
end