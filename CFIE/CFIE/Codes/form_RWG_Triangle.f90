subroutine form_RWG_Triangle
  use Common_array
  implicit none
  integer i,j,k
  integer node1_index,node2_index
  integer*8 RWG_index,key,temp
  integer last,key1,temp1
  integer*4,allocatable:: sidetemp(:,:)
  integer*8,allocatable:: heap(:)
  integer*4,allocatable:: heap1(:)
  integer*8,parameter:: maxnum=2**31-1
   
  last=0
  allocate(heap(3*Num_triangle),heap1(3*Num_triangle))
  
  
  do i=1,Num_triangle 
    do j=1,3
      k=j+1
      if(k==4) k=1
      
      node1_index=triangle_point(j,i)
      node2_index=triangle_point(k,i)  
      if(node1_index<node2_index) then
        RWG_index=node1_index*maxnum+node2_index     
      else
        RWG_index=node2_index*maxnum+node1_index
      endif   
    
      call insert_int_heap1(3*Num_triangle,i,RWG_index,last,heap,heap1)
      
    enddo  
  enddo
  
  
  allocate(sidetemp(4,3*Num_triangle))
  
  
  Num_side=0
  
  temp=-1
 
  do i=1,3*Num_triangle

    call int_heap_get_first1(3*Num_triangle,key,key1,last,heap,heap1)
    
    if(temp/=key) then
      temp=key
      temp1=key1
    else
      
      Num_side=Num_side+1
    
      sidetemp(1,Num_side)=key/maxnum
      sidetemp(2,Num_side)=mod(key,maxnum)  
      sidetemp(3,Num_side)=temp1
      sidetemp(4,Num_side)=key1     
    endif
   
  enddo

  
  deallocate(heap,heap1)
  
  allocate(side_index(4,Num_side))
  side_index=sidetemp(:,1:Num_side)
  deallocate(sidetemp)
  
  
end


subroutine int_heap_get_first1(n,key,key1,last,heap,heap1)
  implicit none
  integer n
  integer*8 key
  integer key1
  integer last
  integer*8 heap(n)
  integer*4 heap1(n)

  integer i,j
  integer*8 temp

  key     = heap(1)
  heap(1) = heap(last)
  
  key1     = heap1(1)
  heap1(1) = heap1(last)
  
  last    = last - 1


  i = 1
  do 
    if(i>(last/2)) exit
    if ((heap(2*i)<heap(2*i+1)).or.(2*i == last)) then 
      j = 2*i
    else
      j = 2*i + 1
    end if

    if(heap(i)>heap(j)) then 
      temp =heap(i)
      heap(i)=heap(j)
      heap(j)=temp
      
      temp =heap1(i)
      heap1(i)=heap1(j)
      heap1(j)=temp
      i=j 
    else
      exit
    endif
  enddo

end



subroutine insert_int_heap1(n,j,key,last,heap,heap1)
  implicit none 
  integer n,j
  integer*8 key
  integer last
  integer*8 heap(n)
  integer*4 heap1(n)
  
  integer i, i2
  integer*8 temp

  last=last+1

  i=last
  heap(i)=key
  heap1(i)=j
  
  do
    if(i<=1) exit
    i2 = i/2
    if (heap(i) < heap(i2)) then 
      temp     = heap(i)
      heap(i)  = heap(i2)
      heap(i2) = temp
      
      temp     = heap1(i)
      heap1(i)  = heap1(i2)
      heap1(i2) = temp     
      
      i        = i2
    else
      exit
    endif
  enddo

end