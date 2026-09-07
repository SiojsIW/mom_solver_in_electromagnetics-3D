subroutine bicgstab(N,x,b,tolerance,itMax)
  implicit none
  integer,intent(in):: N
  complex,intent(out):: x(N)
  complex,intent(in):: b(N)
  integer,intent(in):: itMax
  real,intent(in):: tolerance
  
  integer k
  real errtol,norm_b,zeta
  complex omega,beta,alpha,tau
  complex,allocatable:: v(:),p(:),s(:),t(:),r(:),hatr0(:),rho(:)
  real rel_res(itMax)

  real tarry(2),et(2),etime
  
  write(*,*) "采用BiCGstab求解器进行迭代求解，计算中..."
  write(*,*)  "最大迭代步数：",itMax,"，收敛值：",tolerance
 
  et(1)=etime(tarry)
  
  x=0

  norm_b=sqrt(dot_product(b,b))
  errtol=tolerance*norm_b

 
  allocate(rho(itMax+2))

  allocate(r(N))
  if(dot_product(x,x)/=0.0) then
    call MVP_MoM(x,r)
    r=b-r
  else
    r=b
  endif

  
  allocate(hatr0(N))
  hatr0=r

  k=0

  rho(1)=1.0
  alpha=1.0
  omega=1.0


  allocate(v(N),p(N),t(N),s(N))
  v=(0.0,0.0)
  p=(0.0,0.0)
  t=(0.0,0.0)
  s=(0.0,0.0)
  
  rho(2)=dot_product(hatr0,r)

  zeta=sqrt(dot_product(r,r))

  do while((zeta>errtol).and.(k<itMax))
    k=k+1
    if(omega==0.0) then
      write(*,*) "Bi-CGSTAB求解器崩溃, 1"
      return
    endif

    beta=(rho(k+1)*alpha)/(rho(k)*omega)


    p=r+beta*(p-omega*v)

    
    call MVP_MoM(p,v)


    tau=dot_product(hatr0,v)
    

    if(tau==0.0) then
      write(*,*) "***错误***：Bi-CGSTAB求解器崩溃！2"
      return
    endif 

    alpha=rho(k+1)/tau 
    
    
    s=r-alpha*v
    

    call MVP_MoM(s,t)
    
    tau=dot_product(t,t)


    if(tau==0.0) then
      write(*,*) "***错误***：Bicgstab求解器崩溃！3"
      return
    endif

    omega=dot_product(t,s)/tau
    
    rho(k+2)=-omega*dot_product(hatr0,t)

    x=x+alpha*p+omega*s

    r=s-omega*t
    zeta=sqrt(dot_product(r,r))
    
    rel_res(k)=zeta/norm_b
    
    write(*,*) "第 ",k," 步，残差：",rel_res(k)
    
  enddo

  
  if(k==0) then

  elseif(1<=k.and.k<itMax) then
    write(*,*) "已收敛，共迭代 ",k," 步，最终残差：",rel_res(k)
  elseif(k==itMax) then
    write(*,*) "**警告**：Bicgstab已达到最大迭代步数，但仍未能收敛到真解"
  endif

  deallocate(r,v,p,hatr0,s,t,rho)
    
  et(2)=etime(tarry)
  
  write(*,*) "Bicgstab迭代求解完毕。耗时：",(et(2)-et(1))/60.0," min"
  
end