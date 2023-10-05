module Input
contains
subroutine readFile(d, d1, d2, d3, mass, convergencePar, medimuType, nIts, nLineF, par, fre, weight)
implicit none
integer nPar                !参数个数
integer d                   !多项式阶数
integer nLineC              !计算的频线数
real    mass                !质量
real    convergencePar      !收敛参数
integer nIts                !迭代次数
integer nLineF              !用于拟合的频线数
real    poison              !估计泊松比
real, allocatable:: par(:)  !估计参数
real    d1, d2, d3          !样品尺寸
real, allocatable:: fre(:), weight(:)  !测量所得频率及权重
integer medimuType !介质类型
real fre1, fre2, w
integer i 

!读取输入文件
open(10,file='cylin.dat')
read(10,*)
read(10,*) nPar, d, nLineC, mass, convergencePar, nIts, nLineF, poison
allocate(par(nPar))
read(10,*) par(1:nPar) !估计参数
read(10,*) d1, d2, d3  !样品尺寸
!读取
do 
  read(10,*,iostat=i) fre1, fre2, w
  if(i/=0) exit
  fre = [fre,fre1]
  weight = [weight,w]
end do
close(10)

nLineF = min(nLineF,count(weight>1.0e-5))

!确定介质类型
select case(nPar)
case(2)
medimuType = 1
case(3)
medimuType = 2
case(5)
medimuType = 3
case(6)
medimuType = 5
case(9)
medimuType = 6
end select

nLineC = max(nLineF+4, nLineC)

!由于对称性，样品尺寸减半后计算
d1 = 0.5 * d1; d2 = 0.5 * d2; d3 = 0.5 * d3
end subroutine

!读取反演信息
subroutine readInvPar(nParticle, w, c1, c2, rang)
implicit none
integer nParticle
real w, c1, c2, rang(:,:)
open(10,file='InvPar.txt',status='old')
read(10,*); read(10,*) nParticle; write(*,"('nParticle: ', g0)")  nParticle
read(10,*); read(10,*) w, c1, c2; write(*,"('w, c1, c2: ', *(g0,2x))")  w, c1, c2
read(10,*); read(10,*) rang; 
!write(*,"('rang: ')")  
!write(*,"(2(g0,2x))")  rang
close(10)
end subroutine
end module