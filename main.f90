program ParticleSwarm 
use RUS
use PSO
use Input
use ifport, only: rtc
use fit_mod
use Differential_Evolution
implicit none
real, parameter:: EPS = 1.0e-5
real d1, d2, d3      !样品尺寸(由于对称性，为实际值的一半)
integer:: medimuType !介质类型
integer:: shape=1    !形状: 1圆柱体，2球，3长方体
integer:: d          !拟合多项式的阶数
real mass            !质量
real err0            !截断误差, 终止迭代条件
integer nIts         !最大迭代次数 
integer nLineF       !用于拟合的频线数
real, allocatable:: par0(:)  !估计参数
real, allocatable:: obsFre(:), weight(:) !观测频率及权重
real, allocatable:: par(:), fre(:), rang0(:,:), rang(:,:)
real(8) time0, time1 !计时
integer nPar
real err
integer i
time0 = rtc()
call RANDOM_SEED()

!读取输入文件
call readFile(d, d1, d2, d3, mass, err0, medimuType, nIts, nLineF, par0, obsFre, weight)
nPar = size(par0)
allocate(par(nPar))
par = par0

!初始化: 分配内存, 计算索引表、e
call init(d, d1, d2, d3, mass, shape, medimuType, fre)

!正演
if(nIts<=1) then
  call calFrequence(par0, fre) !计算频线
  open(10,file='fre.txt')
  do i = 1, size(fre)
    if(fre(i)>EPS) write(10,*) fre(i)
  end do
  close(10)
  write(*,"('forward DONE...')")
  read(*,*)
  stop
end if

!读取反演参数
allocate(rang0(2,nPar),rang(2,nPar))
!call readInvPar(nParticle, w, c1, c2, rang0)
call readDEPar(nParticle, F, CR, rang0)
rang = rang0

f_obs = obsFre
f_weight = weight
f_nLineF = nLineF
f_nFre = size(fre)
call DE(F, CR ,nParticle ,nIts, fit, rang0, err0, par, par0)



!最佳拟合的参数
call calFrequence(par, fre)
err = fitnessFun(fre, obsFre, weight, nLineF, .true.)

time1 = rtc()
write(*,"('CPU Time: ', g0, 's')") real(time1-time0,4)
write(30,*)  
write(30,"('CPU Time:, ', g0, 's')") real(time1-time0,4)

close(30)
read(*,*)  
end program ParticleSwarm