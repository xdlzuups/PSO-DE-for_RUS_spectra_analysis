module PSO !Particle Swarm Optimization
real:: w=0.9, c1=2, c2=2
integer:: nParticle = 80
contains
 
subroutine particleSwarmOptimization(w, c1, c2, nParticle, rang, nIts, nLineF, nFre, obs, weight, err0, res, par0)
use RUS, only: calFrequence
!$ use omp_lib
implicit none
real, intent(in):: w, c1, c2 
integer, intent(in):: nParticle, nIts, nLineF, nFre
real,    intent(in):: rang(:,:), obs(:), weight(:), err0
real, intent(out)::   res(:)
real, intent(in), optional:: par0(:)
real, allocatable:: x(:,:), pBest(:,:), gBest(:), minErr(:), fre(:)
real, allocatable:: rand1(:), rand2(:), v(:,:), vmax(:), vLimit(:), err(:)
integer nPar
integer its, i

nPar = size(rang,2) !参数个数
allocate(x(nPar,nParticle), v(nPar,nParticle), pBest(nPar,nParticle), gBest(nPar)) !位置、速度, 最佳位置
allocate(rand1(nParticle), rand2(nParticle))   !随机数
allocate(minErr(nParticle),vmax(nPar),vLimit(nPar),fre(nFre),err(nParticle))

!初始位置及速度
call RANDOM_NUMBER(x); call RANDOM_NUMBER(v(1,1)); call RANDOM_NUMBER(v)
vmax(:) = (rang(2,:)-rang(1,:)) !各个参数的扰动幅度
do i = 1, nParticle
  x(:,i) = rang(1,:) + x(:,i)*vmax(:)
  v(:,i) = (v(:,i)-0.5) * vmax(:)
end do
if(present(par0)) x(:,1) = par0(:)

!初始化误差
minErr = huge(1.0)
gBest = x(:,1)
!迭代
do its = 1, nIts
  write(30,"(a,g0)") 'its = ', its
  !速度更新量递减
  vLimit = vmax * (1.0-(its-1.0)/nIts) * 0.25
  call RANDOM_NUMBER(rand1); call RANDOM_NUMBER(rand2(1)); call RANDOM_NUMBER(rand2)
  do i = 1, nParticle
    !误差
    call calFrequence(x(:,i), fre)
    err(i) = fitnessFun(fre, obs, weight, nLineF, .false.)   
    !更新个体最优
    if(err(i)<minErr(i)) then
      pBest(:,i) = x(:,i)
      minErr(i) = err(i)
    end if 
    !更新速度
    v(:,i) = (1.0-0.5*(its-1.0)/nIts)*w*v(:,i) + c1*rand1(i)*(pBest(:,i)-x(:,i)) + c2*rand2(i)*(gBest(:)-x(:,i))
    v(:,i) = v(:,i) / maxval(abs(v(:,i)/vLimit))
    !where(abs(v(:,i))>vLimit) v(:,i) = sign(vLimit, v(:,i))
  end do
  do i = 1, nParticle
    !输出粒子轨迹, 历史最优
    write(30,"(*(g0,','))") i, pBest(:,i), minErr(i), x(:,i), err(i)
    !新位置
    x(:,i) = x(:,i) + v(:,i)
    where(x(:,i)<rang(1,:)) x(:,i) = x(:,i) + rang(2,:) - rang(1,:)
    where(x(:,i)>rang(2,:)) x(:,i) = x(:,i) - rang(2,:) + rang(1,:) 
    !x(:,i) = min(rang(2,:),x(:,i)); x(:,i) = max(rang(1,:),x(:,i))
  end do
  gBest = pBest(:,minloc(minErr,1))
  !输出全局最优
  write(30,"('gBest,',*(g0,','))") gBest(:), minval(minErr)
  write(*,"('its = ',i3,'/',g0,', minErr = ',g0)") its, nIts, minval(minErr)
  !if(minval(minErr)<err0) exit !满足误差，退出迭代
end do

!结果
res = gBest
end subroutine

!误差函数
function fitnessFun(fre, obs, weight, nLineF, isOutput) result(res)
implicit none
real, parameter:: EPS = 1.0d-6, HUG = huge(1.0)/10.0
real, intent(in):: fre(:), obs(:), weight(:)
integer,intent(in):: nLineF  
logical,intent(in):: isOutput 
real err(size(obs)), fittedFre(size(obs))
real res, e
integer i, j, k, i1
!计算百分比误差
err = huge(1.0)
k = 1 !检索起始位置，obs必须升序排列
do j = 1, size(obs)
  if(weight(j)<EPS .or. obs(j)<EPS) cycle
  i1 = k
  do i = i1, size(fre)
    e = abs(obs(j)-fre(i))
    if(e<err(j)) then
      err(j) = e
      k = i
    else if(e>err(j)+EPS) then
      exit
    end if
  end do
  fittedFre(j) = fre(k)
  err(j) = err(j) / obs(j) * weight(j) * 100.0
end do
do i = nLineF+1, count(err<HUG)
	j = maxloc(err,1,err<HUG)
  err(j) = huge(1.0)
end do
res = 0; e = 0
do i = 1, size(err)
  if(err(i)>HUG) cycle
  res = res + abs(err(i))
  e = e + weight(i)
end do
res = res / e

if(.not.isOutput) return
!输出
!open(20,file='fitted_Fre.csv')
write(30,*)  
write(30,'("obsFre,calFre,err(%)")')
do i = 1, size(err)
	if(err(i)>HUG) cycle
  write(30,"(*(g0,','))") obs(i),fittedFre(i),err(i) 
end do
end function
end module