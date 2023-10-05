program ParticleSwarm 
use RUS
use PSO
use Input
use ifport, only: rtc
implicit none
real, parameter:: EPS = 1.0e-5
real d1, d2, d3      !��Ʒ�ߴ�(���ڶԳ��ԣ�Ϊʵ��ֵ��һ��)
integer:: medimuType !��������
integer:: shape=1    !��״: 1Բ���壬2��3������
integer:: d          !��϶���ʽ�Ľ���
real mass            !����
real err0            !�ض����, ��ֹ��������
integer nIts         !���������� 
integer nLineF       !������ϵ�Ƶ����
real, allocatable:: par0(:)  !���Ʋ���
real, allocatable:: obsFre(:), weight(:) !�۲�Ƶ�ʼ�Ȩ��
real, allocatable:: par(:), fre(:), rang0(:,:), rang(:,:)
real(8) time0, time1 !��ʱ
integer nPar
real err
integer i
time0 = rtc()
call RANDOM_SEED()

!��ȡ�����ļ�
call readFile(d, d1, d2, d3, mass, err0, medimuType, nIts, nLineF, par0, obsFre, weight)
nPar = size(par0)
allocate(par(nPar))

!��ʼ��: �����ڴ�, ����������e
call init(d, d1, d2, d3, mass, shape, medimuType, fre)

!����
if(nIts<=1) then
  call calFrequence(par0, fre) !����Ƶ��
  open(10,file='fre.txt')
  do i = 1, size(fre)
    if(fre(i)>EPS) write(10,*) fre(i)
  end do
  close(10)
  write(*,"('forward DONE...')")
  read(*,*)
  stop
end if

!��ȡ���ݲ���
allocate(rang0(2,nPar),rang(2,nPar))
call readInvPar(nParticle, w, c1, c2, rang0)
rang = rang0

open(30,file='Result.csv')
select case(nPar)
case(2)
!write(10,"(a)") 'its time,minErr,initial C11,initial C12,final C11,final C12,N of lost peaks'
write(30,"('no,C11,C12,err(%),C11,C12,err(%)')") 
write(30,"(',pBest,,,position')")  
case(3)
!write(10,"(a)") 'its time,minErr,initial C11,initial C12,initial C44,final C11,final C12,final C44,N of lost peaks'
write(30,"('no,C11,C12,C44,err(%),C11,C12,C44,err(%)')") 
write(30,"(',pBest,,,,position')") 
end select

!����
write(*,*)  
write(*,"('BEGIN...')")
write(*,"('rang: ')")  
write(*,"(2(g0,2x))")  rang
call particleSwarmOptimization(w, c1, c2, nParticle, rang, nIts, nLineF, size(fre), obsFre, weight, err0, par, par0)
write(*,"('moduli: ', *(g0,', '))") par 
write(*,"('DONE...')")  
write(*,*)  

!�����ϵĲ���
call calFrequence(par, fre)
err = fitnessFun(fre, obsFre, weight, nLineF, .true.)

time1 = rtc()
write(*,"('CPU Time: ', g0, 's')") real(time1-time0,4)
write(30,*)  
write(30,"('CPU Time:, ', g0, 's')") real(time1-time0,4)

close(30)
read(*,*)  
end program ParticleSwarm