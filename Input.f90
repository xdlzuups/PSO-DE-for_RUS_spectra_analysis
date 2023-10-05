module Input
contains
subroutine readFile(d, d1, d2, d3, mass, convergencePar, medimuType, nIts, nLineF, par, fre, weight)
implicit none
integer nPar                !��������
integer d                   !����ʽ����
integer nLineC              !�����Ƶ����
real    mass                !����
real    convergencePar      !��������
integer nIts                !��������
integer nLineF              !������ϵ�Ƶ����
real    poison              !���Ʋ��ɱ�
real, allocatable:: par(:)  !���Ʋ���
real    d1, d2, d3          !��Ʒ�ߴ�
real, allocatable:: fre(:), weight(:)  !��������Ƶ�ʼ�Ȩ��
integer medimuType !��������
real fre1, fre2, w
integer i 

!��ȡ�����ļ�
open(10,file='cylin.dat')
read(10,*)
read(10,*) nPar, d, nLineC, mass, convergencePar, nIts, nLineF, poison
allocate(par(nPar))
read(10,*) par(1:nPar) !���Ʋ���
read(10,*) d1, d2, d3  !��Ʒ�ߴ�
!��ȡ
do 
  read(10,*,iostat=i) fre1, fre2, w
  if(i/=0) exit
  fre = [fre,fre1]
  weight = [weight,w]
end do
close(10)

nLineF = min(nLineF,count(weight>1.0e-5))

!ȷ����������
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

!���ڶԳ��ԣ���Ʒ�ߴ��������
d1 = 0.5 * d1; d2 = 0.5 * d2; d3 = 0.5 * d3
end subroutine

!��ȡ������Ϣ
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