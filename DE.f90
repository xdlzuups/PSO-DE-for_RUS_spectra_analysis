module Differential_Evolution
real :: F = 0.8, CR = 0.2
contains
subroutine DE(F, CR ,np ,nit, fit, bound, error, c, c0)
implicit none
real, intent(in):: F      !��������
real, intent(in):: CR     !��������
integer, intent(in):: np  !��Ⱥ��
integer, intent(in):: nit !��������
real, intent(in):: bound(:,:) !�߽�
real, intent(out):: c(:)   !����
real, intent(in):: error  !���
real, intent(in), optional::c0(:) !��ʼֵ
interface
real function fit(c)
real, intent(in):: c(:)
end function
end interface

integer nd  !��������
real,allocatable:: Xi(: ,:) !����
real,allocatable:: Vi(: ,:) !�м���
real,allocatable:: Ui(: ,:) !�Ӵ�
real,allocatable:: err0(: ) !���
logical, save:: isInit = .true.
real,allocatable:: Rand(:)
real rand1
integer r1 ,r2 ,r3, it, i, j, k

if(isInit) then
  isInit = .false.
  call RANDOM_SEED()
end if
open(10 ,file = 'result.csv') !���
open(11 ,file = 'err.csv')

nd = size(c)
allocate(Xi(nd ,np))
allocate(Vi(nd ,np))
allocate(Ui(nd ,np))
Xi = 0.0
Vi = 0.0
Ui = 0.0
allocate(err0(np))
allocate(Rand(nd))

!��ʼ��
do i = 1 ,np
  call RANDOM_NUMBER(rand)
  if(i==1.and.present(c0)) then
    Xi(:,1) = c0
  else
    Xi(: ,i) = bound(1 ,:) + rand(:) * (bound(2 ,:) - bound(1 ,:))
  end if
  err0(i) = Fit(Xi(:,i))
enddo

do it = 1, nit
  !�������
  do i = 1 ,np
    call MY_INT_RAND(i ,np ,r1 ,r2 ,r3)
    Vi(:,i) = Xi(:,r1) + F * (Xi(:,r2) - Xi(:,r3))
  enddo

  do i = 1 ,np
    !�������
    call RANDOM_NUMBER(rand1)
    k = int(rand1*nd-0.001) + 1
    Ui(k,i) = Vi(k,i)
    do j = 1, nd
      if (j /= k)then
        call RANDOM_NUMBER(rand1)
        if (rand1 <= CR)then
          Ui(j ,i) = Vi(j ,i)
        else
          Ui(j ,i) = Xi(j ,i)
        endif
      endif
      !�߽紦��
      if (Ui(j ,i) > bound(2 ,j) .or. Ui(j ,i) < bound(1 ,j))then
        call RANDOM_NUMBER(rand1)
        Ui(j ,i) = bound(1 ,j) + rand1 * (bound(2 ,j) - bound(1 ,j))
      endif
    enddo
    !�������
    rand1 = Fit(Ui(: ,i))
    if(rand1 < err0(i))then
      Xi(: ,i) = Ui(: ,i)
      err0(i) = rand1
    endif
  end do

  !��С������
  i = minloc(err0,1)
  c = Xi(:,i)
  write(10,"(*(g0,:,','))") it, i, err0(i)
  write(11,"(*(g0,:,','))") it, err0(i)
  write(*,"(*(g0,:,'  '))") it, err0(i), Xi(:,i)
  do j = 1, np
    write(10,"(*(g0,:,','))") j, err0(j), Xi(:,j)
  end do
  if(err0(i) < error) exit
end do
!���ŵ�10����
write(10,*)  
do i = 1, min(10,np)
  j = minloc(err0,1)
  write(10,"(*(g0,:,','))") j, err0(j), Xi(:,j)
  err0(j) = huge(1.0)
end do
close(10)
close(11)
end subroutine

subroutine MY_INT_RAND(i ,np ,r1 ,r2 ,r3)
implicit none
integer, intent(in):: i ,np
integer, intent(out):: r1 ,r2 ,r3
real rand1(3)
integer i1, i2, i3

call RANDOM_NUMBER(rand1)
r1 = int(rand1(1)*(np-2)+0.0001) + 1 ! 1:np-1
if(r1>=i) r1 = r1 + 1
if(r1>i) then
  i1 = i
  i2 = r1
else
  i1 = r1
  i2 = i
end if

r2 = int(rand1(2)*(np-3)+0.0001) + 1! 1:np-2
if(r2>=i2) then
  r2 = r2 + 2
else if(r2>=i1) then
  r2 = r2 + 1
end if
if(r2>i2) then
  i3 = r2
else if(r2>i1) then
  i3 = i2
  i2 = r2
else
  i3 = i2
  i2 = i1
  i1 = r2
end if

r3 = int(rand1(3)*(np-4)+0.0001) + 1! 1:np-3
if(r3>=i3) then
  r3 = r3 + 3
else if(r3>=i2) then
  r3 = r3 + 2
else if(r3>=i1) then
  r3 = r3 + 1
end if
end subroutine

subroutine readDEpar(nParticle, F, CR, rang)
implicit none
integer nParticle
real F, CR, rang(:,:)
open(10,file='DEpar.txt',status='old')
read(10,*); read(10,*) nParticle; write(*,"('nParticle: ', g0)")  nParticle
read(10,*); read(10,*) F, CR; write(*,"('F, CR: ', *(g0,2x))")  F, CR
read(10,*); read(10,*) rang; 
close(10)
end subroutine
end module
