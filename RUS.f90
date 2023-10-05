module RUS !Resonant Ultrasound Spectroscopy
implicit none
real d1, d2, d3, rho !��Ʒ�ߴ�(���ڶԳ��ԣ�Ϊʵ��ֵ��һ��), �ܶ�
integer:: shape=1    !��״: 1Բ���壬2��3������
integer:: medimuType !��������
integer:: d          !��϶���ʽ�Ľ���
integer irk(8)       !�洢ÿ���Ӿ��������ֵ����
integer, allocatable:: itab(:), ltab(:), mtab(:), ntab(:) !������
type array
  real,  allocatable:: dat(:,:)
end type
type(array) e0(8)
private
public init, calFrequence
contains 

!****************************************************************************
! �����ڴ�, ����������e
!****************************************************************************
subroutine init(d_in, d1_in, d2_in, d3_in, mass, shape_in, medimuType_in, wsort)
integer, intent(in):: d_in, shape_in, medimuType_in
real,    intent(in):: d1_in, d2_in, d3_in, mass
real, allocatable, intent(out):: wsort(:)
integer i
!��ֵ
d = d_in; d1 = d1_in; d2 = d2_in; d3 = d3_in
rho = mass / volume(d1, d2, d3, shape) !�ܶ�
shape = shape_in; medimuType = medimuType_in
!����������
i = (3*(d+1)*(d+2)*(d+3))/6
allocate(itab(i), ltab(i), mtab(i), ntab(i))
call index_relationship(itab, ltab, mtab, ntab, d, irk)
!���ٿռ�洢Ƶ��
allocate(wsort(sum(irk))); wsort = 0
!����e
do i = 1, 8
  allocate(e0(i)%dat(irk(i),irk(i)))
  call e_fill(e0(i)%dat, itab, ltab, mtab, ntab, d1, d2, d3, shape, i, irk)
end do
end subroutine

!****************************************************************************
! ����Ƶ��
! ����:
! par(:)     ... ��Ʒ����
! ���:
! fre(:)     ... Ƶ�� MHz
!****************************************************************************
pure subroutine calFrequence(par, fre)
use lapack95
implicit none
real, intent(in):: par(:)
real, intent(out):: fre(:)
real, parameter:: PI = acos(-1.0), EPS = 1.0d-5
integer:: info
real c(3,3,3,3)
type(array) e(8), g(8) 
integer i
e = e0; g = e0 !��ʱ�ռ�, �ᱻsygv��Ⱦ
!תΪ4������
c = trans2Cijkl( medimuType, par/rho/100.0 )
do i = 1, 8
  !���� gamma
  call gamma_fill(g(i)%dat, itab, ltab, mtab, ntab, d1, d2, d3, c, shape, i, irk)
  !������ֵ
  call sygv(a=g(i)%dat, b=e(i)%dat, w=fre(sum(irk(:i-1))+1:sum(irk(:i))), itype=1, info=info)
end do
!����
call sort(fre)
!��Ƶ��תΪƵ��
where(fre<EPS) fre = 0.0
fre = sqrt(fre) / (PI*2.0) 
end subroutine

!****************************************************************************
! ���㵯�Բ���4������
! ����:
! mediumType   ... ��������
! c(:)         ... ���Բ���
! isotropic    c11, c44
! cubic        c11, c12, c44
! VTI          c33, c23, c12, c44, c66
! HTI          c11, c33, c12, c44, c66
! tetragonal   c11, c33, c23, c12, c44, c66
! orthorhombic c11, c22, c33, c23, c13, c12, c44, c55, c66
! ����ֵ:
! res(3,3,3,3)
!****************************************************************************
pure function trans2Cijkl( medimuType, c ) result(res)
  implicit none
  integer, intent(in):: medimuType
  real,    intent(in):: c(:)  
  integer,  parameter:: IND(3,3) = reshape([1,6,5,6,2,4,5,4,3],[3,3]) !convert Aijkl into Amn
  real res(3,3,3,3), a(6,6)
  integer i, j, k, L, m, n
  
  res = 0; a = 0
  select case(medimuType)
  case(1) !isotropic
  a(1,1) = c(1);   a(4,4) = c(2)
  a(2,2) = a(1,1); a(3,3) = a(1,1); a(5,5) = a(4,4); a(6,6) = a(4,4)
  a(1,2) = a(1,1) - 2.0*a(4,4)
  a(1,3) = a(1,2); a(2,3) = a(1,2)
  a(2,1) = a(1,2); a(3,1) = a(1,2); a(3,2) = a(1,2) 
  case(2) !cubic
  a(1,1) = c(1);   a(1,2) = c(2);   a(4,4) = c(3)
  a(2,2) = a(1,1); a(3,3) = a(1,1); a(5,5) = a(4,4); a(6,6) = a(4,4)
  a(1,3) = a(1,2); a(2,3) = a(1,2); a(3,1) = a(1,2); a(3,2) = a(1,2)
  a(2,1) = a(1,2)
  case(3) !VTI
  a(3,3) = c(1);   a(2,3) = c(2);   a(1,2) = c(3);   a(4,4) = c(4);   a(6,6) = c(5)
  a(1,1) = 2.0*a(6,6) + a(1,2);     a(2,2) = a(1,1); a(1,3) = a(2,3); a(3,1) = a(2,3)
  a(3,2) = a(2,3); a(2,1) = a(1,2); a(5,5) = a(4,4)
  case(4) !HTI
  a(1,1) = c(1);   a(3,3) = c(2);   a(1,2) = c(3);   a(4,4) = c(4);   a(6,6) = c(5)
  a(2,3) = a(3,3) - 2.0*a(4,4);     a(3,2) = a(2,3); a(1,3) = a(1,2); a(2,1) = a(1,2)
  a(3,1) = a(1,2); a(2,2) = a(3,3); a(5,5) = a(6,6)  
  case(5) !tetragonal
  a(1,1) = c(1);   a(3,3) = c(2);   a(2,3) = c(3);   a(4,4) = c(4);   a(1,2) = c(5);   a(6,6) = c(6)
  a(2,2) = a(1,1); a(1,3) = a(2,3); a(3,1) = a(2,3); a(2,1) = a(1,2)
  a(3,2) = a(2,3); a(5,5) = a(4,4);
  case(6) !orthorhombic
  a(1,1) = c(1);   a(2,2) = c(2);   a(3,3) = c(3);   a(2,3) = c(4);   a(1,3) = c(5)
  a(1,2) = c(6);   a(4,4) = c(7);   a(5,5) = c(8);   a(6,6) = c(9)
  a(3,1) = a(1,3); a(2,1) = a(1,2); a(3,2) = a(2,3)
  end select
  !תΪ4������
  do i = 1, 3
  do j = 1, 3
    m = IND(i,j)
    do k = 1, 3
    do L = 1, 3
      n = IND(k,L)
      res(i,j,k,L) = a(m,n)
    end do; end do
  end do; end do
end function 

!****************************************************************************
! �������
! ����:
! d1, d2, d3     ... ��Ʒ�ߴ��һ��
! l, m, n        ... ����
! shape          ... ��״
! ���:
! res            ... ���ֽ��
!****************************************************************************
elemental function volIntegral(d1, d2, d3, l, m, n, shape) result(res)
implicit none
real,    intent(in):: d1, d2, d3
integer, intent(in):: l, m, n
integer, intent(in):: shape
real, parameter:: PI = acos(-1.0)
real res
res = 0
if(any(mod([l,m,n],2)==1)) return !�溯��, ���Ϊ��
select case(shape)
  case(1) !Բ����
  res = 4.0*PI * d1**(l+1) * d2**(m+1) * d3**(n+1) / (n+1)
  res = res * doubleFact(l-1) * doubleFact(m-1) / doubleFact(l+m+2)
  case(2) !����
  res = 4.0*PI * d1**(l+1) * d2**(m+1) * d3**(n+1) 
  res = res * doubleFact(l-1) * doubleFact(m-1) * doubleFact(n-1) / doubleFact(l+m+n+3)
  case default !������
  res = 8.0 / ((l+1)*(m+1)*(n+1))
  res = res * d1**(l+1) * d2**(m+1) * d3**(n+1)
end select

contains
elemental recursive function doubleFact(n) result(res)
implicit none
integer, intent(in):: n
real res
if(n<2) then
  res = 1
else
  res = n * doubleFact(n-2)
end if
end function
end function

!****************************************************************************
! �������
! ����:
! d1, d2, d3   ... ��Ʒ�ߴ�(һ��)
! shape        ... ��״
! ����ֵ:
! res
!****************************************************************************
elemental function volume(d1, d2, d3, shape) result(res)
implicit none
real,    intent(in):: d1, d2, d3
integer, intent(in):: shape
real, parameter:: PI = acos(-1.0)
real res
res = 0
select case(shape)
  case(1) !Բ����
  res = 2.0*PI * d1 * d2 * d3
  case(2) !����
  res = 4.0/3.0*PI * d1 * d2 * d3
  case default !������
  res = 8.0 * d1 * d2 * d3
end select
end function

!****************************************************************************
! ����������
! ����:
! d            ... ����ʽ����
! ���:
! itab(:), ltab(:), mtab(:), ntab(:)  ... ������
! irk(8)       ... ÿ���Ӿ��������ֵ����
!****************************************************************************
pure subroutine index_relationship(itab, ltab, mtab, ntab, d, irk)
implicit none
integer,  intent(in):: d
integer, intent(out):: itab(:), ltab(:), mtab(:), ntab(:), irk(8)
integer, parameter:: p(3,3,8) = reshape([0,0,0,1,1,0,1,0,1,&
                                         0,0,1,1,1,1,1,0,0,&
                                         0,1,0,1,0,0,1,1,1,&
                                         0,1,1,1,0,1,1,1,0,&
                                         1,0,0,0,1,0,0,0,1,&
                                         1,0,1,0,1,1,0,0,0,&
                                         1,1,0,0,0,0,0,1,1,&
                                         1,1,1,0,0,1,0,1,0 ],[3,3,8])
integer ir, i, k, l, m, n
ir = 0; irk = 0
do k = 1, 8
do i = 1, 3
do l = 0, d
do m = 0, d-l
do n = 0, d-l-m
if( all(mod([l,m,n],2)==p(:,i,k)) ) then
  ir = ir + 1;  irk(k) = irk(k) + 1
  itab(ir) = i; ltab(ir) = l
  mtab(ir) = m; ntab(ir) = n
end if
end do; end do; end do; end do; end do
end subroutine

!****************************************************************************
! ����նȾ��� e
! ����:
! itab(:), ltab(:), mtab(:), ntab(:)  ... ������
! d1, d2, d3   ... ��Ʒ�ߴ�(һ��)
! shape        ... ��״
! irk(8)       ... ÿ���Ӿ��������ֵ����
! k            ... ��k���Ӿ���
! ���:
! e(:,:)       ... �նȾ���
!****************************************************************************
pure subroutine e_fill(e, itab, ltab, mtab, ntab, d1, d2, d3, shape, k, irk)
implicit none
real,   intent(out):: e(:,:)
real,    intent(in):: d1, d2, d3
integer, intent(in):: itab(:), ltab(:), mtab(:), ntab(:), shape, k, irk(8)
integer irs, irf, irv, irh, ir1, ir2, l, m, n
e = 0
irs = sum(irk(:k-1)); irf = irs + irk(k)
irv = 0
do ir1 = irs+1, irf
  irv = irv + 1
  irh = 0
  do ir2 = irs+1, irf
    irh = irh + 1
    if(itab(ir1)==itab(ir2)) then
      l = ltab(ir1) + ltab(ir2)
      m = mtab(ir1) + mtab(ir2)
      n = ntab(ir1) + ntab(ir2)
      e(irh,irv) = volIntegral(d1, d2, d3, l, m, n, shape)
    end if
  end do
end do
end subroutine

!****************************************************************************
! ���� gamma
! ����:
! itab(:), ltab(:), mtab(:), ntab(:)  ... ������
! d1, d2, d3   ... ��Ʒ�ߴ�(һ��)
! c(3,3,3,3)   ... �Ľ׵��Բ���
! shape        ... ��״
! irk(8)       ... ÿ���Ӿ��������ֵ����
! k            ... ��k���Ӿ���
! ���:
! g(:,:)       ... 
!****************************************************************************
pure subroutine gamma_fill(g, itab, ltab, mtab, ntab, d1, d2, d3, c, shape, k, irk)
implicit none
real,   intent(out):: g(:,:)
real,    intent(in):: d1, d2, d3, c(3,3,3,3)
integer, intent(in):: itab(:), ltab(:), mtab(:), ntab(:), shape, k, irk(8)
integer irs, irf, irv, irh, ir1, ir2, j1, j2
integer i1, i2, l1, l2, m1, m2, n1, n2, l, m, n
g = 0
irs = sum(irk(:k-1)); irf = irs + irk(k)
irv = 0
do ir1 = irs+1, irf
  irv = irv + 1
  irh = 0
  do ir2 = irs+1, irf
    irh = irh + 1
    i1=itab(ir1); i2=itab(ir2)
    l1=ltab(ir1); l2=ltab(ir2)
    m1=mtab(ir1); m2=mtab(ir2)
    n1=ntab(ir1); n2=ntab(ir2)
    !L
    if(l1>0) then
      j1 = 1
      if(l2>0) then
        j2 = 1; l = l1+l2-2; m = m1+m2; n = n1+n2
        g(irh,irv) = g(irh,irv) + c(i1,j1,i2,j2)*l1*l2*volIntegral(d1, d2, d3, l, m, n, shape)
      end if
      if(m2>0) then
        j2 = 2; l = l1+l2-1; m = m1+m2-1; n = n1+n2
        g(irh,irv) = g(irh,irv) + c(i1,j1,i2,j2)*l1*m2*volIntegral(d1, d2, d3, l, m, n, shape)
      end if
      if(n2>0) then
        j2 = 3; l = l1+l2-1; m = m1+m2; n = n1+n2-1
        g(irh,irv) = g(irh,irv) + c(i1,j1,i2,j2)*l1*n2*volIntegral(d1, d2, d3, l, m, n, shape)
      end if 
    end if
    !m
    if(m1>0) then
      j1 = 2
      if(l2>0) then
        j2 = 1; l = l1+l2-1; m = m1+m2-1; n = n1+n2
        g(irh,irv) = g(irh,irv) + c(i1,j1,i2,j2)*m1*l2*volIntegral(d1, d2, d3, l, m, n, shape)
      end if
      if(m2>0) then
        j2 = 2; l = l1+l2; m = m1+m2-2; n = n1+n2
        g(irh,irv) = g(irh,irv) + c(i1,j1,i2,j2)*m1*m2*volIntegral(d1, d2, d3, l, m, n, shape)
      end if
      if(n2>0) then
        j2 = 3; l = l1+l2; m = m1+m2-1; n = n1+n2-1
        g(irh,irv) = g(irh,irv) + c(i1,j1,i2,j2)*m1*n2*volIntegral(d1, d2, d3, l, m, n, shape)
      end if 
    end if
    !n
    if(n1>0) then
      j1 = 3
      if(l2>0) then
        j2 = 1; l = l1+l2-1; m = m1+m2; n = n1+n2-1
        g(irh,irv) = g(irh,irv) + c(i1,j1,i2,j2)*n1*l2*volIntegral(d1, d2, d3, l, m, n, shape)
      end if
      if(m2>0) then
        j2 = 2; l = l1+l2; m = m1+m2-1; n = n1+n2-1
        g(irh,irv) = g(irh,irv) + c(i1,j1,i2,j2)*n1*m2*volIntegral(d1, d2, d3, l, m, n, shape)
      end if
      if(n2>0) then
        j2 = 3; l = l1+l2; m = m1+m2; n = n1+n2-2
        g(irh,irv) = g(irh,irv) + c(i1,j1,i2,j2)*n1*n2*volIntegral(d1, d2, d3, l, m, n, shape)
      end if 
    end if  
  end do
end do
end subroutine

!****************************************************************************
! ��������
! �������:
! a(:) 
!****************************************************************************
pure subroutine sort(a)
implicit none
real, intent(inout):: a(:)
real c
integer i, j, k, n
n = size(a)
do i = 1, n-1
  k = i
  do j = i+1, n
    if(a(k)>a(j)) k = j
  end do
  if(k/=i) then
    c = a(i); a(i) = a(k); a(k) = c
  end if
end do
end subroutine
end module

