module source2_partialwave_api
  use iso_c_binding
  use source2_api
  use source2_projection_math
  implicit none
  private
  public :: source2_partialwaves
  integer, parameter :: nq=128,maxj=20
  real(8), save :: x(nq),qw(nq),weights(nq,4,4,6,0:maxj)
  logical, save :: ready=.false.
contains
  subroutine source2_partialwaves(pout,pin,j,inn,variant,v) bind(C,name='source2_partialwaves')
    real(c_double), value :: pout,pin
    integer(c_int), value :: j,inn,variant
    real(c_double), intent(out) :: v(6)
    real(8) :: kernel(nq,4,4,4),pi
    integer :: jj,c,s,lo,li,t,z
    if(j<0.or.j>maxj) error stop 'source2: J must lie between 0 and 20'
    if(inn<1.or.inn>3) error stop 'source2: inn must be 1,2,3'
    if(pout<0.or.pin<0) error stop 'source2: momenta must be nonnegative'
    pi=acos(-1.d0)
    if(.not.ready) then
      call gauss_legendre(nq,x,qw)
      weights=0.d0
      do jj=0,maxj
        do c=1,6
          call quantum_numbers(jj,c,s,lo,li)
          if(lo<0.or.li<0.or.(jj==0.and.c==2)) cycle
          call projection_weights(jj,s,lo,li,nq,x,qw,weights(:,:,:,c,jj))
          ! Radial-wave phase i**(Lin-Lout), then continuum normalization.
          weights(:,:,:,c,jj)=weights(:,:,:,c,jj)*(-1.d0)**((lo-li)/2)/(2.d0*pi)**3
        enddo
      enddo
      ready=.true.
    endif
    call source2_initialize(variant)
    call source2_kernel(pout,pin,nq,x,kernel)
    v=0.d0
    do c=1,6
      call quantum_numbers(j,c,s,lo,li)
      if(lo<0.or.li<0.or.(j==0.and.c==2)) cycle
      t=mod(lo+s+1,2)
      if(inn/=2.and.t/=1) cycle
      select case(inn)
      case(1);z=1
      case(2);z=3-t
      case(3);z=4
      end select
      v(c)=sum(kernel(:,:,:,z)*weights(:,:,:,c,j))
    enddo
  end subroutine

  subroutine quantum_numbers(j,c,s,lo,li)
    integer,intent(in) :: j,c
    integer,intent(out) :: s,lo,li
    s=1
    select case(c)
    case(1);s=0;lo=j;li=j
    case(2);lo=j;li=j
    case(3);lo=j+1;li=j+1
    case(4);lo=j-1;li=j-1
    case(5);lo=j+1;li=j-1
    case(6);lo=j-1;li=j+1
    end select
  end subroutine
end module

! Legacy CENS calling interface. Coulomb is absent.
subroutine n2loopt_cens(variant)
  use source2_partialwave_api
  implicit none
  integer,intent(in) :: variant
  integer :: j,inn
  real(8) :: v(6),xmev,ymev
  logical :: heform,sing,trip,coup,endep
  character(4) :: label
  common /cpot/ v,xmev,ymev
  common /cstate/ j,heform,sing,trip,coup,endep,label
  common /cnn/ inn
  if(heform) error stop 'source2 CENS adapter supports LSJ, not helicity'
  call source2_partialwaves(xmev,ymev,j,inn,variant,v)
  if(.not.sing) v(1)=0.d0
  if(.not.trip) v(2)=0.d0
  if(.not.coup) v(3:6)=0.d0
  endep=.false.
  label='N2op'
end subroutine
