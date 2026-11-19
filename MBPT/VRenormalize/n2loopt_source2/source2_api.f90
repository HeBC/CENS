! Continuous-momentum interface to the corrected thesis implementation.
module source2_api
  use iso_c_binding
  use chiral_constants_nnlo_opt
  use chiral_potentials_nnlo_opt
  implicit none
  private
  public :: source2_initialize, source2_kernel, source2_raw
  integer :: initialized_variant=-1
contains
  subroutine source2_initialize(variant) bind(C,name='source2_initialize')
    integer(c_int), value :: variant
    real(8) :: energy, common_cs, common_ct
    if(variant/=0.and.variant/=1) error stop 'variant must be 0 or 1'
    if(initialized_variant==variant) return
    CE=0.d0; CD=0.d0
    call init_chp_constants_nnlo_opt(1.d0,1,energy)
    ! Remove finite-box scales: input vectors now MeV, result MeV**(-2).
    mom_scale=1.d0/hbarc
    iescale=1.d0
    if(variant==1) then
      mnuc(:)=0.5d0*(proton_mass+neutron_mass)
      mnuc2=mnuc**2
      mpi(:)=(pi_m_mass+pi_0_mass+pi_p_mass)/3.d0
      mpi2=mpi**2; mpi4=mpi**4; mpi5=mpi**5
      twompi=2.d0*mpi; fourmpi2=4.d0*mpi2
      loop_s=chp_NLO_sfr_two_pion_exchange_loop_s(2)
      ! Average charge-dependent singlet contacts, preserving triplet contact.
      common_cs=sum(CS)/3.d0; common_ct=sum(CT)/3.d0
      CS=common_cs; CT=common_ct
    endif
    initialized_variant=variant
  end subroutine

  subroutine source2_raw(pout,pin,m1,m2,m3,m4,t1,t2,t3,t4,value) bind(C,name='source2_raw')
    real(c_double), intent(in) :: pout(3),pin(3)
    integer(c_int), value :: m1,m2,m3,m4,t1,t2,t3,t4
    complex(c_double_complex), intent(out) :: value
    type(t_vector) :: a,b,c,d
    real(8) :: saved_average_mass
    if(initialized_variant<0) call source2_initialize(0_c_int)
    a%k=pout; b%k=-pout; c%k=pin; d%k=-pin
    a%s=m1; b%s=m2; c%s=m3; d%s=m4
    a%t=t1; b%t=t2; c%t=t3; d%t=t4
    ! The upstream two-pion 1/m functions read mnuc(0). For the physical
    ! variant use the pp, pn, or nn nucleon mass during this evaluation.
    ! Preserve the common mass for the symmetric variant and restore state
    ! so later pn calls cannot inherit a preceding pp/nn mass.
    saved_average_mass=mnuc(0)
    if(initialized_variant==0) mnuc(0)=mnuc((t1+t2)/2)
    value=chiral_pot_nnlo_opt_np(a,b,c,d)
    mnuc(0)=saved_average_mass
  end subroutine

  subroutine source2_kernel(pout,pin,n,x,kernel) bind(C,name='source2_kernel')
    real(c_double), value :: pout,pin
    integer(c_int), value :: n
    real(c_double), intent(in) :: x(n)
    ! basis: singlet, triplet M=-1,0,+1; charge slots pp,pn T1,pn T0,nn.
    real(c_double), intent(out) :: kernel(n,4,4,4)
    real(8) :: u(4,4),a(3),b(3),mat(4,4,4),sq
    complex(8) :: direct,exchange
    integer :: ix,i,k,z,ma(4),mb(4)
    if(initialized_variant<0) call source2_initialize(0_c_int)
    sq=1.d0/sqrt(2.d0)
    ! Product basis ++,+-,-+,--.
    ma=[1,1,-1,-1]; mb=[1,-1,1,-1]
    u=0.d0
    u(2,1)=sq;u(3,1)=-sq
    u(4,2)=1.d0;u(2,3)=sq;u(3,3)=sq;u(1,4)=1.d0
    b=[0.d0,0.d0,pin]
    do ix=1,n
      a=[pout*sqrt(max(0.d0,1.d0-x(ix)**2)),0.d0,pout*x(ix)]
      do k=1,4
        do i=1,4
          call source2_raw(a,b,ma(i),mb(i),ma(k),mb(k),-1,-1,-1,-1,direct)
          mat(i,k,1)=real(direct)
          call source2_raw(a,b,ma(i),mb(i),ma(k),mb(k),1,1,1,1,direct)
          mat(i,k,4)=real(direct)
          call source2_raw(a,b,ma(i),mb(i),ma(k),mb(k),1,-1,1,-1,direct)
          call source2_raw(a,b,ma(i),mb(i),ma(k),mb(k),1,-1,-1,1,exchange)
          mat(i,k,2)=real(direct+exchange)
          mat(i,k,3)=real(direct-exchange)
        enddo
      enddo
      do z=1,4
        kernel(ix,:,:,z)=matmul(transpose(u),matmul(mat(:,:,z),u))
      enddo
    enddo
  end subroutine
end module
