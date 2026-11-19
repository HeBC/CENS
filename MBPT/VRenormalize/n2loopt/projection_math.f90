! Angular projection only: no radial measure, units, or (2*pi)**(-3).
! Incoming momentum is along +z; outgoing momentum lies in the xz plane.
! Spin basis: 1=(S=0,M=0), 2=(1,-1), 3=(1,0), 4=(1,+1).
! Condon-Shortley spherical harmonics and Clebsch-Gordan convention.
! Raw momentum angular states: the adapter, not this module, must apply
! i**(Lin-Lout) for the real radial-wave convention of NuHamil and the
! source's named partial-wave contacts. For equal parities this factor is
! (-1)**((Lout-Lin)/2); it changes only coupled off-diagonal signs here.
module source2_projection_math
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: gauss_legendre, projection_weights, source2_cg
  real(real64), parameter :: pi = acos(-1.0_real64)
contains

  pure real(real64) function fact(n) result(value)
    integer, intent(in) :: n
    integer :: k
    value = 1.0_real64
    do k=2,n
      value = value*real(k,real64)
    end do
  end function fact

  pure real(real64) function source2_cg(j1,m1,j2,m2,j,m) result(value)
    integer, intent(in) :: j1,m1,j2,m2,j,m
    integer :: k,low,high
    real(real64) :: prefactor,series
    value = 0.0_real64
    if (m1+m2 /= m .or. abs(m1)>j1 .or. abs(m2)>j2 .or. abs(m)>j) return
    if (j<abs(j1-j2) .or. j>j1+j2 .or. min(j1,j2,j)<0) return
    prefactor = sqrt(real(2*j+1,real64)*fact(j+j1-j2)*fact(j-j1+j2) &
                     *fact(j1+j2-j)/fact(j1+j2+j+1))
    prefactor = prefactor*sqrt(fact(j+m)*fact(j-m)*fact(j1-m1)*fact(j1+m1) &
                               *fact(j2-m2)*fact(j2+m2))
    low = max(0,j2-j-m1,j1+m2-j)
    high = min(j1+j2-j,j1-m1,j2+m2)
    series = 0.0_real64
    do k=low,high
      series = series + (-1.0_real64)**k &
           /(fact(k)*fact(j1+j2-j-k)*fact(j1-m1-k)*fact(j2+m2-k) &
             *fact(j-j2+m1+k)*fact(j-j1-m2+k))
    end do
    value = prefactor*series
  end function source2_cg

  pure real(real64) function y_phi_zero(l,m,x) result(value)
    integer, intent(in) :: l,m
    real(real64), intent(in) :: x
    integer :: a,i,ll
    real(real64) :: pmm,pmmp1,pll,p_previous,p_current
    value = 0.0_real64
    a = abs(m)
    if (a>l) return
    pmm = 1.0_real64
    do i=1,a
      pmm = -pmm*real(2*i-1,real64)*sqrt(max(0.0_real64,1.0_real64-x*x))
    end do
    if (l==a) then
      pll = pmm
    else
      pmmp1 = x*real(2*a+1,real64)*pmm
      pll = pmmp1
      p_previous = pmm
      p_current = pmmp1
      do ll=a+2,l
        pll = (real(2*ll-1,real64)*x*p_current &
               -real(ll+a-1,real64)*p_previous)/real(ll-a,real64)
        p_previous = p_current
        p_current = pll
      end do
    end if
    value = sqrt(real(2*l+1,real64)/(4.0_real64*pi)*fact(l-a)/fact(l+a))*pll
    if (m<0) value = value*(-1.0_real64)**a
  end function y_phi_zero

  subroutine gauss_legendre(n,x,w)
    integer, intent(in) :: n
    real(real64), intent(out) :: x(n),w(n)
    integer :: i,j,iteration
    real(real64) :: z,z_previous,p1,p2,p3,derivative
    if (n<1) error stop 'gauss_legendre: n must be positive'
    do i=1,(n+1)/2
      z = cos(pi*(real(i,real64)-0.25_real64)/(real(n,real64)+0.5_real64))
      do iteration=1,100
        p1 = 1.0_real64
        p2 = 0.0_real64
        do j=1,n
          p3 = p2
          p2 = p1
          p1 = (real(2*j-1,real64)*z*p2-real(j-1,real64)*p3)/real(j,real64)
        end do
        derivative = real(n,real64)*(z*p1-p2)/(z*z-1.0_real64)
        z_previous = z
        z = z_previous-p1/derivative
        if (abs(z-z_previous)<4.0_real64*epsilon(z)) exit
      end do
      if (iteration>100) error stop 'gauss_legendre: root did not converge'
      ! Reevaluate derivative at the converged root, including the last step.
      p1 = 1.0_real64
      p2 = 0.0_real64
      do j=1,n
        p3 = p2
        p2 = p1
        p1 = (real(2*j-1,real64)*z*p2-real(j-1,real64)*p3)/real(j,real64)
      end do
      derivative = real(n,real64)*(z*p1-p2)/(z*z-1.0_real64)
      x(i) = -z
      x(n+1-i) = z
      w(i) = 2.0_real64/((1.0_real64-z*z)*derivative*derivative)
      w(n+1-i) = w(i)
    end do
  end subroutine gauss_legendre

  subroutine projection_weights(j,s,lo,li,n,x,w,weights)
    integer, intent(in) :: j,s,lo,li,n
    real(real64), intent(in) :: x(n),w(n)
    real(real64), intent(out) :: weights(n,4,4)
    integer :: mi,mo,ml,ii,oi,k
    real(real64) :: prefactor,incoming,outgoing
    if (s<0 .or. s>1) error stop 'projection_weights: spin must be 0 or 1'
    if (min(lo,li)<0 .or. j<abs(lo-s) .or. j>lo+s &
        .or. j<abs(li-s) .or. j>li+s) error stop 'projection_weights: invalid LSJ triangle'
    weights = 0.0_real64
    prefactor = 8.0_real64*pi*pi/real(2*j+1,real64) &
                *sqrt(real(2*li+1,real64)/(4.0_real64*pi))
    do mi=-s,s
      incoming = source2_cg(li,0,s,mi,j,mi)
      if (incoming==0.0_real64) cycle
      ii = 1
      if (s==1) ii=mi+3
      do mo=-s,s
        ml = mi-mo
        outgoing = source2_cg(lo,ml,s,mo,j,mi)
        oi = 1
        if (s==1) oi=mo+3
        do k=1,n
          weights(k,oi,ii) = prefactor*incoming*outgoing*y_phi_zero(lo,ml,x(k))*w(k)
        end do
      end do
    end do
  end subroutine projection_weights
end module source2_projection_math
