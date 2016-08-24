module Pseudospec
  use constants
  implicit none
  
  real(dl), dimension(XRANGE) :: xgrid, wgrid
  real(dl), dimension(XRANGE,BRANGE) :: b, bp, bpp, mmt
  real(dl), dimension(XRANGE) :: freal, del
  real(dl), dimension(BRANGE) :: fspec

contains

!!!!!!!!!!!!!!!!!!!!!!!!
! Chebyshev Polynomials
!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gauss-Chebychev "zeros" grid on interval [-1:1]
  subroutine gauss_cheby(xvals, wvals, order)
    real(dl), intent(out) :: xvals(0:order), wvals(0:order)
    integer, intent(in) :: order

    integer :: i
    real(dl) :: dkcol

    dkcol = twopi / 4. /dble(order+1)
    do i=0,order
       xvals(i) = -cos((2*i+1)*dkcol)  ! negative sign just puts the points in increasing order from [-1:1]
    enddo
    wvals = 2. / dble(order+1)
  end subroutine gauss_cheby

  !> Gauss-Lobatto endpoints and extrema grid for Chebychev polynomials
  subroutine gauss_lobatto(xvals, wvals, order)
    real(dl), intent(out) :: xvals(XRANGE), wvals(XRANGE)
    integer, intent(in) :: order

    integer :: i
    real(dl) :: dkcol
    
    dkcol = twopi / 2. / dble(order)
    do i=0,order
       xvals(i) = -cos(i*dkcol)
    enddo
    wvals = 2. / dble(order)
    wvals(0) = 1. / dble(order); wvals(order) = 1. / dble(order)
  end subroutine gauss_lobatto

  !> Compute Chebychev polynomials up to order nmax at position x.
  !> Also computes first and second derivatives.
  subroutine chebychev(nmax,x,T,Tp,Tpp)
    integer, intent(in) :: nmax
    real(dl), intent(in) :: x
    real(dl), dimension(0:nmax), intent(out) :: T, Tp, Tpp

    real(dl) :: xt, cn, sn
    real(dl) :: pt, ptt
    integer :: i

    T(0) = 1.
    Tp(0) = 0.
    Tpp(0) = 0.
    T(1) = x
    Tp(1) = 1.
    Tpp(1) = 0.
    do i=2,nmax
       T(i) = 2.*x*T(i-1) - T(i-2)
       Tp(i) = 2.*x*Tp(i-1) - Tp(i-2) + 2.*T(i-1)
       Tpp(i) = 2.*x*Tpp(i-1) - Tpp(i-2) + 4.*Tp(i-1)
    enddo
  end subroutine chebychev


!!!!!!!!!!!!!!!!!!!!!!!
! Legendre Polynomials
!!!!!!!!!!!!!!!!!!!!!!!

  !> Gaussian quadrature for Legendre polynomials
  !> @todo Write this
  subroutine gauss_lobatto_legendre(xvals, order)
    real, intent(out) :: xvals(0:order)
    integer, intent(in) :: order

    xvals = 0.
  end subroutine gauss_lobatto_legendre

  !> Gauss-Lobatto quadrature for Legendre polynomials
  !> @todo Write this
  subroutine gauss_legendre(xvals, order)
    real, intent(out) :: xvals(0:order)
    integer, intent(in) :: order

    xvals = 0.
  end subroutine gauss_legendre

  !> Compute Legendre polynomials up to order lmax at position x,
  subroutine legendre(P, lmax, x)
    integer :: l, lmax
    real :: P(0:lmax), x

    P(0) = 1.
    P(1) = x
    do l=2,lmax
       P(l) = ((2*l-1)*x*P(l-1) - (l-1)*P(l-2))/l
    enddo
  end subroutine legendre


end module Pseudospec
