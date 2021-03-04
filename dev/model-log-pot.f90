!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> University College London
!>
!>@brief
!> Store the equations of motion to be solved by our nonlinear boundary value solver
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Preprocessors for inlining.
#define POTENTIAL(f) ( 0.5_dl*f**2*(1._dl-log(f**2 + eps**2)) )
#define VPRIME(f) ( -f*(log(f**2+eps**2) - eps**2/(f**2+eps**2)) )
#define VDPRIME(f) ( -log(f**2+eps**2) + (eps**4-3._dl*eps**2*f**2-2._dl*f**4)/(eps**2+f**2)**2 )


module Model
  use constants
  use Cheby
  implicit none
  private :: eps, i_
  
  real(dl) :: eps
  real(dl), parameter :: eps_min = 0._dl, eps_max = exp(0.5_dl)
  integer :: i_
  integer, parameter :: nEps = 200
  real(dl), parameter :: log_e_min = -7._dl
  real(dl), parameter :: scanVals(nEps) = 10.**(/ ( log_e_min + (i_-1)*(log10(eps_max)-log_e_min)/dble(100) , i_=1,nEps ) /)  ! Scan parameters
  
contains
  
  subroutine set_model_params(params,dim)
    real(dl), dimension(1), intent(in) :: params
    integer, intent(in) :: dim
    eps = params(1)
  end subroutine set_model_params
  
  subroutine grid_params_(w,len,r0,w0)
    real(dl), intent(out) :: w, len
    real(dl), intent(in) :: r0, w0

    len = 2._dl; w = 1._dl
  end subroutine grid_params_

  subroutine bubble_parameters_nd_(delta,dim,r0,meff)
    real(dl), intent(in) :: delta, dim
    real(dl), intent(out) :: r0, meff

    meff = 1._dl
    r0 = 1._dl/meff
  end subroutine bubble_parameters_nd_

  subroutine get_minima(phif,phit)
    real(dl), intent(out) :: phif, phit
    phif = 0._dl; phit = 10._dl
  end subroutine get_minima

  elemental function potential(phi)
    real(dl) :: potential
    real(dl), intent(in) :: phi
    potential = POTENTIAL(phi)
  end function potential

  elemental function vprime(phi)
    real(dl) :: vprime
    real(dl), intent(in) :: phi
    vprime =  VPRIME(phi)
  end function vprime

  elemental function vdprime(phi)
    real(dl) :: vdprime
    real(dl), intent(in) :: phi
    vdprime =  VDPRIME(phi)
  end function vdprime

  !>@brief
  !> The unperturbed part of the potential used in computing the thin-wall domain wall
  elemental function potential_tw(phi)
    real(dl) :: potential_tw
    real(dl), intent(in) :: phi
    potential_tw = 0._dl
  end function potential_tw

  !>@brief
  !> A subroutine to set general Robin boundary conditions on our fields
  !>  \f[
  !>    \alpha_L f(x_L) + \beta_L f'(x_L) = c_L
  !>  \f]
  !>  \f[
  !>    \alpha_R f(x_R) + \beta_R f'(x_R) = c_R
  !>  \f]
  subroutine boundaries(L,S,c,bc)
    real(dl), intent(inout) :: L(1:,1:), S(1:)
    real(dl), dimension(1:3,1:2), intent(in) :: c
    logical, dimension(1:2), intent(in) :: bc
    integer :: sz
    sz = size(S)

    if (bc(1)) then
       L(:,:) = c(1,1) + c(2,1)
       S(1) = c(3,1)
    endif
    if (bc(2)) then
       L(:,:) = c(1,2) + c(2,2)
       S(sz) = c(3,2)
    endif
  end subroutine boundaries

end module Model
