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

! Preprocessors for inlining.  Could also move this to another file
#define POTENTIAL(f) ( 2._dl*beta*f**2*( (2.*beta+1._dl-d_space)/(2._dl*beta+1._dl)*abs(f)**(1._dl/beta) - abs(f)**(2._dl/beta) ) )
#define VPRIME(f) ( 2._dl*f*( (2._dl*beta+1._dl-d_space)*abs(f)**(1._dl/beta) - 2._dl*(beta+1._dl)*abs(f)**(2._dl/beta) ) )
#define VDPRIME(f) ( 2._dl*(beta+1._dl)/(beta)*( (2._dl*beta+1._dl-d_space)*abs(f)**(1._dl/beta) - 2._dl*(beta+2._dl)*abs(f)**(2._dl/beta) ) )

module Model
  use constants
  use Cheby
  implicit none
  private ::  ndim, beta, d_space
  
  real(dl) :: beta, d_space
  integer :: ndim
  logical :: init = .false.

contains

  subroutine set_model_params(params,dim)
    real(dl), dimension(:), intent(in) :: params
    integer, intent(in) :: dim
    beta = params(1); d_space=dble(dim)
    ndim = dim
    init = .true.
  end subroutine set_model_params

  subroutine get_minima(phif,phit)
    real(dl), intent(out) :: phif, phit
    phif = 0._dl; phit = 1._dl
  end subroutine get_minima

  real(dl) elemental function potential(phi)
    real(dl), intent(in) :: phi
    potential = POTENTIAL(phi)
  end function potential

  real(dl) elemental function vprime(phi)
    real(dl), intent(in) :: phi
    vprime =  VPRIME(phi)
  end function vprime

  real(dl) elemental function vdprime(phi)
    real(dl), intent(in) :: phi
    vdprime =  VDPRIME(phi)
  end function vdprime
  
  !>@brief
  !> Given specified radius and width of a bubble profile, adjust grid mapping parameters.
  !>
  !> The relationship between the radius and mapping length are fixed by choice of polynomials
  !> Should probably be moved into the chebyshev class
  subroutine grid_params_(w,len,r0,w0)
    real(dl), intent(out) :: w, len
    real(dl), intent(in) :: r0, w0
    
    len = 1._dl; w = 2._dl
  end subroutine grid_params_

  ! These need to be adjusted for every model.  Might be worth moving it
  ! Change delta to parameters for the model
  subroutine bubble_parameters_nd_(delta,dim,r0,meff)
    real(dl), intent(in) :: delta, dim
    real(dl), intent(out) :: r0, meff

    meff = 1._dl 
    r0 = delta  ! this is beta
  end subroutine bubble_parameters_nd_

  !>@brief
  !> The unperturbed part of the potential used in computing the thin-wall domain wall
  elemental function potential_tw(phi)
    real(dl) :: potential_tw
    real(dl), intent(in) :: phi
    potential_tw = 0._dl  !del*sin(phi)**2  This is clearly wrong, fix it
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
