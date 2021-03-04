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
#define POTENTIAL(f) ( cos(f) + del*sin(f)**2 + 1._dl )
#define VPRIME(f) ( -sin(f) + del*sin(2._dl*f) )
#define VDPRIME(f) ( -cos(f) + 2._dl*del*cos(2._dl*f) )

module Model
  use constants
  use Cheby
  implicit none
  private :: del
  
  real(dl) :: del

contains

  subroutine set_model_params(params,dim)
    real(dl), dimension(1), intent(in) :: params
    integer, intent(in) :: dim
    del = params(1)
  end subroutine set_model_params
  
  subroutine get_minima(phif,phit)
    real(dl), intent(out) :: phif, phit
    phif = 0._dl; phit = pi
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
  !> Given specified radius and width of a bubble profile, adjust grid mapping parameters.
  !>
  !> The relationship between the radius and mapping length are fixed by choice of polynomials
  !> Should probably be moved into the chebyshev class
  subroutine grid_params_(w,len,r0,w0)
    real(dl), intent(out) :: w, len
    real(dl), intent(in) :: r0, w0
    real(dl), parameter :: wscl = 6.3_dl   ! TWEAK THIS!!!!!
    
    len = r0*3._dl**0.5
    w = wscl * w0 / r0
    if (w0 > r0) then
       len = 3._dl**0.5*w0
       w = 1._dl
    endif
  end subroutine grid_params_
  
  ! Change delta to parameters for the model
  subroutine bubble_parameters_nd_(delta,dim,r0,meff)
    real(dl), intent(in) :: delta, dim
    real(dl), intent(out) :: r0, meff

    meff = sqrt(2._dl*delta)*(1._dl-0.25_dl/delta**2)**0.5
    r0 = dim * sqrt(2._dl*delta)
  end subroutine bubble_parameters_nd_

  !>@brief
  !> The unperturbed part of the potential used in computing the thin-wall domain wall
  elemental function potential_tw(phi)
    real(dl) :: potential_tw
    real(dl), intent(in) :: phi
    potential_tw = del*sin(phi)**2
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
