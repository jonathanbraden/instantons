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
  private ::  beta, d_space, i_
  
  real(dl) :: beta, d_space
  integer, parameter :: nPar = 1  ! number of parameters
  logical :: init = .false.

  integer :: i_
  integer, parameter :: nBeta = 100
  real(dl), parameter :: log_b_min=-5., log_b_max=3.
  real(dl), dimension(1:nBeta), parameter :: scanVals = 1._dl + 10.**(/ (log_b_min+(i_-1)*(log_b_max-log_b_min)/dble(nBeta),i_=1,nBeta) /)

  integer, parameter :: p_guess = 5
contains
  
  subroutine set_model_params(params,dim)
    real(dl), intent(in) :: params(1:nPar), dim
    beta = params(1); d_space=dble(dim)
    init = .true.
  end subroutine set_model_params

  function get_model_params() result(params)
    real(dl), dimension(1:nPar) :: params
    params(1) = beta
  end function get_model_params
  
  subroutine get_minima(phif,phit)
    real(dl), intent(out) :: phif, phit
    phif = 0._dl; phit = 1._dl
  end subroutine get_minima

  ! Adjust this if I change my potential definitions
  real(dl) function phi_max()
    phi_max = ( (2._dl*beta+1._dl-d_space)/(2._dl*(1._dl+beta)) )**beta
  end function phi_max
  
  ! Write this
  real(dl) function phi_out()
    phi_out = 0._dl
  end function phi_out
  
  real(dl) function meff_max() result(meff)
    meff = (2._dl*beta+1._dl-d_space)/sqrt(beta*(1._dl+beta)) 
  end function meff_max
  
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

  subroutine get_grid_params(params,dim,len,scl,pow)
    real(dl), intent(in) :: params(1:nPar), dim
    real(dl), intent(out) :: len, scl
    real(dl), intent(out), optional :: pow

    real(dl) :: beta_loc, meff
    beta_loc = params(1)
    meff = meff_max()

    print*,"meff is ",meff
    len = 2. / meff
    scl = 1._dl
    if (beta_loc < 2.) len = 1._dl
    if (present(pow)) then
       pow = 1._dl
       if (beta_loc < 2.) pow = 2./beta_loc
    endif
  end subroutine get_grid_params

  subroutine get_guess_params(params,dim,params_ic,prof_type)
    real(dl), intent(in) :: params(1:nPar), dim
    real(dl), dimension(1:4), intent(out) :: params_ic
    integer, intent(out), optional :: prof_type

    params_ic(1) = params(1)
    params_ic(2) = 1._dl
    params_ic(3) = 0._dl
    params_ic(4) = 1._dl
    if (present(prof_type)) prof_type = 5
  end subroutine get_guess_params

  logical function use_previous(params,dim) result(prev)
    real(dl), intent(in) :: params(1:nPar), dim
    prev = .false.
  end function use_previous
  
  !>@brief
  !> Given specified radius and width of a bubble profile, adjust grid mapping parameters.
  !>
  !> The relationship between the radius and mapping length are fixed by choice of polynomials
  !> Should probably be moved into the chebyshev class
  subroutine grid_params_(w,len,r0,w0)
    real(dl), intent(out) :: w, len
    real(dl), intent(in) :: r0, w0
    
    len = 3._dl**0.5 ; w = 0.2_dl
  end subroutine grid_params_

  subroutine profile_params_(p_type,params)
    integer, intent(out) :: p_type
    real(dl), dimension(2) :: params

    p_type = p_guess
    params(1) = beta   ! power
    params(2) = 1._dl  ! effective mass
  end subroutine profile_params_
  
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
