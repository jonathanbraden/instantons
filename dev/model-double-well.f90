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
#define POTENTIAL(f) ( 0.25_dl*(f**2-1._dl)**2 + del*(f**3/3._dl-f-2._dl/3._dl) )
#define VPRIME(f) ( (f+del)*(f**2-1._dl) )
#define VDPRIME(f) ( 3._dl*f**2 - 1._dl + 2._dl*del*f )

module Model
  use constants
  use Cheby
  implicit none
  private :: ndim, del, i_
  
  real(dl) :: del, ndim
  integer, parameter :: nPar = 1
  
  integer :: i_
  integer, parameter :: nEps = 160
  real(dl), parameter :: eps_min=0._dl, eps_max=1._dl
  real(dl), parameter :: log_a_min = -9._dl, log_a_max = 7._dl
  real(dl), parameter, dimension(nEps) :: aVals = 10.**(/ (log_a_min+(i_-1)*(log_a_max-log_a_min)/dble(nEps), i_=nEps,1,-1) /)
  real(dl), parameter, dimension(nEps) :: scanVals = 1._dl/(1._dl+aVals)

  integer, parameter :: p_guess = 2 ! Using 3 gives a weird nonconvergence at delta=0.11
  integer, parameter :: p_guess_2 = 1
  
contains

  subroutine set_model_params(params,dim)
    real(dl), dimension(1:nPar), intent(in) :: params
    real(dl), intent(in) :: dim
    del = params(1); ndim = dim
  end subroutine set_model_params

  function get_model_params() result(params)
    real(dl), dimension(1:nPar) :: params
    params(1) = del
  end function get_model_params
  
  subroutine get_minima(phif,phit)
    real(dl), intent(out) :: phif, phit
    phif = -1._dl; phit = 1._dl
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

  real(dl) function m2eff_max(params) result(m2)
    real(dl), dimension(1:nPar), intent(in), optional :: params

    m2 = -1._dl + del**2
    if (present(params)) m2 = -1._dl + params(1)**2
  end function m2eff_max

  real(dl) function m2eff_fv(params) result(m2)
    real(dl), dimension(1:nPar), intent(in), optional :: params

    m2 = 2._dl*(1._dl+del)
    if (present(params)) m2 = 2._dl*(1._dl+params(1)) 
  end function m2eff_fv

  real(dl) function phi_max(params) result(phi)
    real(dl), dimension(1:nPar), intent(in), optional :: params

    phi = -del
    if (present(params)) phi = -params(1)
  end function phi_max

!!!!
!!!! TO DO : This must be implemented
!!!!
  real(dl) function phi_out(params) result(phi)
    real(dl), dimension(1:nPar), intent(in), optional :: params

    phi = 0._dl
  end function phi_out
  
  function get_profile_type(params,dim) result(p)
    real(dl), dimension(1:nPar), intent(in) :: params
    real(dl), intent(in) :: dim
    integer :: p

    real(dl) :: r0, meff
    call bubble_parameters_nd_(params(1),dim,r0,meff)
    if (r0 > 100.) then
       p = 2
    else
       p = 1
    endif
  end function get_profile_type
  
  subroutine get_guess_params(params,dim,params_ic,prof_type)
    real(dl), intent(in) :: params(1:nPar), dim
    real(dl), dimension(1:4), intent(out) :: params_ic
    integer, intent(out), optional :: prof_type  ! Remove optional nature
    real(dl) :: delta, n_dim

    n_dim = dim; delta = params(1)
    call bubble_parameters_nd_(delta,n_dim,params_ic(1),params_ic(2))
    call get_minima(params_ic(3),params_ic(4))
    if (present(prof_type)) then
       prof_type = p_guess
       if (params_ic(2)*params_ic(1) < 10._dl) prof_type = p_guess_2
    endif
  end subroutine get_guess_params
  
  subroutine get_grid_params(params,dim,len,scl)
    real(dl), dimension(1:nPar), intent(in) :: params
    real(dl), intent(in) :: dim
    real(dl), intent(out) :: len, scl
    real(dl) :: meff_max, r0
    real(dl) :: rad_scl
    real(dl), parameter :: wscl = 0.5_dl*twopi !5.8_dl !5.8_dl*3.**0.5

    rad_scl = 1._dl ! sqrt(3._dl) 
    call bubble_parameters_nd_(params(1),dim,r0,meff_max)
    len = r0*rad_scl; scl = wscl / (meff_max*len) ! *2.**0.5
    if (meff_max*r0 < 1._dl) then
       len = rad_scl/meff_max
       scl = 1._dl
    endif
  end subroutine get_grid_params

  logical function use_previous(params,dim) result(prev)
    real(dl), intent(in) :: params(1:nPar), dim
    real(dl) :: r0, meff
    call bubble_parameters_nd_(params(1),dim,r0,meff)
    prev = (r0*meff < 2.)
  end function use_previous
  
  !>@brief
  !> Given specified radius and width of a bubble profile, adjust grid mapping parameters.
  !>
  !> The relationship between the radius and mapping length are fixed by choice of polynomials
  !> Should probably be moved into the chebyshev class
  subroutine grid_params_(w,len,r0,w0)
    real(dl), intent(out) :: w, len
    real(dl), intent(in) :: r0, w0
    real(dl), parameter :: wscl = 8.96_dl   ! decent for cubic, need to tweak delta -> 1 part
    
    len = r0*3._dl**0.5
    w = wscl * w0 / r0
    if (w0 > r0) then
       len = w0*3._dl**0.5
       w = 1._dl
    endif
  end subroutine grid_params_
  
  ! These need to be adjusted for very model.  Might be worth moving it
  ! Change delta to parameters for the model
  subroutine bubble_parameters_nd_(delta,dim,r0,meff)
    real(dl), intent(in) :: delta, dim
    real(dl), intent(out) :: r0, meff

    meff = sqrt(1._dl - delta**2) ! effective mass at maximum
    r0 = dim / (sqrt(2._dl)*delta)
  end subroutine bubble_parameters_nd_

  !>@brief
  !> The unperturbed part of the potential used in computing the thin-wall domain wall
  elemental function potential_tw(phi)
    real(dl) :: potential_tw
    real(dl), intent(in) :: phi
    potential_tw = 0.25_dl*(phi**2-1._dl)**2
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
