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
  private :: eps, ndim, i_

  integer, parameter :: nPar = 1
  real(dl) :: eps, ndim
  real(dl), parameter :: eps_min = 0._dl, eps_max = exp(0.5_dl)

  integer :: i_
  integer, parameter :: nEps = 60
  real(dl), parameter :: log_e_min = -6._dl, log_e_split = -1._dl
  real(dl), parameter :: scanVals(nEps) = 10.**(/ ( log_e_min + (i_-1)*(log_e_split-log_e_min)/dble(nEps), i_=1,nEps ) /)  ! Scan parameters piling up at zero

  integer, parameter :: nEps_r = 50
  real(dl), parameter :: log_de_min = -5._dl, log_de_max = log10(eps_max - 10.**(log_e_split))
  real(dl), parameter :: scanVals_r(nEps_r) = eps_max - 10.**(/ (log_de_min + (i_-1)*(log_de_max-log_de_min)/dble(nEps_r-1), i_=nEps_r,1,-1) /) 
  
  integer, parameter :: p_guess = 4
  
contains
  
  subroutine set_model_params(params,dim)
    real(dl), intent(in) :: params(1:nPar), dim
    eps = params(1); ndim = dim
  end subroutine set_model_params
  
  subroutine grid_params_(w,len,r0,w0)
    real(dl), intent(out) :: w, len
    real(dl), intent(in) :: r0, w0

    len = 10._dl; w = 1._dl  ! Muck around with these
  end subroutine grid_params_

  subroutine bubble_parameters_nd_(delta,dim,r0,meff)
    real(dl), intent(in) :: delta, dim
    real(dl), intent(out) :: r0, meff

    !meff = sqrt(2._dl)  ! epsilon = 0 result
    !meff = sqrt(-vdprime(get_maximum()))
    print*,"meff is ",meff
    r0 = 1._dl/meff
  end subroutine bubble_parameters_nd_

  subroutine get_minima(phif,phit)
    real(dl), intent(out) :: phif, phit
    phif = 0._dl; phit = exp(0.5_dl*2._dl)  ! Adjust this
  end subroutine get_minima

  !>@brief
  !> Return point on other side of the barrier where the energy is zero
  real(dl) function phi_out()
    phi_out = sqrt(exp(1.)-eps**2)
  end function phi_out
  
  !>@brief
  !> Returns the maximum of the potential by solving
  !> \f[
  !>  (x^2+e^2)log(x^2+e^2) = e^2
  !> \f]
  real(dl) function phi_max(params,phi0,tol) result(phi)
    real(dl), intent(in) :: params(1:nPar)
    real(dl), intent(in), optional :: phi0
    real(dl), intent(in), optional :: tol

    real(dl) :: tol_, eps, dp

    eps = params(1)
    tol_ = 1.e-12; if (present(tol)) tol_ = tol
    
    phi = 1._dl; if (present(phi0)) phi = phi0
    dp = -( (phi+eps**2)*log(phi+eps**2) - eps**2 ) / ( log(phi+eps**2) + 1._dl)

    do while ( (abs(dp) >tol_) .or. (abs((phi+eps**2)*log(phi+eps**2)-eps**2) > tol_) )
       phi = phi + dp
       dp = -( (phi+eps**2)*log(phi+eps**2) - eps**2 ) / ( log(phi+eps**2) + 1._dl)
    enddo
    phi = sqrt(phi)
  end function phi_max

  real(dl) function meff_max(params) result(meff)
    real(dl), intent(in) :: params(1:nPar)
    real(dl) :: phi

    phi = phi_max(params,phi0=1.,tol=1.e-12)
    meff = sqrt(-vdprime(phi))
  end function meff_max
  
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

  subroutine get_grid_params(params,dim,len,scl)
    real(dl), intent(in) :: params(1:nPar), dim
    real(dl), intent(out) :: len, scl

    real(dl) :: rad_scl, meff

    rad_scl = sqrt(3._dl)
    meff = meff_max(params)
    len = 2._dl*1._dl*rad_scl/meff
    scl = 1._dl
  end subroutine get_grid_params

  subroutine get_guess_params(params,dim,params_ic,prof_type)
    real(dl), dimension(1:nPar), intent(in) :: params
    real(dl), intent(in) :: dim
    real(dl), dimension(1:4), intent(out) :: params_ic
    integer, intent(out), optional :: prof_type

    params_ic(1) = 1._dl
    params_ic(2) = (0.5_dl)**0.5
    params_ic(3) = 0._dl
    params_ic(4) = exp(0.5_dl*(dim+1))
    if (present(prof_type)) prof_type = p_guess
  end subroutine get_guess_params

  logical function use_previous(params,dim) result(prev)
    real(dl), intent(in) :: params(1:nPar), dim

    prev = (params(1) > 0.2_dl)
  end function use_previous
  
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
