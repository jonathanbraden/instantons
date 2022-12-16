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
#define POTENTIAL(f) ( 0.5_dl*( sin(f)**2 + (1._dl/del)*(cos(f)-1._dl)) )
#define VPRIME(f) ( 0.5_dl*( sin(2._dl*f) - (1._dl/del)*sin(f)) )
#define VDPRIME(f) ( 0.5_dl*( 2._dl*cos(2._dl*f) - (1._dl/del)*cos(f)) )

module Model
  use constants
  use Cheby
  implicit none
  private :: del, i_

  integer, parameter :: nPar=1
  real(dl) :: del
  
  real(dl), parameter :: del_min = 0.5_dl
  integer :: i_
  integer, parameter :: nPar_ = 200
  real(dl), parameter :: log_del_min = -15._dl, log_del_max = 10._dl
  real(dl), parameter :: scanVals(1:nPar_) = del_min +  10.**([ (log_del_min+(log_del_max-log_del_min)*(i_-1)/dble(nPar_), i_=nPar_,1,-1) ] )
  integer, parameter :: p_guess = 3  ! Check this
  
contains

  subroutine set_model_params(params,dim)
    real(dl), dimension(1), intent(in) :: params
    real(dl), intent(in) :: dim

    real(dl) :: lam, phiout, phimax, m2max, norm
    del = params(1)
    lam = sqrt(2.*del)
    phiout = acos((2._dl-lam**2)/lam**2); phimax = acos(1./lam**2)
    m2max = (lam**2-1._dl)*(lam**2+1._dl)/lam

    norm = (1._dl/(m2max*phiout**2))
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

  subroutine get_grid_params(params,dim,len,scl)
    real(dl), dimension(1:nPar), intent(in) :: params
    real(dl), intent(in) :: dim
    real(dl), intent(out) :: len, scl

    real(dl), parameter :: wscl = 6.3*1.74
    real(dl) :: meff_max, r0, w0, delta, rad_fac

    delta = params(1); rad_fac = sqrt(3._dl)
    call bubble_parameters_nd_(delta,1._dl*dim,r0,meff_max)
    len = r0*rad_fac  ! Fix this for semi-infinite grid, when there's no
    w0 = 1._dl/meff_max
    scl = wscl * (w0/len) ! Change to (w0/len) and adjust wscl
    if (w0 > r0) then
       len = rad_fac*w0
       scl = 1._dl
    endif
  end subroutine get_grid_params

  subroutine get_guess_params(params,dim,params_ic,prof_type)
    real(dl), dimension(1:nPar), intent(in) :: params
    real(dl), intent(in) :: dim
    real(dl), dimension(1:4), intent(out) :: params_ic
    integer, intent(out), optional :: prof_type

    real(dl) :: delta, n_dim
    delta = params(1); n_dim = dim
    call bubble_parameters_nd_(delta,n_dim,params_ic(1),params_ic(2))
    call get_minima(params_ic(3),params_ic(4))
    if (present(prof_type)) then
       prof_type = p_guess
       !if (params_ic(2)*params_ic(1) < 10._dl) prof_type = 1
    endif
  end subroutine get_guess_params

  logical function use_previous(params,dim) result(prev)
    real(dl), intent(in) :: params(1:nPar), dim
    real(dl) :: r0, meff
    call bubble_parameters_nd_(params(1),dim,r0,meff)
    prev = (r0*meff < 2.)
  end function use_previous
  
  subroutine get_profile_params(params,dim,prof_type,p_params)
    real(dl), dimension(1:nPar), intent(in) :: params
    real(dl), intent(in) :: dim
    integer, intent(out) :: prof_type
    real(dl), dimension(1:4), intent(out) :: p_params
    real(dl) :: phi_inf, phi_out

    call get_minima(p_params(1),p_params(2))
    call bubble_parameters_nd_(params(1),dim,p_params(3),p_params(4))
    prof_type = 1
  end subroutine get_profile_params
  
  !>@brief
  !> Given specified radius and width of a bubble profile, adjust grid mapping parameters.
  !>
  !> The relationship between the radius and mapping length are fixed by choice of polynomials
  !> Should probably be moved into the chebyshev class
  subroutine grid_params_(w,len,r0,w0)
    real(dl), intent(out) :: w, len
    real(dl), intent(in) :: r0, w0
    real(dl), parameter :: wscl = 6.3_dl   ! TWEAK THIS!!!!! ! was 6.3
    
    len = r0*3._dl**0.5
    w = wscl * (w0 / r0)
    if (w0 > r0) then
       len = 3._dl**0.5*w0
       w = 1._dl
    endif
  end subroutine grid_params_
  
  ! Change delta to parameters for the model
  subroutine bubble_parameters_nd_(delta,dim,r0,meff)
    real(dl), intent(in) :: delta, dim
    real(dl), intent(out) :: r0, meff

    ! Fix this (check it now)
    meff = (1._dl-0.25_dl/delta**2)**0.5
    ! Fix this (it's definitely wrong
    r0 = dim * 2._dl*delta
  end subroutine bubble_parameters_nd_

  !>@brief
  !> The unperturbed part of the potential used in computing the thin-wall domain wall
  elemental function potential_tw(phi)
    real(dl) :: potential_tw
    real(dl), intent(in) :: phi
    potential_tw = sin(phi)**2
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
