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
#define POTENTIAL(f) ( cos(f) + del*sin(f)**2 - 1._dl + 0.5_dl*m2*f**2 + 0.25_dl*lam*f**4 )
#define VPRIME(f) ( -sin(f) + del*sin(2._dl*f) + m2*f + lam*f**3 )
#define VDPRIME(f) ( -cos(f) + 2._dl*del*cos(2._dl*f) +m2 + 3._dl*lam*f**2)

module Model
  use constants
  use Cheby
  implicit none
  private :: del, m2, lam, i_

  integer, parameter :: nPar=1
  real(dl) :: del, m2, lam
  
  real(dl), parameter :: del_min = 0.5_dl
  integer :: i_
  integer, parameter :: nPar_ = 200
  real(dl), parameter :: log_del_min = -15._dl, log_del_max = 10._dl
  real(dl), parameter :: scanVals(1:nPar_) = del_min +  10.**([ (log_del_min+(log_del_max-log_del_min)*(i_-1)/dble(nPar_), i_=nPar_,1,-1) ] )
  integer, parameter :: p_guess = 3  ! Check this

  real(dl), parameter :: m2Vals(1:nPar_) = ([ (-(i_-1)*0.42/dble(nPar_), i_=1,nPar_) ])
  
contains

  subroutine set_model_params(params,dim)
    real(dl), dimension(1:nPar), intent(in) :: params
    real(dl), intent(in) :: dim
    del = 0.5_dl*1.2**2
    m2 = params(1)
    lam = -m2/(0.5_dl*twopi)**2
  end subroutine set_model_params

  function get_model_params() result(params)
    real(dl), dimension(1:nPar) :: params
    params(1) = m2
  end function get_model_params

  subroutine write_model_header(fNum)
    integer, intent(in) :: fNum

    write(fNum,*) "# Model potential is:"
    write(fNum,*) "# V(phi) = cos(phi) - 1 + del*sin(phi)**2 + 0.5*m2*phi**2 - 0.25*m2/pi**2*phi**4"
    write(fNum,*) "#"
    write(fNum,*) "# del     m2      action variables"
  end subroutine write_model_header
  
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

    real(dl), parameter :: wscl = 0.5_dl*twopi !6.3*1.74
    real(dl) :: meff_max, r0, w0, delta
    real(dl) :: rad_scl
    
    delta = del; rad_scl = 1._dl !sqrt(3._dl)
    call bubble_parameters_nd_(delta,dim,r0,meff_max)
    len = r0*rad_scl
    w0 = 1._dl/meff_max
    scl = wscl * (w0/len) ! Change to (w0/len) and adjust wscl
    if (w0 > r0) then
       len = rad_scl*w0
       scl = 1._dl
    endif
  end subroutine get_grid_params

  subroutine get_guess_params(params,dim,params_ic,prof_type)
    real(dl), dimension(1:nPar), intent(in) :: params
    real(dl), intent(in) :: dim
    real(dl), dimension(1:4), intent(out) :: params_ic
    integer, intent(out), optional :: prof_type

    real(dl) :: delta, n_dim
    !    delta = params(1); n_dim = dim
    delta = del; n_dim = dim
    call bubble_parameters_nd_(delta,n_dim,params_ic(1),params_ic(2))
    call get_minima(params_ic(3),params_ic(4))
    if (present(prof_type)) then
       prof_type = p_guess
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

end module Model
