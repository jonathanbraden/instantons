module Instanton_Class
  use constants, only : dl
  use Cheby
  use Model  ! Try to remove this dependence
  use Nonlinear_Solver
  implicit none
  
  type Instanton
     type(Chebyshev) :: tForm
     real(dl), dimension(:), allocatable :: phi
     integer :: ord
     integer :: dim
     real(dl) :: r0, meff, w
     logical :: exists = .false.
  end type Instanton
  
contains

  subroutine create_instanton(this,ord,d)
    type(Instanton), intent(out) :: this
    integer, intent(in) :: ord,d
    this%dim = d; this%ord = ord
    if (allocated(this%phi)) deallocate(this%phi) ! Remove this to only allocate if size has changed
    allocate(this%phi(0:ord))
    this%exists = .true.
  end subroutine create_instanton

  subroutine destroy_instanton(this)
    type(Instanton), intent(inout) :: this
    if (.not.this%exists) return
    deallocate(this%phi)
    this%exists = .false.
  end subroutine destroy_instanton
  
  function interpolate_instanton_(this,r_new) result(f_int)
    type(Instanton), intent(in) :: this
    real(dl), dimension(:), intent(in) :: r_new
    real(dl), dimension(1:size(r_new)) :: f_int

    real(dl) :: L,w
    real(dl) :: xvals(1:size(r_new)), spec(1:this%tForm%ord+1), bVals(1:this%tForm%ord+1,0:2)
    integer :: i; real(dl) :: winv
    
    w = this%tForm%scl; L = this%tForm%len
    winv = 1._dl/w
    xvals = r_new / sqrt(r_new**2+L**2)
    xvals = atan(winv*tan(pi*(xvals-0.5_dl)))/pi + 0.5_dl
    xvals = 2._dl*xvals**2 - 1._dl
#ifdef USEBLAS
    spec = 
#else
    spec = matmul(this%tForm%fTrans,this%phi)
#endif
    do i=1,size(r_new)
       call evaluate_chebyshev(this%tForm%ord,xvals(i),bVals,2)
       f_int(i) = sum(spec*bVals(:,0))
    enddo
  end function interpolate_instanton_
  
  !TO DO: phif is extracted from the solution, not from input.
  subroutine output_instanton(this)
    type(Instanton), intent(in) :: this
    
    real(dl), dimension(:), allocatable :: phi_spec, dphi, d2phi
    logical :: o
    integer :: i, sz
    real(dl) :: phif

    integer, parameter :: u = 56
    
    inquire(opened=o,unit=u)
    if (.not.o) open(unit=u,file='instanton_.dat')

    sz = size(this%phi); phif = this%phi(sz-1)
    allocate(phi_spec(0:sz-1),dphi(0:sz-1),d2phi(0:sz-1))
    
    phi_spec = matmul(this%tForm%fTrans,this%phi) 
    dphi = matmul(this%tForm%derivs(:,:,1),this%phi)
    d2phi = matmul(this%tForm%derivs(:,:,2),this%phi)
    do i=0,sz-1
       write(u,*) this%tForm%xGrid(i), this%tForm%weights(i), this%phi(i), potential(this%phi(i)) - potential(phif), vprime(this%phi(i)), vdprime(this%phi(i)), dphi(i), d2phi(i), phi_spec(i), this%tForm%wFunc(i)
    enddo
    write(u,*)
    
    deallocate(phi_spec,dphi,d2phi)
  end subroutine output_instanton

  function compute_action_(this) result(action)
    type(Instanton), intent(in) :: this
    real(dl), dimension(1:7) :: action
    real(dl), allocatable, dimension(:) :: dfld, lag
    real(dl) :: r0, d, phif
    integer :: ord

    ord = size(this%phi)-1
    d = dble(this%dim); r0 = this%r0
    allocate( dfld(0:ord), lag(0:ord) )
    phif = this%phi(ord)  ! fix this
    
    dfld = matmul(this%tForm%derivs(:,:,1),this%phi)
    lag = 0.5_dl*dfld**2 + potential(this%phi) - potential(phif)
    action(1) = quadrature(this%tForm,lag*this%tForm%xGrid(:)**d)
    action(2) = quadrature(this%tForm,dfld**2)
    action(3) = quadrature(this%tForm,0.5_dl*dfld**2*this%tForm%xGrid(:)**d)
    action(4) = quadrature(this%tForm,(potential(this%phi)-potential(phif))*this%tForm%xGrid(:)**d)
    action(5) = 0.5_dl*action(2)*r0**d
    action(6) = quadrature(this%tForm,(0.5_dl*dfld**2+potential_tw(this%phi))*this%tForm%xGrid(:)**d)
    action(7) = quadrature(this%tForm,this%phi*vprime(this%phi)*this%tForm%xGrid(:)**d)

    deallocate(dfld,lag)
  end function compute_action_

  subroutine compute_profile_(this,delta,phi_init,out)
    type(Instanton), intent(inout) :: this
    real(dl), intent(in) :: delta
    real(dl), intent(in), optional :: phi_init(0:this%ord)
    logical, intent(in), optional :: out

    logical :: outLoc; integer :: order, dim, n
    type(Solver) :: solv

    ! Clean up all this extraneous crap
    real(dl) :: w, len      ! These seem extraneous
    real(dl) :: r0, meff    ! These seem extraneous
    real(dl) :: phif, phit  ! These also do

    dim = this%dim; order = this%ord
    outLoc = .false.; if (present(out)) outLoc = out
    n = order+1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FIX THIS
    !   call get_minima(phif,phit)
    phif = 0._dl; phit = pi  ! Replace this with the function call above
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call bubble_parameters_nd_(delta,dim*1._dl,r0,meff)
    call grid_params_(w,len,r0,1._dl/meff)

    call create_grid_(this%tForm,order,w,len) ! Replace this with the library call
    call create_solver(solv,n,100,0.1_dl)
    call initialise_equations(this%tForm,delta,dim)
    
    ! Modularise this part
    if (present(phi_init)) then
       this%phi(0:order) = phi_init(0:order)
    else
!!! Need to fix this now, since I have a different meff in the subroutine
       call initialise_fields_(this%phi,delta,dim,.false.,this%tForm)  ! FIX THIS
    endif
    call solve(solv,this%phi)

    if (outLoc) call output_instanton(this)
  end subroutine compute_profile_

  !>@brief
  !> Given specified radius and width of a bubble profile, adjust grid mapping parameters.
  !>
  !> The relationship between the radius and mapping length are fixed by choice of polynomials
  !> Should probably be moved into the chebyshev class
  subroutine grid_params_(w,len,r0,w0)
    real(dl), intent(out) :: w, len
    real(dl), intent(in) :: r0, w0
    real(dl), parameter :: wscl = 5._dl

    len = r0*3._dl**0.5
    w = wscl * w0 / len
    if (w0 > r0) then
       len = w0*3._dl**0.5
       w = 1._dl
    endif
  end subroutine grid_params_

!!!! This functionality should be moved into the chebyshev code
!!! I'm pretty sure it's in there already, so just kill this an use the call in the library
  subroutine create_grid_(tForm,ord,w,l)
    type(Chebyshev), intent(out) :: tForm
    integer, intent(in) :: ord
    real(dl), intent(in) :: w,l

    call create_chebyshev(tForm,ord,2,.false.,.true.)
    call transform_to_evens(tForm)
    call cluster_points(tForm,w,.true.)
    call transform_double_infinite(tForm,l)
  end subroutine create_grid_

  !!!! Model dependent
  subroutine bubble_parameters_nd_(delta,dim,r0,meff)
    real(dl), intent(in) :: delta, dim
    real(dl), intent(out) :: r0, meff

    meff = (2._dl*delta)**0.5*(1._dl-0.25_dl/delta**2)**0.5
    r0 = dim*(2._dl*delta)**0.5
  end subroutine bubble_parameters_nd_

!!!!!!!!!!!!!!
  ! Clean these next few up so they aren't so ugly
  !!!!!!!!!!!!!
  !>@brief
  !> Initialise our initial guess for the instanton profile
  !>
  !>@param[in] prev  If True - initialise profile using previous numerical profile
  !>                 If False - initialise using an analytic formula
  subroutine initialise_fields_(phi_i,delta,dim,prev,tForm,w,len)
    real(dl), dimension(:), intent(out) :: phi_i
    real(dl), intent(in) :: delta
    integer, intent(in) :: dim
    logical, intent(in) :: prev
    type(Chebyshev), intent(in) :: tForm
    real(dl), intent(in), optional :: w,len
    real(dl) :: r0, meff

    ! This is all fucked
    if (prev) then
    else
       call bubble_parameters_nd_(delta,dim*1._dl,r0,meff)
       if (meff*r0 < 100.) then
!          phi_i = -2._dl*atan(-0.5_dl*exp(meff*r0)/cosh(meff*tForm%xGrid))  ! other minima
          phi_i = 2._dl*atan(-0.5_dl*exp(meff*r0)/cosh(meff*tForm%xGrid)) + pi
       else
          phi_i = 2._dl*atan(exp((tForm%xGrid-r0)*meff))
!          phi_i = -2._dl*atan(exp((tForm%xGrid-r0)*meff)) + pi  ! other minima
       endif
    endif
  end subroutine initialise_fields_

  !>@brief
  !> Initialise our initial profile guess based on the given radius and width.
  !> A choice to use tanh, arctan or breather initial profiles are given
  subroutine profile_guess(this,r0,meff,s)
    type(Instanton), intent(inout) :: this
    real(dl), intent(in) :: r0,meff
    integer, intent(in) :: s ! Select type of profile to use

    !call breather_profile(this%tForm%xGrid,this%phi,r0,w)
    ! this%phi(:) = 2._dl*atan(-0.5_dl*exp(meff*r0)/cosh(meff*this%xGrid(:)) + pi
    !call tanh_profile(this%tForm%xGrid,this%phi,r0,w)
    ! this%phi(:) = 0._dl
    !call atan_profile(this%tForm%xGrid,this%phi,r0,w)
    ! this%phi(:) = 2._dl*atan(exp(meff*(this%tForm%xGrid(:)-r0)))
  end subroutine profile_guess
  
end module Instanton_Class
