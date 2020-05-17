module Instanton_Class
  use constants, only : dl
  use Utils, only : newunit
  use Cheby
  use Model  ! Try to remove this dependence
  use Nonlinear_Solver
  implicit none
  
  type Instanton
     type(Chebyshev) :: tForm
     real(dl), dimension(:), allocatable :: phi
     integer :: ord  
     integer :: dim  ! Generalize for this to be noninteger
     real(dl) :: l, w  ! Specifies grid
     real(dl) :: r0, meff, phif, phit  ! Specifies initial bubble
     logical :: exists = .false.
     integer :: unit = -1
  end type Instanton
  
contains

  ! To do in here : Pass in the potential, its derivative and it's second derivative as functions inside the type
  
  ! Expand to noninteger dimensions
  subroutine create_instanton(this,ord,d)
    type(Instanton), intent(out) :: this
    integer, intent(in) :: ord,d
    this%dim = d; this%ord = ord
    if (allocated(this%phi)) deallocate(this%phi) ! Remove this to only allocate if size has changed
    allocate(this%phi(0:ord))
    this%exists = .true.
    this%unit = 56    !! Fix this to be automated
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

    ! This should all be moved somewhere more modular where it doesn't rely on knowledge of how grid was transformed
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
    
    real(dl), dimension(0:this%ord) :: phi_spec, dphi, d2phi
    logical :: o
    integer :: i, sz
    real(dl) :: phif

    integer :: u  ! Fix this to automatically select an open unit

    ! Move this to autoselect the unit (can I just set u = -1 initially and this will work?)
    u = this%unit
    inquire(opened=o,unit=u)
    if (.not.o) open(unit=u,file='instanton_.dat')

    sz = size(this%phi); phif = this%phi(sz-1)
    
    phi_spec = matmul(this%tForm%fTrans,this%phi) 
    dphi = matmul(this%tForm%derivs(:,:,1),this%phi)
    d2phi = matmul(this%tForm%derivs(:,:,2),this%phi)
    do i=0,sz-1
       write(u,*) this%tForm%xGrid(i), this%tForm%weights(i), this%phi(i), potential(this%phi(i)) - potential(phif), vprime(this%phi(i)), vdprime(this%phi(i)), dphi(i), d2phi(i), phi_spec(i), this%tForm%wFunc(i)
    enddo
    write(u,*)    
  end subroutine output_instanton

  function compute_action_(this) result(action)
    type(Instanton), intent(in) :: this
    real(dl), dimension(1:8) :: action
    real(dl), dimension(0:this%ord) :: dfld, lag
    real(dl) :: r0, d, phif
    integer :: ord

    ord = size(this%phi)-1
    d = dble(this%dim); r0 = this%r0  ! Fix r0 not being set
    phif = this%phi(ord)  ! fix this
    
    dfld = matmul(this%tForm%derivs(:,:,1),this%phi)
    lag = 0.5_dl*dfld**2 + potential(this%phi) - potential(phif)

    action(1) = quadrature(this%tForm,lag*this%tForm%xGrid(:)**d)
    action(2) = quadrature(this%tForm,0.5_dl*dfld**2*this%tForm%xGrid(:)**d)
    action(3) = quadrature(this%tForm,(potential(this%phi)-potential(phif))*this%tForm%xGrid(:)**d)
    action(4) = quadrature(this%tForm,dfld**2)
    action(5) = 0.5_dl*action(4)*r0**d
    action(6) = quadrature(this%tForm,(0.5_dl*dfld**2+potential_tw(this%phi))*this%tForm%xGrid(:)**d)
    action(7) = quadrature(this%tForm,this%phi*vprime(this%phi)*this%tForm%xGrid(:)**d)
    action(8) = quadrature(this%tForm, this%tForm%xGrid(:)**(d-1)*lag ) ! Why this index?
  end function compute_action_

  !>@brief
  !> Compute the instanton solution for the symmetry breaking parameter delta
  !>
  !>@param[inout] this
  !>@param[in] delta
  !>@param[in] (optional) phi_init
  !>@param[in] (optional) out  - Boolean to write result to file or not
  !>@param[in] (optional) p_i  - Integer choice of initial analytic profile guess
  subroutine compute_profile_(this,params,phi_init,out,p_i)
    type(Instanton), intent(inout) :: this
    real(dl), dimension(:), intent(in) :: params
    real(dl), intent(in), optional :: phi_init(0:this%ord)
    logical, intent(in), optional :: out
    integer, intent(in), optional :: p_i

    logical :: outLoc; integer :: order, dim, n, p_loc
    type(Solver) :: solv

    ! Clean up all this extraneous crap
    real(dl) :: w, len      ! These seem extraneous
    real(dl) :: r0, meff    ! These seem extraneous
    real(dl) :: phif, phit  ! These also do
    integer :: i 
    integer :: u
    
    dim = this%dim; order = this%ord
    outLoc = .false.; if (present(out)) outLoc = out
    p_loc = 1; if (present(p_i)) p_loc = p_i
    n = order+1

    ! The first two calls here are just getting params for the initial guess
    call get_minima(phif,phit)  ! Change this to pass in params
    call bubble_parameters_nd_(params(1),dim*1._dl,r0,meff)  ! Change this to pass in params
    call grid_params_(w,len,r0,1._dl/meff)

    call create_grid_(this%tForm,order,w,len) ! Replace this with the library call
    call create_solver(solv,n,100,0.1_dl)
    call initialise_equations(this%tForm,params,dim)
    
    if (present(phi_init)) then
       this%phi(0:order) = phi_init(0:order)
    else
       call profile_guess(this,r0,meff,phif,phit,p_loc)
       !this%phi = 1._dl / (1._dl + this%tForm%xGrid**2/params(1))**(0.9*params(1))  ! This is useful for the the Fubini
    endif
    
    open(unit=newunit(u),file='ic.dat')
    do i=0,this%ord
       write(u,*) this%tForm%xGrid(i), this%phi(i), potential(this%phi(i)), vprime(this%phi(i)), vdprime(this%phi(i))
    enddo
    close(unit=u)
    
    call solve(solv,this%phi)

    if (outLoc) call output_instanton(this)
  end subroutine compute_profile_

!!!! This functionality should be moved into the chebyshev code
!!! I'm pretty sure it's in there already, so just kill this and use the call in the library
  subroutine create_grid_(tForm,ord,w,l)
    type(Chebyshev), intent(out) :: tForm
    integer, intent(in) :: ord
    real(dl), intent(in) :: w,l

    call create_chebyshev(tForm,ord,2,.false.,.true.)
    call transform_to_evens(tForm)
    call cluster_points(tForm,w,.true.)
    call transform_double_infinite(tForm,l)
  end subroutine create_grid_
  
  !>@brief
  !> Initialise our initial profile guess based on the given radius and width.
  !> A choice to use tanh, arctan, breather, Gaussian or generalised witchhat initial profiles are given
  subroutine profile_guess(this,r0,meff,phif,phit,s)
    type(Instanton), intent(inout) :: this
    real(dl), intent(in) :: r0,meff,phif, phit
    integer, intent(in) :: s ! Select type of profile to use
    
    select case (s)
    case (1)
       this%phi(:) = breather_p(this%tForm%xGrid(:),r0,meff,phif,phit)
    case (2)
       this%phi(:) = tanh_p(this%tForm%xGrid(:),r0,meff,phif,phit)
    case (3)
       this%phi(:) = atan_p(this%tForm%xGrid(:),r0,meff,phif,phit)
    case (4)
       this%phi(:) = gaussian_p(this%tForm%xGrid(:),r0,meff,phif,phit)
    case (5)
       this%phi(:) = witchhat_p(this%tForm%xGrid(:),r0,meff,phif,phit)
    case default
       this%phi(:) = breather_p(this%tForm%xGrid(:),r0,meff,phif,phit)
    end select
  end subroutine profile_guess

  elemental function breather_p(x,r0,m,phif,phit) result(f)
    real(dl), intent(in) :: x
    real(dl), intent(in) :: r0,m,phif,phit
    real(dl) :: f
    f = (phif-phit)*(2._dl/pi)*atan(-0.5_dl*exp(m*r0)/cosh(m*x)) + phif
  end function breather_p

  elemental function tanh_p(x,r0,m,phif,phit) result(f)
    real(dl), intent(in) :: x
    real(dl), intent(in) :: r0,m,phif,phit
    real(dl)  :: f
    f = 0.5_dl*(phif-phit)*tanh(m*(x-r0)) + 0.5_dl*(phit+phif)
  end function tanh_p

  elemental function atan_p(x,r0,m,phif,phit) result(f)
    real(dl), intent(in) :: x
    real(dl), intent(in) :: r0,m,phif,phit
    real(dl) :: f

!    real(dl) :: mscl
!    mscl = 1.5_dl*m
    f = (phif-phit)*(2._dl/pi)*atan(exp(m*(x-r0))) + phit
  end function atan_p

  elemental function gaussian_p(x,r0,m,phif,phit) result(f)
    real(dl), intent(in) :: x
    real(dl), intent(in) :: r0,m,phif,phit
    real(dl) :: f
    f = (phit-phif)*exp(-(m*x)**2) + phif
  end function gaussian_p

  elemental function witchhat_p(x,r0,m,phif,phit) result(f)
    real(dl), intent(in) :: x
    real(dl), intent(in) :: r0,m,phif,phit
    real(dl) :: f
    f = (phit-phif)/(1.+m*x**2/r0)**r0 + phif
  end function witchhat_p
    
  subroutine breather_profile(x,f,r0,m,phif,phit)
    real(dl), dimension(:), intent(in) :: x
    real(dl), dimension(1:size(x)), intent(out) :: f
    real(dl), intent(in) :: r0,m,phif,phit
    f(:) =  (phif-phit)*(2._dl/pi)*atan(-0.5_dl*exp(m*r0)/cosh(m*x)) + phif
  end subroutine breather_profile

  subroutine tanh_profile(x,f,r0,m,phif,phit)
    real(dl), dimension(:), intent(in) :: x
    real(dl), dimension(:), intent(out) :: f
    real(dl), intent(in) :: r0,m,phif,phit

    f(:) = 0.5_dl*(phif-phit)*tanh(m*(x(:)-r0)) + 0.5_dl*(phit+phif)
  end subroutine tanh_profile

  subroutine atan_profile(x,f,r0,m,phif,phit)
    real(dl), dimension(:), intent(in) :: x
    real(dl), intent(in) :: r0,m,phif,phit
    real(dl), dimension(1:size(x)), intent(out) :: f

    f(:) = (phif-phit)*(2._dl/pi)*atan(exp(1.5*m*(x(:)-r0))) + phit
  end subroutine atan_profile

  !>@brief
  !> Use Newton method to find FV minimum.
  !> Required input is initial guess, and function for potential derivative and curvature
  !>
  !> TO DO : Fix this for the case where the curvature at the minimum is very small
  subroutine solve_for_fv(phifv,tol_deriv_,tol_diff_)
    real(dl), intent(inout) :: phifv
    real(dl), intent(in), optional :: tol_diff_,tol_deriv_
    
    real(dl) :: phi_prev, tol_diff, tol_deriv
    integer :: i
    
    tol_deriv = 1.e-12; if (present(tol_deriv_)) tol_deriv = tol_deriv_
    tol_diff = tol_deriv; if (present(tol_diff_)) tol_diff = tol_diff_
    phi_prev = phifv + 10.*tol_diff

    i = 0
    do while ( abs(phi_prev-phifv) > tol_diff .and. vdprime(phifv) > tol_deriv )
       phi_prev = phifv
       phifv = phifv - vprime(phifv)/vdprime(phifv)
       i = i+1  ! Add error checking on max number of iterations
    enddo
  end subroutine solve_for_fv
  
end module Instanton_Class
