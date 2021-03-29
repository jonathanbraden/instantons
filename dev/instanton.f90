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
     real(dl) :: dim
     logical, dimension(1:2) :: bc
     character(8) :: grid_type
     real(dl) :: l, w  ! Specifies grid
     real(dl) :: r0, meff, phif, phit  ! Specifies initial bubble (can remove these)
     logical :: exists = .false.
     integer :: fNum = -1
  end type Instanton
  
contains

  ! To do in here : Pass in the potential, its derivative and it's second derivative as functions inside the type
  
  subroutine create_instanton(this,ord,dim)
    type(Instanton), intent(out) :: this
    integer, intent(in) :: ord
    real(dl), intent(in) :: dim
    
    this%dim = dim; this%ord = ord
    if (allocated(this%phi)) deallocate(this%phi)
    allocate(this%phi(0:ord))
    this%exists = .true.
    this%fNum = 56   !! Fix this to be automated
  end subroutine create_instanton

  subroutine destroy_instanton(this)
    type(Instanton), intent(inout) :: this
    if (.not.this%exists) return
    deallocate(this%phi)
    this%exists = .false.
  end subroutine destroy_instanton

!!!!!
!!!! This needs to be updated with the new options for the grid mappings
!!!!
  function interpolate_instanton_(this,r_new) result(f_int)
    type(Instanton), intent(in) :: this
    real(dl), dimension(:), intent(in) :: r_new
    real(dl), dimension(1:size(r_new)) :: f_int

    real(dl) :: L,w
    real(dl) :: xvals(1:size(r_new)), spec(1:this%tForm%ord+1), bVals(1:this%tForm%ord+1,0:2)
    integer :: i; real(dl) :: winv

    ! This should all be moved somewhere more modular where it doesn't rely on knowledge of how grid was transformed (I've started doing this)

    ! call invert_mapping(this,r_new,xvals)  ! Fix to actual call, add rational mappings
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
    integer :: u

    u = this%fNum
    inquire(opened=o,file='instanton.dat')
    if (.not.o) open(unit=u,file='instanton.dat')

    sz = size(this%phi); phif = this%phi(sz-1) ! This isn't happy for the log potential, possibly from huge error in calculation of potential at origin.
    
    phi_spec = matmul(this%tForm%fTrans,this%phi) 
    dphi = matmul(this%tForm%derivs(:,:,1),this%phi)
    d2phi = matmul(this%tForm%derivs(:,:,2),this%phi)
    do i=0,sz-1
       write(u,*) this%tForm%xGrid(i), this%phi(i), dphi(i), d2phi(i), potential(this%phi(i)) - potential(phif), vprime(this%phi(i)), vdprime(this%phi(i)), phi_spec(i), this%tForm%weights(i), this%tForm%wFunc(i)
    enddo
    write(u,*)    
  end subroutine output_instanton

  function compute_action(this) result(action)
    type(Instanton), intent(in) :: this
    real(dl), dimension(1:8) :: action
    real(dl), dimension(0:this%ord) :: dfld, lag
    real(dl) :: r0, d, phif
    integer :: ord

    ord = size(this%phi)-1
    d = this%dim
    phif = this%phi(ord)  ! fix this
    r0 = this%r0  ! Fix this in case r0 isn't specified
    
    dfld = matmul(this%tForm%derivs(:,:,1),this%phi)
    lag = 0.5_dl*dfld**2 + potential(this%phi) - potential(phif)

    action(1) = quadrature(this%tForm,lag*this%tForm%xGrid(:)**d)
    action(2) = quadrature(this%tForm,0.5_dl*dfld**2*this%tForm%xGrid(:)**d)
    action(3) = quadrature(this%tForm,(potential(this%phi)-potential(phif))*this%tForm%xGrid(:)**d)
    action(4) = quadrature(this%tForm,dfld**2)
    action(5) = 0.5_dl*action(4)*r0**d
    action(6) = quadrature(this%tForm,(0.5_dl*dfld**2+potential_tw(this%phi))*this%tForm%xGrid(:)**d)
    action(7) = quadrature(this%tForm,this%phi*vprime(this%phi)*this%tForm%xGrid(:)**d)
    action(8) = quadrature(this%tForm, this%tForm%xGrid(:)**(d-1)*lag )
  end function compute_action

  
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
  
  subroutine compute_profile(this,params,get_grid,get_ic,phi_init,grid_type,out)
    type(Instanton), intent(inout) :: this
    real(dl), dimension(:), intent(in) :: params
    procedure(get_grid_params) :: get_grid
    procedure(get_guess_params) :: get_ic
    real(dl), intent(in), optional :: phi_init(0:this%ord)
    character(8), intent(in), optional :: grid_type
    logical, intent(in), optional :: out

    character(8) :: grid_; logical :: out_
    type(Solver) :: solv; integer, parameter :: maxit = 100; real(dl), parameter :: kick = 0.1_dl
    integer :: order, n
    real(dl) :: dim, len, w, pow
    real(dl), dimension(1:4) :: params_ic; integer :: prof_type
    
    grid_ = 'FULL_MID'; if (present(grid_type)) grid_=grid_type
    out_ = .true.; if (present(out)) out_ = out
    
    dim = this%dim; order = this%ord; n = order+1
    call set_model_params(params,dim)
    call get_grid(params,dim,len,w)
!    call create_instanton_grid(this,grid_,len,w,1.)
    call create_instanton_grid(this,grid_,len,w)

    ! Now initialise the fields
    ! Clean this up so it isn't so ugly.  Pass in the guess_params function
    ! call initial_guess
    if (present(phi_init)) then
       this%phi = phi_init
    else
       call get_ic(params,dim,params_ic,prof_type)
       call profile_guess(this,params_ic(1),params_ic(2),params_ic(3),params_ic(4),prof_type)
    endif
    
    call create_solver(solv,n,maxit,kick)
    call initialise_equations(this%tForm,params,dim,this%bc)
    call solve(solv,this%phi)  
    if (out_) call output_instanton(this)
  end subroutine compute_profile
  
  !>@brief
  !> Compute the instanton solution for the symmetry breaking parameter delta
  !>
  !>@param[inout] this
  !>@param[in] delta
  !>@param[in] (optional) phi_init
  !>@param[in] (optional) (len,w) parameters for the grid.  Otherwise uses a function call
  !>@param[in] (optional) out  - Boolean to write result to file or not
  !>@param[in] (optional) p_i  - Integer choice of initial analytic profile guess
  subroutine compute_profile_(this,params,phi_init,grid,out,p_i,prof_params)
    type(Instanton), intent(inout) :: this
    real(dl), dimension(:), intent(in) :: params
    real(dl), intent(in), optional :: phi_init(0:this%ord)
    real(dl), dimension(1:2), intent(in), optional :: grid
    logical, intent(in), optional :: out
    integer, intent(in), optional :: p_i
    real(dl), dimension(:), intent(in), optional :: prof_params

    logical :: outLoc; integer :: order, n, p_loc
    real(dl) :: dim
    type(Solver) :: solv
    real(dl) :: w, len, params_ic(1:4)
    integer :: i 
    
    dim = this%dim; order = this%ord
    outLoc = .false.; if (present(out)) outLoc = out
    p_loc = 1; if (present(p_i)) p_loc = p_i
    n = order+1
    
    call set_model_params(params,dim)
    if (present(grid)) then
       len = grid(1); w = grid(2)
    else
       call get_grid_params(params,dim,len,w)
    endif
    print*,"grid params ",params(1),len,w
    call create_instanton_grid(this,'FULL_MID',len,w)
    ! Better IR behaviour for Fubini (add this functionality to new call)
    !call create_grid_stretch(this%tForm,order,len,20._dl)

    if (present(phi_init)) then
       this%phi(0:order) = phi_init(0:order)
    else
       print*,"Making guess"
       call get_guess_params(params,dim,params_ic)
       print*,params_ic
       call profile_guess(this,params_ic(1),params_ic(2),params_ic(3),params_ic(4),p_loc)
    endif

    print*,"making solver"
    call create_solver(solv,n,100,0.1_dl)
    call initialise_equations(this%tForm,params,dim,this%bc)
    call solve(solv,this%phi)

    ! Add option to try a new profile if the solver doesn't converge
    ! if ( (.not.present(phi_init)) .and. (.not.convg) )
    !   call profile_guess(this,r0,meff,phif,phit,p_loc_2
    !
    
    if (outLoc) call output_instanton(this)
    call delete_solver(solv)
  end subroutine compute_profile_

  subroutine solve_profile(this,params)
    type(Instanton), intent(inout) :: this
    real(dl), dimension(1:nPar), intent(in) :: params
    
    type(Solver) :: solv
    integer :: n
    integer, parameter :: maxIt = 100
    
    n = this%ord+1
    call create_solver(solv,n,maxIt,0.1_dl)
    call initialise_equations(this%tForm,params,this%dim,this%bc) ! Need to pre call this one to simplify calling signature
    call solve(solv,this%phi)
  end subroutine solve_profile
    
  subroutine create_instanton_grid(this,type,len,w,pow)
    type(Instanton), intent(inout) :: this
    character(8), intent(in) :: type
    real(dl), intent(in) :: len, w
    real(dl), intent(in), optional :: pow

    logical :: rational
    integer :: n_deriv
    logical :: num_inv
    
    n_deriv = 2; num_inv = .false.
    rational=.false.; if (present(pow)) rational=.true.
    
    select case (type)
    case ('FULL_MID') ! Gauss grid with even rationals on double infinite
       this%bc = (/.false.,.false./) 
       call create_chebyshev_full(this%tForm,this%ord,n_deriv,'MIDS',num_inv)
       call transform_to_evens(this%tForm)
       call cluster_points(this%tForm,w,.true.)
       if (rational) then
          call transform_double_infinite_rational(this%tForm,len,pow)
          print*,"Warning, rational mapping not fully implemented for interpolation"
       else
          call transform_double_infinite(this%tForm,len)
          this%grid_type = 'FULL_MID'
       endif
       
    case ('FULL_RHT') ! Radau grid with even rationals on double infinite
       this%bc = (/.false.,.true./) 
       print*,"Warning, FULL_RHT not fully tested"
       call create_chebyshev_full(this%tForm,this%ord,n_deriv,'RGHT',num_inv)
       call transform_to_evens(this%tForm)
       call cluster_points(this%tForm,w,.true.)
       if (rational) then
          call transform_double_infinite_rational(this%tForm,len,pow)
          print*,"Warning, rational mapping not fully implemented for interpolation"
       else
          call transform_double_infinite(this%tForm,len)
          this%grid_type = 'FULL_RHT'
       endif

    case ('HALF_END') ! Lobatto grid with rationals on half infinite
       this%bc = (/.true.,.true./) 
       call create_chebyshev_full(this%tForm,this%ord,n_deriv,'ENDS',num_inv)
       call cluster_points(this%tForm,w,.false.) 
       call transform_semi_infinite(this%tForm,len)
       this%grid_type = 'HALF_END'

    case ('HALF_LFT') ! Radau grid with left endpoint on half infinite
       this%bc = (/.true.,.false./) 
       print*,"Warning, HALF_LFT not fully tested"
       call create_chebyshev_full(this%tForm,this%ord,n_deriv,'LEFT',num_inv)
       call cluster_points(this%tForm,w,.false.)
       call transform_semi_infinite(this%tForm,len)
       this%grid_type = 'HALF_LFT'
       
    case default
       print*,"Invalid choice of collocation grid, defaulting to FULL_MID"
       this%bc = (/.false.,.false./) 
       call create_chebyshev_full(this%tForm,this%ord,n_deriv,'MIDS',num_inv)
       call transform_to_evens(this%tForm)
       call cluster_points(this%tForm,w,.true.)
       call transform_double_infinite(this%tForm,len)
       this%grid_type = 'FULL_MID'
    end select
  end subroutine create_instanton_grid

  !>@brief
  !> Invert the mapping of the grid to map a collection of radial grid points
  !> r_new into collocation points of the original Chebyshev polynomials x_new
  subroutine invert_grid_mapping(this,r_new,xvals)
    type(Instanton), intent(inout) :: this
    real(dl), dimension(:), intent(in) :: r_new
    real(dl), dimension(1:size(r_new)), intent(out) :: xvals
    
    real(dl) :: L, w, pow, winv

    select case (this%grid_type)
    case ('FULL_MID')
       w = this%tForm%scl; L = this%tForm%len
       winv = 1._dl/w
       xvals = r_new / sqrt(r_new**2+L**2)
       xvals = atan(winv*tan(pi*(xvals-0.5_dl)))/pi + 0.5_dl
       xvals = 2._dl*xvals**2 - 1._dl

    case ('FULL_RHT')
       w = this%tForm%scl; L = this%tForm%len
       winv = 1._dl/w
       xvals = r_new / sqrt(r_new**2+L**2)
       xvals = atan(winv*tan(pi*(xvals-0.5_dl)))/pi + 0.5_dl
       xvals = 2._dl*xvals**2 - 1._dl
              
    case ('HALF_END')
       w = this%tForm%scl; L = this%tForm%len
       winv = 1._dl/w
       xvals = (r_new - L)/(r_new + L) !!! Check this
       xvals = atan(winv*tan(pi*xvals))/pi !!! Check this
       print*,"Still testing HALF_END interpolation"

    case ('HALF_LFT')
       print*,"HALF_LFT interpolation not yet tested"
       
    case default
       print*,"Warning, attempting to invert a nonexistant grid mapping"
    end select
  end subroutine invert_grid_mapping
  
  ! As with the above, delete this once I've fully tested the create_instanton_grid subroutine
  subroutine create_grid_stretch(tForm,ord,l,p)
    type(Chebyshev), intent(out) :: tForm
    integer, intent(in) :: ord
    real(dl), intent(in) :: l,p

    call create_chebyshev(tForm,ord,2,.false.,.true.)
    call transform_to_evens(tForm)
    call transform_double_infinite_rational(tForm,l,p)
  end subroutine create_grid_stretch

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
    do while ( abs(phi_prev-phifv) > tol_diff .or. vprime(phifv) > tol_deriv )
       phi_prev = phifv
       phifv = phifv - vprime(phifv)/vdprime(phifv)
       i = i+1  ! Add error checking on max number of iterations
    enddo
  end subroutine solve_for_fv

  subroutine solve_for_phi_out(phi_out,tol_deriv_,tol_diff_)
    real(dl), intent(inout) :: phi_out
    real(dl), intent(in), optional :: tol_diff_, tol_deriv_

    real(dl) :: phi_prev, tol_diff, tol_deriv
    integer :: i
    
    tol_deriv = 1.e-12; if (present(tol_deriv_)) tol_deriv = tol_deriv_
    tol_diff = tol_deriv; if (present(tol_diff_)) tol_diff = tol_diff_

    i = 0; phi_prev = phi_out + 10.*tol_diff
    do while ( abs(phi_prev-phi_out) > tol_diff .or. potential(phi_out) > tol_deriv )
       phi_prev = phi_out
       phi_out = phi_out - potential(phi_out)/vprime(phi_out)
       i = i+1  ! Add error checking on max number of iterations
    enddo
  end subroutine solve_for_phi_out

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
    f = 0.5_dl*(phif-phit)*tanh((m/2.**0.5)*(x-r0)) + 0.5_dl*(phit+phif)
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

  elemental function witchhat_p(x,pow,m,phif,phit) result(f)
    real(dl), intent(in) :: x
    real(dl), intent(in) :: pow,m
    real(dl), intent(in) :: phif,phit
    real(dl) :: f
    f = (phit-phif)/(1.+m**2*x**2/pow)**pow + phif
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

    f(:) = 0.5_dl*(phif-phit)*tanh((m/2.**0.5)*(x(:)-r0)) + 0.5_dl*(phit+phif)
  end subroutine tanh_profile

  subroutine atan_profile(x,f,r0,m,phif,phit)
    real(dl), dimension(:), intent(in) :: x
    real(dl), intent(in) :: r0,m,phif,phit
    real(dl), dimension(1:size(x)), intent(out) :: f

    f(:) = (phif-phit)*(2._dl/pi)*atan(exp(0.5_dl*pi*m*(x(:)-r0))) + phit
  end subroutine atan_profile
  
end module Instanton_Class
