module Instanton_Class
  use constants, only : dl
  use Cheby
  use Model  ! Try to remove this dependence
  use Nonlinear_Solver
  implicit none
  
  type Instanton_Multi
     type(Chebyshev) :: tForm
     integer :: nfld, ord
     real(dl), dimension(:,:), pointer :: phi
     real(dl), dimesnion(:), allocatable, target :: phi_dat
     integer :: dim
     real(dl) :: r0, meff, phif, phit
     logical :: exists = .false.
  end type Instanton_Multi
  
contains

  subroutine create_instanton_multi(this,nf,ord,d)
    type(Instanton_multi), intent(out) :: this
    integer, intent(in) :: ord,d, nf
    this%dim = d; this%ord = ord; this%nfld = nf
    if (allocated(this%phi)) deallocate(this%phi) ! Remove this to only allocate if size has changed
    allocate(this%phi(0:ord,1:nf))  ! This is the wrong way
    this%exists = .true.
  end subroutine create_instanton_multi

  subroutine destroy_instanton_multi(this)
    type(Instanton_multi), intent(inout) :: this
    if (.not.this%exists) return
    deallocate(this%phi)
    this%exists = .false.
  end subroutine destroy_instanton_mulit
  
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

    integer, parameter :: u = 56 ! Fix this to automatically select an open unit
    
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
    real(dl), dimension(1:8) :: action
    real(dl), allocatable, dimension(:) :: dfld, lag
    real(dl) :: r0, d, phif
    integer :: ord

    ord = size(this%phi)-1
    d = dble(this%dim); r0 = this%r0
    allocate( dfld(0:ord), lag(0:ord) )
    phif = this%phi(ord)  ! fix this
    
    dfld = matmul(this%tForm%derivs(:,:,1),this%phi)
    lag = 0.5_dl*dfld**2 + potential(this%phi) - potential(phif)

    action(1) = quadrature(this%tForm,dfld**2)
    action(2) = quadrature(this%tForm,lag*this%tForm%xGrid(:)**d)
    action(3) = quadrature(this%tForm,0.5_dl*dfld**2*this%tForm%xGrid(:)**d)
    action(4) = quadrature(this%tForm,(potential(this%phi)-potential(phif))*this%tForm%xGrid(:)**d)
    action(5) = 0.5_dl*action(2)*r0**d
    action(6) = quadrature(this%tForm,(0.5_dl*dfld**2+potential_tw(this%phi))*this%tForm%xGrid(:)**d)
    action(7) = quadrature(this%tForm,this%phi*vprime(this%phi)*this%tForm%xGrid(:)**d)
    action(8) = quadrature(this%tForm, this%tForm%xGrid(:)**(d-1)*lag )
    
    deallocate(dfld,lag)
  end function compute_action_

  !>@brief
  !> Compute the instanton solution for the symmetry breaking parameter delta
  !>
  !>@param[inout] this
  !>@param[in] delta
  !>@param[in] (optional) phi_init
  !>@param[in] (optional) out  - Boolean to write result to file or not
  !>@param[in] (optional) p_i  - Integer choice of initial analytic profile guess
  subroutine compute_profile_(this,delta,phi_init,out,p_i)
    type(Instanton), intent(inout) :: this
    real(dl), intent(in) :: delta
    real(dl), intent(in), optional :: phi_init(0:this%ord)
    logical, intent(in), optional :: out
    integer, intent(in), optional :: p_i

    logical :: outLoc; integer :: order, dim, n, p_loc
    type(Solver) :: solv

    ! Clean up all this extraneous crap
    real(dl) :: w, len      ! These seem extraneous
    real(dl) :: r0, meff    ! These seem extraneous
    real(dl) :: phif, phit  ! These also do

    dim = this%dim; order = this%ord
    outLoc = .false.; if (present(out)) outLoc = out
    p_loc = 3; if (present(p_i)) p_loc = p_i
    n = order+1

    call get_minima(phif,phit)
    call bubble_parameters_nd_(delta,dim*1._dl,r0,meff)
    call grid_params_(w,len,r0,1._dl/meff)

    call create_grid_(this%tForm,order,w,len) ! Replace this with the library call
    call create_solver(solv,n,100,0.1_dl)
    call initialise_equations(this%tForm,delta,dim)
    
    ! Modularise this part
    if (present(phi_init)) then
       this%phi(0:order) = phi_init(0:order)
    else
       call profile_guess(this,r0,meff,phif,phit,p_loc)
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
  !> A choice to use tanh, arctan or breather initial profiles are given
  subroutine profile_guess(this,r0,meff,phif,phit,s)
    type(Instanton), intent(inout) :: this
    real(dl), intent(in) :: r0,meff,phif, phit
    integer, intent(in) :: s ! Select type of profile to use

    print*,s
    
    select case (s)
    case (1)
       !call breather_profile(this%tForm%xGrid,this%phi,r0,meff,phif,phit)
       this%phi(:) = breather_p(this%tForm%xGrid(:),r0,meff,phif,phit)
    case (2)
       !call tanh_profile(this%tForm%xGrid,this%phi,r0,w)
       this%phi(:) = tanh_p(this%tForm%xGrid(:),r0,meff,phif,phit)
    case (3)
       !call atan_profile(this%tForm%xGrid,this%phi,r0,meff,phif,phit)
       this%phi(:) = atan_p(this%tForm%xGrid(:),r0,meff,phif,phit)
    case default
       ! call breather_profile(this%tForm%xGrid,this%phi,r0,meff,phif,phit)
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

!!!! Model dependent
!!! Move this
  subroutine bubble_parameters_nd_(delta,dim,r0,meff)
    real(dl), intent(in) :: delta, dim
    real(dl), intent(out) :: r0, meff

    meff = (2._dl*delta)**0.5*(1._dl-0.25_dl/delta**2)**0.5
    r0 = dim*(2._dl*delta)**0.5
  end subroutine bubble_parameters_nd_

  subroutine get_minima(phif,phit)
    real(dl), intent(out) :: phif, phit
    phif = pi; phit = 0._dl
  end subroutine get_minima
  
end module Instanton_Class
