program Instanton
  use constants
  use Cheby
  use Model
  use Nonlinear_Solver
  implicit none

! Nonlinear Solver Parameters and storage.  To be moved to a separate module upon code cleanup
!  type(Field_Model) :: model

  real(dl), dimension(:), allocatable :: xNew
  
  type Solution
     type(Chebyshev) :: tForm
     real(dl), dimension(:), allocatable :: phi
  end type Solution
  
  type(Solver) :: solv
  real(dl), dimension(:), allocatable :: phi, phi_prev
  real(dl), dimension(:), allocatable :: model_params

  integer :: order, n
  real(dl) :: w, len  ! parameters controlling mapping of collocation grid
  real(dl), parameter :: wbase = 10._dl !20._dl
  real(dl) :: meff
  integer :: i

  ! Move these parameters somewhere else
!  type(Bubble_Params) :: bubble
  real(dl) :: phit, phif,r0,w0
  real(dl) :: w_old, l_old
  
  real(dl) :: delta
!  real(dl), dimension(1:18) :: deltas = [ 500., 200., 100., 50., 20., 10., 5., 2., 1., 0.8, 0.6, 0.55, 0.52, 0.51, 0.505, 0.501, 0.5001, 0.50001 ]
!  real(dl), dimension(1:2) :: deltas = [ 0.6, 0.8 ]
  real(dl), dimension(1:25) :: deltas
  
  deltas = [ (-4.+0.25*(i-1), i=25,1,-1) ]
  deltas = 0.5_dl + 10.**deltas
  print*,deltas
  
  ! Values for double well
!  delta = 0.5_dl
!  r0 = 1.5_dl*2._dl**0.5/delta; w0=2.**0.5
!  phif=-1.; phit=1.

!  delta = 0.54_dl
  delta = 0.56_dl
  
  phif = 0._dl; phit = pi
  order = 100; n=order+1
  r0 = 3._dl*(2._dl*delta)**0.5
  meff = (2.*delta)**0.5!*(1._dl-0.25_dl/delta**2)**0.5
  len = r0*3._dl**0.5; w = wbase / len / meff
  
  ! Initialise our derivatives and set up collocation grid
  call create_grid(transform,order,w,len)
!  call create_solver(solv,n,100,0.1_dl)
!  call initialise_equations(transform,delta,3)
!  call create_model(model)

  allocate(phi(0:order),phi_prev(0:order))  ! This is ugly, fix it to allow for varying orders
!  call initialise_fields(phi,delta,.false.)
  ! Move this vacuum solving somewhere else
!  call get_vacuum(phif); call get_vacuum(phit)
!  print*,"vacua are ", phit, phif
  
!  call solve(solv, phi)
!  call output_simple(transform,phi,.true.)
!  call get_action(phi,transform)
!  print*,eval_action(phi,transform,3,r0)

  call compute_profile(delta,order,phi)

  delta = 0.52
  r0 = 3._dl*(2._dl*delta)**0.5
  meff = (2.*delta)**0.5!*(1._dl-0.25_dl/delta**2)**0.5
  call grid_params(w,len,r0,1._dl/meff)

!  len = r0*3._dl**0.5; w = wbase / len / meff

  allocate(xNew(0:order))
  xNew = chebyshev_grid(order,len,w)
  phi_prev = interpolate_instanton(xNew,phi,transform)
  call compute_profile(delta,order,phi,phi_prev)

  delta = 0.51
  r0 = 3._dl*(2._dl*delta)**0.5
  meff = (2.*delta)**0.5!*(1._dl-0.25_dl/delta**2)**0.5
  len = r0*3._dl**0.5; w = wbase / len / meff

  xNew = chebyshev_grid(order,len,w)
  phi_prev = interpolate_instanton(xNew,phi,transform)
  call compute_profile(delta,order,phi,phi_prev)

  delta = 0.505
  r0 = 3._dl*(2._dl*delta)**0.5
  meff = (2.*delta)**0.5!*(1._dl-0.25_dl/delta**2)**0.5
  len = r0*3._dl**0.5; w = wbase / len / meff

  xNew = chebyshev_grid(order,len,w)
  phi_prev = interpolate_instanton(xNew,phi,transform)
  call compute_profile(delta,order,phi,phi_prev)

  delta = 0.501
  r0 = 3._dl*(2._dl*delta)**0.5
  meff = (2.*delta)**0.5!*(1._dl-0.25_dl/delta**2)**0.5
  len = r0*3._dl**0.5; w = wbase / len / meff

  xNew = chebyshev_grid(order,len,w)
  phi_prev = interpolate_instanton(xNew,phi,transform)
  call compute_profile(delta,order,phi,phi_prev)

  delta = 0.5001
  r0 = 3._dl*(2._dl*delta)**0.5
  meff = (2.*delta)**0.5!*(1._dl-0.25_dl/delta**2)**0.5
  len = r0*3._dl**0.5; w = wbase / len / meff

  xNew = chebyshev_grid(order,len,w)
  phi_prev = interpolate_instanton(xNew,phi,transform)
  call compute_profile(delta,order,phi,phi_prev)

 
  delta = 0.50005
  r0 = 3._dl*(2._dl*delta)**0.5
  meff = (2.*delta)**0.5!*(1._dl-0.25_dl/delta**2)**0.5
  len = r0*3._dl**0.5; w = wbase / len / meff

  xNew = chebyshev_grid(order,len,w)
  phi_prev = interpolate_instanton(xNew,phi,transform)
  call compute_profile(delta,order,phi,phi_prev)
  
  delta = 0.50001
  r0 = 3._dl*(2._dl*delta)**0.5
  meff = (2.*delta)**0.5!*(1._dl-0.25_dl/delta**2)**0.5
  len = r0*3._dl**0.5; w = wbase / len / meff

  xNew = chebyshev_grid(order,len,w)
  phi_prev = interpolate_instanton(xNew,phi,transform)
  call compute_profile(delta,order,phi,phi_prev)

   delta = 0.500005
  r0 = 3._dl*(2._dl*delta)**0.5
  meff = (2.*delta)**0.5!*(1._dl-0.25_dl/delta**2)**0.5
  len = r0*3._dl**0.5; w = wbase / len / meff

  xNew = chebyshev_grid(order,len,w)
  phi_prev = interpolate_instanton(xNew,phi,transform)
  call compute_profile(delta,order,phi,phi_prev)

  delta = 0.500003
  r0 = 3._dl*(2._dl*delta)**0.5
  meff = (2.*delta)**0.5!*(1._dl-0.25_dl/delta**2)**0.5
  len = r0*3._dl**0.5; w = wbase / len / meff

  xNew = chebyshev_grid(order,len,w)
  phi_prev = interpolate_instanton(xNew,phi,transform)
  call compute_profile(delta,order,phi,phi_prev)

 
!  call scan_profiles(deltas)
!  call delete_solver(solv)
  
contains

  subroutine interpolate_new_field(order,phi_old,phi_new,tForm)
    integer, intent(in) :: order
    real(dl), dimension(0:order), intent(in) :: phi_old
    real(dl), dimension(0:order), intent(out) :: phi_new
    type(Chebyshev), intent(in) :: tForm
    real(dl), dimension(0:order) :: xNew
    xNew = chebyshev_grid(order,len,w)
    
  end subroutine interpolate_new_field
  
  ! GRRR, this thing is all fucked because of the order I need to compute things in.  Design a useful way to modularise it
  
  !>@brief
  !> Initialise our initial guess for the instanton profile
  !>
  !>@param[in] prev  If True - initialise profile using previous numerical profile
  !>                 If False - initialise using an analytic formula
  subroutine initialise_fields(phi_i,delta,prev,w,len,phi_prev)
    real(dl), dimension(:), intent(out) :: phi_i
    real(dl), intent(in) :: delta
    logical, intent(in) :: prev
    real(dl), intent(in), optional :: w,len
    real(dl), intent(in), optional :: phi_prev
    real(dl) :: r0, meff
    real(dl) :: phi_new(1:size(phi_i))
    real(dl), allocatable :: xNew(:)

    ! This is all fucked
    if (prev) then
!       allocate(xNew(0:order))
!       xNew = chebyshev_grid()
!       phi_new = interpolate_instanton(xNew,phi_prev,transform,len,w)
!       phi_i = phi_new
       ! deallocate(xNew)
    else
       r0 = 3._dl*(2.*delta)**0.5
       meff = (2.*delta)**0.5*(1._dl-0.25_dl/delta**2)**0.5
!       if (delta > 0.7) then
!          phi = -2._dl*atan(exp((transform%xGrid-r0)*meff)) + pi ! For Drummond
!       else
!          phi = -0.5_dl*(phit-phif)*tanh((transform%xGrid-r0)/w0) + 0.5_dl*(phit+phif)
!       endif
       !       call thin_wall_profile()
       if (meff*r0 < 100.) then
          phi_i = -2._dl*atan(-0.5_dl*exp(meff*r0)/cosh(meff*transform%xGrid)) ! Brutal nonlocality
       else
          phi_i = -2._dl*atan(exp((transform%xGrid-r0)*meff)) + pi
       endif
    endif
  end subroutine initialise_fields

  
  subroutine initial_profile(phi,delta,phi_prev)
    type(Solution), intent(out) :: phi, phi_prev
    real(dl), intent(in) :: delta
    real(dl) :: r0, meff

    r0 = 3._dl*(delta*2._dl)**0.5
    meff = (delta*2._dl)**0.5*(1._dl-0.25/delta**2)**0.5
!    if (meff*r0 < 1.) then
    if (.false.) then
    ! Add interpolation and condition above
    else
       call breather_profile(phi%tForm%xGrid,phi%phi,r0,1./meff)  ! ugly nonlocality with transform
    endif
       
  end subroutine initial_profile
  
  subroutine scan_profiles(delVals)
    real(dl), dimension(:), intent(in) :: delVals
    
    open(unit=80,file='integrals.dat')
    do i=1,size(delVals)
       delta = delVals(i)
       r0 = 3._dl*2._dl**0.5*delta**0.5; w0=2.**0.5;
       meff = (2._dl*delta)**0.5*(1._dl-0.25_dl/delta**2)**0.5
       len=r0*3._dl**0.5; w = wbase/len/meff
       call destroy_chebyshev(transform)
       call create_grid(transform,order,w,len)
       call initialise_equations(transform,delta,3)
       !     if (meff*r0 < 2._dl) then
       if (.false.) then
          print*,"using prev, for delta = ", delta
          phi = interpolate_instanton(transform%xGrid,phi_prev,transform) ! Need to store old transform
       else
          call initialise_fields(phi,delta,.false.)
       endif
       call solve(solv, phi)
       !call print_solver(solv)
       !call output_simple(transform,phi,.false.)
       print*,delta, eval_action(phi,transform,3,r0)
       write(80,*) delta, eval_action(phi,transform,3,r0)
       w_old = w; l_old = len; phi_prev = phi
    enddo
  end subroutine scan_profiles
  
  !>@todo
  !>@arg Write this.  It will have identical functionality to the compute_profile subroutine, just with a function return.  For easiest maintainability, probably easiest to just call the subroutine
  function instanton_profile(delta,order,phi_init) result(phi)
    real(dl), intent(in) :: delta
    integer, intent(in) :: order
    real(dl), dimension(0:order), optional :: phi_init
    real(dl), dimension(0:order) :: phi
  end function instanton_profile
  
  !>@brief
  !> Solve for a single bubble profile.  Optionally, store the output field values, derivatives, potential, etc. in a file.
  subroutine compute_profile(delta, order, phi, phi_init)
    real(dl), intent(in) :: delta
    integer, intent(in) :: order
    real(dl), dimension(0:order) :: phi
    real(dl), intent(in), optional :: phi_init(0:order)
    
    type(Solver) :: solv
    real(dl) :: w,len         ! Collocation grid parameters
    real(dl) :: r0, w0, meff  ! Bubble Parameters
    real(dl) :: phif, phit    ! Move this elsewhere
    integer :: n
    
    n = order+1

!   call get_minima(phif,phit)
    phif = 0._dl; phit = pi  ! Replace this with the function call above

!    r0 = 3._dl*2._dl**0.5*delta**0.5
!    meff = (2.*delta)**0.5!*(1._dl-0.25_dl/delta**2)**0.5
    call bubble_parameters(delta,r0,meff)
    call grid_params(w,len,r0,1._dl/meff)

    ! The way transform is used here is ugly and nonlocal.  Fix it!!!
    call destroy_chebyshev(transform)
    call create_grid(transform,order,w,len)
    call create_solver(solv,n,100,0.1_dl)
    call initialise_equations(transform,delta,3)

    ! Modularise this part
    if (present(phi_init)) then
       phi(0:order) = phi_init(0:order)
    else
       call initialise_fields(phi,delta,.false.)  ! Replace this with a call to thin-wall profile or something
    endif
    call solve(solv,phi)
    print*,eval_action(phi,transform,3,r0)
    call output_simple(transform,phi,.false.)
  end subroutine compute_profile

  !>@brief
  !> Compute the thin wall profile.  This will be more easily stored in the model subroutine
  subroutine thin_wall_profile(rvals,phi,r0,w0,phit,phif)
    real(dl), dimension(:), intent(in) :: rvals
    real(dl), dimension(:), intent(out) :: phi
    real(dl), intent(in) :: r0, w0, phit, phif
    
    !phi = -0.5_dl*(phit-phif)*tanh((rvals-r0)/w0) + 0.5_dl*(phit+phif)
    phi = -2._dl*atan(exp((rvals-r0)/w0)) + pi
  end subroutine thin_wall_profile

  subroutine breather_profile(rvals,phi,r0,w0)
    real(dl), dimension(:), intent(in) :: rvals
    real(dl), dimension(:), intent(out) :: phi
    real(dl), intent(in) :: r0,w0
    real(dl) :: meff

    meff = 1._dl/20
    phi = -2._dl*atan(-0.5_dl*exp(meff*r0)/cosh(meff*rvals))
  end subroutine breather_profile
  
  !>@brief
  !> Compute the radius and width of the thin-wall using the tension, vacuum splitting and number of dimensions
  subroutine thin_wall_params(rinit,width,sigma,drho,dim)
    real(dl), intent(out) :: rinit, width
    real(dl), intent(in) :: sigma, drho
    integer, intent(in) :: dim

    rinit = dble(dim)*sigma / drho
    width = 0._dl  ! Make this actually work
  end subroutine thin_wall_params

  subroutine bubble_parameters(delta,r0,meff)
    real(dl), intent(in) :: delta
    real(dl), intent(out) :: r0, meff

    meff = (2._dl*delta)**0.5!*(1._dl-0.25/delta**2)**0.5
    r0 = 3._dl*(2._dl*delta)**0.5
!    w0 = 1._dl/(2._dl*delta)**0.5  ! Adjust as necessary
  end subroutine bubble_parameters
  
  subroutine grid_params(w,len,r0,w0)
    real(dl), intent(out) :: w, len
    real(dl), intent(in) :: r0,w0
    real(dl), parameter :: wscl = 10._dl
    len = r0*3._dl**0.5
    w = wscl * w0 / len 
  end subroutine grid_params

    
  subroutine bubble_params(delta,r0,w0,len,w)
    real(dl), intent(in) :: delta
    real(dl), intent(out) :: r0, w0, len, w
    real(dl) :: meff
    meff = (2._dl*delta)**0.5*(1._dl-0.25_dl/delta**2)**0.5
    r0 = 3._dl*(2._dl*delta)**0.5
    w0 = 1._dl/meff  ! fix this.  Depends on delta in my writing of the potential

    len = r0*3._dl**0.5
    w = w0/len  ! Adjust this
  end subroutine bubble_params
  
  subroutine create_grid(tForm,ord,w,l)
    type(Chebyshev), intent(out) :: tForm
    integer, intent(in) :: ord
    real(dl), intent(in) :: w,l

    call create_chebyshev(tForm,ord,2,.false.,.true.)
    call transform_to_evens(tForm)
    call cluster_points(tForm,w,.true.)
    call transform_double_infinite(tForm,l)
  end subroutine create_grid

  subroutine interpolate_phi(phi_new, xNew, phi_old, xOld)
    real(dl), dimension(1:), intent(out) :: phi_new
    real(dl), dimension(1:), intent(in) :: xNew, phi_old, xOld
    real(dl), dimension(:), allocatable :: phi_spec
    real(dl), dimension(:,:), allocatable :: basis
    integer :: i, ord

    ord = size(phi_old)-1
    allocate( phi_spec(0:ord) ); allocate(basis(0:ord,0:2))
    phi_spec = matmul(transform%fTrans, phi_old)
    do i=1,size(xNew)
       call evaluate_chebyshev(ord,xNew(i),basis,0)
       phi_new(i) = sum(phi_spec(:)*basis(:,0))
    enddo

    deallocate(phi_spec)
  end subroutine interpolate_phi

  !>@brief
  !> Get the eigenvalues around the instanton solution on the original collocation grid
  subroutine get_evalues(phi,n,init)
    real(dl), dimension(:), intent(in) :: phi
    integer, intent(in) :: n
    logical, intent(in) :: init

    double precision, dimension(1:n,1:n) :: L
    double precision, dimension(1:n) :: eval_real, eval_imag
    integer :: ierror
    double precision :: dummy(1:1)
    real, allocatable, dimension(:), save :: work
    integer, save :: iwork
    integer :: i, imax(1)

    ! Allocate workspace
    if (init) then
       call DGEEV('N','N',n,L,n,eval_real,eval_imag,dummy,1,dummy,1,dummy,-1,ierror)
       if (ierror == 0) then
          iwork = int(dummy(1))
          allocate(work(1:iwork))
       else
          print*,"Error allocating workspace for eigenvalue problem, exiting"
          stop
       endif
    endif

    call variation(phi,L)  ! Change this if different boundary conditions are needed.  Currently is will barf since it's getting the wrong boundaries for the eigenvalue problem

    call DGEEV('N','N',n,L,n,eval_real,eval_imag,dummy,1,dummy,1,work,iwork,ierror)
    if (ierror /= 0) then
       print*,"Error in eigenvalue solution"
       stop
    endif
   
  end subroutine get_evalues

  subroutine get_evectors()
  end subroutine get_evectors
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Delete the stuff up until the end of delete
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_vacuum(fld)
    real(dl), intent(inout) :: fld

    integer, parameter :: maxit=16
    real(dl), parameter :: min_tol=1.e-14
    real(dl) :: vp,vpp,dfld
    integer :: l

    print*,"initial field is ",fld
    do l=1,maxit
       vpp = vdprime(fld)
       vp = vprime(fld)
       dfld = -vp/vpp
       fld = fld + dfld
       print*,"new field ",fld
       if (abs(dfld) < min_tol) exit
    enddo
    
    if (l.eq.maxit) then
       print*,"Failed to find local minimum of potential. Adust guess"
       stop
    endif
    print*,"Vacuum is ",fld," derivative is ",vprime(fld)
  end subroutine get_vacuum


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End Delete
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_simple(tForm,phi,init)
    type(Chebyshev), intent(in) :: tForm
    real(dl), dimension(1:tForm%ord+1), intent(in) :: phi
    logical, intent(in) :: init
    integer :: i, sz
    integer, parameter :: u = 60  ! Fix this ugliness
    real(dl), dimension(1:size(phi)) :: phi_spec, dphi, d2phi
    logical :: o
    
    inquire(opened=o,unit=u)
    if (.not.o) open(unit=u,file='instanton.dat')

    sz = size(phi)
    phi_spec = matmul(tForm%fTrans,phi) 
    dphi = matmul(tForm%derivs(:,:,1),phi)
    d2phi = matmul(tForm%derivs(:,:,2),phi)
    do i=1,sz
       write(u,*) tForm%xGrid(i-1), tForm%weights(i-1), phi(i), potential(phi(i)) - potential(phif), vprime(phi(i)), vdprime(phi(i)), dphi(i), d2phi(i), phi_spec(i), tForm%wFunc(i-1)
    enddo
    write(u,*)
  end subroutine output_simple

  
  !>@todo
  !>@arg Modify this to allow for more than one field being solved for
  subroutine output(phi,tForm,sol,init)
    real(dl), dimension(1:), intent(in) :: phi
    type(Chebyshev), intent(in) :: tForm
    type(Solver), intent(in) :: sol
    logical, intent(in) :: init
    integer :: u, i, sz
    real(dl), dimension(1:size(phi)) :: phi_spec, dphi, d2phi

    sz = size(phi)
    phi_spec = matmul(tForm%fTrans,phi)
    dphi = matmul(tForm%derivs(:,:,1),phi)
    d2phi = matmul(tForm%derivs(:,:,2),phi)
    if (init) then
       u = 98
       open(unit=u,file='instanton-full.dat')
    endif

    do i=1,sz
       write(u,*) tForm%xGrid(i-1), phi(i), potential(phi(i))-potential(phif), vprime(phi(i)), vdprime(phi(i)), dphi(i), d2phi(i), phi_spec(i), tForm%weights(i-1), tForm%wFunc(i-1)
    enddo
    write(u,*)
  end subroutine output
    
  !>@brief
  !> Compute the instanton action given the transform grid and field values
  subroutine get_action(fld,transform)
    real(dl), dimension(:), intent(in) :: fld
    type(Chebyshev), intent(in) :: transform
    real(dl), dimension(1:size(fld)) :: dfld, lag
    real(dl) :: action, tension, KE, PE
    
#ifdef USE_BLAS
    call DGEMV()
#else
    dfld = matmul(transform%derivs(:,:,1),fld)
#endif
    lag = 0.5_dl*dfld**2 + potential(fld) - potential(phif)
    action = quadrature(transform,lag*transform%xGrid(:)**3)
    tension = quadrature(transform,dfld**2)
    KE = quadrature(transform,0.5_dl*dfld**2*transform%xGrid(:)**3)
    PE = quadrature(transform,(potential(fld)-potential(phif))*transform%xGrid(:)**3)
    
    print*,"integrals are",action,tension,KE,PE,0.5_dl*tension*r0**3
    write(80,*) delta, action, tension, KE, PE, tension*r0**3
  end subroutine get_action

#define V_TW(f) ( del*sin(f)**2 )  ! This isn't going to work anymore since del is private
  function new_action(sol,d,r0) result(action)
    real(dl), dimension(1:7) :: action
    type(Solution), intent(in) :: sol
    integer, intent(in) :: d
    real(dl), intent(in) :: r0
    real(dl), dimension(size(sol%phi)) :: dfld, lag

    dfld = matmul(sol%tForm%derivs(:,:,1),sol%phi)
    lag = 0.5_dl*dfld**2 + potential(sol%phi) - potential(phif)
    action = 0._dl
  end function new_action
   
  function eval_action(fld,transform,d,r0) result(action)
    real(dl), dimension(1:7) :: action
    real(dl), dimension(:), intent(in) :: fld
    type(Chebyshev), intent(in) :: transform
    integer, intent(in) :: d
    real(dl), intent(in) :: r0
    real(dl), dimension(1:size(fld)) :: dfld, lag

    dfld = matmul(transform%derivs(:,:,1),fld)
    lag = 0.5_dl*dfld**2 + potential(fld) - potential(phif)
    action(1) = quadrature(transform,lag*transform%xGrid(:)**d)
    action(2) = quadrature(transform,dfld**2)
    action(3) = quadrature(transform,0.5_dl*dfld**2*transform%xGrid(:)**d)
    action(4) = quadrature(transform,(potential(fld)-potential(phif))*transform%xGrid(:)**d)
    action(5) = 0.5_dl*action(2)*r0**d
    action(6) = quadrature(transform,(0.5_dl*dfld**2+potential_tw(fld))*transform%xGrid(:)**d)
    action(7) = quadrature(transform,fld*vprime(fld)*transform%xGrid(:)**d)
    ! Add in the "thin-wall potential" here
  end function eval_action

  ! This doesn't seem to be working yet
  function interpolate_instanton(r_new,f_cur,tForm) result(f_int)
    real(dl), dimension(:), intent(in) :: r_new, f_cur
    type(Chebyshev), intent(in) :: tForm
    real(dl) :: L,w
    real(dl), dimension(1:size(r_new)) :: f_int
    real(dl) :: xvals(1:size(r_new)), spec(1:size(f_cur)),  bVals(1:size(f_cur),0:2)
    integer :: i; real(dl) :: winv

    w = tForm%scl; L = tForm%len
    winv = 1._dl/w
    xvals = r_new / sqrt(r_new**2+L**2)
    xvals = atan(winv*tan(pi*(xvals-0.5_dl)))/pi + 0.5_dl
    xvals = 2._dl*xvals**2-1._dl

#ifdef USEBLAS

#else
    spec = matmul(tForm%fTrans,f_cur)
#endif
    do i=1,size(r_new)
       call evaluate_chebyshev(tForm%ord,xvals(i),bVals,2)
       f_int(i) = sum(spec*bVals(:,0))
    enddo
  end function interpolate_instanton

  !>@brief
  !> Return the collocation grid for even Chebyshevs on the doubly-infinite interval with clustering
  function chebyshev_grid(order,L,w) result(xVals)
    integer, intent(in) :: order
    real(dl), intent(in) :: L,w
    real(dl), dimension(0:order) :: xVals
    integer :: i; real(dl) :: dkcol

    ! Replace with a function call
    dkcol = 0.5_dl*pi / dble(order+1)
    do i=0,order
       xVals(i) = -cos( dble(2*i+1)*dkcol )
    enddo
    xVals = sqrt(0.5_dl*xVals + 0.5_dl)
    xVals = (1._dl/pi)*atan(w*tan(pi*(xVals-0.5_dl))) + 0.5_dl
    xVals = L*xVals / sqrt(1._dl - xVals**2)
  end function chebyshev_grid
  
end program instanton
