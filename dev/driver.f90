program Instanton_Solve
  use constants
  use Utils, only : newunit
  use Cheby
  use Model
  use Nonlinear_Solver
  implicit none
  
  real(dl), dimension(:), allocatable :: phi
  real(dl), dimension(:), allocatable :: model_params
  
  real(dl), dimension(:), allocatable :: x_uniform, phi_int
  integer, parameter :: n_un = 2048
  real(dl) :: len_un
  
  integer :: order_
  integer :: i, u

  real(dl) :: phit, phif  ! kill these to restore locality, currently used in output subroutine
  
  real(dl) :: delta
  real(dl), allocatable :: deltas(:); integer :: nDel
  integer :: dim
  type(Chebyshev) :: tForm

! Values for double well
!  delta = 0.5_dl
!  r0 = 1.5_dl*2._dl**0.5/delta; w0=2.**0.5
!  phif=-1.; phit=1.
  
  dim = 1
  nDel = 55; allocate(deltas(1:nDel))
  deltas = 0.5_dl + 10.**([ (-7.+0.2*(i-1), i=size(deltas),1,-1) ])

  ! For BEC paper
!  nDel = 6; allocate(deltas(1:nDel))
!  deltas = 0.5_dl*([ 1.2, 1.3, 1.4, 1.5, 1.55, 1.6 ])**2
  
  call get_minima(phif,phit)
  order_=100
  allocate(phi(0:order_))
  call compute_profile(0.5*1.2_dl**2,order_,dim,phi,tForm,out=.true.)
  
  len_un = 50._dl*2._dl**0.5  ! Fix this
  allocate(x_uniform(1:n_un), phi_int(1:n_un))
  x_uniform = [ ( (i-1)*(len_un/n_un), i=1,n_un ) ]
  
  phi_int = interpolate_instanton(x_uniform,phi,tForm)
  open(unit=newunit(u),file='instanton_interp.dat')
  do i=n_un,2,-1
     write(u,*) -x_uniform(i), phi_int(i)
  enddo
  do i=1,n_un
     write(u,*) x_uniform(i), phi_int(i)
  enddo
  close(u)
    
  !call scan_profiles(deltas,dim,200)
  
contains

  ! TO DO : Finish writing this
  subroutine sim_bubble_profile(phi,tForm,n,len,rc)
    real(dl), dimension(:), intent(in) :: phi
    type(Chebyshev), intent(in) :: tForm
    real(dl), intent(in) :: len, rc
    integer, intent(in) :: n

    real(dl) :: dx
    real(dl), dimension(1:n) :: xNew, rad, phi_new
    integer :: i, u

    ! Fix for odd number of grid sites
    dx = len / dble(n)
    xNew = (/ ( (i-1-n/2)*dx, i=1,n) /);
    rad = (/ ( abs(xNew(i)-rc), i=1,n ) /)
    
    open(unit=newunit(u),file='init_bubble.dat')
    phi_new = interpolate_instanton(rad,phi,tForm)
    do i=1,n
       write(u,*) phi_new(i), 0._dl
    enddo
    close(unit=u)
  end subroutine sim_bubble_profile
  
  subroutine get_minima(phif,phit)
    real(dl), intent(out) :: phif, phit
    phif = pi; phit = 0._dl
  end subroutine get_minima

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

    meff = 1._dl/20  ! Why is this hardcoded?
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

  subroutine bubble_parameters_nd(delta,dim,r0,meff)
    real(dl), intent(in) :: delta, dim
    real(dl), intent(out) :: r0, meff

    meff = (2._dl*delta)**0.5*(1._dl-0.25_dl/delta**2)**0.5
    r0 = dim*(2._dl*delta)**0.5
  end subroutine bubble_parameters_nd
    
  subroutine grid_params(w,len,r0,w0)
    real(dl), intent(out) :: w, len
    real(dl), intent(in) :: r0,w0
    real(dl), parameter :: wscl = 5._dl
    len = r0*3._dl**0.5
    w = wscl * w0 / len
    if (w0 > r0) then
       len = w0*3._dl**0.5
       w = 1._dl
    endif
  end subroutine grid_params
      
  subroutine create_grid(tForm,ord,w,l)
    type(Chebyshev), intent(out) :: tForm
    integer, intent(in) :: ord
    real(dl), intent(in) :: w,l

    call create_chebyshev(tForm,ord,2,.false.,.true.)
    call transform_to_evens(tForm)
    call cluster_points(tForm,w,.true.)
    call transform_double_infinite(tForm,l)
  end subroutine create_grid
  
  ! GRRR, this thing is all fucked because of the order I need to compute things in.  Design a useful way to modularise it
  
  !>@brief
  !> Initialise our initial guess for the instanton profile
  !>
  !>@param[in] prev  If True - initialise profile using previous numerical profile
  !>                 If False - initialise using an analytic formula
  subroutine initialise_fields(phi_i,delta,prev,tForm,w,len,phi_prev)
    real(dl), dimension(:), intent(out) :: phi_i
    real(dl), intent(in) :: delta
    logical, intent(in) :: prev
    type(Chebyshev), intent(in) :: tForm
    real(dl), intent(in), optional :: w,len
    real(dl), intent(in), optional :: phi_prev
    real(dl) :: r0, meff
    real(dl) :: phi_new(1:size(phi_i))
    real(dl), allocatable :: xNew(:)

    ! This is all fucked
    if (prev) then
    else
!       r0 = dim*(2.*delta)**0.5
       !       meff = (2.*delta)**0.5*(1._dl-0.25_dl/delta**2)**0.5
       call bubble_parameters_nd(delta,dim*1._dl,r0,meff)
       if (meff*r0 < 100.) then
!          phi_i = -2._dl*atan(-0.5_dl*exp(meff*r0)/cosh(meff*tForm%xGrid))  ! other minima
          phi_i = 2._dl*atan(-0.5_dl*exp(meff*r0)/cosh(meff*tForm%xGrid)) + pi
       else
          phi_i = 2._dl*atan(exp((tForm%xGrid-r0)*meff))
!          phi_i = -2._dl*atan(exp((tForm%xGrid-r0)*meff)) + pi  ! other minima
       endif
    endif
  end subroutine initialise_fields

#ifdef CLASS
  subroutine initial_profile(phi,delta,phi_prev)
    type(Instanton), intent(out) :: phi, phi_prev
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

  !>@brief
  !> Obtain instanton profiles and actions over the prescribed range of potential parameters
  subroutine scan_profiles_(delVals,dim,ord)
    real(dl), dimension(:), intent(in) :: delVals
    integer, intent(in) :: dim, ord

    type(Instanton) :: inst
    real(dl) :: delta

    real(dl) :: w, len, r0, meff
    real(dl) :: phit, phif
    real(dl), dimension(0:ord) :: xNew, phi_prev, phi
    integer :: i,u
    
    open(unit=newunit(u),file='actions.dat')
    call create_instanton(inst,ord,dim)

    delta = delVals(1)
    call compute_profile_(inst,delta,out=.false.)
    write(u,*) delta, compute_action_(inst)
    
    do i=2,size(delVals)
       delta = delVals(i)
!       call get_minima(phif,phit)  
       if ( abs(phi(0)-phif) < 0.9*abs(phif-phit) ) then
          call bubble_parameters_nd(delta,dim*1._dl,r0,meff); call grid_params(w,len,r0,1._dl/meff)
          ! Add initialisatin using previous profile
       else
          call compute_profile_(inst,delta,out=.false.)
       endif
       write(u,*) delta, compute_action_(inst)
    enddo
    close(u)
  end subroutine scan_profiles_
#endif
  
!!!! TO FIX: nonlocality of w,r0,len,meff in this subroutine
  !!!! Extra calls to create chebyshev grids
  subroutine scan_profiles(delVals,dim,ord)
    real(dl), dimension(:), intent(in) :: delVals
    integer, intent(in) :: dim, ord

    type(Chebyshev) :: tForm
    real(dl) :: delta
    real(dl) :: w, len, r0, meff
    real(dl) :: phit, phif
    real(dl), dimension(0:ord) :: xNew, phi_prev, phi
    integer :: actFile

    call get_minima(phif,phit)  ! In general this depends on the potential parameters
    open(unit=newunit(actFile),file='actions.dat')

    delta = delVals(1)
    call compute_profile(delta,ord,dim,phi,tForm,out=.true.)
    write(actFile,*) delta, eval_action(phi,tForm,dim,r0)
    
    do i=2,size(delVals)
       delta = delVals(i)
       print*,i,"delta = ",delta, abs(phi(0)-phif)/pi
       if ( abs(phi(0)-phif) < 0.9*abs(phif-phit)) then
          call bubble_parameters_nd(delta,dim*1._dl,r0,meff); call grid_params(w,len,r0,1._dl/meff)
          print*,i," prev profile", r0, meff
          xNew = chebyshev_grid(ord,len,w)
          phi_prev = interpolate_instanton(xNew,phi,tForm)
          call compute_profile(delta,ord,dim,phi,tForm,phi_prev,out=.true.)
       else
          call compute_profile(delta,ord,dim,phi,tForm,out=.true.)
       endif
       write(actFile,*) delta, eval_action(phi,tForm,dim,r0)
    enddo
    close(actFile)
  end subroutine scan_profiles
  
  !>@brief
  !> Solve for a single bubble profile.  Optionally, store the output field values, derivatives, potential, etc. in a file.
  subroutine compute_profile(delta, order, dim, phi, tForm, phi_init, out)
    real(dl), intent(in) :: delta
    integer, intent(in) :: order, dim
    real(dl), dimension(0:order), intent(out) :: phi
    type(Chebyshev), intent(out) :: tForm
    real(dl), intent(in), optional :: phi_init(0:order)
    logical, intent(in), optional :: out

    logical :: outLoc
    type(Solver) :: solv
    real(dl) :: w,len         ! Collocation grid parameters
    real(dl) :: r0, meff  ! Bubble Parameters
    real(dl) :: phif, phit    ! Move this elsewhere
    integer :: n

    outLoc = .false.; if (present(out)) outLoc = out
    n = order+1

    call get_minima(phif,phit)
    
    call bubble_parameters_nd(delta,dim*1._dl,r0,meff)
    call grid_params(w,len,r0,1._dl/meff)

    ! The way transform is used here is ugly and nonlocal.  Fix it!!!
!    call destroy_chebyshev(transform)
    call create_grid(tForm,order,w,len)
    call create_solver(solv,n,100,0.1_dl)
    call initialise_equations(tForm,delta,dim)
    
    ! Modularise this part
    if (present(phi_init)) then
       phi(0:order) = phi_init(0:order)
    else
!!! Need to fix this now, since I have a different meff in the subroutine
       call initialise_fields(phi,delta,.false.,tForm)  ! Replace this with a call to thin-wall profile or something
    endif
    call solve(solv,phi)

    if (outLoc) call output_simple(tForm,phi,.false.)
  end subroutine compute_profile

  !! This is super wasteful if we are going to interpolate onto the same grid many times
  !! If we are, should store the interpolation matrix separately and just do matrix multiplies
  subroutine interpolate_phi(phi_new, xNew, phi_old, xOld, tForm)
    real(dl), dimension(1:), intent(out) :: phi_new
    real(dl), dimension(1:), intent(in) :: xNew, phi_old, xOld
    type(Chebyshev), intent(in) :: tForm
    real(dl), dimension(:), allocatable :: phi_spec
    real(dl), dimension(:,:), allocatable :: basis
    integer :: i, ord

    ord = size(phi_old)-1
    allocate( phi_spec(0:ord) ); allocate(basis(0:ord,0:2))
    phi_spec = matmul(tForm%fTrans, phi_old)
    do i=1,size(xNew)
       call evaluate_chebyshev(ord,xNew(i),basis,0)
       phi_new(i) = sum(phi_spec(:)*basis(:,0))
    enddo

    deallocate(phi_spec)
  end subroutine interpolate_phi

  ! This doesn't seem to be working yet
  function interpolate_instanton(r_new,f_cur,tForm) result(f_int)
    real(dl), dimension(:), intent(in) :: r_new, f_cur
    type(Chebyshev), intent(in) :: tForm
    real(dl), dimension(1:size(r_new)) :: f_int

    real(dl) :: L,w
    real(dl) :: xvals(1:size(r_new)), spec(1:size(f_cur)),  bVals(1:size(f_cur),0:2)
    integer :: i; real(dl) :: winv

    w = tForm%scl; L = tForm%len
    winv = 1._dl/w
    xvals = r_new / sqrt(r_new**2+L**2)
    xvals = atan(winv*tan(pi*(xvals-0.5_dl)))/pi + 0.5_dl
    xvals = 2._dl*xvals**2-1._dl
#ifdef USEBLAS
    call DGEMV
#else
    spec = matmul(tForm%fTrans,f_cur)
#endif
    do i=1,size(r_new)
       call evaluate_chebyshev(tForm%ord,xvals(i),bVals,2)
       f_int(i) = sum(spec*bVals(:,0))
    enddo
  end function interpolate_instanton
  
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

    call variation(phi,L)  ! Change this if different boundary conditions are needed.  Currently this will barf since it's getting the wrong boundaries for the eigenvalue problem

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
       print*,"Failed to find local minimum of potential. Adjust guess"
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

    print*,"Outputting profile"
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
       
  function eval_action(fld,tForm,d,r0) result(action)
    real(dl), dimension(1:7) :: action
    real(dl), dimension(:), intent(in) :: fld
    type(Chebyshev), intent(in) :: tForm
    integer, intent(in) :: d
    real(dl), intent(in) :: r0
    real(dl), dimension(1:size(fld)) :: dfld, lag

    dfld = matmul(tForm%derivs(:,:,1),fld)
    lag = 0.5_dl*dfld**2 + potential(fld) - potential(phif)
    action(1) = quadrature(tForm,lag*tForm%xGrid(:)**d)
    action(2) = quadrature(tForm,dfld**2)
    action(3) = quadrature(tForm,0.5_dl*dfld**2*tForm%xGrid(:)**d)
    action(4) = quadrature(tForm,(potential(fld)-potential(phif))*tForm%xGrid(:)**d)
    action(5) = 0.5_dl*action(2)*r0**d
    action(6) = quadrature(tForm,(0.5_dl*dfld**2+potential_tw(fld))*tForm%xGrid(:)**d)
    action(7) = quadrature(tForm,fld*vprime(fld)*tForm%xGrid(:)**d)
    ! Add in the "thin-wall potential" here
  end function eval_action

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

  !>@brief
  !> Return Gauss-Lobatto Grid for the even Chebyshevs
  function chebyshev_grid_lobatto(order,L,w) result(xVals)
    integer, intent(in) :: order
    real(dl), intent(in) :: L,w
    real(dl), dimension(0:order) :: xVals
    integer :: i; real(dl) :: dkcol

    dkcol = pi / dble(order)
    do i=0,order
       xVals(i) = -dcos(dble(i)*dkcol)
    enddo
    xVals = sqrt(0.5_dl*xVals + 0.5_dl)
    xVals = (1._dl/pi)*atan(w*tan(pi*(xVals-0.5_dl))) + 0.5_dl
    xVals = L*xVals / sqrt(1._dl - xVals**2)
  end function chebyshev_grid_lobatto
  
end program Instanton_Solve
