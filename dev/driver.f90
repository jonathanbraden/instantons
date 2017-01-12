program Instanton
  use constants
  use Cheby
  use Model
  use Nonlinear_Solver
  implicit none

! Nonlinear Solver Parameters and storage.  To be moved to a separate module upon code cleanup
!  type(Field_Model) :: model

!  type(Chebyshev) :: transform
  type(Solver) :: solv
  real(dl), dimension(:), allocatable :: phi, phi_prev
  real(dl), dimension(:), allocatable :: model_params

  integer :: order, n
  real(dl) :: w, len  ! parameters controlling mapping of collocation grid
  real(dl), parameter :: wbase = 25._dl
  real(dl) :: meff
  integer :: i

  ! Move these parameters somewhere else
  real(dl) :: phit, phif,r0,w0
  
  real(dl) :: delta
  real(dl), dimension(1:20) :: deltas = (/ 0.50001, 0.5001, 0.501, 0.505, 0.51, 0.52, 0.55, 0.6, 0.7, 0.8, 0.9, 1., 2. , 4., 10., 20., 50., 100., 200., 500. /)

  ! Values of double well
!  delta = 0.5_dl
!  r0 = 1.5_dl*2._dl**0.5/delta; w0=2.**0.5
!  phif=-1.; phit=1.

  delta = 4._dl; meff = (2.*delta)**0.5
  r0 = 3._dl*2._dl**0.5*delta**0.5; w0=2.**0.5 ! This w0 should depend on delta since the mass does.  Fix it!!!!  This will also require adjusting the w used on my collocation grid.  easiest to just adjust wbase?
  phif = 0._dl; phit = pi

  len = r0*3._dl**0.5
  
  order = 200; n=order+1
  if (delta > 0.7) then
     w = wbase/len/meff
  else
     w = wbase / len
  endif
  ! Initialise our derivatives and set up collocation grid
  print*,"w is ",w
  call create_grid(transform,order,w,len)

  call create_solver(solv,n,100,0.1_dl)
  call initialise_equations(transform,delta)
!  call create_model(model)

  allocate(phi(0:order),phi_prev(0:order))
  call initialise_fields(.false.)
  call get_vacuum(phif); call get_vacuum(phit)
  print*,"vacua are ", phit, phif
     
  call solve(solv, phi)
  call output_simple(transform%xGrid,phi,.true.)
  call get_action(phi,transform)
  !  call output()

#ifdef FULL
  open(unit=80,file='integrals.dat')
  do i=1,size(deltas)
     ! call bubble_params()
     ! call create_grid()
     ! call initialise_field() ! Add interpolation in here
     delta = deltas(i)
     r0 = 3._dl*2._dl**0.5*delta**0.5; w0=2.**0.5; meff = (2._dl*delta)**0.5
     len=r0*3._dl**0.5; w=wbase / len
     if (delta > 0.7) then
        w = wbase/len/meff
     else
        w = wbase / len
     endif
     print*,"w is ",w
     call destroy_chebyshev(transform)
     call create_grid(transform,order,w,len)
     deallocate(l0) ! This is really ugly.  Fix it up somehow
     call initialise_equations(transform,deltas(i))
     call initialise_fields(.false.)
     call solve(solv, phi)
     !call print_solver(solv)
     !call output_simple(transform%xGrid,phi,.false.)
     !instanton_prev(:) = instanton(:)
     call get_action(phi,transform)
  enddo
#endif
  
  call delete_solver(solv)
  
contains

  subroutine bubble_params(r0,w0,len,w)
    real(dl), intent(out) :: r0, w0, len, w
    r0 = 3._dl*2._dl**0.5*delta**0.5
    w0 = 2.**0.5
    len = r0*3._dl**0.5
    w = 0.5_dl  ! Adjust this
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

    call variation(phi,L)  ! Change this if different boundary conditions are needed

    call DGEEV('N','N',n,L,n,eval_real,eval_imag,dummy,1,dummy,1,work,iwork,ierror)
    if (ierror /= 0) then
       print*,"Error in eigenvalue solution"
       stop
    endif
   
  end subroutine get_evalues

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

  subroutine output_simple(xvals,phi,init)
    real(dl), dimension(1:), intent(in) :: xvals, phi
    logical, intent(in) :: init
    integer :: i, sz
    integer, parameter :: u = 60
    real(dl), dimension(:), allocatable :: phi_spec, dphi, d2phi

    if (init) open(unit=u,file='instanton.dat')

    sz = size(xvals)
    allocate(phi_spec(sz), dphi(sz), d2phi(sz) )
    phi_spec = matmul(transform%fTrans,phi) 
    dphi = matmul(transform%derivs(:,:,1),phi)
    d2phi = matmul(transform%derivs(:,:,2),phi)
    do i=1,sz
       write(u,*) xvals(i), transform%weights(i-1), phi(i), potential(phi(i)) - potential(phif), vprime(phi(i)), vdprime(phi(i)), dphi(i), d2phi(i), phi_spec(i), transform%wFunc(i-1)
    enddo
    write(u,*)
  end subroutine output_simple

  !>@brief
  !> Initialise our initial guess for the instanton profile
  !>
  !>@param[in] prev  If True - initialise profile using previous numerical profile
  !>                 If False - initialise using an analytic formula
  subroutine initialise_fields(prev)
    logical, intent(in) :: prev

    if (prev) then
       phi(:) = phi_prev(:)
    else
       if (delta > 0.7) then
          phi = -2._dl*atan(exp((transform%xGrid-r0)*meff)) + pi ! For Drummond
       else
          phi = -0.5_dl*(phit-phif)*tanh((transform%xGrid-r0)/w0) + 0.5_dl*(phit+phif)
       endif
!       call thin_wall_profile()
    endif
  end subroutine initialise_fields

  !>@todo
  !>@arg Modify this to allow for more than one field being solved for
  subroutine output(phi,tForm,sol,init)
    real(dl), dimension(1:), intent(in) :: phi
    type(Chebyshev), intent(in) :: tForm
    type(Solver), intent(in) :: sol
    logical, intent(in) :: init
    integer :: u, i, sz
    real(dl), dimension(:), allocatable :: phi_spec

    sz = size(phi)
    allocate(phi_spec(1:sz))
    phi_spec = matmul(tForm%fTrans,phi)
    !!!! Make this safer
    if (init) then
       u = 98
       open(unit=u,file='instanton-full.dat')
    endif

    do i=1,sz
       write(u,*) tForm%xGrid(i-1), phi(i), potential(phi(i)), vdprime(phi(i)), phi_spec(i)
    enddo
  end subroutine output

  !>@brief
  !> Compute the radius and width of the thin-wall using the tension, vacuum splitting and number of dimensions
  subroutine thin_wall_params(rinit,width,sigma,drho,dim)
    real(dl), intent(out) :: rinit, width
    real(dl), intent(in) :: sigma, drho
    integer, intent(in) :: dim

    rinit = dble(dim)*sigma / drho
    width = 1._dl  ! Make this actually work
  end subroutine thin_wall_params

  !>@brief
  !> Compute the thin wall profile.  This will be more easily stored in the model subroutine
  subroutine thin_wall_profile(rvals,phi,r0,w0,phit,phif)
    real(dl), dimension(:), intent(in) :: rvals
    real(dl), dimension(:), intent(out) :: phi
    real(dl), intent(in) :: r0, w0, phit, phif
    
    phi = -0.5_dl*(phit-phif)*tanh((rvals-r0)/w0) + 0.5_dl*(phit+phif)
  end subroutine thin_wall_profile

  !>@brief
  !> Given the previous instanton profile, adjust the coordinate mapping to the new grid
  subroutine adapt_grid()
  end subroutine adapt_grid

  !>@brief
  !> Interpolate our previous instanton solution to the new collocation grid
  subroutine interpolate_to_new_grid()
  end subroutine interpolate_to_new_grid
  
  subroutine get_evectors()
  end subroutine get_evectors

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
  
end program instanton
