program Instanton
  use constants
  use Cheby
  use Nonlinear_Solver
  implicit none

! Nonlinear Solver Parameters and storage.  To be moved to a separate module upon code cleanup
!  type(Field_Model) :: model
  type(Chebyshev) :: transform
  type(Solver) :: solv
  real(dl), dimension(:), allocatable :: phi, phi_prev
  real(dl), dimension(:), allocatable :: model_params

  integer :: order, n
  real(dl) :: w, len  ! parameters controlling mapping of collocation grid
  integer :: i

  ! Move these parameters somewhere else
  real(dl) :: phit, phif,r0,w0

  real(dl) :: delta

  delta = 0.1_dl
  r0 = 1.5_dl*2._dl**0.5/delta; w0=2.**0.5
  phif=-1.; phit=1.
  len = r0*3._dl**0.5
  
  order = 100; n=order+1
  w = 0.5_dl
  ! Initialise our derivatives and set up collocation grid
  call create_chebyshev(transform,order,2,.false.,.true.)
  call transform_to_evens(transform)
  call cluster_points(transform,w,.true.)
  call transform_double_infinite(transform,len)

!  print*,""
!  print*,"X Grid is "
!  print*,transform%xGrid
!  print*,""

  call create_solver(solv,n,100,0.1_dl)
  call initialise_equations(transform,delta)
!  call create_model(model)
 
!  print*,""
!  print*,"L0 is "
!  do i=1,n
!     print*,L0(i,:)
!  enddo
!  print*,""

  allocate(phi(0:order),phi_prev(0:order))
  call initialise_fields(.false.)
  
  call get_vacuum(phif); call get_vacuum(phit)
  print*,"vacua are ", phit, phif
  phi = -0.5_dl*(phit-phif)*tanh((transform%xGrid-r0)/w0) + 0.5_dl*(phit+phif)

!  print*,"phi is "
!  print*,phi
!  print*,""

  call solve(solv, phi)
  call output_simple(transform%xGrid,phi,.true.)

!  do i=1,nparam
!     instanton(:) = instanton_prev(:)
!     call solve(solv, instanton)
!     instanton_prev(:) = instanton(:)
!  enddo

  call delete_solver(solv)
  
contains

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

!
! Define the model through it's potential
!
  elemental function potential(phi)
    real(dl) :: potential
    real(dl), intent(in) :: phi
    potential = 0.25_dl*(phi**2-1._dl)**2 + delta*(phi**3/3._dl - phi)
  end function potential

  elemental function vprime(phi)
    real(dl) :: vprime
    real(dl), intent(in) :: phi
    vprime =  (phi+delta)*(phi**2 - 1.)
  end function vprime

  elemental function vdprime(phi)
    real(dl) :: vdprime
    real(dl), intent(in) :: phi
    vdprime =  3.*phi**2 - 1. + 2.*delta*phi
  end function vdprime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End Delete
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_simple(xvals,phi,init)
    real(dl), dimension(1:), intent(in) :: xvals, phi
    logical, intent(in) :: init
    integer :: i, sz
    integer, parameter :: u = 60

    if (init) open(unit=u,file='instanton.dat')

    sz = size(xvals)
    do i=1,sz
       write(u,*) xvals(i), phi(i)
    enddo
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
!       call thin_wall_profile()
    endif
  end subroutine initialise_fields

  !>@todo
  !>@arg Modify this to allow for more than one field being solved for
  subroutine output(phi,tForm,sol,init)
    real(dl), dimension(:), intent(in) :: phi
    type(Chebyshev), intent(in) :: tForm
    type(Solver), intent(in) :: sol
    logical, intent(in) :: init
    integer :: u, i, sz

    sz = size(phi)
    !!!! Make this safer
    if (init) then
       u = 98
       open(unit=u,file='instanton.dat')
    endif

    do i=0,sz-1
       write(u,*) tForm%xGrid(i), phi(i)
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

  subroutine get_evalues()
  end subroutine get_evalues
  
  subroutine get_evectors()
  end subroutine get_evectors
  
  subroutine get_action()
  end subroutine get_action
  
end program instanton
