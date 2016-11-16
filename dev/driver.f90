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
  
  order = 3; n=order+1
  w = 1._dl; len = 20._dl
  ! Initialise our derivatives and set up collocation grid
  call create_chebyshev(transform,order,2,.false.,.true.)
  call transform_to_evens(transform)
  call cluster_points(transform,w,.true.)
  call transform_double_infinite(transform,len)

  call create_solver(solv,n,5,0.1_dl)
!  call create_model(model)
 
  allocate(phi(0:order),phi_prev(0:order))
  call initialise_fields(.false.)
  call solve(solv, phi)
  
!  do i=1,nparam
!     instanton(:) = instanton_prev(:)
!     call solve(solv, instanton)
!     instanton_prev(:) = instanton(:)
!  enddo

  call delete_solver(solv)
  
contains

  !>@brief
  !> Initialise our initial guess for the instanton profile
  !>
  !>@param[in] prev  If True - initialise profile using previous numerical profile
  !>                 If False - initialise using an analytic formula
  subroutine initialise_fields(prev)
    logical, intent(in) :: prev

    if (prev) then
       
    else
!       call thin_wall_profile()
    endif
  end subroutine initialise_fields

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
