!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> University College London
!>
!>@brief
!> Contains subroutines for solution of the non-linear boundary value problem
!> required to obtain instanton profiles
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define SOLVER_DEBUG 1

module Nonlinear_Solver
  use constants
  use Cheby
  use Model
  implicit none
  
  type Solver
     integer :: nVar, u
     integer :: maxIter; real(dl) :: kick_param  ! Adjustable solver parameters
     real(dl), dimension(:,:), allocatable :: L, L_const  ! Linear operator
     real(dl), dimension(:), allocatable :: S    ! Source term
     real(dl), dimension(:), allocatable :: del  ! Residual and perturbation
     ! Make these ones private
     integer, dimension(:), allocatable :: ipiv  ! Needed for DGESV in the solver
     real(dl), dimension(:), allocatable :: f_prev
     logical :: spec_space  ! whether or not we work in spectral space (not yet used)
  end type Solver

!  abstract interface
!     subroutine src(fld,eom)
!       real(dl), dimension(:), intent(in) :: fld
!       real(dl), dimension(:), intent(out) :: eom
!     end subroutine src
!  end interface

contains

  subroutine create_solver(this,n, nit, kick)
    type(Solver), intent(inout) :: this
    integer, intent(in) :: n,nit
    real(dl), intent(in) :: kick
    this%nVar = n; this%maxIter = nit; this%kick_param = kick
    allocate(this%L(1:n,1:n),this%L_const(1:n,1:n)); allocate(this%S(1:n))
    allocate(this%del(1:n)); allocate(this%f_prev(1:n))
    allocate(this%ipiv(1:n))
    
    this%u = 99  ! Change this to a call to the next open solver
    open(unit=this%u,file='solver-output.dat')
  end subroutine create_solver

  subroutine delete_solver(this)
    type(Solver), intent(inout) :: this
    this%nVar = -1
    deallocate(this%L, this%L_const, this%S)
    deallocate(this%del, this%f_prev)
    deallocate(this%ipiv)
  end subroutine delete_solver

  subroutine solve(this,f_cur)
    type(Solver), intent(inout) :: this
    real(dl), dimension(:), intent(inout) :: f_cur
    integer :: i; logical :: test
    
    do i=1,this%maxIter
       print*,"Iterating solver ",i
       call line_iteration(this,f_cur)
       call output_solver(this)
       if (stop_solver(this)) exit
    enddo
  end subroutine solve
  
  !>@brief
  !> Stopping condition for our nonlinear solver
  function stop_solver(this) result(test)
    type(Solver), intent(in) :: this
    logical :: test

    test = (maxval(abs(this%S(:))) < 1.e-10) .and. (maxval(abs(this%del(:))) < 1.e-10)
  end function stop_solver

  !>@brief
  !> Nonlinear solver based on Newton's method, with the extension to consider
  !> variable length paths along the Newton iteration to improve convergence properties.
  !> This version requires as input a subroutine to fill the linear and source terms
!  subroutine line_iteration(this, hessian, source)
!    type(Solver), intent(inout) :: this
!    procedure (lin), nopass :: hessian
!    procedure (src), nopass :: source

!    call eqn(this%L,this%S,phi_cur)
!  end subroutine line_iteration

  !>@brief
  !> Nonlinear solver based on Newton's method, with the extension to consider variable
  !> distances along the Newton iteration path to ensure the solution is converging
  !>
  !> TODO: figure out what to do with the b1 and b0
  !>       Add in the source, and hessian calls somewhere
  subroutine line_iteration(this,f_cur)
    type(Solver), intent(inout) :: this
    real(dl), dimension(:), intent(inout) :: f_cur
    integer :: i, info, n
    real(dl) :: alpha
    real(dl) :: err_rms, err_max, res
    real(dl) :: b1  ! remove this later
    
    n = this%nVar
    
    call variation(f_cur, this%L)
    call source(f_cur, this%S)
    
    res = sqrt(sum(this%S**2))  ! Current residual
    print*,"Residual is ",res
    ! Compute the required perturbation
    call DGESV(n,1,this%L,n,this%ipiv,this%S,n,info)
    this%del = this%S
    if (info /= 0) then
       print*,"Error inverting linear matrix in solver"
       stop  ! Improve this error handling
    endif

!    print*,"del vector is "; print*,this%del

    ! Why am I ever setting this thing here?
    b1 = res + this%kick_param ! If we're not converging, this allows us to kick ourselves
    alpha = 2._dl
    do i=1,8
       alpha = alpha/2._dl
       this%f_prev = f_cur + alpha*this%del
       call source(this%f_prev, this%S) ! replace with a call to source
       err_max = maxval(abs(this%S))
       err_rms = sqrt(sum(this%S**2))
       print*,"b1 is ",err_rms**2," b0 is ",res**2," alpha is ",alpha
       if (err_rms < res) exit
    enddo

    print*,"Done line search with alpha = ",alpha
    print*,""
    this%f_prev = f_cur            ! Store previous iteration
    f_cur = f_cur + alpha*this%del  ! Update function
    call source(f_cur, this%S)          ! Store violation of the equation of motion
  end subroutine line_iteration

  !>@brief
  !> Look for solutions using combined Newton and gradient descent method
  !> This will presumably be more robust, although somewhat slower.
  subroutine levenberg_marquadt(this)
    type(Solver), intent(in) :: this
  end subroutine levenberg_marquadt

  !>@brief
  !> Look for solutions to the differential equation by performing a gradient flow.
  !>
  !> Solutions for the differential equation are found by minimising the scalar function defined by
  !>  \f[
  !>    E = \sum_i\left(\frac{\partial S}{\partial\phi_i}\right)^2 = \sum_i S_i
  !>  \f]
  !> where \f$S\f$ is the action from which the equations of motion are derived.
  !> This method follows the gradient of \f$E\f$ with respect to variations in the field values.
  !> The corresponding gradient vector is given by
  !>  \f[
  !>    \delta\phi_l \propto \sum_i S_i\frac{\delta S_i}{\delta \phi_l} = \sum_i S_iL_{ij}
  !>  \f]
  !> The solution is then updated to
  !>  \f[
  !>    \phi^{(i+1)} = \phi^{(i)} + \alpha\delta\phi^{(i)}
  !>  \f]
  !> where the parameter \f$\alpha\f$ is determined by doing a line-search along the gradient vector direction
  !> such that it minimises \f$E(\phi^{(i+1)})\f$
  !>@todo
  !>@arg Fill in some more explanation of what this is doing
  subroutine gradient_solve(this, f_cur)
    type(Solver), intent(in) :: this
    real(dl), dimension(:), intent(inout) :: f_cur

    real(dl) :: alpha
    integer :: i
  end subroutine gradient_solve

  !>@todo
  !>@arg Implement this thing
  subroutine conjugate_gradient(this, f_cur)
    type(Solver), intent(in) :: this
    real(dl), dimension(:), intent(inout) :: f_cur
  end subroutine conjugate_gradient

  !>@brief
  !> Output information from the nonlinear solver.  This can be useful for debugging, or to assess
  !> convergence of our iterator
  subroutine output_solver(this)
    type(Solver), intent(in) :: this
    integer :: i

    do i=1,this%nVar
       write(this%u,*) this%f_prev(i), this%del(i), this%S(i)
    enddo
    write(this%u,*)
  end subroutine output_solver

end module Nonlinear_Solver
