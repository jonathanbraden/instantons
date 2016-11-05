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
module Nonlinear_Solver
  use constants
  use Model
  implicit none
  
  type Solver
     integer :: nVar
     real, dimension(:,:), allocatable :: L  ! Linear operator
     real, dimension(:), allocatable :: S    ! Source term
     real, dimension(:), allocatable :: res, del  ! Residual and perturbation
     ! Whatever else I need
  end type Solver
  
contains

  subroutine intialise_solver(this,n)
    type(Solver), intent(inout) :: this
    this%nVar = n
    allocate(this%L(1:n,1:n)); allocate(this%S(1:n))
    allocate(this%res(1:n), this%del(1:n))
  end subroutine intialise_solver

  subroutine delete_solver(this)
    type(Solver), intent(inout) :: this
    this%nVar = -1
    deallocate(this%L, this%S)
    deallocate(this%res, this%del)
  end subroutine delete_solver

  subroutine solve(this)
    type(Solver), intent(in) :: this
    integer :: i
    integer, parameter :: imax = 100
    
    do i=1,imax
       call line_iteration(this)
    enddo
  end subroutine solve
  
  !>@brief
  !> Nonlinear solver based on Newton's method, with the extension to consider variable
  !> distances along the Newton iteration path to ensure the solution is converging
  subroutine line_iteration(this)
    type(Solver), intent(in) :: this
    integer :: i
    
    ! Use the equations of motion to fill the Linear and Source term
    call fill_equations(this%L,this%S) ! write this.  Currently in model file
    
    this%res = sum(this%S*this%S)
    ! Fix this one up.  Where does ipiv come from
    call DGESV(this%nVar,1,this%L,ipiv,this%S,this%nVar,info)
    this%del = S
    if (info /= 0) then
       print*,"Error inverting linear matrix in solver"
       stop  ! Improve this error handling
    endif

    b1 = b0 + 0.1_dl ! Why the hell do I need this?
    alpha = 2._dl
    do i=1,8
       alpha = alpha/2._dl
       tmp = freal + alpha*del
       S = -matmul(l0,tmp) + vprime(tmp)
       b1 = sum(S*S)
       if (b1 < b0) exit
    enddo

    freal = freal + alpha*del

    this%res = source()
  end subroutine line_iteration

  !>@brief
  !> Look for solutions using combined Newton and gradient descent method
  subroutine newton_gradient(this)
    type(Solver), intent(in) :: this
    integer ::

  end subroutine newton_gradient
  
end module Nonlinear_Solver
