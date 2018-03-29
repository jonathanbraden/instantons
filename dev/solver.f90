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

! Currently uncommenting this breaks the solver
!#define DEBUG_SOLVER 1

module Nonlinear_Solver
  use constants
  use Cheby
  use Model
  implicit none

  real(dl), external :: DLANGE
  
  type Solver_Storage
     integer, dimension(:), allocatable :: iwork
     real(dl), dimension(:), allocatable :: rwork
     real(dl), dimension(:,:), allocatable :: lu_factor
     real(dl), dimension(:), allocatable :: row_scale, col_scale
  end type Solver_Storage
  
  type Solver
     integer :: nVar, u
     real(dl) :: alpha
     integer :: maxIter; real(dl) :: kick_param  ! Adjustable solver parameters
     real(dl), dimension(:,:), allocatable :: L, L_const  ! Linear operator
     real(dl), dimension(:), allocatable :: S    ! Source term
     real(dl), dimension(:), allocatable :: del  ! Residual and perturbation
     ! Make these ones private
     integer, dimension(:), allocatable :: ipiv  ! Needed for DGESV in the solver
     real(dl), dimension(:), allocatable :: f_prev, S_prev
     logical :: spec_space  ! whether or not we work in spectral space (not yet used)
#ifdef DEBUG_SOLVER
     type(Solver_Storage) :: mat_store
#endif
  end type Solver

!  abstract interface
!     subroutine src(fld,eom)
!       real(dl), dimension(:), intent(in) :: fld
!       real(dl), dimension(:), intent(out) :: eom
!     end subroutine src
!  end interface

contains

  !>@brief
  !> Set aside storage space and necessary variables for our nonlinear iterations used
  !> to solve a boundary value problem
  !>
  !>@param[in,out] this
  !>@param[in] n  Number of variables to solve for
  !>@param[in] nit Maximum number of iterations for the solver before quitting
  !>@param[in] kick Parameter controlling the ability to leave a local minimum (currently not functional)
  subroutine create_solver(this,n, nit, kick)
    type(Solver), intent(inout) :: this
    integer, intent(in) :: n,nit
    real(dl), intent(in) :: kick
    this%nVar = n; this%maxIter = nit; this%kick_param = kick
    allocate(this%L(1:n,1:n),this%L_const(1:n,1:n)); allocate(this%S(1:n))
    allocate(this%del(1:n)); allocate(this%f_prev(1:n))
    allocate(this%ipiv(1:n))
    allocate(this%S_prev(1:n))
#ifdef DEBUG_SOLVER
    call create_solver_storage(this%mat_store,n)
#endif
    this%u = 99  ! Change this to a call to the next open solver
    open(unit=this%u,file='solver-output.dat')
  end subroutine create_solver

  subroutine create_solver_storage(this, n)
    type(Solver_Storage), intent(out) :: this
    integer, intent(in) :: n

    allocate(this%iwork(1:n)); allocate(this%rwork(1:4*n))
    allocate(this%lu_factor(1:n,1:n)); allocate(this%row_scale(1:n),this%col_scale(1:n))
  end subroutine create_solver_storage

  subroutine delete_solver(this)
    type(Solver), intent(inout) :: this
    this%nVar = -1
    deallocate(this%L, this%L_const, this%S)
    deallocate(this%del, this%f_prev)
    deallocate(this%ipiv)
  end subroutine delete_solver

  subroutine solve(this,f_cur,test)
    type(Solver), intent(inout) :: this
    real(dl), dimension(:), intent(inout) :: f_cur
    logical, intent(out), optional :: test
    integer :: i; logical :: convg

    convg = .false.
    do i=1,this%maxIter
       call line_iteration(this,f_cur)
!       call output_solver(this)
!       call print_solver(this)
       if (stop_solver(this)) then; print*,"Converged in ",i," steps"; convg=.true.; exit; endif
    enddo
    if (.not.convg) print*,"Failed to converge "
    call print_solver(this)   
    write(this%u,*) ""

    if (present(test)) test = convg
  end subroutine solve
  
  !>@brief
  !> Stopping condition for our nonlinear solver
  function stop_solver(this) result(test)
    type(Solver), intent(in) :: this
    logical :: test
    real(dl), parameter :: eps = 1.e-10
    
    test = (maxval(abs(this%S(:))) < eps) .and. (maxval(abs(this%del(:))) < eps)
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

  subroutine solve_radius(this,f_cur)
    type(Solver), intent(inout) :: this
    real(dl), dimension(:), intent(inout) :: f_cur
    real(dl), dimension(:), allocatable :: f_prime

    allocate(f_prime(1:size(f_cur)))
    f_prime = matmul(transform%derivs(:,:,1),f_cur)
    
  end subroutine solve_radius

  !>@brief
  !> Decompose the proposed perturbation to the field into eigenmodes of the linear operator.
  !> In particular, exctract the dependence of the current derivative of the solution,
  !> corresponding to the collective mode associated with the overall bubble radius
  subroutine decompose_perturbation(this)
    type(Solver), intent(in) :: this

    ! 1. Compute derivative of the current solution
    ! 2. Get the eigenmodes of the linear operator
    ! 3. Decompose the source into the eigenmodes (find the weight for orthogonality)
    ! 4. Figure out what these modes correspond to
  end subroutine decompose_perturbation
  
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
#ifdef DEBUG_SOLVER
    real(dl) :: cond_num, mat_norm(1:2)
    character :: norm
#endif

    n = this%nVar
    
    call variation(f_cur, this%L)
    call source(f_cur, this%S)
    
    res = sqrt(sum(this%S**2))  ! Current residual
    this%S_prev = this%S
    ! Compute the required perturbation
    call DGESV(n,1,this%L,n,this%ipiv,this%S,n,info)
#ifdef DEBUG_SOLVER
    !mat_norm = F06RAF('1',n,n,this%L,n, ) ! NAG version of DLANGE
    mat_norm(1) = DLANGE('1',n,n,this%L,n,this%mat_store%rwork(1:n))
    mat_norm(2) = DLANGE('I',n,n,this%L,n,this%mat_store%rwork(1:n))
    call DGECON('1',n,this%L,n,mat_norm,cond_num,this%mat_store%rwork(1:4*n),this%mat_store%iwork(1:n),info)
    print*,"One-norm condition number of matrix is ", cond_num
    call DGECON('I',n,this%L,n,mat_norm,cond_num,this%mat_store%rwork(1:4*n),this%mat_store%iwork(1:n),info)
    print*,"Infinity-norm condition number of matrix is ", cond_num
    call DGETRF(n,n,this%L,n,this%ipiv,info) ! Get LU factorisation
    call DGETRS('N',n,1,this%L,n,this%ipiv,this%S,n,info)
    ! Do eigenvalue decomposition

!    call DGESVX(,'N',n,1,this%L,n,this%lu_factor,n,this%ipiv,'N', , ,this%S,n, ,
#endif

    this%del = this%S
    if (info /= 0) then
       print*,"Error inverting linear matrix in solver"
       stop  ! Improve this error handling
    endif

    ! Why am I ever setting this thing here?
    b1 = res + this%kick_param ! If we're not converging, this allows us to kick ourselves
    alpha = 2._dl
    do i=1,8
       alpha = alpha/2._dl
       this%f_prev = f_cur + alpha*this%del
       call source(this%f_prev, this%S) ! replace with a call to source
       err_max = maxval(abs(this%S))
       err_rms = sqrt(sum(this%S**2))
#ifdef DEBUG_SOLVER
       print*,"RMS error is ",err_rms**2," Residual is ",res**2," alpha is ",alpha
#endif
       if (err_rms < res) exit
    enddo
#ifdef DEBUG_SOLVER
    print*,""
#endif
    this%f_prev = f_cur            ! Store previous iteration
    f_cur = f_cur + alpha*this%del  ! Update function
    call source(f_cur, this%S)
    this%alpha = alpha          ! Store violation of the equation of motion
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
    real(dl), dimension(:), allocatable :: f_deriv

    allocate(f_deriv(1:size(this%f_prev)))
    f_deriv = matmul(transform%derivs(:,:,1),this%f_prev)
    do i=1,this%nVar
      write(this%u,*) this%f_prev(i), this%del(i), this%S(i), this%S_prev(i), f_deriv(i)
    enddo
    write(this%u,*)
  end subroutine output_solver

  !>@brief
  !> Write a brief summary of the last state of the nonlinear solver to the screen
  subroutine print_solver(this)
    type(Solver), intent(in) :: this

    print*,"RMS violation of the EOM is ",sqrt(sum(this%S**2))
    print*,"Alpha on previous step is ",this%alpha
  end subroutine print_solver

end module Nonlinear_Solver
