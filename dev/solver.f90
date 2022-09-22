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
  use Utils, only : newunit
  use Cheby
  use Equations
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
     integer :: nVar, u, fNum_summary
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
     logical :: exists = .false.
     logical :: output_result, output_iterations
     real(dl) :: tol_dphi, tol_eom, tol_action

     logical :: use_pre = .false.
     real(dl), dimension(:,:), allocatable :: precond
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
    if (this%exists) call delete_solver(this)
    this%nVar = n; this%maxIter = nit; this%kick_param = kick
    allocate(this%L(1:n,1:n),this%L_const(1:n,1:n)); allocate(this%S(1:n))
    allocate(this%del(1:n)); allocate(this%f_prev(1:n))
    allocate(this%ipiv(1:n))
    allocate(this%S_prev(1:n))
#ifdef DEBUG_SOLVER
    call create_solver_storage(this%mat_store,n)  ! Need to delete this as needed
#endif
    this%exists = .true.

    ! Need to pass in options for the rest of these
    open(unit=newunit(this%u),file='solver-output.dat')
    open(unit=newunit(this%fNum_summary),file='solver-summary.dat')
    this%output_result = .true.
    this%output_iterations = .true.
    this%tol_dphi = 1.e-10; this%tol_eom = 1.e-10; this%tol_action = 1.e-10

    if (this%use_pre) allocate(this%precond(1:n,1:n))
    
  end subroutine create_solver

  subroutine create_solver_storage(this, n)
    type(Solver_Storage), intent(out) :: this
    integer, intent(in) :: n

    allocate(this%iwork(1:n)); allocate(this%rwork(1:4*n))
    allocate(this%lu_factor(1:n,1:n)); allocate(this%row_scale(1:n),this%col_scale(1:n))
  end subroutine create_solver_storage

  subroutine delete_solver(this)
    type(Solver), intent(inout) :: this
    logical :: o
    
    this%nVar = -1
    deallocate(this%L, this%L_const, this%S)
    deallocate(this%del, this%f_prev)
    deallocate(this%ipiv)
    deallocate(this%S_prev)
    inquire(opened=o,unit=this%u); if (o) close(this%u) 
    inquire(opened=o,unit=this%fNum_summary); if (o) close(this%fNum_summary) 
    this%exists = .false.
  end subroutine delete_solver

  !>@brief
  !> Set the tolerances for the nonlinear solver.
  !
  !>@param[in] (optional) tol_global - A global value for all tolerances
  !>@param[in] (optional) tol_dphi - Tolerance on magnitude of change in solution
  !>@param[in] (optional) tol_eom  - Tolerance of EOM violation
  !>@param[in] (optional) tol_action - Tolerance on change in action
  subroutine set_tolerances(this,tol_global,tol_dphi,tol_eom,tol_action)
    type(Solver), intent(inout) :: this
    real(dl), intent(in), optional :: tol_global, tol_dphi, tol_eom, tol_action

    if (present(tol_global)) then
       this%tol_dphi = tol_global; this%tol_eom = tol_global; this%tol_action = tol_global
    endif
    if (present(tol_dphi)) this%tol_dphi = tol_dphi
    if (present(tol_eom)) this%tol_eom = tol_eom
    if (present(tol_action)) this%tol_action = tol_action
  end subroutine set_tolerances
  
  !>@brief
  !> Output information from the nonlinear solver.  This can be useful for debugging, or to assess convergence of our iterations
  subroutine output_solver(this)!,tPair)
    type(Solver), intent(in) :: this
!    type(Chebyshev), intent(in) :: tForm
    integer :: i
!    real(dl), dimension(1:size(this%f_prev)) :: f_deriv

!    f_deriv = matmul(tForm%derivs(:,:,1),this%f_prev)
    do i=1,this%nVar
      write(this%u,*) this%f_prev(i), this%del(i), this%S(i), this%S_prev(i)!, f_deriv(i)
    enddo
    write(this%u,*)
  end subroutine output_solver

  subroutine output_summary(this)
    type(Solver), intent(in) :: this
    write(this%fNum_summary,*) this%alpha, sqrt(sum(this%S**2)), maxval(abs(this%S)), this%S(1), sqrt(sum(this%del**2))
  end subroutine output_summary
  
  !>@brief
  !> Write a brief summary of the last state of the nonlinear solver to the screen
  subroutine print_solver(this)
    type(Solver), intent(in) :: this
    print*,"RMS violation of the EOM is ",sqrt(sum(this%S**2))
    print*,"Maximal violation of the EOM is ",maxval(abs(this%S))
    print*,"RMS change in field is ",sqrt(sum(this%del**2))
    print*,"Maximum change in field is ",maxval(abs(this%del))
    print*,"Alpha on previous step is ",this%alpha
  end subroutine print_solver
  
  subroutine solve(this,f_cur,test)
    type(Solver), intent(inout) :: this
    real(dl), dimension(:), intent(inout) :: f_cur
    logical, intent(out), optional :: test
    integer :: i; logical :: convg

    convg = .false.
    do i=1,this%maxIter
       call line_iteration(this,f_cur)
       if (this%output_iterations) call output_solver(this)
       if (this%output_iterations) call output_summary(this)
       if (stop_solver(this)) then
          if (this%output_result) print*,"Converged in ",i," steps"
          convg=.true.
          exit
       endif
    enddo
    if (.not.convg) then
       print*,"Failed to converge "
       call print_solver(this)
    endif
    write(this%u,*) ""

    if (present(test)) test = convg
  end subroutine solve
  
  !>@brief
  !> Stopping condition for our nonlinear solver
  function stop_solver(this) result(test)
    type(Solver), intent(in) :: this
    logical :: test
    
    test = (maxval(abs(this%S(:))) < this%tol_eom) .and. (maxval(abs(this%del(:))) < this%tol_dphi)
    !test = (maxval(abs(this%S(:))) < this%tol_eom)
  end function stop_solver

  !>@brief
  !> Nonlinear solver based on a variable step length Newton's method
  !>
  !> TODO: Add in the source, and hessian calls somewhere
  !>       Figure out what to do with b1
  subroutine line_iteration(this,f_cur)
    type(Solver), intent(inout) :: this
    real(dl), dimension(:), intent(inout) :: f_cur
    integer :: i, info, n
    real(dl) :: alpha
    real(dl) :: err_rms, err_max, res
#ifdef DEBUG_SOLVER
    real(dl) :: cond_num, mat_norm(1:2)
    character :: norm
#endif

    n = this%nVar
    call variation(f_cur, this%L)  ! pass this in
    call source(f_cur, this%S)     ! pass this in
    
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
       print*,"Error inverting linear matrix in solver, terminating"
       stop 
    endif

    alpha = 2._dl
    do i=1,8
       alpha = alpha/2._dl
       this%f_prev = f_cur + alpha*this%del
       call source(this%f_prev, this%S) 
       err_max = maxval(abs(this%S))
       err_rms = sqrt(sum(this%S**2))
#ifdef DEBUG_SOLVER
       print*,"RMS error is ",err_rms**2," Residual is ",res**2," alpha is ",alpha
#endif
       if (err_rms < res) exit  ! Maybe this isn't the best condition
    enddo
#ifdef DEBUG_SOLVER
    print*,""
#endif

    this%f_prev = f_cur             ! Store previous iteration
    f_cur = f_cur + alpha*this%del  ! Update function
    call source(f_cur, this%S)      ! This seems redundant with a call above
    this%alpha = alpha              ! Store violation of the equation of motion
  end subroutine line_iteration
  
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

end module Nonlinear_Solver
