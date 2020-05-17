module Fluctuations
  use constants, only : dl
  use Utils, only : newunit
  use Cheby
  use Instanton_Class
  use equations, only : variation

!  real(dl), dimension(:), allocatable :: work
!  real(dl), dimension(:), allocatable :: evalreal, evalimag
!  integer :: asize, ierror, iwork
!  real(dl) :: dummy(1:1)
  
contains

  subroutine init_evals()
  end subroutine init_evals
  
  !>@brief
  !> Compute the eigenvalues for spherical harmonic l around the instanton profile
  subroutine get_eigenvalues_neumann(this,ev_r,ev_i,l,d)
    type(Instanton), intent(in) :: this
    real(dl), dimension(:), intent(out) :: ev_r,ev_i
    integer, intent(in) :: l,d

    integer :: npt, ierror
    real(dl) :: norm
    real(dl), dimension(0:this%ord,0:this%ord) :: op
    real(dl), dimension(:), allocatable :: work
    integer :: iwork
    real(dl), dimension(1) :: dummy
    
    norm = l*(l+d-1._dl)  ! Check this normalisation
    
    ! Step 1: Transform interpolation grid to lobatto if necessary
    ! Step 2: Interpolate phi onto lobatto grid
    ! Step 3: Call variation matrix (on appropriate grid)
    ! Step 4: Add angular momentum barrier to variation matrix (-l(l+d-2)/r^2)
    ! Step 5: Fix boundar conditions on variation matrix (if Lobatto grid)
    ! Step 6: Solve eigenvalue problem
    !   - allocate solver storage for LAPACK
    !   - call eigenvalue solver

    
    ! Assign workspace for solver (move this to a separate place so I don't have to keep redoing it
    npt = this%ord + 1
    call DGEEV('N','N',npt,op,npt,ev_r,ev_i,dummy,1,dummy,1,dummy,-1,ierror)
    if (ierror == 0) then
       iwork = int(dummy(1))
       allocate(work(1:iwork))
    endif

    ! Make sure I set b.c. in here properly
    ! 5. Set the boundary conditions appropriately here
    call variation(this%phi,op)
    op = -op
    do i=0,this%ord
       op(i,i) = op(i,i) + norm/this%tForm%xGrid(i)**2
    enddo

    ! Step 6: obtain eigenvalues
    call DGEEV('N','N',npt,op,npt,ev_r,ev_i,dummy,1,dummy,1,work,iwork,ierror)
  end subroutine get_eigenvalues_neumann

  subroutine get_eigenvalues_dirichlet(this,ev,l,d)
    type(Instanton), intent(in) :: this
    real(dl), dimension(:), intent(out) :: ev
    integer, intent(in) :: l, d
  end subroutine get_eigenvalues_dirichlet
  
  subroutine get_eigenvectors(this,ev_r,ev_i,evec)
    type(Instanton), intent(in) :: this
    real(dl), dimension(:), intent(out) :: ev_r,ev_i
    real(dl), dimension(:,:), intent(out) :: evec

    ! Allocation
!    call DGEEV('N','V',npt,op,npt,ev_r,ev_i,dummy_,npt,v_r,npt,dummy_,-1,ierror)
    ! Compuationa
!    call DGEEV('N','V',npt,op,npt,ev_r,ev_i,dummy_,npt,v_r,npt,work,iwork,ierror)


  end subroutine get_eigenvectors

  !>@brief
  !> Integrate Ricatti equation to get partial wave determinant
  subroutine determinant_l(this,l,d)
    type(Instanton), intent(in) :: this
    integer, intent(in) :: l,d

  end subroutine determinant_l

  subroutine compute_determinant(this,lmax,d)
    type(Instanton), intent(in) :: this
    integer, intent(in) :: lmax, d

    integer :: i

    ! Compute l=0 determinant here
    ! Compute l=1 determinant here, make sure it's zero (compute zero removed determinant
    do i=2,lmax
    enddo
  end subroutine compute_determinant
  
end module Fluctuations
