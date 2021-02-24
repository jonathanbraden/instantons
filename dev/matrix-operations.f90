!#define OUTPUT T

module Matrix
  use constants
  implicit none
  
contains

  !>@brief
  !> Invert an NxN matrix
  !>
  !>@param[in]   A: real(dl) - nxn matrix to invert
  !>@param[out]  Ainv: real(dl) - nxn array to store the inverse of A
  !>@param[in]   n : size of arrays
  subroutine invert_matrix(A, Ainv, n)
    integer, intent(in) :: n
    real(dl), dimension(n,n), intent(in) :: A
    real(dl), dimension(n,n), intent(out) :: Ainv

    ! Workspace for LAPACK
    real(dl), allocatable, dimension(:) :: work
    integer :: lwork, lda, info
    integer, allocatable, dimension(:) :: ipiv
    integer :: deallocatestatus

  ! check how much space LAPACK needs here
    lda = n  ! should really assign separately as input
    lwork = n*n
    allocate (work(lwork)); allocate (ipiv(n))

! Perform LU decomposition
    Ainv = A
    call DGETRF( n, n, Ainv, lda, ipiv, info )

! Print an error if the LU decomposition failed
    if (info.eq.0) then
#ifdef OUTPUT
       print*, "LU decomposition successful"
#endif
    elseif (info < 0) then
       print*, "LU decomposition: illegal value at ", info
       stop
    elseif (info > 0) then
       print*, "Singular Matrix U = 0 at ",info
    endif

! USE LU decomposition to compute inverse
! To do: finish workspace query to determine size of work array
!  call DGETRI(N, Ainv, LDA, IPIV, WORK, -1, INFO)
!  lwork = work(1)
!  deallocate(work, 
!  allocate(work(lwork))
    call DGETRI( n, Ainv, lda, IPIV, work, lwork, info )

    if (info.ne.0) then
       stop "Matrix inversion failed!"
    else
#ifdef OUTPUT
       print*, "Inverse Successful"
#endif
    endif

    ! clean up temporary workspace
    deallocate(ipiv, STAT=deallocatestatus)
    deallocate(work, STAT=deallocatestatus)

  end subroutine invert_matrix

#ifdef WRAPPERS
  !>@brief
  !> Wrapper for DGEEV that automatically assigns workspace
  subroutine compute_eigenvalues(A,ev_r,ev_i)
    real(dl), dimension(:,:), intent(in) :: A
    real(dl), dimension(:), intent(out) :: ev_r,ev_i

    real(dl), dimension(:), allocatable :: work
    integer :: iwork, npt
    real(dl), dimension(1) :: dummy

    npt = size(ev_r)  ! fix this
    call DGEEV('N','N', ,A, ,ev_r,ev_i,dummy,1,dummy,1,dummy,-1,ierror)
    if (ierror == 0) then
       iwork = int(dummy(1))
       allocate(work(1:iwork))
    endif
    call DGEEV('N','N',npt,A,npt,ev_r,ev_i,dummy,1,dummy,1,work,iwork,ierror)
    deallocate(work)
  end subroutine compute_eigenvalues

  ! Combine this with the above by making evec and optional 
  subroutine compute_eigenvalues_eigenvectors(A,ev_r,ev_i,evec)
    real(dl), dimension(:,:), intent(in) :: A
    real(dl), dimension(:), intent(out) :: ev_r,ev_i
    real(dl), dimension(:,:), intent(out), optional :: evec
    
    integer :: npt, ierror
    real(dl), dimension(:), allocatable :: work
    integer :: iwork
    real(dl), dimension(1) :: dummy
    logical :: get_vec

    get_vec = .false.; if (present(evec)) get_vec = .true.

    npt = size(ev_r)

    if (get_vec) then
       call DGEEV('N','V',npt,A,npt,ev_r,ev_i,dummy,1,evec,npt,dummy,-1,ierror)
    else
       call DGEEV('N','N',npt,A,npt,ev_r,ev_i,dummy,1,dummy,1,dummy,-1,ierror)
    endif
    
    if (ierror == 0) then
       iwork = int(dummy(1))
       allocate(work(1:iwork))
    endif

    if (get_vec) then
       call DGEEV('N','V',npt,op,npt,ev_r,ev_i,dummy,1,evec,npt,work,iwork,ierror)
    else
       call DGEEV('N','V',npt,op,npt,ev_r,ev_i,dummy,1,evec,npt,work,iwork,ierror)
    endif
    deallocate(work)
  end subroutine compute_eigenvalues_eigenvectors
#endif
  
end module Matrix
