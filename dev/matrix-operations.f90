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

end module Matrix
