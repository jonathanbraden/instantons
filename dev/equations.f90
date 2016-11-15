!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> University College London
!>
!>@brief
!> Store generic equations for solving the instanton equation for both canonical
!> scalar fields and non-linear sigma model scalar fields
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Equations
  use constants, only: dl
  use Cheby  ! Do I actually need this
  implicit none

  real(dl), dimension(:,:), allocatable :: L0

contains
  
  subroutine initialise_canonical(nlat,nfld)
    integer, intent(in) :: nlat, nfld
    integer :: nvar
    nvar = nlat*nfld

    allocate(L0(1:nvar,1:nvar))
  end subroutine initialise_canonical

  !>@brief
  !> Compute the result of the equations of motion (which should be 0) on the proposed
  !> field solution.  This acts as a source term for a Newton method iteration
  subroutine source_canonical(S,f_cur,nlat,nfld)
    real(dl), dimension(:), intent(out) :: S
    real(dl), dimension(1:n,1:nfld), intent(in) :: f_cur ! Nope, wrong dimensions
    integer, intent(in) :: nlat, nfld
    integer :: i,l; real(dl) :: tmp

    do l=1,nfld
       S((l-1)*n+1:l*n) = matmul( L0,f_cur(:,l) )
       do i=1,n
          S(i+(l-1)*n) = S(i+(l-1)*n) + vprime(f_cur,l)
       enddo
    enddo
  end subroutine source_canonical

  !>@brief
  !> Compute the first variation of the equations of motion (second variation of the action)
  !> in order to do a Newton method iteration on the nonlinear solver
  subroutine variation_canonical(L,f_cur,nlat,nfld)
    real(dl), dimension(:,:), intent(out) :: L
    real(dl), dimension(1:n), intent(in) :: f_cur
    integer, intent(in) ::  n,nfld
    integer :: i,l; real(dl) :: tmp

    L = 0._dl  ! Not the most efficient, but the easiest to understand
    do l=1,nfld
       L(l*n+1:(l+1)*n,l*n+1:(l+1)*n) = L0
    enddo
    do i=1,n
       tmp = vdprime(f_cur(i)) ! Nope, this doesn't work
       do l=1,nfld
          L(i+l,i+l) = L(i+l,i+l) + tmp  ! Turn this into a simple array slicing
       enddo
    enddo
    ! Fill in the off-diagonal blocks here
  end subroutine variation_canonical

end module Equations
