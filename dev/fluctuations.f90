module Fluctuations
  use constants, only : dl
  use Utils, only : newunit
  use Cheby
  use Instanton_Class
  use equations, only : variation, set_bc

!  real(dl), dimension(:), allocatable :: work
!  real(dl), dimension(:), allocatable :: evalreal, evalimag
!  integer :: asize, ierror, iwork
!  real(dl) :: dummy(1:1)

  ! Create an appropriate type here
  
contains


  !>@brief
  !> Compute the eigenvalues for spherical harmonic l around the instanton profile
  subroutine get_eigenvalues(this,ev_r,ev_i,l)
    type(Instanton), intent(in) :: this
    real(dl), dimension(:), intent(out) :: ev_r,ev_i
    integer, intent(in) :: l

    integer :: npt, ierror
    real(dl) :: norm
    real(dl), dimension(0:this%ord,0:this%ord) :: op
    real(dl), dimension(:), allocatable :: work
    integer :: iwork
    real(dl), dimension(1) :: dummy
    
    ! Assign workspace for solver
    npt = this%ord + 1
    call DGEEV('N','N',npt,op,npt,ev_r,ev_i,dummy,1,dummy,1,dummy,-1,ierror)
    if (ierror == 0) then
       iwork = int(dummy(1))
       allocate(work(1:iwork))
    endif

    call fluc_op(this,op,l,this%dim)
    call DGEEV('N','N',npt,op,npt,ev_r,ev_i,dummy,1,dummy,1,work,iwork,ierror)
  end subroutine get_eigenvalues

  subroutine get_eigenvalues_()
  end subroutine get_eigenvalues_
  
  subroutine get_eigenvalues_dirichlet(this,ev,l,d)
    type(Instanton), intent(in) :: this
    real(dl), dimension(:), intent(out) :: ev
    integer, intent(in) :: l, d
  end subroutine get_eigenvalues_dirichlet
  
  subroutine get_eigenvectors(this,ev_r,ev_i,evec,l)
    type(Instanton), intent(in) :: this
    real(dl), dimension(:), intent(out) :: ev_r,ev_i
    real(dl), dimension(:,:), intent(out) :: evec
    integer, intent(in) :: l

    integer :: npt, ierror
    real(dl) :: norm
    real(dl), dimension(0:this%ord,0:this%ord) :: op
    real(dl), dimension(:), allocatable :: work
    integer :: iwork
    real(dl), dimension(1) :: dummy
      
    ! Assign workspace for solver
    npt = this%ord + 1
    call DGEEV('N','V',npt,op,npt,ev_r,ev_i,dummy,1,evec,npt,dummy,-1,ierror)
    if (ierror == 0) then
       iwork = int(dummy(1))
       allocate(work(1:iwork))
    endif

    call fluc_op(this,op,l,this%dim)
    call DGEEV('N','V',npt,op,npt,ev_r,ev_i,dummy,1,evec,npt,work,iwork,ierror)
  end subroutine get_eigenvectors

  subroutine fluc_op(this,op,l,d)
    type(Instanton), intent(in) :: this
    real(dl), dimension(0:this%ord,0:this%ord), intent(out) :: op
    integer, intent(in) :: l,d
    real(dl) :: norm; integer :: i

    norm = l*(l+this%dim-1._dl)
    call set_bc( (/.false.,.false./) )  ! Check this is correct
    call variation(this%phi,op)
    op = -op
    do i=0,this%ord
       op(i,i) = op(i,i) + norm/this%tForm%xGrid(i)**2
    enddo
  end subroutine fluc_op
  
  subroutine write_eval_evec()
  end subroutine write_eval_evec
  
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
