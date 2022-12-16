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

contains

  subroutine write_eigenvalues(ev_r,ev_i)
    real(dl), dimension(1:), intent(in) :: ev_r, ev_i
    integer :: u, i

    open(unit=newunit(u),file='eval.dat')
    write(u,*) "# l = "
    write(u,*) "# Real  Imaginary"
    do i=1,size(ev_r)
       write(u,*) ev_r(i), ev_i(i)
    enddo
    close(u)
  end subroutine write_eigenvalues

  subroutine write_eigenvectors(evec,xv)
    real(dl), dimension(1:,1:), intent(in) :: evec
    real(dl), dimension(1:), intent(in) :: xv
    integer :: u, i, j, ne, nl

    open(unit=newunit(u),file='evec.dat')
    ne = size(evec(1,:)); nl = size(evec(:,1))
    write(u,*) "# R     F(R)"
    do i=1,ne
       do j=1,nl
          write(u,*) xv(j), evec(j,i)
       enddo
       write(u,*)
    enddo
    close(u)
  end subroutine write_eigenvectors
  
  !>@brief
  !> Compute the eigenvalues for spherical harmonic l around the instanton profile
  subroutine get_eigenvalues(this,ev_r,ev_i,l,out)
    type(Instanton), intent(in) :: this
    real(dl), dimension(1:this%ord+1), intent(out) :: ev_r,ev_i
    integer, intent(in) :: l
    logical, intent(in), optional :: out

    logical :: out_
    integer :: npt, ierror
    real(dl) :: norm
    real(dl), dimension(0:this%ord,0:this%ord) :: op
    real(dl), dimension(:), allocatable :: work
    integer :: iwork
    real(dl), dimension(1) :: dummy

    out_ = .false.; if (present(out)) out_ = out
    ! Assign workspace for solver
    npt = this%ord + 1
    call DGEEV('N','N',npt,op,npt,ev_r,ev_i,dummy,1,dummy,1,dummy,-1,ierror)
    if (ierror == 0) then
       iwork = int(dummy(1))
       allocate(work(1:iwork))
    endif

    call fluc_op(this,op,l,this%dim)
    call DGEEV('N','N',npt,op,npt,ev_r,ev_i,dummy,1,dummy,1,work,iwork,ierror)
    if (out_) call write_eigenvalues(ev_r,ev_i)
  end subroutine get_eigenvalues
  
  subroutine get_eigenvalues_dirichlet(this,ev,l,d)
    type(Instanton), intent(in) :: this
    real(dl), dimension(:), intent(out) :: ev
    integer, intent(in) :: l, d
  end subroutine get_eigenvalues_dirichlet
  
  subroutine get_eigenvectors(this,ev_r,ev_i,evec,l,out)
    type(Instanton), intent(in) :: this
    real(dl), dimension(:), intent(out) :: ev_r,ev_i
    real(dl), dimension(:,:), intent(out) :: evec
    integer, intent(in) :: l
    logical, dimension(1:2), intent(in), optional :: out

    logical, dimension(1:2) :: out_
    integer :: npt, ierror
    real(dl) :: norm
    real(dl), dimension(0:this%ord,0:this%ord) :: op
    real(dl), dimension(:), allocatable :: work
    integer :: iwork
    real(dl), dimension(1) :: dummy

    out_ = (/.false.,.false./); if (present(out)) out_ = out
    
    ! Assign workspace for solver
    npt = this%ord + 1
    call DGEEV('N','V',npt,op,npt,ev_r,ev_i,dummy,1,evec,npt,dummy,-1,ierror)
    if (ierror == 0) then
       iwork = int(dummy(1))
       allocate(work(1:iwork))
    endif

    call fluc_op(this,op,l,this%dim)
    call DGEEV('N','V',npt,op,npt,ev_r,ev_i,dummy,1,evec,npt,work,iwork,ierror)

    if (out_(1)) call write_eigenvalues(ev_r,ev_i)
    if (out_(2)) call write_eigenvectors(evec,this%tForm%xGrid)
  end subroutine get_eigenvectors

!!!!!
!!!!  Update this now that I've included the boundary conditions in the definition of the grid
!!!!!
  subroutine fluc_op(this,op,l,d)
    type(Instanton), intent(in) :: this
    real(dl), dimension(0:this%ord,0:this%ord), intent(out) :: op
    integer, intent(in) :: l
    real(dl), intent(in) :: d  ! Do I actually use this?
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
