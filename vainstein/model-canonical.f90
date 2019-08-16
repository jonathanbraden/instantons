!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> University College London
!>
!>@brief
!> Store the equations of motion to be solved by our nonlinear boundary value solver
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Model
  use constants
  use Cheby
  implicit none
  private :: L0, S0, L4
  private :: del, m2, npow, ndim
  
  real(dl), dimension(:,:), allocatable :: L0, L4
  real(dl), dimension(:), allocatable :: S0
  real(dl) :: del
  real(dl) :: m2
  integer :: npow
  integer :: ndim

contains
  subroutine initialise_equations(tForm, delta, dim)
    type(Chebyshev), intent(in) :: tForm
    real(dl), intent(in) :: delta
    integer, intent(in) :: dim
    integer :: i, sz
    
    sz = size(tForm%xGrid); ndim = dim
    if (allocated(L0)) deallocate(L0); allocate( L0(1:sz,1:sz) )
    do i=0,sz-1
       L0(i+1,:) = tForm%derivs(i,:,2) + dble(dim)*tForm%derivs(i,:,1)/tForm%xGrid(i)
    enddo

    if (allocated(L4)) deallocate(L4); allocate( L4(1:sz,1:sz) )
    L4 = matmul(L0,L0)  ! Check accuracy of this

    if (allocated(S0)) deallocate(S0); allocate( S0(1:sz) )
    S0(:) = external_source(tForm%xGrid(:))

    del = delta; m2 = 1.e-8; npow = 3
  end subroutine initialise_equations
  
  subroutine source(fld,src)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:), intent(out) :: src
    integer :: sz
    real(dl), dimension(1:size(fld)) :: lap

    sz = size(fld)
    lap(:) = matmul(L0,fld); lap = lap**npow
    src(:) = -matmul(L0,fld) + m2*fld(:) + del*matmul(L0,lap)
    src(:) = src(:) + S0(:)
    src(sz) = 0._dl  ! Set boundary condition at infinity
  end subroutine source

  !>@brief
  !> The second variation of our nonlinear equation
  !
  !> Notes on storage.  The first index is ? and the second index is?
  !> The boundary conditions are designed to pick out phi(r_max),
  !> which is then set to 0 with the source term
  subroutine variation(fld,var)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:size(fld),1:size(fld)), intent(out) :: var
    real(dl), dimension(1:size(fld)) :: lap
    integer :: i, sz

    sz = size(fld)
    var(1:sz,1:sz) = L0(1:sz,1:sz)
    lap = matmul(L0,fld); lap = lap**(npow-1)
    do i=1,sz
       var(:,i) = var(:,i) - del*lap(i)*L4(:,i)  ! check ordering
    enddo
    lap = matmul(L0,lap)
    do i=1,sz
       var(:,i) = var(:,i) - del*lap(i)*L0(:,i)  ! check ordering
    enddo
    do i=1,sz
       var(i,i) = var(i,i) - m2
    enddo
    ! boundary condition at infinity
    var(sz,:) = 0._dl; var(sz,sz) = 1._dl
  end subroutine variation

  elemental function external_source(r) result(src)
    real(dl), intent(in) :: r
    real(dl) :: src
    if (r >= 1._dl) then
       src = 0._dl
    else
       src = 0.75_dl*0.5_dl*twopi/(0.25_dl*twopi**2-6._dl)*(cos(0.5*twopi*r)+1._dl)
    endif
  end function external_source
  
end module Model
