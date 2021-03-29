!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author Jonathan Braden
!>        Canadian Institute for Theoretical Astrophysics
!
! TO DO: - Fix initialisation of L0 if the grid includes the origin or infinity
!        - Decide if I only want to initialise neumann, etc if said bcs are selected
!        - Either include a header with the potentials for inlining, or call model
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Equations
  use constants, only : dl
  use Cheby
  use Model
  implicit none
  private :: L0, neumann, bc
  
  real(dl), dimension(:,:), allocatable :: L0
  real(dl), dimension(:), allocatable :: neumann
  logical, dimension(1:2) :: bc
  
contains

  subroutine set_bc(bc_)
    logical, dimension(1:2), intent(in) :: bc_
    bc = bc_
  end subroutine set_bc
  
  subroutine initialise_equations(tForm, params, dim, bc_)
    type(Chebyshev), intent(in) :: tForm
    real(dl), dimension(:), intent(in) :: params
    real(dl), intent(in) :: dim
    logical, dimension(1:2), intent(in), optional :: bc_
    integer :: i, sz, imin

    ! This line is awful, set it somewhere else
    ! call set_model_params(params,dim)  Moved into compute_profile_
    
    sz = size(tForm%xGrid)
    if (allocated(L0)) deallocate(L0); allocate( L0(1:sz,1:sz) )
    if (allocated(neumann)) deallocate(neumann); allocate( neumann(1:sz) )
    bc = .false.; if (present(bc_)) bc = bc_

    ! Bad if this contains the origin  !!! FIX IT
    imin = 0; if (bc_(1)) imin=1
    do i=imin,sz-1
       L0(i+1,:) = tForm%derivs(i,:,2) + dim*tForm%derivs(i,:,1)/tForm%xGrid(i)
    enddo
    neumann(:) = tForm%derivs(0,:,1)
    if (bc_(1)) L0(1,:) = neumann(:)
    if (bc_(2)) then; L0(sz,:) = 0._dl; L0(sz,sz) = 1._dl; endif
  end subroutine initialise_equations

  subroutine source(fld,src)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:), intent(out) :: src
    integer :: sz

    sz = size(fld)
    src(:) = -matmul(L0,fld) + vprime(fld(:))
    if ( bc(2) ) src(sz) = 0._dl
    if ( bc(1) ) src(1) = 0._dl
  end subroutine source

  subroutine variation(fld,var)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:,1:), intent(out) :: var
    integer :: i,sz
    
    sz = size(fld)
    var(1:sz,1:sz) = L0(1:sz,1:sz)
    do i=1,sz
       var(i,i) = var(i,i) - vdprime(fld(i))
    enddo
    if ( bc(1) ) then; var(1,:) = neumann(:); endif
    if ( bc(2) ) then; var(sz,:) = 0._dl; var(sz,sz) = 1._dl; endif
  end subroutine variation
    
end module Equations

 
