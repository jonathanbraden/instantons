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

!!!!! Don't repeat this.  Use either a header or just call the function
#define POTENTIAL(f) ( 2._dl*beta*f**2*( (2.*beta+1._dl-d_space)/(2._dl*beta+1._dl)*abs(f)**(1._dl/beta) - abs(f)**(2._dl/beta) ) )
#define VPRIME(f) ( 2._dl*f*( (2._dl*beta+1._dl-d_space)*abs(f)**(1._dl/beta) - 2._dl*(beta+1._dl)*abs(f)**(2._dl/beta) ) )
#define VDPRIME(f) ( 2._dl*(beta+1._dl)/(beta)*( (2._dl*beta+1._dl-d_space)*abs(f)**(1._dl/beta) - 2._dl*(beta+2._dl)*abs(f)**(2._dl/beta) ) )

module Equations
  use constants, only : dl
  use Cheby
  use Model
  implicit none
  private :: L0, neumann, bc  ! streamline this
  
  real(dl), dimension(:,:), allocatable :: L0
  real(dl), dimension(:), allocatable :: neumann
  logical, dimension(1:2) :: bc
  
contains

  subroutine initialise_equations(tForm, params, dim, bc_)
    type(Chebyshev), intent(in) :: tForm
    real(dl), dimension(:), intent(in) :: params
    integer, intent(in) :: dim
    logical, dimension(1:2), intent(in), optional :: bc_
    integer :: i, sz

    ! This line is awful, set it somewhere else
    call set_model_params(params,dim)
    
    sz = size(tForm%xGrid); ndim = dim  ! why am I storing ndim?
    if (allocated(L0)) deallocate(L0); allocate( L0(1:sz,1:sz) )
    if (allocated(neumann)) deallocate(neumann); allocate( neumann(1:sz) )
    bc = .false.; if (present(bc_)) bc = bc_

    ! Bad if this contains the origin  !!! FIX IT
    do i=0,sz-1
       L0(i+1,:) = tForm%derivs(i,:,2) + dble(dim)*tForm%derivs(i,:,1)/tForm%xGrid(i)
    enddo
    neumann(:) = tForm%derivs(0,:,1)
  end subroutine initialise_equations

  subroutine source(fld,src)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:), intent(out) :: src
    integer :: sz

    sz = size(fld)
    src(:) = -matmul(L0,fld)
    src(:) = src(:) + VPRIME(fld(:))  ! Inlining problem
    if ( bc(2) ) src(sz) = 1._dl
    if ( bc(1) ) src(1) = 1._dl
  end subroutine source

  subroutine variation(fld,var)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:,1:), intent(out) :: var
    integer :: i,sz
    
    sz = size(fld)
    var(1:sz,1:sz) = L0(1:sz,1:sz)
    do i=1,sz
       var(i,i) = var(i,i) - ( VDPRIME(fld(i)) )
    enddo
    if ( bc(1) ) then; var(1,:) = neumann(:); endif
    if ( bc(2) ) then; var(sz,:) = 0._dl; var(sz,sz) = 1._dl; endif
  end subroutine variation
    
end module Equations

 
