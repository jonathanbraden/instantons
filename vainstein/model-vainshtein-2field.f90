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
  private :: L0, S0, L2, neumann
  private :: m2_lt, m2_hvy, ndim, lam, alpha
  
  real(dl), dimension(:,:), allocatable :: L0, L2
  real(dl), dimension(:), allocatable :: neumann
  real(dl), dimension(:), allocatable :: S0
  real(dl) :: alpha
  real(dl) :: m2_lt, m2_hvy, lam
  integer :: ndim
  integer, parameter :: nfld = 2
  
contains

  subroutine initialise_equations(tForm, params, dim)
    type(Chebyshev), intent(in) :: tForm
    real(dl), dimension(1:4), intent(in) :: params
    integer, intent(in) :: dim
    call initialise_equations_gauss(tForm,params,dim)
  end subroutine initialise_equations
  
  subroutine initialise_equations_gauss(tForm, params, dim)
    type(Chebyshev), intent(in) :: tForm
    real(dl), dimension(1:4), intent(in) :: params
    integer, intent(in) :: dim
    integer :: sz, n
    integer :: i

    n = size(tForm%xGrid); sz = nfld*size(tForm%xGrid)
    alpha = params(1); lam = params(2); m2_lt = params(3); m2_hvy = params(4)
    ndim = dim ! This should just be 2 here
    
    ! Compute radial Laplacian
    if (allocated(L2)) deallocate(L2); allocate( L2(1:n,1:n) )
    do i=0,n-1
       L2(i+1,:) = tForm%derivs(i,:,2) + dble(dim)*tForm%derivs(i,:,1)/tForm%xGrid(i)
    enddo

    ! Compute linear part of variation operator
    if (allocated(L0)) deallocate(L0); allocate( L0(1:sz,1:sz) )
    do i=1,nfld
       L0((i-1)*n+1:i*n,(i-1)*n+1:i*n) = L2(:,:)
    enddo
    L0(1:n,n+1:2*n) = -alpha*L2(:,:)
    L0(n+1:2*n,1:n) = -alpha*L2(:,:)
    do i=1,n
       L0(i,i) = L0(i,i) - m2_lt
       L0(i+n,i+n) = L0(i+n,i+n) - m2_hvy
    enddo
    
    if (allocated(S0)) deallocate(S0); allocate( S0(1:sz) )
    S0 = 0._dl
    S0(1:n) = external_source(tForm%xGrid(:))
  end subroutine initialise_equations_gauss

  subroutine initialise_equations_lobatto(tForm, params, dim)
    type(Chebyshev), intent(in) :: tForm
    real(dl), dimension(1:4), intent(in) :: params
    integer, intent(in) :: dim
    integer :: i, sz, n

    n = size(tForm%xGrid); sz = nfld*size(tForm%xGrid)
    alpha = params(1); lam = params(2); m2_lt = params(3); m2_hvy = params(4)
    ndim = dim
    
    ! Compute radial Laplacian
    if (allocated(L2)) deallocate(L2); allocate( L2(1:n,1:n) )
    do i=1,n-1
       L2(i+1,:) = tForm%derivs(i,:,2) + dble(dim)*tForm%derivs(i,:,1)/tForm%xGrid(i)
    enddo
    L2(1,:) = 0._dl

    ! For Neumann boundary condition at origin
    if (allocated(neumann)) deallocate(neumann); allocate( neumann(1:n) )
    neumann = tForm%derivs(0,:,1)
    
    if (allocated(L0)) deallocate(L0); allocate( L0(1:sz,1:sz) )
    do i=1,nfld
       L0((i-1)*n+1:i*n,(i-1)*n+1:i*n) = L2(:,:)
    enddo
    ! Now add cross terms
    L0(1:n,n+1:2*n) = -alpha*L2(:,:)
    L0(n+1:2*n,1:n) = -alpha*L2(:,:)
    do i=1,n
       L0(i,i) = L0(i,i) - m2_lt
       L0(i+n,i+n) = L0(i+n,i+n) - m2_hvy
    enddo
    ! Add Neumann boundary condition
    ! Write this
    
    if (allocated(S0)) deallocate(S0); allocate( S0(1:sz) )
    S0 = 0._dl
    S0(1:n) = external_source(tForm%xGrid(:))
  end subroutine initialise_equations_lobatto
  
  subroutine source(fld,src)
    real(dl), dimension(1:), intent(in) :: fld  ! Fix shape
    real(dl), dimension(1:), intent(out) :: src
    integer :: i, sz, n

    sz = size(fld); n = sz / nfld

    src = -matmul(L0,fld) + S0(:)
    ! Nonlinear term
    src(n+1:2*n) = src(n+1:2*n) + (1._dl/6._dl)*lam*fld(n+1:2*n)**3

    ! Dirichlet boundary condition at infinity
    do i=1,nfld
       src(i*n) = 0._dl
    enddo
 !   do i=1,nfld
 !      src((i-1)*n+1) = 0._dl  ! Neumann at origin
 !   enddo
  end subroutine source

  !>@brief
  !> The second variation of our nonlinear equation
  subroutine variation(fld,var)
    real(dl), dimension(1:), intent(in) :: fld 
    real(dl), dimension(1:size(fld),1:size(fld)), intent(out) :: var 
    integer :: i, sz, n

    sz = size(fld); n = size(fld)/nfld 
    var(1:sz,1:sz) = L0(1:sz,1:sz)

    ! Now put in the nonlinear terms
    do i=1,n
       var(i+n,i+n) = var(i+n,i+n) - 0.5_dl*lam*fld(i+n)**2
    enddo
       
    ! Dirichlet boundary condition at infinity
    do i=1,nfld
       var(i*n,:) = 0._dl; var(i*n,i*n) = 1._dl
    enddo
    ! Neumann boundary condition at origin
!    do i=1,nfld
!       var((i-1)*n+1,:) = 0._dl
!       var((i-1)*n+1,(i-1)*n+1:i*n) = neumann(:)
!    enddo
  end subroutine variation

  elemental function external_source(r) result(src)
    real(dl), intent(in) :: r
    real(dl) :: src
    if (r >= 1._dl) then
       src = 0._dl
    else
       src = 0.75_dl*0.5_dl*twopi/(0.25_dl*twopi**2-6._dl)*(cos(0.5_dl*twopi*r)+1._dl)
    endif
  end function external_source

  subroutine boundary_conditions(var,src)
    real(dl), dimension(:,:), intent(inout) :: var
    real(dl), dimension(:), intent(inout) :: src

    integer :: i,n
    ! Find a nicer way to smoothly determine grid size
    n = size(src)/nfld
    
    ! Dirichlet boundary condition at infinity
    do i=1,nfld
       var(i*n,:) = 0._dl; var(i*n,i*n) = 1._dl
       src(i*n) = 0._dl
    enddo
    ! Neumann boundary condition
    do i=1,nfld
       var((i-1)*n+1,:) = 0._dl; var((i-1)*n+1,(i-1)*n+1:i*n) = neumann(:)
       src((i-1)*n+1) = 0._dl
    enddo
  end subroutine boundary_conditions
  
end module Model
