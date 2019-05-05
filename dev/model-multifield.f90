!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> University College London
!>
!>@brief
!> Store the equations of motion to be solved by our nonlinear boundary value solver
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Model
  use constants
  use Cheby
  implicit none
  private :: L0, L0_base, dim, nfld, params
  
  real(dl), dimension(:,:), allocatable :: L0, L0_base
  real(dl), dimension(:), allocatable :: params
  integer :: nfld
  integer :: dim
  
contains

  subroutine initialise_equations(tForm, par, dim_, nf)
    type(Chebyshev), intent(in) :: tForm
    real(dl), dimension(:), intent(in) :: par
    integer, intent(in) :: dim  ! Change to a real number
    integer, intent(in) :: nf
    integer :: i,j, sz, nx

    dim = dim_; params = par; nfld = nf
    nx = size(tForm%xGrid); sz = nfld*nx
    allocate( L0(1:sz,1:sz) ); allocate( L0_base(1:nx,1:nx) )

    L0_base(:,:) = tForm%derivs(:,:,2)
    do i=0,nx-1
       L0_base(i+1,:) = L0_base(i+1,:) + dble(dim)*tForm%derivs(i,:,1)/tForm%xGrid(i)
    enddo
    L0 = 0._dl
    do i=1,nfld
       L0((i-1)*nx+1:i*nx,(i-1)*nx+1:i*nx) = L0_base(:,:)
    enddo
  end subroutine initialise_equations

  subroutine source(fld,src)
    real(dl), dimension(1:,1:nfld), intent(in) :: fld
    real(dl), dimension(1:nx*nfld), intent(out) :: src
    integer :: i, m, nx
    real(dl), dimension(1:nfld) :: dv
    nx = size(fld(:,1))
    do i=1,nfld
       src((i-1)*nx+1,i*nx) = -matmul(L0_base,fld(:,i))
    enddo
    do m=1,nx
       dv = vprime(m)
       do i=1,nfld
          src((i-1)*nx+m) = src((i-1)*nx+m) + vp(i)  ! replace with stride as below
       enddo
       ! Test this to replace the loop above
       !src((i-1)*nx+m:m+nx*(nfld-1):nx) = src(m:m+nx*(nfld-1):nx) + dv(1:nfld)
    enddo
    ! Boundary condition at infinity
    do i=1,nfld
       src(i*nx) = 0._dl
    enddo
    ! Test this to replace loop as well
    !src(nx:nx*nfld:nx) = 0._dl
  end subroutine source

  subroutine variation(fld,var)
    real(dl), dimension(1:,1:nfld), intent(in) :: fld
    real(dl), dimension(1:,1:), intent(out) :: var
    integer :: i,j, nx, m
    real(dl), dimension(1:nfld,1:nfld) :: d2v
    
    nx = size(fld(:,1))
    var = L0
    do m=1,nx
       d2v = vdprime(fld(m,:))
       ! Replace this with a nicer slicing operator
       do j=1,nfld; do i=1,nfld
          var( (i-1)*nx+1:i*nx,(j-1)*nx+1:j*nx ) = var( (i-1)*nx+1:i*nx,(j-1)*nx+1:j*nx ) + d2v(i,j)
       enddo; enddo
    enddo
    ! Set boundary conditions at infinity (try merging with above)
    ! Also replace this with a nice slicing operator if possible
    do i=1,nfld
       var(i*nx,:) = 0._dl
       var(i*nx,i*nx) = 1._dl
    enddo
  end subroutine variation

  ! Factor these out into a separate file (since the the functions above are generic as long as the kinetic term is canonical
  function potential(phi)
    real(dl) :: potential
    real(dl), dimension(:), intent(in) :: phi
    potential = 0.25_dl*(phi(1)**2-1._dl)**2 + 0.25_dl*(phi(2)**2-1._dl)**2
  end function potential
  
  function vprime(phi)
    real(dl), dimension(1:nfld) :: vprime
    real(dl), dimension(1:nfld), intent(in) :: phi

    vprime(1) = (phi(1)**2-1._dl)*phi(1)
    vprime(2) = (phi(2)**2-1._dl)*phi(2)
  end function vprime

  function vdprime(phi,i,j)
    real(dl), dimension(1:nfld,1:nfld) :: vdprime
    real(dl), dimension(1:nfld), intent(in) :: phi

    vdprime(1,1) = 3._dl*phi(1)**2-1._dl
    vdprime(2,2) = 3._dl*phi(2)**2-1._dl
    vdprime(2,1) = 0._dl
    vdprime(1,2) = 0._dl
  end function vdprime
  
end module Model
