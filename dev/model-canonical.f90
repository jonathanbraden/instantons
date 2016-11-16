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

  real(dl), dimension(:,:), allocatable :: L0
  real(dl) :: del

contains
  subroutine initialise_equations(tForm, delta)
    type(Chebyshev), intent(in) :: tForm
    real(dl), intent(in) :: delta
    integer :: i, sz
    sz = size(tForm%xGrid)
    allocate( L0(1:sz,1:sz) )
    do i=0,sz-1
       L0(i+1,:) = tForm%derivs(i,:,2) + 3._dl*tForm%derivs(i,:,1)/tForm%xGrid(i)
    enddo
    del = delta
    ! In the original code, there's a multiplication by the MMT.  Make sure that this is already included in my definition of derivs
  end subroutine initialise_equations

!
! Define the model through it's potential
!
  elemental function potential(phi)
    real(dl) :: potential
    real(dl), intent(in) :: phi
    potential = 0.25_dl*(phi**2-1._dl)**2 + del*(phi**3/3._dl - phi)
  end function potential

  elemental function vprime(phi)
    real(dl) :: vprime
    real(dl), intent(in) :: phi
    vprime =  (phi+del)*(phi**2 - 1.)
  end function vprime

  elemental function vdprime(phi)
    real(dl) :: vdprime
    real(dl), intent(in) :: phi
    vdprime =  3.*phi**2 - 1. + 2.*del*phi
  end function vdprime
  
! Preprocessor for inlining
#define VPRIME(f) ((f+del)*(f**2-1._dl))
  subroutine source(fld,src)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:), intent(out) :: src
    integer :: sz

    sz = size(fld)
    src(:) = -matmul(L0,fld)
    src(:) = src(:) + (fld(:)+del)*(fld(:)**2-1._dl)  !VPRIME(fld(:))
    src(sz) = 0._dl  ! Set boundary condition at infinity
  end subroutine source

#define VDPRIME(f) (3._dl*f**2 - 1._dl + 2._dl*del*f)
  subroutine variation(fld,var)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:,1:), intent(out) :: var
    integer :: i, sz

    sz = size(fld)
    var(1:sz,1:sz) = L0(1:sz,1:sz)
    do i=1,sz
       var(i,i) = var(i,i) - (3._dl*fld(i)**2-1._dl +2._dl*del*fld(i))  !VDPRIME(fld(i))
    enddo
    ! boundary condition at infinity
    var(sz,:) = 0._dl
    var(sz,sz) = 1._dl
  end subroutine variation

  !>@brief
  !> A subroutine to set general Robin boundary conditions on our fields
  !>  \f[
  !>    \alpha_L f(x_L) + \beta_L f'(x_L) = c_L
  !>  \f]
  !>  \f[
  !>    \alpha_R f(x_R) + \beta_R f'(x_R) = c_R
  !>  \f]
  subroutine boundaries(L,S,c,bc)
    real(dl), intent(inout) :: L(1:,1:), S(1:)
    real(dl), dimension(1:3,1:2), intent(in) :: c
    logical, dimension(1:2), intent(in) :: bc
    integer :: sz
    sz = size(S)

    if (bc(1)) then
       L(:,:) = c(1,1) + c(2,1)
       S(1) = c(3,1)
    endif
    if (bc(2)) then
       L(:,:) = c(1,2) + c(2,2)
       S(sz) = c(3,2)
    endif
  end subroutine boundaries

end module Model
