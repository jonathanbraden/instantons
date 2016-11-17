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
  end subroutine initialise_equations

! Preprocessors for inlining
!#define POTENTIAL(f) ( 0.25_dl*(f**2-1._dl)**2 + del*(f**3/3._dl-f) )
!#define VPRIME(f) ( (f+del)*(f**2-1._dl) )
!#define VDPRIME(f) ( 3._dl*f**2 - 1._dl + 2._dl*del*f )
#define POTENTIAL(f) ( cos(f) + del*sin(f)**2 + 1._dl )
#define VPRIME(f) ( -sin(f) + del*sin(2._dl*f) )
#define VDPRIME(f) ( -cos(f) + 2._dl*del*cos(2._dl*f) )


  elemental function potential(phi)
    real(dl) :: potential
    real(dl), intent(in) :: phi
    potential = POTENTIAL(phi) !0.25_dl*(phi**2-1._dl)**2 + del*(phi**3/3._dl - phi)
  end function potential

  elemental function vprime(phi)
    real(dl) :: vprime
    real(dl), intent(in) :: phi
    vprime =  VPRIME(phi) !(phi+del)*(phi**2 - 1.)
  end function vprime

  elemental function vdprime(phi)
    real(dl) :: vdprime
    real(dl), intent(in) :: phi
    vdprime =  VDPRIME(phi) !3.*phi**2 - 1. + 2.*del*phi
  end function vdprime
  
  subroutine source(fld,src)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:), intent(out) :: src
    integer :: sz

    sz = size(fld)
    src(:) = -matmul(L0,fld)
    src(:) = src(:) + VPRIME(fld(:)) !+ (fld(:)+del)*(fld(:)**2-1._dl)
    src(sz) = 0._dl  ! Set boundary condition at infinity
  end subroutine source

  subroutine variation(fld,var)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:,1:), intent(out) :: var
    integer :: i, sz

    sz = size(fld)
    var(1:sz,1:sz) = L0(1:sz,1:sz)
    do i=1,sz
       var(i,i) = var(i,i) - VDPRIME(fld(i)) !- (3._dl*fld(i)**2-1._dl +2._dl*del*fld(i))
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