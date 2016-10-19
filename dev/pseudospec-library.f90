!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! MODULE: PSpec
!> @author Jonathan Braden
!>         University College London
!>
!> This provides an object oriented approach to pseudospectral differentiation.
!> With this version of the code, only Chebyshev expansions are supported,
!> with the necessary machinery included in the external module chebyshev.f90.
!> In a future release, the functionality will be extended to include 
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module PSpec
  use, intrinsic :: iso_c_binding
#ifdef USEOMP
  use omp_lib
#endif
  use constants
  implicit none

  type TransformSpec1D
     integer :: nx, ord
     real(C_DOUBLE), allocatable :: realSpace(:)
     real(C_DOUBLE), allocatable :: specSpace(:)
     real(C_DOUBLE), allocatable :: fTrans(:,:)
     real(C_DOUBLE), allocatable :: invTrans(:,:)
     real(C_DOUBLE), allocatable :: derivs(:,:,:)
  end type TransformSpec1D

contains
  
  subroutine initialize_transform_1d(this, n, nd)
    type(TransformSpec1D), intent(out) :: this
    integer, intent(in) :: n, nd

    this%nx = n; this%ord=n-1
    allocate( realSpace(1:nx),specSpace(0:ord) )
    ! Now initialize the transforms, etc.
    allocate( fTrans(1:n,1:n) )
    allocate( invTrans(1:n,1:n) )
    allocate( derivs(1:n,1:n,1:nd) )
  end subroutine initialize_transform_1d

  subroutine destroy_transform_1d(this)
    type(TransformPair1D), intent(inout) :: this

    this%nx = -1; this%ord = -1
    deallocate(realSpace,specSpace,fTrans,invTrans,derivs)
  end subroutine destroy_transform_1d

  subroutine derivativeNth(this)
    type(TransformPair1D), intent(inout) :: this
  end subroutine derivativeNth

! Implement Clenshaw-Curtis algorithm for performance testing
  subroutine clenshawCurtis()
  end subroutine clenshawCurtis
  
end module PSpec
