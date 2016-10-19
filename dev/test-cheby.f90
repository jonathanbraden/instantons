!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> @author Jonathan Braden
!>         University College London
!>
!> @brief
!> Test my Chebyshev expansion based pseudospectral interpolation and differentiation module.
!> This also functions as a unit testing module for my code as well as a compilation verification
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! TODO: Combine all of the ugly separate subroutines into a single subroutine calling specialised
!       separate functions to reduce code duplication
!       Basically use my Test_FFT module as a blueprint

program test_chebyshev
  use, intrinsic :: iso_c_binding
  use constants
  use Chebyshev
  implicit none

  integer :: order
  type(Chebyshev) :: transform

  order = 50
  call check_gaussian(order,.false.)
  call check_gaussian(order,.true.)
contains

  subroutine check_gaussian(o,endpoints)
    integer, intent(in) :: o
    logical, intent(in) :: endpoints
    real(dl), dimension(0:o) :: test_func
    real(dl), dimension(0:o) :: deriv, deriv_n
    real(dl), parameter :: sigma2_inv = 1._dl/4._dl

    call create_chebyshev(transform,o,2,endpoints,.false.)
    test_func = exp(-0.5_dl*sigma2_inv*transform%xGrid(:)**2)

    ! Test first derivative
    deriv = -sigma2_inv*transform%xGrid(:)*test_func(:)
    deriv_n = matmul(transform%derivs(:,:,1),test_func(:))
    print*,"Maximal error in first derivative is ", maxval(abs(deriv-deriv_n))

    ! Test second derivative
    deriv = -sigma2_inv*test_func(:) + sigma2_inv**2*transform%xGrid(:)**2*test_func(:)
    deriv_n = matmul(transform%derivs(:,:,2),test_func(:))
    print*,"Maximal error in second derivative is ", maxval(abs(deriv-deriv_n))

    call destroy_chebyshev(transform)
  end subroutine check_gaussian

! Use the above subroutine as a template
  subroutine check_sinewave(o,endpoints)
    integer, intent(in) :: o
    logical, intent(in) :: endpoints
  end subroutine check_sinewave

end program test_chebyshev
