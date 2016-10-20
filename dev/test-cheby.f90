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
!
!       Add a subroutine to test the collocation points are correct
!       Test the use of multiple applications of numerical derivatives vs. analytic calculation of second derivative matrix
!

program test_chebyshev
  use, intrinsic :: iso_c_binding
  use constants
  use Cheby
  implicit none

  integer :: order
  type(Chebyshev) :: transform

  order = 150
  open(unit=99,file='tests.dat') ! fix this ugly non-locality
  open(unit=98,file='basis.dat')
  call check_basis(order,.true.)
!  call check_basis_derivs(order,.false.)
  
  call check_gaussian(order,.false.)
  call check_gaussian(order,.true.)
contains
  
  subroutine check_basis(o,endpoints)
    integer, intent(in) :: o
    logical, intent(in) :: endpoints
    integer :: i, j
    real(dl), dimension(0:o) :: test_func, tmp, tmp_c

    open(unit=97,file='basis-test.dat')
    call create_chebyshev(transform,o,2,endpoints,.false.)
    print*,transform%xGrid
    !transform%fTrans(o,:) = 0.5_dl*transform%fTrans(o,:) ! This line fixed the Lobatto and breaks the Gauss
    do i=0,o
       test_func = transform%invTrans(:,i)
!       tmp_c = matmul(transform%fTrans,test_func)
!       tmp = matmul(transform%invTrans,tmp_c)
       call DGEMV('N',(o+1),(o+1),1.d0, transform%fTrans, (o+1), test_func, 1, 0.d0, tmp_c, 1)
       call DGEMV('N',(o+1),(o+1),1.d0, transform%invTrans, (o+1), tmp_c, 1, 0.d0, tmp, 1)
       do j=0,o
          write(97,*) transform%xGrid(j), test_func(j), tmp(j), tmp_c(j)
       enddo
       write(97,*)
       print*,"Maximal error transforming basis function ",i," is ",maxval(abs(tmp-test_func))
       print*,"Maximal analytic error for ",i," is ",maxval(abs(tmp-cos(dble(i)*acos(transform%xGrid))))
       print*,"Maximal error in original func is ",i," is ",maxval(abs(test_func-cos(dble(i)*acos(transform%xGrid))))
    enddo

    call destroy_chebyshev(transform)
    close(unit=97)
  end subroutine check_basis

  !>@todo
  !> Initialise functions analytically to see how the error on derivative matrix propagates
  subroutine check_basis_derivs(o,endpoints)
    integer, intent(in) :: o
    logical, intent(in) :: endpoints
    integer :: i,j
    real(dl), dimension(0:o) :: test_func, tmp, tmp_c
    real(dl), dimension(0:o,0:o) :: dfunc
    integer :: nd

    nd = 1
    open(unit=97,file='basis-deriv.dat')
    call create_chebyshev(transform,o,2,endpoints,.false.)
    dfunc = matmul(transform%derivs(:,:,nd),transform%invTrans)
    if (.not. endpoints) then
       do i=0,o
          test_func = transform%invTrans(:,i)
          tmp = matmul(transform%derivs(:,:,nd),test_func)
          tmp_c = matmul(transform%fTrans,tmp)
          print*,"Maximal ",nd," derivative error is ",maxval(abs(tmp(:)-dfunc(:,i)))
          print*,"Maximal error compared to analytic formula is ",maxval(abs(tmp(:)-dble(i)*sin(i*acos(transform%xGrid(:)))/sqrt(1._dl-transform%xGrid(:)**2)))
          print*,"Maximal error from analytic formula for recurrence is ",maxval(abs(dfunc(:,i)-dble(i)*sin(i*acos(transform%xGrid(:)))/sqrt(1._dl-transform%xGrid(:)**2)))
          do j=0,o
             write(97,*) transform%xGrid(j), dfunc(j,i), tmp(j), tmp_c(j)
          enddo
          write(97,*)
       enddo
    endif

    ! Now test against the analytically initialised basis functions
    if (.not. endpoints) then
       do i=0,o
          test_func = dcos(dble(i)*dacos(transform%xGrid))
          tmp = matmul(transform%derivs(:,:,nd),test_func)
          tmp_c = matmul(transform%fTrans,tmp)
          do j=0,o
             write(97,*) transform%xGrid(j), dfunc(j,i), tmp(j)
          enddo
          write(97,*)
       enddo
    endif
    
    call destroy_chebyshev(transform)
    close(unit=97)
  end subroutine check_basis_derivs
  
  subroutine check_gaussian(o,endpoints)
    integer, intent(in) :: o
    logical, intent(in) :: endpoints
    real(dl), dimension(0:o) :: test_func
    real(dl), dimension(0:o) :: deriv, deriv_n
    real(dl), parameter :: sigma2_inv = 4.*64._dl
    integer :: i,j
    real(dl) :: c(0:1)
    c(0) = 1.; c(1) = 0.5
    
    call create_chebyshev(transform,o,2,endpoints,.false.)
    if (.not. endpoints) then
       do i=0,o
          do j=0,o
             write(98,*) transform%xGrid(j), transform%invTrans(j,i)
          enddo
          write(98,*)
       enddo
    endif
    
    test_func = exp(-0.5_dl*sigma2_inv*transform%xGrid(:)**2)
    
    deriv = matmul(transform%fTrans,test_func)
    deriv_n = matmul(transform%invTrans,deriv)
    do i=0,o
       write(99,*) transform%xGrid(i), test_func(i), deriv(i), deriv_n(i)
    enddo
    write(99,*)
    print*,"Maximal error in double transformation is ", maxval(abs(deriv_n(:)-test_func(:)))
    
    ! Test first derivative
    deriv = -sigma2_inv*transform%xGrid(:)*test_func(:)
!    deriv_n = matmul(transform%fTrans(:,:),test_func(:))
!    deriv_n = matmul(transform%derivs(:,:,1),deriv_n)
    deriv_n = matmul(transform%derivs(:,:,1),test_func(:))
    print*,"Maximal error in first derivative is ", maxval(abs(deriv(:)-deriv_n(:)))
    do i=0,o
       write(99,*) transform%xGrid(i), deriv_n(i), deriv(i)
    enddo
    write(99,*)
    
    ! Test second derivative
    deriv = -sigma2_inv*test_func(:) + sigma2_inv**2*transform%xGrid(:)**2*test_func(:)
    deriv_n = matmul(transform%derivs(:,:,2),test_func(:))
    print*,"Maximal error in second derivative is ", maxval(abs(deriv-deriv_n))
    do i=0,o
       write(99,*) transform%xGrid(i), deriv_n(i), deriv(i)
    enddo
    write(99,*)
    call destroy_chebyshev(transform)
  end subroutine check_gaussian

! Use the above subroutine as a template
  subroutine check_sinewave(o,endpoints)
    integer, intent(in) :: o
    logical, intent(in) :: endpoints
  end subroutine check_sinewave

  subroutine check_polynomial(o,endpoints,coeffs)
    integer, intent(in) :: o
    logical, intent(in) :: endpoints
    real(dl), intent(in) :: coeffs(0:o)
  end subroutine check_polynomial
  
end program test_chebyshev
