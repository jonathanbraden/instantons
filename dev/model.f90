module Model
  use constants
  implicit none
  
  integer, parameter :: nfield = 1
  real(dl), parameter :: lambda = 2.
  
  integer, parameter :: ndim = 3

! These hardcoded parameters should be replaced by a computational subroutine
  real(dl), parameter :: tension = 2._dl**1.5
  real(dl), parameter :: delrho = 2._dl
  real(dl), parameter :: rinit = dble(ndim)*tension / delrho

contains
  elemental function potential(phi)
    real(dl), intent(in) :: phi
    real(dl) :: potential
    potential = cos(phi) + lambda*sin(phi)**2
  end function potential

  elemental function vprime(phi)
    real(dl), intent(in) :: phi
    real(dl) :: vprime
    vprime = -sin(phi) + lambda*sin(2._dl*phi)
  end function vprime

  elemental function vdprime(phi)
    real(dl), intent(in) :: phi
    real(dl) :: vdprime
    vdprime = -cos(phi) + 2._dl*lambda*cos(2._dl*phi)
  end function vdprime

  subroutine get_vacuum(fld)
    real(dl), intent(inout) :: fld

    integer, parameter :: maxit=16
    real(dl), parameter :: min_tol=1.e-14
    real(dl) :: vp,vpp,dfld

    print*,"initial field is ",fld
    do l=1,maxit
       vpp = vdprime(fld)
       vp = vprime(fld)
       dfld = -vp/vpp
       fld = fld + dfld
       print*,"new field ",fld
       if (abs(dfld) < min_tol) exit
    enddo
    
    if (l.eq.maxit) then
       print*,"Failed to find local minimum of potential. Adust guess"
       stop
    endif
    print*,"Vacuum is ",fld," derivative is ",vprime(fld)
  end subroutine get_vacuum

  subroutine initial_profile()

  end subroutine initial_profile

end module Model
