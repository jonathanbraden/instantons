!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! INSTANTON SOLVER
!
!>@author
!> Jonathan Braden, University College London
!
!> Driver program to solve instanton profiles for
!> vacuum decay in both the zero-temperature and
!> high temperature limit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program Instanton_Solver
  use constants, only : dl, pi, twopi
  use Utils, only : newunit
  use Cheby
  use Instanton_Class
  implicit none
  
  type(Instanton) :: inst
  real(dl), dimension(:), allocatable :: phi_prev
  integer :: ndel
  real(dl), dimension(:), allocatable :: dVals
  integer :: i
  
  ndel = 100
  allocate(dVals(1:ndel))
  dVals = [( i*0.0001, i=1,ndel )]

  
  call create_instanton(inst,100,2)
  allocate(phi_prev(1:101))
  call compute_profile_(inst,0._dl,out=.true.)
  phi_prev = inst%phi

!  do i=1,ndel
!     call compute_profile_(inst,dVals(i),out=.true.,phi_init=phi_prev)
!     phi_prev = inst%phi
!  enddo
  
contains

  !>@brief
  !> Interpolate the instanton profile onto a uniform grid for external simulation
  subroutine interp_uniform(this,nlat,len)
    type(Instanton), intent(in) :: this
    integer, intent(in) :: nlat
    real(dl), intent(in) :: len

    integer :: u, i, n
    real(dl), dimension(1:nlat) :: r, phi_i
    n = (nlat+1)/2
    
    r = [ ( (i-1)*(len/nlat), i=1,nlat) ]
    phi_i = interpolate_instanton_(this,r)
    
    open(unit=newunit(u), file='instanton_interp_.dat')
    do i=nlat,2,-1
       write(u,*) -r(i), phi_i(i)
    enddo
    do i=1,nlat
       write(u,*) r(i), phi_i(i)
    enddo
    close(u)
  end subroutine interp_uniform

  function prev_test(inst) result(prev)
    type(Instanton), intent(in) :: inst
    logical :: prev

    real(dl) :: phif, phit
    call get_minima(phif,phit)
    prev = abs(inst%phi(0)-phif) < 0.9_dl*abs(phif-phit)
  end function prev_test

  !>@brief
  !> Return the collocation grid for even Chebyshevs on the doubly-infinite interval with clustering.
  !> Used by my scan program
  function chebyshev_grid(order,L,w) result(xVals)
    integer, intent(in) :: order
    real(dl), intent(in) :: L,w
    real(dl), dimension(0:order) :: xVals
    integer :: i; real(dl) :: dkcol

    ! Replace with a function call
    dkcol = 0.5_dl*pi / dble(order+1)
    do i=0,order
       xVals(i) = -cos( dble(2*i+1)*dkcol )
    enddo
    xVals = sqrt(0.5_dl*xVals + 0.5_dl)
    xVals = (1._dl/pi)*atan(w*tan(pi*(xVals-0.5_dl))) + 0.5_dl
    xVals = L*xVals / sqrt(1._dl - xVals**2)
  end function chebyshev_grid

  ! Move this somewhere else
  
  ! Compute initial guess
  subroutine profile_guess_(this,delta,npow)
    type(Instanton), intent(inout) :: this
    real(dl), intent(in) :: delta
    integer, intent(in) :: npow

    ! First get W from Poisson equation lap W = -rho
    ! call DGESV()
    this%phi(:) = this%phi(:)**(1._dl/npow)
    ! Solve Poisson equation to get initial Pi guess lap Pi = W
    ! call DGESV()
    this%phi(:) = this%phi(:) / delta**(1._dl/npow)
  end subroutine profile_guess_

end program Instanton_Solver
