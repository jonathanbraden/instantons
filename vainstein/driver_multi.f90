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
  
  type(Instanton_multi) :: inst
  real(dl), dimension(:,:), allocatable :: phi_prev
  integer :: ndel
  real(dl), dimension(:), allocatable :: dVals
  integer :: i, ord
  
  ndel = 100
  allocate(dVals(1:ndel))
  dVals = [( i*0.0001, i=1,ndel )]

  ord = 100
  call create_instanton_multi(inst,2,ord,2)
  allocate(phi_prev(1:ord+1,1:2))
  call compute_profile_multi(inst,(/0.1_dl,1.e20,1.e-6,1.e-2/),out=.true.)
  phi_prev = inst%phi

!  do i=1,ndel
!     call compute_profile_(inst,dVals(i),out=.true.,phi_init=phi_prev)
!     phi_prev = inst%phi
!  enddo
  
contains

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

end program Instanton_Solver
