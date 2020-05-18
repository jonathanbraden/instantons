!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! INSTANTON SOLVER
!
!>@author
!> Jonathan Braden, Canadian Institute for Theoretical Astrophysics
!
!> Driver program to solve instanton profiles for
!> vacuum decay in both the zero-temperature and
!> high temperature limit
!
!> Important TO DO: Write output to some sort of binary format
!> Fits looks easy enough ...
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program Instanton_Solver
  use constants, only : dl, pi, twopi
  use Utils, only : newunit
  use Cheby
  use Instanton_Class
  use Fluctuations
  implicit none
  
  type(Instanton) :: inst
  real(dl), dimension(:), allocatable :: dVals
  integer :: nDel, i

  ! Eigenvalue and eigenvector storage.  Move this to fluctuations.f90
  real(dl), dimension(:), allocatable :: ev_r, ev_i
  real(dl), dimension(:,:), allocatable :: v_r
  integer :: u, j
  complex(dl), parameter :: iImag = (0._dl,1._dl)
  
  call create_instanton(inst,25,2)
!  call compute_profile_(inst,0.5*10._dl**2,out=.true.)  
!  call interp_uniform(inst,2048,50._dl*2._dl**0.5)

  nDel = 65; allocate(dVals(1:nDel))
!  dVals = 0.5_dl +  10.**([ (-7.+0.2*(i-1), i=size(dVals),1,-1) ] )
  dVals = 10.**([ (-7.+0.2*(i-1), i=1,size(dVals)) ] )
  
  call compute_profile_(inst,(/ 0.5*100._dl**2 /),out=.true.,p_i=3)
!  call compute_profile_(inst,(/0.001/),out=.true.,p_i=2)    ! Cubic double well
!  call compute_profile_(inst,(/ 1.+1.e1 /),out=.true.,p_i=5)    ! Fubini Potential
!  call compute_profile_(inst,(/ 1.5 /),out=.true.,p_i=4)  ! Logarithmic potential
  print*, compute_action(inst)
  
  allocate(ev_r(0:inst%ord),ev_i(0:inst%ord))
  allocate(v_r(0:inst%ord,0:inst%ord))
  call get_eigenvectors(inst,ev_r,ev_i,v_r,0)
  print*,"l = 0 : ",minval(ev_r)
  open(unit=newunit(u),file='eval.dat')
  do i=0,inst%ord
     write(u,*) ev_r(i), ev_i(i)
  enddo
  close(u)
  open(unit=newunit(u),file='evec.dat')
  do i=0,inst%ord
     do j=0,inst%ord
        write(u,*) inst%tForm%xGrid(j), v_r(j,i)
     enddo
     write(u,*)
  enddo
  close(u)
  call get_eigenvalues(inst,ev_r,ev_i,1)
  print*,"l = 1 : ",minval(ev_r)
  call get_eigenvalues(inst,ev_r,ev_i,2)
  print*,"l = 2 : ",minval(ev_r)
  
  ! This seems broken at the moment.  Figure out why.  3D vomits and dies horribly
!  call scan_profiles_(dVals,1,10ls0,.false.)
  
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
  
  !>@todo Write this to include the adjustment of the initial bubble profile
  !>@todo Remove need for chebyshev_grid, or put it in the pseudospec module
  subroutine scan_profiles_(deltas,dim,ord,out)
    real(dl), dimension(:), intent(in) :: deltas
    integer, intent(in) :: dim, ord
    logical, intent(in), optional :: out

    type(Instanton) :: inst
    real(dl) :: dCur
    integer :: i, u
    logical :: outL

    real(dl) :: r0, meff, phif, phit
    real(dl) :: len, w
    real(dl), dimension(0:ord) :: xNew, phi_prev
    
    outL = .false.; if (present(out)) outL = out
    open(unit=newunit(u),file='actions.dat')
    call create_instanton(inst,ord,dim)

    dCur = deltas(1)
    call compute_profile_(inst,(/dCur/),out=outL)
    write(u,*) dCur, compute_action(inst)
    
    do i=2,size(deltas)
       dCur = deltas(i)
       if ( prev_test(inst) ) then
          call bubble_parameters_nd_(dCur,dim*1._dl,r0,meff); call grid_params_(w,len,r0,1._dl/meff)
          xNew = chebyshev_grid(inst%ord,len,w)
          phi_prev = interpolate_instanton_(inst,xNew)
          call compute_profile_(inst,(/dCur/),phi_prev)
          !call compute_profile_(inst,dCur,interpolate_instanton_(inst,chebyshev_grid(ord,len,w))) ! This avoids declaring some arrays, but requires more memory assignment
       else
          call compute_profile_(inst,(/dCur/),out=outL)
       endif
       write(u,*) dCur, compute_action(inst)
    enddo
    close(u)
  end subroutine scan_profiles_

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
  
end program Instanton_Solver
