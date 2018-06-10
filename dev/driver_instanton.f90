!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! INSTANTON SOLVER
!
!>@author
!> Jonathan Braden, University College London
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program Instanton_Solver
  use constants, only : dl, pi, twopi
  use Utils, only : newunit
  use Cheby
  use Instanton_Class

  type(Instanton) :: inst
  real(dl), dimension(:), allocatable :: dVals
  integer :: nDel, i
  
  call create_instanton(inst,100,1)
  call compute_profile_(inst,0.5*1.05_dl**2,out=.true.)

!  call interp_uniform(inst,2048,50._dl*2._dl**0.5)

  nDel = 55; allocate(dVals(1:nDel))
  dVals = 0.5_dl +  10.**([ (-7.+0.2*(i-1), i=size(dVals),1,-1) ] )

  call compute_profile_(inst,dVals(1),out=.true.)
  
  call scan_profiles_(dVals,1,200,.false.)
  
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
    call compute_profile_(inst,dCur,out=outL)
    write(u,*) dCur, compute_action_(inst)
    
    do i=2,size(deltas)
       dCur = deltas(i)
       if ( prev_test(inst) ) then
          call bubble_parameters_nd_(dCur,dim*1._dl,r0,meff); call grid_params_(w,len,r0,1._dl/meff)
          xNew = 0._dl
          !xNew = chebyshev_grid(inst%ord,len,w)
          phi_prev = interpolate_instanton_(inst,xNew)
          call compute_profile_(inst,dCur,phi_prev)
          !call compute_profile_(inst,dCur,interpolate_instanton_(inst,chebyshev_grid(ord,len,w))) ! This avoids declaring some arrays, but requires more memory assignment
       else
          call compute_profile_(inst,dCur,out=outL)
       endif
       write(u,*) dCur, compute_action_(inst)
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
  
end program Instanton_Solver
