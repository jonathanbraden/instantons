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
  
  call create_instanton(inst,100,1)
  call compute_profile_(inst,0.5*1.2_dl**2,out=.true.)

  call interp_uniform(inst,2048,50._dl*2._dl**0.5)
  
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
  
  !>@todo Write this an include the adjustment of the initial bubble profile
  subroutine scan_profiles_(deltas,dim,ord,out)
    real(dl), dimension(:), intent(in) :: deltas
    integer, intent(in) :: dim, ord
    logical, intent(in), optional :: out

    type(Instanton) :: inst
    real(dl) :: dCur
    integer :: i, u
    logical :: outL
    
    outL = .false.; if (present(out)) outL = out
    open(unit=newunit(u),file='actions.dat')
    call create_instanton(inst,ord,dim)

    dCur = deltas(1)
    call compute_profile_(inst,dCur,out=outL)
    write(u,*) dCur, compute_action_(inst)
    
    do i=2,size(deltas)
       dCur = deltas(i)
       if ( .false. ) then
          ! Add code to use previous profile
       else
          call compute_profile_(inst,dCur,out=outL)
       endif
       write(u,*) dCur, compute_action_(inst)
    enddo
    close(u)
  end subroutine scan_profiles_
  
end program Instanton_Solver
