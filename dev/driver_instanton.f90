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

  integer, dimension(1) :: i1
  
  call create_instanton(inst,200,1)

  nDel = 65; allocate(dVals(1:nDel))
  dVals = 0.5_dl +  10.**([ (-7.+0.2*(i-1), i=size(dVals),1,-1) ] )
!  call scan_profiles_(dVals,1,200,.false.)
!  dVals = 10.**([ (-7.+0.2*(i-1), i=1,size(dVals)) ] )

!  call compute_profile_(inst,(/ 2./10._dl**2 /), out=.true.,p_i=3) ! Drummond new normalisation
!  call compute_profile_(inst,(/ 0.5*1.2_dl**2 /),out=.true.,p_i=3) ! Drummond old normalisation
!  call compute_profile_(inst,(/0.001/),out=.true.,p_i=2)    ! Cubic double well
!  call compute_profile_(inst,(/ 1.+1.e1 /),out=.true.,p_i=5)    ! Fubini Potential
  !  call compute_profile_(inst,(/ 1.5 /),out=.true.,p_i=4)  ! Logarithmic potential
!  print*,"Action Components are :"
!  print*, compute_action(inst)
  !call interp_uniform(inst,2048,50._dl*2._dl**0.5)

  call compute_profile_(inst,(/ 0.5*1000._dl**2 /),out=.true.,p_i=3)
  print*,compute_action(inst)
  call interp_uniform(inst,256,50.*2.**0.5)
  call interp_simulation(inst,2048,2.*50.*2.**0.5)
  
  allocate(ev_r(0:inst%ord),ev_i(0:inst%ord))
  allocate(v_r(0:inst%ord,0:inst%ord))
  call get_eigenvectors(inst,ev_r,ev_i,v_r,0,(/.true.,.true./))
  print*,"l = 0 : ",minval(ev_r), minloc(ev_r)
  i1 = minloc(ev_r); ev_r(i1(1)-1) = 10._dl
  print*,"l = 0 : ",minval(ev_r),minloc(ev_r)
  call get_eigenvalues(inst,ev_r,ev_i,1)
  print*,"l = 1 : ",minval(ev_r)
  call get_eigenvalues(inst,ev_r,ev_i,2)
  print*,"l = 2 : ",minval(ev_r)
  
  ! This seems broken at the moment.  Figure out why.  3D vomits and dies horribly
!  call scan_profiles_(dVals,1,10ls0,.false.)
!  call scan_resolutions((/ 2./2._dl**2 /),2)

  dVals = ([ (0.5*(1.02+0.02*(i-1))**2 ,i=1,size(dVals)) ])
!  call scan_eigens(dVals,63,2)
  
contains

  subroutine scan_eigens(par,ord,dim)
    real(dl), dimension(:), intent(in) :: par
    integer, intent(in) :: ord,dim

    type(Instanton) :: inst
    integer :: i,j
    integer :: u1, u2, u3
    real(dl), dimension(0:ord) :: ev_i, ev_r
    real(dl), dimension(0:ord,0:ord) :: v_r
    integer, dimension(1) :: ev_ind
    real(dl) :: mv
    
    open(unit=newunit(u1),file='eigenvalues.dat')
    open(unit=newunit(u2),file='eigenvectors.dat')
    open(unit=newunit(u3),file='actions.dat')
    
    call create_instanton(inst,ord,dim)
    do i=1,size(par)
       call compute_profile_(inst,(/ par(i) /),out=.false.,p_i=3)
       write(u3,*) compute_action(inst)
       call get_eigenvectors(inst,ev_r,ev_i,v_r,0,(/.true.,.true./))
       mv = minval(ev_r); ev_ind = minloc(ev_r)
       ev_r(ev_ind(1)-1) = 1.e8
       write(u1,*) par(i), mv, minval(ev_r)
       do j=0,ord
          write(u2,*) inst%tForm%xGrid(j), v_r(j,ev_ind(1)-1)
       enddo
       write(u2,*) ""
    enddo
  end subroutine scan_eigens
  
  subroutine scan_resolutions(par,d)
    real(dl), dimension(:), intent(in) :: par
    integer, intent(in) :: d
    type(Instanton) :: inst
    integer :: u1,u2
    integer :: i,o

    open(unit=newunit(u1),file="actions.dat"); open(unit=newunit(u2),file="profiles.dat")
    do i=2,200,4
       !o = 2**i-1
       o = i
       call create_instanton(inst,o,d)
       call compute_profile_(inst,par,out=.false.,p_i=3)
       write(u1,*) o, compute_action(inst)
       ! Add interpolation and output of profiles here
    enddo
    close(u1); close(u2)
  end subroutine scan_resolutions

  !>@brief
  !> Interpolate the instanton profile onto a uniform grid fo
  !>
  !>@param[in] this - Instanton profile to interpolate (type(Instanton))
  !>@param[in] nlat - Number of lattice sites with r>=0
  !>@param[in] len  - Half side length of interpolation grid
  subroutine interp_uniform(this,nlat,len)
    type(Instanton), intent(in) :: this
    integer, intent(in) :: nlat
    real(dl), intent(in) :: len

    integer :: u, i, nHalf
    real(dl), dimension(1:nlat) :: r, phi_i
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

  !!! Point nlat/2 is the origin, means I need nlat/2-1 to the left (for even), nlat/2+1 is origin for odd
  subroutine interp_simulation(this,nlat,len)
    type(Instanton), intent(in) :: this
    integer, intent(in) :: nlat
    real(dl), intent(in) :: len

    integer :: u,i,nHalf
    real(dl) :: dx
    real(dl), dimension(1:nlat/2+1) :: phi_i,r

    open(unit=newunit(u), file='instanton_sim.dat')
    dx = len/dble(nlat)

    nHalf = nlat/2
    r = [( (i-1)*dx, i=1,nHalf+1 )]
    phi_i = interpolate_instanton_(this,r)
    
    if (modulo(nlat,2).eq.0) then
       do i=nHalf,2,-1
          !          write(u,*) -r(i), phi_i(i)
          write(u,*) phi_i(i), 0._dl
       enddo
       do i=1,nHalf+1
          !          write(u,*) r(i), phi_i(i)
          write(u,*) phi_i(i), 0._dl
       enddo
    else
       do i=nHalf+1,2,-1
          !          write(u,*) -r(i), phi_i(i)
          write(u,*) phi_i(i), 0._dl
       enddo
       do i=1,nHalf+1
          !          write(u,*) r(i), phi_i(i)
          write(u,*) phi_i(i), 0._dl
       enddo
    endif
  end subroutine interp_simulation

#ifdef SCALE
  !>@brief
  !> 
  subroutine scale_solution(fld,scl)
    real(dl), dimension(:), intent(in) :: fld
    real(dl), intent(in) :: scl
    ! In here write code to first put the scaled instanton on a new grid
    ! Then go ahead and compute the eigenvalues
    ! Is this just a matter of defining a new grid and doing the interpolation?
    ! Can I just scale the coordinates (but what happens to derivatives?)
    ! I think they just scale ...
    ! deriv => scale times deriv

    real(dl), dimension(1:size(fld),1:size(fld)) :: var
    real(dl), dimension(1:size(fld)) :: src
    integer :: i
    real(dl) :: action_scl

    ! action_scl = 1/scl**(dim-1)
    var = L0
    do i=1,sz
       var(i,i) = var(i,i) - vdprime(fld(i))/scl**2
    enddo
    src = -matmul(L0,fld) + vprime(fld(:))/scl**2

    ! Add boundary conditions here (careful that Neumann is scaled correctly
    ! Figure out correct overall scale
    ! Or, I could just scale the potential pieces, then put the overall norm on.
    ! Action scales as 1/s^{d-1} if I scale potential term as 1/s^2, or 1/s^d, KE by s and PE by 1/s
  end subroutine scale_solution
#endif
  
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
    call compute_profile_(inst,(/dCur/),out=outL,p_i=3)
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
