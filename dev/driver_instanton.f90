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
  real(dl), dimension(:,:), allocatable :: L_lin, L_inv
  integer :: u, j, uu
  complex(dl), parameter :: iImag = (0._dl,1._dl)

  integer, dimension(1) :: i1

  ! Only for log-potential, comment if not using it
!  real(dl), dimension(1:nEps+nEPS_r) :: epsScan
!  epsScan(1:nEps) = scanVals(:); epsScan(nEps+1:nEps+nEps_r) = scanVals_r(:)
  
!  call create_instanton(inst,100,3.)
!  call compute_profile_(inst,(/ 0.5 /),out=.true.,p_i=5)
!  call compute_profile(inst,(/1.+10./),get_grid_params,get_guess_params,grid_type='FULL_MID',out=.true.)
!  print*,compute_action(inst)

  call scan_dimensions('bec-norm',100,scanVals, get_grid_params, get_guess_params, use_previous)
!  call scan_profiles(epsScan,2.,100,get_grid_params,get_guess_params,use_previous,out=.true.)
!  call scan_profiles(scanVals,2.,100,get_grid_params,get_guess_params,use_previous,out=.true.)
!  call scan_profiles_(scanVals,1.,100,.true.,p_guess)
  
!  call compute_profile_(inst,(/ 2./10._dl**2 /), out=.true.,p_i=3) ! Drummond new normalisation
!  call compute_profile_(inst,(/ 0.5*1.2_dl**2 /),out=.true.,p_i=3) ! Drummond old normalisation
!  call compute_profile_(inst,(/0.001/),out=.true.,p_i=2)    ! Cubic double well
!  call compute_profile_(inst,(/ 1.+1.e1 /),out=.true.,p_i=5)    ! Fubini Potential
!  call compute_profile_(inst,(/ 1.5 /),out=.true.,p_i=4)  ! Logarithmic potential

#ifdef MORE
  call compute_profile_(inst,(/ 0.5*1.5_dl**2 /),out=.true.,p_i=3)
  call compute_profile_(inst,(/ 0.5*1.6_dl**2 /),out=.true.,p_i=3)
  call compute_profile_(inst,(/ 0.5*6._dl**2 /),out=.true.,p_i=3)
!  call interp_uniform(inst,256,50.*2.**0.5)
!  call interp_simulation(inst,2048,16._dl*4.)
  
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
#endif
#ifdef LINOP
  allocate(L_lin(0:inst%ord,0:inst%ord))
  call linear_operator(L_lin,inst,100,2,vdprime(0._dl))  ! improve this
  allocate(L_inv(0:inst%ord,0:inst%ord))
  call invert_matrix(L_lin,L_inv,inst%ord+1)
  open(unit=newunit(uu),file='lin-inv.dat')
  do j=0,inst%ord
     write(uu,*) L_inv(j,:)
  enddo
  write(uu,*)
  do j=0,inst%ord
     write(uu,*) L_lin(j,:)
  enddo
  write(uu,*)
  call variation(inst%phi,L_lin)
  do j=0,inst%ord
     write(uu,*) L_lin(j,:)
  enddo
#endif
  
!  call scan_resolutions((/ 2./2._dl**2 /),2)

!  dVals = ([ (0.5*(1.02+0.02*(i-1))**2 ,i=1,size(dVals)) ])
!  call scan_eigens(dVals,63,2)
  
contains
  !!
  ! Add the three functions as parameters to pass in
  !!
  !>@brief
  !> Do full scans over parameter values for varying numbers of spatial dimensions.
  !>
  !>@param[in] modName (String) - name of model to use to name files
  !>@param[in] ord (Integer) - Order of interpolating polynomials
  !>@param[in] scanParams (real array) - Values of model parameter to scan over
  subroutine scan_dimensions(modName,ord,scanParams,get_grid,get_ic_params,prev_test)
    character(*), intent(in) :: modName
    integer, intent(in) :: ord
    real(dl), dimension(:), intent(in) :: scanParams
    procedure(get_grid_params) :: get_grid
    procedure(get_guess_params) :: get_ic_params
    procedure(use_previous) :: prev_test
    
    integer :: i
    real(dl) :: dim
    character(6) :: dimStr; character(3) :: ordStr

    do i=0,30
       dim = 0.5+i*0.1
       write(dimStr,'(A1,F3.1,A1)') '-',dim,'d'
       write(ordStr,'(I3)') ord

       call scan_profiles(scanParams,dim,ord,get_grid,get_ic_params,prev_test,out=.false.)
       call execute_command_line('mv actions.dat actions-'//trim(adjustl(modName))//trim(adjustl(dimStr))//'-o'//trim(adjustl(ordStr))//'.dat')
    enddo
  end subroutine scan_dimensions
  
  !! TO DO: Put the interface for get_grid_params and get_guess_params somewhere besides the module
  subroutine scan_profiles(deltas,dim,ord,get_grid,get_ic_params,prev_test,grid_type,out)
    real(dl), dimension(:), intent(in) :: deltas
    real(dl), intent(in) :: dim
    integer, intent(in) :: ord
    procedure(get_grid_params) :: get_grid
    procedure(get_guess_params) :: get_ic_params
    procedure(use_previous) :: prev_test
    character(8), intent(in), optional :: grid_type
    logical, intent(in), optional :: out
    
    type(Instanton) :: inst
    type(Solver) :: solv
    real(dl) :: dCur
    logical :: outL
    integer :: i, fNum
    real(dl) :: len, w
    real(dl), dimension(1:4) :: params_ic
    integer :: p_i
    real(dl), dimension(0:ord) :: xNew, phi_prev
    character(8) :: grid_
    real(dl) :: phizero
    
    grid_= 'FULL_MID'; if (present(grid_type)) grid_ = grid_type
    outL = .false.; if (present(out)) outL = out
    open(unit=newunit(fNum),file='actions.dat')
    call create_instanton(inst,ord,dim)
    call create_solver(solv,ord+1,100,0.1_dl)

    do i=1,size(deltas)
       dCur = deltas(i)
       call set_model_params((/dCur/),dim)

       call get_grid( (/dCur/),dim,len,w )
       call get_ic_params( (/dCur/),dim,params_ic,p_i )
       ! TO DO: Fix the ugly call to get the new chebyshev grid
       ! Optional : Don't even bother to interpolate solution onto the new grid, then I don't have to do anything
       !
       ! Problem : The interpolation seems broken for the Drummond potential
       if ( prev_test( (/dCur/),dim) ) then
          xNew = chebyshev_grid(inst%ord,len,w) ! Fix this to allow for different grids
          phi_prev = interpolate_instanton_(inst,xNew)
          call create_instanton_grid(inst,grid_,len,w)
          inst%phi(:) = phi_prev(:)
       else
          call create_instanton_grid(inst,grid_,len,w)
          call profile_guess(inst,params_ic(1),params_ic(2),params_ic(3),params_ic(4),p_i)
       endif       
       
       call initialise_equations(inst%tForm,(/dCur/),dim,inst%bc) ! Can I add this into the grid creation?

       ! Now solve
       call solve(solv,inst%phi)
       write(fNum,*) dCur, compute_action(inst), interpolate_instanton_(inst,(/0._dl/)), inst%phi(0)
       if (outL) call output_instanton(inst)
    enddo

    close(fNum); call delete_solver(solv); call destroy_instanton(inst)
  end subroutine scan_profiles

  subroutine scan_resolutions(par,d)
    real(dl), dimension(:), intent(in) :: par
    real(dl), intent(in) :: d
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
  
!  subroutine scan_grid_params(inst,ic_params)
!    type(Instanton), intent(inout) :: inst
!    procedure(get_ic_params) :: ic_params
!  end subroutine scan_grid_params

  !
  ! This is old, I can probably delete it
  !
  !>@todo Write this to include the adjustment of the initial bubble profile
  !>@todo Remove need for chebyshev_grid, or put it in the pseudospec module
  !>@todo Allow passing in of check for whether or not to interpolate profile as a function
  subroutine scan_profiles_(deltas,dim,ord,out,prof_type)
    real(dl), dimension(:), intent(in) :: deltas
    real(dl), intent(in) :: dim
    integer, intent(in) :: ord
    logical, intent(in), optional :: out
    integer, intent(in), optional :: prof_type

    type(Instanton) :: inst
    real(dl) :: dCur
    integer :: i, u
    logical :: outL

    real(dl) :: r0, meff, phif, phit
    real(dl) :: len, w
    real(dl), dimension(0:ord) :: xNew, phi_prev
    integer :: p_i

    p_i = 3; if (present(prof_type)) p_i = prof_type
    outL = .false.; if (present(out)) outL = out
    open(unit=newunit(u),file='actions.dat')
    call create_instanton(inst,ord,dim)

    dCur = deltas(1)
    call compute_profile_(inst,(/dCur/),out=outL,p_i=p_i)
    write(u,*) dCur, compute_action(inst)
    
    do i=2,size(deltas)
       dCur = deltas(i)
       if ( prev_test(inst) ) then
       !if ( prev_test_dw(inst) ) then
          call bubble_parameters_nd_(dCur,dim*1._dl,r0,meff); call grid_params_(w,len,r0,1._dl/meff)  ! Remove this nastiness
          call get_grid_params((/dCur/),dim,len,w)
          xNew = chebyshev_grid(inst%ord,len,w)
          phi_prev = interpolate_instanton_(inst,xNew)
          call compute_profile_(inst,(/dCur/),phi_prev,out=outL)
       else
          call compute_profile_(inst,(/dCur/),out=outL,p_i=p_i)
       endif
       write(u,*) dCur, compute_action(inst)
    enddo
    close(u)
  end subroutine scan_profiles_
  
  ! Super ugly, fix it
  subroutine linear_operator(L_lin,this,ord,dim,m2eff)
    type(Instanton), intent(in) :: this
    real(dl), intent(out) :: L_lin(0:ord,0:ord)
    integer, intent(in) :: ord, dim
    real(dl), intent(in) :: m2eff
    
    integer :: i

    do i=0,ord
       L_lin(i,:) = this%tForm%derivs(i,:,2) + dble(dim)*this%tForm%derivs(i,:,1)/this%tForm%xGrid(i)
    enddo
    do i=0,ord
       L_lin(i,i) = L_lin(i,i) - m2eff
    enddo
  end subroutine linear_operator
  
  subroutine scan_eigens(par,ord,dim)
    real(dl), dimension(:), intent(in) :: par
    integer, intent(in) :: ord
    real(dl), intent(in) :: dim
    
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
  !>@brief
  !> Interpolate the instanton onto a specified simulation grid.
  !> The output is stored in instanton_sim.dat, with the first column being the field values and the second column the time derivative.
  !> Each row corresponds to a single lattice site (so that the x-grid is in units of the lattice spacing)
  !> The instanton is centered in the middle of the simulation grid.
  !>
  !>@param[in] this - the instanton solution to be interpolated (type(Instanton))
  !>@param[in] nlat - the number of lattice sites in the simulation
  !>@param[in] len  - length of the simulation
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

#ifdef IC_FACTOR
  subroutine deform_initial_guess(this)
    type(Instanton), intent(in) :: this
    real(dl), dimension(0:this%ord) :: xNew, phi_prev

    call bubble_parameters_nd_(dCur,dim*1._dl,r0,meff)
    call grid_params_(w,len,r0,1._dl/meff)
    
    xNew = chebyshev_grid(this%ord,len,w) ! Fix this
    phi_prev = interpolate_instanton_(inst,xNew)
    call compute_profile_(this,(/dCur/),phi_prev,out=outL)  ! fix this
  end subroutine deform_initial_guess
#endif
  
  function prev_test(inst) result(prev)
    type(Instanton), intent(in) :: inst
    logical :: prev

    real(dl) :: phif, phit
    call get_minima(phif,phit)
    prev = abs(inst%phi(0)-phif) < 0.99_dl*abs(phif-phit)
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
