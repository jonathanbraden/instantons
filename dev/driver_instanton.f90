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

  real(dl), dimension(1:101) :: del_vals_
  integer, parameter :: n_params = 51 ! 5001  ! adjust this
  real(dl), dimension(1:nPar,1:n_params) :: param_vals
  real(dl), dimension(:,:), allocatable :: params
  real(dl) :: drho
  
 
  integer :: i_, istart

  ! This is for strict Drummond model
  !param_vals = 0.5*1.2_dl**2
  ! call scan_profiles(param_vals(:1,:2),1.,100,get_grid_params,get_guess_params,use_previous,out=.true.)
  ! param_vals(1,1) = 2._dl; param_vals(2,1) = -2._dl; param_vals(3,1) = 0._dl

  ! testing lambda = 1.2 sims
  !param_vals(1,:) = 1._dl; param_vals(2,:) = -0.25*1.2**2; param_vals(3,:) = 0.
  !call scan_profiles(param_vals(:,:2),1.,100,get_grid_params,get_guess_params,use_previous,out=.true.)
  
  ! Only for log-potential, comment if not using it
!  real(dl), dimension(1:nEps+nEPS_r) :: epsScan
!  epsScan(1:nEps) = scanVals(:); epsScan(nEps+1:nEps+nEps_r) = scanVals_r(:)

!  call scan_s2_bec_vals(1.2,0.99,101,del_vals_)
!  call read_param_file('eff_params_l1.2.txt',params,4)
  
!  call read_param_file('eff_params_cosine_l1.2.txt',params,8)
  !  call read_param_file('eff_params_cosine.txt',params,8)
  call read_param_file('params.txt',params,8,n_head=4)
  istart = 0
  do i_ = 1,n_params
     !call convert_masses_to_model(params(5,istart+i_),params(6,istart+i_),param_vals(:,i_)) ! 2 term version
     drho=0.25
     call convert_masses_and_rho_to_model(params(5,istart+i_),params(6,istart+i_),drho,param_vals(:,i_)) ! 3 term fixed rho version
     !param_vals(1,i_) = exp(-0.5*params(7,i_)); param_vals(3,i_) = 0.; param_vals(2,i_) = -0.25*2._dl**2*exp(-2.*params(7,i_))
  enddo
  call scan_profiles(param_vals(:,:),1.,100,get_grid_params,get_guess_params,use_previous,out=.true.)
  ! Using lambda-effective in Gaussian resum
  !param_vals(1,:) = 0.5_dl*params(4,:)

  
!  call scan_profiles(param_vals(:,:1),1.,100,get_grid_params,get_guess_params,use_previous,out=.true.)
  
  !call scan_profiles(sigVals,1.,100,get_grid_params,get_guess_params,use_previous,grid_type='FULL_MID',out=.true.)
  
!  call compute_profile(inst,(/0.09_dl/),get_grid_params,get_guess_params,use_previous,grid_type='FULL_MID',out=.true.)
  
!  call scan_dimensions('dw',25,scanVals,get_grid_params,get_guess_params,use_previous)
!  call scan_dimensions('bec-norm',100,scanVals, get_grid_params, get_guess_params, use_previous)
  
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

  subroutine read_param_file(fName, params,line_length,n_head)
    character(*), intent(in) :: fName
    real(dl), dimension(:,:), allocatable, intent(out) :: params
    integer, intent(in) :: line_length
    integer, intent(in), optional :: n_head

    integer :: n_head_
    integer :: u, i, iostat, nl
    character(100) :: burn
    
    n_head_ = 0; if (present(n_head)) n_head_ = n_head
    
    open(unit=newunit(u),file=fName)
    nl = 0
    do i=1,n_head_; read(u,*); enddo
    do
       read(u,*,iostat=iostat)
       if (iostat /= 0) exit
       nl = nl+1
    enddo
    close(u)

    allocate(params(1:line_length,1:nl)) ! Fix this
    open(unit=newunit(u), file=fName)
    do i=1,n_head_; read(u,*); enddo
    do i=1,nl
       read(u,*) params(:,i)
    enddo
    close(u)
  end subroutine read_param_file
  
  subroutine scan_s2_bec_vals(l_base,s2_frac,ns,del_vals)
    real(dl), intent(in) :: l_base, s2_frac
    integer, intent(in) :: ns
    real(dl), dimension(1:ns), intent(out) :: del_vals
    
    real(dl), dimension(1:ns) :: sig_vals
    real(dl) :: sig_min, sig_max
    integer :: i_

    sig_max = (2._dl/3._dl)*log(l_base**2)*s2_frac
    sig_min = 0._dl
    sig_vals = (/ (sig_min + (i_-1)*(sig_max-sig_min)/dble(ns-1), i_=1,ns) /)
    del_vals = 0.5_dl*exp(-1.5*sig_vals)*l_base**2
  end subroutine scan_s2_bec_vals
  
  !>@brief
  !> Do full scans over parameter values for varying numbers of spatial dimensions.
  !>
  !>@param[in] modName (String) - name of model to use to name files
  !>@param[in] ord (Integer) - Order of interpolating polynomials
  !>@param[in] scanParams (real array) - Values of model parameter to scan over
  subroutine scan_dimensions(modName,ord,scanParams,get_grid,get_ic_params,prev_test)
    character(*), intent(in) :: modName
    integer, intent(in) :: ord
    real(dl), dimension(:,:), intent(in) :: scanParams
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
       !call execute_command_line('mv instantons.dat instantons-'//trim(adjustl(modName))//trim(adjustl(dimStr))//'-o'//trim(adjustl(ordStr))//'.dat')
    enddo
  end subroutine scan_dimensions

  ! Write this one
  subroutine scan_length_param(deltas,dim,ord,get_grid,get_ic_params,prev_test,grid_type,out)
    real(dl), dimension(:,:), intent(in) :: deltas
    real(dl), intent(in) :: dim
    integer, intent(in) :: ord
    procedure(get_grid_params) :: get_grid
    procedure(get_guess_params) :: get_ic_params
    procedure(use_previous) :: prev_test
    character(8), intent(in), optional :: grid_type
    logical, intent(in), optional :: out

    character(8) :: grid_
    logical :: out_
    integer :: i, nLen
    real(dl) :: lenFac

    out_ = .false.; if (present(out)) out_ = out
    grid_ = 'FULL_MID'; if (present(grid_type)) grid_ = grid_type
    
    nLen = 50
    do i=0,nLen
       lenFac = 1._dl  ! Adjust this here
       call scan_profiles(deltas,dim,ord,get_grid,get_ic_params,prev_test,grid_,out_) ! this won't work since I need to adjust the length parameter.  Redefine the get_grid function call
    enddo
  end subroutine scan_length_param
  
  !! TO DO: Put the interface for get_grid_params and get_guess_params somewhere besides the module
  !>@brief
  !> Solve for bounce profiles by scanning over the model parameters specified in the array deltas.
  !> Outputs the actions as a function of parameter, and optionally the full field profiles
  !>
  !>@params[in] deltas - [nPar,nProfile] Array of model parameters.  First index is simply the number of parameters needed to specify the model.  Second index scans over models to solve for
  !>@params[in] dim - Spatial dimension (possible float)
  !>@params[in] ord - Order of interpolation to use
  !>@params[in] get_grid - Name of subroutine that provides the grid parameters for the model
  !>@params[in] get_ic_params - Name of subrouting that provides parameters for initial profile guess
  !>@params[in] prev_test  - Boolean function that determines whether to use previous profile or not
  !>@params[in] grid_typ   - Optional -  Type of collocation grid to use.  Default is "FULL_MID"
  !>@params[in] out - Optional - (Boolean) Whether to output profiles or not.  Default is false
  subroutine scan_profiles(param_grid,dim,ord,get_grid,get_ic_params,prev_test,grid_type,out)
    real(dl), dimension(:,:), intent(in) :: param_grid
    real(dl), intent(in) :: dim
    integer, intent(in) :: ord
    procedure(get_grid_params) :: get_grid
    procedure(get_guess_params) :: get_ic_params
    procedure(use_previous) :: prev_test
    character(8), intent(in), optional :: grid_type
    logical, intent(in), optional :: out
    
    type(Instanton) :: inst
    type(Solver) :: solv
    real(dl) :: par_cur(1:nPar)
    logical :: outL
    integer :: i, fNum
    real(dl) :: len, w
    real(dl), dimension(1:4) :: params_ic
    integer :: p_i
    real(dl), dimension(0:ord) :: xNew, phi_prev
    character(8) :: grid_
    real(dl) :: phizero
    real(dl) :: pow
    
    grid_= 'FULL_MID'; if (present(grid_type)) grid_ = grid_type
    outL = .false.; if (present(out)) outL = out
    open(unit=newunit(fNum),file='actions.dat')

    call create_instanton(inst,ord,dim)
    call create_solver(solv,ord+1,100,0.1_dl)

    do i=1,size(param_grid(1,:))
       par_cur = param_grid(:,i)
       call set_model_params(par_cur,dim)

       ! Check if I can comment these things
       !call get_grid( (/dCur/),dim,len,w,pow )  ! Where is this pow thing defined?  Only for some models?
       call get_grid(get_model_params(),dim,len,w)
       call get_ic_params(get_model_params(),dim,params_ic,p_i)
       
       ! This is where I decide whether or not to use the previous solution
       ! TO DO: Fix the ugly call to get the new chebyshev grid
       ! Optional : Don't even bother to interpolate solution onto the new grid, then I don't have to do anything
       ! Problem : The interpolation seems broken for the Drummond potential
       if ( prev_test( get_model_params(),dim) .and. i.ne.1 ) then  ! Fix this to not have the 1 explicitly
          xNew = chebyshev_grid(inst%ord,len,w) ! Fix this to allow for different grids
          phi_prev = interpolate_instanton_(inst,xNew)
          call create_instanton_grid(inst,grid_,len,w)
          inst%phi(:) = phi_prev(:)
          print*,par_cur,"interpolate"
       else
          !call create_instanton_grid(inst,grid_,len,w,pow)
          call create_instanton_grid(inst,grid_,len,w)
          call profile_guess(inst,params_ic(1),params_ic(2),params_ic(3),params_ic(4),p_i)
       endif       
       
       call initialise_equations(inst%tForm,dim,inst%bc) ! Can I add this into the grid creation?

       ! Now solve
       call solve(solv,inst%phi)

       write(fNum,*) par_cur, compute_action(inst), interpolate_instanton_(inst,(/0._dl/)), inst%phi(0), 1./sqrt(abs(phi_prev(0))), sqrt(maxval(abs(phi_prev))), len, 1./w
       if (outL) call output_instanton(inst)
    enddo

    close(fNum); call delete_solver(solv); call destroy_instanton(inst)
  end subroutine scan_profiles

  subroutine extract_grid_params(inst,len,w,pow)
    type(Instanton), intent(in) :: inst
    real(dl), intent(out) :: len, w
    real(dl), intent(out), optional :: pow
    
    real(dl), dimension(0:inst%ord) :: dphi
    real(dl) :: m2_0, m2_max, r_max
    integer :: ind(1)
    
    dphi = matmul(inst%tForm%derivs(:,:,2),inst%phi)
    m2_0 = abs(dphi(0)) 
    m2_max = maxval(abs(dphi))
    ind = maxloc(abs(dphi))
    r_max = inst%tForm%xGrid(ind(1))

    if (sqrt(m2_0)*r_max < 1._dl) then
       len = 0.5_dl*twopi/sqrt(m2_0); w = 1._dl
    else
       len = r_max; w = 0.5_dl*twopi/(len*sqrt(m2_max))
    endif

    if (present(pow)) then
       ! Now fit power law to tail
    endif
  end subroutine extract_grid_params
    
  subroutine scan_resolutions(params,dim,ord,get_grid,get_ic_params,grid_type)
    real(dl), dimension(:), intent(in) :: params
    real(dl), intent(in) :: dim
    procedure(get_grid_params) :: get_grid
    procedure(get_guess_params) :: get_ic_params
    character(8), intent(in), optional :: grid_type

    type(Instanton) :: inst
    integer :: fNum_a, fnum_p
    integer :: ord
    character(8) :: grid_
    type(Solver) :: solv
    real(dl) :: len, w, pow, params_ic(1:4)
    integer :: p_i
    
    grid_ = 'FULL_MID'; if (present(grid_type)) grid_ = grid_type
    
    open(unit=newunit(fNum_a),file="actions.dat")
    open(unit=newunit(fNum_p),file="profiles.dat")

    ! Set up a grid interpolation here
    
    call set_model_params(params,dim)
    !call get_grid(params,dim,len,w,pow)
    !call get_grid(params,dim,len,w)  ! old version that doesn't allow for multiple params
    
    call get_grid(get_model_params(),dim,len,w)
    !    call get_ic_params(params,dim,params_ic,p_i)
    call get_ic_params(get_model_params(),dim,params_ic,p_i)
    do ord=10,200,4
       call create_instanton(inst,ord,dim); call create_solver(solv,ord+1,100,0.1_dl)
       !call create_instanton_grid(inst,grid_,len,w,pow)
       call create_instanton_grid(inst,grid_,len,w)
       call profile_guess(inst,params_ic(1),params_ic(2),params_ic(3),params_ic(4),p_i)
       call initialise_equations(inst%tForm,dim,inst%bc)
       call solve(solv,inst%phi)
       
       ! Interpolate onto a fixed grid then output so I can directly compare the instantons
       ! phi_int = matmul(interp,inst%phi)
       write(fNum_a,*) ord, compute_action(inst)
       call delete_solver(solv); call destroy_instanton(inst)
    enddo

    close(fNum_a); close(fNum_p)
  end subroutine scan_resolutions
    
  ! Super ugly, fix it
  subroutine linear_operator(L_lin,this,ord,dim,m2eff)
    type(Instanton), intent(in) :: this
    integer, intent(in) :: ord, dim
    real(dl), intent(out) :: L_lin(0:ord,0:ord)
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

  subroutine make_interpolation_matrix(tForm,xNew,interp)
    real(dl), dimension(1:), intent(in) :: xNew
    type(Chebyshev), intent(in) :: tForm
    real(dl), dimension(1:size(xNew),0:tForm%ord), intent(out) :: interp

    integer :: i,j,sz
    real(dl) :: xc, xBase(1:size(xNew))
    
    sz = size(xNew)
    ! First, convert the xNew coordinates to the base coordinate system
    ! Now evaluate Chebyshev polynomials at those xValues
    ! This gives my interpolation grid
    do i=1,sz
       xc = xBase(i)  ! This conversion needs to be done
       do j=0,tForm%ord
          interp(i,j) = 0._dl  ! evaluate Chebyshevs
       enddo
    enddo
  end subroutine make_interpolation_matrix

  logical function new_guess(params,dim) result(prev)
    real(dl), intent(in) :: params(1:nPar), dim
    prev = .false.
  end function new_guess
  
end program Instanton_Solver
