module Instanton_Class
  use constants, only : dl
  use Cheby
  use Model  ! Try to remove this dependence
  use Nonlinear_Solver
  implicit none
  
  type Instanton
     type(Chebyshev) :: tForm
     real(dl), dimension(:), allocatable :: phi
     integer :: ord  
     integer :: dim  ! Generalize for this to be noninteger
     real(dl) :: r0, meff, phif, phit  ! Do I really need these?
     logical :: exists = .false.
  end type Instanton
  
contains

  ! Expand to noninteger dimensions
  subroutine create_instanton(this,ord,d)
    type(Instanton), intent(out) :: this
    integer, intent(in) :: ord,d
    this%dim = d; this%ord = ord
    if (allocated(this%phi)) deallocate(this%phi) ! Remove this to only allocate if size has changed
    allocate(this%phi(0:ord))
    this%exists = .true.
  end subroutine create_instanton

  subroutine destroy_instanton(this)
    type(Instanton), intent(inout) :: this
    if (.not.this%exists) return
    deallocate(this%phi)
    this%exists = .false.
  end subroutine destroy_instanton
  
  function interpolate_instanton_(this,r_new) result(f_int)
    type(Instanton), intent(in) :: this
    real(dl), dimension(:), intent(in) :: r_new
    real(dl), dimension(1:size(r_new)) :: f_int

    real(dl) :: L,w
    real(dl) :: xvals(1:size(r_new)), spec(1:this%tForm%ord+1), bVals(1:this%tForm%ord+1,0:2)
    integer :: i; real(dl) :: winv
    
    w = this%tForm%scl; L = this%tForm%len
    winv = 1._dl/w
    xvals = r_new / sqrt(r_new**2+L**2)
    xvals = atan(winv*tan(pi*(xvals-0.5_dl)))/pi + 0.5_dl
    xvals = 2._dl*xvals**2 - 1._dl
#ifdef USEBLAS
    spec = 
#else
    spec = matmul(this%tForm%fTrans,this%phi)
#endif
    do i=1,size(r_new)
       call evaluate_chebyshev(this%tForm%ord,xvals(i),bVals,2)
       f_int(i) = sum(spec*bVals(:,0))
    enddo
  end function interpolate_instanton_
  
  !TO DO: phif is extracted from the solution, not from input.
  subroutine output_instanton(this)
    type(Instanton), intent(in) :: this
    
    real(dl), dimension(:), allocatable :: phi_spec, dphi, d2phi, phi_spec_y
    logical :: o
    integer :: i, sz
    real(dl) :: phif

    integer, parameter :: u = 56 ! Fix this to automatically select an open unit
    
    inquire(opened=o,unit=u)
    if (.not.o) open(unit=u,file='instanton_.dat')

    sz = size(this%phi); phif = this%phi(sz-1)
    allocate(phi_spec(0:sz-1),dphi(0:sz-1),d2phi(0:sz-1),phi_spec_y(0:sz-1))
    
    phi_spec = matmul(this%tForm%fTrans,this%phi) 
    dphi = matmul(this%tForm%derivs(:,:,1),this%phi)
    d2phi = matmul(this%tForm%derivs(:,:,2),this%phi)
    phi_spec_y = matmul(this%tForm%fTrans,d2phi)
    do i=0,sz-1
       write(u,*) this%tForm%xGrid(i), this%tForm%weights(i), this%phi(i), dphi(i), d2phi(i), phi_spec(i), phi_spec_y(i), external_source(this%tForm%xGrid(i))
    enddo
    write(u,*)
    
    deallocate(phi_spec,dphi,d2phi)
  end subroutine output_instanton

  !>@brief
  !> Compute the instanton solution for the symmetry breaking parameter delta
  !>
  !>@param[inout] this
  !>@param[in] delta
  !>@param[in] (optional) phi_init
  !>@param[in] (optional) out  - Boolean to write result to file or not
  !>@param[in] (optional) p_i  - Integer choice of initial analytic profile guess
  subroutine compute_profile_(this,delta,phi_init,out)
    type(Instanton), intent(inout) :: this
    real(dl), intent(in) :: delta
    real(dl), intent(in), optional :: phi_init(0:this%ord)
    logical, intent(in), optional :: out

    logical :: outLoc; integer :: order, dim, n, p_loc
    type(Solver) :: solv

    ! Fix up this initialisation of the grid, etc.
    real(dl) :: len, w
    
    dim = this%dim; order = this%ord
    outLoc = .false.; if (present(out)) outLoc = out
    n = order+1

    len = 3._dl**0.5  ! r0 is 1
    w = 0.8_dl   ! no stretching as a test
    call create_grid_(this%tForm,order,w,len)
    call create_solver(solv,n,100,0.1_dl)
    call initialise_equations(this%tForm,delta,dim)
    
    ! Modularise this part (Fix this for Vainstein
    if (present(phi_init)) then
       this%phi(0:order) = phi_init(0:order)
    else
       !      call profile_guess(this,r0,meff,phif,phit,p_loc)
       this%phi = 0._dl ! Fix this
    endif
    
    call solve(solv,this%phi)
    if (outLoc) call output_instanton(this)
  end subroutine compute_profile_

!!!! This functionality should be moved into the chebyshev code
!!! I'm pretty sure it's in there already, so just kill this and use the call in the library
  subroutine create_grid_(tForm,ord,w,l)
    type(Chebyshev), intent(out) :: tForm
    integer, intent(in) :: ord
    real(dl), intent(in) :: w,l

    call create_chebyshev(tForm,ord,2,.false.,.true.)
    call transform_to_evens(tForm)
    call cluster_points(tForm,w,.true.)
    call transform_double_infinite(tForm,l)
  end subroutine create_grid_

  !>@brief
  !> Initialise our initial profile guess based on the given radius and width.
  subroutine profile_guess(this)
    type(Instanton), intent(inout) :: this
  end subroutine profile_guess
  
end module Instanton_Class
