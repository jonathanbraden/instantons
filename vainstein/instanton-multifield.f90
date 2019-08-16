module Instanton_Class
  use constants, only : dl
  use Utils, only : newunit
  use Cheby
  use Model  ! Try to remove this dependence
  use Nonlinear_Solver
  implicit none
  
  type Instanton_Multi
     type(Chebyshev) :: tForm
     integer :: nfld, ord
     real(dl), dimension(:,:), allocatable :: phi
     integer :: dim  ! Generalize this to noninteger
     logical :: exists = .false.
     integer :: unit = -1
  end type Instanton_Multi
  
contains

  subroutine create_instanton_multi(this,nf,ord,d)
    type(Instanton_multi), intent(out) :: this
    integer, intent(in) :: ord,d,nf
    
    this%dim = d; this%ord = ord; this%nfld = nf
    if (allocated(this%phi)) deallocate(this%phi)
    allocate( this%phi(0:ord,1:nf) )
    this%exists = .true.
    this%unit = 56  ! automate this to open an available unit
  end subroutine create_instanton_multi

  subroutine destroy_instanton_multi(this)
    type(Instanton_multi), intent(inout) :: this
    if (.not.this%exists) return
    deallocate(this%phi)
    this%exists = .false.
    this%unit = -1
  end subroutine destroy_instanton_multi

! Add interpolation subroutine here
  
  subroutine output_instanton_multi(this)
    type(Instanton_multi), intent(in) :: this

    real(dl), dimension(:,:), allocatable :: phi_spec, dphi, d2phi
    integer :: u; logical :: o
    integer :: ord, nf
    integer :: i
    
    ord = this%ord; nf = this%nfld
    ! File opening, fix it up and write here
    u = this%unit
    inquire(opened=o,unit=u) ! Move this and next line into initialization
    if (.not.o) open(unit=u,file='instanton_multi.dat')

    allocate(phi_spec(0:ord,1:nf),dphi(0:ord,1:nf),d2phi(0:ord,1:nf))
    
    phi_spec = matmul(this%tForm%fTrans,this%phi) 
    dphi = matmul(this%tForm%derivs(:,:,1),this%phi)
    d2phi = matmul(this%tForm%derivs(:,:,2),this%phi)
    do i=0,ord
       write(u,*) this%tForm%xGrid(i), this%phi(i,:), dphi(i,:), d2phi(i,:), phi_spec(i,:)
    enddo
    write(u,*)
    
    deallocate(phi_spec,dphi,d2phi)
  end subroutine output_instanton_multi

  !>@brief
  !> Compute the instanton solution for the symmetry breaking parameter delta
  !>
  !>@param[inout] this
  !>@param[in] delta
  !>@param[in] (optional) phi_init
  !>@param[in] (optional) out  - Boolean to write result to file or not
  !>@param[in] (optional) p_i  - Integer choice of initial analytic profile guess
  subroutine compute_profile_multi(this,delta,phi_init,out)
    type(Instanton_multi), intent(inout) :: this
    real(dl), dimension(:), intent(in) :: delta  ! Adjust this to number of model parameters
    real(dl), intent(in), optional :: phi_init(0:this%ord,1:this%nfld)
    logical, intent(in), optional :: out

    logical :: outLoc; integer :: order, dim, nv, nf
    type(Solver) :: solv

    ! Clean up all this extraneous crap
    real(dl) :: w, len      ! These seem extraneous

    real(dl), dimension(1:this%nfld*(this%ord+1)) :: data; integer :: i, n
    
    dim = this%dim; order = this%ord; nf = this%nfld
    outLoc = .false.; if (present(out)) outLoc = out
    nv = this%nfld*(this%ord+1)

    ! Need to rewrite this stuff
    len = 3._dl**0.5; w = 4._dl
    !    call create_grid_gauss(this%tForm,order,w,len) ! Replace this with the library call
    call create_grid_gauss(this%tForm,order,w,len)
    call create_solver(solv,nv,100,0.1_dl)
    call initialise_equations(this%tForm,delta,dim)
    
    ! Modularise this part
    if (present(phi_init)) then
       this%phi(0:order,1:nf) = phi_init(0:order,1:nf)
    else
       !call profile_guess(this,r0,meff,phif,phit,p_loc)
       this%phi = 0._dl
    endif

    ! As a hack, copy the fields into a 1D array to pass to solve
    n = order + 1
    do i=1,nf
       data((i-1)*n+1:i*n) = this%phi(:,i)
    enddo
    print*,"nv = ",nv," size data is ",size(data)
    call solve(solv,data)
    do i=1,nf
       this%phi(:,i) = data((i-1)*n+1:i*n)
    enddo
       
    if (outLoc) call output_instanton_multi(this)
  end subroutine compute_profile_multi

!!!! This functionality should be moved into the chebyshev code
!!! I'm pretty sure it's in there already, so just kill this and use the call in the library
  subroutine create_grid_gauss(tForm,ord,w,l)
    type(Chebyshev), intent(out) :: tForm
    integer, intent(in) :: ord
    real(dl), intent(in) :: w,l

    call create_chebyshev(tForm,ord,2,.false.,.true.)
    call transform_to_evens(tForm)
    call cluster_points(tForm,w,.true.)
    call transform_double_infinite(tForm,l)
    call IR_stretch(tForm,1.e-2)
  end subroutine create_grid_gauss

  subroutine create_grid_lobatto(tForm,ord,w,l)
    type(Chebyshev), intent(out) :: tForm
    integer, intent(in) :: ord
    real(dl), intent(in) :: w,l
    logical :: evens

    call create_chebyshev(tForm,ord,2,.true.,.true.)
    call cluster_points(tForm,w,.false.)
    call transform_semi_infinite(tForm,l)
    !call IR_stretch(tForm,1.e-2)
  end subroutine create_grid_lobatto
  
end module Instanton_Class
