!#def PRINT_WARNINGS T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! MODULE: Chebyshev
!> @author Jonathan Braden
!>         University College London
!>
!> @brief
!> This provides the functionality for doing pseudospectral differentiation,
!> interpolation, (and Gaussian quadrature) using Chebyshev polynomial expansion
!>
!> This module provides the functionality for pseudospectral differential, interpolation and quadrature
!> based on a Chebyshev polynomial expansion.
!> The function we approximating is assumed to be expanded as
!>  \f[ f(x) = \sum_{i=0}^N c_iB_i(x)  \f]
!> where \f$B_i(x)\f$ is the \f$i\f$th Chebyshev polynomial and \f$N\f$ is the order of the expansion.
!> The \f$B_i\f$'s satisfy the following recurrence relation
!>  \f[ 
!>     B_{i+1}(x) = 2xB_{i}(x)  - B_{i-1}(x) 
!>   \f]
!> with the initial conditions
!>  \f[ 
!>    B_0(x) = 1 \qquad B_1(x) = x \, .
!>  \f]
!> They have the following integral normalisation
!>  \f[  
!>    \int dx \frac{B_i(x)B_j(x)}{\sqrt{1-x^2}} = \begin{cases} 
!>                                                 \pi  &\quad : i=j=0 \\
!>                                                 \frac{\pi}{2}  &\quad : i=j\neq 0\\
!>                                                 0  &\quad : i \neq j
!>                                                \end{cases}
!>  \f]
!> and are suitable for quadruature integration with weight function \f$w(x) = \left(1-x^2\right)^{1/2}\f$
!>  \f[
!>    \int dx \frac{f(x)}{\sqrt{1-x^2}} \approx \sum_i w_if(x_i)
!>  \f]
!> with collocation points given by
!>  \f[
!>    x_i = \cos\left(\frac{\pi}{2(N+1)}(2i+1)\right) \qquad i=0,\dots,N
!>  \f]
!> for the Gaussian quadrature grid and
!>  \f[
!>    x_i = \cos\left(\frac{i\pi}{N}\right) \qquad i=0,\dots,N
!>  \f]
!> for the Gauss-Lobatto (i.e. extrema and endpoints) grid.
!> The corresponding weight functions are
!>  \f[
!>    w_i = 
!>  \f]
!> and
!>  \f[
!>    w_i = 
!>  \f]
!>  respectively.
!>
!> Implementation Details
!>
!> Example usage
!>
!>@todo
!> Factor out all of the setup from the create_chebyshev call.  This will allow easier modularisation later
!> Add some inheritance properties so that we can easily implement other polynomial expansions
!> Add all of the stretching and compressing transformations.
!>
!>@todo
!>@arg Add nice description of what this module is doing for the purposes of documentation
!>@arg When I'm doing the mappings, I currently have to first undo the transformation from spectral
!> to real space, apply the change in derivatives, then transform back.  Fix this.
!>@arg Flags for error checking of derivative orders for development purposes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#define USEBLAS  ! put this in a header

module Cheby
  use, intrinsic :: iso_c_binding
  use utils, only : newunit
#ifdef USEOMP
  use omp_lib
#endif
  use constants
  use Matrix
  implicit none

  !>@brief
  !> Stores information necessary for performing transformations between
  !> Chebyshev pseudospectral basis and real space.
  !>
  !>@todo
  !> norm - Stores normalisation of orthogonal polynomials (needed to generate transformation matrix)
  !> wVals - Weights to be used in quadrature integrations
  !> inv_coord_trans - Store the coordinate transform matrix for future use
  type Chebyshev
     integer :: nx, ord, nDeriv
     real(dl), allocatable :: xGrid(:), weights(:), norm(:)
     real(dl), allocatable :: fTrans(:,:)
     real(dl), allocatable :: invTrans(:,:)
     real(dl), allocatable :: derivs(:,:,:)
     real(dl), allocatable :: wFunc(:), w_quad(:)

     character(4) :: grid  ! Add this functionality
     character(4) :: transforms  ! Add this functionality
     real(dl), dimension(:), allocatable :: transform_params
     real(dl) :: len, scl
     real(dl) :: pow, ir_st
     logical :: evens=.false.
     logical :: do_cluster=.false., do_stretch=.false.
     
     logical :: init=.false.
  end type Chebyshev

contains

  ! These all use the old interface, which lacked the choice of collocation grid
  
  type(Chebyshev) function new_chebyshev_finite(ord,evens,w) result(tForm)
    integer, intent(in) :: ord
    logical, intent(in) :: evens
    real(dl), intent(in), optional :: w

    call create_chebyshev(tForm,ord,2,.false.,.false.)
    if (evens) call transform_to_evens(tForm)
    if (present(w)) call cluster_points(tForm,w,evens)
  end function new_chebyshev_finite

  type(Chebyshev) function new_chebyshev_semi_infinite(ord,len,w) result(tForm)
    integer, intent(in) :: ord
    real(dl), intent(in) :: len
    real(dl), intent(in), optional :: w

    tForm = new_chebyshev_finite(ord,.false.)
    if (present(w)) call cluster_points(tForm,w,.false.)
    call transform_semi_infinite(tForm,len)
  end function new_chebyshev_semi_infinite

  !>@brief
  !> Rational Chebyshev basis on double infinite interval.  Clustering applied after mapping to evens.
  type(Chebyshev) function new_chebyshev_double_infinite_(ord,len,evens,w) result(tForm)
    integer, intent(in) :: ord
    real(dl), intent(in) :: len
    logical, intent(in) :: evens
    real(dl), intent(in), optional :: w

    tForm = new_chebyshev_finite(ord,evens)
    if (present(w)) call cluster_points(tForm,w,evens)
    call transform_double_infinite(tForm,len)
  end function new_chebyshev_double_infinite_

  !>@brief
  !> Rational Chebyshev basis on double infinite interval.  Clustering applied before mapping to evens.
  type(Chebyshev) function new_chebyshev_double_infinite(ord,len,w) result(tForm)
    integer, intent(in) :: ord
    real(dl), intent(in) :: len
    real(dl), intent(in), optional :: w

    tForm = new_chebyshev_finite(ord,.false.)
    if (present(w)) call cluster_points(tForm,w,.false.)
    call transform_to_evens(tForm); tForm%evens = .true.
    call transform_double_infinite(tForm,len)
  end function new_chebyshev_double_infinite
  
  type(Chebyshev) function new_chebyshev_double_infinite_rational(ord,len,pow,evens,w) result(tForm)
    integer, intent(in) :: ord
    real(dl), intent(in) :: len, pow
    logical, intent(in) :: evens
    real(dl), intent(in), optional :: w

    tForm = new_chebyshev_finite(ord,evens)
    if (present(w)) call cluster_points(tForm,w,evens)
    call transform_double_infinite_rational(tForm,len,pow)
  end function new_chebyshev_double_infinite_rational
  
  !>@brief
  !> Create a Chebyshev transformation object to do polynomial approximation of order ord
  !>
  !>@param[out] this
  !>@param[in] ord
  !>@param[in] nd
  !>@param[in] end
  !>@param[in] num_inv
  subroutine create_chebyshev(this, ord, nd, end, num_inv)
    type(Chebyshev), intent(out) :: this
    integer, intent(in) :: ord, nd
    logical, intent(in) :: end, num_inv

    real(dl) :: x, BVals(0:ord,0:nd) ! figure out how to remove the ugly temporary BVals
    integer :: i

    this%len = -1._dl; this%scl = 1._dl
    if (this%init) call destroy_chebyshev(this)
    call allocate_chebyshev(this, ord, nd)

    ! Get the collocation points
    if (end) then
       call chebyshev_lobatto_nodes(this%xGrid,this%weights,this%ord)
    else
       call chebyshev_gauss_nodes(this%xGrid,this%weights,this%ord)
    endif

    ! Fix this to include the Radau grid
    ! Is this being overwritten when I call compute_basis norms?
    ! What if I comment this out?  Does anything break?
    this%norm = pi
    this%norm(0) = 0.5_dl*this%norm(0)
    if (end) this%norm(ord) = 0.5_dl*this%norm(ord) ! Why the hell is this one different?
    ! Even more annoyance with the above.  It looks like the factors of 1/2 are actually wrong ...
    
    ! Evaluate the basis functions on the collocation grid
    do i=0,ord
       x = this%xGrid(i)
!          call evaluate_chebyshev(ord,x,BVals,nd)
       call evaluate_chebyshev_trig(ord,x,BVals,nd)  ! Will break for endpoints
       this%invTrans(i,:) = BVals(:,0)
       this%derivs(i,:,1:nd) = BVals(:,1:nd)
    enddo
    ! Find a more elegant way to remove the dependence on end.
    Call compute_basis_norms(this)
    call make_mmt(this,num_inv,end)
 
    !>@todo
    !> Rather than doing a numerical matrix multiplication here, use explicit formulas in the above loop
    !> Add DGEMM support
    do i=1,nd
#ifdef USEBLAS
       call DGEMM( , )
#else
       this%derivs(:,:,i) = matmul(this%derivs(:,:,i),this%fTrans(:,:))
#endif
    enddo

    this%wFunc(:) = sqrt(1._dl-this%xGrid(:)**2)  ! Move this somewhere else
    this%w_quad(:) = this%weights(:)*this%wFunc(:)
    this%init = .true.
  end subroutine create_chebyshev

  subroutine create_chebyshev_full(this,ord,nd,grid,num_inv,recurse)
    type(Chebyshev), intent(out) :: this
    integer, intent(in) :: ord, nd
    character(4), intent(in) :: grid
    logical, intent(in) :: num_inv
    logical, intent(in), optional :: recurse
    
    real(dl) :: x, BVals(0:ord,0:nd) ! figure out how to remove BVals
    integer :: npt
    integer :: i

    logical :: use_recursion, use_num_inv

    use_num_inv = num_inv
    use_recursion = .false.; if (present(recurse)) use_recursion = recurse
    
    this%len = -1._dl; this%scl = 1._dl  ! See if I need this (better to store parameters for the transform)
    if (this%init) call destroy_chebyshev(this)
    call allocate_chebyshev(this, ord, nd)

    select case (grid)
    case ('MIDS')
       call chebyshev_gauss_nodes(this%xGrid,this%weights,this%ord)
       this%grid = 'MIDS'
    case ('ENDS')
       call chebyshev_lobatto_nodes(this%xGrid,this%weights,this%ord)
       this%grid = 'ENDS'
       use_recursion = .true.
    case ('LEFT')
       call chebyshev_radau_L_nodes(this%xGrid,this%weights,this%ord)
       this%grid = 'LEFT'
       use_recursion = .true.
    case ('RGHT')
       call chebyshev_radau_R_nodes(this%xGrid,this%weights,this%ord)
       this%grid = 'RGHT'
       use_recursion = .true.
    case default
       print*,"Invalid Collocation Grid"
       this%grid = 'NONE'
    end select
    
    ! Evaluate the basis functions on the collocation grid
    do i=0,ord
       x = this%xGrid(i)
       if (use_recursion) then
          call evaluate_chebyshev(ord,x,BVals,nd)      ! works for endpoints?
       else
          call evaluate_chebyshev_trig(ord,x,BVals,nd) ! breaks for endpoints?
       endif
       this%invTrans(i,:) = BVals(:,0)
       this%derivs(i,:,1:nd) = BVals(:,1:nd)
    enddo
    
    ! Find a more elegant way to remove the dependence on end.
    call compute_basis_norms(this,num_sum=.true.)
    call make_mmt(this,num_inv,.false.) ! This shouldn't be hardcoded in evens parameter
    
    !>@todo
    !> Rather than doing a numerical matrix multiplication here, use explicit formulas in the above loop
    !> Add DGEMM support
    npt = ord + 1
    do i=1,nd
#ifdef USEBLAS
!       call DGEMM('N','N',npt,npt,npt,1._dl,this%derivs(:,:,i),npt,this%fTrans(:,:),npt,)
#else
       this%derivs(:,:,i) = matmul(this%derivs(:,:,i),this%fTrans(:,:))
#endif
    enddo

    this%wFunc(:) = sqrt(1._dl-this%xGrid(:)**2)  ! Move this somewhere else
    this%w_quad(:) = this%weights(:)*this%wFunc(:)
    this%init = .true.
  end subroutine create_chebyshev_full
  
  !>@brief
  !> Free the memory stored in the input Chebyshev object
  !>
  !>@param[inout] this
  subroutine destroy_chebyshev(this)
    type(Chebyshev), intent(inout) :: this

    if (this%init) then
       this%nx = -1; this%ord = -1; this%nDeriv = -1
       deallocate(this%xGrid,this%weights,this%norm)
       deallocate(this%fTrans,this%invTrans)
       deallocate(this%derivs)
       deallocate(this%wFunc,this%w_quad)
    endif
    this%init = .false.
  end subroutine destroy_chebyshev
  
  !>@brief
  !> Allocate space in our transformation to store grid points, weights, derivative matrices, etc.
  !>
  !>@param[out] this
  !>@param[in]  ord   The maximal polynomial order
  !>@param[in]  nd    Number of derivative matrices to compute and store
  subroutine allocate_chebyshev(this,ord,nd)
    type(Chebyshev), intent(out) :: this
    integer, intent(in) :: ord, nd

    this%ord = ord; this%nDeriv = nd; this%nx = ord+1
    allocate( this%xGrid(0:ord), this%weights(0:ord), this%norm(0:ord) )
    allocate( this%fTrans(0:ord,0:ord) )
    allocate( this%invTrans(0:ord,0:ord) )
    allocate( this%derivs(0:ord,0:ord,1:nd) )

    allocate( this%wFunc(0:ord), this%w_quad(0:ord) )
  end subroutine allocate_chebyshev
  
  !>@brief
  !> Copy the contents of an existing transform
  !>
  !>@param[in]  old
  !>@param[out] new_t
  !>
  !>@todo
  !>@arg Test this
  subroutine copy_chebyshev(old,new_t)
    type(Chebyshev), intent(in) :: old
    type(Chebyshev), intent(inout) :: new_t

    call allocate_chebyshev(new_t,old%ord,old%nDeriv)
    new_t%xGrid = old%xGrid; new_t%weights = old%weights; new_t%norm = old%norm
    new_t%fTrans = old%fTrans; new_t%invTrans = old%invTrans
    new_t%derivs = old%derivs
    new_t%wFunc = old%wFunc; new_t%w_quad = old%w_quad
    new_t%init = .true.
  end subroutine copy_chebyshev

  !>@brief
  !> Write details about the quadrature grid, including the weights, continuum weights,
  !> and weights needed to do a direct sum
  !>
  !>@todo
  !>@arg Write this thing properly
  subroutine output_quadrature_grid(this)
    type(Chebyshev), intent(in) :: this
    integer :: u, i
    open(unit=newunit(u),file='chebyshev-grid.dat')
    do i=0,this%ord
       write(u,*) this%xGrid(i), this%weights(i), this%w_quad(i)
    enddo
    write(u,*)
    close(u)
  end subroutine output_quadrature_grid
  
  !>@brief
  !> Write the basis functions used by the transform to file
  !>
  !>@todo
  !>@arg Add some options to output a specified number of derivatives, or only a finite
  !>     number of the basis functions
  !>@arg Test this
  subroutine output_basis_functions(this)
    type(Chebyshev), intent(in) :: this
    integer :: u,i,j
    open(unit=newunit(u),file='basis-functions.dat')
    do i=0,this%ord
       do j=0,this%ord
          write(99,*) this%xGrid(j), this%invTrans(j,i)
       enddo
       write(99,*)
    enddo
    close(u)
  end subroutine output_basis_functions

  !>@brief
  !> Compute the norms of the basis polynomials based on the numerical quadrature
  subroutine compute_basis_norms(this,num_sum)
    type(Chebyshev), intent(inout) :: this
    logical, intent(in), optional :: num_sum

    integer :: i
    logical :: num_

    num_ = .true.; if (present(num_sum)) num_ = num_sum

    if (num_) then
       do i=0,this%ord
          this%norm(i) = sum(this%weights(:)*this%invTrans(:,i)**2)
       enddo
    else  ! This looks questionable, like it's off by a factor of 0.5
       this%norm(0) = twopi; this%norm(1:) = pi
       ! Lobatto isn't exact when computing the square of the final basis function
       if (this%grid=='ENDS') this%norm(this%ord) = twopi 
    endif
  end subroutine compute_basis_norms
  
  !>@brief
  !> Helper routine to evaluate the basis functions and derivatives on the collocation grid
  !>
  !>@param[in,out] this The Chebyshev transform to initialise
  !>@param[in] ord      Polynomial order of the interpolation
  !>@param[in] nd       Number of derivatives to compute
  !>@param[in] recur    Whether evaluate basis functions using recursion of analytic formulae
  subroutine eval_basis_functions(this,ord,nd,recur)
    type(Chebyshev), intent(inout) :: this
    integer, intent(in) :: nd, ord
    logical, intent(in) :: recur
    integer :: i
    real(dl) :: x, BVals(0:ord,0:nd)
    
    do i=0,ord
       x = this%xGrid(i)
       if (recur) then
          call evaluate_chebyshev(ord,x,BVals,nd)
       else
          call evaluate_chebyshev_trig(ord,x,BVals,nd)
       endif
       this%invTrans(i,:) = BVals(:,0)
       this%derivs(i,:,1:nd) = BVals(:,1:nd)
       !this%fTrans(:,i) = this%weights(i)*BVals(:,0) !/ this%norm(:)
    enddo
  end subroutine eval_basis_functions

  !>@brief
  !> Make Matrix Multiplication Transform matrix from real space to spectral space  
  !>
  !> For the case of a transformation from real space to spectral space
  !> based on Gaussian quadrature, the \f$k\f$th spectral coefficient \f$c_k\f$ is given by
  !> \f[ 
  !>   c_k = \sum_i w_i f(x_i) B_k(x_i)
  !> \f]
  !> where $B_k(x_i)$ are the basis functions evaluated at collocation point \f$x_i\f$, \f$w_i\f$ are the collocation weights, and \f$\mathcal{N}_k\f$ is the normalisation of the orthogonal functions 
  !> \f[ \int dx w(x)B_k(x)^2 = \mathcal{N}_k  \, . \f]
  !> The transformation from real to spectral space can be implemented using the matrix multiplication \f$c_k = M_{ki}f_i\f$ with
  !> \f[
  !>   M_{ki} = 
  !> \f]
  !>
  !>@param[in,out] this
  !>@param[in]     num_inv   Whether to compute the MMT using a numerical inversion of the forward transform or not
  !>@param[in]     end       Whether or not we are using the Gauss-Lobatto grid.  Requires modifed normalisation of the highest mode if we are (when doing analytical inverse)
  !>
  !>@todo
  !>@arg  Implement the numerical inverse
  subroutine make_mmt(this,num_inv,end)
    type(Chebyshev), intent(inout) :: this
    logical, intent(in) :: num_inv, end
    integer :: i, ord
    
    ord = this%ord
    if (num_inv) then
       call invert_matrix(this%invTrans,this%fTrans,this%nx)
    else
       do i=0,ord
          this%fTrans(:,i) = this%weights(:)*this%invTrans(:,i)/this%norm(i)
       enddo
       this%fTrans = transpose(this%fTrans)
    endif    
  end subroutine make_mmt
  
  !>@brief
  !> Evaluate the given function via numerical quadrature
  !>  \f[
  !>    I = \int f(x) dx \approx sum_i f(x_i)w_i
  !>  \f]
  !> Returns the value of \f$I\f$ as defined above.
  !>
  !>@todo
  !>@arg For the chebyshev, implement this more efficiently since the weights are trivial
  real(dl) function quadrature(this,fVals) result(quad)
    type(Chebyshev), intent(in) :: this
    real(dl), dimension(0:this%ord), intent(in) :: fVals

    !quad = sum(this%weights(:)*this%wFunc(:)*fVals(:))
    quad = sum(this%w_quad(:)*fVals(:))
#ifdef FAST_QUAD
    quad = sum(fVals(:))
    quad = quad*(pi/dble(order))
#endif
  end function quadrature

  !>@brief
  !> Evaluate the integral of the given function times the continuum weight function
  !> appropriate for the remapped collocation grid using numerical quadrature
  !>  \f[
  !>    I \equiv \int w(x)f(x)dx
  !>  \f]
  !> Returns the value of \f$I\f$ as defined above
  !>
  !>@param[in] this  They Chebyshev object storing the collocation grid information
  !>@param[in] fVals The function values evaluated at the collocation points
  real(dl) function quadrature_w_weight(this,fVals) result(quad)
    type(Chebyshev), intent(in) :: this
    real(dl), dimension(0:this%ord), intent(in) :: fVals
    quad = sum( this%weights(:)*fVals(:) )
  end function quadrature_w_weight
    
  !>@brief
  !> Evaluate the given function (passed as a function pointer) using numerical quadrature
  !>
  !>@todo
  !>@arg Implement this using nicer vectorisation and atomic functions instead of the loop
!  real(dl) function quadrature_function(this,f) result(quad)
!    type(Chebyshev), intent(in) :: this
!    procedure :: f  ! fill this in properly.
!    integer :: i
!    real(dl) :: ftmp

!    quad = 0._dl
!    do i=0,this%ord
!       ftmp = f(this%xGrid(i))
!       quad = quad + ftmp*this%weights(i)
!    enddo
!  end function quadrature_function
  
  !>@brief
  !> Create the nth derivative operator by repeated application of the first derivative
  !>
  !> WARNING: the first derivative must already have been initialised for this to work
  !>@todo
  !>@arg Add DGEMM functionality
  !>@arg Fix this so it is actually correct
  subroutine generate_nth_derivative_matrix(this,nd,d_dx,dn_dx)
    type(Chebyshev), intent(in) :: this
    integer, intent(in) :: nd
    real(dl), dimension(:,:), intent(in) :: d_dx
    real(dl), dimension(:,:), intent(out) :: dn_dx
    integer :: i
    dn_dx = this%fTrans
    do i=1,nd
#ifdef USEBLAS

#else
       dn_dx = matmul(d_dx,dn_dx)
#endif
    enddo
  end subroutine generate_nth_derivative_matrix
  
! Some useful debugging subroutines
! 1. Check that forward and inverse transforms are inverses
! 2. Check how different multiple applications of derivative matrix and direct calculation of higher derivatives is

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Non Object Specific Subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>@brief
  !> Evaluate the Chebyshev polynomials 
  !>
  !> Evaluate the Chebyshev polynomials and their derivatives using the recurrence relation
  !> \f[ 
  !>   T_{i}(x) = 2xT_{i-1}(x) - T_{i-2}
  !> \f]
  !> along with the boundary conditions
  !> \f[ T_0(x) = 1 \qquad T_1(x) = x \f] \, .
  !> This generates the corresponding Chebyshev functions given by
  !> \f[ T_i =  \f]
  !> which have the normalisation
  !>@todo: Fill in normalisation
  !>
  !>@param[in]  ord  The maximal order of the polynomials required
  !>@param[in]  x    The value of \f$x\f$
  !>@param[out] T    Array to store the value of the Chebyshev polynomials and its derivatives
  !>@param[in]  nd   Number of derivatives to compute
  !>
  !>@todo
  !>@arg Modify to allow and arbitrary number of derivatives
  subroutine evaluate_chebyshev(ord,x,T,nd)
    integer, intent(in) :: ord, nd
    real(dl), intent(in) :: x
    real(dl), intent(out) :: T(0:ord,0:2)
    integer :: i

    ! Add some checks for how many derivatives are needed in here
    ! Extend to allow for higher than second order derivatives
    T(0,0) = 1._dl
    do i=1,nd
       T(0,i) = 0._dl
    enddo
    T(1,0) = x
    T(1,1) = 1._dl
    do i=2,nd
       T(1,i) = 0._dl
    enddo
    do i=2,ord
       T(i,0) = 2._dl*x*T(i-1,0) - T(i-2,0)
       T(i,1) = 2._dl*x*T(i-1,1) - T(i-2,1) + 2._dl*T(i-1,0)
       T(i,2) = 2._dl*x*T(i-1,2) - T(i-2,2) + 4._dl*T(i-1,1)
    enddo
  end subroutine evaluate_chebyshev

  !@brief
  !> Evaluate Chebyshev's using the recurrence relation for the individual polynomials
  !> \f[
  !>   B_n(x) = 2xB_{n-1}(x) - B_{n-2}(x)
  !> \f]
  !> and the derivative recurrence for all derivatives
  !>  \f[
  !>    2 B_n(x) = \frac{1}{n+1}B'_{n+1}(x) - \frac{1}{n-1}B'_{n-1}(x)
  !>  \f]
  !>
  !>@todo
  !> Finish writing this thing for arbitrary choices of nd
  subroutine evaluate_chebyshev_recur(ord,x,T,nd)
    integer, intent(in) :: ord, nd
    real(dl), intent(in) :: x
    real(dl), intent(out) :: T(0:ord,0:nd)
    integer :: i,j

    print*,"Warning evaluate_chebyshev_recur not fully tested"

    if (nd < 0) then
       print*,"Error, need to at least evaluate the basis functions"
       stop
    endif
    T(0:1,:) = 0._dl
    T(0,0) = 1._dl; T(1,0) = x
    if (nd >= 1) T(1,1) = 1._dl

    do i=2,ord
       T(i,0) = 2._dl*x*T(i-1,0) - T(i-2,0)
    enddo

    ! Check if I can remove this if statement and just not have the loop execute
    if (nd >=1) then
       do j=1,nd
          do i=2,ord
             T(i,j) = 2._dl*x*T(i-1,j) - T(i-2,j) + 2._dl*j*T(i-1,j-1)
          enddo
       enddo
    endif
  end subroutine evaluate_chebyshev_recur
    
  !>@brief
  !> Evaluate Chebyshev polynomials using trigonometric representation
  !>
  !> The Chebyshev polynomials at grid point x are evaluated using
  !> \f[
  !>   B_n(x) = \cos\left(n\cos^{-1}(x)\right)
  !>  \f]
  !> and the corresponding derivatives
  !> \f[
  !>   B_n'(x) =
  !> \f]
  !>@todo
  !> Finish writing this, allow for Lobatto endpoints
  !> Add appropriate control flow to allow nd to be specified
  subroutine evaluate_chebyshev_trig(ord,x,T,nd)
    integer, intent(in) :: ord,nd
    real(dl), intent(in) :: x
    real(dl), intent(out) :: T(0:ord,0:2)
    
    real(dl) :: icn, dn
    integer, dimension(0:ord) :: iVals
    integer :: i

    icn = dacos(x); dn = 1._dl / sqrt(1._dl-x**2)
    iVals = (/ (i, i=0,ord) /)
    T(:,0) = cos(iVals*icn)
    T(:,1) = dble(iVals) * sin(dble(iVals)*icn) * dn
    T(:,2) = -cos(iVals*icn) * dble(iVals**2) * dn**2 + sin(iVals*icn) * iVals * dn**3 * x
  end subroutine evaluate_chebyshev_trig

  !>@brief
  !> Interpolate the given function off of the collocation grid
  subroutine interpolate(this,fVals,xNew,fNew)
    type(Chebyshev), intent(in) :: this
    real(dl), intent(in) :: fVals(:), xNew(1:)
    real(dl), intent(out) :: fNew(:)
    
    integer :: n,o,i
    real(dl) :: spec(0:this%ord), B_tmp(0:this%ord,0:2)

    ! Add error check to make sure xNew and fNew are the same size
    ! Make sure that the indexing of xNew actually starts at 1 when it's defined implicitly
    ! as it is here

    fNew = 0._dl
    n = size(xNew); o = this%ord
#ifdef USEBLAS
    call DGEMV
#else
    spec = matmul(this%fTrans,fVals)
#endif
    do i=1,n
       call evaluate_chebyshev(o,xNew(i),B_tmp,0)
       fNew(i) = sum(B_tmp(0:o,0)*spec(0:o))  ! use evaluate_chebyshev
    enddo
  end subroutine interpolate
  
! Add Gauss-Cheby and Gauss-Lobatto nodes
! These are all fucked up from the point of view of actually doing quadrature.
  ! There are missing factors due to the normalisation of the Chebyshevs

  !>@brief
  !> Compute the Gauss-Chebyshev abscissa and weights for a given maximal polynomial order
  !> These are located at the zeroes of \f$T_{N+1}\f$ where \f$N\f$ is the interpolation order.
  !>
  !>@param[out] x      The quadrature collocation points
  !>@param[out] w      The quadrature weights
  !>@param[in]  order  Interpolating polynomial order
  subroutine chebyshev_gauss_nodes(x,w,order)
    integer, intent(in) :: order
    real(dl), dimension(0:order), intent(out) :: x,w
    double precision :: dkcol
    integer :: i

    dkcol = 0.5_dl* pi / dble(order+1)
    do i=0,order
       x(i) = -dcos( dble(2*i+1)*dkcol )
    enddo
    w = pi / dble(order+1)
  end subroutine chebyshev_gauss_nodes

  !>@brief
  !> Compute the Gauss-Chebyshev-Lobatto (endpoints and extrema) absicca and weights.
  !> These are given by the solutions of \f$(1-x_i^2)T'_{N}(x_i) = 0\f$
  !>
  !>@param[out] x      The quadrature collocation points
  !>@param[out] w      The quadrature weights
  !>@param[in]  order  Interpolating polynomial order
  subroutine chebyshev_lobatto_nodes(x,w,order)
    integer, intent(in) :: order
    real(dl), dimension(0:order), intent(out) :: x,w
    double precision :: dkcol
    integer :: i

    dkcol = pi  / dble(order)
    do i=0,order
       x(i) = -dcos(dble(i)*dkcol)
    enddo
    w = pi / dble(order)
    w(0) = 0.5_dl*pi / dble(order); w(order) = 0.5_dl*pi / dble(order)
  end subroutine chebyshev_lobatto_nodes

!!!! Need to test both of these
!!!! I'm pretty sure the weights in the textbook on spectral methods for turbulence have a typo, which then got transferred here, so need to check these
  
  !>@brief
  !> Compute the Gauss-Chebyshev-Raud abscissa and weights, with the left endpoint included in the grid.
  subroutine chebyshev_radau_R_nodes(x,w,order)
    integer, intent(in) :: order
    real(dl), dimension(0:order), intent(out) :: x,w
    double precision :: dkcol
    integer :: i

    dkcol = pi/(2._dl*dble(order)+1._dl)
    do i=0,order
       x(i) = -dcos((2._dl*dble(i)+1._dl)*dkcol)
    enddo
    w(0:order-1) = twopi/(2._dl*dble(order)+1._dl)
    w(order) = pi/(2._dl*dble(order)+1._dl)
  end subroutine chebyshev_radau_R_nodes

    !>@brief
  !> Compute the Gauss-Chebyshev-Radau abscissa and weights, with the right endpoint included in the grid.
  subroutine chebyshev_radau_L_nodes(x,w,order)
    integer, intent(in) :: order
    real(dl), dimension(0:order), intent(out) :: x,w
    double precision :: dkcol
    integer :: i

    dkcol = twopi/(2._dl*dble(order)+1._dl)
    do i=0,order
       x(i) = -dcos(dble(i)*dkcol)
    enddo
    w(1:order) = twopi/(2._dl*dble(order)+1._dl)
    w(0) = pi/(2._dl*dble(order)+1._dl)
  end subroutine chebyshev_radau_L_nodes

!!! Finish testing

  
  !>@brief
  !> Make transform matrix from real space to spectral space  
  !>
  !> For the case of a transformation from real space to spectral space
  !> based on Gaussian quadrature, the \f$k\f$th spectral coefficient \f$c_k\f$ is given by
  !> \f[ 
  !>   c_k = \sum_i w_i f(x_i) B_k(x_i)
  !> \f]
  !> where $B_k(x_i)$ are the basis functions evaluated at collocation point \f$x_i\f$, \f$w_i\f$ are the collocation weights, and \f$\mathcal{N}_k\f$ is the normalisation of the orthogonal functions 
  !> \f[ \int dx w(x)B_k(x)^2 = \mathcal{N}_k  \, .\f]
  !> The transformation from real to spectral space can be implemented using the matrix multiplication \f$c_k = M_{ki}f_i\f$ with
  !> \f[
  !>   M_{ki} = 
  !> \f]
  subroutine make_forward_transform(ord, trans)
    integer, intent(in) :: ord
    real(C_DOUBLE), dimension(0:ord,0:ord), intent(out) :: trans
    integer :: i
    
    do i=0,ord
       trans(:,i) = 0._dl !b(:,i)*wgrid(:)  ! Need cheby values and weights
    enddo
    trans(:,0) = 0.5_dl*trans(:,0)
    trans = transpose(trans)
  end subroutine make_forward_transform

  !>@brief
  !> Make transform matrix from spectral space to real space
  !>
  !> To transform from spectral space back to real space with the function evaluated at the collocation points, we simply use
  !> \f[ f(x_i) = \sum_k c_k B_k(x_i)  \f]
  !> which can be implemented as a matrix multiplication
  !> \f[ f_i = B_{ik} c_k \f]
  !> with 
  !> \f[ B_{ik} = B_k(x_i) \f]
  subroutine make_backward_transform(ord,trans)
    integer, intent(in) :: ord
    real(C_DOUBLE), dimension(0:ord,0:ord), intent(out) :: trans
    integer :: i
    real(dl), dimension(0:ord,0:1) :: B
    real(dl) :: x

    ! This loop is all broken right now.
    do i=0,ord
       call evaluate_chebyshev(ord,x,B,1)  ! Insert line to get Chebyshevs evaluated at x
       trans(:,i) = B(:,0)
    enddo
  end subroutine make_backward_transform

  !>@brief
  !> Apply the coordinate transformation
  !> \f$ y(x) = \sqrt{\frac{x}{2}+\frac{1}{2}}\f$
  !> to transform to a basis of only the even Chebyshev polynomials
  !>  \f[
  !>    B_{2i}(y) = 
  !>  \f]
  !>
  !>@todo
  !> Add higher than second derivatives in here
  subroutine transform_to_evens(this)
    type(Chebyshev), intent(inout) :: this
    real(dl), dimension(:,:), allocatable :: dmap
    integer :: i, nd, ord
    integer, parameter :: ndmax = 2

    this%xGrid = sqrt(0.5_dl*this%xGrid + 0.5_dl)
    nd = this%nDeriv; ord = this%ord
    if (nd > ndmax) then
       print*,"Error, only derivatives up to order ",ndmax," implemented in transform_to_evens, defaulting to ",ndmax
       nd = ndmax
    endif
    allocate( dmap(0:this%ord,nd) )
    dmap(:,1) = 4._dl*this%xGrid
    dmap(:,2) = 4._dl
    do i=3,nd
       dmap(:,i) = 0._dl
    enddo
    call transform_derivatives(this,dmap)
    deallocate( dmap )
    this%evens = .true.
  end subroutine transform_to_evens

  subroutine invert_evens(xVals)
    real(dl), dimension(:), intent(inout) :: xVals
    xVals = 2._dl*xVals**2 - 1._dl
  end subroutine invert_evens
  
  !>@brief
  !> Transform to rational Chebyshev functions on the double infinite interval using the transform
  !>  \f[
  !>    y(x) = L\frac{x}{\sqrt{1-x^2}}
  !>  \f]
  subroutine transform_double_infinite(this, len)
    type(Chebyshev), intent(inout) :: this
    real(dl), intent(in) :: len
    real(dl), dimension(:,:), allocatable :: dmap
    integer :: ord

    this%len = len
    ord = this%ord
    this%xGrid(:) = len*this%xGrid(:) / sqrt(1._dl - this%xGrid(:)**2)

    allocate( dmap(0:this%ord,this%nDeriv) )
    dmap(:,1) = len**2 / (len**2 + this%xGrid(:)**2)**1.5
    dmap(:,2) = -3._dl*len**2*this%xGrid / (len**2 + this%xGrid(:)**2)**2.5

    ! Check if I need to deal with endpoints differently
    call transform_derivatives(this,dmap)
    deallocate(dmap)
  end subroutine transform_double_infinite

  subroutine invert_double_infinite(xVals,len)
    real(dl), dimension(:), intent(inout) :: xVals
    real(dl), intent(in) :: len
    xVals = xVals / sqrt(xVals**2+len**2)
  end subroutine invert_double_infinite
  
  subroutine transform_tanh_infinite(this,len)
    type(Chebyshev), intent(inout) :: this
    real(dl), intent(in) :: len
    real(dl), dimension(0:this%ord,1:2) :: dmap

    ! Add a storage of parameter
    dmap(:,1) = (1._dl-this%xGrid**2)/len
    dmap(:,2) = -2._dl*this%xGrid*(1._dl-this%xGrid**2)/len**2
    this%xGrid = atanh(this%xGrid)/len
    call transform_derivatives(this,dmap)
  end subroutine transform_tanh_infinite

  subroutine invert_tanh_infinite(xVals,len)
    real(dl), dimension(:), intent(inout) :: xVals
    real(dl), intent(in) :: len
    xVals = tanh(len*xVals)/len
  end subroutine invert_tanh_infinite
    
  !>@brief
  !> Transform to rational Chebyshev functions on the singly infinite interval using the coordinate transform
  !>  \f[
  !>    y(x) = L\frac{1+x}{1-x}
  !>  \f]
  !> or
  !>  \f[
  !>    x(y) = \frac{r-L}{r+L}
  !>  \f]
  !>@param[in,out] this
  !>@param[in] len  The length parameter \f$L\f$ in the transformation
  subroutine transform_semi_infinite(this, len)
    type(Chebyshev), intent(inout) :: this
    real(dl), intent(in) :: len
    real(dl), dimension(0:this%ord,this%nDeriv) :: dmap

    this%len = len
    this%xGrid(:) = len*( (1._dl+this%xGrid(:))/(1._dl-this%xGrid(:)) )
    dmap(:,1) = 2._dl*len / (len+this%xGrid)**2
    dmap(:,2) = -4._dl*len / (len+this%xGrid)**3

    ! See if I need to treat the endpoints differently
    ! And if there's a more stable numerical approach to the above
    
    call transform_derivatives(this,dmap)
  end subroutine transform_semi_infinite

  subroutine invert_semi_infinite(xVals,len)
    real(dl), dimension(:), intent(inout) :: xVals
    real(dl), intent(in) :: len
    xVals = (xVals - len)/(xVals + len)
  end subroutine invert_semi_infinite
  
  subroutine transform_double_infinite_rational(this,l,p)
    type(Chebyshev), intent(inout) :: this
    real(dl), intent(in) :: l,p
    real(dl), dimension(0:this%ord,this%nDeriv) :: dmap

    print*,"Warning, double infinite rational transform not tested"
    this%len = l; this%pow = p
    
    ! First compute y^(n)(x)
    !dmap(:,1) = l*(1._dl+(p-1.)*this%xGrid(:)**2)/(1._dl-this%xGrid(:)**2)**(0.5*p+1)
    !dmap(:,2) = l*this%xGrid(:)*p*(3.+(p-1.)*this%xGrid(:)**2)/(1._dl-this%xGrid(:)**2)**(0.5*p+2)
    ! Now covert to x^(n)(y)
    ! TO DO: replace two lines below with this call
    ! call transform_derivatives_forward_to_backward(dmap(:,1:2))
    !dmap(:,1) = 1._dl/dmap(:,1)
    !dmap(:,2) = -dmap(:,2)*dmap(:,1)**3  ! Is this correct
    
    ! Combine everything analytically, which avoids some potential roundoff errors
    dmap(:,2) = -p*this%xGrid(:)*(1._dl-this%xGrid(:)**2)**(p+1.)/l**2*(3._dl+(p-1._dl)*this%xGrid(:)**2)/(1._dl+(p-1._dl)*this%xGrid(:)**2)**3
    dmap(:,1) = (1._dl-this%xGrid(:)**2)**(0.5_dl*p+1.) / (1._dl+(p-1._dl)*this%xGrid(:)**2) / l
    this%xGrid(:) = l*this%xGrid(:) / (1._dl-this%xGrid(:)**2)**(0.5*p)
    call transform_derivatives(this,dmap)
  end subroutine transform_double_infinite_rational

  ! Test this one
  subroutine transform_asinh_stretch(this,w)
    type(Chebyshev), intent(inout) :: this
    real(dl), intent(in) :: w
    real(dl), dimension(0:this%ord,1:2) :: dmap
    
    this%xGrid = sinh(w*this%xGrid)/w
    dmap(:,1) = cosh(w*this%xGrid)
    dmap(:,2) = w*sinh(w*this%xGrid)
    call transform_derivatives(this,dmap)
  end subroutine transform_asinh_stretch

  subroutine invert_asinh_stretch(this,w)
    type(Chebyshev), intent(inout) :: this
    real(dl), intent(in) :: w

    this%xGrid = asinh(w*this%xGrid)/w
  end subroutine invert_asinh_stretch
  
  ! This one seems wrong.  Check it because I think I messed it up
  subroutine transform_exponential_stretch(this,m)
    type(Chebyshev), intent(inout) :: this
    real(dl), intent(in) :: m
    real(dl), dimension(0:this%ord,1:2) :: dmap

    ! Add storage of parameter
    dmap(:,1) = exp(-m*this%xGrid)
    dmap(:,2) = -exp(-2._dl*m*this%xGrid)
    this%xGrid = exp(m*this%xGrid)/m-1._dl
    call transform_derivatives(this,dmap)
  end subroutine transform_exponential_stretch

  subroutine invert_exponential_stretch()
  end subroutine invert_exponential_stretch
  
  subroutine transform_ir_stretch(this,s)
    type(Chebyshev), intent(inout) :: this
    real(dl), intent(in) :: s
    real(dl), dimension(0:this%ord,this%nDeriv) :: dmap
    !this%s = s  ! needed if I want to resample

    print*,"Warning, this is unfinished and untested"
    dmap(:,1) = cosh(s*this%xGrid(:))
    dmap(:,2) = s*sinh(s*this%xGrid(:))
    dmap(:,1) = 1./dmap(:,1)
    dmap(:,2) = -dmap(:,2)*dmap(:,1)**3
    this%xGrid(:) = sinh(s*this%xGrid(:))/s

    call transform_derivatives(this,dmap)
  end subroutine transform_ir_stretch

  subroutine invert_ir_stretch(xVals,s)
    real(dl), dimension(:), intent(inout) :: xVals
    real(dl), intent(in) :: s
    xVals = asinh(s*xVals)/s
  end subroutine invert_ir_stretch
  
  !>@brief
  !> Cluster collocation points near a boundary layer using the coordinate mapping
  !>  \f[
  !>    y(x) = \frac{s}{\pi}\tan^{-1}\left(w\tan\left(\frac{\pi}{s}\left[x-x_0\right] \right)\right) + x_0
  !>  \f]
  !> where \f$x_0\f$ is \f$(b-a)/2 + (a+b)/2\f$ and \f$s=b-a\f$
  !> and the collocation points are on the interval \f$ [a,b] \f$. 
  !> The corresponding first and second derivatives are
  !>  \f[
  !>    \frac{dx}{dy}(y) = w^{-1}\frac{1+t(y)^2}{1+w^{-2}t(y)^2} \qquad t(y) \equiv \tan\left(\frac{\pi}{s}(y-x_0)\right)
  !>  \f]
  !>  \f[
  !>    \frac{d^2x}{dy^2} = 2w^{-1}\frac{\pi}{s}\left(1-w^{-2}\right)t(y)\left(\frac{1+t(y)^2}{1+w^{-2}t(y)^2}\right)^2
  !>  \f]
  !>
  !>@param[in,out] this
  !>@param[in]  w      Width parameter \f$w\f$ describing the transformation.
  !>@param[in]  evens  Whether or not the grid is the infinite or semi-infinite interval
  !>
  !>@todo
  !>@arg Implement the choice of collocation point to cluster around
  subroutine cluster_points(this, w, evens)
    type(Chebyshev), intent(inout) :: this
    real(dl), intent(in) :: w
    logical, intent(in) :: evens

    real(dl) :: winv, x0, s, sinv
    real(dl) :: xshift(0:this%ord), dmap(0:this%ord,1:this%nDeriv)

#ifdef PRINT_WARNINGS
    print*,"Point clustering not yet fully tested"
#endif
    this%scl = w; winv = 1._dl / w

    if (evens) then
       s = pi/1._dl; sinv = 1._dl/pi; x0=0.5_dl
    else
       s = pi/2._dl; sinv = 2._dl/pi; x0=0._dl
    endif

    ! This can be automated for any finite grid
    xshift = s*(this%xGrid-x0)
    this%xGrid = sinv*atan(w*tan(xshift)) + x0
    xshift = s*(this%xGrid-x0)
    dmap(:,1) = winv*(1._dl+tan(xshift)**2) / (1._dl + winv**2*tan(xshift)**2)
    dmap(:,2) = 2._dl*s*winv*(1._dl-winv**2)*tan(xshift)*(1._dl+tan(xshift)**2)/(1._dl+winv**2*tan(xshift)**2)**2    

    ! Fix these for the endpoints
    ! dmap(0,1) = w; dmap(this%ord,1) = w
    ! dmap(0,2) = 0._dl; dmap(this%ord,2) = 0._dl

    call transform_derivatives(this,dmap)
  end subroutine cluster_points

  ! Check this one
  subroutine invert_cluster_points(xVals,w,evens)
    real(dl), dimension(:), intent(inout) :: xVals
    real(dl), intent(in) :: w
    logical, intent(in) :: evens

    real(dl) :: winv, s, sinv, x0

    winv = 1._dl/w
    if (evens) then; s = pi; sinv = 1._dl/pi; x0 = 0.5_dl
    else; s = 0.5_dl*pi; sinv = 2._dl/pi; x0 = 0._dl
    endif
    xVals = sinv*atan(winv*tan(s*(xVals-x0))) + x0
  end subroutine invert_cluster_points
  
  !>@brief
  !> Given the derivatives of a coordinate mapping to a new set of coordinates for our spatial grid,
  !> \f[
  !>    y = y(x)
  !> \f]
  !> this subroutine will make the appropriate transformations of the derivative operators so that they act in the new coordinate system.
  !> The appropriate transformations follow from the chain rule, with the first few given by
  !>  \f{align}{
  !>    \frac{df}{dy} &= \frac{dx}{dy}\frac{df}{dx} \\
  !>    \frac{d^2f}{dy^2} &= \frac{d^2x}{dy^2}\frac{df}{dx} + \left(\right)\frac{d^2f}{dx^2} \\
  !>      &\dots
  !>    \frac{d^{(n)}f}{dy^2} &=
  !>  \f}
  !> and the general result given by recursion as
  !>  \f[
  !>    Fill this in
  !>  \f]
  !> Since differentiation via Chebyshev pseudospectral expansion is mildly ill-conditioned,
  !> only low order derivatives should be computed in this manner.
  !>
  !> Furthermore, in order to facilitate numerical quadrature calculations, the volume transformation encoded in the Jacobian is folded into the width calculation, so the new weights are
  !> \f[
  !>   \tilde{w}_i = \frac{dx_{new}}{dx_{old}}(x_{old,i})w_i
  !>  \f]
  !>
  !>@param[in,out] this
  !>@param[in]  dmap  An array storing the derivatives up to the maximal derivative order in this.  These are dx^n/dy^2, where y is the new coordinate and x the old coordinate.
  !>
  !>@todo
  !>@arg Extend to higher than second order derivatives
  !>@arg Implement the weighting recalculation.
  subroutine transform_derivatives(this,dmap)
    type(Chebyshev), intent(inout) :: this
    real(dl), dimension(:,:), intent(in) :: dmap
    integer :: ord, nd
    integer :: j
    integer, parameter :: ndmax = 2

    ord = this%ord; nd = this%nDeriv
    if (nd > ndmax) then
       print*,"Error, cannot transform coordinates for derivatives larger than ",ndmax
       print*,"Reverting to ",ndmax," transformed derivatives"
       nd = ndmax
    endif

    ! Convert back to spectral space to get a matrix of basis function derivatives
    ! Temporally commented.  Se if I need it
!    do j=1,nd
!       this%derivs(:,:,j) = matmul(this%derivs(:,:,j),this%invTrans(:,:))
!    enddo
    ! Transform the derivatives to the new coordinates
    do j=0,ord
       this%derivs(:,j,2) = dmap(:,1)**2*this%derivs(:,j,2) + dmap(:,2)*this%derivs(:,j,1)
       this%derivs(:,j,1) = dmap(:,1)*this%derivs(:,j,1)
    enddo
    ! Transform the derivative matrices to act in real space
    ! Temporarily commented again for testing
    !do j=1,nd
    !   this%derivs(:,:,j) = matmul(this%derivs(:,:,j),this%fTrans(:,:))
    !enddo
    this%weights = this%weights / dmap(:,1)
    this%w_quad = this%w_quad / dmap(:,1)
  end subroutine transform_derivatives

  subroutine transform_derivatives_forward(this,dmap)
    type(Chebyshev), intent(inout) :: this
    real(dl), dimension(:,:), intent(in) :: dmap
    integer :: ord, nd
    integer :: j
    integer, parameter :: ndmax = 2
    
    ord = this%ord; nd = this%nDeriv
    if (nd > ndmax) then
       print*,"Error, derivative transforms for order larger than ",ndmax," not yet supported."
       print*,"Defaulting to ",ndmax," transformed derivatives"
       nd = ndmax
    endif

    do j=0,ord
       this%derivs(:,j,2) = this%derivs(:,j,2)/dmap(:,1)**2 - dmap(:,2)/dmap(:,1)**3*this%derivs(:,j,1)
       this%derivs(:,j,1) = this%derivs(:,j,1)/dmap(:,1)
    enddo
    this%weights = this%weights / dmap(:,1)
    this%w_quad = this%w_quad / dmap(:,1)
  end subroutine transform_derivatives_forward
       
  !>@brief
  !> Given the derivatives of the forward transform \f$y(x)\f$
  !> i.e.
  !>  \f[
  !>    \frac{dy}{dx} \qquad \frac{d^2y}{dx^2} \qquad \dots
  !>  \f]
  !> convert these to derivatives of the backward transform
  !> i.e.
  !> \f[
  !>   \frac{dx}{dy} \qquad \frac{d^2x}{dy^2} \qquad \dots
  !> \f]
  subroutine transform_derivatives_forward_to_backward(this,dmap)
    type(Chebyshev), intent(inout) :: this
    real(dl), dimension(0:this%ord,1:2), intent(inout) :: dmap

    print*,"Derivative transforming not yet tested fully"
    dmap(:,1) = 1._dl/dmap(:,1)
    dmap(:,2) = -dmap(:,2)*dmap(:,1)**3
  end subroutine transform_derivatives_forward_to_backward
  
end module Cheby
