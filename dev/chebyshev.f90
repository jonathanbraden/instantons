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
  type Chebyshev
     integer :: nx, ord, nDeriv
     real(dl), allocatable :: xGrid(:), weights(:), norm(:)
     real(dl), allocatable :: fTrans(:,:)
     real(dl), allocatable :: invTrans(:,:)
     real(dl), allocatable :: derivs(:,:,:)
     real(dl), allocatable :: wFunc(:)
  end type Chebyshev

contains

!#ifdef INCLUDE_FUNCTIONS
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

  type(Chebyshev) function chebyshev_double_infinite(ord,len,evens,w) result(tForm)
    integer, intent(in) :: ord
    real(dl), intent(in) :: len
    logical, intent(in) :: evens
    real(dl), intent(in), optional :: w

    tForm = new_chebyshev_finite(ord,.true.)
    if (evens) call transform_to_evens(tForm)
    if (present(w)) call cluster_points(tForm,w,evens)
    call transform_double_infinite(tForm,len)
  end function chebyshev_double_infinite
!#endif

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

    call allocate_chebyshev(this, ord, nd)
    
    ! Get the collocation points
    if (end) then
       call chebyshev_lobatto_nodes(this%xGrid,this%weights,this%ord)
    else
       call chebyshev_gauss_nodes(this%xGrid,this%weights,this%ord)
    endif

    ! Fix this
    this%norm = pi
    this%norm(0) = 0.5_dl*this%norm(0)
    if (end) this%norm(ord) = 0.5_dl*this%norm(ord)

    ! Evaluate the basis functions on the collocation grid
    do i=0,ord
       x = this%xGrid(i)
!          call evaluate_chebyshev(ord,x,BVals,nd)
       call evaluate_chebyshev_trig(ord,x,BVals,nd)
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

    this%wFunc = sqrt(1._dl-this%xGrid(:)**2)  ! Move this somewhere else
  end subroutine create_chebyshev

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

    allocate( this%wFunc(0:ord) )
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
  end subroutine copy_chebyshev
  
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
    u=99
    open(unit=u,file='basis-functions.dat')
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
  subroutine compute_basis_norms(this)
    type(Chebyshev), intent(inout) :: this
    integer :: i
    do i=0,this%ord
       this%norm(i) = sum(this%weights(:)*this%invTrans(:,i)**2)
    enddo
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
  !> Free the memory stored in the input Chebyshev object
  !>
  !>@param[inout] this
  subroutine destroy_chebyshev(this)
    type(Chebyshev), intent(inout) :: this

    this%nx = -1; this%ord = -1; this%nDeriv = -1
    deallocate(this%xGrid,this%weights,this%norm)
    deallocate(this%fTrans,this%invTrans)
    deallocate(this%derivs)
  end subroutine destroy_chebyshev

  !>@brief
  !> Evaluate the given function via numerical quadrature
  !>
  !>@todo
  !>@arg For the chebyshev, implement this more efficiently since the weights are trivial
  real(dl) function quadrature(this,fVals) result(quad)
    type(Chebyshev), intent(in) :: this
    real(dl), dimension(:), intent(in) :: fVals

    ! Add a check that weights and f are the same size.  Add a preprocessor flag to include
#ifdef ERRORS
    if ( size(f) /= size(this%weights) ) then
       print*,"Error, size of function and order of polynomial expansion are incompatible"
    endif
#endif
    quad = sum(this%weights(:)*this%wFunc(:)*fVals(:))

#ifdef FAST_QUAD
    quad = sum(fVals(:))
    quad = quad*(pi/dble(order))
#endif
  end function quadrature

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
    integer, intent(in) :: ord
    real(dl), intent(in) :: x, nd
    real(dl), intent(out) :: T(0:ord,0:2)
    integer :: i

    T(0,0) = 1._dl; T(1,0) = x
    T(0,1:) = 0._dl; T(1,1) = 1._dl
    T(1,2:) = 0._dl
    T(2,0) = 2._dl*x*T(1,0) - T(0,0)
    do i=3,ord
       T(i,0) = 2._dl*x*T(i-1,0) - T(i-2,0)
       T(i,1) = 0._dl
    enddo
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
    integer, intent(in) :: ord, nd
    real(dl), intent(in) :: x
    real(dl), intent(out) :: T(0:ord,0:nd)
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
    integer :: n,o, i
    real(dl), allocatable :: B_tmp(:,:), spec(:)

    ! Add error check to make sure xNew and fNew are the same size
    ! Make sure that the indexing of xNew actually starts at 1 when it's defined implicitly
    ! as it is here

    fNew = 0._dl
    
    n = size(xNew); o = this%ord
    allocate(B_tmp(0:o,0:2)); allocate(spec(0:o))
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
    print*,""
    !w = 2._dl / dble(order+1)  ! wrong one
    w = 1._dl*pi / dble(order+1)
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
    w = 1._dl*pi / dble(order)
    w(0) = 0.5_dl*pi / dble(order); w(order) = 0.5_dl*pi / dble(order)
  end subroutine chebyshev_lobatto_nodes

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
    if (nd > 2) then
       print*,"Error, only derivatives up to order 2 implemented in transform_to_evens, defaulting to ",ndmax
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
  end subroutine transform_to_evens
  
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

    ord = this%ord
    this%xGrid(:) = len*this%xGrid(:) / sqrt(1._dl - this%xGrid(:)**2)

    allocate( dmap(0:this%ord,this%nDeriv) )
    dmap(:,1) = len**2 / (len**2 + this%xGrid(:)**2)**1.5
    dmap(:,2) = -3._dl*len**2*this%xGrid / (len**2 + this%xGrid(:)**2)**2.5
    call transform_derivatives(this,dmap)
    deallocate(dmap)
  end subroutine transform_double_infinite

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
    integer :: ord
    real(dl), dimension(:,:), allocatable :: dmap
    
    ord = this%ord
    this%xGrid(:) = len*( (1._dl+this%xGrid(:))/(1._dl-this%xGrid(:)) )
    allocate( dmap(0:this%ord,this%nDeriv) )
    dmap(:,1) = 2._dl*len / (len+this%xGrid)**2
    dmap(:,2) = -4._dl*len / (len+this%xGrid)**3
    
    call transform_derivatives(this,dmap)
    deallocate(dmap)
  end subroutine transform_semi_infinite

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
    real(dl), allocatable :: xshift(:), dmap(:,:)

    print*,"Point clustering not yet fully tested"
    
    allocate( xshift(0:this%ord) ); allocate( dmap(0:this%ord,1:this%nDeriv) )
    winv = 1._dl / w

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
    
!    if (dble_inf) then
!       xshift = pi*(this%xGrid-0.5_dl)
!       this%xGrid = atan(w*tan(xshift)) / pi + 0.5_dl
!       xshift = pi*(this%xGrid-0.5_dl)  ! We are now in y-coordinates
!       dmap(:,1) = winv*(1._dl+tan(xshift)**2) / (1._dl + winv**2*tan(xshift)**2 )
!       dmap(:,2) = 2._dl*pi*winv*(1._dl-winv**2)*tan(xshift)*(1._dl+tan(xshift)**2)/(1._dl+winv**2*tan(xshift**2) )
!    else
!       print*,"Error, clustering transform not yet written for half-infinite intervals"
!    endif

    call transform_derivatives(this,dmap)
    deallocate(xshift); deallocate(dmap)
  end subroutine cluster_points

  !>@brief
  !> Given the derivatives of the mapping parameter, make the appropriate transformations of the
  !> derivative matrices
  !>
  !> Given the derivatives of a coordinate mapping to a new set of coordinates for our spatial grid,
  !> this subroutine will make the appropriate transformations of the derivative operators so that
  !> that act in the new coordinate system.
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
  !>@param[in]  dmap  An array storing the derivatives up to the maximal derivative order in this
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
    do j=1,nd
       this%derivs(:,:,j) = matmul(this%derivs(:,:,j),this%invTrans(:,:))
    enddo
    ! Transform the derivatives to the new coordinates
    do j=0,ord
       this%derivs(:,j,2) = dmap(:,1)**2*this%derivs(:,j,2) + dmap(:,2)*this%derivs(:,j,1)
       this%derivs(:,j,1) = dmap(:,1)*this%derivs(:,j,1)
    enddo
    ! This line should be used for higher derivatives
    ! do i=2,nd; do j=0,ord
    !    this%derivs(:,j,i) = *this%derivs(:,j,i-1) + 
    ! enddo; enddo
    ! Transform the derivative matrices to act in real space
    do j=1,nd
       this%derivs(:,:,j) = matmul(this%derivs(:,:,j),this%fTrans(:,:))
    enddo
    this%weights = this%weights / dmap(:,1)
  end subroutine transform_derivatives

  !>@brief
  !> Given the coordinate mapping between two coordinate systems at the given sets of collocation points,
  !> compute the derivatives between the mappings using numerical differentiation.
  !> Due to instabilities of numerical differentiation and other numerical artifacts, it is preferable
  !> to have analytic expressions for these mappings if possible.
  !>
  !> Given an original set of collocation points \f$x_i\f$, and their mapped locations \f$y_i = y(x_i)\f$,
  !> use numerical differentiation to compute
  !>  \f[
  !>    \frac{d^{(m)}x}{dy^{(m)}} = 
  !>  \f]
  !> up to some specified order.
  !>
  !>
  !>@todo
  !>@arg Add a warning label to the use of this automated procedure
  !>@arg Write this procedure
  !>@arg If I want to write a general thing, it might be better to use a refined collocation grid to compute the
  !> numerical derivatives and then transform back
  !>@arg Test this
  subroutine differentiate_mapping(this, xNew, dMap, nd)
    type(Chebyshev), intent(in) :: this
    real(dl), dimension(:), intent(in) :: xNew
    real(dl), dimension(:,:), intent(out) :: dMap
    integer, intent(in) :: nd
    integer :: i, o

    o = this%ord
    do i=1,nd
#ifdef USEBLAS
       call DGEMV('N',(o+1),(o+1),1.d0,this%derivs(:,:,i), (o+1),xNew,1, 0.d0,dMap(:,i),1)
#else
       dMap(:,i) = matmul(this%derivs(:,:,i),xNew(:))
#endif
    enddo
  end subroutine differentiate_mapping

  !>@brief
  !> Perform a coordinate transform specified by the subroutine pointer
  !>
  !>@todo
  !>@arg Write the functionality for this subroutine
  subroutine transform_coordinates(this)
    type(Chebyshev), intent(inout) :: this
    this%xGrid = this%xGrid
  end subroutine transform_coordinates

end module Cheby
