! To do. Add a flag to do error checking for derivative orders, etc. that can be easily turned off when performance is important
!
! Implement integration of a function of the entire interval by inverting the derivative matrix
! Then test that this thing actually works

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! MODULE:Chebyshev
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
!> and are suitable for quadruature integration with weight function \f$w(x) = \left(1-x^2\right^{1/2}\f$
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
!> Add nice description of what this module is doing for the purposes of documentation
!> When I'm doing the mappings, I currently have to first undo the transformation from spectral
!> to real space, apply the change in derivatives, then transform back.  Fix this.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define USEBLAS  ! put this in a header

module Cheby
  use, intrinsic :: iso_c_binding
#ifdef USEOMP
  use omp_lib
#endif
  use constants
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
  end type Chebyshev

contains
  
  !>@brief
  !> Create a Chebyshev transformation object to do polynomial approximation of order ord
  subroutine create_chebyshev(this, ord, nd, end, num_inv)
    type(Chebyshev), intent(out) :: this
    integer, intent(in) :: ord, nd
    logical, intent(in) :: end, num_inv

    real(dl) :: x, BVals(0:ord,0:nd) ! figure out how to remove the ugly temporary BVals
    integer :: i

    this%ord = ord; this%nDeriv = nd; this%nx = ord+1
    allocate( this%xGrid(0:ord), this%weights(0:ord), this%norm(0:ord) )
    allocate( this%fTrans(0:ord,0:ord) )
    allocate( this%invTrans(0:ord,0:ord) )
    allocate( this%derivs(0:ord,0:ord,1:nd) )

    ! Get the collocation points
    if (end) then
       call chebyshev_lobatto_nodes(this%xGrid,this%weights,this%ord)
    else
       call chebyshev_gauss_nodes(this%xGrid,this%weights,this%ord)
    endif
    
    this%norm = pi
    this%norm(0) = 0.5_dl*this%norm(0)
    
    ! Evaluate the basis functions on the collocation grid
    do i=0,ord
       x = this%xGrid(i)
!          call evaluate_chebyshev(ord,x,BVals,nd)
       call evaluate_chebyshev_trig(ord,x,BVals,nd)
       this%invTrans(i,:) = BVals(:,0)
       this%derivs(i,:,1:nd) = BVals(:,1:nd)
       !this%fTrans(:,i) = this%weights(i)*BVals(:,0) !/ this%norm(:)
    enddo

    if (num_inv) then
       ! Fill in code to numerically invert to get the inverse transform
    else
       do i=0,ord
          this%fTrans(:,i) = this%weights(:)*this%invTrans(:,i)
       enddo
       this%fTrans(:,0) = 0.5_dl*this%fTrans(:,0)
       this%fTrans = transpose(this%fTrans)
       if (end) this%fTrans(ord,:) = 0.5_dl*this%fTrans(ord,:)
       ! For the Lobatto grid, this needs to be modified due to the incorrect normalisation of the highest mode on this grid.  Check this carefully with Legendre etc.
    endif

    !>@todo
    !> Rather than doing a numerical matrix multiplication here, use explicit formulas in the above loop
    !> Add DGEMM support
    do i=1,nd
       !call DGEMM()
       this%derivs(:,:,i) = matmul(this%derivs(:,:,i),this%fTrans(:,:))
    enddo
  end subroutine create_chebyshev

  !>@todo
  !> Helper routine to evaluate the basis functions and derivatives on the collocation grid
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

  !>@todo
  !> Create the Matrix Multiplication Transform from real to spectral space
  subroutine make_mmt(this,num_inv)
    type(Chebyshev), intent(inout) :: this
    logical, intent(in) :: num_inv
    integer :: i, ord
    
    if (num_inv) then
       ! Fill in code to numerically invert the inverse transform
    else
       do i=0,ord
          this%fTrans(:,i) = this%weights(:)*this%invTrans(:,i)
       enddo
       this%fTrans(:,0) = 0.5_dl*this%fTrans(:,0)
       this%fTrans = transpose(this%fTrans)
    endif    
  end subroutine make_mmt
  
  !>@brief
  !> Free the memory stored in the input Chebyshev object
  subroutine destroy_chebyshev(this)
    type(Chebyshev), intent(inout) :: this

    this%nx = -1; this%ord = -1
    deallocate(this%xGrid,this%weights,this%norm)
    deallocate(this%fTrans,this%invTrans)
  end subroutine destroy_chebyshev

  !>@brief
  !> Evaluate the given function via numerical quadrature
  real(dl) function quadrature(this,f) result(quad)
    type(Chebyshev), intent(in) :: this
    real(dl), dimension(:), intent(in) :: f

    ! Add a check that weights and f are the same size.  Add a preprocessor flag to include
#if ERRORS
    if ( size(f) /= size(this%weights) ) then
       print*,"Error, size of function and order of polynomial expansion are incompatible"
    endif
#endif
    quad = sum(this%weights(:)*f(:))
  end function quadrature

  !>@brief
  !> Create the nth derivative operator by repeated application of the first derivative
  !>
  !> WARNING: the first derivative must already have been initialised for this to work
  subroutine generate_nth_derivative_matrix(this,nd,d_dx,dn_dx)
    type(Chebyshev), intent(in) :: this
    integer, intent(in) :: nd
    real(dl), dimension(:,:), intent(in) :: d_dx
    real(dl), dimension(:,:), intent(out) :: dn_dx
    integer :: i
    dn_dx = this%fTrans
    do i=1,nd
       dn_dx = matmul(d_dx,dn_dx)
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
  !> \f[ T_0(x) = 1. \qquad T_1(x) = x \f] \, .
  !> This generates the corresponding Chebyshev functions given by
  !> \f[ T_i =  \f]
  !> which have the normalisation
  !>@todo: Fill in normalisation
  !>
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
  !> \f[ B_n(x) = 2xB_{n-1}(x) - B_{n-2}(x)\f]
  !> and the derivative recurrence for all derivatives
  !> \f[ 2B_n(x) = \frac{1}{n+1}B'_{n+1}(x) - \frac{1}{n-1}B'_{n-1}(x)\f]
  !>
  !>@todo
  !> Finish writing this thing
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
  !> \f[ B_n(x) = \cos\left(n\cos^{-1}(x)\right) \f]
  !> and the corresponding derivatives
  !>
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
  subroutine interpolate(fVals,xNew,fNew)
    real(dl), intent(in) :: fvals(:), xNew(:)
    real(dl), intent(out) :: fNew(:)
    integer :: n, i

    ! Add error check to make sure xNew and fNew are the same size
    n = size(xNew)
    do i=1,n
!       fNew(i) =   ! use evaluate_chebyshev
    enddo
  end subroutine interpolate

! Add Gauss-Cheby and Gauss-Lobatto nodes
! These are all fucked up from the point of view of actually doing quadrature.
! There are missing factors due to the normalisation of the Chebyshevs
  subroutine chebyshev_gauss_nodes(x,w,order)
    real(dl), dimension(0:order), intent(out) :: x,w
    integer, intent(in) :: order
    real*8 :: dkcol
    integer :: i

    dkcol = 0.5_dl* pi / dble(order+1)
    do i=0,order
       x(i) = -dcos( dble(2*i+1)*dkcol )
    enddo
    print*,""
    w = 2._dl / dble(order+1)
  end subroutine chebyshev_gauss_nodes

  subroutine chebyshev_lobatto_nodes(x,w,order)
    real(dl), dimension(0:order), intent(out) :: x,w
    integer, intent(in) :: order
    real*8 :: dkcol
    integer :: i

    dkcol = pi  / dble(order)
    do i=0,order
       x(i) = -dcos(dble(i)*dkcol)
!       print*,"Lobatto node is ",x(i),cos(dble(order)*acos(x(i)))
    enddo
    w = 2._dl / dble(order)
    w(0) = 1._dl / dble(order); w(order) = 1._dl / dble(order)
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
    integer :: i,j
    
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
    integer :: i,j
    real(dl), dimension(0:ord,0:1) :: B
    real(dl) :: x

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
    integer :: j, ord

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
  subroutine transform_semi_infinite(this, len)
    type(Chebyshev), intent(inout) :: this
    real(dl), intent(in) :: len
  end subroutine transform_semi_infinite


  !>@brief
  !> Given the derivatives of the mapping parameter, make the appropriate transformations of the
  !> derivative matrices
  !>
  !> Given the derivatives of a coordinate mapping to a new set of coordinates for our spatial grid,
  !> this subroutine will make the appropriate transformations of the derivative operators so that
  !> that act in the new coordinate system.
  subroutine transform_derivatives(this,dmap)
    type(Chebyshev), intent(inout) :: this
    real(dl), dimension(:,:), intent(in) :: dmap
    integer :: ord, nd
    integer :: j
    integer, parameter :: ndmax = 2

    ord = this%ord; nd = this%nDeriv
    if (nd > 2) then
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
    ! Transform the derivative matrices to act in real space
    do j=1,nd
       this%derivs(:,:,j) = matmul(this%derivs(:,:,j),this%fTrans(:,:))
    enddo
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
  !> Add a warning label to the use of this automated procedure
  !> Write this procedure
  !> If I want to write a general thing, it might be better to use a refined collocation grid to compute the
  !> numerical derivatives and then transform back
  subroutine differentiate_mapping(this, xOld,xNew, dMap, nd)
    type(Chebyshev), intent(in) :: this
    real(dl), dimension(:), intent(in) :: xOld, xNew
    real(dl), dimension(:,:), intent(out) :: dMap
    integer, intent(in) :: nd
  end subroutine differentiate_mapping

  !>@brief
  !> Perform a coordinate transform specified by the subroutine pointer
  !>@todo Write this stupid thing
  subroutine transform_coordinates(this)
    type(Chebyshev), intent(inout) :: this
  end subroutine transform_coordinates

end module Cheby
