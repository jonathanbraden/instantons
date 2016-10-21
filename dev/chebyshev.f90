! To do. Add a flag to do error checking for derivative orders, etc. that can be easily turned off when performance is important

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! MODULE:Chebyshev
!> @author Jonathan Braden
!>         University College London
!>
!> This provides the functionality for doing pseudospectral differentiation,
!> interpolation, (and Gaussian quadrature) using Chebyshev polynomial expansion
!>
!>@todo
!> Factor out all of the setup from the create_chebyshev call.  This will allow easier modularisation later
!> Add some inheritance properties so that we can easily implement other polynomial expansions
!> Add all of the stretching and compressing transformations.
!>
!>@todo
!> Add nice description of what this module is doing for the purposes of documentation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
       ! Fill in code to numerically invert the inverse transform
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
       !call DGEMM
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
    ! Make sure that the indexing of xNew actually starts at 1 when it's defined implicitly
    ! as it is here
    n = size(xNew)
    do i=1,n
!       fNew(i) = sum()  ! use evaluate_chebyshev
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

! Add in mapping subroutines here
  
end module Cheby
