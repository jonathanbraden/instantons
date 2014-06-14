!
! Code to obtain instanton profiles (and in the future actions, etc.)
! in scalar field theories
!
! Compile with
!  gfortran -fdefault-real-8 -fdefault-double-8 -O3 -xhost -o instanton instanton.f90 -llapack -lm
!
!#define CUBIC 1
#define LINEAR 1


#define XLOOP do i=0,nx
#define ENDXLOOP enddo
#define XRANGE 0:nx
#define BRANGE 0:nbasis

program instanton
  implicit none

  real, parameter :: twopi =  6.283185307179586476925867665590
!!!!!!!!!!!
! Parameters for the field model
!!!!!!!!!!!
  real, parameter :: g2=1., epsilon=-0.01
  real :: delta
  integer, parameter :: numdelta = 1
  real, parameter, dimension(numdelta) :: deltavals = (/ 0.1 /)
  integer, parameter :: nfield=2
  real, dimension(nfield), parameter :: vac_guess = (/-1., 0./)
  real, parameter :: phiout = 1.
  real :: rinit, len
  real, dimension(nfield) :: falsevac

  integer :: m

!!!!!!!!!!
! Parameters controlling order of interpolation, etc.
!!!!!!!!!!
  integer, parameter :: nmax = 511 ! order of highest interpolating polynomial
  integer, parameter :: nx = nmax,  nbasis=nx
  integer, parameter :: fsize = nfield*(nx+1)  ! needed for DGESV

!!!!!!!!!
! Storage for derivative operators, grids, transform matrices, etc.
! Leave these alone unless you wish to extend functionality
!!!!!!!!! 
  real, dimension(XRANGE) :: xgrid, wgrid
  real, dimension(XRANGE,BRANGE) :: b, bp, bpp, mmt
  real, dimension(XRANGE,nfield) :: freal, del
  real, dimension(BRANGE,nfield) :: fspec
  real, dimension(XRANGE) :: map_p, map_pp

! Storage for various pieces of the differential operators.  These will depend on the precise equation
  real, dimension(XRANGE, BRANGE) :: l0, l1
  real, dimension(XRANGE) :: l2
  real, dimension(XRANGE,nfield) :: error  ! violation of equations of motion
  real :: errtot, deltot
  real, parameter :: errtol = 5.e-10, deltol = 5.e-10
  integer, parameter :: maxit = 1000

  integer :: i,j,l
  real :: lambda = 1.e-3

! Find the false vacuum
  delta = deltavals(1)
  falsevac = vac_guess
  call get_vacuum(falsevac)
  print*,"# Derivative at proposed false vacuum is ", vprime(falsevac,1), vprime(falsevac,2)

! Use overall 1.5 for the double well w/ phi^3, a 1 otherwise
#ifdef LINEAR
  rinit =  2.1**0.5/deltavals(1)
#endif
#ifdef CUBIC
  rinit = 1.5*2.**0.5/deltavals(1)
#endif
  len = rinit  ! do I need this anywhere else?

  call get_collocation_points(xgrid, wgrid, nmax)
  call get_basis_matrix(b, bp, bpp)  ! evaluate basis functions and derivatives at collocation points
  call transform_coords(xgrid, rinit, .true.)
  call make_transform_matrix()

  call initialize_fields()
  freal(:,1) = -0.5*(phiout-falsevac(1))*tanh((xgrid-rinit)/2.**0.5) + 0.5*(phiout+falsevac(1))
  freal(:,2) = falsevac(2)
  call output(.true.)  ! initialize output file


!  call init_operators()
! Initialize the linear part of the differential operator ( think about this second piece )
  do j=1,nbasis  ! starting at 1, I'm not doing one column?
     XLOOP
        l0(i,j) = bpp(i,j) + 3.*bp(i,j)/xgrid(i)
     ENDXLOOP
  enddo
! transform operator to act in real space not spectral space
  l0 = matmul(l0,mmt)

! Iterate guess until it converges
! While I'm debugging, I'll probably want to output the intermediate iterates
  open(unit=98,file='debug.dat')
  open(unit=71,file='errors.dat')
  write(71,*)"# Max error   Max delta   RMS error"
  print*,"# Residuals before,  after,  step size"

  do l=1,maxit
     call linesolve()
! Add stopping condition
     write(71,*) maxval(abs(error)),maxval(abs(del)),sum(error**2)
     if ( (maxval(abs(error(0:nx-1,:))) < errtol) .and. (maxval(abs(del)) < deltol) ) exit
  enddo
  if (l < maxit) then
     print*,"Converged in ", l," iterations"
     print*,"Total violation of EOM is ",sum(error**2)**0.5
     print*,"Total RMS of last step is ",sum(del**2)**0.5
  else
     print*,"Error, failed to converge in ",maxit," steps"
  endif

  call output()
!  call output_resampled()  ! output for lattice code
contains

  subroutine initialize_fields()
    integer :: i, j

    do i=0,nx
       freal(i,1) = 1.
       freal(i,2) = 0.
    enddo
  end subroutine initialize_fields

!Define the potential derivatives and hessian
#ifdef LINEAR
  function vprime(phi,fld)
    real :: vprime
    real, dimension(nfield), intent(in) :: phi
    integer, intent(in) :: fld

    if (fld == 1) then
       vprime = phi(1)*(phi(1)**2-1.) - delta + g2*(phi(1)-1.)*phi(2)**2
    else
       vprime = epsilon + g2*(phi(1)-1.)**2*phi(2)
    endif
  end function vprime

  function vdprime(phi, fld1, fld2)
   real :: vdprime
    real, dimension(nfield), intent(in) :: phi
    integer, intent(in) :: fld1, fld2

    if (fld1==fld2) then
       if (fld1==1) then
          vdprime = 3.*phi(1)**2 - 1. + g2*phi(2)**2
       else
          vdprime = g2*(phi(1)-1.)**2
       endif
    else
       vdprime = 2.*g2*(phi(1)-1.)*phi(2)
    endif
  end function vdprime
#endif

#ifdef CUBIC
!  function potential(phi)
!    real :: potential
!    real, dimension(nfield), intent(in) :: phi

!    potential = 0.25*(phi(1)**2-1.)**2 + delta*((1./3.)*phi(1)**3 - phi(1)) + 0.5*g2*(phi(1)-1)**2*phi(2)**2 - epsilon*phi(2)
!  end function potential

  function vprime(phi, fld)
    real :: vprime
    real, dimension(nfield), intent(in) :: phi
    integer, intent(in) :: fld

    if (fld == 1) then
       vprime =  (phi(1)+delta)*(phi(1)**2 - 1.) + g2*(phi(1)-1.)*phi(2)**2
    else
       vprime = epsilon + g2*(phi(1)-1.)**2*phi(2)
    endif
  end function vprime

  function vdprime(phi, fld1, fld2)
    real :: vdprime
    real, dimension(nfield), intent(in) :: phi
    integer, intent(in) :: fld1,fld2

    if (fld1 == fld2) then
       if (fld1 == 1) then
         vdprime =  3.*phi(1)**2 - 1. + 2.*delta*phi(1) + g2*phi(2)**2
      else
         vdprime = g2*(phi(1)-1.)**2
      endif
   else
      vdprime = 2.*g2*(phi(1)-1.)*phi(2)
   endif
  end function vdprime
#endif

  subroutine output_debug(init)
    logical, optional :: init
    
    integer :: i
    real, dimension(XRANGE) :: derivs

    if (present(init)) then
       if (init) then
          open(unit=98, file="debug.dat")
          del = 0.
       endif
       return
    endif

    derivs = matmul(l0,freal(:,1))
    XLOOP
       write(98,*) xgrid(i), freal(i,:), del(i,:), error(i,:), derivs(i), vprime(freal(i,:),1)
    ENDXLOOP
    write(98,*)
  end subroutine output_debug

  subroutine output(init)
    logical, optional :: init
    integer :: i

    if (present(init)) then
       if (init) then
          open(unit=99,file="instanton.dat")
       endif
       return
    endif

    do i=1,nfield
       fspec(:,i) = matmul(mmt,freal(:,i))
    enddo
    XLOOP
       write(99,*) delta, xgrid(i), freal(i,:), del(i,:), fspec(i,:), error(i,:)
    ENDXLOOP
    write(99,*)
  end subroutine output

! Linear Solver for Frechet Derivative of the equation
! Currently Solves L(x)f = S(x) with L(x) a linear operator and S(x) a source term
! which is simply a matrix inversion of L(x)
!
! If we have numerical instead of behavioural boundary conditions, then we must "boundary border" that L matrix below or do basis recombination
!

! Newton Line-Solver
  subroutine linesolve()
    real, dimension(1:fsize,1:fsize) :: L
    real, dimension(1:fsize) :: S
    real, dimension(XRANGE,nfield) :: tmp, tmpl
    integer :: ipiv(1:fsize), info

    integer :: i
    real :: b0, b1, bb1  ! residual of step
    real :: alpha, alpha1
    real :: resid(-20:20)
    real, parameter :: alphar(-20:20) = (/ (0.1*i, i=-20,20) /)
    integer, parameter :: nxmax = nx+1
    integer :: n1,n2

!
! Compute linear differential operator and source term
!
    L = 0.
! Linear differential operators (diagonal in spectral basis)
    do i=1,nfield
       n1 = (i-1)*nxmax + 1
       n2 = i*nxmax
       L(n1:n2,n1:n2) = l0(XRANGE,BRANGE)    ! Derivative parts of operator
    enddo
! Potential pieces ("diagonal" in real space", w/ two side bands)
    do i=1,nxmax
       L(i,i) = L(i,i) - vdprime(freal(i-1,:),1,1)
       L(i,i+nxmax) = -vdprime(freal(i-1,:),1,2)
       L(i+nxmax,i) = L(i,i+nxmax)
       L(i+nxmax,i+nxmax) = L(i+nxmax,i+nxmax) - vdprime(freal(i-1,:),2,2)
    enddo
! Now compute source term
    do i=1,nfield
       n1 = (i-1)*nxmax + 1
       n2 = i*nxmax
       S(n1:n2) = -matmul(l0,freal(:,i))
    enddo
! I think this is right
    do j=1,nfield
       do i=1,nxmax
          S(i+(j-1)*nxmax) = S(i+(j-1)*nxmax) + vprime(freal(i-1,:),j)
       enddo
    enddo

! This is a hack, doesn't work for arbitrary numbers of fields
!    do i=1,nxmax
!       n1 = (i-1)*nxmax + 1
!       n2 = i*nxmax
!       S(i) = S(i) + vprime(freal(i-1,:),1)
!       S(i+nxmax) = S(i+nxmax) + vprime(freal(i-1,:),2)
!    enddo

! Set boundary condition at infinity (again, this isn't general for an arbitrary number of fields
    L(nxmax,:) = 0.
    L(nxmax,nxmax) = 1.
    S(nxmax) = 0.

    L(fsize,:) = 0.
    L(fsize,fsize) = 1.
    S(fsize) = 0.

! Check residual of previous step
    b0 = sum(S*S)
    call DGESV(fsize,1,L,fsize,ipiv,S,fsize,info)
    if (info /= 0) then
       print*,"Error inverting linear matrix in DGESV"
       stop
    endif
    do i=1,nfield
       n1=(i-1)*nxmax+1
       n2=i*nxmax
       del(XRANGE,i) = S(n1:n2)
    enddo

! This stuff is just for debugging, although it would be interesting to put in a paper
! output the residual as a function of step we take along the proposed solution
    do i=-20,20
       tmp = freal + alphar(i)*del
       call get_residual(error, tmp)
       resid(i) = sum(error**2)
       write(70,*) alphar(i), resid(i)
    enddo
    write(70,*) 

! Check residual of proposed step
! Clean this part up to do a line search

    b1 = b0 + 0.1
    alpha = 2.
! Iterate over a maximum number of steps
! If we're still stuck and we haven't stopped iterating, probably a local min, so kick ourselves

! Fix this stuff for multifield
!
! This is basically the guts of the Newton iteration.  For multiple fields it might prove useful to let different fields vary at
! different rates.  How do I arrange for this? (even better, use eigenmodes in some way, but how?)
    do i=1,8
       alpha = alpha / 2.
       tmp = freal + alpha*del
       call get_residual(error, tmp)
       b1 = sum(error(0:nx-1,:)**2)
       if (b1 < b0) exit
    enddo

    freal = freal + alpha*del
!    errtot = sum(error**2)
!    deltot = sum(del**2)

    call get_residual(error, freal)

    print*, b0, b1,alpha
    call output_debug()
  end subroutine linesolve

!
! Fill out the source term matrix
!
!  subroutine get_source(source)
!    real, dimension(blah:blah), intent(out) :: source
!    integer :: i,j

!    do j=1,nfield
!       source(:) = 0.
!    enddo 
!  end subroutine get_source

  subroutine get_residual(res, fld)
    real, dimension(XRANGE,nfield), intent(in) :: fld
    real, dimension(XRANGE,nfield), intent(out) :: res
    integer :: i,j
    do j=1,nfield
       res(:,j) = -matmul(l0,fld(:,j))
       XLOOP
       res(i,j) = res(i,j) + vprime(fld(i,:),j)
       ENDXLOOP
    enddo
  end subroutine get_residual

!!!!!!!!!!!!!!
! Root-finding to set the false vacuum minimum
!!!!!!!!!!!!!!
  subroutine get_vacuum(fld)
    real, dimension(nfield), intent(inout) :: fld

    real, dimension(nfield,nfield) :: hessian
    integer :: i,j
    real, dimension(nfield) :: dfld, dpot
    integer :: pivs(1:nfield), info  ! temporary storage for LAPACK
    integer :: l
    integer, parameter :: maxit = 16  ! maximum number of iterations
    real, parameter :: min_tol = 1.e-14

    print*,"initial field is ",fld
    do l=1,maxit
! Do Newton iteration to try and find nearest minimum to initial guess
       do j=1,nfield; do i=1,nfield
          hessian(i,j) = vdprime(fld,i,j)
       enddo; enddo
       do i=1,nfield
          dpot(i) = -vprime(fld,i)
       enddo
       dfld = dpot

       call DGESV(nfield, 1, hessian, nfield, pivs, dfld, nfield,info)
! Check for convergence
       if (sum(abs(dfld)) < min_tol) exit
       fld = fld + dfld
       print*,"fld is ",fld," dfld is ",dfld
    enddo

    if (l.eq.maxit) then
       print*,"Failed to find local minimum of potential.  Adjust initial guess.  Quitting."
       stop
    endif
  end subroutine get_vacuum


!!!!!!!!!!!!!!
! General purpose subroutines for setting up the spectral solver
!!!!!!!!!!!!!!

!
! Subroutines to set up basis functions and transformation matrices, etc.
!
! Input: real - length : parameter L that appears to transform to rational chebyshev's on doubly infinite interval
!        logical - x_to_y : if .true. then input coordinates are on interval [-1:1].  Otherwise on intervals [-infinity,infinity]
!
  subroutine transform_coords(xvals, length, x_to_y)
    real, dimension(XRANGE), intent(inout) :: xvals
    real, intent(in) :: length
    logical, intent(in) :: x_to_y

    integer :: j
!!!!!
! Map to only even Chebychev using x(y)=2*y^2-1 formula for transform
!!!!!
    if (x_to_y) then
       xvals = sqrt(0.5*xvals+0.5)
       map_p = 4.*xvals
       map_pp = 4.
       do j=0,nbasis
          bpp(:,j) = map_p**2*bpp(:,j) + map_pp*bp(:,j)
          bp(:,j) = bp(:,j) * map_p
       enddo
    else
!
! This option needs debugging !!!!!
!
! map_p and map_pp are dy/dx evaluated at the original gridpoints with y=y(x) s.t. B_new(y) = B_old(x(y))
!    map_p = 2.*len/(1-xgrid)**2
!    map_pp = 4.*len/(1-xgrid)**3
!    do j=0,nbasis
!       bp(:,j) = bp(:,j) / map_p
!       bpp(:,j) = ( bpp(:,j) - map_pp*bp(:,j) ) / map_p**2
!    enddo
!    xgrid = len*(1+xgrid)/(1-xgrid)
    endif

!!!!!!
! transform for rational chebychev on doubly-infinite intervals
!!!!!!
! Here's the same transform, but using x(y) instead of y(x) to transform
    if (x_to_y) then
       xvals = length*xvals / sqrt(1-xvals**2)
       map_p = length**2 / (length**2 + xvals**2)**1.5
       map_pp = -3.*length**2*xvals / (length**2 + xvals**2)**2.5
       do j=0,nbasis
          bpp(:,j) = map_p**2*bpp(:,j) + map_pp*bp(:,j)
          bp(:,j) = map_p * bp(:,j)
       enddo
    else
!!!!
! Debug this still
!!!!
!    map_p = length / (1-xvals**2)**1.5
!    map_pp = 3.*length * xvals / (1-xvals**2)**2.5
!    do j=0,nbasis
!       bp(:,j) = bp(:,j) / map_p
!       bpp(:,j) = ( bpp(:,j) - map_pp*bp(:,j) ) / map_p**2
!    enddo
!    xvals = length*xvals/sqrt(1-xvals**2)   
    endif
  end subroutine transform_coords

  subroutine make_transform_matrix()
! Matrix to transform between real and spectral basis, includes weights (This thing is dependent on choice of basis functions, so it should be put in a subroutine
! Used as c = mmt*r,  r=b*c, c: spectral vector, r: real space vector, mmt: transform.  For Chebyshev, probably better to use DCT for speed (worry later)
    do i=0,nbasis
       mmt(:,i) = b(:,i)*wgrid(:)
    enddo
    mmt(:,0) = 0.5*mmt(:,0)  ! integral normalization of T_0 is different
    mmt = transpose(mmt)
  end subroutine make_transform_matrix

  subroutine get_basis_matrix(b,bp,bpp)
    real, dimension(XRANGE, BRANGE) :: b,bp,bpp
    real, dimension(BRANGE) :: c,cp,cpp

    real :: xcur
    integer :: i
    XLOOP
       xcur = xgrid(i)
       call chebychev(nbasis, xcur, c, cp, cpp)
       b(i,:) = c(:)
       bp(i,:) = cp(:)
       bpp(i,:) = cpp(:)
    ENDXLOOP
  end subroutine get_basis_matrix

!
! Initialize collocation points for our solver
! Any transforms should be done in here
!
  subroutine get_collocation_points(xvals, wvals, order)
    real :: xvals(0:order), wvals(0:order)
    integer :: order
    
    call gauss_cheby(xvals, wvals, order)
!    call gauss_lobatto(xvals, wvals, order)
  end subroutine get_collocation_points
!
! Here live some generic subroutines for computing collocation points, etc.  If I wanted I could create a separate module to store these things in.
!
!
! Gauss-Chebychev "zeros" grid on interval [-1:1]
!
  subroutine gauss_cheby(xvals, wvals, order)
    real, intent(out) :: xvals(0:order), wvals(0:order)
    integer, intent(in) :: order

    integer :: i
    real :: dkcol

    dkcol = twopi / 4. /dble(order+1)
    do i=0,order
       xvals(i) = -cos((2*i+1)*dkcol)  ! negative sign just puts the points in increasing order from [-1:1]
    enddo
    wvals = 2. / dble(order+1)
  end subroutine gauss_cheby

!
! Gauss-Lobatto endpoints and extrema grid for Chebychev polynomials
!
  subroutine gauss_lobatto(xvals, wvals, order)
    real, intent(out) :: xvals(XRANGE), wvals(XRANGE)
    integer, intent(in) :: order

    integer :: i
    real :: dkcol
    dkcol = twopi / 2. / dble(order)
    do i=0,order
       xvals(i) = -cos(i*dkcol)
    enddo
    wvals = 2. / dble(order)
    wvals(0) = 1. / dble(order)
    wvals(order) = 1. / dble(order)
  end subroutine gauss_lobatto

!
! Gaussian quadrature abscissas for 
!
  subroutine gauss_lobatto_legendre(xvals, order)
    real, intent(out) :: xvals(0:order)
    integer, intent(in) :: order

  end subroutine gauss_lobatto_legendre

  subroutine gauss_legendre(xvals, order)
    real, intent(out) :: xvals(0:order)
    integer, intent(in) :: order

  end subroutine gauss_legendre

! Subroutines to evaluate desired spectral basis functions and derivatives

!
! Compute Chebychev polynomials up to order nmax at position x
!
  subroutine chebychev(nmax,x,T,Tp,Tpp)
    integer, intent(in) :: nmax
    real, intent(in) :: x
    real, dimension(0:nmax), intent(out) :: T, Tp, Tpp

    real :: xt, cn, sn
    real :: pt, ptt
    integer :: i

! The commented stuff is what Boyd uses in terms of transforming via cosines
!    xt = acos(x)
!    sn = sin(xt)
!    cn = cos(xt)  ! why does Boyd need this?

!    do i=0,nmax
!       T(i) = cos(dble(i)*xt)
!       pt = -dble(i) * sin(dble(i)*xt)
!       ptt = -dble(i)*dble(i) * T(i)
!       Tp(i) = -pt / sn
!       Tpp(i) = (sn*ptt - cn*pt) / sn**3
!    enddo

    T(0) = 1.
    Tp(0) = 0.
    Tpp(0) = 0.
    T(1) = x
    Tp(1) = 1.
    Tpp(1) = 0.
    do i=2,nmax
       T(i) = 2.*x*T(i-1) - T(i-2)
       Tp(i) = 2.*x*Tp(i-1) - Tp(i-2) + 2.*T(i-1)
       Tpp(i) = 2.*x*Tpp(i-1) - Tpp(i-2) + 4.*Tp(i-1)
    enddo
  end subroutine chebychev
 
!
! Compute Legendre polynomials up to order lmax at position x
! To do : add in derivatives etc.
!
  subroutine legendre(P, lmax, x)
    integer :: l, lmax
    real :: P(0:lmax), x

    P(0) = 1.
    P(1) = x
    do l=2,lmax
       P(l) = ((2*l-1)*x*P(l-1) - (l-1)*P(l-2))/l
    enddo
  end subroutine legendre

  subroutine output_resampled()
!    integer, intent(in) :: ngrid  ! number of gridpoints on a side
!    real, intent(in) :: length
    integer, parameter  :: sz = 1024
    real, parameter :: dx = 13.*1.5*2.**0.5/ dble(sz) / 0.2
    integer, parameter :: rsize=3.*sz**2+1
    real, dimension(rsize) :: xvals
    real, dimension(rsize, nfield) :: f
    integer :: i,j
    real :: xtmp(1), ftmp(1)

! Output on a uniformly spaced grid in r^2

    do j=1,nfield
       fspec(:,j) = matmul(mmt,freal(:,j))
       do i=1,rsize
          xtmp = dx*sqrt(dble(i-1))
          call evalf(xtmp,fspec(:,j),ftmp)
          xvals(i) = xtmp(1)
          f(i,j) = ftmp(1)
       enddo
    enddo

    open(unit=80, file="instanton_newgrid.dat")
    do i=1, size(xvals)
       write(80,*) xvals(i), f(i,:)
    enddo
    close(80)    
  end subroutine output_resampled
!
! Auxilliary routines aided in outputting on a noncollocation grid
!
! To do : Since I tend to reach the roundoff plateu, include a new cutoff
!  in here
!
! Input : x - coorinate to evaluate basis at
!       : c - vector of spectral coefficients
!
  subroutine evalf(x, c, f)
    real, intent(inout), dimension(:) :: x
    real, intent(in), dimension(BRANGE) :: c
    real, intent(out), dimension(:) :: f

    real, dimension(BRANGE) :: bvals  ! store values of basis functions
    real :: xtrans
    integer :: i

    do i=1, size(x)
       xtrans = x(i)
       xtrans = xtrans / sqrt(xtrans**2+len**2)     !len*x(i) / sqrt(1-x(i)**2)
       xtrans = 2.*xtrans**2-1.
       call evalbasis(xtrans, bvals)
       f(i) = sum(c*bvals)
    enddo
  end subroutine evalf

  subroutine evalbasis(x, b)
    real, intent(in) :: x
    real, intent(out), dimension(:) :: b

    integer :: ncoeff
    real, dimension(BRANGE) :: tmp1, tmp2
    ncoeff = size(b)

    call chebychev(ncoeff, x, b, tmp1,tmp2)
  end subroutine evalbasis

end program instanton
