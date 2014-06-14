!
! Code to obtain instanton profiles (and in the future actions, etc.)
! in scalar field theories
!
! Compile with
!  gfortran -fdefault-real-8 -fdefault-double-8 -O3 -xhost -o instanton instanton.f90 -llapack -lm
!


#define XLOOP do i=0,nx
#define ENDXLOOP enddo
#define XRANGE 0:nx
#define BRANGE 0:nbasis

program instanton
  implicit none

  real, parameter :: twopi =  6.283185307179586476925867665590
  real, parameter :: pi = 0.5*twopi

! Number of fields, basis functions and collocation points
  integer, parameter :: nfield=1
  integer, parameter :: nmax = 511 ! order of highest interpolating polynomial
  integer, parameter :: nx = nmax,  nbasis=nx
  integer, parameter :: lsize = nx+1         ! 
  integer, parameter :: fsize = (nfield+1)*(nx+1) + 1 ! needed for DGESV, total number of variables (nx+1 per field, nx+1 for scale factor, 1 for size of euclidean dS)

  real, dimension(XRANGE) :: xgrid, wgrid
  real, dimension(XRANGE,BRANGE) :: b, bp, bpp, mmt
! Solution in real and spectral space (have to fix for multiple fields)
  real, dimension(XRANGE) :: freal
  real, dimension(BRANGE) :: fspec
  real, dimension(XRANGE) :: areal
  real, dimension(BRANGE) :: aspec
  real :: tmax
  real, dimension(1:fsize) :: del

! Storage for various pieces of the differential operators.  These will depend on the precise equation
  real, dimension(XRANGE, BRANGE) :: D2D2, D1D1
! Store various splittings of our operators
  real, dimension(XRANGE, BRANGE) :: l0, l1
  real, dimension(XRANGE) :: l2
  real, dimension(1:fsize) :: error  ! store violation of equations of motion
  real :: errtot, deltot
  real, parameter :: errtol = 1.e-10, deltol = 1.e-10
  integer, parameter :: maxit = 1000

  real :: delta
  real, parameter :: epsilon = 0.1  ! controls strength of gravity (basically phi_vev / M_pl)
  integer, parameter :: numdelta = 1
  real, parameter, dimension(numdelta) :: deltavals = (/ 0. /)
  integer :: m
  real, parameter :: phifalse = -1.
  real, parameter :: phiout = 1.
  real :: rinit, len, V0
  real,dimension(XRANGE) :: map_p, map_pp

  integer :: i,j,l

  real :: lambda = 1.e-3

  delta = deltavals(1)

! Set up the spectral matrices
  call get_collocation_points(xgrid, wgrid, nmax)
  call get_basis_matrix(b, bp, bpp)
!  call transform_coords() 
  call make_transform_matrix()

  D2D2 = matmul(bpp,mmt)
  D1D1 = matmul(bp,mmt)

! Debugging to check transforms
!  call debug_check()

! Now do all of the initialization of the fields
!  call initialize_fields()
!  call initialize_gravity()
  V0 = potential(0.)
  freal = 0.01*sin(0.5*pi * xgrid)
  tmax = 0.5 * (pi / sqrt(epsilon*V0) ) ! my definition gives tmax as half the value of tau_max
  areal = (2.*tmax/pi)* cos(0.5*pi* xgrid )

  print*,"tmax is ",tmax," V0 is ",V0

  call output(.true.)  ! initialize output file
  print*,"Done initializing"

! Iterate guess until it converges
! While I'm debugging, I'll probably want to output the intermediate iterates
  open(unit=98,file='debug.dat')
  call output_debug()

  open(unit=71,file='errors.dat')
  write(71,*)"# Max error    Max delta   RMS error"
  print*,"# Residuals before,    after,   step size"

  do l=1,maxit
     call linesolve()
! Add stopping condition
     write(71,*) maxval(abs(error)), maxval(abs(del)), sum(error**2)
     print*,"tmax is ",tmax
     if ( (maxval(abs(error(:))) < errtol) .and. (maxval(abs(del(:))) < deltol) ) exit
  enddo
  if (l < maxit) then
     print*,"Converged in ", l," iterations"
     print*,"Total violation of EOM is ",errtot
     print*,"Total RMS of last step is ",deltot
  else
     print*,"Error, failed to converge in ",maxit," steps"
  endif

  call output()
contains

!
! Define the model through it's potential
! Since I'm assuming a finite domain, make sure this has a positive minimum (to get deSitter, deSitter transitions)
!
  elemental function potential(phi)
    real :: potential
    real, intent(in) :: phi

    potential = 0.25*(phi**2-1.)**2 + delta/3.*(phi**3 - 3.*phi + 2.) + 80.
  end function potential

  function vprime(phi)
    real :: vprime
    real, intent(in) :: phi

    vprime =  (phi+delta)*(phi**2 - 1.)
  end function vprime

  function vdprime(phi)
    real :: vdprime
    real, intent(in) :: phi

    vdprime =  3.*phi**2 - 1. + 2.*delta*phi
  end function vdprime

  subroutine output_debug(init)
    logical, optional :: init
    
    integer :: i
    real, dimension(XRANGE) :: aprime, fprime

    if (present(init)) then
       if (init) then
          open(unit=98, file="debug.dat")
          del = 0.
       endif
       return
    endif

!    aspec = matmul(D1D1,areal)
!    fspec = matmul(D1D1,freal)
!    error = aspec**2 - tmax**2 + epsilon*areal**2*(0.5*fspec**2-tmax**2*potential(freal))

    aspec = matmul(D2D2, areal)
    fspec = matmul(D2D2, freal)
    fprime = matmul(D1D1, freal)
    aprime = matmul(D1D1, areal)
    error(1:lsize) = areal*fspec + 3.*aprime*fprime - tmax**2*vprime(freal)
    error(lsize+1:2*lsize) = aspec + areal*epsilon*(0.5*fspec**2 + tmax**2*potential(freal))

    XLOOP
       write(98,*) xgrid(i), freal(i), areal(i), error(i+1), error(lsize+i+1), aspec(i), fprime(i), tmax
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
    endif

    fspec = matmul(mmt,freal)
    aspec = matmul(mmt,areal)
    XLOOP
       write(99,*) delta, xgrid(i), freal(i), areal(i), fspec(i), aspec(i), error(i+1)
    ENDXLOOP
    write(99,*)

  end subroutine output

  subroutine debug_check()
    freal = b(:,6)
    l2 = matmul(mmt,freal)
    do i=0,nx; write(91,*) l2(i); enddo
       
    fspec = 0.
    fspec(5) = 1.
    l2 = matmul(b, fspec)
    do i=0,nx; write(90,*) xgrid(i), l2(i); enddo
          
! Debugging, output the basis vectors
    open(unit=97,file="basis.dat")
    do i=0,nbasis
       do j=0,nx
          write(97,*) xgrid(j), b(j,i), bp(j,i), bpp(j,i)
       enddo
       write(97,*)
    enddo
    close(unit=97)
  end subroutine debug_check

! Linear Solver for Frechet Derivative of the equation
! Currently Solves L(x)f = S(x) with L(x) a linear operator and S(x) a source term
! which is simply a matrix inversion of L(x)

  subroutine linesolve()
    real, dimension(1:fsize,1:fsize) :: L
    real, dimension(1:fsize) :: S

    real, dimension(1:lsize) :: aprime, adprime, fprime, fdprime, hubble
    real, dimension(1:lsize) :: tmpf, tmpa
    real :: tmptmax
    real, dimension(1:lsize) :: eomf, eoma

    integer :: ipiv(1:fsize), info

    integer :: i, j
    real :: b0, b1, bb1  ! residual of step
    real :: alpha, alpha1
    real :: resid(-20:20)
    real, parameter :: alphar(-20:20) = (/ (0.1*i, i=-20,20) /)

    fprime = matmul(D1D1,freal)
    fdprime = matmul(D2D2,freal)
    aprime = matmul(D1D1,areal)
    adprime = matmul(D2D2,areal)
    hubble(1) = 0.
    hubble(2:lsize-1) = aprime(2:lsize-1) / areal(1:nx-1) 
    hubble(lsize) = 0.

!!!!!
! This part needs fixing for multiple fields
!!!!!
! Construct the linearized matrix operator (basically Frechet derivative)
    L = 0.
    L(1:lsize,1:lsize) = D2D2(XRANGE,BRANGE)
! Neuman b.c. at 0
    L(1,1:lsize) = D1D1(0,BRANGE)
    L(1,lsize+1:fsize) = 0.
    do i=2,lsize-1
      L(i,1:lsize) = L(i,1:lsize) + 3.*hubble(i)*D1D1(i-1,:)
      L(i,i) = L(i,i) - tmax**2*vdprime(freal(i-1))
      L(i,lsize+1:2*lsize) = (3./areal(i-1))*fprime(i)*( D1D1(i-1,BRANGE) - hubble(i) )
      L(i,fsize) = -2.*tmax*vprime(freal(i-1))
    enddo
! Neuman b.c. at 1
    L(lsize,1:lsize) = D1D1(nx,BRANGE)
    L(lsize,lsize+1:fsize) = 0.

! Now construct the part acting on scale factor
! b.c.'s are a(0) = a(1) = 0.
    L(lsize+1,:)=0.
    L(lsize+1,lsize+1) = 1.
    do i=lsize+2,2*lsize-1
       j = i - lsize
       L(i,lsize+1:2*lsize) = D2D2(j-1,BRANGE)
       L(i,i) = L(i,i) + epsilon*( 0.5*fprime(j)**2 + tmax**2*potential(freal(j-1)) )
!       L(i,i) = L(i,i) - adprime(j)/areal(j-1)
       L(i,1:lsize) = epsilon*areal(j-1)*( fprime(j)*D1D1(j-1,BRANGE) + tmax**2*vprime(freal(j-1)) )
       L(i,fsize) = 2.*epsilon*areal(j-1)*potential(freal(j-1))
    enddo
    L((nfield+1)*lsize,:) = 0.
    L((nfield+1)*lsize,(nfield+1)*lsize) = 1.

! Finally, impose aprime=t_max at the left endpoint
    L(fsize,1:lsize) = 0.
    L(fsize,lsize+1:2*lsize) = D1D1(0,BRANGE)
    L(fsize,fsize) = -1.

! Now construct source term (this needs fixing)
    S(1) = 0.
    S(2:lsize-1) = -fdprime(2:lsize-1) - 3.*hubble(2:lsize-1)*fprime(1:nx-1) + tmax**2*vprime(freal(1:nx-1))
    S(lsize) = 0.
    S(lsize+1) = 0.
    S(lsize+2:2*lsize-1) = -adprime(2:lsize-1) - epsilon*areal(1:nx-1)*( 0.5*fprime(1:nx-1)**2 + tmax**2*potential(freal(1:nx-1)) )   ! 0.  ! Could be replaced by appropriate equation (since we're solving for the constraint, it won't be satisfied at intermediate steps in general
    S(2*lsize) = 0.
    S(fsize) = -aprime(1) + tmax

! Is b0 actually giving the EOM now?
    b0 = sum(S*S)
    print*,"Final source is ",S(fsize)," eom error is ", b0

    call DGESV(fsize,1,L,fsize,ipiv,S,fsize,info); del=S
    if (info /= 0) then
       print*,"Error inverting DGESV"
       stop
    endif

    fprime = matmul(D1D1, del(1:lsize))
    print*,"derivatives at end ", fprime(1), fprime(lsize)

! Ok, edit this below so I'm actually computing the correct norm here
    alpha = 0.125  ! 2.
! Iterate over a maximum number of steps
! If we're still stuck and we haven't stopped iterating, probably a local min, so kick ourselves
    do i=1,4
       alpha = alpha / 2.
       tmpf = freal + alpha*del(1:lsize)
       tmpa = areal + alpha*del(lsize+1:2*lsize)
       tmptmax = tmax + alpha*del(fsize)

! Now get the two equations I'm trying to solve
       fprime = matmul(D1D1,tmpf)
       fdprime = matmul(D2D2,tmpf)
       aprime = matmul(D1D1,tmpa)
       adprime = matmul(D2D2,tmpa)

       eomf(1) = 0.
       eomf(2:lsize-1) = fdprime(2:lsize-1) +  &
               3.*aprime(2:lsize-1)*fprime(2:lsize-1) / tmpa(2:lsize-1) - &
               tmptmax**2*vprime(tmpf(2:lsize-1))
       eomf(lsize) = 0.

       eoma = adprime(:) + epsilon*tmpa(:)*( 0.5*fprime(:)**2 + tmptmax**2*potential(tmpf(:)) )

       b1 = sum(eomf*eomf) + sum(eoma*eoma) + (aprime(1) - tmptmax)**2

! For debugging, plot these residuals
       write(70,*) alpha, b1, sum(eomf**2), sum(eoma**2)
       if (b1 < b0) exit
    enddo
    write(70,*)

    freal = freal + alpha*del(1:lsize)
    areal = areal + alpha*del(lsize+1:2*lsize)
    tmax = tmax + alpha*del(fsize)

    error(1:lsize) = eomf
    error(lsize+1:2*lsize) = eoma
    error(fsize) = (aprime(1)-tmax)

!    errnorm = sum(error**2)
!    errmax = maxval(abs(error))
    deltot = sum(del**2)

    print*, b0, b1, alpha
    call output_debug()
  end subroutine linesolve

  ! Initial guess for the field profile
  ! Use the thin-wall solution, or else previously computed profile for different adjustable model parameter
  elemental function profile(x)
    real :: profile
    real, intent(in) :: x

    real, parameter :: rinit = 0. !1./delta

    profile = tanh((x-rinit)/2.**0.5)
  end function profile

!
! Subroutines to set up basis functions and transformation matrices, etc.
!
  subroutine transform_coords()
! Map to only even Chebychev using x(y)=2*y^2-1 formula for transform
    xgrid = sqrt(0.5*xgrid+0.5)
    map_p = 4.*xgrid
    map_pp = 4.
    do j=0,nbasis
       bpp(:,j) = map_p**2*bpp(:,j) + map_pp*bp(:,j)
       bp(:,j) = bp(:,j) * map_p
    enddo

! map_p and map_pp are dy/dx evaluated at the original gridpoints with y=y(x) s.t. B_new(y) = B_old(x(y))
!    map_p = 2.*len/(1-xgrid)**2
!    map_pp = 4.*len/(1-xgrid)**3
!    do j=0,nbasis
!       bp(:,j) = bp(:,j) / map_p
!       bpp(:,j) = ( bpp(:,j) - map_pp*bp(:,j) ) / map_p**2
!    enddo
!    xgrid = len*(1+xgrid)/(1-xgrid)

! transform for rational chebychev on doubly-infinite intervals
!    map_p = len / (1-xgrid**2)**1.5
!    map_pp = 3.*len * xgrid / (1-xgrid**2)**2.5
!    do j=0,nbasis
!       bp(:,j) = bp(:,j) / map_p
!       bpp(:,j) = ( bpp(:,j) - map_pp*bp(:,j) ) / map_p**2
!    enddo
!    xgrid = len*xgrid/sqrt(1-xgrid**2)

! Here's the same transform, but using x(y) instead of y(x) to transform
    xgrid = len*xgrid / sqrt(1-xgrid**2)
    map_p = len**2 / (len**2 + xgrid**2)**1.5
    map_pp = -3.*len**2*xgrid / (len**2 + xgrid**2)**2.5
    do j=0,nbasis
       bpp(:,j) = map_p**2*bpp(:,j) + map_pp*bp(:,j)
       bp(:,j) = map_p * bp(:,j)
    enddo
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
    
!    call gauss_cheby(xvals, wvals, order)
    call gauss_lobatto(xvals, wvals, order)
  end subroutine get_collocation_points

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

end program instanton
