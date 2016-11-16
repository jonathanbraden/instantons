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
  real, parameter :: pi = twopi / 2.

! Number of fields, basis functions and collocation points
  integer, parameter :: nfield=1
  integer, parameter :: nmax = 100 !90 ! order of highest interpolating polynomial
  integer, parameter :: nx = nmax,  nbasis=nx
  integer, parameter :: fsize = nx+1  ! needed for DGESV

  real, dimension(XRANGE) :: xgrid, wgrid
  real, dimension(XRANGE,BRANGE) :: b, bp, bpp, mmt
! Solution in real and spectral space
  real, dimension(XRANGE) :: freal, del
  real, dimension(BRANGE) :: fspec

! Storage for various pieces of the differential operators.  These will depend on the precise equation
  real, dimension(XRANGE, BRANGE) :: l0, l1
  real, dimension(XRANGE) :: l2
  real, dimension(XRANGE) :: error  ! store violation of equations of motion
  real :: errtot, deltot
  real, parameter :: errtol = 5.e-10, deltol = 5.e-10
  integer, parameter :: maxit = 10000

  real :: delta
  integer, parameter :: numdelta = 1

  real, parameter, dimension(numdelta) :: deltavals = (/ 0.1 /)
  integer :: m
  real :: phifalse
  real, parameter :: phiout = pi
  real :: rinit, len
  real,dimension(XRANGE) :: map_p, map_pp

  integer :: i,j,l
  real :: lambda = 1.e-3 ! what does this control?

  real :: width,wback
  real :: phif, phit

  ! For Drummond
!  rinit =  3.*2.**0.5*deltavals(1)**0.5
!  phif=0.
!  phit=pi
  
!  len = 10.
!  len=1.6*rinit
  delta = deltavals(1)
  phif = -1.; phit = 1.
  rinit = 1.5*2.**0.5/delta
  len = rinit*3.**0.5
  print*,"length is ",len

  call get_vacuum(phif)
  phifalse = phif
  call get_vacuum(phit)
  call get_collocation_points(xgrid, wgrid, nmax)
  call get_basis_matrix(b, bp, bpp)
  call transform_to_evens(xgrid)
  call cluster_points(xgrid,2.,.true.)
  call transform_double_infinite(xgrid,len)
  call make_transform_matrix()

! Debugging to check transforms
!  call debug_check()

!  call init_guess(freal)
  freal = -0.5*(phit-phif)*tanh((xgrid-rinit)/2.**0.5) + 0.5*(phit+phif)
  print*,""
  print*,"vacua are ",phit,phif
!  print*,"phi is"
!  print*,freal
!  print*,""
  call output(.true.)  ! initialize output file

! Initialize the linear part of the differential operator
  do j=1,nbasis
     XLOOP
        l0(i,j) = bpp(i,j) + 3.*bp(i,j)/xgrid(i)
     ENDXLOOP
  enddo
  l0 = matmul(l0,mmt) ! act in real not spectral space

!  print*,""
!  print*,"X grid is "
!  print*,xgrid
!  print*,""

!  print*,""
!  print*,"L0 is "
!  do i=0,nx
!     print*,l0(i,:)
!  enddo
!  print*,""

! Iterate guess until it converges
! While I'm debugging, I'll probably want to output the intermediate iterates
  open(unit=98,file='debug.dat')

  do l=1,maxit
     print*,"Iterating solver ",l
     call linesolve()
! Stopping condition
     print*,maxval(abs(error)), maxval(abs(del))
     print*,""
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
!  call output_resampled_1d(1024,2.**0.5*13.*10./1024.)
!  call output_resampled_1d(512.,2.**0.5*13.*10./512.)
!  call output_resampled()
!  call get_evalues()
!  call get_evectors()
  
!  call get_action()
  
contains

! Initial guess for the fields
  subroutine initialize_fields()
  end subroutine initialize_fields

  subroutine get_vacuum(fld)
    real, intent(inout) :: fld

    integer, parameter :: maxit=16
    real, parameter :: min_tol=1.e-14
    real :: vp,vpp,dfld

    print*,"initial field is ",fld
    do l=1,maxit
       vpp = vdprime(fld)
       vp = vprime(fld)
       dfld = -vp/vpp
       fld = fld + dfld
       print*,"new field ",fld
       if (abs(dfld) < min_tol) exit
    enddo
    
    if (l.eq.maxit) then
       print*,"Failed to find local minimum of potential. Adust guess"
       stop
    endif
    print*,"Vacuum is ",fld," derivative is ",vprime(fld)
  end subroutine get_vacuum
!
! Define the model through it's potential
  !
  elemental function potential(phi)
    real :: potential
    real, intent(in) :: phi

    potential = 0.25*(phi**2-1.)**2 + delta*(phi**3/3. - phi)
!    potential = cos(phi) + delta*sin(phi)**2 + 1.
  end function potential
  
  elemental function vprime(phi)
    real :: vprime
    real, intent(in) :: phi

    vprime = (phi+delta)*(phi**2-1.)
!    vprime =  -sin(phi) + delta*sin(2.*phi)
  end function vprime

  elemental function vdprime(phi)
    real :: vdprime
    real, intent(in) :: phi

    vdprime = 3.*phi**2 - 1. + 2.*delta*phi
!    vdprime =  -cos(phi) + 2.*delta*cos(2.*phi)
  end function vdprime

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

    derivs = matmul(l0,freal)
    XLOOP
       write(98,*) xgrid(i), freal(i), del(i), error(i), l2(i)
    ENDXLOOP
    write(98,*)
  end subroutine output_debug

  subroutine output(init)
    logical, optional :: init
    integer :: i
    real, dimension(XRANGE) :: fderiv, fderiv2

    if (present(init)) then
       if (init) then
          open(unit=99,file="instanton.dat")
       endif
       return
    endif
    fspec = matmul(mmt,freal)
    fderiv = matmul(bp,fspec)
    fderiv2 = matmul(bpp,fspec)
    XLOOP
       write(99,*) delta, xgrid(i), freal(i), del(i), fspec(i), error(i), 3.*fderiv(i)/xgrid(i), fderiv2(i), vprime(freal(i))
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

  subroutine linesolve()
    real, dimension(1:fsize,1:fsize) :: L
    real, dimension(1:fsize) :: S
    real, dimension(XRANGE) :: tmp, tmpl
    integer :: ipiv(1:fsize), info

    integer :: i
    real :: b0, b1, bb1  ! residual of step
    real :: alpha, alpha1
    real :: resid(-20:20)
    real, parameter :: alphar(-20:20) = (/ (0.1*i, i=-20,20) /)

! Get coefficients for linearized problem Ax=B
    L(1:fsize,1:fsize) = l0(XRANGE,BRANGE)
    do i=1,fsize
      L(i,i) = L(i,i) - vdprime(freal(i-1)) 
    enddo
    l2 =  vprime(freal)
    S(1:fsize) = -matmul(l0,freal) + l2(XRANGE)
    tmpl = S
    print*,"Residual is ",sqrt(sum(S**2))

! Set boundary condition at infinity
    L(fsize,:) = 0.
    L(fsize,fsize) = 1.
    S(fsize) = 0.

! Check residual of previous step
    b0 = sum(S*S)
    call DGESV(fsize,1,L,fsize,ipiv,S,fsize,info); del = S
    if (info /= 0) then
       print*,"Error inverting linear matrix in DGESV"
       stop
    endif

!    print*,"del vector is "; print*,del

! Check residual of proposed step
! Clean this part up to do a line search

    b1 = b0 + 0.1
    alpha = 2.
! Iterate over a maximum number of steps
! If we're still stuck and we haven't stopped iterating, probably a local min, so kick ourselves
    do i=1,8
       alpha = alpha / 2.
       tmp = freal + alpha*del
       S = -matmul(l0,tmp) + vprime(tmp)
       S(fsize) = 0.
       b1 = sum(S*S)
       print*,"b1 is ",b1," b0 is ",b0," alpha is ",alpha
       if (b1 < b0) exit
    enddo

    freal = freal + alpha*del

    l2 = vprime(freal)
    error = -matmul(l0,freal) + vprime(freal)
    errtot = sum(error**2)
    deltot = sum(del**2)

    print*, b0, b1, alpha
    call output_debug()
  end subroutine linesolve

! Do do, pass in function instead of directly using freal
  subroutine get_evalues( )
    real, dimension(1:fsize,1:fsize) :: L
    real, dimension(1:fsize) :: evalreal, evalimag
!    real, dimension(1:fsize) :: evecl, evecr
    integer :: asize
    integer :: ierror
    real :: dummy(1:1)
    real, allocatable, dimension(:) :: work
    integer :: iwork, i, imax(1)

    asize = fsize ! (or nx or something)
    ! Allocate workspace
    call DGEEV('N','N',asize,L,asize, evalreal, evalimag, dummy,1,dummy,1,dummy,-1,ierror)
    if (ierror.eq.0) then
       iwork = int(dummy(1))
       allocate(work(1:iwork))
    else
       print*,"Error allocating workspace for eigenvalue problem, exiting"
       stop
    endif

! Fill second variation of action
    L(1:fsize,1:fsize) = -l0(XRANGE,BRANGE)
    do i=1,fsize
       L(i,i) = L(i,i) + vdprime(freal(i-1))
    enddo
! Hmm, I should add some boundary conditions in here

! Get eigenvalues
    call DGEEV('N','N',asize,L,asize,evalreal, evalimag, dummy,1,dummy,1,work,iwork,ierror)
    if (ierror /= 0) then
       print*,"Error in eigenvalue solution"
       stop
    endif

    ! Print out 10 smallest eigenvalues
    do i=1,10
       imax = minloc(evalreal)
       print*,evalreal(imax(1)), evalimag(imax(1))
       evalreal(imax(1))=1000000.
    enddo

    do i=1,10
       imax = maxloc(evalimag)
       print*,evalreal(imax(1)), evalimag(imax(1))
       evalimag(imax(1)) = -1.
    enddo
  end subroutine get_evalues

  subroutine get_evectors()
    real, dimension(1:fsize,1:fsize) :: L
    real, dimension(1:fsize) :: evalreal, evalimag
!    real, dimension(1:fsize) :: evecl, evecr
    integer :: asize
    integer :: ierror
    real :: dummy(1:1)
    real :: evec(1:fsize,1:fsize)
    real, allocatable, dimension(:) :: work
    integer :: iwork, i, imax(1), j

    asize = fsize ! (or nx or something)
    ! Allocate workspace
    call DGEEV('N','V',asize,L,asize, evalreal, evalimag, dummy,1,evec,fsize,dummy,-1,ierror)
    if (ierror.eq.0) then
       iwork = int(dummy(1))
       allocate(work(1:iwork))
    else
       print*,"Error allocating workspace for eigenvalue problem, exiting"
       stop
    endif

! Fill second variation of action
    L(1:fsize,1:fsize) = -l0(XRANGE,BRANGE)
    do i=1,fsize
       L(i,i) = L(i,i) + vdprime(freal(i-1))
    enddo

! Get eigenvalues
    call DGEEV('N','V',asize,L,asize,evalreal, evalimag, dummy,1,evec,fsize,work,iwork,ierror)
    if (ierror /= 0) then
       print*,"Error in eigenvalue solution"
       stop
    endif

    ! Print out 10 smallest eigenvalues
    open(unit=80,file="eigenvectors.dat")
    do i=1,10
       imax = minloc(evalreal)
       print*,evalreal(imax(1)), evalimag(imax(1))
       evalreal(imax(1))=1000000.
       do j=1,fsize
          write(80,*) evec(j,imax(1))
       enddo
       write(80,*)
    enddo
  end subroutine get_evectors
  
  subroutine get_action()
    real, dimension(XRANGE) :: dfsq, pot, lag
    real :: action
    integer, parameter :: ndim = 3
    real, parameter :: sphere_surface = 2.*pi**2 ! adjust for dimension
    dfsq = matmul(bp,fspec)
    dfsq = 0.5*dfsq**2
    lag = dfsq + potential(freal) ! Euclidean Lagrangian
    
    action = sum(wgrid*lag*xgrid**ndim) ! Modify to allow arbitrary dimensions
    action = sphere_surface * action
    print*,"action is ",action
!    print*,dfsq
  end subroutine get_action
  
  ! Initial guess for the field profile
  ! Use the thin-wall solution, or else previously computed profile for different adjustable model parameter
  elemental function profile(x)
    real :: profile
    real, intent(in) :: x

    real, parameter :: rinit = 0. !1./delta

    profile = tanh((x-rinit)/2.**0.5)
  end function profile

  subroutine transform_to_evens(xvals)
    real, dimension(XRANGE), intent(inout) :: xvals
    integer :: j
    
    xvals = sqrt(0.5*xvals+0.5)
    map_p = 4.*xvals
    map_pp = 4.
    do j=0,nbasis
       bpp(:,j) = map_p**2*bpp(:,j) + map_pp*bp(:,j)
       bp(:,j) = bp(:,j) * map_p
    enddo
  end subroutine transform_to_evens

  subroutine transform_double_infinite(xvals, length)
    real, dimension(XRANGE), intent(inout) :: xvals
    real, intent(in) :: length
    integer :: j

    xvals = length*xvals / sqrt(1.-xvals**2)
    map_p = length**2 / (length**2 + xvals**2)**1.5
    map_pp = -3.*length**2*xvals / (length**2 + xvals**2)**2.5
    do j=0,nbasis
       bpp(:,j) = map_p**2*bpp(:,j) + map_pp*bp(:,j)
       bp(:,j) = map_p * bp(:,j)
    enddo
  end subroutine transform_double_infinite

  subroutine cluster_points(xvals, w, evens)
    real, dimension(XRANGE), intent(inout) :: xvals
    real, intent(in) :: w
    logical, intent(in) :: evens

    integer :: j
    real, dimension(XRANGE) :: xshift
    real :: winv

    winv = 1./w
    if (evens) then
       xshift = pi*(xvals-0.5)
       xvals = atan(w*tan(xshift)) / pi + 0.5
       xshift = pi*(xvals-0.5) ! why did I do this again?
       map_p = winv*(1.+tan(xshift)**2) / (1. + winv**2*tan(xshift)**2)
       map_pp = 2.*pi*winv*(1.-winv**2)*tan(xshift)*(1.+tan(xshift)**2)/(1.+winv**2*tan(xshift)**2)**2
       do j=0, nbasis
          bpp(:,j) = map_p**2*bpp(:,j) + map_pp*bp(:,j)
          bp(:,j) = map_p * bp(:,j)
       enddo
       wback = winv
    else
       print*,"Error, transform currently only works to cluster even Chebyshev points"
    endif
  end subroutine cluster_points
!
! Subroutines to set up basis functions and transformation matrices, etc.
!
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
    integer, intent(in) :: order
    
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

    xvals = 0.
  end subroutine gauss_lobatto_legendre

  subroutine gauss_legendre(xvals, order)
    real, intent(out) :: xvals(0:order)
    integer, intent(in) :: order

    xvals = 0.
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

  subroutine output_resampled()
    integer, parameter :: sz = 1024
    real, parameter :: dx = 6.5*10.*2.**0.5/dble(sz)    !0.125
    integer, parameter :: rsize=3.*sz**2+1
    real, dimension(rsize) :: xvals
    real, dimension(rsize) :: f
    integer :: i
    real :: xtmp(1), ftmp(1)

! Output on a uniformly spaced grid in r^2

    fspec = matmul(mmt,freal)
    do i=1,rsize
       xtmp = dx*sqrt(dble(i-1))
       call evalf(xtmp,fspec,ftmp)
       xvals(i) = xtmp(1)
       f(i) = ftmp(1)
    enddo
!    do i=1,sz
!       xvals(i) = dx**2*(i-1)
!    enddo
!    call evalf(xvals, fspec, f)

    open(unit=80, file="instanton_newgrid.dat")
    do i=1, size(xvals)
       write(80,*) xvals(i), f(i)
    enddo
    close(80)    
  end subroutine output_resampled

  subroutine output_resampled_1d(ngrid, dx)
    integer, intent(in) :: ngrid
    real, intent(in) :: dx
    
    real, dimension(ngrid) :: xvals
    real, dimension(ngrid) :: f, fp, fpp
    integer :: i,j
    real :: xtmp(1)
    real, dimension(3) :: ftmp

    fspec = matmul(mmt,freal)
    do i=1,ngrid
       xtmp = dx*(i-1)
       call evalf(xtmp,fspec(:),ftmp)
       f(i) = ftmp(1)
       fp(i) = ftmp(2)
       fpp(i) = ftmp(3)
       xvals(i) = xtmp(1)
    enddo

    open(unit=80, file="instanton_offgrid.dat")
! Avoid singularity at origin
    write(80,*) xvals(1), f(1), 0., fpp(1), vprime(f(1)), 0.
    do i=2, ngrid
       write(80,*) xvals(i), f(i), 3.*fp(i)/xvals(i), fpp(i), vprime(f(i)), fpp(i) + 3.*fp(i)/xvals(i)-vprime(f(i))
    enddo
  end subroutine output_resampled_1d

!
! Auxilliary routines aided in outputting on a noncollocation grid
!
! To do : Since I tend to reach the roundoff plateu, include a new cutoff
!  in here
!
! Input : x - coorinate to evaluate at
!       : c - vector of spectral coefficients
! Output : f - function, and first and second derivatives
!
  subroutine evalf(x, c, f)
    real, intent(inout), dimension(:) :: x
    real, intent(in), dimension(BRANGE) :: c
    real, intent(out), dimension(3,size(x)) :: f

    real, dimension(3,BRANGE) :: bvals  ! store values of basis functions and derivatives
    real :: xtrans, dx, ddx
    real :: wcur, xshift
    integer :: i

    wcur = 1./wback
    do i=1, size(x)
       xtrans = x(i)
       xtrans = xtrans / sqrt(xtrans**2+len**2)
       xtrans = atan( wback*tan(pi*(xtrans-0.5)) )/pi + 0.5
       xtrans = 2.*xtrans**2-1.
       call evalbasis(xtrans, bvals,size(c)-1)
 ! Now add variable change for derivatives
       xtrans=sqrt(0.5*xtrans+0.5); dx=4.*xtrans; ddx=4.
       bvals(3,:) = dx**2*bvals(3,:) + ddx*bvals(2,:)
       bvals(2,:) = dx*bvals(2,:)
       xtrans = atan(wcur*tan(pi*(xtrans-0.5)))/pi + 0.5
       xshift = pi*(xtrans-0.5)
       dx = wback*(1.+tan(xshift)**2) / (1.+wback**2*tan(xshift)**2)
       ddx = 2.*pi*wback*(1.-wback**2)*tan(xshift)*(1.+tan(xshift)**2)/(1.+wback**2*tan(xshift)**2)**2
       bvals(3,:) = dx**2*bvals(3,:) + ddx*bvals(2,:)
       bvals(2,:) = dx*bvals(2,:)
! Transform to infinite interval
       xtrans=len*xtrans/sqrt(1.-xtrans**2)
       dx=len**2/(len**2+xtrans**2)**1.5
       ddx=-3.*len**2*xtrans / (len**2 + xtrans**2)**2.5
       bvals(3,:) = dx**2*bvals(3,:) + ddx*bvals(2,:)
       bvals(2,:) = dx*bvals(2,:)

       f(1,i) = sum(c*bvals(1,:))
       f(2,i) = sum(c*bvals(2,:))
       f(3,i) = sum(c*bvals(3,:))
    enddo

  end subroutine evalf

  subroutine evalbasis(x, b, ncoeff)
    integer, intent(in) :: ncoeff
    real, intent(in) :: x
    real, intent(out), dimension(3,0:ncoeff) :: b
!    real, intent(out), dimension(3,ncoeff) :: b

    call chebychev(ncoeff, x, b(1,:), b(2,:), b(3,:))
  end subroutine evalbasis

end program instanton
