! Compile with
!  gfortran -fdefault-real-8 -fdefault-double-8 -O3 -xhost -o instanton instanton.f90 -llapack -lm
!

#define CUBIC 1
!#define LINEAR 1

#define XLOOP do i=0,nx
#define ENDXLOOP enddo
#define XRANGE 0:nx
#define BRANGE 0:nbasis

program instanton
  use constants
  use Cheby
  use Nonlinear_Solver
  implicit none

! Number of fields, basis functions and collocation points
  integer, parameter :: nfield=1
  integer, parameter :: nmax = 90 ! order of highest interpolating polynomial
  integer, parameter :: nx = nmax,  nbasis=nx
  integer, parameter :: fsize = nx+1  ! needed for DGESV

! Storage for various pieces of the differential operators.  These will depend on the precise equation
  real, dimension(XRANGE, BRANGE) :: l0, l1
  real, dimension(XRANGE) :: l2
  real, dimension(XRANGE) :: error  ! store violation of equations of motion
  real :: errtot, deltot
  real, parameter :: errtol = 5.e-10, deltol = 5.e-10
  integer, parameter :: maxit = 10000

  real :: delta
  integer, parameter :: numdelta = 1
  real, parameter, dimension(numdelta) :: deltavals = (/ 0.99 /)
  integer :: m
  real :: phifalse
  real, parameter :: phiout = 1.
  real :: rinit, len
  real,dimension(XRANGE) :: map_p, map_pp

  integer :: i,j,l
  real :: lambda = 1.e-3  ! this won't work here

  real :: width,wback ! wback is needed in evalf
  real :: phif, phit

! Added with refactoring.  The stuff above this should eventually disappear
  type(Chebyshev) :: transform
  type(Solver) :: solv
  real(dl), parameter :: wp = 2._dl

#ifdef CUBIC
  rinit =  1.5*2.**0.5/deltavals(1)
  phif=-1.
  phit=1.
#endif
#ifdef LINEAR
  rinit = 2.**0.5/deltavals(1)
  phif=-1.
  phit=1.
#endif
  len = 10.
!  len=1.6*rinit
!  len = rinit *3.**0.5
  print*,"length is ",len

  delta = deltavals(1)
  call get_vacuum(phif)
  phifalse = phif
  call get_vacuum(phit)

  ! Set up our spectral derivatives and transform
  call create_chebyshev(transform,nmax,2,.false.,.true.)  ! Check the last two flags
  call transform_to_evens(transform)
  call cluster_points(transform,w,.true.)
  call transform_double_infinite(transform,len)

  call initialise_solver(solv,nmax,100,0.1_dl)

! Debugging to check transforms
!  call debug_check()

!  call init_guess(freal)
  freal = -0.5*(phit-phif)*tanh((xgrid-rinit)/2.**0.5) + 0.5*(phit+phif)
  call output(.true.)  ! initialize output file

! Initialize the linear part of the differential operator
  do j=1,nbasis
     XLOOP
        l0(i,j) = bpp(i,j) + 3.*bp(i,j)/xgrid(i)
     ENDXLOOP
  enddo
  l0 = matmul(l0,mmt) ! act in real not spectral space

! Iterate guess until it converges
! While I'm debugging, I'll probably want to output the intermediate iterates
  open(unit=98,file='debug.dat')

  do l=1,maxit
     call linesolve()
 
! Stopping condition
     print*,maxval(abs(error)), maxval(abs(del))
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
  call output_resampled_1d(1024,2.**0.5*13.*10./1024.)
!  call output_resampled_1d(512.,2.**0.5*13.*10./512.)
!  call output_resampled()
  call get_evalues()

  call delete_solver(solv)
  call 

contains

! Initial guess for the fields
  subroutine initialize_fields()
  end subroutine initialize_fields
!
! Define the model through it's potential
!
  elemental function vprime(phi)
    real :: vprime
    real, intent(in) :: phi

#ifdef CUBIC
    vprime =  (phi+delta)*(phi**2 - 1.)
#endif
#ifdef LINEAR
    vprime = phi*(phi**2-1.) - delta
#endif
  end function vprime

  elemental function vdprime(phi)
    real :: vdprime
    real, intent(in) :: phi

#ifdef CUBIC
    vdprime =  3.*phi**2 - 1. + 2.*delta*phi
#endif
#ifdef LINEAR
    vdprime = 3.*phi**2 - 1.
#endif
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

#ifdef OLD
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
       b1 = sum(S*S)
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
#endif
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
  end subroutine get_evalues

  ! Initial guess for the field profile
  ! Use the thin-wall solution, or else previously computed profile for different adjustable model parameter
  elemental function profile(x)
    real :: profile
    real, intent(in) :: x

    real, parameter :: rinit = 0. !1./delta

    profile = tanh((x-rinit)/2.**0.5)
  end function profile

! Subroutines to evaluate desired spectral basis functions and derivatives
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

end program instanton
