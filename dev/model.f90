!>@todo
!>@arg Factor out the setting of boundary conditions into a separate part of the code
module Model
  use constants
  implicit none
  
  integer, parameter :: nfield = 1
  real(dl), parameter :: lambda = 2.
  
  integer, parameter :: ndim = 3

! These hardcoded parameters should be replaced by a computational subroutine
  real(dl), parameter :: tension = 2._dl**1.5
  real(dl), parameter :: delrho = 2._dl
  real(dl), parameter :: rinit = dble(ndim)*tension / delrho

  ! Temporary storage variables for derivatives
  real(dl), dimension(:,:), allocatable :: L0, L2

  type Model
     type(Chebyshev) :: tForm
     real(dl) :: lambda, delta
     real(dl), dimension(:), allocatable :: params
     real(dl) :: tension, drho, phit, phif
     real(dl) :: r0, width
  end type Model

contains
  ! Think more carefully about how to do this
  subroutine initialise_model(this,trans)
    type(Model), intent(inout) :: this
    type(Chebyshev), intent(in) :: trans

    this%tForm = trans ! Is this a pointer or new variable?
    
  end subroutine initialise_model

  !>@brief
  !> Fill the source and linear terms in spectral space
  subroutine eval_func_spec(f,L,S)
    real(dl), dimension(1:), intent(in) :: f
    real(dl), dimension(1:,1:), intent(out) :: L
    real(dl), dimension(1:), intent(out) :: S
  end subroutine eval_func_spec
  
  !>@brief
  !> Fill the source term and linear term for our differential equation
  !>  \f[
  !>    L[\phi]\delta\phi = S[\phi]
  !>  \f]
  subroutine eval_func(f,L,S)
    real(dl), dimension(1:), intent(in) :: f
    real(dl), dimension(1:,1:), intent(out) :: L
    real(dl), dimension(1:), intent(out) :: S

    real(dl), dimension(:), allocatable :: tmpS
    integer :: i, n

    n = size(f)
    allocate( tmpS(1:n) )
    
    L = L0  ! need to set L0 initially
    do i=1,
       L(i,i) = L(i,i) - vdprime(f(i))
    enddo

    tmpS = vprime(f)
    S = -matmul(l0,f) + tmpS

    ! Dirichlet BC
    L(n,:) = 0._dl
    L(n,n) = 1._dl
    S(n) = 0._dl
    
  end subroutine eval_func

  !>@brief
  !> Set the boundary conditions in my Newton solver
  subroutine set_bc(L,S)
    real(dl), dimension(:,:), intent(inout) :: L
    real(dl), dimension(:), intent(inout) :: S
  end subroutine set_bc
  
  !>@brief
  !> Compute the source term for my iterative solver.
  !> This is simply the rights hand side of the first order set of equations
  !>  \f[
  !>    \vec{S}[\vec{f}] = \frac{d\vec{f}}{dt}
  !>  \f]
  subroutine source(f,S)
    real(dl), dimension(0:), intent(in) :: f
    real(dl), dimension(0:), intent(out) :: S

    S(:) = -matmul(l0,f) + vprime(f)
  end subroutine source

  !>@brief
  !> Compute the linear part of the perturbed differential equation
  !> (ie. first variation of the differential equation)
  subroutine linear(f,L)
    real(dl), dimension(:), intent(in) :: f
    real(dl), dimension(:,:), intent(out) :: L
  end subroutine linear

  elemental function potential(phi)
    real(dl), intent(in) :: phi
    real(dl) :: potential
    potential = cos(phi) + lambda*sin(phi)**2
  end function potential

  elemental function vprime(phi)
    real(dl), intent(in) :: phi
    real(dl) :: vprime
    vprime = -sin(phi) + lambda*sin(2._dl*phi)
  end function vprime

  elemental function vdprime(phi)
    real(dl), intent(in) :: phi
    real(dl) :: vdprime
    vdprime = -cos(phi) + 2._dl*lambda*cos(2._dl*phi)
  end function vdprime

  subroutine get_vacuum(fld)
    real(dl), intent(inout) :: fld

    integer, parameter :: maxit=16
    real(dl), parameter :: min_tol=1.e-14
    real(dl) :: vp,vpp,dfld

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
       print*,"Failed to find local minimum of potential. Adjust guess"
       stop
    endif
    print*,"Vacuum is ",fld," derivative is ",vprime(fld)
  end subroutine get_vacuum

  subroutine initial_profile()

  end subroutine initial_profile

end module Model
