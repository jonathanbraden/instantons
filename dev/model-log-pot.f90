!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> University College London
!>
!>@brief
!> Store the equations of motion to be solved by our nonlinear boundary value solver
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Preprocessors for inlining.  Could also move this to another file
#define POTENTIAL(f) ( 0.5_dl*f**2*(1._dl-log(f**2+del**2)) )
#define VPRIME(f) ( -f*(log(f**2+del**2)-del**2) )
#define VDPRIME(f) ( -log(f**2+del**2) - 2._dl*f**2/(f**2+del**2) + del**2 )

module Model
  use constants
  use Cheby
  implicit none
  private :: L0, ndim, del
  
  real(dl), dimension(:,:), allocatable :: L0
  real(dl) :: del
  integer :: ndim

contains

  subroutine grid_params_(w,len,r0,w0)
    real(dl), intent(out) :: w, len
    real(dl), intent(in) :: r0, w0

    len = 2._dl; w = 1._dl
  end subroutine grid_params_

  subroutine bubble_parameters_nd_(delta,dim,r0,meff)
    real(dl), intent(in) :: delta, dim
    real(dl), intent(out) :: r0, meff

    meff = 1._dl
    r0 = 1._dl/meff
  end subroutine bubble_parameters_nd_

  subroutine get_minima(phif,phit)
    real(dl), intent(out) :: phif, phit
    phif = 0._dl; phit = 10._dl
  end subroutine get_minima

  subroutine initialise_equations(tForm, delta, dim)
    type(Chebyshev), intent(in) :: tForm
    real(dl), intent(in) :: delta
    integer, intent(in) :: dim
    integer :: i, sz
    
    sz = size(tForm%xGrid); ndim = dim
    if (allocated(L0)) deallocate(L0); allocate( L0(1:sz,1:sz) )
    do i=0,sz-1
       L0(i+1,:) = tForm%derivs(i,:,2) + dble(dim)*tForm%derivs(i,:,1)/tForm%xGrid(i)
    enddo
    del = delta
  end subroutine initialise_equations

  subroutine initialise_equations_scaled(tForm, delta, dim, scl)
    type(Chebyshev), intent(in) :: tForm
    real(dl), intent(in) :: delta, scl
    integer, intent(in) :: dim
    integer :: i, sz
    real(dl) :: d, m

    sz = size(tForm%xGrid)
    if (allocated(L0)) deallocate(L0)
    allocate( L0(1:sz,1:sz) )
    d = dble(dim) - 2._dl*scl
    m = scl*(scl+1._dl-dble(dim))
    do i=0,sz-1
       L0(i+1,:) = tForm%derivs(i,:,2) + d*tForm%derivs(i,:,1)/tForm%xGrid(i)
       L0(i+1,i+1) = L0(i+1,i+1) + m/tForm%xGrid(i)**2
    enddo
  end subroutine initialise_equations_scaled

  ! Move these to another module
  elemental function potential(phi)
    real(dl) :: potential
    real(dl), intent(in) :: phi
    potential = POTENTIAL(phi)
  end function potential

  elemental function vprime(phi)
    real(dl) :: vprime
    real(dl), intent(in) :: phi
    vprime =  VPRIME(phi)
  end function vprime

  elemental function vdprime(phi)
    real(dl) :: vdprime
    real(dl), intent(in) :: phi
    vdprime =  VDPRIME(phi)
  end function vdprime

  !!!!!!!!!!!!!!!!!!! This isn't well-defined for the Logarithmic potential
  !>@brief
  !> The unperturbed part of the potential used in computing the thin-wall domain wall
  elemental function potential_tw(phi)
    real(dl) :: potential_tw
    real(dl), intent(in) :: phi
    potential_tw = del*sin(phi)**2
  end function potential_tw
  
  subroutine source(fld,src)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:), intent(out) :: src
    integer :: sz

    sz = size(fld)
    src(:) = -matmul(L0,fld)
    src(:) = src(:) + VPRIME(fld(:))
    src(sz) = 0._dl  ! Set boundary condition at infinity
  end subroutine source

  !>@brief
  !> The second variation of our nonlinear equation
  !
  !> Notes on storage.  The first index is ? and the second index is ?
  !> The boundary conditions are designed to pick out phi(r_max),
  !> which is then set to 0 with the source term
  subroutine variation(fld,var)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:size(fld),1:size(fld)), intent(out) :: var
    integer :: i, sz

    sz = size(fld)
    var(1:sz,1:sz) = L0(1:sz,1:sz)
    do i=1,sz
       var(i,i) = var(i,i) - ( VDPRIME(fld(i)) ) 
    enddo
    ! boundary condition at infinity
    var(sz,:) = 0._dl
    var(sz,sz) = 1._dl
  end subroutine variation

  subroutine source_scaled(fld,src,scl,rvals)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:), intent(out) :: src
    real(dl), intent(in) :: scl
    real(dl), dimension(1:), intent(in) :: rvals
    real(dl), dimension(:), allocatable :: ftmp
    integer :: sz
    
    sz = size(fld); allocate( ftmp(sz) )
    ftmp = fld / rvals**scl
    src(:) = -matmul(L0,fld)
    src(:) = src(:) + rvals**scl*VPRIME(ftmp)
    deallocate(ftmp)
  end subroutine source_scaled

  subroutine variation_scaled(fld,var,scl,rvals)
    real(dl), dimension(1:), intent(in) :: fld
    real(dl), dimension(1:,1:), intent(out) :: var
    real(dl), intent(in) :: scl
    real(dl), dimension(1:), intent(in) :: rvals
    integer :: i, sz
    real(dl), dimension(:), allocatable :: ftmp

    sz = size(fld); allocate( ftmp(sz) )
    ftmp = fld / rvals**scl

    var(1:sz,1:sz) = L0(1:sz,1:sz)
    do i=1,sz
       ! Hmm, I think this multiplication by rvals on the V'' term is wrong.  Which means scaled and unscaled versions are the same
       var(i,i) = var(i,i) - rvals(i)**scl*( VDPRIME(ftmp(i)) ) !- (3._dl*fld(i)**2-1._dl +2._dl*del*fld(i))
    enddo
    ! boundary condition at infinity
    var(sz,:) = 0._dl
    var(sz,sz) = 1._dl
  end subroutine variation_scaled

  !>@brief
  !> A subroutine to set general Robin boundary conditions on our fields
  !>  \f[
  !>    \alpha_L f(x_L) + \beta_L f'(x_L) = c_L
  !>  \f]
  !>  \f[
  !>    \alpha_R f(x_R) + \beta_R f'(x_R) = c_R
  !>  \f]
  subroutine boundaries(L,S,c,bc)
    real(dl), intent(inout) :: L(1:,1:), S(1:)
    real(dl), dimension(1:3,1:2), intent(in) :: c
    logical, dimension(1:2), intent(in) :: bc
    integer :: sz
    sz = size(S)

    if (bc(1)) then
       L(:,:) = c(1,1) + c(2,1)
       S(1) = c(3,1)
    endif
    if (bc(2)) then
       L(:,:) = c(1,2) + c(2,2)
       S(sz) = c(3,2)
    endif
  end subroutine boundaries

end module Model
