program instanton
  use constants
  use chebyshev
  use field_model
  implicit none

  type(Field_Model) :: model
  type(Chebyshev) :: pspec
  real(dl), dimension(:), allocatable :: instanton, instanton_i
  real(dl), dimension(:), allocatable :: model_params

  integer :: order
  real(dl) :: w, len  ! parameters controlling mapping of collocation grid
  integer :: i
  
  order = 100
  w = 1._dl; len = 20. 
  ! Initialise our derivatives
  call initialise_chebyshev(pspec)
  call create_model(model)

  call initialise_field_analytic()
  do i=1,nparam
     call initialise_field_numerical()
     call nonlinear_solve()
     call output()
  enddo
  
contains

  !>@brief
  !> Initialise our initial guess for the instanton profile
  !>
  !>@param[in] prev  If True - initialise profile using previous numerical profile
  !>                 If False - initialise using an analytic formula
  subroutine initialise_fields(prev)
    logical, intent(in) :: prev

    if (prev) then
       
    else
       call thin_wall_profile()
    endif
  end subroutine initialise_fields

  !>@brief
  !> Given the previous instanton profile, adjust the coordinate mapping to
  subroutine adapt_grid()
  end subroutine adapt_grid

  !>@brief
  !> Interpolate our previous instanton solution to the new collocation grid
  subroutine interpolate_to_new_grid()
  end subroutine interpolate_to_new_grid

  subroutine get_evalues()
  end subroutine get_evalues
  
  subroutine get_evectors()
  end subroutine get_evectors
  
  subroutine get_action()
  end subroutine get_action
  
end program instanton
