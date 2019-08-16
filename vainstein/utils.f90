!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! MODULE: Utils
!
!>@author
!> Jonathan Braden, University College London
!>
!>@brief
!> A collection of useful utility functions for Fortran
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Utils
  implicit none
  
contains

  !>@brief
  !> Search for an available output unit
  integer function newunit(unit)
    integer, intent(out), optional :: unit

    integer, parameter :: LUN_MIN=10, LUN_MAX=1000
    logical :: o; integer :: lun
    newunit = -1
    do lun = LUN_MIN, LUN_MAX
       inquire(unit=lun,opened=o)
       if (.not.o) then
          newunit = lun
          exit
       endif
    enddo
    if (present(unit)) unit=newunit
  end function newunit
  
end module Utils
