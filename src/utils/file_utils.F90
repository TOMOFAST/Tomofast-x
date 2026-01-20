
!========================================================================
!
!                          T o m o f a s t - x
!                        -----------------------
!
!           Authors: Vitaliy Ogarko, Jeremie Giraud, Roland Martin.
!
!               (c) 2021 The University of Western Australia.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

!========================================================================
! Functions work with files.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!========================================================================
module file_utils
  implicit none

  private
  public :: create_directory

contains

!========================================================================
! Creates a folder, works on both Linux and Windows.
!========================================================================
subroutine create_directory(path)
  character(len=*), intent(in) :: path
  integer :: istat

#ifdef WINDOWS
  call execute_command_line('mkdir "'//trim(path)//'" 2>nul', exitstat=istat)
#else
  call execute_command_line('mkdir -p "'//trim(path)//'" 2>/dev/null', exitstat=istat)
#endif
end subroutine create_directory

end module file_utils
