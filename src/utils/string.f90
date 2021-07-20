
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

!==============================================================================================================
! Contains various helper functions to work with strings.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015.
!==============================================================================================================
module string

  use global_typedefs, only: CUSTOM_REAL

  implicit none

  private

  public :: str
  public :: str2int
  public :: str2real
  public :: extract_number
  public :: extract_int_number

contains

!============================================
! Convert an integer to string.
!============================================
character(len=20) function str(k)
  integer, intent(in) :: k

  write (str, *) k
  str = adjustl(str)
end function str

!============================================
! Convert a string to integer.
!============================================
function str2int(str) result(res)
  character(len=*), intent(in) :: str
  integer :: res

  read(str, *) res
end function str2int

!============================================
! Convert a string to real.
!============================================
function str2real(str) result(res)
  character(len=*), intent(in) :: str
  real(kind=CUSTOM_REAL) :: res

  read(str, *) res
end function str2real

!=================================================
! Extracts a real number from the string.
!=================================================
function extract_number(str, symbol) result(res)
  character(len=*), intent(in) :: str
  character(len=1), intent(in) :: symbol
  real(kind=CUSTOM_REAL) :: res
  integer :: where

  where = index(str, symbol)
  res = str2real(str(where + 1:))

end function extract_number

!=================================================
! Extracts an integer number from the string.
!=================================================
function extract_int_number(str, symbol) result(res)
  character(len=*), intent(in) :: str
  character(len=1), intent(in) :: symbol
  real(kind=CUSTOM_REAL) :: res
  integer :: where

  where = index(str, symbol)
  res = str2int(str(where + 1:))

end function extract_int_number

end module string
