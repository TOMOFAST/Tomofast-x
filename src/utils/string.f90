
!========================================================================
!
!                    T O M O F A S T X  Version 1.0
!                  ----------------------------------
!
!              Main authors: Vitaliy Ogarko, Roland Martin,
!                   Jeremie Giraud, Dimitri Komatitsch.
! CNRS, France, and University of Western Australia.
! (c) CNRS, France, and University of Western Australia. January 2018
!
! This software is a computer program whose purpose is to perform
! capacitance, gravity, magnetic, or joint gravity and magnetic tomography.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
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
