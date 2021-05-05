
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

!==========================================================================
! A class to work with 3D vectors.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015.
!==========================================================================
module vector

  implicit none

  private

  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: assignment(=)

  integer, parameter :: CUSTOM_REAL = 8

  !-----------------------------------------------
  ! A main class (floating point vectors).
  type, public :: t_vector
    real(kind=CUSTOM_REAL) :: x, y, z

  contains
    private

    procedure, public, pass :: cross_product => vector_cross_product
    procedure, public, pass :: dot_product => vector_dot_product
    procedure, public, pass :: get_norm => vector_get_norm
  end type t_vector

  interface t_vector
    module procedure vector_constructor
  end interface t_vector

  interface operator(+)
    module procedure vector_add
  end interface

  interface operator(-)
    module procedure vector_subtract
  end interface

  interface operator(*)
    module procedure vector_mult
  end interface

  interface assignment(=)
    module procedure vector_assign
  end interface

  !-----------------------------------
  ! A class for integer vectors.
  type, public :: t_ivector
    integer :: x, y, z
  end type t_ivector

  interface t_ivector
    module procedure ivector_constructor
  end interface t_ivector

contains

!=======================================================================
! Constructor for t_ivector type.
!=======================================================================
function ivector_constructor(i, j, k) result(res)
  integer, intent(in) :: i, j, k
  type(t_ivector) :: res

  res%x = i
  res%y = j
  res%z = k

end function ivector_constructor

!=======================================================================
! Constructor for t_vector type.
!=======================================================================
function vector_constructor(x, y, z) result(res)
  real(kind=CUSTOM_REAL) :: x, y, z
  type(t_vector) :: res

  res%x = x
  res%y = y
  res%z = z

end function vector_constructor

!=======================================================================
! Returns cross-product between vectors.
!=======================================================================
pure function vector_cross_product(this, vec) result(res)
  class(t_vector), intent(in) :: this
  type(t_vector), intent(in) :: vec
  type(t_vector) :: res

  res%x = this%y * vec%z - this%z * vec%y
  res%y = this%z * vec%x - this%x * vec%z
  res%z = this%x * vec%y - this%y * vec%x

end function vector_cross_product

!=======================================================================
! Returns dot-product between vectors.
!=======================================================================
pure function vector_dot_product(this, vec) result(res)
  class(t_vector), intent(in) :: this
  type(t_vector), intent(in) :: vec
  real(kind=CUSTOM_REAL) :: res

  res = this%x * vec%x + this%y * vec%y + this%z * vec%z

end function vector_dot_product

!=======================================================================
! Returns the norm of a vector.
!=======================================================================
pure function vector_get_norm(this) result(res)
  class(t_vector), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = sqrt(this%x**2 + this%y**2 + this%z**2)

end function vector_get_norm

!=======================================================================
! Returns the sum of two vectors.
!=======================================================================
pure function vector_add(v1, v2) result(res)
  type(t_vector), intent(in) :: v1, v2
  type(t_vector) :: res

  res%x = v1%x + v2%x
  res%y = v1%y + v2%y
  res%z = v1%z + v2%z

end function vector_add

!=======================================================================
! Returns the subtraction of two vectors.
!=======================================================================
pure function vector_subtract(v1, v2) result(res)
  type(t_vector), intent(in) :: v1, v2
  type(t_vector) :: res

  res%x = v1%x - v2%x
  res%y = v1%y - v2%y
  res%z = v1%z - v2%z

end function vector_subtract

!=======================================================================
! Assign a vector to a scalar.
!=======================================================================
subroutine vector_assign(lhs, rhs)
  type(t_vector), intent(out) :: lhs
  real(kind=CUSTOM_REAL), intent(in) :: rhs

  lhs%x = rhs
  lhs%y = rhs
  lhs%z = rhs

end subroutine vector_assign

!=======================================================================
! Returns a scalar times vector.
!=======================================================================
pure function vector_mult(const, vec) result(res)
  real(kind=CUSTOM_REAL), intent(in) :: const
  type(t_vector), intent(in) :: vec
  type(t_vector) :: res

  res%x = const * vec%x
  res%y = const * vec%y
  res%z = const * vec%z

end function vector_mult

end module vector
