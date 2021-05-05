
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


module forward_problem_gravmag

  use global_typedefs, only: CUSTOM_REAL
  use data_gravmag
  use sensitivity_gravmag
  use parameters_gravmag
  use inversion_arrays

  implicit none

  private

  public :: solve_forward_problem

contains

!=========================================================================================
! Solves forward problem for gravity or magnetism.
!=========================================================================================
subroutine solve_forward_problem(par, iarr, data, myrank)
  class(t_parameters_base), intent(in) :: par
  type(t_data), intent(inout) :: data
  type(t_inversion_arrays), intent(inout) :: iarr
  integer, intent(in) :: myrank

  type(t_sensitivity_gravmag) :: sens

  if (par%calc_data_directly == 1) then
    ! Calculate data directly without computing sensitivity matrix.
    call sens%calc_data_directly(par, iarr%model, data, myrank)

    return
  endif

  ! Calculate sensitivity kernel (analytically).
  call sens%calculate_sensitivity(par, iarr%model%grid, data, iarr%matrix_sensit, par%distance_threshold, myrank)

  ! Calculate the data using sensitivity (S) and prior model (m) as d = S * m.
  call iarr%model%calculate_data(par%ndata, iarr%matrix_sensit, data%val_calc, myrank)

end subroutine solve_forward_problem

end module forward_problem_gravmag
