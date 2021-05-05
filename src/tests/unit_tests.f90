
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

!===========================================================================================
! A module that runs all unit tests.
!
! Author: Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!===========================================================================================

module unit_tests

  use ftnunit, only: test
  use tests_ect
  use tests_inversion
  use tests_lsqr
  use tests_parallel_tools
  use tests_method_of_weights
  use tests_sparse_matrix

  implicit none

  private

  ! Runs main test suit.
  public :: test_all

contains

!==================================================================================
! Routine that simply runs all unit tests via the general routine "test".
!==================================================================================
subroutine test_all(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  if (nbproc == 1) then
  ! These tests are not adapted for the parallel version.

    call test(test_geometry_all, "Tests of the system geometry (ECT).", myrank, nbproc)
    call test(test_boundary_conditions_all, "Tests of the boundary conditions and right hand side (ECT).", myrank, nbproc)
    !call test(test_analytical_comparison_all, "Compare to analytical solution (ECT).", myrank, nbproc)

    call test(test_trans_mult_matrix, "Test of the matrices multiplication.", myrank, nbproc)
    call test(test_normalize_columns, "Test of the matrix columns normalization.", myrank, nbproc)
    call test(test_get_value, "Test of getting matrix element value.", myrank, nbproc)

    call test(test_method_of_weights_1, "Test of method of weights by Van Loan.", myrank, nbproc)
    call test(test_method_of_weights_2, "Test of method of weights by Per-Ake Wedin.", myrank, nbproc)
  endif

  ! Parallel tests:
  call test(test_get_number_elements_on_other_cpus, &
            "Test of getting number of elements on other cpus (parallel tools).", myrank, nbproc)
  call test(test_get_total_number_elements, "Test of getting total number of elements (parallel tools).", myrank, nbproc)
  call test(test_get_full_array, "Test of getting the full array (parallel tools).", myrank, nbproc)
  call test(test_get_full_array_in_place, "Test of getting the full array in place (parallel tools).", myrank, nbproc)

  call test(test_add_damping_identity_matrix, "Test of adding damping identity matrix (inversion).", myrank, nbproc)
  !call test(test_cross_gradient_calculate_all, "Tests of calculating cross gradient matrix (joint inversion).", myrank, nbproc)

  call test(test_lsqr_determined, "Test of the LSQR solver for a determined system.", myrank, nbproc)
  call test(test_lsqr_overdetermined_1, "Test 1 of the LSQR solver for an overdetermined system.", myrank, nbproc)
  call test(test_lsqr_overdetermined_2, "Test 2 of the LSQR solver for an overdetermined system.", myrank, nbproc)
  call test(test_lsqr_underdetermined_1, "Test 1 of the LSQR solver for an underdetermined system.", myrank, nbproc)
  call test(test_lsqr_underdetermined_2, "Test 2 of the LSQR solver for an underdetermined system.", myrank, nbproc)
  call test(test_lsqr_underdetermined_3, "Test 3 of the LSQR solver for an underdetermined system.", myrank, nbproc)

end subroutine test_all

end module unit_tests
