
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
! Routines to calculate the cost functions (e.g. data misfit).
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!===========================================================================================

module costs

  use global_typedefs

  implicit none

  private

  public :: calculate_cost
  public :: calculate_cost_model

contains

!=====================================================================
! Calculates and prints out the cost functions.
!=====================================================================
subroutine calculate_cost(n, arr1, arr2, cost, myrank)

  integer, intent(in) :: n, myrank
  real(kind=CUSTOM_REAL), intent(in) :: arr1(:)
  real(kind=CUSTOM_REAL), intent(in) :: arr2(:)

  real(kind=CUSTOM_REAL), intent(out) :: cost

  real(kind=CUSTOM_REAL) :: cost1, cost2
  integer :: i

  cost = 0._CUSTOM_REAL
  cost1 = 0._CUSTOM_REAL
  cost2 = 0._CUSTOM_REAL

  do i = 1, n
    cost1 = cost1 + (arr1(i) - arr2(i))**2
    cost2 = cost2 + arr1(i)**2
  enddo

  if (cost2 /= 0) cost = cost1 / cost2

  if (myrank == 0) print *, 'cost =', cost, 'cost1 =', cost1, 'cost2 =', cost2

end subroutine calculate_cost

!===============================================================================================
! Computes the Lp norm of the difference of two models, that are split between CPUs.
!===============================================================================================
subroutine calculate_cost_model(nelements, norm_power, model, model_prior, damping_weight, &
                                cost_model, nbproc)
  ! MPI variables.
  integer, intent(in) :: nbproc
  ! Local dimension of the model vector.
  integer, intent(in) :: nelements

  ! Power of the Lp norm on the model variation x-xprior : Int_V |x-xprior|^p dV.
  real(kind=CUSTOM_REAL), intent(in) :: norm_power
  real(kind=CUSTOM_REAL), intent(in) :: model(:)
  real(kind=CUSTOM_REAL), intent(in) :: model_prior(:)
  ! Inversion weights to scale damping.
  real(kind=CUSTOM_REAL), intent(in) :: damping_weight(:)
  real(kind=CUSTOM_REAL), intent(out) :: cost_model

  ! Local variables.
  real(kind=CUSTOM_REAL) :: cost_model_glob
  integer :: i, ierr

  ! Compute local cost model.
  cost_model = 0._CUSTOM_REAL
  do i = 1, nelements
    ! TODO: Do we need to use norm_power for damping_weight too? Then also need to do this in damping.F90,
    !       and change the cell volume weight in capacitance.
    cost_model = cost_model + abs(model(i) - model_prior(i))**norm_power &
                 * (damping_weight(i))**2
  enddo

  ! Compute global cost model.
  cost_model_glob = 0._CUSTOM_REAL
  if (nbproc > 1) then
    call mpi_allreduce(cost_model, cost_model_glob, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    cost_model = cost_model_glob
  endif

end subroutine calculate_cost_model

end module costs


