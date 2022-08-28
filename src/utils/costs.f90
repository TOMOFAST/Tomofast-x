
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

!===========================================================================================
! Routines to calculate the cost functions (e.g. data misfit).
!
! Vitaliy Ogarko, UWA, CET, Australia.
!===========================================================================================

module costs

  use global_typedefs

  implicit none

  private

  public :: calculate_cost
  public :: calculate_cost_model

contains

!=====================================================================
! Calculates the relative cost.
! Flag 'in_parallel' for when arrays are split between the CPUs.
!=====================================================================
subroutine calculate_cost(arr1, arr2, cost, in_parallel, nbproc)
  integer, intent(in) :: nbproc
  real(kind=CUSTOM_REAL), intent(in) :: arr1(:)
  real(kind=CUSTOM_REAL), intent(in) :: arr2(:)
  logical, intent(in) :: in_parallel

  real(kind=CUSTOM_REAL), intent(out) :: cost

  real(kind=CUSTOM_REAL) :: cost1, cost2
  real(kind=CUSTOM_REAL) :: cost1_glob, cost2_glob
  integer :: ierr

  cost = 0._CUSTOM_REAL
  cost1 = 0._CUSTOM_REAL
  cost2 = 0._CUSTOM_REAL

  cost1 = sum((arr1 - arr2)**2)
  cost2 = sum((arr1)**2)

  if (in_parallel .and. nbproc > 1) then
    call MPI_Allreduce(cost1, cost1_glob, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(cost2, cost2_glob, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  else
    cost1_glob = cost1
    cost2_glob = cost2
  endif

  if (cost2_glob /= 0) cost = cost1_glob / cost2_glob

end subroutine calculate_cost

!==========================================================================================================
! Computes the model cost.
!==========================================================================================================
subroutine calculate_cost_model(nelements_total, norm_power, model, model_prior, column_weight, cost_model)
  integer, intent(in) :: nelements_total
  real(kind=CUSTOM_REAL), intent(in) :: norm_power
  real(kind=CUSTOM_REAL), intent(in) :: model(nelements_total)
  real(kind=CUSTOM_REAL), intent(in) :: model_prior(nelements_total)
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(nelements_total)

  real(kind=CUSTOM_REAL), intent(out) :: cost_model

  ! Local variables.
  real(kind=CUSTOM_REAL) :: model_diff
  integer :: i

  cost_model = 0._CUSTOM_REAL
  do i = 1, nelements_total
    ! Make calculations consistent with damping.F90 where we also scale the model difference.
    model_diff = (model(i) - model_prior(i)) / column_weight(i)
    cost_model = cost_model + (abs(model_diff))**norm_power
  enddo

end subroutine calculate_cost_model

end module costs


