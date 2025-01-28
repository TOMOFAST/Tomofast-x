
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
subroutine calculate_cost(n, arr1, arr2, cost, in_parallel, nbproc)
  integer, intent(in) :: n, nbproc
  real(kind=CUSTOM_REAL), intent(in) :: arr1(n)
  real(kind=CUSTOM_REAL), intent(in) :: arr2(n)
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
    call mpi_allreduce(cost1, cost1_glob, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    call mpi_allreduce(cost2, cost2_glob, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  else
    cost1_glob = cost1
    cost2_glob = cost2
  endif

  if (cost2_glob /= 0) then
    cost = sqrt(cost1_glob / cost2_glob)
  endif

end subroutine calculate_cost

!===============================================================================================
! Computes the Lp norm of the difference of two models, that are split between CPUs.
!===============================================================================================
subroutine calculate_cost_model(nelements, norm_power, model, model_prior, column_weight, &
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
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  real(kind=CUSTOM_REAL), intent(out) :: cost_model

  ! Local variables.
  real(kind=CUSTOM_REAL) :: model_diff
  real(kind=CUSTOM_REAL) :: cost_model_glob
  integer :: i, ierr

  ! Compute local cost model.
  cost_model = 0._CUSTOM_REAL
  do i = 1, nelements
    ! Make calculations consistent with damping.F90 where we also scale the model difference.
    if (column_weight(i) /= 0.d0) then
      model_diff = (model(i) - model_prior(i)) / column_weight(i)
    else
      model_diff = 0.d0
    endif
    cost_model = cost_model + (abs(model_diff))**norm_power
  enddo

  ! Compute global cost model.
  cost_model_glob = 0._CUSTOM_REAL
  if (nbproc > 1) then
    call mpi_allreduce(cost_model, cost_model_glob, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    cost_model = cost_model_glob
  endif

end subroutine calculate_cost_model

end module costs


