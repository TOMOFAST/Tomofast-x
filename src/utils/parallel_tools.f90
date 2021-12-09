
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

!================================================================================================
! Static functions (i.e., no object passed) to work with parallel arrays and scalars.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!================================================================================================
module parallel_tools

  use global_typedefs
  use mpi_tools, only: exit_MPI

  implicit none

  private

  type, public :: t_parallel_tools
    private

  contains
    private

    procedure, public, nopass :: calculate_nelements_at_cpu
    procedure, public, nopass :: get_nsmaller
    procedure, public, nopass :: get_number_elements_on_other_cpus
    procedure, public, nopass :: get_total_number_elements
    procedure, public, nopass :: get_mpi_partitioning
    procedure, public, nopass :: get_full_array
    procedure, public, nopass :: get_full_array_in_place
    procedure, public, nopass :: get_full_array_in_place2

  end type t_parallel_tools

contains

!==================================================================================
! Define model splitting for parallelization:
! calculate local number of elements (at current CPU).
!==================================================================================
function calculate_nelements_at_cpu(nelements_total, myrank, nbproc) result(nelements)
  integer, intent(in) :: nelements_total, myrank, nbproc
  integer :: nelements

  nelements = nelements_total / nbproc

  if (myrank == nbproc - 1 .and. mod(nelements_total, nbproc) /= 0) then
    ! Last rank gets the remaining elements.
    nelements = nelements + mod(nelements_total, nbproc)
  endif
end function calculate_nelements_at_cpu

!=====================================================================================
! Calculates the number of elements on CPUs with rank smaller than myrank.
!=====================================================================================
function get_nsmaller(nelements, myrank, nbproc) result(nsmaller)
  integer, intent(in) :: nelements
  integer, intent(in) :: myrank, nbproc

  integer :: i
  ! The number of elements on every CPU.
  integer :: nelements_at_cpu(nbproc)
  integer :: nsmaller

  ! Calculate the number of elements on every CPU.
  nelements_at_cpu = get_number_elements_on_other_cpus(nelements, myrank, nbproc)

  ! Note that CPU_index = myrank + 1.
  nsmaller = 0
  do i = 1, myrank
    nsmaller = nsmaller + nelements_at_cpu(i)
  enddo

end function get_nsmaller

!======================================================================================================
! Gather the nelements from other CPUs.
!======================================================================================================
function get_number_elements_on_other_cpus(nelements, myrank, nbproc) result(nelements_at_cpu)
  integer, intent(in) :: nelements, myrank, nbproc
  integer :: ierr
  integer :: nelements_at_cpu(nbproc)

  call MPI_Allgather(nelements, 1, MPI_INTEGER, nelements_at_cpu, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI error in get_number_elements_on_other_cpus!", myrank, ierr)

end function get_number_elements_on_other_cpus

!======================================================================================================
! Returns the total number of model elements.
!======================================================================================================
function get_total_number_elements(nelements, myrank, nbproc) result(nelements_total)
  integer, intent(in) :: nelements, myrank, nbproc

  integer :: nelements_at_cpu(nbproc)
  integer :: nelements_total

  ! Calculate the number of elements on every CPU.
  nelements_at_cpu = get_number_elements_on_other_cpus(nelements, myrank, nbproc)

  nelements_total = sum(nelements_at_cpu)
end function get_total_number_elements

!=========================================================================================
! Get index arrays (displacement and number of elements) for mpi_scatterv / mpi_gatherv.
!=========================================================================================
subroutine get_mpi_partitioning(nelements, displs, nelements_at_cpu, myrank, nbproc)
  ! Parameters.
  integer, intent(in) :: nelements
  integer, intent(in) :: myrank, nbproc

  ! Displacement for mpi_scatterv / mpi_gatherv.
  integer, intent(out) :: displs(:)
  ! The number of elements on every CPU for mpi_scatterv / mpi_gatherv.
  integer, intent(out) :: nelements_at_cpu(:)
  integer :: i

  nelements_at_cpu = get_number_elements_on_other_cpus(nelements, myrank, nbproc)

  ! Note, start indexing from 1.
  displs(1) = 0
  do i = 2, nbproc
    displs(i) = displs(i - 1) + nelements_at_cpu(i - 1)
  enddo

end subroutine get_mpi_partitioning

!======================================================================================================
! Returns the full array (that is split between CPUs)
! bcast = false: to only master CPU.
! bcast = true:  to all CPUs (using Bcast).
! This function is unit tested in tests_inversion.f90 in test_get_full_array().
!======================================================================================================
subroutine get_full_array(val, nelements, array, bcast, myrank, nbproc)
  real(kind=CUSTOM_REAL), intent(in) :: val(:)
  integer, intent(in) :: nelements
  logical, intent(in) :: bcast
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(out) :: array(:)

  ! Displacement for MPI_Gatherv.
  integer :: displs(nbproc)
  ! The number of elements on every CPU for MPI_Gatherv.
  integer :: nelements_at_cpu(nbproc)
  integer :: nelements_total
  integer :: ierr

  if (nbproc == 1) then
    array = val

  else
    ! Get partitioning for MPI_Gatherv.
    call get_mpi_partitioning(nelements, displs, nelements_at_cpu, myrank, nbproc)

    ! Gather the full model vector.
    call MPI_Gatherv(val, nelements, CUSTOM_MPI_TYPE, &
                     array, nelements_at_cpu, displs, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

    if (bcast) then
    ! Cast the array to all CPUs.
      nelements_total = sum(nelements_at_cpu)

      call MPI_Bcast(array, nelements_total, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
    endif

    if (ierr /= 0) call exit_MPI("MPI error in get_full_array!", myrank, ierr)
  endif

end subroutine get_full_array

!======================================================================================================
! Same as get_full_array() but uses the same send and recv buffer in MPI_Gatherv.
! This function is unit tested in tests_inversion.f90 in test_get_full_array_in_place().
!======================================================================================================
subroutine get_full_array_in_place(nelements, array, bcast, myrank, nbproc)
  integer, intent(in) :: nelements
  logical, intent(in) :: bcast
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(inout) :: array(:)

  ! Displacement for MPI_Gatherv.
  integer :: displs(nbproc)
  ! The number of elements on every CPU for MPI_Gatherv.
  integer :: nelements_at_cpu(nbproc)
  integer :: nelements_total
  integer :: ierr

  ! Get partitioning for MPI_Gatherv.
  call get_mpi_partitioning(nelements, displs, nelements_at_cpu, myrank, nbproc)

  if (myrank == 0) then
    call MPI_Gatherv(MPI_IN_PLACE, nelements, CUSTOM_MPI_TYPE, &
                     array, nelements_at_cpu, displs, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  else
    call MPI_Gatherv(array, nelements, CUSTOM_MPI_TYPE, &
                     array, nelements_at_cpu, displs, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  endif

  if (bcast) then
  ! Cast the array to all CPUs.
    nelements_total = sum(nelements_at_cpu)

    call MPI_Bcast(array, nelements_total, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  endif

  if (ierr /= 0) call exit_MPI("MPI error in get_full_array_in_place!", myrank, ierr)

end subroutine get_full_array_in_place

!======================================================================================================
! Same as get_full_array_in_place() but uses MPI_Allgatherv to combine Gather with Bcast.
!======================================================================================================
subroutine get_full_array_in_place2(nelements, array, myrank, nbproc)
  integer, intent(in) :: nelements
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(inout) :: array(:)

  ! Displacement for MPI_Gatherv.
  integer :: displs(nbproc)
  ! The number of elements on every CPU for MPI_Gatherv.
  integer :: nelements_at_cpu(nbproc)
  integer :: ierr

  ! Get partitioning for MPI_Gatherv.
  call get_mpi_partitioning(nelements, displs, nelements_at_cpu, myrank, nbproc)

  call MPI_Allgatherv(MPI_IN_PLACE, nelements, CUSTOM_MPI_TYPE, &
                      array, nelements_at_cpu, displs, CUSTOM_MPI_TYPE, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI error in get_full_array_in_place2!", myrank, ierr)

end subroutine get_full_array_in_place2

end module parallel_tools
