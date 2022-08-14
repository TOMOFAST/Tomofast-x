
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

!===============================================================================================
! Unit tests for parallel (MPI) tools.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!===============================================================================================
module tests_parallel_tools

  use global_typedefs
  use ftnunit

  use parallel_tools

  implicit none

  private

  public :: test_get_number_elements_on_other_cpus
  public :: test_get_total_number_elements
  public :: test_get_full_array
  public :: test_get_full_array_in_place

contains

!=============================================================================================
! Testing get_number_elements_on_other_cpus().
!=============================================================================================
subroutine test_get_number_elements_on_other_cpus(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  integer :: nelements_at_cpu(nbproc)
  integer :: i, nelements

  nelements = 5 * (myrank + 1) + 3

  nelements_at_cpu = get_number_elements_on_other_cpus(nelements, myrank, nbproc)

  do i = 1, nbproc
    call assert_equal_int(nelements_at_cpu(i), 5 * i + 3, "test_get_number_elements_on_other_cpus failed.")
  enddo

end subroutine test_get_number_elements_on_other_cpus

!=============================================================================================
! Testing get_total_number_elements().
!=============================================================================================
subroutine test_get_total_number_elements(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  integer :: nelements, nelements_total

  ! Set different number of elements on every CPU, as:
  ! CPU 1: 1 element,
  ! CPU 2: 2 elements,
  ! CPU 3: 3 elements, etc.
  nelements = myrank + 1

  nelements_total = get_total_number_elements(nelements, myrank, nbproc)

  call assert_equal_int(nelements_total, (nbproc + 1) * nbproc / 2, "Wrong nelements_total in test_get_total_number_elements.")

end subroutine test_get_total_number_elements

!=============================================================================================
! Testing get_full_array().
!=============================================================================================
subroutine test_get_full_array(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  integer :: nelements, nelements_total
  real(kind=CUSTOM_REAL), allocatable :: model(:)
  real(kind=CUSTOM_REAL), allocatable :: model_all(:)
  integer :: i

  ! Set different number of elements on every CPU, as:
  ! CPU 1: 1 element,
  ! CPU 2: 2 elements,
  ! CPU 3: 3 elements, etc.
  nelements = myrank + 1

  nelements_total = get_total_number_elements(nelements, myrank, nbproc)

  allocate(model(nelements))
  allocate(model_all(nelements_total))

  ! Form a model vector (1, 2, 3, 4, ...).
  ! To do this, sum nelements on the lower ranks,
  ! assuming the above relation nelements(myrank) = myrank + 1.
  do i = 1, nelements
    model(i) = dble((myrank + 1) * myrank / 2 + i)
  enddo

  call get_full_array(model, nelements, model_all, .true., myrank, nbproc)

  ! Check the result.
  do i = 1, nelements_total
    call assert_equal_int(int(model_all(i)), i, "model_all(i) /= i in test_get_full_array.")
  enddo

  deallocate(model)
  deallocate(model_all)

end subroutine test_get_full_array

!=============================================================================================
! Testing get_full_array_in_place().
!=============================================================================================
subroutine test_get_full_array_in_place(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  integer :: nelements, nelements_total
  real(kind=CUSTOM_REAL), allocatable :: model_all(:)
  integer :: i

  ! Set different number of elements on every CPU, as:
  ! CPU 1: 1 element,
  ! CPU 2: 2 elements,
  ! CPU 3: 3 elements, etc.
  nelements = myrank + 1

  nelements_total = get_total_number_elements(nelements, myrank, nbproc)

  allocate(model_all(nelements_total))

  ! Form a model vector (1, 2, 3, 4, ...).
  ! To do this, sum nelements on the lower ranks,
  ! assuming the above relation nelements(myrank) = myrank + 1.
  do i = 1, nelements
    model_all(i) = dble((myrank + 1) * myrank / 2 + i)
  enddo

  call get_full_array_in_place(nelements, model_all, .true., myrank, nbproc)

  ! Check the result.
  do i = 1, nelements_total
    call assert_equal_int(int(model_all(i)), i, "model_all(i) /= i in test_get_full_array_in_place.")
  enddo

  deallocate(model_all)

end subroutine test_get_full_array_in_place

end module tests_parallel_tools
