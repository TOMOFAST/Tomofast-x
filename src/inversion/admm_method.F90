
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
! A class for adding the ADMM method constraints.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!===========================================================================================
module admm_method

  use global_typedefs
  use mpi_tools, only: exit_MPI

  implicit none

  private

  type, public :: t_admm_method
    integer :: nelements

    real(kind=CUSTOM_REAL), allocatable :: z(:)
    real(kind=CUSTOM_REAL), allocatable :: u(:)

  contains
    private

    procedure, public, pass :: initialize => admm_method_initialize
    procedure, public, pass :: iterate_admm_arrays => admm_method_iterate_admm_arrays

  end type t_admm_method

contains

!===========================================================================================
! Initialization.
!===========================================================================================
subroutine admm_method_initialize(this, nelements, myrank)
  class(t_admm_method), intent(inout) :: this
  integer, intent(in) :: nelements
  integer, intent(in) :: myrank

  integer :: ierr

  this%nelements = nelements

  ierr = 0

  ! Set initial values for u[0] and z[0] to zero.
  allocate(this%z(nelements), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%u(nelements), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in admm_method_initialize!", myrank, ierr)

end subroutine admm_method_initialize

!===========================================================================================
! Calculating the main ADMM arrays needed to add the ADMM constraints.
!===========================================================================================
subroutine admm_method_iterate_admm_arrays(this, nlithos, xmin, xmax, x, x0, myrank)
  class(t_admm_method), intent(inout) :: this
  integer, intent(in) :: nlithos
  real(kind=CUSTOM_REAL), intent(in) :: xmin(nlithos, this%nelements)
  real(kind=CUSTOM_REAL), intent(in) :: xmax(nlithos, this%nelements)
  real(kind=CUSTOM_REAL), intent(in) :: x(this%nelements)
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: x0(this%nelements)

  real(kind=CUSTOM_REAL) :: arg, mindist, val, closest_boundary
  real(kind=CUSTOM_REAL) :: lambda
  integer :: i, j
  logical :: inside

  if (myrank == 0) print *, 'Calculating the ADMM arrays.'

  lambda = 1.d-2

  ! Calculate z[k + 1] = Pc(x[k + 1] + u[k]).
  do i = 1, this%nelements

    ! Calculate the indicator function.
    arg = x(i) + this%u(i)

    !-----------------------------------------
    ! Test 1.
    if (abs(arg) <= lambda) then
      this%z(i) = 0.d0

    else if (arg > lambda) then
      this%z(i) = arg - lambda

    else
      this%z(i) = arg + lambda
    endif

  enddo

  ! Calculate u[k + 1] = u[k] + x[k + 1] - z[k + 1].
  this%u = this%u + x - this%z

  x0 = this%z - this%u

end subroutine admm_method_iterate_admm_arrays

end module admm_method
