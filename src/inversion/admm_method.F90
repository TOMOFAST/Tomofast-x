
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
subroutine admm_method_iterate_admm_arrays(this, nlithos, xmin, xmax, x, x0)
  class(t_admm_method), intent(inout) :: this
  integer, intent(in) :: nlithos
  real(kind=CUSTOM_REAL), intent(in) :: xmin(nlithos, this%nelements)
  real(kind=CUSTOM_REAL), intent(in) :: xmax(nlithos, this%nelements)
  real(kind=CUSTOM_REAL), intent(in) :: x(this%nelements)

  real(kind=CUSTOM_REAL), intent(out) :: x0(this%nelements)

  real(kind=CUSTOM_REAL) :: arg, mindist, val, closest_boundary
  integer :: i, j
  logical :: inside

  ! Calculate z[k + 1] = Pc(x[k + 1] + u[k]).
  do i = 1, this%nelements

    ! Calculate the indicator function.
    arg = x(i) + this%u(i)

    !if (arg < xmin(i)) then
    !  this%z(i) = xmin(i)
    !
    !else if (arg > xmax(i)) then
    !  this%z(i) = xmax(i)
    !
    !else
    !  this%z(i) = arg
    !endif

    inside = .false.
    do j = 1, nlithos
      ! Check if the value lies inside the bounds.
      if (xmin(j, i) <= arg .and. arg <= xmax(j, i)) then
        inside = .true.
        this%z(i) = arg
        exit
      endif
    enddo

    if (.not. inside) then
    ! The value lies outside boundaries, so finding the closest boundary.
      mindist = 1.d30
      do j = 1, nlithos
        val = dabs(xmin(j, i) - arg)
        if (val < mindist) then
          mindist = val
          closest_boundary = xmin(j, i)
        endif

        val = dabs(xmax(j, i) - arg)
        if (val < mindist) then
          mindist = val
          closest_boundary = xmax(j, i)
        endif
      enddo
      this%z(i) = closest_boundary
    endif
  enddo

  ! Calculate u[k + 1] = u[k] + x[k + 1] - z[k + 1].
  this%u = this%u + x - this%z

  x0 = this%z - this%u

end subroutine admm_method_iterate_admm_arrays

end module admm_method
