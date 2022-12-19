
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

!==========================================================================================
! A class that stores allocatable arrays needed for inversion.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!==========================================================================================
module inversion_arrays

  use global_typedefs
  use mpi_tools, only: exit_MPI

  implicit none

  private

  ! Contains allocatable arrays needed for inversion.
  type, public :: t_inversion_arrays

    ! Difference between data measured and data calculated.
    real(kind=CUSTOM_REAL), allocatable :: residuals(:)
    ! Weights to scale the sensitivity matrix columns.
    real(kind=CUSTOM_REAL), allocatable :: column_weight(:)

  contains
    private

    procedure, public, pass :: allocate_aux => inversion_arrays_allocate_aux
    procedure, public, pass :: reallocate_aux => inversion_arrays_reallocate_aux

  end type t_inversion_arrays

contains

!============================================================================================
! Allocates the auxiliarily inversion arrays.
!============================================================================================
subroutine inversion_arrays_allocate_aux(this, nelements, ndata, myrank)
  class(t_inversion_arrays), intent(inout) :: this
  integer, intent(in) :: nelements, ndata
  integer, intent(in) :: myrank

  integer :: ierr

  if (myrank == 0) print *, "Allocating auxiliarily inversion arrays..."

  if (ndata <= 0 .or. nelements <= 0) &
    call exit_MPI("Wrong dimensions in inversion_arrays_allocate_aux!", myrank, 0)

  ierr = 0

  allocate(this%residuals(ndata_components * ndata), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%column_weight(nelements), source=1._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in inversion_arrays_allocate_aux!", myrank, ierr)

end subroutine inversion_arrays_allocate_aux

!============================================================================================
! Reallocates the auxiliarily inversion arrays using new dimensions.
! Note: it reallocates only the arrays which dimenion has changed.
!============================================================================================
subroutine inversion_arrays_reallocate_aux(this, nelements, ndata, myrank)
  class(t_inversion_arrays), intent(inout) :: this
  integer, intent(in) :: nelements, ndata
  integer, intent(in) :: myrank

  integer :: ierr

  if (myrank == 0) print *, "Reallocating auxiliarily inversion arrays..."

  if (ndata <= 0 .or. nelements <= 0) &
    call exit_MPI("Wrong dimensions in inversion_arrays_allocate_aux!", myrank, 0)

  ierr = 0

  if (size(this%residuals) /= ndata * ndata_components) then
    deallocate(this%residuals)
    allocate(this%residuals(ndata_components * ndata), source=0._CUSTOM_REAL, stat=ierr)
  endif

  if (size(this%column_weight) /= nelements) then
    deallocate(this%column_weight)
    allocate(this%column_weight(nelements), source=1._CUSTOM_REAL, stat=ierr)
  endif

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in inversion_arrays_reallocate_aux!", myrank, ierr)

end subroutine inversion_arrays_reallocate_aux

end module inversion_arrays
