
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

!==========================================================================================
! A class that stores allocatable arrays (and types containing them as e.g. t_model) needed for inversion.
! Memory allocation is done here.
!
! Note: we use sourced allocation (with source=0 parameter), which also initializes the arrays.
! This is because allocate() allocates virtual memory,
! and it can be that there is less free physical memory available than virtual memory allocated.
! In this case the program will be killed by OS during initialization.
! When "source=0" parameter is used, then it will invoke an error via "stat=ierr" in such case.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!==========================================================================================
module inversion_arrays

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use model
  use sparse_matrix

  implicit none

  private

  ! Contains allocatable arrays needed for inversion.
  type, public :: t_inversion_arrays

    ! Required parameters to allocate arrays.
    integer :: nelements, ndata
    integer :: nx, ny, nz

    ! Difference between data measured and data calculated.
    real(kind=CUSTOM_REAL), allocatable :: residuals(:)
    ! Prior model values - first initial guess for inversion.
    real(kind=CUSTOM_REAL), allocatable :: model_prior(:)
    ! Weights to scale the sensitivity matrix columns.
    real(kind=CUSTOM_REAL), allocatable :: column_weight(:)
    ! Weights to scale the damping.
    real(kind=CUSTOM_REAL), allocatable :: damping_weight(:)

    ! Sensitivity kernel computed in forward problem.
    real(kind=CUSTOM_REAL), allocatable :: sensitivity(:, :)
    type(t_sparse_matrix) :: matrix_sensit

    ! Solution model (with grid).
    type(t_model) :: model

  contains
    private

    procedure, public, pass :: initialize => inversion_arrays_initialize
    procedure, public, pass :: allocate => inversion_arrays_allocate
    procedure, public, pass :: init_model => inversion_arrays_init_model

    ! Destructor. (Note: bug in gcc 4.9 with warning about 'array final procedure'.)
    final :: inversion_arrays_destructor

  end type t_inversion_arrays

contains

!============================================================================================
! Initialize parameters needed to allocate the inversion arrays.
!============================================================================================
subroutine inversion_arrays_initialize(this, nelements, ndata, nx, ny, nz)
  class(t_inversion_arrays), intent(inout) :: this
  integer, intent(in) :: nelements, ndata, nx, ny, nz

  this%nelements = nelements
  this%ndata = ndata
  this%nx = nx
  this%ny = ny
  this%nz = nz

end subroutine inversion_arrays_initialize

!============================================================================================
! Allocates the inversion arrays.
! NOTE: model objects are allocated separately in init_models().
! USE_LEGACY_SENSIT_MATRIX flag to support old sensitivity matrix allocation (non-sparse) for the ECT problem.
!============================================================================================
subroutine inversion_arrays_allocate(this, USE_LEGACY_SENSIT_MATRIX, nnz, myrank)
  class(t_inversion_arrays), intent(inout) :: this
  integer, intent(in) :: nnz
  logical, intent(in) :: USE_LEGACY_SENSIT_MATRIX
  integer, intent(in) :: myrank

  integer :: nl
  integer :: ierr

  if (myrank == 0) print *, "Allocating inversion arrays..."

  if (this%ndata <= 0 .or. this%nelements <= 0) &
    call exit_MPI("Wrong dimensions in inversion_arrays_initialize!", myrank, 0)

  ierr = 0

  if (.not. allocated(this%residuals)) allocate(this%residuals(this%ndata), source=0._CUSTOM_REAL, stat=ierr)
  if (myrank == 0) print *, "residuals done."

  if (.not. allocated(this%column_weight)) allocate(this%column_weight(this%nelements), source=1._CUSTOM_REAL, stat=ierr)
  if (myrank == 0) print *, "column_weight done."

  if (.not. allocated(this%damping_weight)) allocate(this%damping_weight(this%nelements), source=1._CUSTOM_REAL, stat=ierr)
  if (myrank == 0) print *, "damping_weight done."

  if (.not. allocated(this%model_prior)) allocate(this%model_prior(this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (myrank == 0) print *, "model_prior done."

  if (.not. allocated(this%sensitivity) .and. USE_LEGACY_SENSIT_MATRIX) &
    allocate(this%sensitivity(this%nelements, this%ndata), source=0._CUSTOM_REAL, stat=ierr)
  if (myrank == 0) print *, "sensitivity done."


  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in inversion_arrays_initialize!", myrank, ierr)

  if (.not. USE_LEGACY_SENSIT_MATRIX) then
    nl = this%ndata

    ! Allocating the sparse matrix.
    call this%matrix_sensit%initialize(nl, nnz, myrank)
    if (myrank == 0) print *, "sensitivity (sparse) done."
  endif

  if (myrank == 0) print *, "Inversion arrays allocated!"

end subroutine inversion_arrays_allocate

!==================================================================================
! Allocates memory for the model (and its grid).
!==================================================================================
subroutine inversion_arrays_init_model(this, myrank, nbproc)
  class(t_inversion_arrays), intent(inout) :: this
  integer, intent(in) :: myrank, nbproc

  if (this%nelements <= 0) &
    call exit_MPI("Wrong dimensions in inversion_arrays_init_model!", myrank, 0)

  call this%model%initialize(this%nelements, myrank, nbproc)
  if (myrank == 0) print *, "model done."

  call this%model%init_grid(this%nx, this%ny, this%nz, myrank)
  if (myrank == 0) print *, "grid done."

end subroutine inversion_arrays_init_model

!==================================================================================
! Destructor.
!==================================================================================
subroutine inversion_arrays_destructor(this)
  type(t_inversion_arrays), intent(inout) :: this

  if (allocated(this%sensitivity)) deallocate(this%sensitivity)
  if (allocated(this%column_weight)) deallocate(this%column_weight)
  if (allocated(this%damping_weight)) deallocate(this%damping_weight)
  if (allocated(this%residuals)) deallocate(this%residuals)
  if (allocated(this%model_prior)) deallocate(this%model_prior)

end subroutine inversion_arrays_destructor

end module inversion_arrays
