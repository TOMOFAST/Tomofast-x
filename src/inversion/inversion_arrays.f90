
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
    procedure, public, pass :: allocate_aux => inversion_arrays_allocate_aux
    procedure, public, pass :: allocate_sensit => inversion_arrays_allocate_sensit
    procedure, public, pass :: init_model => inversion_arrays_init_model

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
! Allocates the auxiliarily inversion arrays.
! NOTE: model objects are allocated separately in init_models().
!============================================================================================
subroutine inversion_arrays_allocate_aux(this, myrank)
  class(t_inversion_arrays), intent(inout) :: this
  integer, intent(in) :: myrank

  integer :: ierr

  if (myrank == 0) print *, "Allocating auxiliarily inversion arrays..."

  if (this%ndata <= 0 .or. this%nelements <= 0) &
    call exit_MPI("Wrong dimensions in inversion_arrays_allocate_aux!", myrank, 0)

  ierr = 0

  if (.not. allocated(this%residuals)) allocate(this%residuals(this%ndata), source=0._CUSTOM_REAL, stat=ierr)
  if (myrank == 0) print *, "residuals done."

  if (.not. allocated(this%column_weight)) allocate(this%column_weight(this%nelements), source=1._CUSTOM_REAL, stat=ierr)
  if (myrank == 0) print *, "column_weight done."

  if (.not. allocated(this%damping_weight)) allocate(this%damping_weight(this%nelements), source=1._CUSTOM_REAL, stat=ierr)
  if (myrank == 0) print *, "damping_weight done."

  if (.not. allocated(this%model_prior)) allocate(this%model_prior(this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (myrank == 0) print *, "model_prior done."

end subroutine inversion_arrays_allocate_aux

!==================================================================================
! Allocates memory for the sensitivity kernel.
! USE_LEGACY_SENSIT_MATRIX flag to support old sensitivity matrix allocation (non-sparse) for the ECT problem.
!==================================================================================
subroutine inversion_arrays_allocate_sensit(this, USE_LEGACY_SENSIT_MATRIX, nnz, myrank)
  class(t_inversion_arrays), intent(inout) :: this
  integer, intent(in) :: nnz
  logical, intent(in) :: USE_LEGACY_SENSIT_MATRIX
  integer, intent(in) :: myrank

  integer :: nl
  integer :: ierr

  if (myrank == 0) print *, "Allocating sensitivity kernel..."

  if (this%ndata <= 0 .or. this%nelements <= 0) &
    call exit_MPI("Wrong dimensions in inversion_arrays_allocate_sensit!", myrank, 0)

  ierr = 0

  if (USE_LEGACY_SENSIT_MATRIX) then

    if (.not. allocated(this%sensitivity)) &
      allocate(this%sensitivity(this%nelements, this%ndata), source=0._CUSTOM_REAL, stat=ierr)

    if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in inversion_arrays_allocate_sensit!", myrank, ierr)

    if (myrank == 0) print *, "sensitivity done."


  else
    nl = this%ndata

    ! Allocating the sparse matrix.
    call this%matrix_sensit%initialize(nl, nnz, myrank)
    if (myrank == 0) print *, "sensitivity (sparse) done."
  endif

  if (myrank == 0) print *, "Sensitivity kernel allocated!"

end subroutine inversion_arrays_allocate_sensit

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

end module inversion_arrays
