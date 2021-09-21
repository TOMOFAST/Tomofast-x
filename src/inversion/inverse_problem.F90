
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
! A class to perform parallel inversion of a single data field.
! Works with sparse matrices that are stored using Compressed Sparse Row (CSR) format.
! Calculates the model update (change).
! Uses an object of type t_parameters_inversion to obtain the input parameters and data arrays.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!===============================================================================================
module inverse_problem

  use global_typedefs, only: CUSTOM_REAL
  use mpi_tools, only: exit_MPI
  use parameters_inversion
  use inversion_arrays
  use sparse_matrix
  use model
  use sensitivity_matrix
  use damping
  use lsqr_solver

  implicit none

  private

  ! Main inversion class.
  type, public :: t_inversion
    private

    ! Matrix: sensitivity kernel and damping contribution.
    type(t_sparse_matrix) :: matrix
    ! Right hand side: data residuals with damping extension.
    real(kind=CUSTOM_REAL), allocatable :: b_RHS(:)
    ! Model change.
    real(kind=CUSTOM_REAL), allocatable :: delta_model(:)

    ! Number of elements in the model.
    integer :: nelements

  contains
    private

    procedure, public, pass :: initialize => inversion_initialize
    procedure, public, pass :: reset => inversion_reset
    procedure, public, pass :: solve => inversion_solve
    procedure, public, pass :: get_model_change => inversion_get_model_change

  end type t_inversion

contains

!=======================================================================
! Initialize inversion.
!=======================================================================
subroutine inversion_initialize(this, par, myrank)
  class(t_inversion), intent(inout) :: this
  type(t_parameters_inversion), intent(in) :: par
  integer, intent(in) :: myrank
  integer :: ierr

  if (myrank == 0) print *, 'ndata, nelements =', par%ndata(1), par%nelements

  this%nelements = par%nelements

  call this%matrix%initialize(par%ndata(1) + par%nelements_total, &
                              par%ndata(1) * par%nelements + par%nelements, myrank)

  ierr = 0

  ! Memory allocation.
  if (.not. allocated(this%b_RHS)) &
    allocate(this%b_RHS(this%matrix%get_total_row_number()), source=0._CUSTOM_REAL, stat=ierr)

  if (.not. allocated(this%delta_model)) &
    allocate(this%delta_model(this%nelements), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in inversion_initialize!", myrank, ierr)

end subroutine inversion_initialize

!=======================================================================
! Resets the inversion.
!=======================================================================
pure subroutine inversion_reset(this)
  class(t_inversion), intent(inout) :: this

  call this%matrix%reset()

  this%b_RHS = 0._CUSTOM_REAL
  this%delta_model = 0._CUSTOM_REAL

end subroutine inversion_reset

!========================================================================
! Returns model change, obtained in solve().
!========================================================================
pure function inversion_get_model_change(this) result(res)
  class(t_inversion), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res(this%nelements)

  res = this%delta_model

end function inversion_get_model_change

!========================================================================================
! Inversion of one data set.
! Use legacy sensitivity matrix format to support the ECT inversion.
!========================================================================================
subroutine inversion_solve(this, par, arr, myrank, nbproc)
  class(t_inversion), intent(inout) :: this
  type(t_parameters_inversion), intent(in) :: par
  type(t_inversion_arrays), intent(in) :: arr
  integer, intent(in) :: myrank, nbproc

  type(t_sensitivity_matrix) :: sensit
  type(t_damping) :: damping
  type(t_sparse_matrix) :: matrix_dummy

  real(kind=CUSTOM_REAL), parameter :: problem_weight = 1.d0

  logical, parameter :: USE_LEGACY_SENSIT_MATRIX = .true.

  call sensit%initialize(par%ndata(1), par%nelements, problem_weight)

  ! Add compressed sparse sensitivity matrix in CSR format.
  call sensit%add(this%matrix, this%b_RHS, 0, arr%sensitivity, matrix_dummy, &
                  arr%column_weight, arr%residuals, USE_LEGACY_SENSIT_MATRIX, .true., myrank)

  if (myrank == 0) print *, 'nel = ', this%matrix%get_number_elements()

  ! Add damping.
  call damping%initialize(par%nelements, par%alpha(1), problem_weight, par%norm_power, &
                          par%compression_type, par%nx, par%ny, par%nz)
  call damping%add(this%matrix, this%b_RHS, arr%column_weight, arr%damping_weight, &
                   arr%model, arr%model_prior, 0, myrank, nbproc)

  if (myrank == 0) print *, 'nel (with damping) = ', this%matrix%get_number_elements()

  call this%matrix%finalize(par%nelements, myrank)

  ! Parallel sparse inversion.
  if (par%method == 1) then
    this%delta_model = 0._CUSTOM_REAL
    call lsqr_solve(par%niter, par%rmin, par%gamma, &
                    this%matrix, this%b_RHS, this%delta_model, myrank)

  else if (par%method == 2) then
    ! TODO: BFGS
  endif

! DISABLED.
!  if (par%compress_matrix) then
!    ! Inverse wavelet transform.
!    call iHaar3D(this%delta_model, par%nx, par%ny, par%nz)
!  endif

  ! Scale model back to the original variables.
  call rescale_model(this%delta_model, arr%column_weight, par%nelements)

end subroutine inversion_solve

end module inverse_problem
