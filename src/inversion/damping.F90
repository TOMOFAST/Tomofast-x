
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
! A class to add (model) damping contribution to the matrix and right hand side.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!===========================================================================================
module damping

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use sparse_matrix
  use parallel_tools
  use wavelet_utils

  implicit none

  private

  type, public :: t_damping
    private

    ! Number of model parameters (on current rank).
    integer :: nelements
    ! Total number of elements.
    integer :: nelements_total
    ! Weight of the damping.
    real(kind=CUSTOM_REAL) :: alpha
    ! Weight of the whole problem (damping + misfit) in joint inversion.
    real(kind=CUSTOM_REAL) :: problem_weight
    ! Power p of Lp norm (for LSQR method). Use p=2. for pure LSQR.
    real(kind=CUSTOM_REAL) :: norm_power
    ! Cost of the model objective function phi_m.
    real(kind=CUSTOM_REAL) :: cost

    ! Wavelet compression parameters.
    integer :: compression_type
    integer :: nx, ny, nz

  contains
    private

    procedure, public, pass :: initialize => damping_initialize
    procedure, public, pass :: add => damping_add
    procedure, public, pass :: get_cost => damping_get_cost

    procedure, private, pass :: add_RHS => damping_add_RHS
    procedure, private, pass :: get_norm_multiplier => damping_get_norm_multiplier

  end type t_damping

contains

!===========================================================================================
! Initialization.
!===========================================================================================
subroutine damping_initialize(this, nelements, alpha, problem_weight, norm_power, &
                              compression_type, nx, ny, nz)
  class(t_damping), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: alpha, problem_weight, norm_power
  integer, intent(in) :: nelements
  integer, intent(in) :: compression_type, nx, ny, nz

  this%nelements = nelements
  this%alpha = alpha
  this%problem_weight = problem_weight
  this%norm_power = norm_power
  this%compression_type = compression_type
  this%nx = nx
  this%ny = ny
  this%nz = nz
  this%nelements_total = this%nx * this%ny * this%nz
end subroutine damping_initialize

!===========================================================================================
! 1) Adds damping (below the sensitivity kernel) which is identity matrix times alpha:
!     S_new = (  S )
!             ( aI ),
!    where I is identity matrix, a=alpha is damping parameter,
!    applying the same scaling as applied to the sensitivity matrix.
! 2) Adds the corresponding contribution to the right hand side 'b_RHS'.
!
! Tested in unit_tests.f90 in test_damping_identity_matrix().
!===========================================================================================
subroutine damping_add(this, matrix, nrows, b_RHS, column_weight, &
                       model, model_ref, param_shift, WAVELET_DOMAIN, myrank, nbproc, local_weight)
  class(t_damping), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(this%nelements)
  real(kind=CUSTOM_REAL), intent(in) :: model(this%nelements, 3)
  real(kind=CUSTOM_REAL), intent(in) :: model_ref(this%nelements)
  integer, intent(in) :: nrows, param_shift
  logical, intent(in) :: WAVELET_DOMAIN
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), optional, intent(in) :: local_weight(this%nelements)

  type(t_sparse_matrix), intent(inout) :: matrix
  real(kind=CUSTOM_REAL), intent(inout) :: b_RHS(nrows)

  integer :: i, nsmaller
  integer :: row_beg, row_end
  integer :: ierr
  real(kind=CUSTOM_REAL) :: value
  logical :: SOLVE_PROBLEM(1)

  !---------------------------------------------------------------------
  ! Calculating the model difference: (m - m_ref), which is used in the right-hand side, and in the Lp norm multiplier.
  real(kind=CUSTOM_REAL), allocatable :: model_diff(:)
  real(kind=CUSTOM_REAL), allocatable :: model_diff_full(:)

  real(kind=CUSTOM_REAL) :: model_mag

  allocate(model_diff(this%nelements), source=0._CUSTOM_REAL, stat=ierr)

  !model_diff = model - model_ref

  do i = 1, this%nelements
    ! Vector magnutude (squared).
    model_mag = model(i, 1)**1 + model(i, 2)**1 + model(i, 3)**1

    model_diff(i) = model_mag - model_ref(i)
  enddo

  ! Apply the depth-weighting.
  do i = 1, this%nelements
    model_diff(i) = model_diff(i) / column_weight(i)
  enddo

  if (this%compression_type > 0 .and. WAVELET_DOMAIN) then
    ! Transform the model difference to the wavelet domain.
    SOLVE_PROBLEM = .true.

    if (myrank == 0) then
      ! Allocate a buffer for wavelet transform on the master rank only.
      allocate(model_diff_full(this%nelements_total), source=0._CUSTOM_REAL, stat=ierr)
    else
      allocate(model_diff_full(1), source=0._CUSTOM_REAL, stat=ierr)
    endif

    call apply_wavelet_transform(this%nelements, this%nx, this%ny, this%nz, 1, &
                                 model_diff, model_diff_full, &
                                 .true., this%compression_type, 1, SOLVE_PROBLEM, myrank, nbproc)
  endif

  ! First matrix row (in the big matrix) of the damping matrix that will be added.
  row_beg = matrix%get_current_row_number() + 1

  ! The number of elements on CPUs with rank smaller than myrank.
  nsmaller = get_nsmaller(this%nelements, myrank, nbproc)

  ! Add empty lines.
  call matrix%add_empty_rows(nsmaller, myrank)

  ! Add lines with damping.
  do i = 1, this%nelements
!    value = this%alpha * this%problem_weight
!
!    if (this%norm_power /= 2.d0) then
!      ! Apply the Lp norm.
!      value = value * this%get_norm_multiplier(model_diff(i))
!    endif
!
!    if (present(local_weight)) then
!      ! Apply local weight, which is equivalent to having local alpha.
!      value = value * local_weight(i)
!    endif

    ! X-component
    value = this%alpha * this%problem_weight !* 2.d0 * model(i, 1) / column_weight(i)
    call matrix%add(value, param_shift + i, myrank)

    ! Y-component
    value = this%alpha * this%problem_weight !* 2.d0 * model(i, 2) / column_weight(i)
    call matrix%add(value, param_shift + i + this%nelements, myrank)

    ! Z-component
    value = this%alpha * this%problem_weight !* 2.d0 * model(i, 3) / column_weight(i)
    call matrix%add(value, param_shift + i + 2 * this%nelements, myrank)

    call matrix%new_row(myrank)
  enddo

  ! Add empty lines.
  call matrix%add_empty_rows(this%nelements_total - this%nelements - nsmaller, myrank)

  ! Last matrix row (in the big matrix) of the added damping matrix.
  row_end = matrix%get_current_row_number()

  ! Sanity check.
  if (row_end - row_beg + 1 /= this%nelements_total) &
    call exit_MPI("Sanity check failed in damping_add!", myrank, 0)

  !---------------------------------------------------------------------
  ! Add the damping contribution to the right hand side.
  if (present(local_weight)) then
    call this%add_RHS(b_RHS(row_beg:row_end), model_diff, myrank, nbproc, local_weight)
  else
    call this%add_RHS(b_RHS(row_beg:row_end), model_diff, myrank, nbproc)
  endif

  deallocate(model_diff)

  ! Calculate the damping cost.
  this%cost = sum(b_RHS(row_beg:row_end)**2)

end subroutine damping_add

!=============================================================================================
! Adds damping contribution in the right hand side.
! model_diff - depth weighted (m - m_prior).
!=============================================================================================
subroutine damping_add_RHS(this, b_RHS, model_diff, myrank, nbproc, local_weight)
  class(t_damping), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: model_diff(this%nelements)
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), optional, intent(in) :: local_weight(this%nelements)

  real(kind=CUSTOM_REAL), intent(inout) :: b_RHS(this%nelements_total)

  integer :: i

  do i = 1, this%nelements
    b_RHS(i) = - this%alpha * this%problem_weight * model_diff(i)

!    if (this%norm_power /= 2.d0) then
!      ! Apply the Lp norm.
!      b_RHS(i) = b_RHS(i) * this%get_norm_multiplier(model_diff(i))
!    endif
!
!    if (present(local_weight)) then
!      ! Apply local weight, which is equivalent to having local alpha.
!      b_RHS(i) = b_RHS(i) * local_weight(i)
!    endif
  enddo

  ! Gather full right hand side.
  call get_full_array_in_place(this%nelements, b_RHS, .true., myrank, nbproc)

end subroutine damping_add_RHS

!===========================================================================================
! Returns model objective function cost (norm).
!===========================================================================================
pure function damping_get_cost(this) result(res)
  class(t_damping), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = this%cost

end function damping_get_cost

!===========================================================================================
! Returns a multiplier (for one pixel) to change L2 norm to Lp, in the LSQR method.
!===========================================================================================
pure function damping_get_norm_multiplier(this, model_diff) result(res)
  class(t_damping), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: model_diff
  real(kind=CUSTOM_REAL) :: res

  if (model_diff /= 0.d0) then
    res = (abs(model_diff))**(this%norm_power / 2.d0 - 1.d0)
  else
    res = 1._CUSTOM_REAL
  endif

end function damping_get_norm_multiplier

end module damping
