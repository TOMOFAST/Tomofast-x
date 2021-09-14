
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
! A class to add (inversion) damping contribution to the matrix and right hand side.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!===========================================================================================
module damping

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use sparse_matrix
  use parallel_tools
  use model
  use grid
  use wavelet_transform

  implicit none

  private

  type, public :: t_damping
    private

    ! Number of model parameters.
    integer :: nelements
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
    real(kind=CUSTOM_REAL) :: threshold

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
                              compression_type, nx, ny, nz, threshold)
  class(t_damping), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: alpha, problem_weight, norm_power, threshold
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
  this%threshold = threshold
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
subroutine damping_add(this, matrix, b_RHS, column_weight, local_weight, &
                       model, model_ref, param_shift, myrank, nbproc)
  class(t_damping), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  real(kind=CUSTOM_REAL), intent(in) :: local_weight(:)
  type(t_model), intent(in) :: model
  real(kind=CUSTOM_REAL), intent(in) :: model_ref(:)
  integer, intent(in) :: param_shift
  integer, intent(in) :: myrank, nbproc

  type(t_sparse_matrix), intent(inout) :: matrix
  real(kind=CUSTOM_REAL), intent(inout) :: b_RHS(:)

  integer :: i, nsmaller, nelements_total
  integer :: row_beg, row_end
  integer :: ierr
  real(kind=CUSTOM_REAL) :: value
  type(t_parallel_tools) :: pt

  !---------------------------------------------------------------------
  ! Calculating the model difference: (m - m_ref), which is used in the right-hand side, and in the Lp norm multiplier.
  real(kind=CUSTOM_REAL), allocatable :: model_diff(:)
  real(kind=CUSTOM_REAL), allocatable :: model_diff_full(:)

  allocate(model_diff(this%nelements), source=0._CUSTOM_REAL, stat=ierr)

  model_diff = model%val - model_ref

  ! Apply the depth-weighting.
  do i = 1, this%nelements
    model_diff(i) = model_diff(i) / column_weight(i)
  enddo

  ! The total number of elements.
  nelements_total = pt%get_total_number_elements(this%nelements, myrank, nbproc)

  ! The number of elements on CPUs with rank smaller than myrank.
  nsmaller = pt%get_nsmaller(this%nelements, myrank, nbproc)

  if (this%compression_type == 2) then
  ! Transform the model difference to the wavelet domain.

    if (nbproc > 1) then
    ! Parallel version.

      ! Allocate memory for the full model.
      allocate(model_diff_full(nelements_total), source=0._CUSTOM_REAL, stat=ierr)
      if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in damping_add!", myrank, ierr)

      ! Gather the full model from all processors.
      call pt%get_full_array(model_diff, this%nelements, model_diff_full, .true., myrank, nbproc)

      ! Transform to wavelet domain.
      call Haar3D(model_diff_full, this%nx, this%ny, this%nz)

      ! Extract the local model part.
      model_diff = model_diff_full(nsmaller + 1 : nsmaller + this%nelements)

      deallocate(model_diff_full)

    else
    ! Serial version.
      ! Transform to wavelet domain.
      call Haar3D(model_diff, this%nx, this%ny, this%nz)
    endif
  endif
  !---------------------------------------------------------------------

  ! First matrix row (in the big matrix) of the damping matrix that will be added.
  row_beg = matrix%get_current_row_number() + 1

  ! Add empty lines.
  call matrix%add_empty_rows(nsmaller, myrank)

  ! Add lines with damping.
  do i = 1, this%nelements
    call matrix%new_row(myrank)

    value = this%alpha * this%problem_weight

    ! Apply the Lp norm.
    value = value * this%get_norm_multiplier(model_diff(i))

    ! Apply local weight, which is equivalent to having local alpha.
    value = value * local_weight(i)

    call matrix%add(value, param_shift + i, myrank)

  enddo

  ! Add empty lines.
  call matrix%add_empty_rows(nelements_total - this%nelements - nsmaller, myrank)

  ! Last matrix row (in the big matrix) of the added damping matrix.
  row_end = matrix%get_current_row_number()

  ! Sanity check.
  if (row_end - row_beg + 1 /= nelements_total) &
    call exit_MPI("Sanity check failed in damping_add!", myrank, 0)

  !---------------------------------------------------------------------
  ! Add the damping contribution to the right hand side.
  call this%add_RHS(b_RHS(row_beg:row_end), model_diff, column_weight, local_weight, myrank, nbproc)

  deallocate(model_diff)

  ! Calculate the damping cost.
  this%cost = sum(b_RHS(row_beg:row_end)**2)

end subroutine damping_add

!=============================================================================================
! Adds damping contribution in the right hand side.
! model_ref - reference model.
!=============================================================================================
subroutine damping_add_RHS(this, b_RHS, model_diff, column_weight, local_weight, myrank, nbproc)
  class(t_damping), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: model_diff(:)
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  real(kind=CUSTOM_REAL), intent(in) :: local_weight(:)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), intent(inout) :: b_RHS(:)

  type(t_parallel_tools) :: pt
  integer :: i

  do i = 1, this%nelements
    b_RHS(i) = - this%alpha * this%problem_weight * model_diff(i)

    ! Apply the Lp norm.
    b_RHS(i) = b_RHS(i) * this%get_norm_multiplier(model_diff(i))

    ! Apply local weight, which is equivalent to having local alpha.
    b_RHS(i) = b_RHS(i) * local_weight(i)
  enddo

  ! Gather full right hand side.
  call pt%get_full_array_in_place(this%nelements, b_RHS, .true., myrank, nbproc)

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
