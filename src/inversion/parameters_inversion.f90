
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
! A class that stores input parameters needed for inversion.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!==========================================================================================
module parameters_inversion

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use parallel_tools

  implicit none

  private

  ! Contains scalar parameters needed for inversion (NO allocatable arrays!).
  type, public :: t_parameters_inversion

    ! Original dimensions of the model (i.e., before splitting among CPUs).
    integer :: nx, ny, nz
    ! The local number of elements in the model on every CPU.
    integer :: nelements
    ! The total number of elements in the model (sum of nelements on every CPU).
    integer :: nelements_total
    ! Number of data.
    integer :: ndata(2)
    ! Number of data components.
    integer :: ndata_components(2)
    ! Number of model components.
    integer :: nmodel_components
    ! Number of iterations in inversion solver (e.g. LSQR).
    integer :: niter
    ! Number of inversions (outer inversion loop).
    integer :: ninversions
    ! Write intermediate model every N major iterations.
    integer :: write_model_niter

    ! Global model damping weight.
    real(kind=CUSTOM_REAL) :: alpha(2)
    ! Power p of Lp norm on the model damping term (for LSQR method). Use p=2. for pure LSQR.
    real(kind=CUSTOM_REAL) :: norm_power

    ! Local damping weight.
    integer :: apply_local_damping_weight
    character(len=256) :: damping_weight_file(2)

    ! Damping gradient weight type (1-global, 2-local).
    integer :: damp_grad_weight_type
    ! Damping gradient global weight.
    real(kind=CUSTOM_REAL) :: beta(2)
    ! File with local damping gradient weights.
    character(len=256) :: damping_gradient_file(2)

    ! Minimum relative residual.
    real(kind=CUSTOM_REAL) :: rmin
    ! Target data misfit.
    real(kind=CUSTOM_REAL) :: target_misfit
    ! LSQR threshold to have approx. L1 norm. Use gamma=0 for no thresholding.
    real(kind=CUSTOM_REAL) :: gamma
    ! LSQR=1.
    integer :: method

    ! Wavelet compression parameters.
    integer :: compression_type
    real(kind=CUSTOM_REAL) :: wavelet_threshold

    ! ------ Joint inversion parameters ------------------------------------

    ! Weights of individual problems (grav, mag) in joint inversion.
    real(kind=CUSTOM_REAL) :: problem_weight(2)

    ! Multipliers of the column weight for individual problems (grav, mag) in joint inversion.
    ! To control the influence of the cross-grad term on grav/mag problems.
    real(kind=CUSTOM_REAL) :: column_weight_multiplier(2)

    ! ------ Cross-gradient constraints ------------------------------------

    ! Contribution of the cross-gradient to the cost function.
    real(kind=CUSTOM_REAL) :: cross_grad_weight
    ! The type of derivative use for cross-gradients discretization:
    ! 1 - forward, 2 - central.
    integer :: derivative_type
    ! Flags to define one of the models constant so that it remains unaltered during the inversion.
    integer :: keep_model_constant(2)

    ! A flag for using a vector field for structural constraints.
    integer :: vec_field_type
    ! Vector field file (for structural constraints).
    character(len=256) :: vec_field_file

    ! ------ Clustering constraints ----------------------------------------

    ! Clustering weights.
    real(kind=CUSTOM_REAL) :: clustering_weight_glob(2)
    ! Number of clusters.
    integer :: nclusters
    ! Optimization type (normal or logarithmic).
    integer :: clustering_opt_type
    ! Type of constraints (local vs global).
    integer :: clustering_constraints_type
    ! File name for clustering mixtures (for clustering constraint).
    character(len=256) :: mixture_file
    ! File name for clustering weights per grid cell.
    character(len=256) :: cell_weights_file

    ! ------ ADMM constraints ------------------------------------

    ! 0 - no admm, 1 - with admm.
    integer :: admm_type
    ! 1 - global, 1 - local.
    integer :: admm_bound_type
    ! ADMM bounds for the global bound type.
    type(t_real1d) :: admm_bounds(2)
    ! Number of lithologies.
    integer :: nlithos
    character(len=256) :: bounds_ADMM_file(2)
    real(kind=CUSTOM_REAL) :: rho_ADMM(2)

    ! ADMM parameters for dynamic weight adjustment.
    real(kind=CUSTOM_REAL) :: data_cost_threshold_ADMM
    real(kind=CUSTOM_REAL) :: weight_multiplier_ADMM
    real(kind=CUSTOM_REAL) :: max_weight_ADMM

  contains
    private

    procedure, public, pass :: broadcast => parameters_inversion_broadcast

  end type t_parameters_inversion

contains

!==================================================================================
! MPI broadcast parameters that are read from a Parfile.
!==================================================================================
subroutine parameters_inversion_broadcast(this, myrank)
  class(t_parameters_inversion), intent(inout) :: this
  integer, intent(in) :: myrank
  integer :: ierr

  ierr = 0

  call MPI_Bcast(this%niter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%ninversions, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%write_model_niter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%rmin, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%target_misfit, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%gamma, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%method, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%alpha, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%norm_power, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%apply_local_damping_weight, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%damping_weight_file, 512, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%damp_grad_weight_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%beta, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%damping_gradient_file, 512, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%compression_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%wavelet_threshold, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  ! Joint inversion parameters.
  call MPI_Bcast(this%problem_weight, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%column_weight_multiplier, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  ! Cross-gradient parameters.
  call MPI_Bcast(this%cross_grad_weight, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%derivative_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%keep_model_constant, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%vec_field_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%vec_field_file, 512, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  ! Clustering parameters.
  call MPI_Bcast(this%clustering_weight_glob, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%nclusters, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%clustering_opt_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%clustering_constraints_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  ! ADMM parameters.
  call MPI_Bcast(this%admm_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%admm_bound_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%nlithos, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%bounds_ADMM_file, 512, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%rho_ADMM, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (this%admm_type == 1 .and. this%admm_bound_type == 1) then
    ! They get allocated in init_parameters on the master rank. Allocate for all ranks.
    if (.not. allocated(this%admm_bounds(1)%val)) then
      allocate(this%admm_bounds(1)%val(2 * this%nlithos), source=0._CUSTOM_REAL)
    endif
    if (.not. allocated(this%admm_bounds(2)%val)) then
      allocate(this%admm_bounds(2)%val(2 * this%nlithos), source=0._CUSTOM_REAL)
    endif
    call MPI_Bcast(this%admm_bounds(1)%val, 2 * this%nlithos, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(this%admm_bounds(2)%val, 2 * this%nlithos, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  endif

  ! ADMM parameters (for dynamic weight adjustment).
  call MPI_Bcast(this%data_cost_threshold_ADMM, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%weight_multiplier_ADMM, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%max_weight_ADMM, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI_Bcast error in parameters_inversion_broadcast!", myrank, ierr)

end subroutine parameters_inversion_broadcast

end module parameters_inversion
