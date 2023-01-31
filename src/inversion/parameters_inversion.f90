
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

    ! Error in prior model -- damping for the inverse problem.
    real(kind=CUSTOM_REAL) :: alpha(2)
    ! Power p of Lp norm on the model damping term (for LSQR method). Use p=2. for pure LSQR.
    real(kind=CUSTOM_REAL) :: norm_power

    ! Damping gradient weight type (1-global, 2-local).
    integer :: damp_grad_weight_type
    ! Damping gradient global weight.
    real(kind=CUSTOM_REAL) :: beta(2)

    ! Stopping criterion.
    real(kind=CUSTOM_REAL) :: rmin
    ! LSQR=1, SCA=2.
    integer :: method
    ! LSQR threshold to have kinda L1 norm. Use gamma=0 for no thresholding.
    real(kind=CUSTOM_REAL) :: gamma

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

    ! The number of iterations in the method of weights (equality-constrained LSQR).
    integer :: method_of_weights_niter

    ! The type of derivative use for cross-gradients discretization:
    ! 1 - forward, 2 - central.
    integer :: derivative_type

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
  class(t_parameters_inversion), intent(in) :: this
  integer, intent(in) :: myrank
  integer :: ierr

  ierr = 0

  call MPI_Bcast(this%niter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%ninversions, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%write_model_niter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%rmin, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%method, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%gamma, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%alpha, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%norm_power, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%damp_grad_weight_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%beta, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%compression_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%wavelet_threshold, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  ! Joint inversion parameters.
  call MPI_Bcast(this%problem_weight, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%column_weight_multiplier, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  ! Cross-gradient parameters.
  call MPI_Bcast(this%cross_grad_weight, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%method_of_weights_niter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%derivative_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  ! Clustering parameters.
  call MPI_Bcast(this%clustering_weight_glob, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%nclusters, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%clustering_opt_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%clustering_constraints_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  ! ADMM parameters.
  call MPI_Bcast(this%admm_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%nlithos, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%bounds_ADMM_file, 512, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%rho_ADMM, 2, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  ! ADMM parameters (for dynamic weight adjustment).
  call MPI_Bcast(this%data_cost_threshold_ADMM, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%weight_multiplier_ADMM, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%max_weight_ADMM, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI_Bcast error in parameters_inversion_broadcast!", myrank, ierr)

end subroutine parameters_inversion_broadcast

end module parameters_inversion
