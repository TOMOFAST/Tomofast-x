
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

!========================================================================================
! Parameters needed for forward problems that are common for gravity and magnetism.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!========================================================================================
module parameters_gravmag

  use global_typedefs
  use mpi_tools, only: exit_MPI

  implicit none

  private

  type, public :: t_parameters_base

    ! Problem dimensions.
    integer :: nx, ny, nz
    ! The number of model elements (on current CPU).
    integer :: nelements
    ! Number of data.
    integer :: ndata
    ! Number of data components.
    integer :: ndata_components
    ! Number of model components.
    integer :: nmodel_components

    ! File name for the data.
    character(len=256) :: data_file
    ! File name for the data grid.
    character(len=256) :: data_grid_file

    ! File names for the model (with grid).
    ! (1): Read model.
    ! (2): Prior model.
    ! (3): Starting model.
    character(len=256) :: model_files(3)

    ! Type of prior model.
    integer :: prior_model_type
    ! Number of prior models (for model type 2).
    integer :: number_prior_models
    ! Set prior model to this value.
    real(kind=CUSTOM_REAL) :: prior_model_val

    ! Type of starting model.
    integer :: start_model_type
    ! Set starting model to this value.
    real(kind=CUSTOM_REAL) :: start_model_val

    !------ Depth weighting ------------------------------------------------
    ! Type of the depth weighting (1-depth weight, 2-distance weight).
    integer :: depth_weighting_type
    ! Power constant for depth weighting.
    real(kind=CUSTOM_REAL) :: depth_weighting_power
    ! Power constant for depth weighting strength (external power).
    real(kind=CUSTOM_REAL) :: depth_weighting_beta
    ! Empirical constant (Z-shift) for depth weighting type #1.
    real(kind=CUSTOM_REAL) :: Z0

    ! Local depth weight preconditioning.
    integer :: apply_local_weight
    character(len=256) :: local_weight_file

    !------ Matrix compression ---------------------------------------------
    ! Parameters for reduction of the memory requirements (to store the sensitivity matrix).
    ! 0 -none, 1 - wavelet
    integer :: compression_type
    real(kind=CUSTOM_REAL) :: compression_rate

    !------ Sensitivity kernel ---------------------------------------------
    integer :: sensit_read
    character(len=256) :: sensit_path

    !------ Covariance -----------------------------------------------------
    integer :: use_data_cov
    character(len=256) :: data_cov_file

    !------ Other ----------------------------------------------------------
    ! Model units conversion.
    real(kind=CUSTOM_REAL) :: model_units_mult
    ! Data units conversion.
    real(kind=CUSTOM_REAL) :: data_units_mult
    ! Direction of the Z-axis (1 = down, -1 = up).
    integer :: z_axis_dir
    ! The model output label in the vtk.
    character(len=16) :: vtk_model_label

  contains
    private

    procedure, public, pass :: broadcast => parameters_base_broadcast
  end type t_parameters_base

contains

!=========================================================================
! MPI broadcast parameters that are read from a Parfile.
!=========================================================================
subroutine parameters_base_broadcast(this, myrank)
  class(t_parameters_base), intent(in) :: this
  integer, intent(in) :: myrank
  integer :: ierr

  ierr = 0

  call MPI_Bcast(this%nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%ny, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%nz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%ndata, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%ndata_components, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%nmodel_components, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%prior_model_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%number_prior_models, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%prior_model_val, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%start_model_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%start_model_val, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%depth_weighting_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%depth_weighting_power, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%depth_weighting_beta, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%Z0, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%apply_local_weight, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%local_weight_file, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%compression_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%compression_rate, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%sensit_read, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%sensit_path, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%use_data_cov, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%data_cov_file, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%model_units_mult, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%data_units_mult, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%z_axis_dir, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI_Bcast error in parameters_base_broadcast!", myrank, ierr)

end subroutine parameters_base_broadcast

end module parameters_gravmag
