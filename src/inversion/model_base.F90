!
!!========================================================================
!!
!!                          T o m o f a s t - x
!!                        -----------------------
!!
!!           Authors: Vitaliy Ogarko, Jeremie Giraud, Roland Martin.
!!
!!               (c) 2021 The University of Western Australia.
!!
!! The full text of the license is available in file "LICENSE".
!!
!!========================================================================
!
!!================================================================================================
!! This class contains basic model operations.
!!
!! Vitaliy Ogarko, UWA, CET, Australia.
!!================================================================================================
!module model_base
!
!  use global_typedefs
!  use mpi_tools, only: exit_MPI
!  use paraview
!  use grid
!  use string, only: str
!  use parallel_tools
!  use string
!
!  implicit none
!
!  private
!
!  type, public :: t_model_base
!
!    ! Local model parameters.
!    real(kind=CUSTOM_REAL), allocatable :: val(:)
!
!    ! Full model parameters.
!    ! (Read initial model here, write final model from here, and use also for cross-gradient.)
!    real(kind=CUSTOM_REAL), allocatable :: val_full(:)
!
!    ! TODO: Maybe move this to a separate module and combine with damping_weight?
!    ! Local model covariance (diagonal) matrix.
!    ! This is the weight that will be applied to the model damping term.
!    ! This is equivalent of having local alpha, i.e., changing with model cell.
!    real(kind=CUSTOM_REAL), allocatable :: cov(:) ! Local on one CPU.
!
!    ! Data arrays for local ADMM constraints.
!    integer :: nlithos
!    real(kind=CUSTOM_REAL), allocatable :: min_local_bound(:, :)
!    real(kind=CUSTOM_REAL), allocatable :: max_local_bound(:, :)
!    real(kind=CUSTOM_REAL), allocatable :: local_bound_constraints_weight(:)
!
!    ! Full grid.
!    type(t_grid) :: grid_full
!
!    ! Total number of model parameters.
!    integer :: nelements_total
!    ! Local number of model parameters (on current CPU).
!    integer :: nelements
!
!  contains
!    private
!
!    procedure, public, pass :: initialize => model_initialize
!    procedure, public, pass :: init_grid => model_init_grid
!    procedure, public, pass :: allocate_bound_arrays => model_allocate_bound_arrays
!
!    procedure, public, pass :: distribute => model_distribute
!
!    procedure, public, pass :: get_Xmin => model_get_Xmin
!    procedure, public, pass :: get_Xmax => model_get_Xmax
!    procedure, public, pass :: get_Ymin => model_get_Ymin
!    procedure, public, pass :: get_Ymax => model_get_Ymax
!
!  end type t_model_base
!
!contains
!
!!================================================================================================
!! Initialization.
!!================================================================================================
!subroutine model_initialize(this, nelements, myrank, nbproc)
!  class(t_model_base), intent(inout) :: this
!  integer, intent(in) :: nelements, myrank, nbproc
!
!  integer :: ierr
!  type(t_parallel_tools) :: pt
!
!  this%nelements = nelements
!  this%nelements_total = pt%get_total_number_elements(nelements, myrank, nbproc)
!
!  ierr = 0
!
!  if (.not. allocated(this%val)) allocate(this%val(this%nelements), source=0._CUSTOM_REAL, stat=ierr)
!  if (.not. allocated(this%val_full)) allocate(this%val_full(this%nelements_total), source=0._CUSTOM_REAL, stat=ierr)
!
!  if (.not. allocated(this%cov)) allocate(this%cov(this%nelements), source=1._CUSTOM_REAL, stat=ierr)
!
!  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_initialize!", myrank, ierr)
!
!end subroutine model_initialize
!
!!================================================================================================
!! Allocate bound arrays (used for ADMM).
!!================================================================================================
!subroutine model_allocate_bound_arrays(this, nlithos, myrank)
!  class(t_model_base), intent(inout) :: this
!  integer, intent(in) :: nlithos, myrank
!
!  integer :: ierr
!  ierr = 0
!
!  this%nlithos = nlithos
!
!  if (.not. allocated(this%min_local_bound)) allocate(this%min_local_bound(this%nelements, this%nlithos), &
!    source=0._CUSTOM_REAL, stat=ierr)
!  if (.not. allocated(this%max_local_bound)) allocate(this%max_local_bound(this%nelements, this%nlithos), &
!    source=0._CUSTOM_REAL, stat=ierr)
!  if (.not. allocated(this%local_bound_constraints_weight)) allocate(this%local_bound_constraints_weight(this%nelements), &
!    source=0._CUSTOM_REAL, stat=ierr)
!
!  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_allocate_bound_arrays!", myrank, ierr)
!
!end subroutine model_allocate_bound_arrays
!
!!================================================================================================
!! Initialization of the model grid.
!!================================================================================================
!subroutine model_init_grid(this, nx, ny, nz, myrank)
!  class(t_model_base), intent(inout) :: this
!  integer, intent(in) :: nx, ny, nz
!  integer, intent(in) :: myrank
!
!  call this%grid_full%initialize(this%nelements_total, nx, ny, nz, .true., myrank)
!
!end subroutine model_init_grid
!
!!=================================================================================
!! Distribute the grid and the prior model among CPUs.
!!=================================================================================
!subroutine model_distribute(this, myrank, nbproc)
!  class(t_model_base), intent(inout) :: this
!  integer, intent(in) :: myrank, nbproc
!
!  ! Displacement for mpi_scatterv.
!  integer :: displs(nbproc)
!  ! The number of elements on every CPU for mpi_scatterv.
!  integer :: nelements_at_cpu(nbproc)
!  integer :: ierr
!  type(t_parallel_tools) :: pt
!
!  ! Partitioning for MPI_Scatterv.
!  call pt%get_mpi_partitioning(this%nelements, displs, nelements_at_cpu, myrank, nbproc)
!
!  call MPI_Scatterv(this%val_full, nelements_at_cpu, displs, CUSTOM_MPI_TYPE, &
!                    this%val, this%nelements, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
!
!  if (ierr /= 0) call exit_MPI("Error in MPI_Scatterv in model_distribute!", myrank, ierr)
!
!end subroutine model_distribute
!
!!================================================================================================
!! Get the minimum X-coordinate of the model grid.
!!================================================================================================
!pure function model_get_Xmin(this) result(res)
!  class(t_model_base), intent(in) :: this
!  real(kind=CUSTOM_REAL) :: res
!
!  res = minval(this%grid_full%X1)
!end function model_get_Xmin
!
!!================================================================================================
!! Get the maximum X-coordinate of the model grid.
!!================================================================================================
!pure function model_get_Xmax(this) result(res)
!  class(t_model_base), intent(in) :: this
!  real(kind=CUSTOM_REAL) :: res
!
!  res = maxval(this%grid_full%X2)
!end function model_get_Xmax
!
!!================================================================================================
!! Get the minimum Y-coordinate of the model grid.
!!================================================================================================
!pure function model_get_Ymin(this) result(res)
!  class(t_model_base), intent(in) :: this
!  real(kind=CUSTOM_REAL) :: res
!
!  res = minval(this%grid_full%Y1)
!end function model_get_Ymin
!
!!================================================================================================
!! Get the maximum Y-coordinate of the model grid.
!!================================================================================================
!pure function model_get_Ymax(this) result(res)
!  class(t_model_base), intent(in) :: this
!  real(kind=CUSTOM_REAL) :: res
!
!  res = maxval(this%grid_full%Y2)
!end function model_get_Ymax
!
!end module model_base
