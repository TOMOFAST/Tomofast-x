
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
! A class to perform parallel joint inversion of two data sets.
! Works with sparse matrices that are stored using Compressed Sparse Row (CSR) format.
! Calculates the model update (change).
! Uses an object of type t_parameters_inversion to obtain the input parameters and data arrays.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!===============================================================================================
module joint_inverse_problem

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use parameters_inversion
  use inversion_arrays
  use cross_gradient
  use sparse_matrix
  use vector
  use lsqr_solver
  use damping
  use model
  use admm_method
  use clustering
  use damping_gradient
  use parallel_tools
  use wavelet_utils
  use costs
  use grid
  use paraview

  implicit none

  private

  !---------------------------------------------------------------------
  ! Main inversion class.
  type, public :: t_joint_inversion
    private

    ! Matrix that stores sensitivity kernel(s).
    type(t_sparse_matrix), public :: matrix_sensit

    ! Matrix that stores constraints.
    type(t_sparse_matrix), public :: matrix_cons

    ! Right hand side (corresponding to the matrix).
    real(kind=CUSTOM_REAL), allocatable :: b_RHS(:)

    ! Auxiliary arrays for the ADMM method.
    type(t_real1d) :: x0_ADMM(2)

    integer :: nelements_total

    ! Flags for adding various constraints.
    logical :: add_damping(2)
    logical :: add_damping_gradient(2)
    logical :: add_admm(2)
    logical, public :: add_cross_grad
    logical, public :: add_clustering

    ! Cross gradient data.
    type(t_cross_gradient) :: cross_grad

    ! A compact grid for calculating gradients.
    type(t_grad_grid) :: grad_grid

    ! ADMM method (stores auxiliary arrays).
    type(t_admm_method) :: admm_method(2)

    ! ADMM term cost.
    real(kind=CUSTOM_REAL) :: admm_cost(2)

    ! Damping gradient term cost (in each direction and for each joint problem).
    real(kind=CUSTOM_REAL) :: damping_gradient_cost(3, 2)

    ! Clustering constraints.
    type(t_clustering), public :: clustering

    integer :: ndata_lines

    ! Flag defining if we solve inversion in the wavelet domain for the model vector.
    logical :: WAVELET_DOMAIN

  contains
    private

    procedure, public, pass :: initialize => joint_inversion_initialize
    procedure, public, pass :: initialize2 => joint_inversion_initialize2
    procedure, public, pass :: reset => joint_inversion_reset
    procedure, public, pass :: solve => joint_inversion_solve

    procedure, public, pass :: get_cross_grad => joint_inversion_get_cross_grad
    procedure, public, pass :: get_clustering => joint_inversion_get_clustering
    procedure, public, pass :: get_cross_grad_cost => joint_inversion_get_cross_grad_cost
    procedure, public, pass :: get_clustering_cost => joint_inversion_get_clustering_cost
    procedure, public, pass :: get_admm_cost => joint_inversion_get_admm_cost
    procedure, public, pass :: get_damping_gradient_cost => joint_inversion_get_damping_gradient_cost

    procedure, private, pass :: add_cross_grad_constraints => joint_inversion_add_cross_grad_constraints
    procedure, private, pass :: add_clustering_constraints => joint_inversion_add_clustering_constraints

    procedure, public, nopass :: calculate_matrix_partitioning => joint_inversion_calculate_matrix_partitioning

  end type t_joint_inversion

contains

!=====================================================================================
! Initialize joint inversion.
!=====================================================================================
subroutine joint_inversion_initialize(this, par, nnz_sensit, myrank)
  class(t_joint_inversion), intent(inout) :: this
  type(t_parameters_inversion), intent(in) :: par
  integer(kind=8), intent(in) :: nnz_sensit
  integer, intent(in) :: myrank

  integer :: i

  do i = 1, 2
    if (par%alpha(i) == 0.d0 .or. par%problem_weight(i) == 0.d0) then
      this%add_damping(i) = .false.
    else
      this%add_damping(i) = .true.
    endif
  enddo

  do i = 1, 2
    if (par%beta(i) == 0.d0 .or. par%problem_weight(i) == 0.d0) then
      this%add_damping_gradient(i) = .false.
    else
      this%add_damping_gradient(i) = .true.
    endif
  enddo

  if (par%cross_grad_weight == 0.d0) then
    this%add_cross_grad = .false.
  else
    this%add_cross_grad = .true.
  endif

  if (par%clustering_weight_glob(1) == 0.d0 .and. par%clustering_weight_glob(2) == 0.d0) then
    this%add_clustering = .false.
  else
    this%add_clustering = .true.
  endif

  this%add_admm = .false.
  if (par%admm_type > 0) then
    do i = 1, 2
      if (par%problem_weight(i) /= 0.d0) then
        this%add_admm(i) = .true.
      endif
    enddo
  endif

  this%nelements_total = par%nelements_total

  if (myrank == 0) then
    print *, 'add_damping =', this%add_damping
    print *, 'add_admm =', this%add_admm
    print *, 'add_cross_grad, add_clustering =', this%add_cross_grad, this%add_clustering
    print *, 'myrank, ndata1, nelements1 =', myrank, par%ndata(1), par%nelements
    print *, 'myrank, ndata2, nelements2 =', myrank, par%ndata(2), par%nelements
  endif

  if (this%add_clustering) then
    ! Allocate memory for clustering constraints.
    call this%clustering%initialize(par%nelements_total, par%nelements, par%clustering_weight_glob, &
                                    par%nclusters, par%clustering_opt_type, par%clustering_constraints_type, myrank)
    call this%clustering%read_mixtures(par%mixture_file, par%cell_weights_file, myrank)
  endif

  this%admm_cost = 0.d0
  this%damping_gradient_cost = 0.d0

  ! Determine if we perform inversion using wavelet domain for the model update vector.
  this%WAVELET_DOMAIN = .true.
  if (this%add_cross_grad .or. this%add_clustering .or. &
      this%add_damping_gradient(1) .or. this%add_damping_gradient(2) .or. &
      par%norm_power /= 2.d0 .or. &
      ! Essentially we care only about the local admm weights (i.e., when the bound_weight /= 1).
      par%admm_bound_type /= 1 .or. &
      par%apply_local_damping_weight > 0) then
    this%WAVELET_DOMAIN = .false.
  endif

  if (myrank == 0) print *, 'WAVELET_DOMAIN =', this%WAVELET_DOMAIN

  ! Calculate the number of sensitivity matrix rows.
  this%ndata_lines = 0
  do i = 1, 2
    if (par%problem_weight(i) /= 0.d0) then
      this%ndata_lines = this%ndata_lines + par%ndata(i) * par%ndata_components(i)
    endif
  enddo

  !-------------------------------------------------------------------------------------------
  ! SENSITIVITY MATRIX MEMORY ALLOCATION.
  !-------------------------------------------------------------------------------------------
  call this%matrix_sensit%initialize(this%ndata_lines, &
                                     2 * par%nmodel_components * par%nelements, nnz_sensit, myrank, 0)

end subroutine joint_inversion_initialize

!==================================================================================================
! Remaining memory allocations.
! Split memory allocations to reduce the total memory footprint.
! Call this after reading the sensitivity kernel where some temporary buffer arrays are allocated.
!==================================================================================================
subroutine joint_inversion_initialize2(this, par, arr, model, myrank, nbproc)
  class(t_joint_inversion), intent(inout) :: this
  type(t_parameters_inversion), intent(in) :: par
  type(t_inversion_arrays), intent(in) :: arr(2)
  type(t_model), intent(in) :: model(2)
  integer, intent(in) :: myrank, nbproc

  integer :: ierr
  real(kind=CUSTOM_REAL) :: mem, mem_loc
  integer :: i, k, nl, nl_empty, nl_empty_loc
  integer(kind=8) :: nnz
  real(kind=CUSTOM_REAL), allocatable :: b_dummy(:)
  character(len=256) :: filename
  integer :: nelements_total

  if (this%add_cross_grad) then
    call this%cross_grad%initialize(par%nx, par%ny, par%nz, par%nelements, par%keep_model_constant, &
                                    par%vec_field_type, par%vec_field_file, myrank)

    if (myrank == 0 .and. par%vec_field_type > 0) then
      ! Write vtk model for visualisation of the vector field.
      filename = "cross_grad_vec_field.vtk"
      nelements_total = model(1)%grid_full%nx * model(1)%grid_full%ny * model(1)%grid_full%nz
      call visualisation_paraview_struct_grid(filename, myrank, nelements_total, 3, this%cross_grad%vec_field, &
                                   model(1)%grid_full%X1, model(1)%grid_full%Y1, model(1)%grid_full%Z1, &
                                   model(1)%grid_full%X2, model(1)%grid_full%Y2, model(1)%grid_full%Z2, &
                                   model(1)%grid_full%nx, model(1)%grid_full%ny, model(1)%grid_full%nz, &
                                   1, model(1)%grid_full%nx, 1, model(1)%grid_full%ny, 1, model(1)%grid_full%nz, &
                                   .true., 1.d0, 'Vector_field')
    endif
  endif

  !--------------------------------------------------------------------------------------
  ! Initialize the gradient grid.
  !--------------------------------------------------------------------------------------
  if (this%add_cross_grad .or. &
      this%add_damping_gradient(1) .or. this%add_damping_gradient(2)) then
    call this%grad_grid%init(model(merge(1, 2, par%problem_weight(1) /= 0.d0))%grid_full, myrank)
  endif

  !--------------------------------------------------------------------------------------
  ! Calcualte nl and nnz for the constraints matrix.
  !--------------------------------------------------------------------------------------
  nnz = 0
  nl = 0
  nl_empty = 0
  nl_empty_loc = 0

  do i = 1, 2
    if (this%add_damping(i)) then
      do k = 1, par%nmodel_components
        nl = nl + par%nelements_total
        nnz = nnz + par%nelements
        nl_empty = nl_empty + par%nelements_total - par%nelements
      enddo
    endif
  enddo

  do i = 1, 2
    if (this%add_damping_gradient(i)) then
      do k = 1, par%nmodel_components
        nl = nl + 3 * par%nelements_total
        nnz = nnz + 6 * par%nelements

        ! Note: we increase a multiplier from 3 to 12 to account for elements from other ranks.
        nl_empty_loc = 3 * par%nelements_total - 12 * par%nelements
        if (nl_empty_loc > 0) then
          nl_empty = nl_empty + nl_empty_loc
        endif
      enddo
    endif
  enddo

  if (this%add_cross_grad) then
    allocate(b_dummy(1))
    ! Calculate the cross-gradient without adding it to the matrix and RHS to obtain accurate nnz and nl_nonempty.
    call this%cross_grad%calculate(model(1)%val_full(:, 1), model(2)%val_full(:, 1), &
                                   this%grad_grid, &
                                   arr(1)%column_weight, arr(2)%column_weight, &
                                   this%matrix_cons, b_dummy, .false., par%derivative_type, &
                                   par%cross_grad_weight, myrank, nbproc)

    nl = nl + 3 * par%nelements_total
    nnz = nnz + this%cross_grad%nnz

    nl_empty_loc = 3 * par%nelements_total - this%cross_grad%nl_nonempty
    if (nl_empty_loc > 0) then
      nl_empty = nl_empty + nl_empty_loc
    endif
  endif

  do i = 1, 2
    if (this%add_admm(i)) then
      nl = nl + par%nelements_total
      nnz = nnz + par%nelements
      nl_empty = nl_empty + par%nelements_total - par%nelements

      call this%admm_method(i)%initialize(par%nelements, myrank)
    endif
  enddo

  if (this%add_clustering) then
    nl = nl + 2 * par%nelements_total
    nnz = nnz + 2 * par%nelements
  endif

  !--------------------------------------------------------------------------------------
  ! Constraints matrix allocation.
  !--------------------------------------------------------------------------------------
  call this%matrix_cons%initialize(nl, 2 * par%nmodel_components * par%nelements, nnz, myrank, nl_empty)

  !-------------------------------------------------------------------------------------------------
  ! Allocate the right-hand side.
  !-------------------------------------------------------------------------------------------------
  ! Total number of matrix rows.
  nl = this%matrix_sensit%get_total_row_number() + this%matrix_cons%get_total_row_number()

  mem_loc = kind(this%b_RHS) * nl

  call mpi_allreduce(mem_loc, mem, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (myrank == 0) print *, "Allocating the RHS, memory (GB) =", mem / 1024**3

  allocate(this%b_RHS(nl), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in joint_inversion_initialize2!", myrank, ierr)

  !-------------------------------------------------------------------------------------------------
  do i = 1, 2
    if (this%add_admm(i)) then
      allocate(this%x0_ADMM(i)%val(par%nelements), source=0._CUSTOM_REAL, stat=ierr)
    endif
  enddo

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in joint_inversion_initialize2!", myrank, ierr)

end subroutine joint_inversion_initialize2

!=====================================================================================
! Resets the joint inversion.
!=====================================================================================
subroutine joint_inversion_reset(this, myrank)
  class(t_joint_inversion), intent(inout) :: this
  integer, intent(in) :: myrank

  if (myrank >= 0) continue

  call this%matrix_cons%reset()
  this%b_RHS = 0._CUSTOM_REAL

end subroutine joint_inversion_reset

!=====================================================================================
! Calculate the righ-hand side (corresponding to data misfit).
! Utilize the automatic conversion of the 2D to 1D array.
!=====================================================================================
pure subroutine calculate_b_RHS(line_start, line_end, problem_weight, residuals, b_RHS)
  integer, intent(in) :: line_start, line_end
  real(kind=CUSTOM_REAL), intent(in) :: problem_weight
  real(kind=CUSTOM_REAL), intent(in) :: residuals(line_end - line_start + 1)
  real(kind=CUSTOM_REAL), intent(out) :: b_RHS(:)

  b_RHS(line_start:line_end) = problem_weight * residuals

end subroutine calculate_b_RHS

!================================================================================================
! Joint inversion of two problems.
! It is also used for single inversions.
!================================================================================================
subroutine joint_inversion_solve(this, par, arr, model, delta_model, memory, myrank, nbproc)
  class(t_joint_inversion), intent(inout) :: this
  type(t_parameters_inversion), intent(in) :: par
  type(t_inversion_arrays), intent(in) :: arr(2)
  type(t_model), intent(inout) :: model(2)
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(out) :: memory

  real(kind=CUSTOM_REAL), intent(out) :: delta_model(par%nelements, par%nmodel_components, 2)

  type(t_damping) :: damping
  type(t_damping_gradient) :: damping_gradient
  integer :: i, j, k
  integer :: line_start(2), line_end(2), param_shift(2)
  integer :: damping_param_shift
  real(kind=CUSTOM_REAL) :: cost
  real(kind=CUSTOM_REAL) :: norm_power
  integer :: lc
  logical :: SOLVE_PROBLEM(2)
  integer :: buffer_index

  ! The number of times this subroutine has been called.
  ! This is effectively the major loop iteration number.
  integer, save :: ncalls = 0

  ncalls = ncalls + 1

  ! Store the iteration number.
  this%matrix_sensit%tag = ncalls

  do i = 1, 2
    SOLVE_PROBLEM(i) = (par%problem_weight(i) /= 0.d0)
  enddo

  call this%calculate_matrix_partitioning(par, line_start, line_end, param_shift)

  ! Starting line for constraints.
  lc = this%ndata_lines + 1

  ! ***** Data misfit and damping and damping gradient  *****

  this%damping_gradient_cost = 0.d0

  ! Loop over joint problems.
  do i = 1, 2

    ! Skip the problem with zero problem weight (single inversion).
    if (.not. SOLVE_PROBLEM(i)) then
      cycle
    endif

    if (myrank == 0) print *, 'Adding joint problem #', i, ' weight =', par%problem_weight(i)

    ! Adding the right-hand side (corresponding to data misfit).
    call calculate_b_RHS(line_start(i), line_end(i), par%problem_weight(i), arr(i)%residuals, this%b_RHS)

    if (this%add_damping(i)) then
      if (myrank == 0) print *, 'adding damping with alpha =', par%alpha(i)

      call damping%initialize(par%nelements, par%alpha(i), par%problem_weight(i), par%norm_power, &
                              par%compression_type, par%nx, par%ny, par%nz)

      ! Adding model damping for each component.
      do k = 1, par%nmodel_components
        damping_param_shift = param_shift(i) + (k - 1) * par%nelements
        call damping%add(this%matrix_cons, size(this%b_RHS(lc:)), this%b_RHS(lc:), arr(i)%column_weight, &
                         model(i)%val(:, k), model(i)%val_prior(:, k), damping_param_shift, &
                         this%WAVELET_DOMAIN, myrank, nbproc, model(i)%damping_weight)

        if (myrank == 0) print *, 'model damping cost = ', damping%get_cost()
      enddo
    endif

    if (this%add_damping_gradient(i)) then
      if (myrank == 0) print *, 'adding damping_gradient with beta =', par%beta(i)

      call damping_gradient%initialize(par%beta(i), par%problem_weight(i), par%nx, par%ny, par%nz, par%nelements)

      ! Update the full model.
      call model(i)%update_full(.true., myrank, nbproc)

      ! Adding damping gradient for each component.
      do k = 1, par%nmodel_components
        damping_param_shift = param_shift(i) + (k - 1) * par%nelements
        do j = 1, 3 ! j is direction (1 = x, 2 = y, 3 = z).
          call damping_gradient%add(model(i), this%grad_grid, &
                                    arr(i)%column_weight, model(i)%damping_grad_weight(:, j), &
                                    this%matrix_cons, size(this%b_RHS(lc:)), this%b_RHS(lc:), &
                                    damping_param_shift, j, k, myrank, nbproc)

          cost = damping_gradient%get_cost()
          ! Sum the cost over the model components.
          this%damping_gradient_cost(j, i) = this%damping_gradient_cost(j, i) + cost
        enddo
      enddo
      if (myrank == 0) then
        print *, 'damping gradient cost x =', this%damping_gradient_cost(1, i)
        print *, 'damping gradient cost y =', this%damping_gradient_cost(2, i)
        print *, 'damping gradient cost z =', this%damping_gradient_cost(3, i)
      endif
    endif
  enddo

  ! ***** ADMM method *****
  do i = 1, 2
    if (this%add_admm(i)) then
      call this%admm_method(i)%iterate_admm_arrays(model(i)%nlithos, model(i)%min_bound, model(i)%max_bound, &
                                                   model(i)%val, this%x0_ADMM(i)%val)

      ! Use the L2 norm for the ADMM constraints.
      norm_power = 2.0d0

      call damping%initialize(par%nelements, par%rho_ADMM(i), par%problem_weight(i), norm_power, &
                              par%compression_type, par%nx, par%ny, par%nz)

      call damping%add(this%matrix_cons, size(this%b_RHS(lc:)), this%b_RHS(lc:), arr(i)%column_weight, &
                       model(i)%val(:, 1), this%x0_ADMM(i)%val, param_shift(i), &
                       this%WAVELET_DOMAIN, myrank, nbproc, model(i)%bound_weight)

      ! Calculate the ADMM cost in parallel.
      call calculate_cost(par%nelements, this%admm_method(i)%z, model(i)%val(:, 1), this%admm_cost(i), .true., nbproc)

      if (myrank == 0) print *, "ADMM cost |x - z| / |z| =", this%admm_cost(i)
    endif
  enddo

  ! ***** Cross-gradient *****

  if (this%add_cross_grad) then
    call this%add_cross_grad_constraints(par, arr, model, par%derivative_type, myrank, nbproc)
  endif

  ! ***** Clustering *****

  if (this%add_clustering) then
    call this%add_clustering_constraints(arr, model, myrank, nbproc)
  endif

  !-------------------------------------------------------------------------------------
  ! Parallel sparse inversion.
  !-------------------------------------------------------------------------------------
  call this%matrix_cons%finalize(myrank)

  delta_model = 0._CUSTOM_REAL
  if (par%method == 1) then
    call lsqr_solve_sensit(size(this%b_RHS), size(delta_model), par%niter, par%rmin, par%gamma, par%target_misfit, &
                           this%matrix_sensit, this%matrix_cons, this%b_RHS, delta_model, &
                           SOLVE_PROBLEM, par%nelements, par%nx, par%ny, par%nz, par%nmodel_components, &
                           par%compression_type, this%WAVELET_DOMAIN, memory, myrank, nbproc)
  else
    call exit_MPI("Unknown solver type!", myrank, 0)
  endif

  !-------------------------------------------------------------------------------------
  ! Convert the model update to original variables.
  !-------------------------------------------------------------------------------------
  if (par%compression_type > 0 .and. this%WAVELET_DOMAIN) then
    ! Index for the full model buffer.
    buffer_index = merge(1, 2, SOLVE_PROBLEM(1) .eqv. .true.)

    ! Convert the model update from the wavelet domain.
    call apply_wavelet_transform(par%nelements, par%nx, par%ny, par%nz, par%nmodel_components, &
                                 delta_model, model(buffer_index)%val_full(:, 1), &
                                 .false., par%compression_type, 2, SOLVE_PROBLEM, myrank, nbproc)
  endif

  ! Apply the inverse depth weighting.
  if (SOLVE_PROBLEM(1)) call rescale_model(par%nelements, par%nmodel_components, delta_model(:, :, 1), arr(1)%column_weight)
  if (SOLVE_PROBLEM(2)) call rescale_model(par%nelements, par%nmodel_components, delta_model(:, :, 2), arr(2)%column_weight)

end subroutine joint_inversion_solve

!=======================================================================================================
! Adds the cross-gradient constraints to the least square system.
!=======================================================================================================
subroutine joint_inversion_add_cross_grad_constraints(this, par, arr, model, der_type, myrank, nbproc)
  class(t_joint_inversion), intent(inout) :: this
  type(t_parameters_inversion), intent(in) :: par
  type(t_inversion_arrays), intent(in) :: arr(2)
  type(t_model), intent(inout) :: model(2)
  integer, intent(in) :: der_type
  integer, intent(in) :: myrank, nbproc

  type(t_vector) :: cost
  integer :: lc

  if (myrank == 0) print *, 'Calculating cross gradients. der_type =', der_type

  ! Update the full models.
  call model(1)%update_full(.true., myrank, nbproc)
  call model(2)%update_full(.true., myrank, nbproc)

  ! Starting line for constraints.
  lc = this%ndata_lines + 1

  call this%cross_grad%calculate(model(1)%val_full(:, 1), model(2)%val_full(:, 1), &
                                 this%grad_grid, &
                                 arr(1)%column_weight, arr(2)%column_weight, &
                                 this%matrix_cons, this%b_RHS(lc:), .true., der_type, &
                                 par%cross_grad_weight, myrank, nbproc)

  cost = this%cross_grad%get_cost()

  if (myrank == 0) print *, 'cross-grad cost = ', cost%x, cost%y, cost%z

end subroutine joint_inversion_add_cross_grad_constraints

!=======================================================================================================
! Adds the clustering constraints to the least square system.
!=======================================================================================================
subroutine joint_inversion_add_clustering_constraints(this, arr, model, myrank, nbproc)
  class(t_joint_inversion), intent(inout) :: this
  type(t_inversion_arrays), intent(in) :: arr(2)
  type(t_model), intent(inout) :: model(2)
  integer, intent(in) :: myrank, nbproc

  integer :: i, lc

  ! Starting line for constraints.
  lc = this%ndata_lines + 1

  do i = 1, 2
    call this%clustering%add(model(1), model(2), arr(1)%column_weight, arr(2)%column_weight, &
                             this%matrix_cons, this%b_RHS(lc:), i, myrank, nbproc)

    if (myrank == 0) print *, 'clustering term', i, 'cost = ', this%clustering%get_cost(i)
  enddo

end subroutine joint_inversion_add_clustering_constraints

!==================================================================================================
! Returns the cross-gradient vector magnitude at every model pixel.
!==================================================================================================
pure subroutine joint_inversion_get_cross_grad(this, cross_grad)
  class(t_joint_inversion), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(out) :: cross_grad(this%nelements_total)

  call this%cross_grad%get_magnitude(cross_grad)

end subroutine joint_inversion_get_cross_grad

!==================================================================================================
! Returns the clustering probabilities for every model pixel.
!==================================================================================================
pure subroutine joint_inversion_get_clustering(this, res)
  class(t_joint_inversion), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(out) :: res(this%nelements_total)

  res = this%clustering%get_probabilities()

end subroutine joint_inversion_get_clustering

!==================================================================================================
! Returns the cross-gradient cost (for each direction).
!==================================================================================================
pure function joint_inversion_get_cross_grad_cost(this) result(res)
  class(t_joint_inversion), intent(in) :: this
  type(t_vector) :: res

  if (this%add_cross_grad) then
    res = this%cross_grad%get_cost()
  else
    res = t_vector(0.d0, 0.d0, 0.d0)
  endif

end function joint_inversion_get_cross_grad_cost

!==================================================================================================
! Returns the clustering cost.
!==================================================================================================
pure function joint_inversion_get_clustering_cost(this, problem_type) result(res)
  class(t_joint_inversion), intent(in) :: this
  integer, intent(in) :: problem_type
  real(kind=CUSTOM_REAL) :: res

  if (this%add_clustering) then
    res = this%clustering%get_cost(problem_type)
  else
    res = 0.d0
  endif

end function joint_inversion_get_clustering_cost

!==================================================================================================
! Returns the ADMM cost.
!==================================================================================================
pure function joint_inversion_get_admm_cost(this, problem_type) result(res)
  class(t_joint_inversion), intent(in) :: this
  integer, intent(in) :: problem_type
  real(kind=CUSTOM_REAL) :: res

  res = this%admm_cost(problem_type)

end function joint_inversion_get_admm_cost

!==================================================================================================
! Returns the damping gradient cost.
!==================================================================================================
pure function joint_inversion_get_damping_gradient_cost(this) result(res)
  class(t_joint_inversion), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res(3, 2)

  res = this%damping_gradient_cost

end function joint_inversion_get_damping_gradient_cost

!==================================================================================================
! Calculates the matrix pertitioning.
!==================================================================================================
subroutine joint_inversion_calculate_matrix_partitioning(par, line_start, line_end, param_shift)
  type(t_parameters_inversion), intent(in) :: par
  integer, intent(out) :: line_start(2), line_end(2), param_shift(2)

  logical :: SOLVE_PROBLEM(2)
  integer :: i

  do i = 1, 2
    SOLVE_PROBLEM(i) = (par%problem_weight(i) /= 0.d0)
  enddo

  param_shift(1) = 0
  param_shift(2) = par%nelements * par%nmodel_components

  line_start = 0
  line_end = 0

  if (SOLVE_PROBLEM(1)) then
    line_start(1) = 1
    line_end(1) = par%ndata(1) * par%ndata_components(1)
  endif

  if (SOLVE_PROBLEM(2)) then
    line_start(2) = line_end(1) + 1
    line_end(2) = line_end(1) + par%ndata(2) * par%ndata_components(2)
  endif

end subroutine joint_inversion_calculate_matrix_partitioning

end module joint_inverse_problem
