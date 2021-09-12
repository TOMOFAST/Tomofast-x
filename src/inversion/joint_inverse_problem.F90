
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
  use inverse_problem
  use cross_gradient
  use sparse_matrix
  use vector
  use lsqr_solver
  use sensitivity_matrix
  use damping
  use model
  use method_of_weights
  use sca_solver
  use admm_method
  use clustering
  use damping_gradient
  use parallel_tools
  use wavelet_transform

  implicit none

  private

  logical, parameter :: NORMALIZE_MATRIX_COLUMNS = .false.

  !---------------------------------------------------------------------
  ! Main inversion class.
  type, public :: t_joint_inversion
    private

    ! Matrix that stores both of the joint problems and cross-gradient part.
    type(t_sparse_matrix) :: matrix
    ! Right hand side (corresponding to the matrix).
    real(kind=CUSTOM_REAL), allocatable :: b_RHS(:)

    ! Matrix of the cross-gradient part.
    type(t_sparse_matrix) :: matrix_B
    ! Right hand side for the cross-gradient part.
    real(kind=CUSTOM_REAL), allocatable :: d_RHS(:)

    ! Matrix columns norm.
    real(kind=CUSTOM_REAL), allocatable :: column_norm(:)

    ! Auxiliary arrays for the ADMM method.
    real(kind=CUSTOM_REAL), allocatable :: x0_ADMM(:)
    real(kind=CUSTOM_REAL), allocatable :: weight_ADMM(:)

    integer :: nelements_total

    ! Flags to switch off the use of some terms (for debugging).
    logical :: add_damping(2)
    logical :: add_damping_gradient(2)
    logical, public :: add_cross_grad
    logical, public :: add_clustering

    ! Cross gradient data.
    type(t_cross_gradient) :: cross_grad

    ! ADMM method (stores data arrays).
    ! NOTE: gcc 4.9.2 has problems with this array.
    ! E.g. if the destructor in t_admm_method is commented, the code won't compile (gfc_conv_descriptor_data_get error).
    ! With destructor available it leads to memory segfault during the finalization.
    ! The problem is probably solved in gcc 5.3.0 according to comments from gcc forum.
    ! Avoid using this array for now for back compatibility.
    !type(t_admm_method) :: admm_method(2)
    type(t_admm_method) :: admm_method_1
    type(t_admm_method) :: admm_method_2

    ! Clustering constraints.
    type(t_clustering), public :: clustering

  contains
    private

    procedure, public, pass :: initialize => joint_inversion_initialize
    procedure, public, pass :: reset => joint_inversion_reset
    procedure, public, pass :: solve => joint_inversion_solve

    procedure, public, pass :: get_cross_grad => joint_inversion_get_cross_grad
    procedure, public, pass :: get_clustering => joint_inversion_get_clustering
    procedure, public, pass :: get_cross_grad_cost => joint_inversion_get_cross_grad_cost
    procedure, public, pass :: get_clustering_cost => joint_inversion_get_clustering_cost

    procedure, private, pass :: add_cross_grad_constraints => joint_inversion_add_cross_grad_constraints
    procedure, private, pass :: add_clustering_constraints => joint_inversion_add_clustering_constraints

  end type t_joint_inversion

contains

!=================================================================================
! Initialize joint inversion.
!=================================================================================
subroutine joint_inversion_initialize(this, par, nnz_sensit, myrank)
  class(t_joint_inversion), intent(inout) :: this
  type(t_parameters_inversion), intent(in) :: par
  integer, intent(in) :: nnz_sensit, myrank
  integer :: ierr
  integer :: i, nl, nnz

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

  this%nelements_total = par%nelements_total

  if (myrank == 0) then
    print *, 'add_damping1 =', this%add_damping(1)
    print *, 'add_damping2 =', this%add_damping(2)
    print *, 'add_cross_grad, add_clustering =', this%add_cross_grad, this%add_clustering
    print *, 'myrank, ndata1, nelements1 =', myrank, par%ndata(1), par%nelements
    print *, 'myrank, ndata2, nelements2 =', myrank, par%ndata(2), par%nelements
  endif

  call this%cross_grad%initialize(par%nx, par%ny, par%nz, par%nelements, myrank)

  call this%clustering%initialize(par%nelements_total, par%nelements, par%clustering_weight_glob, &
                                  par%nclusters, par%clustering_opt_type, par%clustering_constraints_type, myrank)

  if (this%add_clustering) then
    call this%clustering%read_mixtures(par%mixture_file, par%cell_weights_file, myrank)
  endif

  ! Calculate number of matrix rows and (approx.) non-zero elements.
  nl = 0
  do i = 1, 2
    if (par%problem_weight(i) /= 0.d0) then
      nl = nl + par%ndata(i)
    endif
  enddo

  nnz = nnz_sensit

  do i = 1, 2
    if (this%add_damping(i)) then
      nl = nl + par%nelements_total
      nnz = nnz + par%nelements
    endif
  enddo

  do i = 1, 2
    if (this%add_damping_gradient(i)) then
      nl = nl + 3 * par%nelements_total
      nnz = nnz + 3 * 2 * par%nelements
    endif
  enddo

  if (this%add_cross_grad) then
    if (par%derivative_type /= 4) then
      nl = nl + 3 * par%nelements_total
      nnz = nnz + this%cross_grad%get_num_elements(par%derivative_type)
    else
      nl = nl + 6 * par%nelements_total
      nnz = nnz + this%cross_grad%get_num_elements(1)
      nnz = nnz + this%cross_grad%get_num_elements(2)
    endif
  endif

  if (par%admm_type > 0) then
    ! Add only one ADMM constraint at a time.
    nl = nl + 1 * par%nelements_total
    nnz = nnz + 1 * par%nelements

    call this%admm_method_1%initialize(par%nelements, myrank)
    call this%admm_method_2%initialize(par%nelements, myrank)
  endif

  if (this%add_clustering) then
    nl = nl + 2 * par%nelements_total
    nnz = nnz + 2 * par%nelements
  endif

  call this%matrix%initialize(nl, nnz, myrank)
  call this%matrix_B%initialize(3 * par%nelements_total, this%cross_grad%get_num_elements(par%derivative_type), myrank)

  call this%matrix%allocate_variance_array(2 * par%nelements, myrank)

  ierr = 0

  ! Memory allocation.
  if (.not. allocated(this%b_RHS)) &
    allocate(this%b_RHS(this%matrix%get_total_row_number()), source=0._CUSTOM_REAL, stat=ierr)

  if (.not. allocated(this%d_RHS)) &
    allocate(this%d_RHS(3 * par%nelements_total), source=0._CUSTOM_REAL, stat=ierr)

  if (.not. allocated(this%column_norm)) &
    allocate(this%column_norm(2 * par%nelements), source=1._CUSTOM_REAL, stat=ierr)

  if (par%admm_type > 0) then
    if (.not. allocated(this%x0_ADMM)) &
      allocate(this%x0_ADMM(par%nelements), source=0._CUSTOM_REAL, stat=ierr)

    if (.not. allocated(this%weight_ADMM)) &
      allocate(this%weight_ADMM(par%nelements), source=1._CUSTOM_REAL, stat=ierr)
  endif

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in joint_inversion_initialize!", myrank, ierr)

end subroutine joint_inversion_initialize

!=====================================================================================
! Resets the joint inversion.
!=====================================================================================
subroutine joint_inversion_reset(this)
  class(t_joint_inversion), intent(inout) :: this

  call this%matrix%reset()
  call this%matrix_B%reset()

  this%b_RHS = 0._CUSTOM_REAL
  this%d_RHS = 0._CUSTOM_REAL

end subroutine joint_inversion_reset

!=====================================================================================
! Joint inversion of two field.
!=====================================================================================
subroutine joint_inversion_solve(this, par, arr, delta_model, matrix_compression_type, myrank, nbproc)
  class(t_joint_inversion), intent(inout) :: this
  type(t_parameters_inversion), intent(in) :: par
  type(t_inversion_arrays), intent(in) :: arr(2)
  integer, intent(in) :: matrix_compression_type
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), intent(out) :: delta_model(:)

  type(t_sensitivity_matrix) :: sensit
  type(t_damping) :: damping
  type(t_damping_gradient) :: damping_gradient
  type(t_parameters_lsqr) :: par_lsqr
  integer :: i, j, k, param_shift(2), nl
  ! Adjusted problem weights for joint inversion.
  real(kind=CUSTOM_REAL) :: problem_weight_adjusted(2), cost
  logical :: solve_gravity_only
  logical :: solve_mag_only
  integer :: der_type
  real(kind=CUSTOM_REAL) :: matrix_dummy(1, 1)
  character(len=32) :: filename
  integer :: nsmaller
  
  ! TODO: Move out with writing the variance.
  integer :: ierr
  real(kind=CUSTOM_REAL), allocatable :: lsqr_var_full1(:)
  real(kind=CUSTOM_REAL), allocatable :: lsqr_var_full2(:)
  real(kind=CUSTOM_REAL), allocatable :: delta_model_full(:)
  type(t_parallel_tools) :: pt

  logical :: SOLVE_PROBLEM(2)

  ! The number of times this subroutine has been called.
  ! This is effectively the major loop iteration number.
  integer, save :: ncalls = 0

  ncalls = ncalls + 1
  
  ! Initialize variance array.
  this%matrix%lsqr_var = 0.d0;
  ! Store iteration number.
  this%matrix%tag = ncalls

  ! Adjust problem weights: solve individual (grav/mag) problems for a few iterations before solving the joint system.
  if (.not. par%single_problem_complete(1, ncalls)) then
    problem_weight_adjusted(1) = par%problem_weight(1)
    problem_weight_adjusted(2) = 0.d0

  else if (.not. par%single_problem_complete(2, ncalls)) then
    problem_weight_adjusted(1) = 0.d0
    problem_weight_adjusted(2) = par%problem_weight(2)

  else
    problem_weight_adjusted = par%problem_weight
  endif

  do i = 1, 2
    SOLVE_PROBLEM(i) = (par%problem_weight(i) /= 0.d0)
  enddo

  param_shift(1) = 0
  param_shift(2) = par%nelements

  ! ***** Data misfit and damping and damping gradient  *****

  ! Loop over joint problems.
  do i = 1, 2

    ! Skip the problem with zero problem weight (single inversion).
    if (.not. SOLVE_PROBLEM(i)) then
      cycle
    endif

    if (myrank == 0 .and. i > 1) print *, '-------------------------------------------------------'

    if (myrank == 0) print *, 'Adding joint problem #', i, ' weight =', problem_weight_adjusted(i)

    call sensit%initialize(par%ndata(i), par%nelements, problem_weight_adjusted(i))

    call sensit%add(this%matrix, this%b_RHS, param_shift(i), matrix_dummy, arr(i)%matrix_sensit, &
                    arr(i)%column_weight, arr(i)%residuals, .false., .false., myrank)

    if (myrank == 0) print *, 'misfit term cost = ', sensit%get_cost()
    if (myrank == 0) print *, 'nel = ', this%matrix%get_number_elements()

    if (this%add_damping(i)) then
      if (myrank == 0) print *, 'adding damping with alpha =', par%alpha(i)

      call damping%initialize(par%nelements, par%alpha(i), problem_weight_adjusted(i), par%norm_power, &
                              par%compression_type, par%nx, par%ny, par%nz, par%wavelet_threshold)

      ! Note: we use model covariance now for the local damping weight, which is equivalent of having local alpha.
      call damping%add(this%matrix, this%b_RHS, arr(i)%column_weight, arr(i)%model%cov, &
                       arr(i)%model, arr(i)%model_prior, param_shift(i), myrank, nbproc)

      if (myrank == 0) print *, 'damping term cost = ', damping%get_cost()
      if (myrank == 0) print *, 'nel (with damping) = ', this%matrix%get_number_elements()
    endif

    if (this%add_damping_gradient(i)) then
      if (myrank == 0) print *, 'adding damping_gradient with beta =', par%beta(i), ' weight type =', par%damp_grad_weight_type

      call damping_gradient%initialize(par%beta(i), problem_weight_adjusted(i), par%nx, par%ny, par%nz, par%nelements, myrank)

      cost = 0.d0

      if (par%damp_grad_weight_type > 1) then
      ! Adding gradient weight based on the second model.
        if (i == 1) k = 2
        if (i == 2) k = 1

        ! Calculate here the local gradient weight.
        if (maxval(arr(k)%model%val_full) /= minval(arr(k)%model%val_full)) then
          ! Normalization.
          damping_gradient%grad_weight = (arr(k)%model%val_full - minval(arr(k)%model%val_full)) &
                                         / (maxval(arr(k)%model%val_full) - minval(arr(k)%model%val_full))

          damping_gradient%grad_weight = damping_gradient%grad_weight**2.d0
        else

          damping_gradient%grad_weight = arr(k)%model%val_full
        endif

      else
        damping_gradient%grad_weight = 1.d0
      endif

      do j = 1, 3 ! j is direction (1 = x, 2 = y, 3 = z).
        call damping_gradient%add(arr(i)%model, damping_gradient%grad_weight, arr(i)%column_weight, &
                                  this%matrix, this%b_RHS, param_shift(i), j, myrank, nbproc)
        cost = cost + damping_gradient%get_cost()

        if (myrank == 0) print *, 'damping_gradient term cost in direction j = ', j, damping_gradient%get_cost()
      enddo

      if (myrank == 0) print *, 'damping_gradient terms total cost = ', cost
      if (myrank == 0) print *, 'nel (with damping_gradient) = ', this%matrix%get_number_elements()
    endif
  enddo

  if (myrank == 0) print *, '-------------------------------------------------------'

  ! ***** ADMM method *****

  ! TODO: Move to a separate function like joint_inversion_add_cross_grad_constraints()

  if (par%admm_type > 0) then

    solve_gravity_only = .false.
    solve_mag_only = .false.

    if (problem_weight_adjusted(1) /= 0.d0 .and. problem_weight_adjusted(2) == 0.d0) then
      solve_gravity_only = .true.

    else if (problem_weight_adjusted(1) == 0.d0 .and. problem_weight_adjusted(2) /= 0.d0) then
      solve_mag_only = .true.
    endif

    if (solve_gravity_only) then
      i = 1
      call this%admm_method_1%iterate_admm_arrays(arr(i)%model%nlithos, &
                                                  arr(i)%model%min_local_bound, arr(i)%model%max_local_bound, &
                                                  arr(i)%model%val, this%x0_ADMM, myrank)
    else if (solve_mag_only) then
      i = 2
      call this%admm_method_2%iterate_admm_arrays(arr(i)%model%nlithos, &
                                                  arr(i)%model%min_local_bound, arr(i)%model%max_local_bound, &
                                                  arr(i)%model%val, this%x0_ADMM, myrank)
    endif

    if (solve_gravity_only .or. solve_mag_only) then
    ! Add ADMM constraints only in separate (single) inversions.

      this%weight_ADMM = 1.d0
      !this%weight_ADMM = arr(i)%damping_weight

      call damping%initialize(par%nelements, par%rho_ADMM(i), problem_weight_adjusted(i), par%norm_power, &
                              par%compression_type, par%nx, par%ny, par%nz, par%wavelet_threshold)

      call damping%add(this%matrix, this%b_RHS, arr(i)%column_weight, this%weight_ADMM, &
                       arr(i)%model, this%x0_ADMM, param_shift(i), myrank, nbproc)

      ! TODO: Calculate the norms with MPI for the parallel case.
      if (norm2(arr(i)%model%val) /= 0.d0) then
        if (myrank == 0) print *, 'ADMM |x - x0| / |x| =', norm2(arr(i)%model%val - this%x0_ADMM) / norm2(arr(i)%model%val)
      endif
      if (myrank == 0) print *, 'ADMM term cost = ', damping%get_cost()
      if (myrank == 0) print *, 'nel (with ADMM) = ', this%matrix%get_number_elements()

    else
      call this%matrix%add_empty_rows(this%nelements_total, myrank)
    endif

  endif

  ! ***** Cross-gradient *****

  if (this%add_cross_grad) then
    if (par%derivative_type /= 4) then
      if (par%derivative_type == 3) then
      ! Mixed derivative - switch derivative type every iteration.
        if (mod(ncalls, 2) == 0) then
          der_type = 1
        else
          der_type = 2
        endif
      else
        der_type = par%derivative_type
      endif

      call this%add_cross_grad_constraints(par, arr, der_type, ncalls, myrank, nbproc)
    else
    ! Adding two cross-grad terms with different derivatives.
      call this%add_cross_grad_constraints(par, arr, 1, ncalls, myrank, nbproc)
      call this%matrix_B%reset()
      call this%add_cross_grad_constraints(par, arr, 2, ncalls, myrank, nbproc)
    endif
  endif

  ! ***** Clustering *****

  if (this%add_clustering) then
    call this%add_clustering_constraints(par, arr, ncalls, myrank, nbproc)
  endif

  !-------------------------------------------------------------------------------------
  call this%matrix%finalize(2 * par%nelements, myrank)

  if (NORMALIZE_MATRIX_COLUMNS) then
    if (myrank == 0) print *, 'Normalizing the matrix columns.'

    call this%matrix%normalize_columns(this%column_norm)
  endif

  ! Parallel sparse inversion.
  delta_model = 0._CUSTOM_REAL
  if (par%method == 1) then
    call lsqr_solve(par%niter, par%rmin, par%gamma, this%matrix, this%b_RHS, delta_model, myrank)
  else
    call sca_solve(par%niter, par%rmin, this%matrix, this%b_RHS, delta_model, myrank, nbproc)
  endif

  ! ***** Method of weights *****

  if (par%method_of_weights_niter > 0 .and. norm2(this%d_RHS) /= 0._CUSTOM_REAL) then
  ! If cross gradient is not zero everywhere.

    par_lsqr%niter = par%niter
    par_lsqr%gamma = par%gamma
    par_lsqr%rmin = par%rmin

    ! Applying weighting method for equality-constrained LSQR.
    call apply_method_of_weights(par_lsqr, par%method_of_weights_niter, this%matrix, this%matrix_B, &
                                 delta_model, this%b_RHS, this%d_RHS, par%cross_grad_weight, 3, myrank)
  endif

  !-------------------------------------------------------------------------------------
  ! Scale models back to the unscaled variables.
  if (NORMALIZE_MATRIX_COLUMNS) then
    ! Calculate the column weight (which is the inverse column norm).
    do i = 1, 2 * par%nelements
      if (this%column_norm(i) /= 0.d0) then
        this%column_norm(i) = 1.d0 / this%column_norm(i)
      else
        this%column_norm(i) = 1.d0
      endif
    enddo
    call rescale_model(delta_model, this%column_norm, 2 * par%nelements)
  endif

  if (matrix_compression_type == 2) then
  ! Applying the Inverse Wavelet Transform.
    if (nbproc > 1) then
    ! Parallel version.
      allocate(delta_model_full(par%nelements_total), source=0._CUSTOM_REAL, stat=ierr)
      nsmaller = pt%get_nsmaller(par%nelements, myrank, nbproc)

      if (SOLVE_PROBLEM(1)) then
        call pt%get_full_array(delta_model(1:par%nelements), par%nelements, delta_model_full, .true., myrank, nbproc)
        call iHaar3D(delta_model_full, par%nx, par%ny, par%nz)

        ! Extract the local model update.
        delta_model(1:par%nelements) = delta_model_full(nsmaller + 1 : nsmaller + par%nelements)
      endif

      if (SOLVE_PROBLEM(2)) then
        call pt%get_full_array(delta_model(par%nelements + 1:), par%nelements, delta_model_full, .true., myrank, nbproc)
        call iHaar3D(delta_model_full, par%nx, par%ny, par%nz)

        ! Extract the local model update.
        delta_model(par%nelements + 1:) = delta_model_full(nsmaller + 1 : nsmaller + par%nelements)
      endif

    else
    ! Serial version.
      if (SOLVE_PROBLEM(1)) call iHaar3D(delta_model(1:par%nelements), par%nx, par%ny, par%nz)
      if (SOLVE_PROBLEM(2)) call iHaar3D(delta_model(par%nelements + 1:), par%nx, par%ny, par%nz)
    endif
  endif

  if (SOLVE_PROBLEM(1)) call rescale_model(delta_model(1:par%nelements), arr(1)%column_weight, par%nelements)
  if (SOLVE_PROBLEM(2)) call rescale_model(delta_model(par%nelements + 1:), arr(2)%column_weight, par%nelements)

  !--------------------------------------------------------------------------------
  ! Writing grav/mag prior and posterior variance.
  ! TODO: Move to a new function.
  !--------------------------------------------------------------------------------
  if (ncalls == 1 .or. ncalls == par%ninversions) then
    ! Calculate standard deviation sigma = sqrt(Var(X)).
    this%matrix%lsqr_var = sqrt(this%matrix%lsqr_var)

    ! Rescale with depth weight.
    if (SOLVE_PROBLEM(1)) call rescale_model(this%matrix%lsqr_var(1:par%nelements), arr(1)%column_weight, par%nelements)
    if (SOLVE_PROBLEM(2)) call rescale_model(this%matrix%lsqr_var(par%nelements + 1:), arr(2)%column_weight, par%nelements)

    if (myrank == 0) then
      ! Allocate array for the whole vector.
      if (SOLVE_PROBLEM(1)) allocate(lsqr_var_full1(par%nelements_total), source=0._CUSTOM_REAL, stat=ierr)
      if (SOLVE_PROBLEM(2)) allocate(lsqr_var_full2(par%nelements_total), source=0._CUSTOM_REAL, stat=ierr)
      if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in joint_inversion_solve!", myrank, ierr)
    endif

    ! Gather full (parallel) vector on the master.
    if (SOLVE_PROBLEM(1)) &
      call pt%get_full_array(this%matrix%lsqr_var(1:par%nelements), par%nelements, lsqr_var_full1, .false., myrank, nbproc)
    if (SOLVE_PROBLEM(2)) &
      call pt%get_full_array(this%matrix%lsqr_var(par%nelements + 1:), par%nelements, lsqr_var_full2, .false., myrank, nbproc)

    if (myrank == 0) then
    ! Writing variance by master CPU.

      if (ncalls == 1) then
        filename = "/lsqr_std_prior"
      else
        filename = "/lsqr_std_posterior"
      endif

      if (SOLVE_PROBLEM(1)) then
      ! Writing grav variance to file.
        open(27, file=trim(trim(path_output)//trim(filename)//"_grav.txt"), access='stream', form='formatted', &
             status='unknown', action='write')

        do i = 1, par%nelements_total
          write (27, *) lsqr_var_full1(i)
        enddo

        close(27)
        deallocate(lsqr_var_full1)
      endif

      if (SOLVE_PROBLEM(2)) then
      ! Writing mag variance to file.
        open(28, file=trim(trim(path_output)//trim(filename)//"_mag.txt"), access='stream', form='formatted', &
             status='unknown', action='write')

        do i = 1, par%nelements_total
          write (28, *) lsqr_var_full2(i)
        enddo

        close(28)
        deallocate(lsqr_var_full2)
      endif

    endif
  endif

end subroutine joint_inversion_solve

!=======================================================================================================
! Adds the cross-gradient constraints to the least square system.
!=======================================================================================================
subroutine joint_inversion_add_cross_grad_constraints(this, par, arr, der_type, ncalls, myrank, nbproc)
  class(t_joint_inversion), intent(inout) :: this
  type(t_parameters_inversion), intent(in) :: par
  type(t_inversion_arrays), intent(in) :: arr(2)
  integer, intent(in) :: der_type, ncalls
  integer, intent(in) :: myrank, nbproc

  type(t_vector) :: cost
  integer :: ibeg, iend

  if (myrank == 0) print *, 'Calculating cross gradients. der_type =', der_type

  call this%cross_grad%calculate(arr(1)%model, arr(2)%model, arr(1)%column_weight, arr(2)%column_weight, &
                                 this%matrix_B, this%d_RHS, this%add_cross_grad, der_type, myrank, nbproc)

  cost = this%cross_grad%get_cost()

  if (myrank == 0) print *, 'cross-grad cost = ', cost%x + cost%y + cost%z

  if (this%add_cross_grad) then
    if (par%single_problem_complete(1, ncalls) .and. par%single_problem_complete(2, ncalls)) then
      ! Adding the corresponding cross-gradient SLAE to the main system.
      call this%matrix_B%finalize(2 * par%nelements, myrank)

      ibeg = this%matrix%get_current_row_number() + 1

      call this%matrix%add_matrix(this%matrix_B, par%cross_grad_weight, myrank)

      iend = this%matrix%get_current_row_number()

      this%b_RHS(ibeg:iend) = par%cross_grad_weight * this%d_RHS

      if (myrank == 0) print *, 'cross-grad term cost = ', sum(this%b_RHS(ibeg:iend)**2)
      if (myrank == 0) print *, 'nel (with cross-grad) = ', this%matrix%get_number_elements()

    else
    ! Finding the initial solution without adding the cross-gradient term, for a few iterations.
      call this%matrix%add_empty_rows(this%matrix_B%get_total_row_number(), myrank)

    endif
  endif

end subroutine joint_inversion_add_cross_grad_constraints

!=======================================================================================================
! Adds the clustering constraints to the least square system.
!=======================================================================================================
subroutine joint_inversion_add_clustering_constraints(this, par, arr, ncalls, myrank, nbproc)
  class(t_joint_inversion), intent(inout) :: this
  type(t_parameters_inversion), intent(in) :: par
  type(t_inversion_arrays), intent(in) :: arr(2)
  integer, intent(in) :: ncalls
  integer, intent(in) :: myrank, nbproc

  integer :: i

  if (par%single_problem_complete(1, ncalls) .and. par%single_problem_complete(2, ncalls)) then
    do i = 1, 2
      call this%clustering%add(arr(1)%model, arr(2)%model, arr(1)%column_weight, arr(2)%column_weight, &
                               this%matrix, this%b_RHS, i, myrank, nbproc)

      if (myrank == 0) print *, 'clustering term', i, 'cost = ', this%clustering%get_cost(i)
      if (myrank == 0) print *, 'nel (with clustering) = ', this%matrix%get_number_elements()
    enddo
  else
  ! Finding the initial solution without adding the clustering term, for a few iterations.
    call this%matrix%add_empty_rows(2 * this%nelements_total, myrank)
  endif

end subroutine joint_inversion_add_clustering_constraints

!==================================================================================================
! Returns the cross-gradient vector magnitude at every model pixel.
!==================================================================================================
pure function joint_inversion_get_cross_grad(this) result(res)
  class(t_joint_inversion), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res(this%nelements_total)

  res = this%cross_grad%get_magnitude()

end function joint_inversion_get_cross_grad

!==================================================================================================
! Returns the clustering probabilities for every model pixel.
!==================================================================================================
pure function joint_inversion_get_clustering(this) result(res)
  class(t_joint_inversion), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res(this%nelements_total)

  res = this%clustering%get_probabilities()

end function joint_inversion_get_clustering

!==================================================================================================
! Returns the cross-gradient cost.
!==================================================================================================
pure function joint_inversion_get_cross_grad_cost(this) result(res)
  class(t_joint_inversion), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res
  type(t_vector) :: cost

  cost = this%cross_grad%get_cost()
  res = cost%get_norm()

end function joint_inversion_get_cross_grad_cost

!==================================================================================================
! Returns the clustering cost.
!==================================================================================================
pure function joint_inversion_get_clustering_cost(this, problem_type) result(res)
  class(t_joint_inversion), intent(in) :: this
  integer, intent(in) :: problem_type
  real(kind=CUSTOM_REAL) :: res

  res = this%clustering%get_cost(problem_type)

end function joint_inversion_get_clustering_cost

end module joint_inverse_problem
