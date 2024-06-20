
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

!==========================================================================
! A class to add cross-gradient joint inversion contribution
! to the System of Linear Algebraic Equations (SLAE) that is stored
! using Compressed Sparse Row (CSR) format.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!==========================================================================
module cross_gradient

  use global_typedefs, only: CUSTOM_REAL
  use mpi_tools, only: exit_MPI
  use sparse_matrix
  use vector
  use parallel_tools
  use string, only: str
  use gradient
  use grid

  implicit none

  private

  !--------------------------------------------------------------------
  ! Data structure to store the cross-gradients data for one element.
  type, private :: t_tau
    ! Cross-gradients vector.
    type(t_vector) :: val

    ! Partial derivatives with respect to the model parameters.
    type(t_vector) :: dm1(5)
    type(t_vector) :: dm2(5)

    ! Indexes in the model corresponding to derivatives.
    type(t_ivector) :: ind(5)
  end type t_tau

  !----------------------------------------------------------------
  ! Main cross-gradient class.
  type, public :: t_cross_gradient
    private

    ! Problem dimensions.
    integer :: nx, ny, nz

    ! Cost (that we want to minimize).
    type(t_vector) :: cost

    ! Total number of parameters.
    integer :: nparams

    ! Local number of parameters.
    integer :: nparams_loc

    ! Flags to define one of the models constant so that it remains unaltered during the inversion.
    logical :: keep_model_constant(2)

    ! Cross gradient vector magnitude (for visualization).
    real(kind=CUSTOM_REAL), allocatable :: cross_grad(:)

    ! A flag for using a vector field for structural constraints.
    ! 0 - don't use, 1 - use for the 1st model, 2 - use for the 2nd model.
    integer :: vec_field_type

    ! Vector field to use as structural constraints.
    real(kind=CUSTOM_REAL), public, allocatable :: vec_field(:, :)

    ! Vector field file.
    character(len=256) :: vec_field_file

    integer(kind=8), public :: nnz
    integer, public :: nl_nonempty

  contains
    private

    procedure, public, pass :: initialize => cross_gradient_initialize
    procedure, public, pass :: calculate => cross_gradient_calculate

    procedure, public, pass :: get_cost => cross_gradient_get_cost
    procedure, public, pass :: get_magnitude => cross_gradient_get_magnitude

    procedure, private :: get_model_gradients => cross_gradient_get_model_gradients

    procedure, private :: calculate_tau => cross_gradient_calculate_tau
    procedure, private, nopass :: normalize_tau => cross_gradient_normalize_tau

    procedure, private :: calculate_tau2 => cross_gradient_calculate_tau2
    procedure, private :: calculate_tau_backward => cross_gradient_calculate_tau_backward

    procedure, private, nopass :: get_num_deriv => cross_gradient_get_num_deriv

  end type t_cross_gradient

  private :: read_vector_field

contains

!============================================================================================================
! Initialize cross-gradients: set dimensions and grid steps.
!============================================================================================================
subroutine cross_gradient_initialize(this, nx, ny, nz, nparams_loc, keep_model_constant, myrank)
  class(t_cross_gradient), intent(inout) :: this
  integer, intent(in) :: nx, ny, nz, nparams_loc
  integer, intent(in) :: keep_model_constant(2)
  integer, intent(in) :: myrank

  integer :: ierr

  if (myrank == 0) print *, "Initializing cross-gradient constraints."

  this%nx = nx
  this%ny = ny
  this%nz = nz

  this%nparams_loc = nparams_loc
  this%nparams = nx * ny * nz

  this%keep_model_constant = (keep_model_constant > 0)

  this%cost = 0._CUSTOM_REAL
  this%nnz = 0
  this%nl_nonempty = 0

  ierr = 0

  if (myrank == 0) then
    allocate(this%cross_grad(this%nparams), stat=ierr)
    if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in cross_gradient_initialize!", myrank, ierr)
  endif

  !----------------------------------------
  ! TODO: Expose to Parfile.
  this%vec_field_type = 2
  this%vec_field_file = 'one_sphere_vector2.txt'

  if (this%vec_field_type > 0) then
    ! Keep dimensions in this order for compatibility with paraview visualisation interface.
    allocate(this%vec_field(this%nparams, 3), stat=ierr)
    if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in cross_gradient_initialize!", myrank, ierr)

    call read_vector_field(this%nparams, this%vec_field, this%vec_field_file, myrank)
  endif

end subroutine cross_gradient_initialize

!==========================================================================================================
! Read the vector field.
!==========================================================================================================
subroutine read_vector_field(nparams, vec_field, file_name, myrank)
  integer, intent(in) :: nparams
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank
  real(kind=CUSTOM_REAL), intent(out) :: vec_field(nparams, 3)

  integer :: ierr, i, nelements_read
  character(len=256) :: msg

  if (myrank == 0) print *, 'Reading the vector field from file ', trim(file_name)

  open(10, file=trim(file_name), status='old', action='read', iostat=ierr, iomsg=msg)
  if (ierr /= 0) call exit_MPI("Error in opening the vector field file! path=" &
                 //file_name//" iomsg="//msg, myrank, ierr)

  read(10, *, iostat=ierr) nelements_read
  if (ierr /= 0) call exit_MPI("Problem while reading the vector field file header!", myrank, ierr)

  ! Sanity check.
  if (nparams /= nelements_read) &
    call exit_MPI("The vector field is not correctly defined!"//new_line('a') &
          //"nelements_read="//str(nelements_read)//new_line('a') &
          //"nelements_total="//str(nparams), myrank, 0)

  ! Reading the vectors.
  do i = 1, nparams
      read(10, *, iostat=ierr) vec_field(i, 1), vec_field(i, 2), vec_field(i, 3)

      if (ierr /= 0) &
        call exit_MPI("Problem with reading the vector field for pixel i = "//str(i), myrank, ierr)
  enddo

  close(10)

end subroutine read_vector_field

!=============================================================================================
! Returns the number of derivatives depending on the finite difference scheme.
!=============================================================================================
pure function cross_gradient_get_num_deriv(der_type) result(num_deriv)
  integer, intent(in) :: der_type
  integer :: num_deriv

  if (der_type == 0 .or. der_type == 1) then
    num_deriv = 3
  else if (der_type == 6) then
    num_deriv = 5
  else
    num_deriv = 4
  endif

end function cross_gradient_get_num_deriv

!=====================================================================================================
! Calculates the cross-gradients of the model1 and model2 and stores in cross_grad(:).
! If the flag add = 'true', then adds cross-gradients to the sparse matrix and right-hand side b_RHS.
!=====================================================================================================
subroutine cross_gradient_calculate(this, model1, model2, grid, column_weight1, column_weight2, &
                                    matrix, b_RHS, add, der_type, glob_weight, myrank, nbproc)
  class(t_cross_gradient), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: model1(:)
  real(kind=CUSTOM_REAL), intent(in) :: model2(:)
  type(t_grad_grid), intent(in) :: grid
  real(kind=CUSTOM_REAL), intent(in) :: column_weight1(:)
  real(kind=CUSTOM_REAL), intent(in) :: column_weight2(:)
  real(kind=CUSTOM_REAL), intent(in) :: glob_weight
  logical, intent(in) :: add
  integer, intent(in) :: der_type
  integer, intent(in) :: myrank, nbproc

  type(t_sparse_matrix), intent(inout) :: matrix
  real(kind=CUSTOM_REAL), intent(inout) :: b_RHS(:)

  ! Cross-gradient data.
  type(t_tau) :: tau
  integer :: i, j, k, l, p
  integer :: ind, nderiv
  integer :: nsmaller
  real(kind=CUSTOM_REAL) :: val1, val2
  logical :: on_left_boundary, on_right_boundary
  logical :: index_included

  ! Set the number of derivatives.
  nderiv = this%get_num_deriv(der_type)

  this%cost = 0._CUSTOM_REAL

  ! Number of parameters on ranks smaller than current one.
  nsmaller = get_nsmaller(this%nparams_loc, myrank, nbproc)

  this%nl_nonempty = 0
  this%nnz = 0
  p = 0
  do k = 1, this%nz
    do j = 1, this%ny
      do i = 1, this%nx
        p = p + 1

        ! Calculate cross-gradient discretization.
        if (der_type == 1 .or. der_type == 2) then
          ! Apply boundary conditions.
          on_left_boundary = .false.
          on_right_boundary = .false.

          if (i == 1 .or. j == 1 .or. k == 1) on_left_boundary = .true.
          if (i == this%nx .or. j == this%ny .or. k == this%nz) on_right_boundary = .true.

          if (on_left_boundary .and. on_right_boundary) then
          ! Skip the cross gradient constraints.
            tau%val = 0.d0
            tau%dm1 = t_vector(0.d0, 0.d0, 0.d0)
            tau%dm2 = t_vector(0.d0, 0.d0, 0.d0)

          else if (der_type == 2 .and. on_left_boundary) then
          ! On the left boundary: use forward difference.
            tau = this%calculate_tau(model1, model2, grid, i, j, k, 1)

          else if (on_right_boundary) then
          ! On the right boundary: use backward difference.
            tau = this%calculate_tau_backward(model1, model2, grid, i, j, k)

          else
          ! Inside the domain.
            tau = this%calculate_tau(model1, model2, grid, i, j, k, der_type)
          endif

        else
          call exit_MPI("Unsupported derivative type!", myrank, der_type)
        endif

        ! Set derivatives to zero to keep the second model constant.
        if (this%keep_model_constant(1)) tau%dm1 = t_vector(0.d0, 0.d0, 0.d0)
        if (this%keep_model_constant(2)) tau%dm2 = t_vector(0.d0, 0.d0, 0.d0)

        ! Normalization.
        !call this%normalize_tau(tau, model1, model2, i, j, k, myrank)

        if (myrank == 0) then
          ! Store cross gradient vector magnitude.
          this%cross_grad(p) = tau%val%get_norm()
        endif

        ! Calculate the cost.
        this%cost%x = this%cost%x + tau%val%x**2
        this%cost%y = this%cost%y + tau%val%y**2
        this%cost%z = this%cost%z + tau%val%z**2

        ! Adding constraints to the matrix and right-hand-side.

        ! Row with x-component.
        index_included = .false.
        do l = 1, nderiv
          ind = tau%ind(l)%x

          if (ind > nsmaller .and. ind <= nsmaller + this%nparams_loc) then
            if (add) then
              ind = ind - nsmaller
              val1 = tau%dm1(l)%x * column_weight1(ind) * glob_weight
              val2 = tau%dm2(l)%x * column_weight2(ind) * glob_weight
              call matrix%add(val1, ind, myrank)
              call matrix%add(val2, ind + this%nparams_loc, myrank)
            endif
            this%nnz = this%nnz + 2
            index_included = .true.
          endif
        enddo
        if (add) then
          call matrix%new_row(myrank)
          b_RHS(matrix%get_current_row_number()) = -tau%val%x * glob_weight
        endif
        if (index_included) this%nl_nonempty = this%nl_nonempty + 1

        ! Row with y-component.
        index_included = .false.
        do l = 1, nderiv
          ind = tau%ind(l)%y

          if (ind > nsmaller .and. ind <= nsmaller + this%nparams_loc) then
            if (add) then
              ind = ind - nsmaller
              val1 = tau%dm1(l)%y * column_weight1(ind) * glob_weight
              val2 = tau%dm2(l)%y * column_weight2(ind) * glob_weight
              call matrix%add(val1, ind, myrank)
              call matrix%add(val2, ind + this%nparams_loc, myrank)
            endif
            this%nnz = this%nnz + 2
            index_included = .true.
          endif
        enddo
        if (add) then
          call matrix%new_row(myrank)
          b_RHS(matrix%get_current_row_number()) = -tau%val%y * glob_weight
        endif
        if (index_included) this%nl_nonempty = this%nl_nonempty + 1

        ! Row with z-component.
        index_included = .false.
        do l = 1, nderiv
          ind = tau%ind(l)%z

          if (ind > nsmaller .and. ind <= nsmaller + this%nparams_loc) then
            if (add) then
              ind = ind - nsmaller
              val1 = tau%dm1(l)%z * column_weight1(ind) * glob_weight
              val2 = tau%dm2(l)%z * column_weight2(ind) * glob_weight
              call matrix%add(val1, ind, myrank)
              call matrix%add(val2, ind + this%nparams_loc, myrank)
            endif
            this%nnz = this%nnz + 2
            index_included = .true.
          endif
        enddo
        if (add) then
          call matrix%new_row(myrank)
          b_RHS(matrix%get_current_row_number()) = -tau%val%z * glob_weight
        endif
        if (index_included) this%nl_nonempty = this%nl_nonempty + 1
      enddo
    enddo
  enddo

  ! Adjust the nnz taking into account the 'keep_model_constant' flag.
  if (this%keep_model_constant(1) .and. this%keep_model_constant(2)) then
    this%nnz = 0
  else if (this%keep_model_constant(1) .or. this%keep_model_constant(2)) then
    this%nnz = this%nnz / 2
  endif

end subroutine cross_gradient_calculate

!==================================================================================================
! Returns the cost for every component (this is what we want to minimize).
!==================================================================================================
pure function cross_gradient_get_cost(this) result(res)
  class(t_cross_gradient), intent(in) :: this
  type(t_vector) :: res

  res = this%cost

end function cross_gradient_get_cost

!==================================================================================================
! Returns the cross-gradient vector at every model pixel.
!==================================================================================================
pure subroutine cross_gradient_get_magnitude(this, cross_grad)
  class(t_cross_gradient), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(out) :: cross_grad(this%nparams)

  cross_grad = this%cross_grad

end subroutine cross_gradient_get_magnitude

!==================================================================================================
! Retrieve the model gradients.
! We either calculate the gradient from the model, or take it from the input vector field.
!==================================================================================================
subroutine cross_gradient_get_model_gradients(this, model1, model2, grid, i, j, k, der_type, m1_grad, m2_grad)
  class(t_cross_gradient), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: model1(:)
  real(kind=CUSTOM_REAL), intent(in) :: model2(:)
  type(t_grad_grid), intent(in) :: grid
  integer, intent(in) :: i, j, k
  character(len=4) :: der_type
  type(t_vector), intent(out) :: m1_grad, m2_grad

  type(t_vector) :: vec
  integer :: ind

  if (this%vec_field_type > 0) then
    ! Use gradient from the provided vector field.
    ind = grid%get_ind(i, j, k)
    vec = t_vector(this%vec_field(ind, 1), this%vec_field(ind, 2), this%vec_field(ind, 3))
  endif

  ! Calculate model gradients.
  if (this%vec_field_type == 1) then
    m1_grad = vec
  else
    m1_grad = get_grad(model1, grid, i, j, k, der_type)
  endif
  if (this%vec_field_type == 2) then
    m2_grad = vec
  else
    m2_grad = get_grad(model2, grid, i, j, k, der_type)
  endif

end subroutine cross_gradient_get_model_gradients

!==================================================================================================
! Calculates the discretized cross-gradients function (and corresponding matrix column indexes)
! between the models model1 and model2, at pixel location (i, j, k).
! der_type = 1:  using a forward difference scheme (Geophys. Res. Lett., Vol. 33, L07303).
! der_type /= 1: using a central difference scheme.
!==================================================================================================
function cross_gradient_calculate_tau(this, model1, model2, grid, i, j, k, der_type) result(tau)
  class(t_cross_gradient), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: model1(:)
  real(kind=CUSTOM_REAL), intent(in) :: model2(:)
  type(t_grad_grid), intent(in) :: grid
  integer, intent(in) :: i, j, k, der_type

  type(t_tau) :: tau
  type(t_vector) :: m1_grad, m2_grad
  type(t_vector) :: step

  ! Retrieve the gradients.
  call this%get_model_gradients(model1, model2, grid, i, j, k, get_der_type(der_type), m1_grad, m2_grad)

  ! Calculate the cross-product between model gradients.
  tau%val = m1_grad%cross_product(m2_grad)

  ! Calculate partial derivatives (with respect to the model parameters),
  ! and corresponding matrix column indexes.

  step = t_vector(grid%dX(i), grid%dY(j), grid%dZ(k))

  if (der_type /= 1) then
    step = 2._CUSTOM_REAL * step
  endif

  !------------------------------------------------------------------
  ! Note: assume both models have the same grid indexes.

  tau%dm1(1)%x =   m2_grad%z / step%y   ! dx / dm1(i, j + 1, k)
  tau%dm2(1)%x = - m1_grad%z / step%y   ! dx / dm2(i, j + 1, k)

  tau%dm1(2)%x = - m2_grad%y / step%z   ! dx / dm1(i, j, k + 1)
  tau%dm2(2)%x =   m1_grad%y / step%z   ! dx / dm2(i, j, k + 1)

  tau%ind(1)%x = grid%get_ind(i, j + 1, k)
  tau%ind(2)%x = grid%get_ind(i, j, k + 1)

  if (der_type == 1) then
    tau%dm1(3)%x = - (m2_grad%z / step%y - m2_grad%y / step%z)  ! dx / dm1(i, j, k)
    tau%dm2(3)%x = - (m1_grad%y / step%z - m1_grad%z / step%y)  ! dx / dm2(i, j, k)

    tau%ind(3)%x = grid%get_ind(i, j, k)

  else
    tau%dm1(3)%x = - tau%dm1(1)%x   ! dx / dm1(i, j - 1, k)
    tau%dm2(3)%x = - tau%dm2(1)%x   ! dx / dm2(i, j - 1, k)

    tau%dm1(4)%x = - tau%dm1(2)%x   ! dx / dm1(i, j, k - 1)
    tau%dm2(4)%x = - tau%dm2(2)%x   ! dx / dm2(i, j, k - 1)

    tau%ind(3)%x = grid%get_ind(i, j - 1, k)
    tau%ind(4)%x = grid%get_ind(i, j, k - 1)
  endif

  !------------
  tau%dm1(1)%y = - m2_grad%z / step%x   ! dy / dm1(i + 1, j, k)
  tau%dm2(1)%y =   m1_grad%z / step%x   ! dy / dm2(i + 1, j, k)

  tau%dm1(2)%y =   m2_grad%x / step%z   ! dy / dm1(i, j, k + 1)
  tau%dm2(2)%y = - m1_grad%x / step%z   ! dy / dm2(i, j, k + 1)

  tau%ind(1)%y = grid%get_ind(i + 1, j, k)
  tau%ind(2)%y = grid%get_ind(i, j, k + 1)

  if (der_type == 1) then
    tau%dm1(3)%y = - (m2_grad%x / step%z - m2_grad%z / step%x)  ! dy / dm1(i, j, k)
    tau%dm2(3)%y = - (m1_grad%z / step%x - m1_grad%x / step%z)  ! dy / dm2(i, j, k)

    tau%ind(3)%y = grid%get_ind(i, j, k)

  else
    tau%dm1(3)%y = - tau%dm1(1)%y   ! dy / dm1(i - 1, j, k)
    tau%dm2(3)%y = - tau%dm2(1)%y   ! dy / dm2(i - 1, j, k)

    tau%dm1(4)%y = - tau%dm1(2)%y   ! dy / dm1(i, j, k - 1)
    tau%dm2(4)%y = - tau%dm2(2)%y   ! dy / dm2(i, j, k - 1)

    tau%ind(3)%y = grid%get_ind(i - 1, j, k)
    tau%ind(4)%y = grid%get_ind(i, j, k - 1)
  endif

  !------------
  tau%dm1(1)%z =   m2_grad%y / step%x   ! dz / dm1(i + 1, j, k)
  tau%dm2(1)%z = - m1_grad%y / step%x   ! dz / dm2(i + 1, j, k)

  tau%dm1(2)%z = - m2_grad%x / step%y   ! dz / dm1(i, j + 1, k)
  tau%dm2(2)%z =   m1_grad%x / step%y   ! dz / dm2(i, j + 1, k)

  tau%ind(1)%z = grid%get_ind(i + 1, j, k)
  tau%ind(2)%z = grid%get_ind(i, j + 1, k)

  if (der_type == 1) then
    tau%dm1(3)%z = - (m2_grad%y / step%x - m2_grad%x / step%y)  ! dz / dm1(i, j, k)
    tau%dm2(3)%z = - (m1_grad%x / step%y - m1_grad%y / step%x)  ! dz / dm2(i, j, k)

    tau%ind(3)%z = grid%get_ind(i, j, k)

  else
    tau%dm1(3)%z = - tau%dm1(1)%z   ! dz / dm1(i - 1, j, k)
    tau%dm2(3)%z = - tau%dm2(1)%z   ! dz / dm2(i - 1, j, k)

    tau%dm1(4)%z = - tau%dm1(2)%z   ! dz / dm1(i, j - 1, k)
    tau%dm2(4)%z = - tau%dm2(2)%z   ! dz / dm2(i, j - 1, k)

    tau%ind(3)%z = grid%get_ind(i - 1, j, k)
    tau%ind(4)%z = grid%get_ind(i, j - 1, k)
  endif

  ! To avoid having non-initialized values.
  if (der_type == 1) then
    tau%ind(4:) = t_ivector(0, 0, 0)
    tau%dm1(4:) = t_vector(0.d0, 0.d0, 0.d0)
    tau%dm2(4:) = t_vector(0.d0, 0.d0, 0.d0)
  else
    tau%ind(5:) = t_ivector(0, 0, 0)
    tau%dm1(5:) = t_vector(0.d0, 0.d0, 0.d0)
    tau%dm2(5:) = t_vector(0.d0, 0.d0, 0.d0)
  endif

end function cross_gradient_calculate_tau

!==================================================================================================
! Calculates the discretized cross-gradients function (and corresponding matrix column indexes)
! between the models model1 and model2, at pixel location (i, j, k).
!
! Same as cross_gradient_calculate_tau but uses three-point forward difference scheme.
!==================================================================================================
function cross_gradient_calculate_tau2(this, model1, model2, grid, i, j, k) result(tau)
  class(t_cross_gradient), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: model1(:)
  real(kind=CUSTOM_REAL), intent(in) :: model2(:)
  type(t_grad_grid), intent(in) :: grid
  integer, intent(in) :: i, j, k

  type(t_tau) :: tau
  type(t_vector) :: m1_grad, m2_grad
  type(t_vector) :: step

  ! Retrieve the gradients.
  call this%get_model_gradients(model1, model2, grid, i, j, k, FWD2_TYPE, m1_grad, m2_grad)

  ! Calculate the cross-product between model gradients.
  tau%val = m1_grad%cross_product(m2_grad)

  ! Calculate partial derivatives (with respect to the model parameters),
  ! and corresponding matrix column indexes.

  step = 2.d0 * t_vector(grid%dX(i), grid%dY(j), grid%dZ(k))

  !------------------------------------------------------------------
  ! Note: assume both models have the same grid indexes.

  ! X-component.
  tau%dm1(1)%x = - 1.d0 * m2_grad%z / step%y                         ! dx / dm1(i, j + 2, k)
  tau%dm1(2)%x = + 4.d0 * m2_grad%z / step%y                         ! dx / dm1(i, j + 1, k)
  tau%dm1(3)%x = + 1.d0 * m2_grad%y / step%z                         ! dx / dm1(i, j, k + 2)
  tau%dm1(4)%x = - 4.d0 * m2_grad%y / step%z                         ! dx / dm1(i, j, k + 1)
  tau%dm1(5)%x = - 3.d0 * (m2_grad%z / step%y - m2_grad%y / step%z)  ! dx / dm1(i, j, k)

  tau%dm2(1)%x = + 1.d0 * m1_grad%z / step%y                         ! dx / dm2(i, j + 2, k)
  tau%dm2(2)%x = - 4.d0 * m1_grad%z / step%y                         ! dx / dm2(i, j + 1, k)
  tau%dm2(3)%x = - 1.d0 * m1_grad%y / step%z                         ! dx / dm2(i, j, k + 2)
  tau%dm2(4)%x = + 4.d0 * m1_grad%y / step%z                         ! dx / dm2(i, j, k + 1)
  tau%dm2(5)%x = - 3.d0 * (m1_grad%y / step%z - m1_grad%z / step%y)  ! dx / dm2(i, j, k)

  tau%ind(1)%x = grid%get_ind(i, j + 2, k)
  tau%ind(2)%x = grid%get_ind(i, j + 1, k)
  tau%ind(3)%x = grid%get_ind(i, j, k + 2)
  tau%ind(4)%x = grid%get_ind(i, j, k + 1)
  tau%ind(5)%x = grid%get_ind(i, j, k)

  ! Y-component.
  tau%dm1(1)%y = - 1.d0 * m2_grad%x / step%z                         ! dy / dm1(i, j, k + 2)
  tau%dm1(2)%y = + 4.d0 * m2_grad%x / step%z                         ! dy / dm1(i, j, k + 1)
  tau%dm1(3)%y = + 1.d0 * m2_grad%z / step%x                         ! dy / dm1(i + 2, j, k)
  tau%dm1(4)%y = - 4.d0 * m2_grad%z / step%x                         ! dy / dm1(i + 1, j, k)
  tau%dm1(5)%y = - 3.d0 * (m2_grad%x / step%z - m2_grad%z / step%x)  ! dy / dm1(i, j, k)

  tau%dm2(1)%y = + 1.d0 * m1_grad%x / step%z                         ! dy / dm2(i, j, k + 2)
  tau%dm2(2)%y = - 4.d0 * m1_grad%x / step%z                         ! dy / dm2(i, j, k + 1)
  tau%dm2(3)%y = - 1.d0 * m1_grad%z / step%x                         ! dy / dm2(i + 2, j, k)
  tau%dm2(4)%y = + 4.d0 * m1_grad%z / step%x                         ! dy / dm2(i + 1, j, k)
  tau%dm2(5)%y = - 3.d0 * (m1_grad%z / step%x - m1_grad%x / step%z)  ! dy / dm2(i, j, k)

  tau%ind(1)%y = grid%get_ind(i, j, k + 2)
  tau%ind(2)%y = grid%get_ind(i, j, k + 1)
  tau%ind(3)%y = grid%get_ind(i + 2, j, k)
  tau%ind(4)%y = grid%get_ind(i + 1, j, k)
  tau%ind(5)%y = grid%get_ind(i, j, k)

  ! Z-component.
  tau%dm1(1)%z = - 1.d0 * m2_grad%y / step%x                         ! dz / dm1(i + 2, j, k)
  tau%dm1(2)%z = + 4.d0 * m2_grad%y / step%x                         ! dz / dm1(i + 1, j, k)
  tau%dm1(3)%z = + 1.d0 * m2_grad%x / step%y                         ! dz / dm1(i, j + 2, k)
  tau%dm1(4)%z = - 4.d0 * m2_grad%x / step%y                         ! dz / dm1(i, j + 1, k)
  tau%dm1(5)%z = - 3.d0 * (m2_grad%y / step%x - m2_grad%x / step%y)  ! dz / dm1(i, j, k)

  tau%dm2(1)%z = + 1.d0 * m1_grad%y / step%x                         ! dz / dm2(i + 2, j, k)
  tau%dm2(2)%z = - 4.d0 * m1_grad%y / step%x                         ! dz / dm2(i + 1, j, k)
  tau%dm2(3)%z = - 1.d0 * m1_grad%x / step%y                         ! dz / dm2(i, j + 2, k)
  tau%dm2(4)%z = + 4.d0 * m1_grad%x / step%y                         ! dz / dm2(i, j + 1, k)
  tau%dm2(5)%z = - 3.d0 * (m1_grad%x / step%y - m1_grad%y / step%x)  ! dz / dm2(i, j, k)

  tau%ind(1)%z = grid%get_ind(i + 2, j, k)
  tau%ind(2)%z = grid%get_ind(i + 1, j, k)
  tau%ind(3)%z = grid%get_ind(i, j + 2, k)
  tau%ind(4)%z = grid%get_ind(i, j + 1, k)
  tau%ind(5)%z = grid%get_ind(i, j, k)

end function cross_gradient_calculate_tau2

!==================================================================================================
! Calculates the discretized cross-gradients function (and corresponding matrix column indexes)
! between the models model1 and model2, at pixel location (i, j, k).
!
! Same as cross_gradient_calculate_tau but uses backward finite difference.
!==================================================================================================
function cross_gradient_calculate_tau_backward(this, model1, model2, grid, i, j, k) result(tau)
  class(t_cross_gradient), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: model1(:)
  real(kind=CUSTOM_REAL), intent(in) :: model2(:)
  type(t_grad_grid), intent(in) :: grid
  integer, intent(in) :: i, j, k

  type(t_tau) :: tau
  type(t_vector) :: m1_grad, m2_grad
  type(t_vector) :: step

  ! Retrieve the gradients.
  call this%get_model_gradients(model1, model2, grid, i, j, k, BWD_TYPE, m1_grad, m2_grad)

  ! Calculate the cross-product between model gradients.
  tau%val = m1_grad%cross_product(m2_grad)

  ! Calculate partial derivatives (with respect to the model parameters),
  ! and corresponding matrix column indexes.

  step = t_vector(grid%dX(i), grid%dY(j), grid%dZ(k))

  !------------------------------------------------------------------
  ! Note: assume both models have the same grid indexes.

  tau%dm1(1)%x = - m2_grad%z / step%y                         ! dx / dm1(i, j - 1, k)
  tau%dm1(2)%x =   m2_grad%y / step%z                         ! dx / dm1(i, j, k - 1)
  tau%dm1(3)%x =   m2_grad%z / step%y - m2_grad%y / step%z    ! dx / dm1(i, j, k)

  tau%dm2(1)%x =   m1_grad%z / step%y                         ! dx / dm2(i, j - 1, k)
  tau%dm2(2)%x = - m1_grad%y / step%z                         ! dx / dm2(i, j, k - 1)
  tau%dm2(3)%x =   m1_grad%y / step%z - m1_grad%z / step%y    ! dx / dm2(i, j, k)

  tau%ind(1)%x = grid%get_ind(i, j - 1, k)
  tau%ind(2)%x = grid%get_ind(i, j, k - 1)
  tau%ind(3)%x = grid%get_ind(i, j, k)

  !------------
  tau%dm1(1)%y =   m2_grad%z / step%x                         ! dy / dm1(i - 1, j, k)
  tau%dm1(2)%y = - m2_grad%x / step%z                         ! dy / dm1(i, j, k - 1)
  tau%dm1(3)%y =   m2_grad%x / step%z - m2_grad%z / step%x    ! dy / dm1(i, j, k)

  tau%dm2(1)%y = - m1_grad%z / step%x                         ! dy / dm2(i - 1, j, k)
  tau%dm2(2)%y =   m1_grad%x / step%z                         ! dy / dm2(i, j, k - 1)
  tau%dm2(3)%y =   m1_grad%z / step%x - m1_grad%x / step%z    ! dy / dm2(i, j, k)

  tau%ind(1)%y = grid%get_ind(i - 1, j, k)
  tau%ind(2)%y = grid%get_ind(i, j, k - 1)
  tau%ind(3)%y = grid%get_ind(i, j, k)

  !------------
  tau%dm1(1)%z = - m2_grad%y / step%x                         ! dz / dm1(i - 1, j, k)
  tau%dm1(2)%z =   m2_grad%x / step%y                         ! dz / dm1(i, j - 1, k)
  tau%dm1(3)%z =   m2_grad%y / step%x - m2_grad%x / step%y    ! dz / dm1(i, j, k)

  tau%dm2(1)%z =   m1_grad%y / step%x                         ! dz / dm2(i - 1, j, k)
  tau%dm2(2)%z = - m1_grad%x / step%y                         ! dz / dm2(i, j - 1, k)
  tau%dm2(3)%z =   m1_grad%x / step%y - m1_grad%y / step%x    ! dz / dm2(i, j, k)

  tau%ind(1)%z = grid%get_ind(i - 1, j, k)
  tau%ind(2)%z = grid%get_ind(i, j - 1, k)
  tau%ind(3)%z = grid%get_ind(i, j, k)

  ! To avoid having non-initialized values.
  tau%dm1(4:) = t_vector(0.d0, 0.d0, 0.d0)
  tau%dm2(4:) = t_vector(0.d0, 0.d0, 0.d0)
  tau%ind(4:) = t_ivector(0, 0, 0)

end function cross_gradient_calculate_tau_backward

!==============================================================================================
! Normalize the cross gradient function.
! (See N Linde, J Doetsch, Joint inversion of crosshole GPR and Seismic traveltime data.)
!==============================================================================================
subroutine cross_gradient_normalize_tau(tau, model1, model2, grid, i, j, k)
  real(kind=CUSTOM_REAL), intent(in) :: model1(:)
  real(kind=CUSTOM_REAL), intent(in) :: model2(:)
  type(t_grad_grid), intent(in) :: grid
  integer, intent(in) :: i, j, k

  type(t_tau), intent(inout) :: tau
  real(kind=CUSTOM_REAL) :: val1, val2, scale
  integer :: l

  val1 = grad_get_par(model1, grid, i, j, k)
  val2 = grad_get_par(model2, grid, i, j, k)

  if (val1 /= 0._CUSTOM_REAL .and. val2 /= 0._CUSTOM_REAL) then
    scale = 1._CUSTOM_REAL / abs(val1 * val2)
  else
  ! Zero denominator!
    scale = 1._CUSTOM_REAL
  endif

  tau%val = scale * tau%val

  do l = 1, 5
    tau%dm1(l) = scale * tau%dm1(l)
    tau%dm2(l) = scale * tau%dm2(l)
  enddo

end subroutine cross_gradient_normalize_tau

end module cross_gradient
