
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
  use model
  use parallel_tools
  use string, only: str
  use gradient

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

    ! Use these when compute different types of derivatives (e.g., fwd and central) for different models.
    type(t_ivector) :: ind1(5)
    type(t_ivector) :: ind2(5)
  end type t_tau

  !----------------------------------------------------------------
  ! Main cross-gradient class.
  type, public :: t_cross_gradient
    private

    ! Problem dimensions.
    type(t_ivector) :: size

    ! Cost (that we want to minimize).
    type(t_vector) :: cost

    ! Total number of parameters.
    integer :: nparams

    ! Local number of parameters.
    integer :: nparams_loc

    ! Cross gradient vector magnitude (for visualization).
    real(kind=CUSTOM_REAL), allocatable :: cross_grad(:)

  contains
    private

    procedure, public, pass :: initialize => cross_gradient_initialize
    procedure, public, pass :: calculate => cross_gradient_calculate

    procedure, public, pass :: get_cost => cross_gradient_get_cost
    procedure, public, pass :: get_magnitude => cross_gradient_get_magnitude
    procedure, public, pass :: get_num_elements => cross_gradient_get_num_elements

    procedure, private, nopass :: calculate_tau => cross_gradient_calculate_tau
    procedure, private, nopass :: normalize_tau => cross_gradient_normalize_tau

    procedure, private, nopass :: calculate_tau_mixed_gradients => cross_gradient_calculate_tau_mixed_gradients
    procedure, private, nopass :: calculate_tau2 => cross_gradient_calculate_tau2
    procedure, private, nopass :: calculate_tau_backward => cross_gradient_calculate_tau_backward

    procedure, private, nopass :: get_num_deriv => cross_gradient_get_num_deriv

  end type t_cross_gradient

contains

!============================================================================================================
! Initialize cross-gradients: set dimensions and grid steps.
!============================================================================================================
subroutine cross_gradient_initialize(this, nx, ny, nz, nparams_loc, myrank)
  class(t_cross_gradient), intent(inout) :: this
  integer, intent(in) :: nx, ny, nz, nparams_loc
  integer, intent(in) :: myrank

  integer :: ierr

  if (myrank == 0) print *, "Initializing cross-gradient constraints."

  this%size%x = nx
  this%size%y = ny
  this%size%z = nz

  this%nparams_loc = nparams_loc
  this%nparams = nx * ny * nz

  this%cost = 0._CUSTOM_REAL

  ierr = 0

  if (.not. allocated(this%cross_grad)) allocate(this%cross_grad(this%nparams), stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in cross_gradient_initialize!", myrank, ierr)

end subroutine cross_gradient_initialize

!=============================================================================================
! Returns the estimated (local) number of elements to be added to the matrix.
!=============================================================================================
pure function cross_gradient_get_num_elements(this, der_type) result(res)
  class(t_cross_gradient), intent(in) :: this
  integer, intent(in) :: der_type
  integer :: res

  res = 3 * this%get_num_deriv(der_type) * 2 * this%nparams_loc

  ! Add this term as derivatives of elements on the boundaries may belong to other CPUs, which leads to more elements.
  ! Note: this can be optimized by adding a smaller number. But I dunno which, e.g., nparams_loc is not enough (tested).
  res = res + this%nparams

end function cross_gradient_get_num_elements

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
subroutine cross_gradient_calculate(this, model1, model2, column_weight1, column_weight2, &
                                    matrix, b_RHS, add, der_type, myrank, nbproc)
  class(t_cross_gradient), intent(inout) :: this
  type(t_model), intent(inout) :: model1
  type(t_model), intent(inout) :: model2
  real(kind=CUSTOM_REAL), intent(in) :: column_weight1(:)
  real(kind=CUSTOM_REAL), intent(in) :: column_weight2(:)
  logical, intent(in) :: add
  integer, intent(in) :: der_type
  integer, intent(in) :: myrank, nbproc

  type(t_sparse_matrix), intent(inout) :: matrix
  real(kind=CUSTOM_REAL), intent(inout) :: b_RHS(:)

  ! Cross-gradient data.
  type(t_tau) :: tau
  integer :: i, j, k, l, p
  integer :: ind1, ind2, nderiv
  integer :: nsmaller
  real(kind=CUSTOM_REAL) :: val1, val2
  logical :: on_left_boundary, on_right_boundary

  ! Update the full models.
  call model1%update_full(.true., myrank, nbproc)
  call model2%update_full(.true., myrank, nbproc)

  ! Set the number of derivatives.
  nderiv = this%get_num_deriv(der_type)

  this%cost = 0._CUSTOM_REAL

  ! Number of parameters on ranks smaller than current one.
  nsmaller = get_nsmaller(this%nparams_loc, myrank, nbproc)

  do p = 1, this%nparams

    ! Grid (i,j,k)-index in the full grid.
    ! (Assume both models have the same grid indexes.)
    i = model1%grid_full%i_(p)
    j = model1%grid_full%j_(p)
    k = model1%grid_full%k_(p)

    ! Sanity check.
    if (i < 1 .or. i > this%size%x .or. &
        j < 1 .or. j > this%size%y .or. &
        k < 1 .or. k > this%size%z) then
        call exit_MPI("Wrong index in cross_gradient_calculate! i, j, k, p =" &
                      //str(i)//" "//str(j)//" "//str(k)//" "//str(p), myrank, 0)
    endif

    ! Calculate cross-gradient discretization.
    if (der_type == 1 .or. der_type == 2) then
      ! Apply boundary conditions.
      on_left_boundary = .false.
      on_right_boundary = .false.

      if (i == 1 .or. j == 1 .or. k == 1) on_left_boundary = .true.
      if (i == this%size%x .or. j == this%size%y .or. k == this%size%z) on_right_boundary = .true.

      if (der_type == 2 .and. on_left_boundary .and. on_right_boundary) then
      ! Skip the cross gradient constraints.
        tau%val = 0.d0
        tau%dm1 = t_vector(0.d0, 0.d0, 0.d0)
        tau%dm2 = t_vector(0.d0, 0.d0, 0.d0)

      else if (der_type == 2 .and. on_left_boundary) then
      ! On the left boundary: use forward difference.
        tau = this%calculate_tau(model1, model2, i, j, k, 1)

      else if (on_right_boundary) then
      ! On the right boundary: use backward difference.
        tau = this%calculate_tau_backward(model1, model2, i, j, k)

      else
      ! Inside the domain.
        tau = this%calculate_tau(model1, model2, i, j, k, der_type)
      endif

    else
      call exit_MPI("Unsupported derivative type!", myrank, der_type)
    endif

!  ! (Start) Skip the cross gradient constraints on the boundaries.
!  if (i == this%size%x .or. j == this%size%y .or. k == this%size%z) then
!    tau%val = 0.d0
!    tau%dm1 = t_vector(0.d0, 0.d0, 0.d0)
!    tau%dm2 = t_vector(0.d0, 0.d0, 0.d0)
!  endif
!
!  if ((der_type == 2) .and. (i == 1 .or. j == 1 .or. k == 1)) then
!    tau%val = 0.d0
!    tau%dm1 = t_vector(0.d0, 0.d0, 0.d0)
!    tau%dm2 = t_vector(0.d0, 0.d0, 0.d0)
!  endif
!  ! (End) Skip the cross gradient constraints on the boundaries.


    ! Normalization.
    !call this%normalize_tau(tau, model1, model2, i, j, k, myrank)

    ! Store cross gradient vector magnitude.
    this%cross_grad(p) = tau%val%get_norm()

    ! Calculate the cost.
    this%cost%x = this%cost%x + tau%val%x**2
    this%cost%y = this%cost%y + tau%val%y**2
    this%cost%z = this%cost%z + tau%val%z**2

    if (add) then
    ! Adding to the matrix and right-hand-side.

      ! Row with x-component.
      call matrix%new_row(myrank)

      b_RHS(matrix%get_current_row_number()) = - tau%val%x

      do l = 1, nderiv
        ind1 = tau%ind1(l)%x
        ind2 = tau%ind2(l)%x

        if (ind1 > nsmaller .and. ind1 <= nsmaller + this%nparams_loc) then
          ind1 = ind1 - nsmaller
          val1 = tau%dm1(l)%x * column_weight1(ind1)
          call matrix%add(val1, ind1, myrank)
        endif

        if (ind2 > nsmaller .and. ind2 <= nsmaller + this%nparams_loc) then
          ind2 = ind2 - nsmaller
          val2 = tau%dm2(l)%x * column_weight2(ind2)
          call matrix%add(val2, ind2 + this%nparams_loc, myrank)
        endif
      enddo

      ! Row with y-component.
      call matrix%new_row(myrank)

      b_RHS(matrix%get_current_row_number()) = - tau%val%y

      do l = 1, nderiv
        ind1 = tau%ind1(l)%y
        ind2 = tau%ind2(l)%y

        if (ind1 > nsmaller .and. ind1 <= nsmaller + this%nparams_loc) then
          ind1 = ind1 - nsmaller
          val1 = tau%dm1(l)%y * column_weight1(ind1)
          call matrix%add(val1, ind1, myrank)
        endif

        if (ind2 > nsmaller .and. ind2 <= nsmaller + this%nparams_loc) then
          ind2 = ind2 - nsmaller
          val2 = tau%dm2(l)%y * column_weight2(ind2)
          call matrix%add(val2, ind2 + this%nparams_loc, myrank)
        endif
      enddo

      ! Row with z-component.
      call matrix%new_row(myrank)

      b_RHS(matrix%get_current_row_number()) = - tau%val%z

      do l = 1, nderiv
        ind1 = tau%ind1(l)%z
        ind2 = tau%ind2(l)%z

        if (ind1 > nsmaller .and. ind1 <= nsmaller + this%nparams_loc) then
          ind1 = ind1 - nsmaller
          val1 = tau%dm1(l)%z * column_weight1(ind1)
          call matrix%add(val1, ind1, myrank)
        endif

        if (ind2 > nsmaller .and. ind2 <= nsmaller + this%nparams_loc) then
          ind2 = ind2 - nsmaller
          val2 = tau%dm2(l)%z * column_weight2(ind2)
          call matrix%add(val2, ind2 + this%nparams_loc, myrank)
        endif
      enddo
    endif

  enddo

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
pure function cross_gradient_get_magnitude(this) result(res)
  class(t_cross_gradient), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res(this%nparams)

  res = this%cross_grad

end function cross_gradient_get_magnitude

!==================================================================================================
! Calculates the discretized cross-gradients function (and corresponding matrix column indexes)
! between the models model1 and model2, at pixel location (i, j, k).
! der_type = 1:  using a forward difference scheme (Geophys. Res. Lett., Vol. 33, L07303).
! der_type /= 1: using a central difference scheme.
!==================================================================================================
function cross_gradient_calculate_tau(model1, model2, i, j, k, der_type) result(tau)
  type(t_model), intent(in) :: model1
  type(t_model), intent(in) :: model2
  integer, intent(in) :: i, j, k, der_type

  type(t_tau) :: tau
  type(t_vector) :: m1_grad, m2_grad
  type(t_vector) :: step
  integer :: ind

  ind = model1%grid_full%get_ind(i, j, k)

  ! Calculate model gradients.
  ! NOTE: These gradients are used to calculate the partial derivatives below,
  ! so if we want to e.g. normalize gradients then expressions for derivatives have to be changed.
  m1_grad = get_grad(model1%val_full(:, 1), model1%grid_full, i, j, k, get_der_type(der_type))
  m2_grad = get_grad(model2%val_full(:, 1), model2%grid_full, i, j, k, get_der_type(der_type))

  ! Calculate the cross-product between model gradients.
  tau%val = m1_grad%cross_product(m2_grad)

  ! Calculate partial derivatives (with respect to the model parameters),
  ! and corresponding matrix column indexes.

  step = t_vector(model1%grid_full%get_hx(ind), model1%grid_full%get_hy(ind), model1%grid_full%get_hz(ind))

  if (der_type /= 1) then
    step = 2._CUSTOM_REAL * step
  endif

  !------------------------------------------------------------------
  ! Note: assume both models have the same grid indexes.

  tau%dm1(1)%x =   m2_grad%z / step%y   ! dx / dm1(i, j + 1, k)
  tau%dm2(1)%x = - m1_grad%z / step%y   ! dx / dm2(i, j + 1, k)

  tau%dm1(2)%x = - m2_grad%y / step%z   ! dx / dm1(i, j, k + 1)
  tau%dm2(2)%x =   m1_grad%y / step%z   ! dx / dm2(i, j, k + 1)

  tau%ind(1)%x = model1%grid_full%get_ind(i, j + 1, k)
  tau%ind(2)%x = model1%grid_full%get_ind(i, j, k + 1)

  if (der_type == 1) then
    tau%dm1(3)%x = - (m2_grad%z / step%y - m2_grad%y / step%z)  ! dx / dm1(i, j, k)
    tau%dm2(3)%x = - (m1_grad%y / step%z - m1_grad%z / step%y)  ! dx / dm2(i, j, k)

    tau%ind(3)%x = model1%grid_full%get_ind(i, j, k)

  else
    tau%dm1(3)%x = - tau%dm1(1)%x   ! dx / dm1(i, j - 1, k)
    tau%dm2(3)%x = - tau%dm2(1)%x   ! dx / dm2(i, j - 1, k)

    tau%dm1(4)%x = - tau%dm1(2)%x   ! dx / dm1(i, j, k - 1)
    tau%dm2(4)%x = - tau%dm2(2)%x   ! dx / dm2(i, j, k - 1)

    tau%ind(3)%x = model1%grid_full%get_ind(i, j - 1, k)
    tau%ind(4)%x = model1%grid_full%get_ind(i, j, k - 1)
  endif

  !------------
  tau%dm1(1)%y = - m2_grad%z / step%x   ! dy / dm1(i + 1, j, k)
  tau%dm2(1)%y =   m1_grad%z / step%x   ! dy / dm2(i + 1, j, k)

  tau%dm1(2)%y =   m2_grad%x / step%z   ! dy / dm1(i, j, k + 1)
  tau%dm2(2)%y = - m1_grad%x / step%z   ! dy / dm2(i, j, k + 1)

  tau%ind(1)%y = model1%grid_full%get_ind(i + 1, j, k)
  tau%ind(2)%y = model1%grid_full%get_ind(i, j, k + 1)

  if (der_type == 1) then
    tau%dm1(3)%y = - (m2_grad%x / step%z - m2_grad%z / step%x)  ! dy / dm1(i, j, k)
    tau%dm2(3)%y = - (m1_grad%z / step%x - m1_grad%x / step%z)  ! dy / dm2(i, j, k)

    tau%ind(3)%y = model1%grid_full%get_ind(i, j, k)

  else
    tau%dm1(3)%y = - tau%dm1(1)%y   ! dy / dm1(i - 1, j, k)
    tau%dm2(3)%y = - tau%dm2(1)%y   ! dy / dm2(i - 1, j, k)

    tau%dm1(4)%y = - tau%dm1(2)%y   ! dy / dm1(i, j, k - 1)
    tau%dm2(4)%y = - tau%dm2(2)%y   ! dy / dm2(i, j, k - 1)

    tau%ind(3)%y = model1%grid_full%get_ind(i - 1, j, k)
    tau%ind(4)%y = model1%grid_full%get_ind(i, j, k - 1)
  endif

  !------------
  tau%dm1(1)%z =   m2_grad%y / step%x   ! dz / dm1(i + 1, j, k)
  tau%dm2(1)%z = - m1_grad%y / step%x   ! dz / dm2(i + 1, j, k)

  tau%dm1(2)%z = - m2_grad%x / step%y   ! dz / dm1(i, j + 1, k)
  tau%dm2(2)%z =   m1_grad%x / step%y   ! dz / dm2(i, j + 1, k)

  tau%ind(1)%z = model1%grid_full%get_ind(i + 1, j, k)
  tau%ind(2)%z = model1%grid_full%get_ind(i, j + 1, k)

  if (der_type == 1) then
    tau%dm1(3)%z = - (m2_grad%y / step%x - m2_grad%x / step%y)  ! dz / dm1(i, j, k)
    tau%dm2(3)%z = - (m1_grad%x / step%y - m1_grad%y / step%x)  ! dz / dm2(i, j, k)

    tau%ind(3)%z = model1%grid_full%get_ind(i, j, k)

  else
    tau%dm1(3)%z = - tau%dm1(1)%z   ! dz / dm1(i - 1, j, k)
    tau%dm2(3)%z = - tau%dm2(1)%z   ! dz / dm2(i - 1, j, k)

    tau%dm1(4)%z = - tau%dm1(2)%z   ! dz / dm1(i, j - 1, k)
    tau%dm2(4)%z = - tau%dm2(2)%z   ! dz / dm2(i, j - 1, k)

    tau%ind(3)%z = model1%grid_full%get_ind(i - 1, j, k)
    tau%ind(4)%z = model1%grid_full%get_ind(i, j - 1, k)
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

  tau%ind1 = tau%ind
  tau%ind2 = tau%ind

end function cross_gradient_calculate_tau

!==================================================================================================
! Calculates the discretized cross-gradients function (and corresponding matrix column indexes)
! between the models model1 and model2, at pixel location (i, j, k).
!
! Same as cross_gradient_calculate_tau but uses three-point forward difference scheme.
!==================================================================================================
function cross_gradient_calculate_tau2(model1, model2, i, j, k) result(tau)
  type(t_model), intent(in) :: model1
  type(t_model), intent(in) :: model2
  integer, intent(in) :: i, j, k

  type(t_tau) :: tau
  type(t_vector) :: m1_grad, m2_grad
  type(t_vector) :: step
  integer :: ind

  ind = model1%grid_full%get_ind(i, j, k)

  ! Calculate model gradients.
  ! NOTE: These gradients are used to calculate the partial derivatives below,
  ! so if we want to e.g. normalize gradients then expressions for derivatives have to be changed.

  m1_grad = get_grad(model1%val_full(:, 1), model1%grid_full, i, j, k, FWD2_TYPE)
  m2_grad = get_grad(model2%val_full(:, 1), model2%grid_full, i, j, k, FWD2_TYPE)

  ! Calculate the cross-product between model gradients.
  tau%val = m1_grad%cross_product(m2_grad)

  ! Calculate partial derivatives (with respect to the model parameters),
  ! and corresponding matrix column indexes.

  step = 2.d0 * t_vector(model1%grid_full%get_hx(ind), model1%grid_full%get_hy(ind), model1%grid_full%get_hz(ind))

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

  tau%ind(1)%x = model1%grid_full%get_ind(i, j + 2, k)
  tau%ind(2)%x = model1%grid_full%get_ind(i, j + 1, k)
  tau%ind(3)%x = model1%grid_full%get_ind(i, j, k + 2)
  tau%ind(4)%x = model1%grid_full%get_ind(i, j, k + 1)
  tau%ind(5)%x = model1%grid_full%get_ind(i, j, k)

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

  tau%ind(1)%y = model1%grid_full%get_ind(i, j, k + 2)
  tau%ind(2)%y = model1%grid_full%get_ind(i, j, k + 1)
  tau%ind(3)%y = model1%grid_full%get_ind(i + 2, j, k)
  tau%ind(4)%y = model1%grid_full%get_ind(i + 1, j, k)
  tau%ind(5)%y = model1%grid_full%get_ind(i, j, k)

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

  tau%ind(1)%z = model1%grid_full%get_ind(i + 2, j, k)
  tau%ind(2)%z = model1%grid_full%get_ind(i + 1, j, k)
  tau%ind(3)%z = model1%grid_full%get_ind(i, j + 2, k)
  tau%ind(4)%z = model1%grid_full%get_ind(i, j + 1, k)
  tau%ind(5)%z = model1%grid_full%get_ind(i, j, k)

  tau%ind1 = tau%ind
  tau%ind2 = tau%ind

end function cross_gradient_calculate_tau2

!==================================================================================================
! Calculates the discretized cross-gradients function (and corresponding matrix column indexes)
! between the models model1 and model2, at pixel location (i, j, k).
!
! Same as cross_gradient_calculate_tau but uses backward finite difference.
!==================================================================================================
function cross_gradient_calculate_tau_backward(model1, model2, i, j, k) result(tau)
  type(t_model), intent(in) :: model1
  type(t_model), intent(in) :: model2
  integer, intent(in) :: i, j, k

  type(t_tau) :: tau
  type(t_vector) :: m1_grad, m2_grad
  type(t_vector) :: step
  integer :: ind

  ind = model1%grid_full%get_ind(i, j, k)

  ! Calculate model gradients.
  ! NOTE: These gradients are used to calculate the partial derivatives below,
  ! so if we want to e.g. normalize gradients then expressions for derivatives have to be changed.

  ! Forward and backward difference.
  m1_grad = get_grad(model1%val_full(:, 1), model1%grid_full, i, j, k, BWD_TYPE)
  m2_grad = get_grad(model2%val_full(:, 1), model2%grid_full, i, j, k, BWD_TYPE)

  ! Calculate the cross-product between model gradients.
  tau%val = m1_grad%cross_product(m2_grad)

  ! Calculate partial derivatives (with respect to the model parameters),
  ! and corresponding matrix column indexes.

  step = t_vector(model1%grid_full%get_hx(ind), model1%grid_full%get_hy(ind), model1%grid_full%get_hz(ind))

  !------------------------------------------------------------------
  ! Note: assume both models have the same grid indexes.

  tau%dm1(1)%x = - m2_grad%z / step%y                         ! dx / dm1(i, j - 1, k)
  tau%dm1(2)%x =   m2_grad%y / step%z                         ! dx / dm1(i, j, k - 1)
  tau%dm1(3)%x =   m2_grad%z / step%y - m2_grad%y / step%z    ! dx / dm1(i, j, k)

  tau%dm2(1)%x =   m1_grad%z / step%y                         ! dx / dm2(i, j - 1, k)
  tau%dm2(2)%x = - m1_grad%y / step%z                         ! dx / dm2(i, j, k - 1)
  tau%dm2(3)%x =   m1_grad%y / step%z - m1_grad%z / step%y    ! dx / dm2(i, j, k)

  tau%ind(1)%x = model1%grid_full%get_ind(i, j - 1, k)
  tau%ind(2)%x = model1%grid_full%get_ind(i, j, k - 1)
  tau%ind(3)%x = model1%grid_full%get_ind(i, j, k)

  !------------
  tau%dm1(1)%y =   m2_grad%z / step%x                         ! dy / dm1(i - 1, j, k)
  tau%dm1(2)%y = - m2_grad%x / step%z                         ! dy / dm1(i, j, k - 1)
  tau%dm1(3)%y =   m2_grad%x / step%z - m2_grad%z / step%x    ! dy / dm1(i, j, k)

  tau%dm2(1)%y = - m1_grad%z / step%x                         ! dy / dm2(i - 1, j, k)
  tau%dm2(2)%y =   m1_grad%x / step%z                         ! dy / dm2(i, j, k - 1)
  tau%dm2(3)%y =   m1_grad%z / step%x - m1_grad%x / step%z    ! dy / dm2(i, j, k)

  tau%ind(1)%y = model1%grid_full%get_ind(i - 1, j, k)
  tau%ind(2)%y = model1%grid_full%get_ind(i, j, k - 1)
  tau%ind(3)%y = model1%grid_full%get_ind(i, j, k)

  !------------
  tau%dm1(1)%z = - m2_grad%y / step%x                         ! dz / dm1(i - 1, j, k)
  tau%dm1(2)%z =   m2_grad%x / step%y                         ! dz / dm1(i, j - 1, k)
  tau%dm1(3)%z =   m2_grad%y / step%x - m2_grad%x / step%y    ! dz / dm1(i, j, k)

  tau%dm2(1)%z =   m1_grad%y / step%x                         ! dz / dm2(i - 1, j, k)
  tau%dm2(2)%z = - m1_grad%x / step%y                         ! dz / dm2(i, j - 1, k)
  tau%dm2(3)%z =   m1_grad%x / step%y - m1_grad%y / step%x    ! dz / dm2(i, j, k)

  tau%ind(1)%z = model1%grid_full%get_ind(i - 1, j, k)
  tau%ind(2)%z = model1%grid_full%get_ind(i, j - 1, k)
  tau%ind(3)%z = model1%grid_full%get_ind(i, j, k)

  ! To avoid having non-initialized values.
  tau%dm1(4:) = t_vector(0.d0, 0.d0, 0.d0)
  tau%dm2(4:) = t_vector(0.d0, 0.d0, 0.d0)
  tau%ind(4:) = t_ivector(0, 0, 0)

  tau%ind1 = tau%ind
  tau%ind2 = tau%ind

end function cross_gradient_calculate_tau_backward


!==================================================================================================
! Calculates the discretized cross-gradients function (and corresponding matrix column indexes)
! between the models model1 and model2, at pixel location (i, j, k).
!
! Same as cross_gradient_calculate_tau but uses forward and central difference at the same time.
!==================================================================================================
function cross_gradient_calculate_tau_mixed_gradients(model1, model2, i, j, k) result(tau)
  type(t_model), intent(in) :: model1
  type(t_model), intent(in) :: model2
  integer, intent(in) :: i, j, k

  type(t_tau) :: tau
  type(t_vector) :: m1_grad, m2_grad
  type(t_vector) :: step
  integer :: ind

  ind = model1%grid_full%get_ind(i, j, k)

  ! Calculate model gradients.
  ! NOTE: These gradients are used to calculate the partial derivatives below,
  ! so if we want to e.g. normalize gradients then expressions for derivatives have to be changed.

  ! Forward and central difference.
!  m1_grad = grad%get_grad(model1%val_full(:, 1), model1%grid_full, i, j, k, FWD_TYPE)
!  m2_grad = grad%get_grad(model2%val_full(:, 1), model2%grid_full, i, j, k, CNT_TYPE)

  ! Forward and backward difference.
  m1_grad = get_grad(model1%val_full(:, 1), model1%grid_full, i, j, k, FWD_TYPE)
  m2_grad = get_grad(model2%val_full(:, 1), model2%grid_full, i, j, k, BWD_TYPE)

  ! Calculate the cross-product between model gradients.
  tau%val = m1_grad%cross_product(m2_grad)

  ! Calculate partial derivatives (with respect to the model parameters),
  ! and corresponding matrix column indexes.

  step = t_vector(model1%grid_full%get_hx(ind), model1%grid_full%get_hy(ind), model1%grid_full%get_hz(ind))

  !------------------------------------------------------------------
  ! Note: assume both models have the same grid indexes.

!-------------------------------------------------------------------------------------------
!  Forward and central difference.
!-------------------------------------------------------------------------------------------

!  tau%dm1(1)%x =   m2_grad%z / step%y                         ! dx / dm1(i, j + 1, k)
!  tau%dm1(2)%x = - m2_grad%y / step%z                         ! dx / dm1(i, j, k + 1)
!  tau%dm1(3)%x = - (m2_grad%z / step%y - m2_grad%y / step%z)  ! dx / dm1(i, j, k)
!
!  tau%dm2(1)%x = - m1_grad%z / step%y / 2.d0                  ! dx / dm2(i, j + 1, k)
!  tau%dm2(2)%x =   m1_grad%y / step%z / 2.d0                  ! dx / dm2(i, j, k + 1)
!  tau%dm2(3)%x = - tau%dm2(1)%x                               ! dx / dm2(i, j - 1, k)
!  tau%dm2(4)%x = - tau%dm2(2)%x                               ! dx / dm2(i, j, k - 1)
!
!  tau%ind1(1)%x = model1%grid_full%get_ind(i, j + 1, k)
!  tau%ind1(2)%x = model1%grid_full%get_ind(i, j, k + 1)
!  tau%ind1(3)%x = model1%grid_full%get_ind(i, j, k)
!
!  tau%ind2(1)%x = model1%grid_full%get_ind(i, j + 1, k)
!  tau%ind2(2)%x = model1%grid_full%get_ind(i, j, k + 1)
!  tau%ind2(3)%x = model1%grid_full%get_ind(i, j - 1, k)
!  tau%ind2(4)%x = model1%grid_full%get_ind(i, j, k - 1)
!
!  !------------
!  tau%dm1(1)%y = - m2_grad%z / step%x                         ! dy / dm1(i + 1, j, k)
!  tau%dm1(2)%y =   m2_grad%x / step%z                         ! dy / dm1(i, j, k + 1)
!  tau%dm1(3)%y = - (m2_grad%x / step%z - m2_grad%z / step%x)  ! dy / dm1(i, j, k)
!
!  tau%dm2(1)%y =   m1_grad%z / step%x / 2.d0                  ! dy / dm2(i + 1, j, k)
!  tau%dm2(2)%y = - m1_grad%x / step%z / 2.d0                  ! dy / dm2(i, j, k + 1)
!  tau%dm2(3)%y = - tau%dm2(1)%y                               ! dy / dm2(i - 1, j, k)
!  tau%dm2(4)%y = - tau%dm2(2)%y                               ! dy / dm2(i, j, k - 1)
!
!  tau%ind1(1)%y = model1%grid_full%get_ind(i + 1, j, k)
!  tau%ind1(2)%y = model1%grid_full%get_ind(i, j, k + 1)
!  tau%ind1(3)%y = model1%grid_full%get_ind(i, j, k)
!
!  tau%ind2(1)%y = model1%grid_full%get_ind(i + 1, j, k)
!  tau%ind2(2)%y = model1%grid_full%get_ind(i, j, k + 1)
!  tau%ind2(3)%y = model1%grid_full%get_ind(i - 1, j, k)
!  tau%ind2(4)%y = model1%grid_full%get_ind(i, j, k - 1)
!
!  !------------
!  tau%dm1(1)%z =   m2_grad%y / step%x                         ! dz / dm1(i + 1, j, k)
!  tau%dm1(2)%z = - m2_grad%x / step%y                         ! dz / dm1(i, j + 1, k)
!  tau%dm1(3)%z = - (m2_grad%y / step%x - m2_grad%x / step%y)  ! dz / dm1(i, j, k)
!
!  tau%dm2(1)%z = - m1_grad%y / step%x / 2.d0                  ! dz / dm2(i + 1, j, k)
!  tau%dm2(2)%z =   m1_grad%x / step%y / 2.d0                  ! dz / dm2(i, j + 1, k)
!  tau%dm2(3)%z = - tau%dm2(1)%z                               ! dz / dm2(i - 1, j, k)
!  tau%dm2(4)%z = - tau%dm2(2)%z                               ! dz / dm2(i, j - 1, k)
!
!  tau%ind1(1)%z = model1%grid_full%get_ind(i + 1, j, k)
!  tau%ind1(2)%z = model1%grid_full%get_ind(i, j + 1, k)
!  tau%ind1(3)%z = model1%grid_full%get_ind(i, j, k)
!
!  tau%ind2(1)%z = model1%grid_full%get_ind(i + 1, j, k)
!  tau%ind2(2)%z = model1%grid_full%get_ind(i, j + 1, k)
!  tau%ind2(3)%z = model1%grid_full%get_ind(i - 1, j, k)
!  tau%ind2(4)%z = model1%grid_full%get_ind(i, j - 1, k)
!
!  ! To avoid having non-initialized values.
!  tau%ind1(4)%x = 0
!  tau%ind1(4)%y = 0
!  tau%ind1(4)%z = 0

!  tau%dm1(4) = 0.d0

!-------------------------------------------------------------------------------------------
!  Forward and backward difference.
!-------------------------------------------------------------------------------------------

  tau%dm1(1)%x =   m2_grad%z / step%y                         ! dx / dm1(i, j + 1, k)
  tau%dm1(2)%x = - m2_grad%y / step%z                         ! dx / dm1(i, j, k + 1)
  tau%dm1(3)%x = - (m2_grad%z / step%y - m2_grad%y / step%z)  ! dx / dm1(i, j, k)

  tau%dm2(1)%x =   m1_grad%z / step%y                         ! dx / dm2(i, j - 1, k)
  tau%dm2(2)%x = - m1_grad%y / step%z                         ! dx / dm2(i, j, k - 1)
  tau%dm2(3)%x =   m1_grad%y / step%z - m1_grad%z / step%y    ! dx / dm2(i, j, k)

  tau%ind1(1)%x = model1%grid_full%get_ind(i, j + 1, k)
  tau%ind1(2)%x = model1%grid_full%get_ind(i, j, k + 1)
  tau%ind1(3)%x = model1%grid_full%get_ind(i, j, k)

  tau%ind2(1)%x = model1%grid_full%get_ind(i, j - 1, k)
  tau%ind2(2)%x = model1%grid_full%get_ind(i, j, k - 1)
  tau%ind2(3)%x = model1%grid_full%get_ind(i, j, k)

  !------------
  tau%dm1(1)%y = - m2_grad%z / step%x                         ! dy / dm1(i + 1, j, k)
  tau%dm1(2)%y =   m2_grad%x / step%z                         ! dy / dm1(i, j, k + 1)
  tau%dm1(3)%y = - (m2_grad%x / step%z - m2_grad%z / step%x)  ! dy / dm1(i, j, k)

  tau%dm2(1)%y = - m1_grad%z / step%x                         ! dy / dm2(i - 1, j, k)
  tau%dm2(2)%y =   m1_grad%x / step%z                         ! dy / dm2(i, j, k - 1)
  tau%dm2(3)%y =   m1_grad%z / step%x - m1_grad%x / step%z    ! dy / dm2(i, j, k)

  tau%ind1(1)%y = model1%grid_full%get_ind(i + 1, j, k)
  tau%ind1(2)%y = model1%grid_full%get_ind(i, j, k + 1)
  tau%ind1(3)%y = model1%grid_full%get_ind(i, j, k)

  tau%ind2(1)%y = model1%grid_full%get_ind(i - 1, j, k)
  tau%ind2(2)%y = model1%grid_full%get_ind(i, j, k - 1)
  tau%ind2(3)%y = model1%grid_full%get_ind(i, j, k)

  !------------
  tau%dm1(1)%z =   m2_grad%y / step%x                         ! dz / dm1(i + 1, j, k)
  tau%dm1(2)%z = - m2_grad%x / step%y                         ! dz / dm1(i, j + 1, k)
  tau%dm1(3)%z = - (m2_grad%y / step%x - m2_grad%x / step%y)  ! dz / dm1(i, j, k)

  tau%dm2(1)%z =   m1_grad%y / step%x                         ! dz / dm2(i - 1, j, k)
  tau%dm2(2)%z = - m1_grad%x / step%y                         ! dz / dm2(i, j - 1, k)
  tau%dm2(3)%z =   m1_grad%x / step%y - m1_grad%y / step%x    ! dz / dm2(i, j, k)

  tau%ind1(1)%z = model1%grid_full%get_ind(i + 1, j, k)
  tau%ind1(2)%z = model1%grid_full%get_ind(i, j + 1, k)
  tau%ind1(3)%z = model1%grid_full%get_ind(i, j, k)

  tau%ind2(1)%z = model1%grid_full%get_ind(i - 1, j, k)
  tau%ind2(2)%z = model1%grid_full%get_ind(i, j - 1, k)
  tau%ind2(3)%z = model1%grid_full%get_ind(i, j, k)

  ! To avoid having non-initialized values.
  tau%ind1(4)%x = 0
  tau%ind1(4)%y = 0
  tau%ind1(4)%z = 0

  tau%ind2(4)%x = 0
  tau%ind2(4)%y = 0
  tau%ind2(4)%z = 0

  tau%dm1(4) = 0.d0
  tau%dm2(4) = 0.d0

end function cross_gradient_calculate_tau_mixed_gradients

!==============================================================================================
! Normalize the cross gradient function.
! (See N Linde, J Doetsch, Joint inversion of crosshole GPR and Seismic traveltime data.)
!==============================================================================================
subroutine cross_gradient_normalize_tau(tau, model1, model2, i, j, k)
  type(t_model), intent(in) :: model1
  type(t_model), intent(in) :: model2
  integer, intent(in) :: i, j, k

  type(t_tau), intent(inout) :: tau
  real(kind=CUSTOM_REAL) :: val1, val2, scale
  integer :: l

  val1 = grad_get_par(model1%val_full(:, 1), model1%grid_full, i, j, k)
  val2 = grad_get_par(model2%val_full(:, 1), model2%grid_full, i, j, k)

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
