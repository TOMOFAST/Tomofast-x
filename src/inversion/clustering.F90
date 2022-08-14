
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
! A class to add clustering constraints to joint inversion.
!
! Vitaliy Ogarko and Jeremie Giraud, UWA, CET, Australia, 2016-2017.
!===========================================================================================
module clustering

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use sparse_matrix
  use model
  use grid
  use string, only: str
  use parallel_tools

  implicit none

  private

  type, public :: t_clustering_data
    ! Likelihood value P(m).
    real(kind=CUSTOM_REAL) :: value
    ! Derivatives for grav and mag problems.
    real(kind=CUSTOM_REAL) :: derivative(2)
    ! The value of each mixture element (Gaussian), i.e., one per cluster.
    real(kind=CUSTOM_REAL), allocatable :: component(:)
  end type

  type, public :: t_clustering
    private

    ! Clustering weight (on the clustering term in the cost function).
    real(kind=CUSTOM_REAL) :: weight_glob(2)

    ! Clustering weights (for individual problems).
    real(kind=CUSTOM_REAL) :: weight_loc(2)

    ! Number of mixtures for clustering constraints.
    integer :: nclusters

    ! Type of optimization (normal or logarithmic constraints).
    integer :: optimization_type

    ! Type of constraints (local - different per cell vs global same for all cells).
    integer :: constraints_type

    ! Mixture parameters (Gaussian).
    ! First dimension corresponds to model type,
    ! i.e., mixture_mu(1, :) for sigma_xx, and mixture_mu(2, :) for sigma_yy, and mixture_mu(3, :) for sigma_xy.
    ! Second dimension corresponds to different clusters.
    real(kind=CUSTOM_REAL), allocatable :: mixture_mu(:, :)
    real(kind=CUSTOM_REAL), allocatable :: mixture_sigma(:, :)
    ! Weights per cluster per cell.
    real(kind=CUSTOM_REAL), allocatable :: cell_weight(:, :)

    ! The maximum value of the mixture (at every cell).
    real(kind=CUSTOM_REAL), allocatable :: mixture_max(:)

    ! Model dimensions.
    integer :: nx, ny, nz
    ! Number of model parameters (local per CPU and total).
    integer :: nelements, nelements_total

    ! Clustering data.
    type(t_clustering_data), allocatable :: data_all(:)

    ! Cost.
    real(kind=CUSTOM_REAL) :: cost(2)

  contains
    private

    procedure, public, pass :: initialize => clustering_initialize
    procedure, public, pass :: read_mixtures => clustering_read_mixtures
    procedure, public, pass :: write_mixtures => clustering_write_mixtures
    procedure, public, pass :: write_data => clustering_write_data

    procedure, public, pass :: add => clustering_add
    procedure, public, pass :: get_cost => clustering_get_cost
    procedure, public, pass :: get_probabilities => clustering_get_probabilities

    procedure, private, nopass :: calculate_Gaussian => clustering_calculate_Gaussian
    procedure, private, nopass :: calculate_Gaussian_1D => clustering_calculate_Gaussian_1D
    procedure, private, nopass :: calculate_Gaussian_2D => clustering_calculate_Gaussian_2D

    procedure, private, pass :: calculate_Gaussian_mixture => clustering_calculate_Gaussian_mixture
    procedure, private, pass :: calculate_Gaussian_mixture_max => clustering_calculate_Gaussian_mixture_max

  end type t_clustering

contains

!===========================================================================================
! Initialization.
!===========================================================================================
subroutine clustering_initialize(this, nelements_total, nelements, weight_glob, &
                                 nclusters, optimization_type, constraints_type, myrank)
  class(t_clustering), intent(inout) :: this
  integer, intent(in) :: nelements_total, nelements
  real(kind=CUSTOM_REAL), intent(in) :: weight_glob(2)
  integer, intent(in) :: nclusters, optimization_type, constraints_type, myrank

  integer :: ierr, i

  if (myrank == 0) print *, "Initializing clustering constraints."

  this%nelements_total = nelements_total
  this%nelements = nelements
  this%weight_glob = weight_glob
  this%nclusters = nclusters
  this%optimization_type = optimization_type
  this%constraints_type = constraints_type

  this%cost = 0.d0

  do i = 1, 2
  ! This way we use 1D Gaussians when one of the global clustering weight is zero.
    if (weight_glob(i) == 0.d0) then
      this%weight_loc(i) = 0.d0
    else
      this%weight_loc(i) = 1.d0
    endif
  enddo

  ierr = 0

  if (.not. allocated(this%mixture_mu)) allocate(this%mixture_mu(2, this%nclusters), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%mixture_sigma)) allocate(this%mixture_sigma(3, this%nclusters), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%cell_weight)) &
    allocate(this%cell_weight(this%nelements_total, this%nclusters), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%mixture_max)) allocate(this%mixture_max(this%nelements_total), source=0._CUSTOM_REAL, stat=ierr)

  if (.not. allocated(this%data_all)) then
    allocate(this%data_all(this%nelements_total), stat=ierr)

    do i = 1, this%nelements_total
      ! Allocate memory for component array at every cell.
      allocate(this%data_all(i)%component(this%nclusters), source=0._CUSTOM_REAL, stat=ierr)
    enddo
  endif

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in clustering_initialize!", myrank, ierr)

end subroutine clustering_initialize

!================================================================================================
! Read clustering mixture parameters.
!================================================================================================
subroutine clustering_read_mixtures(this, file_clusters, file_weights, myrank)
  class(t_clustering), intent(inout) :: this
  character(len=*), intent(in) :: file_clusters, file_weights
  integer, intent(in) :: myrank

  integer :: i, nclusters_read
  integer :: p, nelements_total_read
  integer :: ierr
  character(len=256) :: msg
  real(kind=CUSTOM_REAL), allocatable :: cluster_weight(:)
  real(kind=CUSTOM_REAL), allocatable :: weight(:)

  if (myrank == 0) then
  ! Reading by master CPU only.

    if (.not. allocated(cluster_weight)) allocate(cluster_weight(this%nclusters), source=0._CUSTOM_REAL, stat=ierr)

    print *, 'Reading clustering parameters from file ', trim(file_clusters)

    open(10, file=trim(file_clusters), status='old', action='read', iostat=ierr, iomsg=msg)
    if (ierr /= 0) call exit_MPI("Error in opening the mixture file! path=" &
                               //trim(file_clusters)//" iomsg="//msg, myrank, ierr)

    read(10, *, iostat=ierr) nclusters_read

    ! Sanity check.
    if (this%nclusters /= nclusters_read) &
      call exit_MPI("The number of clusters is inconsistent!"//new_line('a') &
                    //"nclusters="//str(this%nclusters)//new_line('a') &
                    //"nclusters_read="//str(nclusters_read), myrank, 0)

    ! Reading the mixture parameters.
    do i = 1, this%nclusters
      read(10, *, iostat=ierr) cluster_weight(i), &
                               this%mixture_mu(1, i), this%mixture_sigma(1, i), &
                               this%mixture_mu(2, i), this%mixture_sigma(2, i), &
                               this%mixture_sigma(3, i)

      print *, i, cluster_weight(i), &
                  this%mixture_mu(1, i), this%mixture_sigma(1, i), &
                  this%mixture_mu(2, i), this%mixture_sigma(2, i), &
                  this%mixture_sigma(3, i)

      if (ierr /= 0) call exit_MPI("Problem while reading the mixture file in clustering_read_mixtures!", myrank, ierr)
    enddo
    close(10)

    ! Normalize weights.
    cluster_weight = cluster_weight / sum(cluster_weight)

    !------------------------------------------------------------
    if (this%constraints_type /= 1) then
    ! Read local cluster weights (per cell).

      print *, 'Reading clustering cell-weights:'

      open(20, file=trim(file_weights), status='old', action='read', iostat=ierr, iomsg=msg)
      if (ierr /= 0) call exit_MPI("Error in opening the cell-weights file! path=" &
                                 //file_weights//" iomsg="//msg, myrank, ierr)

      read(20, *, iostat=ierr) nelements_total_read, nclusters_read

      ! Sanity checks.
      if (this%nclusters /= nclusters_read) &
        call exit_MPI("The number of clusters is inconsistent!"//new_line('a') &
                      //"nclusters="//str(this%nclusters)//new_line('a') &
                      //"nclusters_read="//str(nclusters_read), myrank, 0)

      if (this%nelements_total /= nelements_total_read) &
        call exit_MPI("The number of cells is inconsistent!"//new_line('a') &
                      //"nelements_total_read="//str(nelements_total_read)//new_line('a'), myrank, 0)

      allocate(weight(this%nclusters))

      ! Reading the mixture parameters.
      do p = 1, this%nelements_total
        read(20, *) weight

        this%cell_weight(p, :) = weight
      enddo
      close(20)

      deallocate(weight)

    else
    ! Set global cluster weights (same for all cells).

      do i = 1, this%nclusters
        this%cell_weight(:, i) = cluster_weight(i)
      enddo
    endif

    deallocate(cluster_weight)

  endif ! myrank == 0

  !-----------------------------------------------------------------------------------------
  ! Broadcast the mixture parameter arrays to all CPUs.
  ierr = 0

  call MPI_Bcast(this%mixture_mu, 2 * this%nclusters, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%mixture_sigma, 3 * this%nclusters, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%cell_weight, this%nelements_total * this%nclusters, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("Error in MPI_Bcast in clustering_read_mixtures!", myrank, ierr)

  !--------------------------------------------------------

  ! Calculating the maximum of Gaussian mixture.
  do p = 1, this%nelements_total
    this%mixture_max(p) = this%calculate_Gaussian_mixture_max(this%cell_weight(p, :))
  enddo

  if (myrank == 0) print *, 'Clustering mixture_max =', maxval(this%mixture_max)

end subroutine clustering_read_mixtures

!==================================================================================================
! Returns the cost for every component (this is what we want to minimize).
!==================================================================================================
pure function clustering_get_cost(this, problem_type) result(res)
  class(t_clustering), intent(in) :: this
  integer, intent(in) :: problem_type
  real(kind=CUSTOM_REAL) :: res

  res = this%cost(problem_type)

end function clustering_get_cost

!==================================================================================================
! Returns the probabilities P(m) for every cell.
!==================================================================================================
pure function clustering_get_probabilities(this) result(res)
  class(t_clustering), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res(this%nelements_total)

  res = this%data_all(:)%value

end function clustering_get_probabilities

!================================================================================================
! Write clustering mixtures values (to plot the function).
! Also calculate the integral of the function.
!================================================================================================
subroutine clustering_write_mixtures(this, file_name, myrank)
  class(t_clustering), intent(in) :: this
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL) :: gauss, deriv(2)
  real(kind=CUSTOM_REAL) :: gauss_loc(this%nclusters)
  real(kind=CUSTOM_REAL) :: dx(2), xmin(2), xmax(2), x(2)
  real(kind=CUSTOM_REAL) :: integral

  ! Skipped.
  return

  dx(1) = 1.d0
  xmin(1) = - 200.d0
  xmax(1) = 700.d0

  dx(2) = 1.d-4
  xmin(2) = - 0.02d0
  xmax(2) = 0.1d0

  integral = 0.d0

  if (myrank == 0) then
  ! Writing by master CPU only.
    open(10, file=trim(file_name), access='stream', form='formatted', status='replace', action='write')

    x = xmin

    do while (x(1) < xmax(1))
      do while (x(2) < xmax(2))

        call this%calculate_Gaussian_mixture(x, gauss, deriv, gauss_loc, this%cell_weight(1, :))

        write(10, *) x(1), x(2), gauss

        integral = integral + dx(1) * dx(2) * gauss

        x(2) = x(2) + dx(2)
      enddo
      x(1) = x(1) + dx(1)
      x(2) = xmin(2)
    enddo

    close(10)

    print *, 'Gaussian mixture 2D integral =', integral
  endif

end subroutine clustering_write_mixtures

!================================================================================================
! Write clustering data (probabilities and derivatives).
!================================================================================================
subroutine clustering_write_data(this, file_name, grid, myrank)
  class(t_clustering), intent(in) :: this
  character(len=*), intent(in) :: file_name
  type(t_grid), intent(in) :: grid
  integer, intent(in) :: myrank

  integer :: p
  character(len=256) :: filename_full

  if (myrank == 0) then
  ! Writing by master CPU only.

    call execute_command_line('mkdir -p '//trim(path_output)//"/Voxet/")

    filename_full = trim(path_output)//"/Voxet/"//file_name

    print *, 'Writing the clustering data to file ', trim(filename_full)

    open(37, file=trim(filename_full), access='stream', form='formatted', status='replace', action='write')

    write (37, *) this%nelements_total

    do p = 1, this%nelements_total
      write(37, *) grid%X1(p), grid%X2(p), grid%Y1(p), grid%Y2(p), grid%Z1(p), grid%Z2(p), &
                   grid%i_(p), grid%j_(p), grid%k_(p), &
                   this%data_all(p)%value, this%data_all(p)%derivative, this%data_all(p)%component
    enddo

    close(37)
  endif

end subroutine clustering_write_data

!===========================================================================================
! Adding clustering constraints to the matrix and right-hand-side.
!===========================================================================================
subroutine clustering_add(this, model1, model2, column_weight1, column_weight2, &
                          matrix, b_RHS, problem_type, myrank, nbproc)
  class(t_clustering), intent(inout) :: this
  type(t_model), intent(in) :: model1
  type(t_model), intent(in) :: model2
  real(kind=CUSTOM_REAL), intent(in) :: column_weight1(:)
  real(kind=CUSTOM_REAL), intent(in) :: column_weight2(:)
  integer, intent(in) :: problem_type
  integer, intent(in) :: myrank, nbproc

  type(t_sparse_matrix), intent(inout) :: matrix
  real(kind=CUSTOM_REAL), intent(inout) :: b_RHS(:)

  real(kind=CUSTOM_REAL) :: matrix_val(2)
  real(kind=CUSTOM_REAL) :: model_val(2)
  real(kind=CUSTOM_REAL) :: gauss, deriv(2)
  real(kind=CUSTOM_REAL) :: gauss_loc(this%nclusters)
  real(kind=CUSTOM_REAL) :: Cp_weight(2)
  real(kind=CUSTOM_REAL) :: func_val
  integer :: row_beg, row_end, nsmaller
  integer :: i, p, ind

  ! Calculate 'Cp-weights'.
  do i = 1, 2
    !Cp_weight(i) = maxval(this%mixture_sigma(i, :))
    Cp_weight(i) = 1.d0

    if (this%weight_loc(i) == 0.d0) then
    ! Single domain inversion.
      Cp_weight(i) = 0.d0
    endif
  enddo

  ! Number of parameters on ranks smaller than current one.
  nsmaller = get_nsmaller(this%nelements, myrank, nbproc)

  ! First matrix row (in the big matrix).
  row_beg = matrix%get_current_row_number() + 1

  ! Add matrix lines (Jacobian).
  do p = 1, this%nelements_total

    call matrix%new_row(myrank)

    model_val(1) = model1%val_full(p)
    model_val(2) = model2%val_full(p)

    call this%calculate_Gaussian_mixture(model_val, gauss, deriv, gauss_loc, this%cell_weight(p, :))

    if (this%optimization_type == 2) then
    ! Minimizing g(x) = - log(f(x)), so derivative  g'(x) = - f'(x) / f(x).
      if (gauss /= 0.d0) then
        deriv = - deriv / gauss
      else
        deriv = 0.d0
      endif
    endif

    if (p > nsmaller .and. p <= nsmaller + this%nelements) then
    ! The current element 'p' belongs to the current CPU, so adding it into the matrix.
      ! Calculate element local index.
      ind = p - nsmaller

      matrix_val(1) = this%weight_glob(1) * column_weight1(ind) * deriv(1) * Cp_weight(1)
      matrix_val(2) = this%weight_glob(2) * column_weight2(ind) * deriv(2) * Cp_weight(2)

      if (problem_type == 1) then
        call matrix%add(matrix_val(1), ind, myrank)

      else if (problem_type == 2) then
        call matrix%add(matrix_val(2), ind + this%nelements, myrank)
      endif
    endif

    ! Calculate function value.
    if (this%optimization_type == 1) then
    ! Normal case: minimizing ||A_max - f(x)||^2.

      ! Note: for least-squares we can switch the sign inside breakets.
      func_val = (gauss - this%mixture_max(p))

    else if (this%optimization_type == 2) then
    ! Logarithmic case: minimizing ||- log(f(x))||^2

      if (gauss > 0.d0) then
        func_val = - log(gauss) + log(this%mixture_max(p))
      else
      ! Exclude the clustering constraints for this pixel.
        func_val =  0.d0
      endif

    else
      call exit_MPI("Wrong optimization type in clustering_add!", myrank, this%optimization_type)
    endif

    ! Setting the right-hand side.
    b_RHS(matrix%get_current_row_number()) = - this%weight_glob(problem_type) * func_val * Cp_weight(problem_type)

    ! Store data for post-processing.
    this%data_all(p)%value = gauss
    this%data_all(p)%derivative = deriv
    this%data_all(p)%component = gauss_loc
  enddo

  ! Last matrix row (in the big matrix).
  row_end = matrix%get_current_row_number()

  ! Calculate the clustering cost.
  this%cost(problem_type) = sum(b_RHS(row_beg:row_end)**2)

  if (myrank == 0) print *, 'Clustering added lines: ', row_end - row_beg + 1

end subroutine clustering_add

!=================================================================================================
! Calculate the 1D Gaussian.
! Returns argument of the exponential and the norm: G = exp(arg) / norm.
!=================================================================================================
subroutine clustering_calculate_Gaussian_1D(x, mu, sigma, arg, norm)
  real(kind=CUSTOM_REAL), intent(in) :: x, mu, sigma
  real(kind=CUSTOM_REAL), intent(out) :: arg, norm

  arg = (- (x - mu)**2 / sigma**2 / 2.d0)
  norm = sqrt(2.d0 * PI * sigma**2)

end subroutine clustering_calculate_Gaussian_1D


!=================================================================================================
! Calculate the 2D Gaussian.
! Returns argument of the exponential and the norm: G = exp(arg) / norm.
!=================================================================================================
subroutine clustering_calculate_Gaussian_2D(x, y, mu1, mu2, sigma11, sigma22, sigma12, arg, norm)
  real(kind=CUSTOM_REAL), intent(in) :: x, y, mu1, mu2, sigma11, sigma22, sigma12
  real(kind=CUSTOM_REAL), intent(out) :: arg, norm

  ! Got these expressions from Mathematica using
  ! PDF[MultinormalDistribution[{mu1, mu2}, {{sigma11^2, sigma12^2}, {sigma12^2, sigma22^2}}], {x, y}]
  ! and then FortranForm[%].

  arg = (- ((-mu2 + y)*(mu2*sigma11**2 - mu1*sigma12**2 + sigma12**2*x - sigma11**2*y))/( sigma12**4 - sigma11**2*sigma22**2) &
         - ((-mu1 + x)*(mu2*sigma12**2 - mu1*sigma22**2 + sigma22**2*x - sigma12**2*y))/(-sigma12**4 + sigma11**2*sigma22**2))/2.

  norm = (2.*PI*sqrt(-sigma12**4 + sigma11**2*sigma22**2))

end subroutine clustering_calculate_Gaussian_2D

!===========================================================================================
! Calculate the joint non-normalized Gaussian for a given parameter.
!===========================================================================================
function clustering_calculate_Gaussian(mu, sigma, weight, val) result(res)
  real(kind=CUSTOM_REAL), intent(in) :: mu(2), sigma(3), weight(2)
  real(kind=CUSTOM_REAL), intent(in) :: val(2)

  real(kind=CUSTOM_REAL) :: res

  real(kind=CUSTOM_REAL) :: arg, norm
  real(kind=CUSTOM_REAL) :: x, y
  real(kind=CUSTOM_REAL) :: mu1, mu2
  real(kind=CUSTOM_REAL) :: sigma11, sigma22, sigma12

  x = val(1)
  y = val(2)
  mu1 = mu(1)
  mu2 = mu(2)
  sigma11 = sigma(1)
  sigma22 = sigma(2)
  sigma12 = sigma(3)

  if (weight(1) /= 0.d0 .and. weight(2) /= 0.d0) then
  ! 2D Gaussian.
    call clustering_calculate_Gaussian_2D(x, y, mu1, mu2, sigma11, sigma22, sigma12, arg, norm)
  else
  ! 1D Gaussian.
    if (weight(2) == 0.d0) then
      call clustering_calculate_Gaussian_1D(val(1), mu(1), sigma(1), arg, norm)

    else if (weight(1) == 0.d0) then
      call clustering_calculate_Gaussian_1D(val(2), mu(2), sigma(2), arg, norm)
    endif
  endif

  ! Sanity check.
  if (norm == 0.d0) then
    print *, 'Error: Zero norm in clustering_calculate_Gaussian! Exiting.'
    stop
  endif

  if (arg < - 100.d0) then
    ! To avoid SIGFPE Floating-point exception.
    ! Chosen some small value.
    res = exp(- 100.d0)
  else
    res = exp(arg) / norm
  endif

end function clustering_calculate_Gaussian

!===========================================================================================================
! Calculate the Gaussian mixture for a given model value.
! Return the Gaussian and its derivative.
!===========================================================================================================
subroutine clustering_calculate_Gaussian_mixture(this, val, gauss, deriv, gauss_loc, cluster_weight)
  class(t_clustering), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: val(2)
  real(kind=CUSTOM_REAL), intent(in) :: cluster_weight(:)

  real(kind=CUSTOM_REAL), intent(out) :: gauss
  real(kind=CUSTOM_REAL), intent(out) :: deriv(2)
  real(kind=CUSTOM_REAL), intent(out) :: gauss_loc(this%nclusters)

  real(kind=CUSTOM_REAL) :: coef(2)
  integer :: i, j

  real(kind=CUSTOM_REAL) :: x, y, mu1, mu2, sigma11, sigma22, sigma12

  gauss = 0.d0
  deriv = 0.d0

  do i = 1, this%nclusters
    gauss_loc(i) = cluster_weight(i) * &
                   this%calculate_Gaussian(this%mixture_mu(:, i), this%mixture_sigma(:, i), this%weight_loc, val)

    gauss = gauss + gauss_loc(i)

    ! ------ Calculate derivatives -------------------------------

    ! To calc. derivative coefficients run in Mathematica:
    ! ClearAll["Global`*"]
    ! f = PDF[MultinormalDistribution[{mu1, mu2}, {{sigma11^2, sigma12^2}, {sigma12^2, sigma22^2}}], {x, y}]
    ! FullSimplify[D[f, x]/f]
    ! FortranForm[%]

    x = val(1)
    y = val(2)
    mu1 = this%mixture_mu(1, i)
    mu2 = this%mixture_mu(2, i)
    sigma11 = this%mixture_sigma(1, i)
    sigma22 = this%mixture_sigma(2, i)
    sigma12 = this%mixture_sigma(3, i)

    ! d/dx
    coef(1) = (sigma22**2*(-mu1 + x) + sigma12**2*(mu2 - y))/(sigma12**4 - sigma11**2*sigma22**2)

    ! d/dy
    coef(2) = (sigma12**2*(mu1 - x) + sigma11**2*(-mu2 + y))/(sigma12**4 - sigma11**2*sigma22**2)

    do j = 1, 2
      deriv(j) = deriv(j) + coef(j) * gauss_loc(i)
    enddo

  enddo

end subroutine clustering_calculate_Gaussian_mixture

!==============================================================================================
! Calculate the maximum of the Gaussian mixture.
!==============================================================================================
function clustering_calculate_Gaussian_mixture_max(this, cluster_weight) result(res)
  class(t_clustering), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: cluster_weight(:)

  real(kind=CUSTOM_REAL) :: res

  real(kind=CUSTOM_REAL) :: gauss, deriv(2)
  real(kind=CUSTOM_REAL) :: gauss_loc(this%nclusters)
  integer :: i

  real(kind=CUSTOM_REAL) :: gauss_max

  gauss_max = 0.d0

  ! NOTE: assuming the maximum is located at one of the cluster center.
  do i = 1, this%nclusters
    ! Calculate the value of Gaussian mixture at the center (mean of Gaussian) of cluster i.
    call this%calculate_Gaussian_mixture(this%mixture_mu(:, i), gauss, deriv, gauss_loc, cluster_weight)

    if (gauss > gauss_max) gauss_max = gauss
  enddo

  res = gauss_max

end function clustering_calculate_Gaussian_mixture_max

end module clustering
