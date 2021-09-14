
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
! A class to calculate sensitivity values for gravity or magnetic field.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!========================================================================================
module sensitivity_gravmag

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use parameters_mag
  use parameters_grav
  use magnetic_field
  use gravity_field
  use grid
  use model
  use data_gravmag
  use sparse_matrix
  use vector
  use wavelet_transform
  use parallel_tools

  implicit none

  private

  type, public :: t_sensitivity_gravmag
    private

  contains
    private

    procedure, public, nopass :: calculate_sensitivity

    procedure, public, nopass :: compress_matrix_line
    procedure, public, nopass :: compress_matrix_line_wavelet

    procedure, private, nopass :: apply_column_weight

  end type t_sensitivity_gravmag

contains

!=============================================================================================
! Calculates the sensitivity kernel (matrix).
!=============================================================================================
subroutine calculate_sensitivity(par, grid, data, column_weight, sensit_matrix, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  type(t_grid), intent(in) :: grid
  type(t_data), intent(in) :: data
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  integer, intent(in) :: myrank, nbproc

  ! Sensitivity matrix.
  type(t_sparse_matrix), intent(inout) :: sensit_matrix

  type(t_magnetic_field) :: mag_field
  integer :: i, p, ierr
  real(kind=CUSTOM_REAL) :: comp_rate, comp_rate_min, comp_rate_max, comp_rate_tot
  integer :: nsmaller, nelements_total
  type(t_parallel_tools) :: pt
  integer :: problem_type
  integer :: nnz_line, nnz_total

  ! Sensitivity matrix row.
  real(kind=CUSTOM_REAL), allocatable :: sensit_line(:)
  real(kind=CUSTOM_REAL), allocatable :: sensit_line2(:)
  real(kind=CUSTOM_REAL), allocatable :: sensit_line3(:)
  real(kind=CUSTOM_REAL), allocatable :: sensit_line_full(:)

  allocate(sensit_line(par%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_sensitivity!", myrank, ierr)

  allocate(sensit_line2(par%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_sensitivity!", myrank, ierr)

  allocate(sensit_line3(par%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_sensitivity!", myrank, ierr)

  if (par%compression_type == 2) then
    if (nbproc > 1) then
      ! Number of parameters on ranks smaller than current one.
      nsmaller = pt%get_nsmaller(par%nelements, myrank, nbproc)

      ! Total number of elements.
      nelements_total = par%nx * par%ny * par%nz

      allocate(sensit_line_full(nelements_total), source=0._CUSTOM_REAL, stat=ierr)
      if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_sensitivity!", myrank, ierr)
    endif
  endif

  select type(par)
  class is (t_parameters_grav)
    if (myrank == 0) print *, 'Calculating GRAVITY sensitivity kernel...'
    problem_type = 1

  class is (t_parameters_mag)
    if (myrank == 0) print *, 'Calculating MAGNETIC sensitivity kernel...'
    problem_type = 2

    ! Precompute common magnetic parameters.
    call mag_field%initialize(par%mi, par%md, par%fi, par%fd, par%theta, par%intensity)
  end select

  !--------------------------------------------------------------------------------------------
  ! Calculating sensitivity and adding to the sparse matrix.
  ! Loop on all the data lines.
  do i = 1, par%ndata
    if (problem_type == 1) then
    ! Gravity problem.
      call graviprism_full(par%nelements, par%ncomponents, grid, data%X(i), data%Y(i), data%Z(i), &
                           sensit_line3, sensit_line2, sensit_line, myrank)
    else if (problem_type == 2) then
    ! Magnetic problem.
      call mag_field%magprism(par%nelements, i, grid, data%X, data%Y, data%Z, sensit_line)
    endif

    ! Applying the depth weight.
    call apply_column_weight(par%nelements, sensit_line, column_weight)

    if (par%compression_type == 1) then
    ! Distance cut-off compression.
      call compress_matrix_line(par%nelements, sensit_line, data, i, grid, par%distance_threshold, comp_rate)

    else if (par%compression_type == 2) then
    ! Wavelet compression.

      if (nbproc > 1) then
      ! Parallel wavelet copression.
        call pt%get_full_array(sensit_line, par%nelements, sensit_line_full, .true., myrank, nbproc)
        call compress_matrix_line_wavelet(par%nx, par%ny, par%nz, sensit_line_full, par%wavelet_threshold, comp_rate)

        ! Extract the local sensitivity part.
        sensit_line = sensit_line_full(nsmaller + 1 : nsmaller + par%nelements)
      else
      ! Serial.
        call compress_matrix_line_wavelet(par%nx, par%ny, par%nz, sensit_line, par%wavelet_threshold, comp_rate)
      endif

      ! Check if we have enough space in the matrix for new elemements, or we need to adjust the compression rate.
      nnz_line = count(sensit_line /= 0.d0)
      if (sensit_matrix%get_number_elements() + nnz_line > sensit_matrix%get_nnz()) then
        call exit_MPI("The matrix size is too small, adjust the compression rate!", myrank, ierr)
      endif
    endif

    call sensit_matrix%new_row(myrank)

    do p = 1, par%nelements
      ! Adding the Z-component only.
      call sensit_matrix%add(sensit_line(p), p, myrank)
    enddo
  enddo

  call sensit_matrix%finalize(par%nelements, myrank)

  ! Calculate the matrix compression rate.
  comp_rate = dble(sensit_matrix%get_number_elements()) / dble(par%nelements) / dble(par%ndata)
  if (nbproc > 1) then
    call mpi_allreduce(comp_rate, comp_rate_min, 1, CUSTOM_MPI_TYPE, MPI_MIN, MPI_COMM_WORLD, ierr)
    call mpi_allreduce(comp_rate, comp_rate_max, 1, CUSTOM_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD, ierr)

    call mpi_allreduce(sensit_matrix%get_number_elements(), nnz_total, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    comp_rate_tot = dble(nnz_total) / dble(nelements_total) / dble(par%ndata)

    if (myrank == 0) print *, 'Compression rate min = ', comp_rate_min
    if (myrank == 0) print *, 'Compression rate max = ', comp_rate_max
    if (myrank == 0) print *, 'Compression rate tot = ', comp_rate_tot
  else
    if (myrank == 0) print *, 'Compression rate = ', comp_rate
  endif

  deallocate(sensit_line)
  deallocate(sensit_line2)
  deallocate(sensit_line3)
  if (allocated(sensit_line_full)) deallocate(sensit_line_full)

  if (myrank == 0) print *, 'Finished calculating the sensitivity kernel.'

end subroutine calculate_sensitivity

!==========================================================================================================
! Compresses matrix line - distance based cutoff.
!==========================================================================================================
subroutine compress_matrix_line(nelements, line, data, idata, grid, distance_threshold, comp_rate)
  integer, intent(in) :: nelements, idata
  type(t_data), intent(in) :: data
  type(t_grid), intent(in) :: grid
  real(kind=CUSTOM_REAL), intent(in) :: distance_threshold

  real(kind=CUSTOM_REAL), intent(inout) :: line(:)
  real(kind=CUSTOM_REAL), intent(out) :: comp_rate

  integer :: i, num_zeros
  type(t_vector) :: dist_vec
  real(kind=CUSTOM_REAL) :: dist

  num_zeros = 0

  do i = 1, nelements
    ! Sphere based cutoff.
!    dist_vec = t_vector(data%X(idata) - grid%get_X_cell_center(i), &
!                        data%Y(idata) - grid%get_Y_cell_center(i), &
!                        data%Z(idata) - grid%get_Z_cell_center(i))

    ! Cylinder based cutoff.
    dist_vec = t_vector(data%X(idata) - grid%get_X_cell_center(i), &
                        data%Y(idata) - grid%get_Y_cell_center(i), &
                        0.d0)

    ! The distance between the source and cell center.
    dist = dist_vec%get_norm()

    if (dist > distance_threshold) then
      line(i) = 0.d0
      num_zeros = num_zeros + 1
    endif
  enddo

  comp_rate = dble(nelements - num_zeros) / dble(nelements)

end subroutine compress_matrix_line

!==========================================================================================================
! Compresses matrix line - wavelet compression
!==========================================================================================================
subroutine compress_matrix_line_wavelet(nx, ny, nz, line, threshold, comp_rate)
  integer, intent(in) :: nx, ny, nz
  real(kind=CUSTOM_REAL), intent(in) :: threshold

  real(kind=CUSTOM_REAL), intent(inout) :: line(:)
  real(kind=CUSTOM_REAL), intent(out) :: comp_rate

  integer :: nelements, i, num_zeros

  nelements = nx * ny * nz
  num_zeros = 0

  call Haar3D(line, nx, ny, nz)

  do i = 1, nelements
    if (abs(line(i)) < threshold) then
      line(i) = 0.d0
      num_zeros = num_zeros + 1
    endif
  enddo

  comp_rate = dble(nelements - num_zeros) / dble(nelements)

end subroutine compress_matrix_line_wavelet

!==========================================================================================================
! Applying the column weight to sensitivity line.
!==========================================================================================================
subroutine apply_column_weight(nelements, sensit_line, column_weight)
  integer, intent(in) :: nelements
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)

  real(kind=CUSTOM_REAL), intent(inout) :: sensit_line(:)

  integer :: i

  do i = 1, nelements
    sensit_line(i) = sensit_line(i) * column_weight(i)
  enddo

end subroutine apply_column_weight

end module sensitivity_gravmag
