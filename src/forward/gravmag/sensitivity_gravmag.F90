
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

  implicit none

  private

  type, public :: t_sensitivity_gravmag
    private

  contains
    private

    procedure, public, nopass :: calculate_sensitivity
    procedure, public, nopass :: calc_data_directly

    procedure, public, nopass :: compress_matrix_line
    procedure, public, nopass :: compress_matrix_line_wavelet

  end type t_sensitivity_gravmag

contains

!=============================================================================================
! Calculates the sensitivity kernel (matrix).
!=============================================================================================
subroutine calculate_sensitivity(par, grid, data, sensit_matrix, myrank)
  class(t_parameters_base), intent(in) :: par
  type(t_grid), intent(in) :: grid
  type(t_data), intent(in) :: data
  integer, intent(in) :: myrank

  ! Sensitivity matrix.
  type(t_sparse_matrix), intent(inout) :: sensit_matrix

  type(t_magnetic_field) :: mag_field
  integer :: i, p, ierr
  real(kind=CUSTOM_REAL) :: comp_rate, comp_rate_max, comp_rate_min

  ! Sensitivity matrix row.
  real(kind=CUSTOM_REAL), allocatable :: sensit_line(:)
  real(kind=CUSTOM_REAL), allocatable :: sensit_line2(:)
  real(kind=CUSTOM_REAL), allocatable :: sensit_line3(:)

  allocate(sensit_line(par%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_sensitivity!", myrank, ierr)

  allocate(sensit_line2(par%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_sensitivity!", myrank, ierr)

  allocate(sensit_line3(par%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_sensitivity!", myrank, ierr)

  comp_rate_max = 0.d0
  comp_rate_min = 1.d0

  select type(par)
  class is (t_parameters_mag)

    if (myrank == 0) print *, 'Calculating MAGNETIC sensitivity kernel...'

    ! Precompute common parameters.
    call mag_field%initialize(par%mi, par%md, par%fi, par%fd, par%theta, par%intensity)

    ! Calculate sensitivity and adding to sparse matrix.
    ! Loop on all the data lines.
    do i = 1, par%ndata
      call mag_field%magprism(par%nelements, i, grid, data%X, data%Y, data%Z, sensit_line)

      call compress_matrix_line(par%nelements, sensit_line, data, i, grid, par%distance_threshold, comp_rate)

      if (comp_rate < comp_rate_min) comp_rate_min = comp_rate
      if (comp_rate > comp_rate_max) comp_rate_max = comp_rate

      if (i == 1 .and. myrank == 0) print *, 'Compression rate (for the 1st matrix line) =', comp_rate

      call sensit_matrix%new_row(myrank)

      do p = 1, par%nelements
        call sensit_matrix%add(sensit_line(p), p, myrank)
      enddo
    enddo

    call sensit_matrix%finalize(par%nelements, myrank)

  !----------- Gravity -------------------------------------------------------------------------
  class is (t_parameters_grav)

    if (myrank == 0) print *, 'Calculating GRAVITY sensitivity kernel...'

    ! Sanity check.
    if (mod(par%ndata, par%ncomponents) /= 0) &
      call exit_MPI("Number of data does not match the number of components!", myrank, 0)

    if (par%ncomponents == 1) then

      ! Loop on all the data lines.
      do i = 1, par%ndata
        call graviprism_full(par, grid, data%X(i), data%Y(i), data%Z(i), &
                             sensit_line, sensit_line2, sensit_line3, myrank)

        if (par%compression_type == 1) then
        ! Distance cut-off compression.
            call compress_matrix_line(par%nelements, sensit_line3, data, i, grid, par%distance_threshold, comp_rate)
        else if (par%compression_type == 2) then
        ! Wavelet compression.
            call compress_matrix_line_wavelet(par%nx, par%ny, par%nz, sensit_line3, par%wavelet_threshold, comp_rate)
        endif

        if (comp_rate < comp_rate_min) comp_rate_min = comp_rate
        if (comp_rate > comp_rate_max) comp_rate_max = comp_rate

        if (i == 1 .and. myrank == 0) print *, 'Compression rate (for the 1st matrix line) =', comp_rate

        call sensit_matrix%new_row(myrank)

        do p = 1, par%nelements
          ! Adding the Z-component only.
          call sensit_matrix%add(sensit_line3(p), p, myrank)
        enddo
      enddo

      call sensit_matrix%finalize(par%nelements, myrank)

    else if (par%ncomponents == 3) then

      if (myrank == 0) print *, 'Error: Not supported case!'
      stop
    endif

  end select

  if (myrank == 0) print *, 'Compression rate (min/max) = ', comp_rate_min, comp_rate_max

  deallocate(sensit_line)
  deallocate(sensit_line2)
  deallocate(sensit_line3)

  if (myrank == 0) print *, 'Finished calculating the sensitivity kernel.'

end subroutine calculate_sensitivity

!======================================================================================
! Calculates data directly from prior model without storing the sensitivity matrix.
! This way we can calculate data for big models/data sets using much less memory.
!======================================================================================
subroutine calc_data_directly(par, model, data, myrank)
  class(t_parameters_base), intent(in) :: par
  type(t_model), intent(in) :: model
  type(t_data), intent(inout) :: data
  integer, intent(in) :: myrank

  type(t_magnetic_field) :: mag_field
  integer :: i, ierr
  real(kind=CUSTOM_REAL), allocatable :: line(:)

  if (myrank == 0) print *, 'Error: Not supported case!'
  stop

  ! Allocate memory for one sensitivity matrix row.
  allocate(line(par%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in write_data_directly!", myrank, ierr)

  select type(par)
  class is (t_parameters_mag)

    if (myrank == 0) print *, 'Computing MAGNETIC data directly...'

    ! Precompute common parameters.
    call mag_field%initialize(par%mi, par%md, par%fi, par%fd, par%theta, par%intensity)

    do i = 1, par%ndata

      ! Calculate one sensitivity matrix line.
      call mag_field%magprism(par%nelements, 1, model%grid, data%X(i:i), data%Y(i:i), data%Z(i:i), line)

      ! Calculate one data using sensitivity line and prior model (m) as d = line * m.
      !call model%calculate_data(1, line, data%val_calc(i:i), myrank)

    enddo

  class is (t_parameters_grav)
    ! TODO.

  end select

end subroutine calc_data_directly

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

end module sensitivity_gravmag
