
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
  use string, only: str
  use sort

  implicit none

  private

  public :: calculate_and_write_sensit
  public :: read_sensitivity_kernel
  public :: read_sensitivity_metadata
  public :: calculate_new_partitioning
  public :: write_depth_weight
  public :: read_depth_weight
  public :: partition_sensitivity_kernel

  private :: apply_column_weight
  private :: get_load_balancing_nelements
  private :: get_nel_compressed
  private :: read_depth_weight_file
  private :: read_sensit_nnz

  character(len=4) :: SUFFIX(2) = ["grav", "magn"]

contains

!===============================================================================================================
! Calculates the number of elements in the compressed sensitivity line.
!===============================================================================================================
pure function get_nel_compressed(par) result(nel_compressed)
  class(t_parameters_base), intent(in) :: par
  integer :: nel_compressed
  integer :: nelements_total

  nelements_total = par%nx * par%ny * par%nz

  if (par%compression_type > 0) then
    nel_compressed = int(par%compression_rate * nelements_total)
  else
    nel_compressed = nelements_total
  endif

end function get_nel_compressed

!===============================================================================================================
! Calculates the sensitivity kernel (parallelized by data) and writes it to files.
!===============================================================================================================
subroutine calculate_and_write_sensit(par, grid_full, data, column_weight, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  type(t_grid), intent(in) :: grid_full
  type(t_data), intent(in) :: data
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(par%nelements)
  integer, intent(in) :: myrank, nbproc

  type(t_magnetic_field) :: mag_field
  integer :: i, k, p, d, nel, ierr
  integer :: nelements_total, nel_compressed
  integer(kind=8) :: nnz_data
  integer :: problem_type
  integer :: idata, ndata_loc, ndata_smaller
  character(len=256) :: filename, filename_full
  real(kind=CUSTOM_REAL) :: comp_rate, comp_error
  integer(kind=8) :: nnz_total
  real(kind=CUSTOM_REAL) :: threshold
  character(len=256) :: msg

  ! Sensitivity matrix row.
  real(kind=CUSTOM_REAL), allocatable :: sensit_line_full(:, :, :)
  real(kind=CUSTOM_REAL), allocatable :: sensit_line_sorted(:)

  ! Arrays for storing the compressed sensitivity line.
  integer, allocatable :: sensit_columns(:)
  real(kind=MATRIX_PRECISION), allocatable :: sensit_compressed(:)

  ! To calculate the partitioning for balanced memory loading among CPUs.
  integer, allocatable :: sensit_nnz(:)
  integer(kind=8) :: sensit_nnz_sum

  real(kind=CUSTOM_REAL) :: cost_full, cost_compressed
  real(kind=CUSTOM_REAL) :: error_r_i, error_r_sum, error_r_sum_loc
  real(kind=CUSTOM_REAL) :: relative_threshold, relative_threshold_sum, relative_threshold_sum_loc

  ! The full column weight.
  real(kind=CUSTOM_REAL), allocatable :: column_weight_full(:)

  ! Sanity check.
  if (par%compression_rate < 0 .or. par%compression_rate > 1) then
    call exit_MPI("Wrong compression rate! It must be between 0 and 1.", myrank, 0)
  endif

  ! Define the problem type.
  select type(par)
  class is (t_parameters_grav)
    if (myrank == 0) print *, 'Calculating GRAVITY sensitivity kernel...'
    problem_type = 1

  class is (t_parameters_mag)
    if (myrank == 0) print *, 'Calculating MAGNETIC sensitivity kernel...'
    problem_type = 2

    ! Precompute common magnetic parameters.
    call mag_field%initialize(par%mi, par%md, par%theta, par%intensity)
  end select

  !---------------------------------------------------------------------------------------------
  ! Define the output file.
  !---------------------------------------------------------------------------------------------
  call execute_command_line('mkdir -p "'//trim(path_output)//'/SENSIT"')

  filename = "sensit_"//SUFFIX(problem_type)//"_"//trim(str(nbproc))//"_"//trim(str(myrank))
  filename_full = trim(path_output)//"/SENSIT/"//filename

  if (myrank == 0) print *, 'Writing the sensitivity to file ', trim(filename_full)

  open(77, file=trim(filename_full), status='replace', access='stream', form='unformatted', action='write', &
       iostat=ierr, iomsg=msg)

  if (ierr /= 0) call exit_MPI("Error in creating the sensitivity file! path=" &
                               //trim(filename_full)//", iomsg="//msg, myrank, ierr)

  !---------------------------------------------------------------------------------------------
  ! Allocate memory.
  !---------------------------------------------------------------------------------------------
  nelements_total = par%nx * par%ny * par%nz

  nel_compressed = get_nel_compressed(par)
  if (myrank == 0) print *, 'nel_compressed =', nel_compressed

  allocate(sensit_line_full(nelements_total, par%nmodel_components, par%ndata_components), source=0._CUSTOM_REAL, stat=ierr)
  allocate(sensit_line_sorted(nelements_total), source=0._CUSTOM_REAL, stat=ierr)
  allocate(column_weight_full(nelements_total), source=0._CUSTOM_REAL, stat=ierr)

  allocate(sensit_columns(nel_compressed), source=0, stat=ierr)
  allocate(sensit_compressed(nel_compressed), source=0._MATRIX_PRECISION, stat=ierr)

  allocate(sensit_nnz(nelements_total), source=0, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_and_write_sensit!", myrank, ierr)

  !---------------------------------------------------------------------------------------------
  ! Calculate sensitivity lines.
  !---------------------------------------------------------------------------------------------
  call get_full_array(column_weight, par%nelements, column_weight_full, .true., myrank, nbproc)

  ndata_loc = calculate_nelements_at_cpu(par%ndata, myrank, nbproc)
  ndata_smaller = get_nsmaller(ndata_loc, myrank, nbproc)

  ! File header.
  write(77) ndata_loc, par%ndata, nelements_total, myrank, nbproc

  nnz_data = 0
  error_r_sum_loc = 0.d0
  relative_threshold_sum_loc = 0.d0

  ! Loop over the local data lines.
  do i = 1, ndata_loc
    ! Global data index.
    idata = ndata_smaller + i

    if (problem_type == 1) then
      ! Calcualte gravity kernel.
      if (par%data_type == 1) then
        call graviprism_z(nelements_total, grid_full, data%X(idata), data%Y(idata), data%Z(idata), &
                          sensit_line_full, myrank)
      else if (par%data_type == 2) then
        ! Gradiometry.
        if (par%ndata_components == 1) then
          ! Only Gzz component.
          call gradiprism_zz(nelements_total, grid_full, data%X(idata), data%Y(idata), data%Z(idata), &
                             sensit_line_full)

        else if (par%ndata_components == 6) then
          ! Full tensor.
          call gradiprism_full(nelements_total, grid_full, data%X(idata), data%Y(idata), data%Z(idata), &
                               sensit_line_full(:, 1, 1), sensit_line_full(:, 1, 2), sensit_line_full(:, 1, 3), &
                               sensit_line_full(:, 1, 4), sensit_line_full(:, 1, 5), sensit_line_full(:, 1, 6), myrank)
        else
          call exit_MPI("Wrong number of gravity gradiometry data components!", myrank, 0)
        endif
      endif

    else if (problem_type == 2) then
      ! Calculate magnetic kernel.
      call mag_field%magprism(nelements_total, par%nmodel_components, par%ndata_components, grid_full, &
                              data%X(idata), data%Y(idata), data%Z(idata), sensit_line_full)
    endif

    ! Loop over the data components.
    do d = 1, par%ndata_components

      ! Loop over the model components.
      do k = 1, par%nmodel_components

        ! Applying the depth weight.
        call apply_column_weight(nelements_total, sensit_line_full(:, k, d), column_weight_full)

        if (par%compression_type > 0) then
        ! Wavelet compression.

          ! The uncompressed line cost.
          cost_full = sum(sensit_line_full(:, k, d)**2)

          ! Apply the wavelet transform.
          call forward_wavelet(sensit_line_full(:, k, d), par%nx, par%ny, par%nz, par%compression_type)

          ! Perform sorting (to determine the wavelet threshold corresponding to the desired compression rate).
          sensit_line_sorted = abs(sensit_line_full(:, k, d))
          call quicksort(sensit_line_sorted, 1, nelements_total)

          ! Calculate the wavelet threshold corresponding to the desired compression rate.
          if (nel_compressed >= nelements_total) then
            ! Taking all elements.
            threshold = -1.d0
          else
            p = nelements_total - nel_compressed
            threshold = abs(sensit_line_sorted(p))
          endif

          if (threshold < 1.d-30) then
            ! Keep small threshold to avoid extremely small values, because when MATRIX_PRECISION is 4 bytes real
            ! they lead to SIGFPE: Floating-point exception - erroneous arithmetic operation.
            threshold = 1.d-30
          endif

          nel = 0
          cost_compressed = 0.d0
          do p = 1, nelements_total
            if (abs(sensit_line_full(p, k, d)) > threshold) then
              ! Store sensitivity elements greater than the wavelet threshold.
              nel = nel + 1
              sensit_columns(nel) = p
              sensit_compressed(nel) = real(sensit_line_full(p, k, d), MATRIX_PRECISION)

              sensit_nnz(p) = sensit_nnz(p) + 1

              cost_compressed = cost_compressed + sensit_line_full(p, k, d)**2
            endif
          enddo

          ! Sanity check.
          if (nel > nel_compressed) then
            call exit_MPI("Wrong number of elements in calculate_and_write_sensit!", myrank, 0)
          endif

          !--------------------------------------------------------------------------------------
          ! Calculate compression statistics.
          !--------------------------------------------------------------------------------------
          if (nel_compressed == nelements_total) then
            ! Assign the error directly as the formula is numerically unstable for this case.
            error_r_i = 0.d0
          else
            ! Compression error for this row. See Eq.(19) in Li & Oldenburg, GJI (2003) 152, 251–265.
            error_r_i = sqrt(1.d0 - cost_compressed / cost_full)
          endif

          ! Relative threshold.
          relative_threshold = threshold / abs(sensit_line_sorted(nelements_total))

          error_r_sum_loc = error_r_sum_loc + error_r_i
          relative_threshold_sum_loc = relative_threshold_sum_loc + relative_threshold

        else
        ! No compression.
          nel = nelements_total
          sensit_compressed = real(sensit_line_full(:, k, d), MATRIX_PRECISION)
          do p = 1, nelements_total
            sensit_columns(p) = p
            sensit_nnz(p) = sensit_nnz(p) + 1
          enddo
        endif

        ! The sensitivity kernel size (when parallelized by data).
        nnz_data = nnz_data + nel

        ! Sanity check.
        if (nnz_data < 0) then
          call exit_MPI("Integer overflow in nnz_data!", myrank, 0)
        endif

        ! Write the sensitivity line to file.
        write(77) idata, nel, k, d
        if (nel > 0) then
          write(77) sensit_columns(1:nel)
          write(77) sensit_compressed(1:nel)
        endif
      enddo ! nmodel_components loop
    enddo ! ndata_components loop

    ! Print the progress.
    if (myrank == 0 .and. mod(i, max(int(0.1d0 * ndata_loc), 1)) == 0) then
      print *, 'Percent completed: ', int(dble(i) / dble(ndata_loc) * 100.d0)
    endif

  enddo ! data loop

  close(77)

  call MPI_Allreduce(MPI_IN_PLACE, sensit_nnz, nelements_total, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  !---------------------------------------------------------------------------------------------
  ! Calculate the kernel compression rate.
  !---------------------------------------------------------------------------------------------
  call MPI_Allreduce(nnz_data, nnz_total, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
  comp_rate = dble(nnz_total) / dble(nelements_total) / dble(par%ndata) / dble(par%nmodel_components) / dble(par%ndata_components)

  if (myrank == 0) print *, 'nnz_total = ', nnz_total
  if (myrank == 0) print *, 'COMPRESSION RATE = ', comp_rate

  ! Sanity check.
  sensit_nnz_sum = 0
  do p = 1, nelements_total
    ! Calculate sum using loop because sum() returns the same type as an array hence leads to overflow.
    sensit_nnz_sum = sensit_nnz_sum + sensit_nnz(p)
  enddo
  if (sensit_nnz_sum /= nnz_total) then
    call exit_MPI("Wrong nnz_total in calculate_and_write_sensit!", myrank, 0)
  endif

  !---------------------------------------------------------------------------------------------
  ! Calculate the kernel compression error.
  !---------------------------------------------------------------------------------------------
  if (par%compression_type > 0) then
    call MPI_Allreduce(error_r_sum_loc, error_r_sum, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(relative_threshold_sum_loc, relative_threshold_sum, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate the average over all rows (compressions).
    comp_error = error_r_sum / dble(par%ndata * par%ndata_components * par%nmodel_components)
    relative_threshold = relative_threshold_sum / dble(par%ndata * par%ndata_components * par%nmodel_components)
  else
    comp_error = 0.d0
    relative_threshold = 0.d0
  endif

  if (myrank == 0) print *, 'COMPRESSION ERROR, r = ', comp_error
  if (myrank == 0) print *, 'Relative threshold, e = ', relative_threshold

  !---------------------------------------------------------------------------------------------
  ! Write the metadata file.
  !---------------------------------------------------------------------------------------------
  if (myrank == 0) then
    filename = "sensit_"//SUFFIX(problem_type)//"_meta.txt"
    filename_full = trim(path_output)//"/SENSIT/"//filename

    print *, 'Writing the sensitivity metadata to file ', trim(filename_full)

    open(77, file=trim(filename_full), form='formatted', status='replace', action='write')

    write(77, *) par%nx, par%ny, par%nz, par%ndata
    write(77, *) nbproc, MATRIX_PRECISION, par%depth_weighting_type
    write(77, *) par%compression_type, comp_error
    write(77, *) par%nmodel_components, par%ndata_components
    write(77, *) nnz_total

    close(77)
  endif

  !---------------------------------------------------------------------------------------------
  ! Write the sensit_nnz (to calculate load balancing for any number of CPUs).
  !---------------------------------------------------------------------------------------------
  if (myrank == 0) then
    filename = "sensit_"//SUFFIX(problem_type)//"_nnz"
    filename_full = trim(path_output)//"/SENSIT/"//filename

    print *, 'Writing the sensit_nnz to file ', trim(filename_full)

    open(77, file=trim(filename_full), form='unformatted', status='replace', action='write', access='stream')

    write(77) nelements_total
    write(77) sensit_nnz

    close(77)
  endif

  !---------------------------------------------------------------------------------------------
  deallocate(sensit_line_full)
  deallocate(sensit_line_sorted)
  deallocate(column_weight_full)
  deallocate(sensit_columns)
  deallocate(sensit_compressed)
  deallocate(sensit_nnz)

  if (myrank == 0) print *, 'Finished calculating the sensitivity kernel.'

end subroutine calculate_and_write_sensit

!=======================================================================================================================
! Write the depth weight to file.
!=======================================================================================================================
subroutine write_depth_weight(par, column_weight, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(par%nelements)
  integer, intent(in) :: myrank, nbproc

  character(len=256) :: filename, filename_full
  integer :: problem_type, nelements_total
  integer :: ierr
  ! The full column weight.
  real(kind=CUSTOM_REAL), allocatable :: column_weight_full(:)

  nelements_total = par%nx * par%ny * par%nz

  ! Allocate the full array on master rank only.
  if (myrank == 0) then
    allocate(column_weight_full(nelements_total), source=0._CUSTOM_REAL, stat=ierr)
  else
    allocate(column_weight_full(1), source=0._CUSTOM_REAL, stat=ierr)
  endif

  call get_full_array(column_weight, par%nelements, column_weight_full, .false., myrank, nbproc)

  if (myrank == 0) then
    ! Define the problem type.
    select type(par)
    class is (t_parameters_grav)
      problem_type = 1
    class is (t_parameters_mag)
      problem_type = 2
    end select

    ! Create the sensit folder.
    call execute_command_line('mkdir -p "'//trim(path_output)//'/SENSIT"')

    ! The output file name.
    filename = "sensit_"//SUFFIX(problem_type)//"_weight"
    filename_full = trim(path_output)//"/SENSIT/"//filename

    print *, 'Writing the depth weight to file ', trim(filename_full)

    open(77, file=trim(filename_full), form='unformatted', status='replace', action='write', access='stream')

    write(77) nelements_total
    write(77) column_weight_full

    close(77)
  endif

  deallocate(column_weight_full)
end subroutine write_depth_weight

!=======================================================================================================================
! Calculate new partitioning for the nnz and nelements at every CPU for the load balancing.
! The load balancing aims to have the same number of nnz at every CPU, which require using different nelements at CPUs.
!=======================================================================================================================
subroutine get_load_balancing_nelements(nelements_total, sensit_nnz, &
                                        nnz_at_cpu_new, nelements_at_cpu_new, myrank, nbproc)
  integer, intent(in) :: nelements_total
  integer, intent(in) :: sensit_nnz(nelements_total)
  integer, intent(in) :: myrank, nbproc

  integer, intent(out) :: nelements_at_cpu_new(nbproc)
  integer(kind=8), intent(out) :: nnz_at_cpu_new(nbproc)

  integer :: cpu, p
  integer :: nelements_new
  integer(kind=8) :: nnz_total, nnz_new, sum_sensit_nnz
  integer(kind=8) :: nnz_at_cpu_best(nbproc)

  ! Calculate nnz_total.
  nnz_total = 0
  do p = 1, nelements_total
    ! Calculate sum using loop because sum() returns the same type as an array hence leads to overflow.
    nnz_total = nnz_total + sensit_nnz(p)
  enddo

  nnz_at_cpu_best(:) = nnz_total / int(nbproc, 8)
  ! Last rank gets the remaining elements.
  nnz_at_cpu_best(nbproc) = nnz_at_cpu_best(nbproc) + mod(nnz_total, int(nbproc, 8))

  cpu = 1
  nnz_new = 0
  nelements_new = 0
  nnz_at_cpu_new = 0
  sum_sensit_nnz = 0

  do p = 1, nelements_total
    nnz_new = nnz_new + sensit_nnz(p)
    sum_sensit_nnz = sum_sensit_nnz + sensit_nnz(p)
    nelements_new = nelements_new + 1

    if ((sum_sensit_nnz >= sum(nnz_at_cpu_best(1:cpu)) .and. (cpu < nbproc)) .or. p == nelements_total) then
      nnz_at_cpu_new(cpu) = nnz_new
      nelements_at_cpu_new(cpu) = nelements_new

      nnz_new = 0
      nelements_new = 0
      cpu = cpu + 1
    endif
  enddo

  ! Sanity check: test that all cpus got elements allocated.
  if (cpu /= nbproc + 1) then
    call exit_MPI("Wrong cpu in get_load_balancing_nelements!", myrank, 0)
  endif

  if (sum(nnz_at_cpu_new) /= nnz_total) then
    call exit_MPI("Wrong nnz_at_cpu_new in get_load_balancing_nelements!", myrank, 0)
  endif
end subroutine get_load_balancing_nelements

!==============================================================================================
! Reads sensit_nnz from file.
!==============================================================================================
subroutine read_sensit_nnz(par, problem_type, nelements_total, sensit_nnz, myrank)
  class(t_parameters_base), intent(in) :: par
  integer, intent(in) :: problem_type, nelements_total, myrank
  integer, intent(out) :: sensit_nnz(nelements_total)

  character(len=256) :: filename, filename_full
  character(len=256) :: msg
  integer :: ierr
  integer :: nelements_total_read

  ! Form the file name.
  filename = "sensit_"//SUFFIX(problem_type)//"_nnz"

  if (par%sensit_read /= 0) then
    filename_full = trim(par%sensit_path)//filename
  else
    filename_full = trim(path_output)//"/SENSIT/"//filename
  endif

  if (myrank == 0) print *, 'Reading the sensit_nnz file ', trim(filename_full)

  ! Open the file.
  open(78, file=trim(filename_full), status='old', access='stream', form='unformatted', action='read', &
       iostat=ierr, iomsg=msg)

  if (ierr /= 0) call exit_MPI("Error in opening the sensitivity file! path=" &
                               //trim(filename_full)//", iomsg="//msg, myrank, ierr)

  ! Reading the file header.
  read(78) nelements_total_read

  ! Consistency check.
  if (nelements_total_read /= nelements_total) then
    call exit_MPI("Wrong file header in calculate_new_partitioning!", myrank, 0)
  endif

  read(78) sensit_nnz

  close(78)
end subroutine read_sensit_nnz

!==============================================================================================
! Calculate new partitioning for the load balancing using the sensit_nnz from file.
!==============================================================================================
subroutine calculate_new_partitioning(par, nnz, nelements_at_cpu, problem_type, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  integer, intent(in) :: problem_type
  integer, intent(in) :: myrank, nbproc

  integer(kind=8), intent(out) :: nnz
  integer, intent(out) :: nelements_at_cpu(nbproc)

  integer(kind=8) :: nnz_at_cpu(nbproc)

  integer :: ierr
  integer :: nelements_total

  integer, allocatable :: sensit_nnz(:)
  integer, allocatable :: sensit_nnz2(:)

  if (myrank == 0) then

    nelements_total = par%nx * par%ny * par%nz

    ! Memory allocation.
    allocate(sensit_nnz(nelements_total), source=0, stat=ierr)

    if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_new_partitioning!", myrank, ierr)

    if (problem_type == 3) then
    ! Joint inversion. Calculate the total sennsit_nnz for both problems.

      ! Memory allocation.
      allocate(sensit_nnz2(nelements_total), source=0, stat=ierr)

      if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_new_partitioning!", myrank, ierr)

      ! Read sensit_nnz from file for problem 1.
      call read_sensit_nnz(par, 1, nelements_total, sensit_nnz2, myrank)

      sensit_nnz = sensit_nnz + sensit_nnz2

      ! Read sensit_nnz from file for problem 2.
      call read_sensit_nnz(par, 2, nelements_total, sensit_nnz2, myrank)

      ! The total sensit_nnz for both problems
      sensit_nnz = sensit_nnz + sensit_nnz2

      deallocate(sensit_nnz2)

    else
      ! Read sensit_nnz from file.
      call read_sensit_nnz(par, problem_type, nelements_total, sensit_nnz, myrank)
    endif

    ! Calculate the load balancing.
    call get_load_balancing_nelements(nelements_total, sensit_nnz, nnz_at_cpu, nelements_at_cpu, myrank, nbproc)

    deallocate(sensit_nnz)

  endif ! myrank == 0

  call MPI_Bcast(nnz_at_cpu, nbproc, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(nelements_at_cpu, nbproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  ! Return the parameter for the current rank.
  nnz = nnz_at_cpu(myrank + 1)

  if (myrank == 0) then
    print *, "nelements_at_cpu =", nelements_at_cpu
    print *, "nnz_at_cpu =", nnz_at_cpu
  endif

end subroutine calculate_new_partitioning

!============================================================================================================================
! Reads the sensitivity kernel from files,
! and stores it in the sparse matrix parallelized by model.
!============================================================================================================================
subroutine read_sensitivity_kernel(par, sensit_matrix, column_weight, problem_weight, data_weight, problem_type, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  real(kind=CUSTOM_REAL), intent(in) :: problem_weight
  real(kind=CUSTOM_REAL), intent(in) :: data_weight(par%ndata_components, par%ndata)
  integer, intent(in) :: problem_type
  integer, intent(in) :: myrank, nbproc

  ! Sensitivity matrix.
  type(t_sparse_matrix), intent(inout) :: sensit_matrix
  ! Depth weight.
  real(kind=CUSTOM_REAL), intent(out) :: column_weight(par%nelements)

  ! Arrays for storing the compressed sensitivity line.
  integer, allocatable :: sensit_columns(:)
  real(kind=MATRIX_PRECISION), allocatable :: sensit_compressed(:)

  real(kind=CUSTOM_REAL) :: comp_rate
  integer(kind=8) :: nnz_total
  integer :: i, j, k, d, nsmaller, ierr
  integer :: rank, nelements_total, nel_compressed
  character(len=256) :: filename, filename_full
  character(len=256) :: msg

  integer :: ndata_loc, ndata_read, nelements_total_read, myrank_read, nbproc_read
  integer :: nbproc_sensit
  integer :: idata, nel, idata_glob, model_component, data_component
  integer(kind=8) :: nnz
  integer :: param_shift(2)
  integer :: istart, iend
  integer :: index_shift
  integer(kind=8) :: pos0
  integer :: ncol

  param_shift(1) = 0
  param_shift(2) = par%nelements * par%nmodel_components

  !---------------------------------------------------------------------------------------------
  ! Allocate memory.
  !---------------------------------------------------------------------------------------------
  nel_compressed = get_nel_compressed(par)

  allocate(sensit_columns(nel_compressed), source=0, stat=ierr)
  allocate(sensit_compressed(nel_compressed), source=0._MATRIX_PRECISION, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in read_sensitivity_kernel!", myrank, ierr)

  !---------------------------------------------------------------------------------------------
  ! Read sensitivity metadata and retrieve nbproc_sensit.
  !---------------------------------------------------------------------------------------------
  call read_sensitivity_metadata(par, nbproc_sensit, problem_type, myrank)

  if (myrank == 0) print *, "nbproc_sensit =", nbproc_sensit

  !---------------------------------------------------------------------------------------------
  ! Reading the sensitivity kernel files.
  !---------------------------------------------------------------------------------------------

  nelements_total = par%nx * par%ny * par%nz

  ! The number of elements on CPUs with rank smaller than myrank.
  nsmaller = get_nsmaller(par%nelements, myrank, nbproc)

  idata_glob = 0
  nnz = 0

  ! Loop over the MPI ranks (as the sensitivity kernel is stored parallelized by data in a separate file for each rank).
  do rank = 0, nbproc_sensit - 1
    ! Form the file name (containing the MPI rank).
    filename = "sensit_"//SUFFIX(problem_type)//"_"//trim(str(nbproc_sensit))//"_"//trim(str(rank))

    if (par%sensit_read /= 0) then
      filename_full = trim(par%sensit_path)//filename
    else
      filename_full = trim(path_output)//"/SENSIT/"//filename
    endif

    if (myrank == 0 .and. rank == 0) print *, 'Reading the sensitivity file (new) ', trim(filename_full)

    open(78, file=trim(filename_full), status='old', access='stream', form='unformatted', action='read', &
         iostat=ierr, iomsg=msg)

    if (ierr /= 0) call exit_MPI("Error in opening the sensitivity file! path=" &
                                 //trim(filename_full)//", iomsg="//msg, myrank, ierr)

    ! Reading the file header.
    read(78) ndata_loc, ndata_read, nelements_total_read, myrank_read, nbproc_read

    ! The current position of file pointer.
    pos0 = 1 + 5 * 4

    ! Consistency check.
    if (ndata_read /= par%ndata .or. nelements_total_read /= nelements_total &
        .or. myrank_read /= rank .or. nbproc_read /= nbproc_sensit) then
      call exit_MPI("Wrong file header in read_sensitivity_kernel!", myrank, 0)
    endif

    ! Loop over the local data chunk within a file.
    do i = 1, ndata_loc
      idata_glob = idata_glob + 1

      ! Loop over data components.
      do d = 1, par%ndata_components

        ! Loop over model components.
        do k = 1, par%nmodel_components

          ! Reading the data descriptor.
          read(78, pos=pos0) idata, nel, model_component, data_component
          pos0 = pos0 + 4 * 4

          ! Make sure the data is stored in the correct order.
          if (idata /= idata_glob) then
            call exit_MPI("Wrong data index in read_sensitivity_kernel!", myrank, 0)
          endif

          ! Sanity check for the nel.
          if (nel > nel_compressed) then
            call exit_MPI("Wrong number of elements in read_sensitivity_kernel!", myrank, 0)
          endif

          ! Sanity check for the model component index.
          if (model_component /= k) then
            call exit_MPI("Wrong model component index in read_sensitivity_kernel!", myrank, 0)
          endif

          ! Sanity check for the data component index.
          if (data_component /= d) then
            call exit_MPI("Wrong data component index in read_sensitivity_kernel!", myrank, 0)
          endif

          if (nel > 0) then
            ! Reading the data.
            read(78) sensit_columns(1:nel)
            pos0 = pos0 + nel * 4
          endif

          ! Determine istart index (for the matrix elements corresponding to the current cpu).
          istart = 0
          do j = 1, nel
            if (sensit_columns(j) > nsmaller) then
              istart = j
              exit
            endif
          enddo

          if (istart > 0) then
          ! Found matrix elements.

            ! Determine iend index (for the matrix elements corresponding to the current cpu).
            iend = nel
            do j = istart, nel
              if (sensit_columns(j) > nsmaller + par%nelements) then
                iend = j - 1
                exit
              endif
            enddo

            ! The number of non-zero matrix elements on this rank.
            ncol = iend - istart + 1

            ! Read only the sensitivity matrix elements corresponding to the current rank.
            read(78, pos=pos0 + (istart - 1) * MATRIX_PRECISION) sensit_compressed(1:ncol)

            ! Column index shift.
            index_shift = param_shift(problem_type) + (k - 1) * par%nelements - nsmaller

            ! The column index in a big (joint) parallel sparse matrix.
            sensit_columns(istart:iend) = sensit_columns(istart:iend) + index_shift

            ! Apply the problem weight.
            sensit_compressed(1:ncol) = sensit_compressed(1:ncol) * real(problem_weight, MATRIX_PRECISION)

            ! Apply the data error.
            sensit_compressed(1:ncol) = sensit_compressed(1:ncol) * real(data_weight(d, idata), MATRIX_PRECISION)

            ! Add a complete matrix row.
            call sensit_matrix%add_row(ncol, sensit_compressed(1:ncol), sensit_columns(istart:iend), myrank)

            ! Number of nonzero elements.
            nnz = nnz + ncol

          endif

          ! Shift the file pointer to the next matrix row.
          pos0 = pos0 + nel * MATRIX_PRECISION

        enddo ! model components loop

        ! Adding the matrix row.
        call sensit_matrix%new_row(myrank)

      enddo ! data components loop

    enddo ! data loop

    close(78)
  enddo

  deallocate(sensit_columns)
  deallocate(sensit_compressed)

  !---------------------------------------------------------------------------------------------
  ! Calculate the read kernel compression rate.
  !---------------------------------------------------------------------------------------------
  call MPI_Allreduce(nnz, nnz_total, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
  comp_rate = dble(nnz_total) / dble(nelements_total) / dble(par%ndata) / dble(par%nmodel_components) / dble(par%ndata_components)

  if (myrank == 0) print *, 'nnz_total (of the read kernel)  = ', nnz_total
  if (myrank == 0) print *, 'COMPRESSION RATE (of the read kernel)  = ', comp_rate

  !---------------------------------------------------------------------------------------------
  ! Read the depth weight
  !---------------------------------------------------------------------------------------------
  call read_depth_weight(par, column_weight, myrank, nbproc)

  if (myrank == 0) print *, 'Finished reading the sensitivity kernel.'

end subroutine read_sensitivity_kernel

!============================================================================================================================
! Reads the sensitivity column indexes and performs partitioning by columns needed
! for efficient reading of the sensitivity parallelised by model parameters.
!============================================================================================================================
subroutine partition_sensitivity_kernel(par, problem_type, myrank, nbproc, nelements_at_cpu)
  class(t_parameters_base), intent(in) :: par
  integer, intent(in) :: problem_type
  integer, intent(in) :: myrank, nbproc
  integer, intent(in) :: nelements_at_cpu(nbproc)

  ! Arrays for storing the compressed sensitivity line.
  integer, allocatable :: sensit_columns(:)

  integer :: i, k, d, ierr
  integer :: nelements_total, nel_compressed, ndata_smaller
  character(len=256) :: filename, filename_full
  character(len=256) :: filename_partit, filename_full_partit
  character(len=256) :: msg

  integer :: ndata_loc, ndata_read, nelements_total_read, myrank_read, nbproc_read
  integer :: idata, nel, idata_glob, model_component, data_component
  integer(kind=8) :: pos0
  integer :: cpu, p
  integer :: nel_at_cpu(nbproc)

  !---------------------------------------------------------------------------------------------
  ! Allocate memory.
  !---------------------------------------------------------------------------------------------
  nel_compressed = get_nel_compressed(par)

  allocate(sensit_columns(nel_compressed), source=0, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in partition_sensitivity_columns!", myrank, ierr)

  !---------------------------------------------------------------------------------------------
  ! Reading the sensitivity kernel files.
  !---------------------------------------------------------------------------------------------
  nelements_total = par%nx * par%ny * par%nz

  ! Form the file name (containing the MPI rank).
  filename = "sensit_"//SUFFIX(problem_type)//"_"//trim(str(nbproc))//"_"//trim(str(myrank))
  filename_partit = "partit_"//SUFFIX(problem_type)//"_"//trim(str(nbproc))//"_"//trim(str(myrank))

  if (par%sensit_read /= 0) then
    filename_full = trim(par%sensit_path)//filename
    filename_full_partit = trim(par%sensit_path)//filename_partit
  else
    filename_full = trim(path_output)//"/SENSIT/"//filename
    filename_full_partit = trim(path_output)//"/SENSIT/"//filename_partit
  endif

  if (myrank == 0) print *, 'Reading the sensitivity file ', trim(filename_full)

  open(78, file=trim(filename_full), status='old', access='stream', form='unformatted', action='read', &
       iostat=ierr, iomsg=msg)

  if (ierr /= 0) call exit_MPI("Error in opening the sensitivity file! path=" &
                               //trim(filename_full)//", iomsg="//msg, myrank, ierr)

  open(77, file=trim(filename_full_partit), form='unformatted', status='replace', action='write', access='stream')

  ! Reading the file header.
  read(78) ndata_loc, ndata_read, nelements_total_read, myrank_read, nbproc_read

  ! The current position of file pointer.
  pos0 = 1 + 5 * 4

  ! Consistency check.
  if (ndata_read /= par%ndata .or. nelements_total_read /= nelements_total &
      .or. myrank_read /= myrank .or. nbproc_read /= nbproc) then
    call exit_MPI("Wrong file header in partition_sensitivity_columns!", myrank, 0)
  endif

  ndata_smaller = get_nsmaller(ndata_loc, myrank, nbproc)

  ! Loop over the local data chunk within a file.
  do i = 1, ndata_loc
    idata_glob = ndata_smaller + i

    ! Loop over data components.
    do d = 1, par%ndata_components

      ! Loop over model components.
      do k = 1, par%nmodel_components

        ! Reading the data descriptor.
        read(78, pos=pos0) idata, nel, model_component, data_component
        pos0 = pos0 + 4 * 4

        ! Make sure the data is stored in the correct order.
        if (idata /= idata_glob) then
          call exit_MPI("Wrong data index in partition_sensitivity_columns!", myrank, 0)
        endif

        ! Sanity check for the nel.
        if (nel > nel_compressed) then
          call exit_MPI("Wrong number of elements in partition_sensitivity_columns!", myrank, 0)
        endif

        ! Sanity check for the model component index.
        if (model_component /= k) then
          call exit_MPI("Wrong model component index in partition_sensitivity_columns!", myrank, 0)
        endif

        ! Sanity check for the data component index.
        if (data_component /= d) then
          call exit_MPI("Wrong data component index in partition_sensitivity_columns!", myrank, 0)
        endif

        if (nel > 0) then
          ! Reading the data.
          read(78) sensit_columns(1:nel)
          pos0 = pos0 + nel * 4
        endif

        ! Shift the file pointer to the next matrix row.
        pos0 = pos0 + nel * MATRIX_PRECISION

        ! Calculate the partitioning.
        nel_at_cpu = 0
        cpu = 1
        do p = 1, nel
          do while (sensit_columns(p) > sum(nelements_at_cpu(1:cpu)))
            cpu = cpu + 1
          enddo
          nel_at_cpu(cpu) = nel_at_cpu(cpu) + 1
        enddo

        if (sum(nel_at_cpu) /= nel) then
          call exit_MPI("Wrong nel_at_cpu in read_sensitivity_kernel!", myrank, 0)
        endif

        write(77) nel_at_cpu

      enddo ! model components loop
    enddo ! data components loop
  enddo ! data loop

  close(78)
  close(77)

  deallocate(sensit_columns)

  if (myrank == 0) print *, 'Finished partitioning of the sensitivity kernel.'

end subroutine partition_sensitivity_kernel

!=============================================================================================
! Reads the depth weight.
!=============================================================================================
subroutine read_depth_weight(par, column_weight, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(out) :: column_weight(par%nelements)

  integer :: problem_type
  character(len=256) :: filename, filename_full

  ! Define the problem type.
  select type(par)
  class is (t_parameters_grav)
    problem_type = 1
  class is (t_parameters_mag)
    problem_type = 2
  end select

  ! Define the file name.
  filename = "sensit_"//SUFFIX(problem_type)//"_weight"

  if (par%sensit_read /= 0) then
    filename_full = trim(par%sensit_path)//filename
  else
    filename_full = trim(path_output)//"/SENSIT/"//filename
  endif

  ! Read the depth weight from file.
  call read_depth_weight_file(par, filename_full, column_weight, myrank, nbproc)
end subroutine read_depth_weight

!=============================================================================================
! Reads the depth weight from file.
!=============================================================================================
subroutine read_depth_weight_file(par, filename, column_weight, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  character(len=*), intent(in) :: filename
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(out) :: column_weight(par%nelements)

  integer :: nelements_total, nelements_total_read
  integer :: ierr
  character(len=256) :: msg

  ! The full column weight.
  real(kind=CUSTOM_REAL), allocatable :: column_weight_full(:)

  if (myrank == 0) print *, "Reading the depth weight file ", trim(filename)

  nelements_total = par%nx * par%ny * par%nz

  ! Allocate the full array on master rank only.
  if (myrank == 0) then
    allocate(column_weight_full(nelements_total), source=0._CUSTOM_REAL, stat=ierr)
  else
    ! Fortran standard requires that allocatable array is allocated when passing by argument.
    allocate(column_weight_full(1), source=0._CUSTOM_REAL, stat=ierr)
  endif

  ! Read the full array by master rank only.
  if (myrank == 0) then
    ! Open the file.
    open(78, file=trim(filename), form='unformatted', status='old', action='read', access='stream', &
         iostat=ierr, iomsg=msg)

    if (ierr /= 0) call exit_MPI("Error in opening the depth weight file! path=" &
                                  //trim(filename)//", iomsg="//msg, myrank, ierr)

    read(78) nelements_total_read
    read(78) column_weight_full

    close(78)

    ! Consistency check.
    if (nelements_total_read /= nelements_total) then
      call exit_MPI("Depth weight file header does not match the Parfile!", myrank, 0)
    endif
  endif

  ! Scatter the depth weight to all CPUs.
  call scatter_full_array(par%nelements, column_weight_full, column_weight, myrank, nbproc)

  deallocate(column_weight_full)
end subroutine read_depth_weight_file

!=====================================================================================================
! Reads the sensitivity kernel metadata and retrieve nbproc used to calculate the sensitivity kernel.
!=====================================================================================================
subroutine read_sensitivity_metadata(par, nbproc_sensit, problem_type, myrank)
  class(t_parameters_base), intent(in) :: par
  integer, intent(in) :: problem_type
  integer, intent(in) :: myrank

  integer, intent(out) :: nbproc_sensit

  integer :: ierr
  character(len=256) :: filename, filename_full
  character(len=256) :: msg
  integer :: nx_read, ny_read, nz_read, ndata_read
  integer :: nmodel_components_read, ndata_components_read
  integer :: compression_type_read, weight_type_read
  integer :: precision_read
  real(kind=CUSTOM_REAL) :: comp_error

  ! Define the file name.
  filename = "sensit_"//SUFFIX(problem_type)//"_meta.txt"

  if (par%sensit_read /= 0) then
    filename_full = trim(par%sensit_path)//filename
  else
    filename_full = trim(path_output)//"/SENSIT/"//filename
  endif

  if (myrank == 0) then

    print *, "Reading the sensitivity metadata file ", trim(filename_full)

    ! Open the file.
    open(78, file=trim(filename_full), form='formatted', status='old', action='read', iostat=ierr, iomsg=msg)

    if (ierr /= 0) call exit_MPI("Error in opening the sensitivity metadata file! path=" &
                                  //trim(filename_full)//", iomsg="//msg, myrank, ierr)

    read(78, *) nx_read, ny_read, nz_read, ndata_read
    read(78, *) nbproc_sensit, precision_read, weight_type_read
    read(78, *) compression_type_read, comp_error
    read(78, *) nmodel_components_read, ndata_components_read

    if (myrank == 0) print *, "COMPRESSION ERROR (read) =", comp_error

    ! Consistency check.
    if (nx_read /= par%nx .or. ny_read /= par%ny .or. nz_read /= par%nz .or. &
        ndata_read /= par%ndata .or. weight_type_read /= par%depth_weighting_type .or. &
        nmodel_components_read /= par%nmodel_components .or. ndata_components_read /= par%ndata_components) then
      call exit_MPI("Sensitivity metadata file info does not match the Parfile!", myrank, 0)
    endif

    if (compression_type_read /= par%compression_type) then
      call exit_MPI("Compression type is inconsistent!", myrank, 0)
    endif

    ! Matrix precision consistency check.
    if (precision_read /= MATRIX_PRECISION) then
      call exit_MPI("Matrix precision is not consistent!", myrank, 0)
    endif

    close(78)
  endif

  call MPI_Bcast(nbproc_sensit, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

end subroutine read_sensitivity_metadata

!==========================================================================================================
! Applying the column weight to sensitivity line.
!==========================================================================================================
subroutine apply_column_weight(nelements, sensit_line, column_weight)
  integer, intent(in) :: nelements
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(nelements)

  real(kind=CUSTOM_REAL), intent(inout) :: sensit_line(nelements)

  integer :: i

  do i = 1, nelements
    sensit_line(i) = sensit_line(i) * column_weight(i)
  enddo

end subroutine apply_column_weight

end module sensitivity_gravmag
