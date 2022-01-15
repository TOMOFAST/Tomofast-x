
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

  implicit none

  private

  public :: calculate_sensitivity_kernel
  public :: predict_sensit_kernel_size

  private :: calculate_sensitivity
  private :: apply_column_weight

  public :: calculate_and_write_sensit
  public :: read_sensitivity_kernel

contains

!=============================================================================================
! Calculates the sensitivity kernel and adds it to a sparse matrix.
!=============================================================================================
subroutine calculate_sensitivity_kernel(par, grid, data, column_weight, sensit_matrix, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  type(t_grid), intent(in) :: grid
  type(t_data), intent(in) :: data
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  integer, intent(in) :: myrank, nbproc

  ! Sensitivity matrix.
  type(t_sparse_matrix), intent(inout) :: sensit_matrix

  integer :: nnz
  logical :: STORE_KERNEL

  STORE_KERNEL = .true.

  call calculate_sensitivity(par, grid, data, column_weight, sensit_matrix, &
                             STORE_KERNEL, nnz, myrank, nbproc)

end subroutine calculate_sensitivity_kernel

!==================================================================================================
! Calculates the compressed sensitivity kernel size.
!==================================================================================================
function predict_sensit_kernel_size(par, grid, data, column_weight, myrank, nbproc) result (nnz)
  class(t_parameters_base), intent(in) :: par
  type(t_grid), intent(in) :: grid
  type(t_data), intent(in) :: data
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  integer, intent(in) :: myrank, nbproc

  integer :: nnz

  type(t_sparse_matrix) :: dummy_matrix
  logical :: STORE_KERNEL

  STORE_KERNEL = .false.

  call calculate_sensitivity(par, grid, data, column_weight, dummy_matrix, &
                             STORE_KERNEL, nnz, myrank, nbproc)

end function predict_sensit_kernel_size

!=============================================================================================
! Calculates the sensitivity kernel (parallelized by data) and writes it to files.
!=============================================================================================
subroutine calculate_and_write_sensit(par, grid_full, data, column_weight, nnz, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  type(t_grid), intent(in) :: grid_full
  type(t_data), intent(in) :: data
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  integer, intent(in) :: myrank, nbproc

  ! The number of non-zero elements (on current CPU) in the sensitivity kernel parallelized by model.
  ! We need this number for reading the sensitivity later from files, for invere problem.
  ! As the inverse problem is parallelized by model, and calculations here are parallelized by data.
  integer, intent(out) :: nnz

  type(t_magnetic_field) :: mag_field
  type(t_parallel_tools) :: pt
  integer :: i, p, nel, ierr
  integer :: nelements_total, nnz_data
  integer :: problem_type
  integer :: idata, ndata_loc, ndata_smaller
  character(len=256) :: filename, filename_full
  real(kind=CUSTOM_REAL) :: comp_rate, nnz_total_dbl, nnz_total_dbl2

  integer :: nnz_model_loc(nbproc), nnz_model(nbproc)
  ! The number of elements on every CPU.
  integer :: nelements_at_cpu(nbproc)
  integer :: nsmaller_at_cpu(nbproc)
  integer :: cpu

  ! Sensitivity matrix row.
  real(kind=CUSTOM_REAL), allocatable :: sensit_line_full(:)
  real(kind=CUSTOM_REAL), allocatable :: dummy1(:), dummy2(:)

  ! Arrays for storing the compressed sensitivity line.
  integer, allocatable :: sensit_columns(:)
  real(kind=CUSTOM_REAL), allocatable :: sensit_compressed(:)

  ! The full column weight.
  real(kind=CUSTOM_REAL), allocatable :: column_weight_full(:)

  ! Define the problem type.
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

  !---------------------------------------------------------------------------------------------
  ! Define the output file.
  !---------------------------------------------------------------------------------------------
  call system('mkdir -p '//trim(path_output)//"/SENSIT/")

  filename = "sensit_"//trim(str(nbproc))//"_"//trim(str(myrank))
  filename_full = trim(path_output)//"/SENSIT/"//filename

  print *, 'Writing the sensitivity to file ', trim(filename_full)

  open(77, file=trim(filename_full), access='stream', form='formatted', status='unknown', action='write')

  !---------------------------------------------------------------------------------------------
  ! Allocate memory.
  !---------------------------------------------------------------------------------------------
  nelements_total = par%nx * par%ny * par%nz

  allocate(sensit_line_full(nelements_total), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_and_write_sensit!", myrank, ierr)

  allocate(column_weight_full(nelements_total), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_and_write_sensit!", myrank, ierr)

  allocate(sensit_columns(nelements_total), source=0, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_and_write_sensit!", myrank, ierr)

  allocate(sensit_compressed(nelements_total), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_and_write_sensit!", myrank, ierr)

  !---------------------------------------------------------------------------------------------
  ! Calculate sensitivity lines.
  !---------------------------------------------------------------------------------------------
  call pt%get_full_array(column_weight, par%nelements, column_weight_full, .true., myrank, nbproc)

  ndata_loc = pt%calculate_nelements_at_cpu(par%ndata, myrank, nbproc)
  ndata_smaller = pt%get_nsmaller(ndata_loc, myrank, nbproc)

  print *, "Calculating sensitivity for myrank, ndata_loc =", myrank, ndata_loc

  ! File header.
  write(77, *) ndata_loc, par%ndata, nelements_total, myrank, nbproc, par%wavelet_threshold

  ! Calculate the number of elements on every CPU.
  nelements_at_cpu = pt%get_number_elements_on_other_cpus(par%nelements, myrank, nbproc)

  ! Calculate the number of elements on ranks smaller than current.
  do i = 1, nbproc
    nsmaller_at_cpu(i) = sum(nelements_at_cpu(1:i-1))
  enddo

  nnz_data = 0
  nnz_model_loc = 0

  ! Loop over the local data lines.
  do i = 1, ndata_loc
    ! Global data index.
    idata = ndata_smaller + i

    if (problem_type == 1) then
    ! Gravity problem.
      call graviprism_full(nelements_total, par%ncomponents, grid_full, data%X(idata), data%Y(idata), data%Z(idata), &
                           dummy1, dummy2, sensit_line_full, myrank)

    else if (problem_type == 2) then
    ! Magnetic problem.
      call mag_field%magprism(nelements_total, idata, grid_full, data%X, data%Y, data%Z, sensit_line_full)
    endif

    ! Applying the depth weight.
    call apply_column_weight(nelements_total, sensit_line_full, column_weight_full)

    if (par%compression_type > 0) then
    ! Wavelet compression.
      call Haar3D_serial(sensit_line_full, par%nx, par%ny, par%nz)

      nel = 0
      cpu = 1
      do p = 1, nelements_total
        ! Partitioning for parallelization by model.
        if (cpu < nbproc .and. p > nsmaller_at_cpu(cpu + 1)) then
          cpu = cpu + 1
        endif

        if (abs(sensit_line_full(p)) >= par%wavelet_threshold) then
        ! Store sensitivity elements greater than the wavelet threshold.
          nel = nel + 1
          sensit_columns(nel) = p
          sensit_compressed(nel) = sensit_line_full(p)

          nnz_model_loc(cpu) = nnz_model_loc(cpu) + 1
        endif
      enddo

    else
    ! No compression.
      nel = nelements_total
      sensit_compressed = sensit_line_full
      do p = 1, nelements_total
        sensit_columns(p) = p
      enddo
      nnz_model_loc = nnz_model_loc + nelements_at_cpu
    endif

    ! The sensitivity kernel size (when parallelized by data).
    nnz_data = nnz_data + nel

    ! Sanity check.
    if (nnz_data < 0) then
      call exit_MPI("Integer overflow in nnz_data! Increase the wavelet threshold or the number of CPUs.", myrank, nnz_data)
    endif

    ! Write the sensitivity line to file.
    write(77, *) idata, nel
    if (nel > 0) then
      write(77, *) sensit_columns(1:nel)
      write(77, *) sensit_compressed(1:nel)
    endif

    ! Print the progress.
    if (myrank == 0 .and. mod(i, int(0.1d0 * ndata_loc)) == 0) then
      print *, 'Percents completed: ', (i / int(0.1d0 * ndata_loc)) * 10 ! Approximate percents.
    endif
  enddo ! data loop

  close(77)

  !---------------------------------------------------------------------------------------------
  ! Calculate the kernel compression rate.
  !---------------------------------------------------------------------------------------------
  call mpi_allreduce(dble(nnz_data), nnz_total_dbl, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  comp_rate = nnz_total_dbl / dble(nelements_total) / dble(par%ndata)

  if (myrank == 0) print *, 'COMPRESSION RATE = ', comp_rate

  !---------------------------------------------------------------------------------------------
  ! Calculate the nnz for the sensitivity kernel parallelized by model.
  !---------------------------------------------------------------------------------------------
  call mpi_allreduce(nnz_model_loc, nnz_model, nbproc, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (myrank == 0)  print *, 'nnz_model = ', nnz_model

  ! Sanity check.
  do i = 1, nbproc
    if (nnz_model(i) < 0) then
      call exit_MPI("Integer overflow in nnz_model! Increase the wavelet threshold or the number of CPUs.", myrank, nnz_model(i))
    endif
  enddo

  ! Sanity check.
  nnz_total_dbl2 = sum(dble(nnz_model))
  if (nnz_total_dbl /= nnz_total_dbl2) then
    print *, myrank, nnz_total_dbl, nnz_total_dbl2
    call exit_MPI("Wrong nnz_model in calculate_and_write_sensit!", myrank, 0)
  endif

  ! Return the nnz for the current CPU.
  nnz = nnz_model(myrank + 1)

  !---------------------------------------------------------------------------------------------
  deallocate(sensit_line_full)
  deallocate(column_weight_full)
  deallocate(sensit_columns)
  deallocate(sensit_compressed)

  if (myrank == 0) print *, 'Finished calculating the sensitivity kernel.'

end subroutine calculate_and_write_sensit

!=============================================================================================
! Reads the sensitivity kernel from files,
! and stores it in the sparse matrix parallelized by model.
!=============================================================================================
subroutine read_sensitivity_kernel(par, sensit_matrix, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  integer, intent(in) :: myrank, nbproc

  ! Sensitivity matrix.
  type(t_sparse_matrix), intent(inout) :: sensit_matrix

  type(t_parallel_tools) :: pt

  ! Arrays for storing the compressed sensitivity line.
  integer, allocatable :: sensit_columns(:)
  real(kind=CUSTOM_REAL), allocatable :: sensit_compressed(:)

  integer :: i, nsmaller, ierr
  integer :: rank, nelements_total
  character(len=256) :: filename, filename_full
  character(len=256) :: msg

  integer :: ndata_loc, ndata_read, nelements_total_read, myrank_read, nbproc_read
  real(kind=CUSTOM_REAL) :: wavelet_threshold_read
  integer :: idata, nel

  !---------------------------------------------------------------------------------------------
  ! Allocate memory.
  !---------------------------------------------------------------------------------------------
  nelements_total = par%nx * par%ny * par%nz

  allocate(sensit_columns(nelements_total), source=0, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_and_write_sensit!", myrank, ierr)

  allocate(sensit_compressed(nelements_total), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in calculate_and_write_sensit!", myrank, ierr)

  !---------------------------------------------------------------------------------------------
  ! Reading the files.
  !---------------------------------------------------------------------------------------------

  ! The number of elements on CPUs with rank smaller than myrank.
  nsmaller = pt%get_nsmaller(par%nelements, myrank, nbproc)

  do rank = 0, nbproc - 1
    ! Form the file name.
    filename = "sensit_"//trim(str(nbproc))//"_"//trim(str(rank))
    filename_full = trim(path_output)//"/SENSIT/"//filename

    if (myrank == 0) print *, 'Reading the sensitivity file ', trim(filename_full)

    open(78, file=trim(filename_full), status='old', action='read', iostat=ierr, iomsg=msg)
    if (ierr /= 0) call exit_MPI("Error in opening the model file! path=" &
                                 //filename_full//" iomsg="//msg, myrank, ierr)

    read(78, *) ndata_loc, ndata_read, nelements_total_read, myrank_read, nbproc_read, wavelet_threshold_read

    ! Sanity check.
    if (ndata_read /= par%ndata .or. nelements_total_read /= nelements_total &
        .or. myrank_read /= rank .or. nbproc_read /= nbproc) then
      call exit_MPI("Wrong file header in read_sensitivity_kernel!", myrank, 0)
    endif

    do i = 1, ndata_loc
      read(78, *) idata, nel
      if (nel > 0) then
        read(78, *) sensit_columns(1:nel)
        read(78, *) sensit_compressed(1:nel)
      endif
      if (myrank == 0) print *, idata, nel
    enddo

    close(78)
  enddo

  !---------------------------------------------------------------------------------------------
  deallocate(sensit_columns)
  deallocate(sensit_compressed)

  if (myrank == 0) print *, 'Finished reading the sensitivity kernel.'

end subroutine read_sensitivity_kernel

!=============================================================================================
! Calculates the sensitivity kernel / or predicts its size without storing the kernel,
! depending on the flag STORE_KERNEL.
!=============================================================================================
subroutine calculate_sensitivity(par, grid, data, column_weight, sensit_matrix, &
                                 STORE_KERNEL, nnz_local, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  type(t_grid), intent(in) :: grid
  type(t_data), intent(in) :: data
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  logical, intent(in) :: STORE_KERNEL
  integer, intent(in) :: myrank, nbproc

  ! The number of non-zero elements in the compressed sensitivity kernel on current CPU.
  integer, intent(out) :: nnz_local

  ! Sensitivity matrix.
  type(t_sparse_matrix), intent(inout) :: sensit_matrix

  type(t_magnetic_field) :: mag_field
  integer :: i, p, ierr
  real(kind=CUSTOM_REAL) :: comp_rate
  integer :: nsmaller, nelements_total
  type(t_parallel_tools) :: pt
  integer :: problem_type
  integer :: nnz_line
  real(kind=CUSTOM_REAL) :: nnz_total_dbl

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

  if (par%compression_type > 0) then
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
  ! Calculating sensitivity and adding to the sparse matrix / or calculating nnz for current CPU.
  nnz_local = 0

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

    if (par%compression_type > 0) then
    ! Wavelet compression.
      if (nbproc > 1) then
      ! Parallel wavelet copression.
        call pt%get_full_array(sensit_line, par%nelements, sensit_line_full, .true., myrank, nbproc)
        call Haar3D(sensit_line_full, par%nx, par%ny, par%nz, myrank, nbproc)

        ! Extract the local sensitivity part.
        sensit_line = sensit_line_full(nsmaller + 1 : nsmaller + par%nelements)
      else
      ! Serial.
        call Haar3D(sensit_line, par%nx, par%ny, par%nz, myrank, nbproc)
      endif

      ! Set values below the threshold to zero.
      do p = 1, par%nelements
        if (abs(sensit_line(p)) < par%wavelet_threshold) then
          sensit_line(p) = 0.d0
        endif
      enddo
    endif

    if (STORE_KERNEL) then
    ! Adding the sensitivity kernel the a sparse matrix.

      ! Sanity check: check if we have enough space in the matrix for new elemements.
      nnz_line = count(sensit_line /= 0.d0)
      if (sensit_matrix%get_number_elements() + nnz_line > sensit_matrix%get_nnz()) then
        call exit_MPI("The matrix size is too small, exiting!", myrank, ierr)
      endif

      call sensit_matrix%new_row(myrank)

      do p = 1, par%nelements
        ! Adding the Z-component only.
        call sensit_matrix%add(sensit_line(p), p, myrank)
      enddo
    endif

    ! The sensitivity kernel size.
    nnz_local = nnz_local + count(sensit_line /= 0.d0)

    ! Printing the progress.
    if (mod(i, int(0.1d0 * par%ndata)) == 0) then
      if (myrank == 0) print *, 'Percents completed: ', (i / int(0.1d0 * par%ndata)) * 10 ! Approximate percents.
    endif

  enddo ! data loop

  ! Sanity check.
  if (nnz_local < 0) then
    call exit_MPI("Integer overflow in nnz_local! Increase the wavelet threshold or the number of CPUs.", myrank, nnz_local)
  endif

  if (STORE_KERNEL) then
    call sensit_matrix%finalize(par%nelements, myrank)
  endif

  ! Calculate the kernel compression rate.
  if (nbproc > 1) then
    call mpi_allreduce(dble(nnz_local), nnz_total_dbl, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    comp_rate = nnz_total_dbl / dble(nelements_total) / dble(par%ndata)
  else
    comp_rate = dble(nnz_local) / dble(par%nelements) / dble(par%ndata)
  endif

  if (STORE_KERNEL) then
    if (myrank == 0) print *, 'COMPRESSION RATE = ', comp_rate
  else
    if (myrank == 0) print *, 'COMPRESSION RATE (estim) = ', comp_rate
  endif

  deallocate(sensit_line)
  deallocate(sensit_line2)
  deallocate(sensit_line3)
  if (allocated(sensit_line_full)) deallocate(sensit_line_full)

  if (myrank == 0) print *, 'Finished calculating the sensitivity kernel.'

end subroutine calculate_sensitivity

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
