
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
! A class to work with data for parallel inversion.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!===============================================================================================
module data_gravmag

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use string
  use paraview

  implicit none

  private

  type, public :: t_data

    ! Number of data points.
    integer :: ndata
    integer :: ncomponents

    ! Units multiplier.
    real(kind=CUSTOM_REAL) :: units_mult

    ! Data positions.
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: X, Y, Z

    ! Data values measured and calculated (using a model from inversion).
    real(kind=CUSTOM_REAL), allocatable :: val_meas(:, :)
    real(kind=CUSTOM_REAL), allocatable :: val_calc(:, :)

  contains
    private

    procedure, public, pass :: initialize => data_initialize
    procedure, public, pass :: read => data_read
    procedure, public, pass :: read_grid => data_read_grid
    procedure, public, pass :: write => data_write

    procedure, pass :: broadcast => data_broadcast
    procedure, pass :: read_points_format => data_read_points_format

  end type t_data

contains

!============================================================================================================
! Initialize data object.
!============================================================================================================
subroutine data_initialize(this, ndata, ncomponents, units_mult, myrank)
  class(t_data), intent(inout) :: this
  integer, intent(in) :: ndata, ncomponents, myrank
  real(kind=CUSTOM_REAL) :: units_mult
  integer :: ierr

  this%ndata = ndata
  this%ncomponents = ncomponents
  this%units_mult = units_mult

  ierr = 0

  if (.not. allocated(this%X)) allocate(this%X(this%ndata), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%Y)) allocate(this%Y(this%ndata), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%Z)) allocate(this%Z(this%ndata), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%val_meas)) allocate(this%val_meas(this%ncomponents, this%ndata), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%val_calc)) allocate(this%val_calc(this%ncomponents, this%ndata), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in data_initialize!", myrank, ierr)

end subroutine data_initialize

!============================================================================================================
! Broadcasts data arrays.
!============================================================================================================
subroutine data_broadcast(this, myrank)
  class(t_data), intent(inout) :: this
  integer, intent(in) :: myrank
  integer :: ierr

  ierr = 0

  call MPI_Bcast(this%X, this%ndata, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%Y, this%ndata, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%Z, this%ndata, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%val_meas, size(this%val_meas), CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI_Bcast error in data_broadcast!", myrank, ierr)

end subroutine data_broadcast

!============================================================================================================
! Read data (coordinates and values) in Universal Transverse Mercator (UTM)
! geographic map coordinate system.
!============================================================================================================
subroutine data_read(this, file_name, myrank)
  class(t_data), intent(inout) :: this
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank

  if (myrank == 0) then
    print *, 'Reading data from file '//trim(file_name)
    call this%read_points_format(file_name, .false., myrank)
  endif

  ! MPI broadcast data arrays.
  call this%broadcast(myrank)

end subroutine data_read

!============================================================================================================
! Read data grid.
!============================================================================================================
subroutine data_read_grid(this, file_name, myrank)
  class(t_data), intent(inout) :: this
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank

  if (myrank == 0) then
    print *, 'Reading data grid from file '//trim(file_name)
    call this%read_points_format(file_name, .true., myrank)
  endif

  ! MPI broadcast data arrays.
  call this%broadcast(myrank)

end subroutine data_read_grid

!============================================================================================================
! Read data in points format.
!============================================================================================================
subroutine data_read_points_format(this, file_name, grid_only, myrank)
  class(t_data), intent(inout) :: this
  character(len=*), intent(in) :: file_name
  logical, intent(in) :: grid_only
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL) :: dummy(this%ncomponents)
  integer :: i, ierr
  integer :: ndata_in_file

  ! Reading my master CPU only.
  if (myrank /= 0) return

  open(unit=10, file=file_name, status='old', form='formatted', action='read', iostat=ierr)
  if (ierr /= 0) call exit_MPI("Error in opening the data file!", myrank, ierr)

  read(10, *) ndata_in_file

  if (ndata_in_file /= this%ndata) &
    call exit_MPI("The number of data in Parfile differs from the number of data in data file!", myrank, ndata_in_file)

  do i = 1, this%ndata
    if (grid_only) then
      ! Read dummy values too to check for file format consistency.
      read(10, *, iostat=ierr) this%X(i), this%Y(i), this%Z(i), dummy(:)
    else
      read(10, *, iostat=ierr) this%X(i), this%Y(i), this%Z(i), this%val_meas(:, i)
   endif

    if (ierr /= 0) call exit_MPI("Problem while reading the data file! Verify the number of data components.", myrank, 0)
  enddo

  ! Convert input data units.
  this%val_meas = this%units_mult * this%val_meas

  close(10)

end subroutine data_read_points_format

!================================================================================================
! Writes the data in two formats:
!   (1) ACSII columns: same as the input data format in Tomofastx.
!   (2) Binary VTK for Paraview visualization.
!
! which=1 - for measured data,
! which=2 - for calculated data.

! name_prefix = prefix for the file name.
!================================================================================================
subroutine data_write(this, name_prefix, which, myrank)
  class(t_data), intent(in) :: this
  character(len=*), intent(in) :: name_prefix
  integer, intent(in) :: which, myrank

  integer :: i
  character(len=512) :: file_name

  ! Write file by master CPU only.
  if (myrank /= 0) return

  file_name  = trim(path_output)//'/'//name_prefix//'data.txt'

  print *, 'Writing data to file '//trim(file_name)

  open(10, file=trim(file_name), access='stream', form='formatted', status='replace', action='write')

  ! Writing a header line.
  write(10, *) this%ndata

  ! Write data.
  if (which == 1) then
    write(10, *) (this%X(i), this%Y(i), this%Z(i), this%val_meas(:, i) / this%units_mult, new_line('a'), i = 1, this%ndata)
  else
    write(10, *) (this%X(i), this%Y(i), this%Z(i), this%val_calc(:, i) / this%units_mult, new_line('a'), i = 1, this%ndata)
  endif

  close(10)

  !------------------------------------------------------------------------------------
  ! Write data in VTK format for Paraview.
  !------------------------------------------------------------------------------------
  file_name  = 'data_'//name_prefix(1:len(name_prefix) - 1)//'.vtk'

  if (which == 1) then
    call visualisation_paraview_points(file_name, myrank, this%ndata, this%ncomponents, &
                                       this%val_meas, this%X, this%Y, this%Z, .true., this%units_mult)
  else
    call visualisation_paraview_points(file_name, myrank, this%ndata, this%ncomponents, &
                                       this%val_calc, this%X, this%Y, this%Z, .true., this%units_mult)
  endif

end subroutine data_write

end module data_gravmag
