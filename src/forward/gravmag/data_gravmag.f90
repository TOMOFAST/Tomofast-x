
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

  implicit none

  private

  type, public :: t_data

    ! Number of data points.
    integer :: ndata

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
subroutine data_initialize(this, ndata, myrank)
  class(t_data), intent(inout) :: this
  integer, intent(in) :: ndata, myrank
  integer :: ierr

  this%ndata = ndata

  ierr = 0

  if (.not. allocated(this%X)) allocate(this%X(this%ndata), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%Y)) allocate(this%Y(this%ndata), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%Z)) allocate(this%Z(this%ndata), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%val_meas)) allocate(this%val_meas(ndata_components, this%ndata), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%val_calc)) allocate(this%val_calc(ndata_components, this%ndata), source=0._CUSTOM_REAL, stat=ierr)

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
      read(10, *, iostat=ierr) this%X(i), this%Y(i), this%Z(i)
    else
      read(10, *, iostat=ierr) this%X(i), this%Y(i), this%Z(i), this%val_meas(:, i)
   endif

    if (ierr /= 0) call exit_MPI("Problem while reading file in data_read_points_format!", myrank, 0)
  enddo

end subroutine data_read_points_format

!================================================================================================
! Writes the data in two formats:
!   (1) to read by read_data();
!   (2) for Paraview visualization.
!
! which=1 - for measured data,
! which=2 - for calculated data.

! name_prefix = prefix for the file name.
!================================================================================================
subroutine data_write(this, name_prefix, which, myrank)
  class(t_data), intent(in) :: this
  character(len=*), intent(in) :: name_prefix
  integer, intent(in) :: which, myrank

  real(kind=CUSTOM_REAL) :: X, Y, Z
  real(kind=CUSTOM_REAL) :: val(ndata_components)
  integer :: i
  character(len=512) :: file_name, file_name2

  ! Write file by master CPU only.
  if (myrank /= 0) return

  file_name  = trim(path_output)//'/'//name_prefix//'data.txt'
  file_name2 = trim(path_output)//'/'//name_prefix//'data_csv.txt'

  print *, 'Writing data to file '//trim(file_name)

  open(10, file=trim(file_name), access='stream', form='formatted', status='replace', action='write')
  ! For Paraview.
  open(20, file=trim(file_name2), access='stream', form='formatted', status='replace', action='write')

  ! Writing a header line.
  write(10, *) this%ndata
  write(20, *) "x,y,z,f"

  ! Write data.
  do i = 1, this%ndata
    X = this%X(i)
    Y = this%Y(i)
    Z = this%Z(i)

    if (which == 1) then
      val = this%val_meas(:, i)
    else
      val = this%val_calc(:, i)
    endif

    write(10, *) X, Y, Z, val
    ! Note: revert Z-axis for Paraview.
    write(20, *) X, ", ", Y, ", ", -Z, ", ", val
  enddo

  close(10)
  close(20)

end subroutine data_write

end module data_gravmag
