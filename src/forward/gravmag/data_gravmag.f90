
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
  use parallel_tools

  implicit none

  private

  type, public :: t_data

    ! Total number of data points.
    integer :: ndata
    ! Local number of data points.
    integer :: ndata_loc

    ! Data positions.
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: X, Y, Z

    ! Data values measured and calculated (using a model from inversion).
    real(kind=CUSTOM_REAL), allocatable :: val_meas(:)
    real(kind=CUSTOM_REAL), allocatable :: val_calc(:)

  contains
    private

    procedure, public, pass :: initialize => data_initialize
    procedure, public, pass :: read => data_read
    procedure, public, pass :: write => data_write

    procedure, pass :: read_points_format => data_read_points_format

  end type t_data

contains

!============================================================================================================
! Initialize data object.
!============================================================================================================
subroutine data_initialize(this, ndata, ndata_loc, myrank)
  class(t_data), intent(inout) :: this
  integer, intent(in) :: ndata, ndata_loc, myrank
  integer :: ierr

  this%ndata = ndata
  this%ndata_loc = ndata_loc

  ierr = 0

  allocate(this%X(ndata_loc), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%Y(ndata_loc), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%Z(ndata_loc), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%val_meas(ndata_loc), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%val_calc(ndata_loc), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in data_initialize!", myrank, ierr)

end subroutine data_initialize

!============================================================================================================
! Read data (coordinates and values) in Universal Transverse Mercator (UTM)
! geographic map coordinate system.
!============================================================================================================
subroutine data_read(this, file_name, myrank, nbproc)
  class(t_data), intent(inout) :: this
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank, nbproc

  if (myrank == 0) print *, 'Reading data from file '//trim(file_name)

  call this%read_points_format(file_name, myrank, nbproc)

end subroutine data_read

!============================================================================================================
! Read data in points format.
!============================================================================================================
subroutine data_read_points_format(this, file_name, myrank, nbproc)
  class(t_data), intent(inout) :: this
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank, nbproc

  integer :: i, i_loc, ierr
  integer :: ndata_read, ndata_smaller

  open(unit=10, file=file_name, status='old', form='formatted', action='read', iostat=ierr)
  if (ierr /= 0) call exit_MPI("Error in opening the data file!", myrank, ierr)

  read(10, *) ndata_read

  if (ndata_read /= this%ndata) &
    call exit_MPI("The number of data in Parfile differs from the number of data in data file!", myrank, ndata_read)

  ndata_smaller = get_nsmaller(this%ndata_loc, myrank, nbproc)

  do i = 1, this%ndata
    if (i > ndata_smaller .and. i <= ndata_smaller + this%ndata_loc) then
      i_loc = i - ndata_smaller
      read(10, *, end=20, err=11) this%X(i_loc), this%Y(i_loc), this%Z(i_loc), this%val_meas(i_loc)
    else
      ! Skip this data on current rank.
      read(10, *)
    endif
  enddo

20 close(unit=10)

  return

11 call exit_MPI("Problem while reading the data!", myrank, 0)

end subroutine data_read_points_format

!================================================================================================
! Writes the data with grid in two formats:
!   (1) To read by read_data().
!   (2) For Paraview visualization.
!
! which=1 - for measured data.
! which=2 - for calculated data.

! name_prefix = prefix for the file name.
!================================================================================================
subroutine data_write(this, name_prefix, which, myrank, nbproc)
  class(t_data), intent(in) :: this
  character(len=*), intent(in) :: name_prefix
  integer, intent(in) :: which, myrank, nbproc

  integer :: i, ierr
  real(kind=CUSTOM_REAL) :: X, Y, Z, val
  character(len=512) :: file_name, file_name2
  integer :: alloc_size

  real(kind=CUSTOM_REAL), allocatable :: x_full(:)
  real(kind=CUSTOM_REAL), allocatable :: y_full(:)
  real(kind=CUSTOM_REAL), allocatable :: z_full(:)
  real(kind=CUSTOM_REAL), allocatable :: val_full(:)

  !----------------------------------------------------------------------------------------
  ! Gather the full daata grid and values.
  !----------------------------------------------------------------------------------------
  if (myrank == 0) then
    ! Need full arrays only on the master. Allocate to zero size on other ranks for correctness.
    alloc_size = this%ndata
  else
    alloc_size = 0
  endif

  allocate(x_full(alloc_size), source=0._CUSTOM_REAL, stat=ierr)
  allocate(y_full(alloc_size), source=0._CUSTOM_REAL, stat=ierr)
  allocate(z_full(alloc_size), source=0._CUSTOM_REAL, stat=ierr)
  allocate(val_full(alloc_size), source=0._CUSTOM_REAL, stat=ierr)

  ! Gather the full data grid.
  call get_full_array(this%X, this%ndata_loc, x_full, .false., myrank, nbproc)
  call get_full_array(this%Y, this%ndata_loc, y_full, .false., myrank, nbproc)
  call get_full_array(this%Z, this%ndata_loc, z_full, .false., myrank, nbproc)

  ! Gather the full data values array.
  if (which == 1) then
    call get_full_array(this%val_meas, this%ndata_loc, val_full, .false., myrank, nbproc)
  else
    call get_full_array(this%val_calc, this%ndata_loc, val_full, .false., myrank, nbproc)
  endif

  !----------------------------------------------------------------------------------------
  ! Write data to file.
  !----------------------------------------------------------------------------------------

  ! Write file by master CPU only.
  if (myrank == 0) then

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
      X = x_full(i)
      Y = y_full(i)
      Z = z_full(i)
      val = val_full(i)

      write(10, *) X, Y, Z, val
      ! Note: revert Z-axis for Paraview.
      write(20, *) X, ", ", Y, ", ", -Z, ", ", val
    enddo

    close(10)
    close(20)

  endif

  deallocate(x_full)
  deallocate(y_full)
  deallocate(z_full)
  deallocate(val_full)

end subroutine data_write

end module data_gravmag
