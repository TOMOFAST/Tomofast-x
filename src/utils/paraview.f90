
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
! Routines for exporting the model data to the vtk format for visualization in Paraview.
!
! Vitaliy Ogarko, UWA, Australia.
!========================================================================================

module paraview

  use global_typedefs
  use mpi_tools, only: exit_MPI

  implicit none

  private

  public :: visualisation_paraview_legogrid
  public :: visualisation_paraview_struct_grid
  public :: visualisation_paraview_points

  private :: index_included

contains

!============================================================================================================
function index_included(p, i_index, j_index, k_index, i1, i2, j1, j2, k1, k2) result(included)
  integer, intent(in) :: p
  integer, intent(in) :: i_index(:), j_index(:), k_index(:)
  integer, intent(in) :: i1, i2, j1, j2, k1, k2
  logical :: included

  integer :: i, j, k

  i = i_index(p)
  j = j_index(p)
  k = k_index(p)

  if ((i1 <= i .and. i <= i2) .and. &
      (j1 <= j .and. j <= j2) .and. &
      (k1 <= k .and. k <= k2)) then
    included = .true.
  else
    included = .false.
  endif
end function index_included

!=======================================================================================================================
! Writes the model in binary VTK format for Paraview visualization.
! A slice of the full model can be specified via i1, i2, j1, j2, k1, k2 grid indexes.
! The Paraview structured grid datastructure is used, which allows additional analysis in Paraview (like contours).
!=======================================================================================================================
subroutine visualisation_paraview_struct_grid(filename, myrank, nelements, ncomponents, val, X1, Y1, Z1, X2, Y2, Z2, &
                                           i_index, j_index, k_index, &
                                           i1, i2, j1, j2, k1, k2, &
                                           INVERT_Z_AXIS, units_mult)
  ! MPI rank of this process.
  integer, intent(in) :: myrank
  ! Total number of cells.
  integer, intent(in) :: nelements, ncomponents
  ! Values for visualization.
  real(kind=CUSTOM_REAL), intent(in) :: val(nelements, ncomponents)
  ! Coordinates of points in the grid.
  real(kind=CUSTOM_REAL), intent(in) :: X1(nelements), Y1(nelements), Z1(nelements)
  real(kind=CUSTOM_REAL), intent(in) :: X2(nelements), Y2(nelements), Z2(nelements)
  integer, intent(in) :: i_index(nelements), j_index(nelements), k_index(nelements)
  integer, intent(in) :: i1, i2, j1, j2, k1, k2
  logical, intent(in) :: INVERT_Z_AXIS
  ! Units multiplier.
  real(kind=CUSTOM_REAL), intent(in) :: units_mult
  ! Output file name.
  character(len=*), intent(in) :: filename

  character(len=256) :: filename_full
  character(len=256) :: msg
  ! I/O error code.
  integer :: ierr
  integer :: npoints, nelements_slice
  integer :: j, p

  real(kind=4), allocatable :: point_centers(:, :)
  real(kind=4), allocatable :: point_data(:, :)
  real(kind=4) :: z_sign

  character :: lf*1, str1*8, str2*8, str3*8
  ! Line feed character.
  lf = char(10)

  ! Create a directory.
  call execute_command_line('mkdir -p "'//trim(path_output)//'/Paraview"')

  filename_full = trim(path_output)//"/Paraview/"//filename

  open(333, file=filename_full, status='replace', access='stream', form='unformatted', action='write', &
       iostat=ierr, iomsg=msg)

  if (ierr /= 0) call exit_MPI("Error with writing the VTK file! path="&
                               //trim(filename_full)//", iomsg="//msg, myrank, ierr)

  !-----------------------------------------------------------------
  ! Write a header.
  !-----------------------------------------------------------------
  write(333) '# vtk DataFile Version 3.0'//lf
  write(333) 'Tomofast-x'//lf
  write(333) 'BINARY'//lf
  write(333) 'DATASET STRUCTURED_GRID'//lf

  write(str1(1:8),'(i8)') (i2 - i1 + 1)
  write(str2(1:8),'(i8)') (j2 - j1 + 1)
  write(str3(1:8),'(i8)') (k2 - k1 + 1)

  write(333) 'DIMENSIONS '//str1//' '//str2//' '//str3//lf

  ! Number of elements in the requested slice.
  nelements_slice = (i2 - i1 + 1) * (j2 - j1 + 1) * (k2 - k1 + 1)

  npoints = nelements_slice

  write(str1(1:8),'(i8)') npoints
  write(333) lf//lf//'POINTS '//str1//' FLOAT'//lf

  !-----------------------------------------------------------------
  ! Allocate memory.
  !-----------------------------------------------------------------
  allocate(point_centers(3, nelements_slice), stat=ierr)
  ! Note we need the first dimension equal to the number of components.
  allocate(point_data(ncomponents, nelements_slice), stat=ierr)

  !-------------------------------------------------------------------
  ! Build the grid.
  !-------------------------------------------------------------------
  if (INVERT_Z_AXIS) then
    z_sign = -1.0
  else
    z_sign = 1.0
  endif

  j = 0
  do p = 1, nelements
    if (index_included(p, i_index, j_index, k_index, i1, i2, j1, j2, k1, k2)) then
      j = j + 1

      point_centers(1, j) = real(0.5 * (X1(p) + X2(p)))
      point_centers(2, j) = real(0.5 * (Y1(p) + Y2(p)))
      point_centers(3, j) = real(0.5 * (Z1(p) + Z2(p)))

      ! Flip the Z-axis of the grid.
      point_centers(3, j) = z_sign * point_centers(3, j)

      point_data(:, j) = real(val(p, :), 4)

      ! Units conversion.
      point_data(:, j) = point_data(:, j) / real(units_mult, 4)

      if (ncomponents == 3) then
        ! Flip the Z-axis of a vector.
        point_data(3, j) = z_sign * point_data(3, j)
      endif
    endif
  enddo

  ! Write the grid to a file.
  write(333) point_centers

  !-------------------------------------------------------------------
  ! Generate element data values.
  !-------------------------------------------------------------------
  write(str1(1:8),'(i8)') nelements_slice

  write(333) lf//lf//'POINT_DATA '//str1//lf

  if (ncomponents == 1) then
  ! Scalar data.
    write(333) 'SCALARS F FLOAT'//lf
    write(333) 'LOOKUP_TABLE default'//lf

  else if (ncomponents == 3) then
  ! Vector data.
    write(333) 'VECTORS vectors FLOAT'//lf
  endif

  write(333) point_data

  close(333)

  deallocate(point_centers)
  deallocate(point_data)

end subroutine visualisation_paraview_struct_grid

!====================================================================================================================
! Writes the model in binary VTK format for Paraview visualization.
! A slice of the full model can be specified via i1, i2, j1, j2, k1, k2 grid indexes.
! An arbitrary "lego" grid is considered.
!====================================================================================================================
subroutine visualisation_paraview_legogrid(filename, myrank, nelements, ncomponents, val, X1, Y1, Z1, X2, Y2, Z2, &
                                           i_index, j_index, k_index, &
                                           i1, i2, j1, j2, k1, k2, &
                                           INVERT_Z_AXIS, units_mult)
  ! MPI rank of this process.
  integer, intent(in) :: myrank
  ! Total number of cells.
  integer, intent(in) :: nelements, ncomponents
  ! Values for visualization.
  real(kind=CUSTOM_REAL), intent(in) :: val(nelements, ncomponents)
  ! Coordinates of points in the grid.
  real(kind=CUSTOM_REAL), intent(in) :: X1(nelements), Y1(nelements), Z1(nelements)
  real(kind=CUSTOM_REAL), intent(in) :: X2(nelements), Y2(nelements), Z2(nelements)
  integer, intent(in) :: i_index(nelements), j_index(nelements), k_index(nelements)
  integer, intent(in) :: i1, i2, j1, j2, k1, k2
  logical, intent(in) :: INVERT_Z_AXIS
  ! Units multiplier.
  real(kind=CUSTOM_REAL), intent(in) :: units_mult
  ! Output file name.
  character(len=*), intent(in) :: filename

  character(len=256) :: filename_full
  character(len=256) :: msg
  ! I/O error code.
  integer :: ierr
  integer :: npoints, nelements_slice
  integer :: i, j, p
  real(kind=CUSTOM_REAL) :: xyzgrid(3, 8)
  real(kind=CUSTOM_REAL) :: Z1_p, Z2_p, z_sign

  real(kind=4), allocatable :: xyzgrid_all(:, :, :)
  real(kind=4), allocatable :: cell_data(:, :)
  integer, allocatable :: cell_indexes(:, :)
  integer, allocatable :: cell_type(:)

  character :: lf*1, str1*8, str2*8
  lf = char(10) ! line feed character

  ! Create a directory.
  call execute_command_line('mkdir -p "'//trim(path_output)//'/Paraview"')

  filename_full = trim(path_output)//"/Paraview/"//filename

  open(333, file=filename_full, status='replace', access='stream', form='unformatted', action='write', &
       iostat=ierr, iomsg=msg)

  if (ierr /= 0) call exit_MPI("Error with writing the VTK file! path="&
                               //trim(filename_full)//", iomsg="//msg, myrank, ierr)

  !-----------------------------------------------------------------
  ! Write a header.
  !-----------------------------------------------------------------
  write(333) '# vtk DataFile Version 3.0'//lf
  write(333) 'TOMOFAST-X'//lf
  write(333) 'BINARY'//lf
  write(333) 'DATASET UNSTRUCTURED_GRID'//lf//lf

  ! Number of elements in the requested slice.
  nelements_slice = (i2 - i1 + 1) * (j2 - j1 + 1) * (k2 - k1 + 1)

  npoints = 8 * nelements_slice

  write(str1(1:8),'(i8)') npoints
  write(333) 'POINTS '//str1//' FLOAT'//lf

  !-----------------------------------------------------------------
  ! Allocate memory.
  !-----------------------------------------------------------------
  allocate(xyzgrid_all(3, 8, nelements_slice), stat=ierr)
  ! Note that we need the first dimension equal to the number of components.
  allocate(cell_data(ncomponents, nelements_slice), stat=ierr)
  allocate(cell_indexes(9, nelements_slice), stat=ierr)
  allocate(cell_type(nelements_slice), stat=ierr)

  !-----------------------------------------------------------------
  ! Build lego-grid.
  !-----------------------------------------------------------------
  if (INVERT_Z_AXIS) then
    z_sign = - 1._CUSTOM_REAL
  else
    z_sign = 1._CUSTOM_REAL
  endif

  j = 0
  do p = 1, nelements
    if (index_included(p, i_index, j_index, k_index, i1, i2, j1, j2, k1, k2)) then
      j = j + 1

      ! Flip the Z-axis of the grid.
      Z1_p = z_sign * Z1(p)
      Z2_p = z_sign * Z2(p)

      ! z = 1
      xyzgrid(1, 1) = X1(p)
      xyzgrid(2, 1) = Y1(p)
      xyzgrid(3, 1) = Z1_p

      xyzgrid(1, 2) = X2(p)
      xyzgrid(2, 2) = Y1(p)
      xyzgrid(3, 2) = Z1_p

      xyzgrid(1, 3) = X1(p)
      xyzgrid(2, 3) = Y2(p)
      xyzgrid(3, 3) = Z1_p

      xyzgrid(1, 4) = X2(p)
      xyzgrid(2, 4) = Y2(p)
      xyzgrid(3, 4) = Z1_p

      ! z = 2
      xyzgrid(1, 5) = X1(p)
      xyzgrid(2, 5) = Y1(p)
      xyzgrid(3, 5) = Z2_p

      xyzgrid(1, 6) = X2(p)
      xyzgrid(2, 6) = Y1(p)
      xyzgrid(3, 6) = Z2_p

      xyzgrid(1, 7) = X1(p)
      xyzgrid(2, 7) = Y2(p)
      xyzgrid(3, 7) = Z2_p

      xyzgrid(1, 8) = X2(p)
      xyzgrid(2, 8) = Y2(p)
      xyzgrid(3, 8) = Z2_p

      ! Store the values.
      xyzgrid_all(:, :, j) = real(xyzgrid, 4)

      cell_data(:, j) = real(val(p, :), 4)

      ! Units conversion.
      cell_data(:, j) = cell_data(:, j) / real(units_mult, 4)

      if (ncomponents == 3) then
        ! Flip the Z-axis of a vector.
        cell_data(3, j) = real(z_sign, 4) * cell_data(3, j)
      endif
    endif
  enddo

  ! Write the grid to a file.
  write(333) xyzgrid_all

  !-------------------------------------------------------------------
  ! Generate elements.
  !-------------------------------------------------------------------
  do p = 1, nelements_slice
    cell_indexes(1, p) = 8
    do i = 1, 8
      cell_indexes(i + 1, p) = 8 * (p - 1) + (i - 1)
    enddo
  enddo

  write(str1(1:8),'(i8)') nelements_slice
  write(str2(1:8),'(i8)') (8 + 1) * nelements_slice
  write(333) lf//lf//'CELLS '//str1//' '//str2//lf

  write(333) cell_indexes

  !-------------------------------------------------------------------
  ! Write cell types.
  !-------------------------------------------------------------------
  write(str1(1:8),'(i8)') nelements_slice
  write(333) lf//lf//'CELL_TYPES '//str1//lf

  ! VTK_VOXEL = 11
  cell_type = 11

  write(333) cell_type

  !-------------------------------------------------------------------
  ! Generate element data values.
  !-------------------------------------------------------------------
  write(str1(1:8),'(i8)') nelements_slice

  write(333) lf//lf//'CELL_DATA '//str1//lf

  if (ncomponents == 1) then
  ! Scalar data.
    write(333) 'SCALARS F FLOAT'//lf
    write(333) 'LOOKUP_TABLE default'//lf

  else if (ncomponents == 3) then
  ! Vector data.
    write(333) 'VECTORS vectors FLOAT'//lf
  endif

  write(333) cell_data

  close(333)

  deallocate(xyzgrid_all)
  deallocate(cell_data)
  deallocate(cell_indexes)
  deallocate(cell_type)

end subroutine visualisation_paraview_legogrid

!=======================================================================================================================
! Writes the data points in binary VTK format for Paraview visualization.
!=======================================================================================================================
subroutine visualisation_paraview_points(filename, myrank, ndata, ncomponents, val, X, Y, Z, INVERT_Z_AXIS, units_mult)
  integer, intent(in) :: myrank
  integer, intent(in) :: ndata, ncomponents
  real(kind=CUSTOM_REAL), intent(in) :: val(ncomponents, ndata)
  real(kind=CUSTOM_REAL), intent(in) :: X(ndata), Y(ndata), Z(ndata)
  logical, intent(in) :: INVERT_Z_AXIS
  ! Units multiplier.
  real(kind=CUSTOM_REAL), intent(in) :: units_mult
  ! Output file name.
  character(len=*), intent(in) :: filename

  character(len=256) :: filename_full
  character(len=256) :: msg
  integer :: ierr
  integer :: p

  real(kind=4), allocatable :: xyzgrid(:, :)
  real(kind=4), allocatable :: point_data(:, :)
  integer, allocatable :: cell_indexes(:, :)
  integer, allocatable :: cell_type(:)

  character :: lf*1, str1*8, str2*8
  lf = char(10) ! line feed character

  ! Create a directory.
  call execute_command_line('mkdir -p "'//trim(path_output)//'/Paraview"')

  filename_full = trim(path_output)//"/Paraview/"//filename

  open(333, file=filename_full, status='replace', access='stream', form='unformatted', action='write', &
       iostat=ierr, iomsg=msg)

  if (ierr /= 0) call exit_MPI("Error with writing the VTK file! path="&
                               //trim(filename_full)//", iomsg="//msg, myrank, ierr)

  !-----------------------------------------------------------------
  ! Write a header.
  !-----------------------------------------------------------------
  write(333) '# vtk DataFile Version 3.0'//lf
  write(333) 'TOMOFAST-X'//lf
  write(333) 'BINARY'//lf
  write(333) 'DATASET UNSTRUCTURED_GRID'//lf//lf

  write(str1(1:8),'(i8)') ndata
  write(333) 'POINTS '//str1//' FLOAT'//lf

  !-----------------------------------------------------------------
  ! Allocate memory.
  !-----------------------------------------------------------------
  allocate(xyzgrid(3, ndata), stat=ierr)
  allocate(point_data(ncomponents, ndata), stat=ierr)
  allocate(cell_indexes(2, ndata), stat=ierr)
  allocate(cell_type(ndata), stat=ierr)

  !-----------------------------------------------------------------
  ! Build lego-grid.
  !-----------------------------------------------------------------
  xyzgrid(1, :) = real(X, 4)
  xyzgrid(2, :) = real(Y, 4)
  xyzgrid(3, :) = real(Z, 4)

  if (INVERT_Z_AXIS) then
    xyzgrid(3, :) = - xyzgrid(3, :)
  endif

  ! Write the grid to a file.
  write(333) xyzgrid

  !-------------------------------------------------------------------
  ! Generate elements.
  !-------------------------------------------------------------------
  do p = 1, ndata
    cell_indexes(1, p) = 1
    cell_indexes(2, p) = p - 1
  enddo

  write(str1(1:8),'(i8)') ndata
  write(str2(1:8),'(i8)') 2 * ndata
  write(333) lf//lf//'CELLS '//str1//' '//str2//lf

  write(333) cell_indexes

  !-------------------------------------------------------------------
  ! Write cell types.
  !-------------------------------------------------------------------
  write(str1(1:8),'(i8)') ndata
  write(333) lf//lf//'CELL_TYPES '//str1//lf

  ! VTK_VERTEX = 1
  cell_type = 1

  write(333) cell_type

  !-------------------------------------------------------------------
  ! Generate element data values.
  !-------------------------------------------------------------------
  write(str1(1:8),'(i8)') ndata

  write(333) lf//lf//'POINT_DATA '//str1//lf

  if (ncomponents == 1) then
  ! Scalar data.
    write(333) 'SCALARS F FLOAT'//lf
    write(333) 'LOOKUP_TABLE default'//lf

  else if (ncomponents == 3) then
  ! Vector data.
    write(333) 'VECTORS vectors FLOAT'//lf
  endif

  point_data = real(val, 4)

  ! Convert data units.
  point_data = point_data / real(units_mult, 4)

  if (ncomponents == 3) then
    if (INVERT_Z_AXIS) then
      point_data(3, :) = - point_data(3, :)
    endif
  endif

  write(333) point_data

  close(333)

  deallocate(xyzgrid)
  deallocate(point_data)
  deallocate(cell_indexes)
  deallocate(cell_type)

end subroutine visualisation_paraview_points

end module paraview
