
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

  use global_typedefs, only: CUSTOM_REAL, path_output
  use mpi_tools, only: exit_MPI

  implicit none

  private

  public :: visualisation_paraview_legogrid
  public :: visualisation_paraview_struct_grid

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

!============================================================================================================
! This subroutine writes the file in (legacy) VTK format used for Paraview visualization.
! A slice of the full model can be specified via i1, i2, j1, j2, k1, k2 grid indexes.
! The Paraview structured grid datastructure is used, which allows additional analysis in Paraview (like contours).
!============================================================================================================
subroutine visualisation_paraview_struct_grid(filename, myrank, nelements, val, X1, Y1, Z1, X2, Y2, Z2, &
                                           i_index, j_index, k_index, &
                                           i1, i2, j1, j2, k1, k2, &
                                           INVERT_Z_AXIS)
  ! MPI rank of this process.
  integer, intent(in) :: myrank
  ! Total number of cells.
  integer, intent(in) :: nelements
  ! Values for visualization.
  real(kind=CUSTOM_REAL), intent(in) :: val(nelements)
  ! Coordinates of points in the grid.
  real(kind=CUSTOM_REAL), intent(in) :: X1(nelements), Y1(nelements), Z1(nelements)
  real(kind=CUSTOM_REAL), intent(in) :: X2(nelements), Y2(nelements), Z2(nelements)
  integer, intent(in) :: i_index(nelements), j_index(nelements), k_index(nelements)
  integer, intent(in) :: i1, i2, j1, j2, k1, k2
  logical, intent(in) :: INVERT_Z_AXIS
  ! Output file name.
  character(len=*), intent(in) :: filename

  character(len=256) :: filename_full
  character(len=256) :: msg
  ! I/O error code.
  integer :: ierr
  integer :: npoints, nelements_slice
  integer :: j, p

  real(kind=4), allocatable :: cell_centers(:, :)
  real(kind=4), allocatable :: cell_data(:)

  character :: lf*1, str1*8, str2*8, str3*8
  ! Line feed character.
  lf = char(10)

  ! Create a directory.
  call execute_command_line('mkdir -p '//trim(path_output)//"/Paraview/")

  filename_full = trim(path_output)//"/Paraview/"//filename

  open(333, file=filename_full, status='replace', access='stream', form='unformatted', action='write', &
       iostat=ierr, iomsg=msg)

  if (ierr /= 0) call exit_MPI("Error with writing the VTK file! path="&
                               //trim(filename_full)//", iomsg="//msg, myrank, ierr)

  ! ************* Generate points ******************

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

  ! Allocate memory.
  allocate(cell_centers(3, nelements_slice), stat=ierr)
  allocate(cell_data(nelements_slice), stat=ierr)

  !-------------------------------------------------------------------
  ! Build the grid.
  !-------------------------------------------------------------------
  j = 0
  do p = 1, nelements
    if (index_included(p, i_index, j_index, k_index, i1, i2, j1, j2, k1, k2)) then
      j = j + 1

      cell_centers(1, j) = real(0.5 * (X1(p) + X2(p)))
      cell_centers(2, j) = real(0.5 * (Y1(p) + Y2(p)))
      cell_centers(3, j) = real(0.5 * (Z1(p) + Z2(p)))

      cell_data(j) = real(val(p), 4)
    endif
  enddo

  if (INVERT_Z_AXIS) then
    cell_centers(3, :) = - cell_centers(3, :)
  endif

  ! Write the grid to a file.
  write(333) cell_centers

  !-------------------------------------------------------------------
  ! Generate element data values.
  !-------------------------------------------------------------------
  write(str1(1:8),'(i8)') nelements_slice

  write(333) lf//lf//'POINT_DATA '//str1//lf
  write(333) 'SCALARS F FLOAT'//lf
  write(333) 'LOOKUP_TABLE default'//lf

  write(333) cell_data

  close(333)

  deallocate(cell_centers)
  deallocate(cell_data)

end subroutine visualisation_paraview_struct_grid

!============================================================================================================
! This subroutine writes the file in (legacy) VTK format used for Paraview visualization.
! A slice of the full model can be specified via i1, i2, j1, j2, k1, k2 grid indexes.
! An arbitrary "lego" grid is considered.
!
! VO: Converted the ASCII version to binary using logics of the code of Renato N. Elias:
!     http://www.nacad.ufrj.br/~rnelias/paraview/VTKFormats.f90
!============================================================================================================
subroutine visualisation_paraview_legogrid(filename, myrank, nelements, val, X1, Y1, Z1, X2, Y2, Z2, &
                                           i_index, j_index, k_index, &
                                           i1, i2, j1, j2, k1, k2, &
                                           INVERT_Z_AXIS)
  ! MPI rank of this process.
  integer, intent(in) :: myrank
  ! Total number of cells.
  integer, intent(in) :: nelements
  ! Values for visualization.
  real(kind=CUSTOM_REAL), intent(in) :: val(nelements)
  ! Coordinates of points in the grid.
  real(kind=CUSTOM_REAL), intent(in) :: X1(nelements), Y1(nelements), Z1(nelements)
  real(kind=CUSTOM_REAL), intent(in) :: X2(nelements), Y2(nelements), Z2(nelements)
  integer, intent(in) :: i_index(nelements), j_index(nelements), k_index(nelements)
  integer, intent(in) :: i1, i2, j1, j2, k1, k2
  logical, intent(in) :: INVERT_Z_AXIS
  ! Output file name.
  character(len=*), intent(in) :: filename

  character(len=256) :: filename_full
  character(len=256) :: msg
  ! I/O error code.
  integer :: ierr
  integer :: npoints, nelements_slice
  integer :: i, j, p
  real(kind=CUSTOM_REAL) :: xyzgrid(3, 8)
  real(kind=CUSTOM_REAL) :: Z1_p, Z2_p

  real(kind=4), allocatable :: xyzgrid_all(:, :, :)
  real(kind=4), allocatable :: cell_data(:)
  integer, allocatable :: cell_indexes(:, :)
  integer, allocatable :: cell_type(:)

  character :: lf*1, str1*8, str2*8
  lf = char(10) ! line feed character

  ! Create a directory.
  call execute_command_line('mkdir -p '//trim(path_output)//"/Paraview/")

  filename_full = trim(path_output)//"/Paraview/"//filename

  open(333, file=filename_full, status='replace', access='stream', form='unformatted', action='write', &
       iostat=ierr, iomsg=msg)

  if (ierr /= 0) call exit_MPI("Error with writing the VTK file! path="&
                               //trim(filename_full)//", iomsg="//msg, myrank, ierr)

  ! ************* Generate points ******************

  write(333) '# vtk DataFile Version 3.0'//lf
  write(333) 'TOMOFAST-X'//lf
  write(333) 'BINARY'//lf
  write(333) 'DATASET UNSTRUCTURED_GRID'//lf//lf

  ! Number of elements in the requested slice.
  nelements_slice = (i2 - i1 + 1) * (j2 - j1 + 1) * (k2 - k1 + 1)

  npoints = 8 * nelements_slice

  write(str1(1:8),'(i8)') npoints
  write(333) 'POINTS '//str1//' FLOAT'//lf

  !=================================================================
  ! Allocate memory.
  !=================================================================
  allocate(xyzgrid_all(3, 8, nelements_slice), stat=ierr)
  allocate(cell_data(nelements_slice), stat=ierr)
  allocate(cell_indexes(9, nelements_slice), stat=ierr)
  allocate(cell_type(nelements_slice), stat=ierr)

  !====================================
  ! Build lego-grid.
  !====================================
  j = 0
  do p = 1, nelements
    if (index_included(p, i_index, j_index, k_index, i1, i2, j1, j2, k1, k2)) then
      j = j + 1

      if (INVERT_Z_AXIS) then
        Z1_p = - Z1(p)
        Z2_p = - Z2(p)
      else
        Z1_p = Z1(p)
        Z2_p = Z2(p)
      endif

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

      cell_data(j) = real(val(p), 4)
    endif
  enddo

  ! Write the grid to a file.
  write(333) xyzgrid_all

  ! ************* Generate elements ******************

  ! See documentation here http://dunne.uni-hd.de/VisuSimple/documents/vtkfileformat.html

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

  ! Write cell types --------------------------------------------
  write(str1(1:8),'(i8)') nelements_slice
  write(333) lf//lf//'CELL_TYPES '//str1//lf

  ! VTK_VOXEL = 11
  cell_type = 11

  write(333) cell_type

  ! ************* Generate element data values ******************

  write(str1(1:8),'(i8)') nelements_slice

  write(333) lf//lf//'CELL_DATA '//str1//lf
  write(333) 'SCALARS F FLOAT'//lf
  write(333) 'LOOKUP_TABLE default'//lf

  write(333) cell_data

  close(333)

  deallocate(xyzgrid_all)
  deallocate(cell_data)
  deallocate(cell_indexes)
  deallocate(cell_type)

end subroutine visualisation_paraview_legogrid

end module paraview
