
!========================================================================
!
!                    T O M O F A S T X  Version 1.0
!                  ----------------------------------
!
!              Main authors: Vitaliy Ogarko, Roland Martin,
!                   Jeremie Giraud, Dimitri Komatitsch.
! CNRS, France, and University of Western Australia.
! (c) CNRS, France, and University of Western Australia. January 2018
!
! This software is a computer program whose purpose is to perform
! capacitance, gravity, magnetic, or joint gravity and magnetic tomography.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

!========================================================================================
! Routines for exporting the model data to the vtk format for visualization in Paraview.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!========================================================================================

module paraview

  use global_typedefs, only: CUSTOM_REAL, path_output
  use mpi_tools, only: exit_MPI

  implicit none

  private

  public :: visualisation_paraview
  public :: visualisation_paraview_legogrid
  public :: visualisation_qgis

contains

!=================================================================================================================
! TODO: Move to another module.
! This subroutine writes model snapshots for visualization in QGIS.
! The file can bu opened in QGIS using "Add Raster Layer" option.
! Tested using QGIS 2.12.3-Lyon.
!=================================================================================================================
subroutine visualisation_qgis(filename, myrank, val, ncols, nrows, cell_size)
  ! MPI rank of this process.
  integer, intent(in) :: myrank
  ! Values for visualization (2D slice).
  real(kind=CUSTOM_REAL), intent(in) :: val(:, :)
  integer, intent(in) :: ncols, nrows
  ! The cell size.
  real(kind=CUSTOM_REAL), intent(in) :: cell_size
  ! Output file name.
  character(len=*), intent(in) :: filename

  ! I/O error code.
  integer :: ierr
  character(len=256) :: filename_full
  character(len=256) :: msg
  integer :: i, j

  ! ************** Create data file **********************
  ! TODO: move this part to a function (and use it also in visualisation_paraview())
  call system('mkdir -p '//trim(path_output)//"/QGIS/")

  filename_full = trim(path_output)//"/QGIS/"//filename

  open(unit=333, file=filename_full, status='unknown', action='write', iostat=ierr, iomsg=msg)

  if (ierr /= 0) call exit_MPI("Error with writing the QGIS file! path="&
                               //filename_full//" iomsg="//msg, myrank, ierr)

  ! ************** Write data into file **********************

  write (333, '(''ncols '',i9)') ncols
  write (333, '(''nrows '',i9)') nrows
  write (333, '(''xllcenter '',f16.8)') 0.d0
  write (333, '(''yllcenter '',f16.8)') 0.d0
  write (333, '(''cellsize '',f16.8)') cell_size

  do j = 1, nrows
    do i = 1, ncols
      write(333, '(f16.8)', advance='no') val(i, j)
    enddo
    write(333, *)
  enddo

end subroutine visualisation_qgis

!============================================================================================================
! This subroutine writes the file in (legacy) VTK format used for Paraview visualization.
! An arbitrary lego grid is considered.
!============================================================================================================
subroutine visualisation_paraview_legogrid(filename, myrank, nelements, val, xgrid, ygrid, zgrid)
  ! MPI rank of this process.
  integer, intent(in) :: myrank
  ! Total number of cells.
  integer, intent(in) :: nelements
  ! Values for visualization.
  real(kind=CUSTOM_REAL), intent(in) :: val(1:nelements)
  ! Coordinates of points in the grid.
  real(kind=CUSTOM_REAL), intent(in) :: xgrid(1:(8 * nelements))
  real(kind=CUSTOM_REAL), intent(in) :: ygrid(1:(8 * nelements))
  real(kind=CUSTOM_REAL), intent(in) :: zgrid(1:(8 * nelements))
  ! Output file name.
  character(len=*), intent(in) :: filename

  character(len=256) :: filename_full
  character(len=256) :: msg
  ! I/O error code.
  integer :: ierr
  integer :: npoints
  integer :: i, ind

  ! (+) Copied from visualisation_paraview -----------------------------------
  call system('mkdir -p '//trim(path_output)//"/Paraview/")

  filename_full = trim(path_output)//"/Paraview/"//filename

  open(unit=333, file=filename_full, status='unknown', action='write', iostat=ierr, iomsg=msg)

  if (ierr /= 0) call exit_MPI("Error with writing the VTK file! path="&
                               //filename_full//" iomsg="//msg, myrank, ierr)
  ! (-) ----------------------------------------------------------------------

  ! ************* generate points ******************

  write (333, '(''# vtk DataFile Version 3.1'')')
  write (333, '(''TOMOFAST3D'')')
  write (333, '(''ASCII'')')

  write (333, '(''DATASET UNSTRUCTURED_GRID'')')

  npoints = 8 * nelements

  write (333, '(''POINTS '',i9,'' FLOAT'')') npoints

  do i = 1, npoints
    write (333, *) xgrid(i), ygrid(i), zgrid(i)
  enddo

  write (333, *)

  ! ************* generate elements ******************

  ! See documentation here http://dunne.uni-hd.de/VisuSimple/documents/vtkfileformat.html
  write (333, '(''CELLS '',2(i9,1x))') nelements, (8 + 1) * nelements

  do i = 1, nelements

    ind = 8 * (i - 1)

    ! Define a cell.
    write (333, "(i1,8(1x,i9))") 8, ind + 0, ind + 1, ind + 2, ind + 3, &
                                    ind + 4, ind + 5, ind + 6, ind + 7
  enddo

  write (333, *)

  write (333, '(''CELL_TYPES '',2(i9,1x))') nelements

  do i = 1, nelements
    if (i < nelements) then
      write (333, "(i2,1x)", advance="no") 11  ! VTK_VOXEL = 11
    else
      write (333, "(i2)") 11
    endif
  enddo

  write (333, *)

  ! ************* generate element data values ******************

  write (333,'(''CELL_DATA '',i9)') nelements
  write (333,'(''SCALARS F FLOAT'')')
  write (333,'(''LOOKUP_TABLE default'')')

  do i = 1, nelements
    write (333, *) val(i)
  enddo

  write (333, *)

  close(333)

end subroutine visualisation_paraview_legogrid

!============================================================================================================
! This subroutine writes the file in (legacy) VTK format used for Paraview visualization.
! data_type = 'POINT_DATA' or 'CELL_DATA', depending on whether we set the data
! at the grid vertex or at the cell center. In the 'POINT_DATA' case Paraview performs nice smoothing of the data.
!============================================================================================================
subroutine visualisation_paraview(filename, myrank, nx, ny, nz, val, xgrid, ygrid, zgrid, &
                                  imin, imax, jmin, jmax, kmin, kmax, &
                                  i_step, j_step, k_step, data_type)
  ! MPI rank of this process.
  integer, intent(in) :: myrank
  ! Dimensions of the problem.
  integer, intent(in) :: nx, ny, nz
  ! Values for visualization.
  real(kind=CUSTOM_REAL), intent(in) :: val(0:, 0:, 0:)
  ! Coordinates of points in the grid.
  real(kind=CUSTOM_REAL), intent(in) :: xgrid(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: ygrid(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: zgrid(0:, 0:, 0:)
  ! The beginning and ending radius and angle grid nodes.
  integer, intent(in) :: imin, imax, jmin, jmax, kmin, kmax
  ! The steps (for loops) in every dimension.
  integer, intent(in) :: i_step, j_step, k_step
  ! Output file name.
  character(len=*), intent(in) :: filename
  ! Type of the Paraview data.
  character(len=*), intent(in) :: data_type

  ! Loop variables.
  integer :: i, j, k, ii, jj, kk, counter
  ! I/O error code.
  integer :: ierr
  ! The grid dimensions for visualization.
  integer :: xdim, ydim, zdim
  ! The number of grid points, the number of volume elements.
  integer :: npoints, nelements
  ! The nodes of the 3D volume element.
  integer :: ival(8)
  character(len=256) :: filename_full
  character(len=256) :: msg
  !integer :: cmd_s

  ! The global index array.
  integer, allocatable :: ibool(:, :, :)

  if (data_type /= 'POINT_DATA' .and. data_type /= 'CELL_DATA') then
    print *, 'Unknown data type in visualisation_paraview()!'
    stop
  endif

  ! ************** create data file **********************
  ! TODO: Put this into a separate function in utils and also use for creating the path_output.
  ! If a folder Paraview does not exist, then it will be created.
  !call execute_command_line('mkdir -p '//trim(path_output)//"/Paraview/",CMDSTAT=cmd_s,CMDMSG=msg)

  ! TODO: ifort with debug flags complains about CMDSTAT and CMDMSG:
  ! "error #6632: Keyword arguments are invalid without an explicit interface".
  !call execute_command_line('mkdir -p '//trim(path_output)//"/Paraview/")

  ! TODO: In some cases on the cluster an error appears:
  ! TODO: "CMDMSG=Termination status of the command-language interpreter cannot be obtained"
  ! TODO: Find why this happens?
  !if (cmd_s /= 0) call exit_MPI("Error with writing the Paraview folder! "&
  !                             //" CMDMSG="//msg, myrank, cmd_s)

  ! ifort 14.0 does not support execute_command_line()
  call system('mkdir -p '//trim(path_output)//"/Paraview/")

  filename_full = trim(path_output)//"/Paraview/"//filename

  open(unit=333, file=filename_full, status='unknown', action='write', iostat=ierr, iomsg=msg)

  if (ierr /= 0) call exit_MPI("Error with writing the VTK file! path="&
                               //filename_full//" iomsg="//msg, myrank, ierr)

  ! ************* allocate memory ************************

  ! TODO: reserve only memory needed for visualization, not for the whole model.
  ! TODO: and remove nx, ny, nz from the interface.
  allocate(ibool(0:nx + 1, 0:ny + 1, 0:nz + 1))

  ! ***********************************************
  ! Grid dimensions.
  xdim = (imax - imin) / i_step + 1
  ydim = (jmax - jmin) / j_step + 1
  zdim = (kmax - kmin) / k_step + 1

  npoints = xdim * ydim * zdim
  nelements = (xdim - 1) * (ydim - 1) * (zdim - 1)

  !if (myrank == 0) print *, 'npoints =', npoints
  !if (myrank == 0) print *, 'nelements =', nelements
  !if (myrank == 0) print *, 'PARAVIEW xdim, ydim, zdim = ', xdim, ydim, zdim

  ! ************* generate points ******************

  write (333, '(''# vtk DataFile Version 3.1'')')
  write (333, '(''TOMOFAST3D'')')
  write (333, '(''ASCII'')')

  if (data_type == 'POINT_DATA') then
    write (333, '(''DATASET POLYDATA'')')
  else ! CELL_DATA
    write (333, '(''DATASET UNSTRUCTURED_GRID'')')
  endif

  write (333, '(''POINTS '',i9,'' FLOAT'')') npoints

  counter = 0
  kk = 0
  do k = kmin, kmax, k_step
    jj = 0
    do j = jmin, jmax, j_step
      ii = 0
      do i = imin, imax, i_step
        ! Write the grid points coordinates.
        write (333, *) xgrid(i, j, k), ygrid(i, j, k), zgrid(i, j, k)

        ! Store the point index.
        ibool(ii, jj, kk) = counter

        ! Start counting points from zero in the VTK format.
        counter = counter + 1
        ii = ii + 1
      enddo
      jj = jj + 1
    enddo
    kk = kk + 1
  enddo

  write (333, *)

  ! ************* generate elements ******************

  if (data_type == 'POINT_DATA') then

    write (333, '(''POLYGONS '',2(i9,1x))') 6 * nelements, 6 * (4 + 1) * nelements

    do k = 0, zdim - 2
      do j = 0, ydim - 2
        do i = 0, xdim - 2
          ! Define the corners of the 3D volume element.
          ival(1) = ibool(i, j, k)        ! 000
          ival(2) = ibool(i, j + 1, k)      ! 010
          ival(3) = ibool(i + 1, j + 1, k)    ! 110
          ival(4) = ibool(i + 1, j, k)      ! 100
          ival(5) = ibool(i + 1, j, k + 1)    ! 101
          ival(6) = ibool(i + 1, j + 1, k + 1)  ! 111
          ival(7) = ibool(i, j + 1, k + 1)    ! 011
          ival(8) = ibool(i, j, k + 1)      ! 001

          ! Define cells.
          write (333,"(i1,4(1x,i9))") 4, ival(1), ival(2), ival(3), ival(4) ! bottom
          write (333,"(i1,4(1x,i9))") 4, ival(5), ival(6), ival(7), ival(8) ! top
          write (333,"(i1,4(1x,i9))") 4, ival(1), ival(2), ival(7), ival(8) ! left
          write (333,"(i1,4(1x,i9))") 4, ival(4), ival(3), ival(6), ival(5) ! right
          write (333,"(i1,4(1x,i9))") 4, ival(1), ival(4), ival(5), ival(8) ! front
          write (333,"(i1,4(1x,i9))") 4, ival(2), ival(3), ival(6), ival(7) ! rear
        enddo
      enddo
    enddo

  else ! CELL_DATA

    ! See documentation here http://dunne.uni-hd.de/VisuSimple/documents/vtkfileformat.html
    write (333, '(''CELLS '',2(i9,1x))') nelements, (8 + 1) * nelements

    counter = 0
    do k = 0, zdim - 2
      do j = 0, ydim - 2
        do i = 0, xdim - 2
          ! Define the corners of the 3D volume element.
          ival(1) = ibool(i, j, k)
          ival(2) = ibool(i + 1, j, k)
          ival(3) = ibool(i, j + 1, k)
          ival(4) = ibool(i + 1, j + 1, k)

          ival(5) = ibool(i, j, k + 1)
          ival(6) = ibool(i + 1, j, k + 1)
          ival(7) = ibool(i, j + 1, k + 1)
          ival(8) = ibool(i + 1, j + 1, k + 1)

          ! Define a cell.
          write (333, "(i1,8(1x,i9))") 8, ival(1), ival(2), ival(3), ival(4), &
                                          ival(5), ival(6), ival(7), ival(8)

          counter = counter + 1
        enddo
      enddo
    enddo

    write (333, *)

    !print *, 'counter1 = ', counter

    write (333, '(''CELL_TYPES '',2(i9,1x))') nelements

    do i = 1, nelements
      if (i < nelements) then
        write (333, "(i2,1x)", advance="no") 11  ! VTK_VOXEL = 11
      else
        write (333, "(i2)") 11
      endif
    enddo

  endif

  write (333, *)

  ! ************* generate element data values ******************

  if (data_type == 'POINT_DATA') then

    write (333,'(''POINT_DATA '',i9)') npoints
    write (333,'(''SCALARS F float'')')
    write (333,'(''LOOKUP_TABLE default'')')

    do k = kmin, kmax, k_step
      do j = jmin, jmax, j_step
        do i = imin, imax, i_step
          write (333, *) val(i, j, k)
        enddo
      enddo
    enddo

  else ! CELL_DATA

    write (333,'(''CELL_DATA '',i9)') nelements
    write (333,'(''SCALARS F FLOAT'')')
    write (333,'(''LOOKUP_TABLE default'')')

    counter = 0
    do k = kmin, kmax - 1, k_step
      do j = jmin, jmax - 1, j_step
        do i = imin, imax - 1, i_step
          write (333, *) val(i, j, k)

          counter = counter + 1
        enddo
      enddo
    enddo

  endif

  write (333, *)

  !print *, 'counter2 = ', counter

  close(333)

  deallocate(ibool)

end subroutine visualisation_paraview

end module paraview
