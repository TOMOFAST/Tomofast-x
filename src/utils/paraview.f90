
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

  public :: visualisation_paraview
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
  real(kind=CUSTOM_REAL), intent(in) :: val(:)
  ! Coordinates of points in the grid.
  real(kind=CUSTOM_REAL), intent(in) :: X1(:), Y1(:), Z1(:)
  real(kind=CUSTOM_REAL), intent(in) :: X2(:), Y2(:), Z2(:)
  integer, intent(in) :: i_index(:), j_index(:), k_index(:)
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
  real(kind=CUSTOM_REAL), intent(in) :: val(:)
  ! Coordinates of points in the grid.
  real(kind=CUSTOM_REAL), intent(in) :: X1(:), Y1(:), Z1(:)
  real(kind=CUSTOM_REAL), intent(in) :: X2(:), Y2(:), Z2(:)
  integer, intent(in) :: i_index(:), j_index(:), k_index(:)
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
  real(kind=CUSTOM_REAL) :: xgrid(8)
  real(kind=CUSTOM_REAL) :: ygrid(8)
  real(kind=CUSTOM_REAL) :: zgrid(8)
  integer :: cell_type

  real(kind=4), allocatable :: xgrid_all(:, :)
  real(kind=4), allocatable :: ygrid_all(:, :)
  real(kind=4), allocatable :: zgrid_all(:, :)
  real(kind=4), allocatable :: cell_data(:)

  character :: lf*1, str1*8, str2*8
  lf = char(10) ! line feed character

  ! Create a directory.
  call execute_command_line('mkdir -p '//trim(path_output)//"/Paraview/")

  filename_full = trim(path_output)//"/Paraview/"//filename

  open(unit=333, file=filename_full, status='replace', access='stream', form='unformatted', &
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
  allocate(xgrid_all(8, nelements_slice), stat=ierr)
  allocate(ygrid_all(8, nelements_slice), stat=ierr)
  allocate(zgrid_all(8, nelements_slice), stat=ierr)
  allocate(cell_data(nelements_slice), stat=ierr)

  !====================================
  ! Build lego-grid.
  !====================================
  j = 0
  do p = 1, nelements
    if (index_included(p, i_index, j_index, k_index, i1, i2, j1, j2, k1, k2)) then
      j = j + 1

      ! z = 1
      xgrid(1) = X1(p)
      ygrid(1) = Y1(p)
      zgrid(1) = Z1(p)

      xgrid(2) = X2(p)
      ygrid(2) = Y1(p)
      zgrid(2) = Z1(p)

      xgrid(3) = X1(p)
      ygrid(3) = Y2(p)
      zgrid(3) = Z1(p)

      xgrid(4) = X2(p)
      ygrid(4) = Y2(p)
      zgrid(4) = Z1(p)

      ! z = 2
      xgrid(5) = X1(p)
      ygrid(5) = Y1(p)
      zgrid(5) = Z2(p)

      xgrid(6) = X2(p)
      ygrid(6) = Y1(p)
      zgrid(6) = Z2(p)

      xgrid(7) = X1(p)
      ygrid(7) = Y2(p)
      zgrid(7) = Z2(p)

      xgrid(8) = X2(p)
      ygrid(8) = Y2(p)
      zgrid(8) = Z2(p)

      if (INVERT_Z_AXIS) then
        zgrid = -zgrid
      endif

      ! Store the values.
      xgrid_all(:, j) = real(xgrid, 4)
      ygrid_all(:, j) = real(ygrid, 4)
      zgrid_all(:, j) = real(zgrid, 4)

      cell_data(j) = real(val(p), 4)

    endif
  enddo

  ! Write the grid to a file.
  write(333) ((xgrid_all(i, j), ygrid_all(i, j), zgrid_all(i, j), i = 1, 8), j = 1, nelements_slice)

  ! ************* Generate elements ******************

  ! See documentation here http://dunne.uni-hd.de/VisuSimple/documents/vtkfileformat.html

  write(str1(1:8),'(i8)') nelements_slice
  write(str2(1:8),'(i8)') (8 + 1) * nelements_slice
  write(333) lf//lf//'CELLS '//str1//' '//str2//lf

  write(333) (8, 8 * (p - 1) + 0, 8 * (p - 1) + 1, 8 * (p - 1) + 2, 8 * (p - 1) + 3, &
                 8 * (p - 1) + 4, 8 * (p - 1) + 5, 8 * (p - 1) + 6, 8 * (p - 1) + 7, p = 1, nelements_slice)


  write(str1(1:8),'(i8)') nelements_slice
  write(333) lf//lf//'CELL_TYPES '//str1//lf

  ! VTK_VOXEL = 11
  cell_type = 11

  write(333) (cell_type, p = 1, nelements_slice)

  ! ************* Generate element data values ******************

  write(str1(1:8),'(i8)') nelements_slice

  write(333) lf//lf//'CELL_DATA '//str1//lf
  write(333) 'SCALARS F FLOAT'//lf
  write(333) 'LOOKUP_TABLE default'//lf

  write(333) (cell_data(p), p = 1, nelements_slice)

  close(333)

  deallocate(xgrid_all)
  deallocate(ygrid_all)
  deallocate(zgrid_all)
  deallocate(cell_data)

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

  ! The global index array.
  integer, allocatable :: ibool(:, :, :)

  if (data_type /= 'POINT_DATA' .and. data_type /= 'CELL_DATA') then
    print *, 'Unknown data type in visualisation_paraview()!'
    stop
  endif

  ! ************** Create data file **********************
  ! If a folder Paraview does not exist, then it will be created.
  call execute_command_line('mkdir -p '//trim(path_output)//"/Paraview/")

  filename_full = trim(path_output)//"/Paraview/"//filename

  open(unit=333, file=filename_full, status='replace', action='write', iostat=ierr, iomsg=msg)

  if (ierr /= 0) call exit_MPI("Error with writing the VTK file! path="&
                               //filename_full//" iomsg="//msg, myrank, ierr)

  ! ************* Allocate memory ************************

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

  ! ************* Generate points ******************

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

  ! ************* Generate elements ******************

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

  ! ************* Generate element data values ******************

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

  close(333)

  deallocate(ibool)

end subroutine visualisation_paraview

end module paraview
