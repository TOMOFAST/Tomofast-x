
!========================================================================
!
!                          T o m o f a s t - x
!                        -----------------------
!
!              Main authors: Vitaliy Ogarko, Roland Martin,
!                   Jeremie Giraud, Dimitri Komatitsch.
!
! (c) CNRS, France, and University of Western Australia. January 2018
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

module paraview_ect

  use global_typedefs, only: CUSTOM_REAL, path_output
  use mpi_tools, only: exit_MPI

  implicit none

  private

  public :: paraview_write_2d_profiles
  public :: paraview_write_3d_profiles

  private :: visualisation_paraview

contains

!============================================================================================================
! This subroutine writes the 2d horizontal profiles (cuts at different hight Z) for Paraview visualization.
!============================================================================================================
subroutine paraview_write_2d_profiles(myrank, name_prefix, ielectrode, nr, ntheta, nz, nzlocal, &
                                      sol, xgrid, ygrid, zgrid)
  ! MPI rank of this process.
  integer, intent(in) :: myrank
  ! Prefix for the file name.
  character(len=*), intent(in) :: name_prefix
  ! Current electrode.
  integer, intent(in) :: ielectrode
  ! Dimensions of the problem (r, theta and z direction, the model is split in the z direction for MPI).
  integer, intent(in) :: nr, ntheta, nz, nzlocal
  ! Solution (e.g. electric potential or permittivity).
  real(kind=CUSTOM_REAL), intent(in) :: sol(0:, 0:, 0:)
  ! Coordinates of points in the grid.
  real(kind=CUSTOM_REAL), intent(in) :: xgrid(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: ygrid(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: zgrid(0:, 0:, 0:)

  character(len=60) :: name
  integer :: k, k1

  do k=1,nzlocal
    k1 = myrank*nzlocal+k

    !if (k1 == 0 .or. k1 == nz/4 .or. k1 == nz/2 .or. k1 == 3*nz/4 .or. k1 == nz) then
    if (k1 == nz/4 .or. k1 == nz/2 .or. k1 == 3*nz/4 .or. &
        k1 == 3*nz/8 .or. k1 == 5*nz/8) then ! Added these heights for two electrode rings.

      ! Make a filename.
      write(name,'("nz",i3.3,"_k",i3.3,"_e",i2.2,"_r",i2.2,".vtk")') nz, k1, ielectrode, myrank

      name = name_prefix//name

      call visualisation_paraview(name, myrank, nr, ntheta, nzlocal, sol, xgrid, ygrid, zgrid, &
                                  0, nr+1, 0, ntheta, k, k+1, 1, 1, 1, 'POINT_DATA')
    endif
  enddo

end subroutine paraview_write_2d_profiles

!================================================================================================
! This subroutine writes the 3d horizontal profiles for Paraview visualization.
!================================================================================================
subroutine paraview_write_3d_profiles(myrank, name_prefix, ielectrode, nr, ntheta, nzlocal, &
                                      sol, xgrid, ygrid, zgrid)
  ! MPI rank of this process.
  integer, intent(in) :: myrank
  ! Prefix for the file name.
  character(len=*), intent(in) :: name_prefix
  ! Current electrode.
  integer, intent(in) :: ielectrode
  ! Dimensions of the problem (r, theta and z direction, the model is split in the z direction for MPI).
  integer, intent(in) :: nr, ntheta, nzlocal
  ! Solution (e.g. electric potential or permittivity).
  real(kind=CUSTOM_REAL), intent(in) :: sol(0:, 0:, 0:)
  ! Coordinates of points in the grid.
  real(kind=CUSTOM_REAL), intent(in) :: xgrid(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: ygrid(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: zgrid(0:, 0:, 0:)

  character(len=60) :: name
  integer :: kmin, kmax

  ! To visualize the model.
  kmin = 1
  kmax = nzlocal + 1

  ! Visualize a vertical profile of the cylinder.
  write(name, '("half_nzlocal",i3.3,"_e",i2.2,"_r",i2.2,".vtk")') nzlocal, ielectrode, myrank

  name = name_prefix//name

  call visualisation_paraview(name, myrank, nr, ntheta, nzlocal, &
                              sol, xgrid, ygrid, zgrid(:, :, (myrank * nzlocal + kmin - 1):(myrank * nzlocal + kmax)), &
                              0, nr + 1, ntheta / 2 + 1, ntheta + 1, kmin, kmax, 1, ntheta / 4, 1, 'POINT_DATA')

end subroutine paraview_write_3d_profiles

!============================================================================================================
! LEGACY function used in the ECT problem.
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
          ival(1) = ibool(i, j, k)              ! 000
          ival(2) = ibool(i, j + 1, k)          ! 010
          ival(3) = ibool(i + 1, j + 1, k)      ! 110
          ival(4) = ibool(i + 1, j, k)          ! 100
          ival(5) = ibool(i + 1, j, k + 1)      ! 101
          ival(6) = ibool(i + 1, j + 1, k + 1)  ! 111
          ival(7) = ibool(i, j + 1, k + 1)      ! 011
          ival(8) = ibool(i, j, k + 1)          ! 001

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

end module paraview_ect
