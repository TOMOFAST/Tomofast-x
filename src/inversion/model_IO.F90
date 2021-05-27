
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

!================================================================================================
! This class contains ONLY I/O operations (reading/writing) with the model.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015-2017.
!================================================================================================
module model_IO

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use paraview
  use grid
  use string, only: str
  use parallel_tools
  use string
  use model_base

  implicit none

  private

  type, extends(t_model_base), public :: t_model_IO


  contains
    private

    procedure, public, pass :: read => model_read
    procedure, public, pass :: write => model_write
    procedure, public, pass :: read_bound_constraints => read_bound_constraints

    procedure, private, pass :: read_voxels_format => model_read_voxels_format
    procedure, private, pass :: write_voxels_format => model_write_voxels_format
    procedure, private, pass :: write_qgis => model_write_qgis
    procedure, private, pass :: write_paraview => model_write_paraview

  end type t_model_IO

contains

!==========================================================================================================
! Read the full model and grid, and then broadcast to all CPUs.
!==========================================================================================================
subroutine model_read(this, file_name, read_grid, myrank)
  class(t_model_IO), intent(inout) :: this
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank
  logical, intent(in) :: read_grid

  integer :: ierr

  if (myrank == 0) then
    call this%read_voxels_format(file_name, read_grid, myrank)
  endif

  if (read_grid) then
    ! Broadcast full grid to all CPUs.
    call this%grid_full%broadcast(myrank)
  endif

  ! Broadcast full model to all CPUs.
  call MPI_Bcast(this%val_full, this%nelements_total, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("Error in MPI_Bcast in model_read!", myrank, ierr)

  if (myrank == 0) then
    print *, 'Xmin, Xmax =', this%get_Xmin(), this%get_Xmax()
    print *, 'Ymin, Ymax =', this%get_Ymin(), this%get_Ymax()
  endif

end subroutine model_read

!==========================================================================================================
! Read the local bound constraints.
!==========================================================================================================
subroutine read_bound_constraints(this, file_name, myrank, nbproc)
  class(t_model_IO), intent(inout) :: this
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank, nbproc
  
  integer :: ierr, nsmaller, ind, i, nelements_read, nlithos_read
  character(len=256) :: msg
  type(t_parallel_tools) :: pt
  character(len=200) :: dummy_line

  if (myrank == 0) print *, 'Reading local bound constraints from file ', trim(file_name)

  open(10, file=trim(file_name), status='old', action='read', iostat=ierr, iomsg=msg)
  if (ierr /= 0) call exit_MPI("Error in opening the bound constraints file! path=" &
                 //file_name//" iomsg="//msg, myrank, ierr)

  read(10, *, iostat=ierr) nelements_read, nlithos_read
  if (ierr /= 0) call exit_MPI("Problem while reading the bound constraints file!", myrank, ierr)

  ! Sanity check.
  if (this%nelements_total /= nelements_read) &
    call exit_MPI("The constraints are not correctly defined!"//new_line('a') &
          //"nelements="//str(this%nelements)//new_line('a') &
          //"nelements_read="//str(nelements_read)//new_line('a') &
          //"nelements_total="//str(this%nelements_total), myrank, 0)

  if (this%nlithos /= nlithos_read) &
    call exit_MPI("The constraints are not correctly defined!"//new_line('a') &
          //"nlithos="//str(this%nlithos)//new_line('a') &
          //"nlithos_read="//str(nlithos_read), myrank, 0)

  if (myrank == 0) print *, 'Read nelements, nlithos = ', nelements_read, nlithos_read

  ! The number of elements on CPUs with rank smaller than myrank.
  nsmaller = pt%get_nsmaller(this%nelements, myrank, nbproc)

  ! Reading.
  do i = 1, this%nelements_total
    if (i > nsmaller .and. i <= nsmaller + this%nelements) then
      ind = i - nsmaller
      read(10, *, iostat=ierr) this%min_local_bound(ind, :), this%max_local_bound(ind, :), this%local_bound_constraints_weight(ind)
    else
      read(10, '(A)', iostat=ierr) dummy_line
    endif

    if (i > nsmaller + this%nelements) then
      exit
    endif
  enddo

  close(10)

end subroutine read_bound_constraints

!================================================================================================
! Read the full model and grid in voxels format.
!================================================================================================
subroutine model_read_voxels_format(this, file_name, read_grid, myrank)
  class(t_model_IO), intent(inout) :: this
  character(len=*), intent(in) :: file_name
  logical, intent(in) :: read_grid
  integer, intent(in) :: myrank

  integer :: i, nelements_read
  integer :: ierr
  character(len=256) :: msg
  real(kind=CUSTOM_REAL) :: dummy, val, cov
  integer :: i_, j_, k_, ind

  if (myrank == 0) then
  ! Reading the full model and grid by master CPU only.
    print *, 'Reading model from file ', trim(file_name)

    open(10, file=trim(file_name), status='old', action='read', iostat=ierr, iomsg=msg)
    if (ierr /= 0) call exit_MPI("Error in opening the model file! path=" &
                               //file_name//" iomsg="//msg, myrank, ierr)

    read(10, *, iostat=ierr) nelements_read
    if (ierr /= 0) call exit_MPI("Problem while reading the model file!", myrank, ierr)

    ! Sanity check.
    if (this%nelements_total /= nelements_read) &
      call exit_MPI("The grid is not correctly defined!"//new_line('a') &
                    //"nelements="//str(this%nelements)//new_line('a') &
                    //"nelements_read="//str(nelements_read)//new_line('a') &
                    //"nelements_total="//str(this%nelements_total), myrank, 0)

    if (read_grid) then
    ! Reading the full grid and the model.
      do i = 1, this%nelements_total
        read(10, *, iostat=ierr) this%grid_full%X1(i), this%grid_full%X2(i), &
                                 this%grid_full%Y1(i), this%grid_full%Y2(i), &
                                 this%grid_full%Z1(i), this%grid_full%Z2(i), &
                                 this%val_full(i), &
                                 this%grid_full%i_(i), this%grid_full%j_(i), this%grid_full%k_(i), &
                                 this%cov_full(i)

        ! Sanity check.
        if (this%grid_full%i_(i) < 1 .or. &
            this%grid_full%j_(i) < 1 .or. &
            this%grid_full%k_(i) < 1 .or. &
            this%grid_full%i_(i) > this%grid_full%nx .or. &
            this%grid_full%j_(i) > this%grid_full%ny .or. &
            this%grid_full%k_(i) > this%grid_full%nz) then

          call exit_MPI("The model grid dimensions in the Parfile are inconsistent with the model 3D indexes!"//new_line('a') &
                    //"i="//str(this%grid_full%i_(i))//new_line('a') &
                    //"j="//str(this%grid_full%j_(i))//new_line('a') &
                    //"k="//str(this%grid_full%k_(i)), myrank, 0)
        endif

        ! Store 1D grid index of the model parameter.
        this%grid_full%ind(this%grid_full%i_(i), this%grid_full%j_(i), this%grid_full%k_(i)) = i

        ! Sanity check.
        if (this%grid_full%X1(i) > this%grid_full%X2(i) .or. &
            this%grid_full%Y1(i) > this%grid_full%Y2(i) .or. &
            this%grid_full%Z1(i) > this%grid_full%Z2(i)) then
          call exit_MPI("The grid is not correctly defined (X1>X2 or Y1>Y2 or Z1>Z2)!", myrank, 0)
        endif

        if (ierr /= 0) call exit_MPI("Problem while reading the model file in model_read_voxels!", myrank, ierr)
      enddo

    else
    ! Reading the model only (without grid). Note: grid should already exist.
      do i = 1, this%nelements_total
        read(10, *, iostat=ierr) dummy, dummy, dummy, dummy, dummy, dummy, val, i_, j_, k_, cov

        ! Obtain 1D grid index of the model parameter.
        ind = this%grid_full%ind(i_, j_, k_)

        ! Set the model value.
        this%val_full(ind) = val

        ! Set the covariance value.
        this%cov_full(ind) = cov

        if (ierr /= 0) call exit_MPI("Problem while reading the model file in model_read_voxels!", myrank, ierr)
      enddo
    endif

    close(10)
  endif

end subroutine model_read_voxels_format

!================================================================================================
! Write the full model and grid in voxels format.
! Using the same format as in model_read_voxels_format subroutine.
!================================================================================================
subroutine model_write_voxels_format(this, file_name, myrank)
  class(t_model_IO), intent(in) :: this
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank

  integer :: i
  character(len=256) :: filename_full

  if (myrank == 0) then

    call system('mkdir -p '//trim(path_output)//"/Voxet/")

    filename_full = trim(path_output)//"/Voxet/"//file_name

    print *, 'Writing the full model to file ', trim(filename_full)

    open(27, file=trim(filename_full), access='stream', form='formatted', status='unknown', action='write')

    write (27, *) this%nelements_total

    ! Writing the full grid and the model.
    do i = 1, this%nelements_total
      write (27, *) this%grid_full%X1(i), this%grid_full%X2(i), &
                    this%grid_full%Y1(i), this%grid_full%Y2(i), &
                    this%grid_full%Z1(i), this%grid_full%Z2(i), &
                    this%val_full(i), &
                    this%grid_full%i_(i), this%grid_full%j_(i), this%grid_full%k_(i)
    enddo

    close(27)
  endif

end subroutine model_write_voxels_format

!======================================================================================================
! Write the model snapshots for visualization.
!======================================================================================================
subroutine model_write(this, name_prefix, gather_full_model, myrank, nbproc)
  class(t_model_IO), intent(inout) :: this
  character(len=*), intent(in) :: name_prefix
  logical, intent(in) :: gather_full_model
  integer, intent(in) :: myrank, nbproc

  type(t_parallel_tools) :: pt

  ! Note: model_write function uses values from val_full array.

  if (gather_full_model) then
    call pt%get_full_array(this%val, this%nelements, this%val_full, .false., myrank, nbproc)
  endif

  call this%write_paraview(name_prefix, myrank)
!  call this%write_qgis(name_prefix, myrank)

  ! Write the full model on voxels format.
  call this%write_voxels_format(name_prefix//"voxet_full.txt", myrank)

end subroutine model_write

!======================================================================================================
! Write the model snapshots in Paraview format for visualization.
!======================================================================================================
subroutine model_write_paraview(this, name_prefix, myrank)
  class(t_model_IO), intent(in) :: this
  character(len=*), intent(in) :: name_prefix
  integer, intent(in) :: myrank

  integer :: i, j, k, p, ierr

  ! For Paraview visualization.
  real(kind=CUSTOM_REAL), allocatable :: xgrid(:, :, :)
  real(kind=CUSTOM_REAL), allocatable :: ygrid(:, :, :)
  real(kind=CUSTOM_REAL), allocatable :: zgrid(:, :, :)
  real(kind=CUSTOM_REAL), allocatable :: model3d(:, :, :)

  real(kind=CUSTOM_REAL), allocatable :: xgrid_lego(:)
  real(kind=CUSTOM_REAL), allocatable :: ygrid_lego(:)
  real(kind=CUSTOM_REAL), allocatable :: zgrid_lego(:)

  real(kind=CUSTOM_REAL) :: Xmax, Ymax, Zmax
  !real(kind=CUSTOM_REAL) :: zshift
  integer :: step_x, step_y, step_z
  integer :: nx, ny, nz, mx, my, mz, ind
  !integer :: istep, jstep, kstep
  !integer :: x, y, z
  character(len=256) :: filename

  ! Write files my master CPU only.
  if (myrank /= 0) return

  print *, 'Writing models for Paraview visualization for ', trim(name_prefix)

  nx = this%grid_full%nx
  ny = this%grid_full%ny
  nz = this%grid_full%nz

  ! Start arrays with 0 as in the (general) interface of visualisation_paraview().
  allocate(xgrid(0:nx + 1, 0:ny + 1, 0:nz + 1), source=0._CUSTOM_REAL, stat=ierr)
  allocate(ygrid(0:nx + 1, 0:ny + 1, 0:nz + 1), source=0._CUSTOM_REAL, stat=ierr)
  allocate(zgrid(0:nx + 1, 0:ny + 1, 0:nz + 1), source=0._CUSTOM_REAL, stat=ierr)
  allocate(model3d(0:nx, 0:ny, 0:nz), source=0._CUSTOM_REAL, stat=ierr)

  allocate(xgrid_lego(1:(8 * this%nelements_total)), source=0._CUSTOM_REAL, stat=ierr)
  allocate(ygrid_lego(1:(8 * this%nelements_total)), source=0._CUSTOM_REAL, stat=ierr)
  allocate(zgrid_lego(1:(8 * this%nelements_total)), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_write!", myrank, ierr)

  ! Calculate 3D grid and model.
  do p = 1, this%nelements_total

    i = this%grid_full%i_(p)
    j = this%grid_full%j_(p)
    k = this%grid_full%k_(p)

    xgrid(i, j, k) = this%grid_full%X1(p)
    ygrid(i, j, k) = this%grid_full%Y1(p)
    zgrid(i, j, k) = this%grid_full%Z1(p)

    model3d(i, j, k) = this%val_full(p)

    !====================================
    ! Build lego-grid.
    !====================================
    ind = 8 * (p - 1) + 1
    ! z = 1
    xgrid_lego(ind + 0) = this%grid_full%X1(p)
    ygrid_lego(ind + 0) = this%grid_full%Y1(p)
    zgrid_lego(ind + 0) = this%grid_full%Z1(p)

    xgrid_lego(ind + 1) = this%grid_full%X2(p)
    ygrid_lego(ind + 1) = this%grid_full%Y1(p)
    zgrid_lego(ind + 1) = this%grid_full%Z1(p)

    xgrid_lego(ind + 2) = this%grid_full%X1(p)
    ygrid_lego(ind + 2) = this%grid_full%Y2(p)
    zgrid_lego(ind + 2) = this%grid_full%Z1(p)

    xgrid_lego(ind + 3) = this%grid_full%X2(p)
    ygrid_lego(ind + 3) = this%grid_full%Y2(p)
    zgrid_lego(ind + 3) = this%grid_full%Z1(p)

    ! z = 2
    xgrid_lego(ind + 4) = this%grid_full%X1(p)
    ygrid_lego(ind + 4) = this%grid_full%Y1(p)
    zgrid_lego(ind + 4) = this%grid_full%Z2(p)

    xgrid_lego(ind + 5) = this%grid_full%X2(p)
    ygrid_lego(ind + 5) = this%grid_full%Y1(p)
    zgrid_lego(ind + 5) = this%grid_full%Z2(p)

    xgrid_lego(ind + 6) = this%grid_full%X1(p)
    ygrid_lego(ind + 6) = this%grid_full%Y2(p)
    zgrid_lego(ind + 6) = this%grid_full%Z2(p)

    xgrid_lego(ind + 7) = this%grid_full%X2(p)
    ygrid_lego(ind + 7) = this%grid_full%Y2(p)
    zgrid_lego(ind + 7) = this%grid_full%Z2(p)

  enddo

  Xmax = this%grid_full%X2(this%grid_full%ind(nx, 1, 1))
  Ymax = this%grid_full%Y2(this%grid_full%ind(1, ny, 1))
  Zmax = this%grid_full%Z2(this%grid_full%ind(1, 1, nz))

  ! Adding last grid points plane:
  ! For example for nx = 3, there are nx + 1 = 4 grid planes:
  !   |          |          |          |
  ! X1(1)      X1(2)      X1(3)      X2(3) = Xmax
  !            X2(1)      X2(2)
  !
  do k = 1, nz + 1
    do j = 1, ny + 1
      do i = 1, nx + 1
        if (i == nx + 1 .or. j == ny + 1 .or. k == nz + 1) then
        ! A point lies on the outer boundary.

          xgrid(i, j, k) = Xmax
          ygrid(i, j, k) = Ymax
          zgrid(i, j, k) = Zmax

          if (i < nx + 1) xgrid(i, j, k) = xgrid(i, ny, nz)
          if (j < ny + 1) ygrid(i, j, k) = ygrid(nx, j, nz)
          if (k < nz + 1) zgrid(i, j, k) = zgrid(nx, ny, k)
        endif
      enddo
    enddo
  enddo

  ! Invert the Z-axis direction.
  zgrid = - zgrid
  zgrid_lego = - zgrid_lego

  step_x = 1
  step_y = 1
  step_z = 1

  ! To show the central slice when there are odd number of elements.
  mx = mod(nx, 2)
  my = mod(ny, 2)
  mz = mod(nz, 2)

  filename = trim(name_prefix)//"model3D_half_x.vtk"
  call visualisation_paraview(filename, myrank, nx, ny, nz, &
                              model3d, xgrid, ygrid, zgrid, &
                              nx / 2 + mx, nx / 2 + mx + step_x, 1, ny + 1, 1, nz + 1, &
                              step_x, step_y, step_z, 'CELL_DATA')

  filename = trim(name_prefix)//"model3D_half_y.vtk"
  call visualisation_paraview(filename, myrank, nx, ny, nz, &
                              model3d, xgrid, ygrid, zgrid, &
                              1, nx + 1, ny / 2 + my, ny / 2 + my + step_y, 1, nz + 1, &
                              step_x, step_y, step_z, 'CELL_DATA')

  filename = trim(name_prefix)//"model3D_full.vtk"
  call visualisation_paraview(filename, myrank, nx, ny, nz, &
                              model3d, xgrid, ygrid, zgrid, &
                              1, nx + 1, 1, ny + 1, 1, nz + 1, &
                              step_x, step_y, step_z, 'CELL_DATA')

  ! Write the full model using lego-grid.
  filename = trim(name_prefix)//"model3D_full_lego.vtk"
  call visualisation_paraview_legogrid(filename, myrank, this%nelements_total, this%val_full, xgrid_lego, ygrid_lego, zgrid_lego)

!  istep = 2
!  jstep = max(ny / 14, 1)
!  kstep = max(nz / 10, 1)
!
!  ! Write horizontal layers.
!  do k = 1, nz, kstep
!    ! Calculate depth to the middle of the voxel.
!    ind = this%grid_full%ind(1, 1, k)
!    z = int(this%grid_full%get_Z_cell_center(ind))
!
!    filename = trim(name_prefix)//"model3D_k="//trim(str(k))//"_z="//trim(str(z))//".vtk"
!
!    call visualisation_paraview(filename, myrank, nx, ny, nz, &
!                                model3d, xgrid, ygrid, zgrid, &
!                                1, nx + 1, 1, ny + 1, k, k + 1, &
!                                1, 1, 1, 'CELL_DATA')
!  enddo
!
!  ! Generate a set of (shifted) profiles that can be seen all on one screen when opened in Paraview.
!  zshift = 0.d0
!  do i = 2, nx, istep
!    ! Calculate depth to the middle of the voxel.
!    ind = this%grid_full%ind(i, 1, 1)
!    x = int(this%grid_full%get_X_cell_center(ind))
!
!    filename = trim(name_prefix)//"model3D_x_profiles_i="//trim(str(i))//"_x="//trim(str(x))//".vtk"
!
!    call visualisation_paraview(filename, myrank, nx, ny, nz, &
!                                model3d, xgrid, ygrid, zgrid, &
!                                i, i + 1, 1, ny + 1, 1, nz + 1, &
!                                step_x, step_y, step_z, 'CELL_DATA')
!
!    zgrid = zgrid + this%grid_full%get_hz() * (dble(nz) + 0.2d0)
!    zshift = zshift + this%grid_full%get_hz() * (dble(nz) + 0.2d0)
!  enddo
!
!  zgrid = zgrid - zshift
!
!  do j = 1, ny, jstep
!    ! Calculate depth to the middle of the voxel.
!    ind = this%grid_full%ind(1, j, 1)
!    y = int(this%grid_full%get_Y_cell_center(ind))
!
!    filename = trim(name_prefix)//"model3D_y_profiles_j="//trim(str(j))//"_y="//trim(str(y))//".vtk"
!
!    call visualisation_paraview(filename, myrank, nx, ny, nz, &
!                                model3d, xgrid, ygrid, zgrid, &
!                                1, nx + 1, j, j + 1, 1, nz + 1, &
!                                step_x, step_y, step_z, 'CELL_DATA')
!
!    zgrid = zgrid + this%grid_full%get_hz() * (dble(nz) + 0.2d0)
!  enddo

  ! Deallocate local arrays.
  deallocate(xgrid)
  deallocate(ygrid)
  deallocate(zgrid)
  deallocate(model3d)

end subroutine model_write_paraview

!======================================================================================================
! Write model snapshots in QGIS raster format for visualization.
!======================================================================================================
subroutine model_write_qgis(this, name_prefix, myrank)
  class(t_model_IO), intent(in) :: this
  character(len=*), intent(in) :: name_prefix
  integer, intent(in) :: myrank

  integer :: i, j, k, p, ierr

  real(kind=CUSTOM_REAL), allocatable :: model3d(:, :, :)

  real(kind=CUSTOM_REAL) :: cell_size
  integer :: istep, jstep, kstep
  integer :: nx, ny, nz
  integer :: nsections_x, nsections_y, nsections_z
  integer :: ind, x, y, z
  character(len=256) :: filename

  ! Write files my master CPU only.
  if (myrank /= 0) return

  nx = this%grid_full%nx
  ny = this%grid_full%ny
  nz = this%grid_full%nz

  allocate(model3d(1:nx, 1:ny, 1:nz), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_write_qgis!", myrank, ierr)

  ! Calculate 3D grid and model.
  do p = 1, this%nelements_total
    i = this%grid_full%i_(p)
    j = this%grid_full%j_(p)
    k = this%grid_full%k_(p)

    model3d(i, j, k) = this%val_full(p)
  enddo

  nsections_x = 10
  nsections_y = 10
  nsections_z = 10

  istep = max(nx / nsections_x, 1)
  jstep = max(ny / nsections_y, 1)
  kstep = max(nz / nsections_z, 1)

  ! QGIS does not support non-cubical cells. Use arbitrary cell size.
  cell_size = 100.d0

  ! Write YZ profiles.
  do i = 1, nx, istep
    ind = this%grid_full%ind(i, 1, 1)
    x = int(this%grid_full%get_X_cell_center(ind))

    filename = trim(name_prefix)//"model_YZ_x="//trim(str(x))//"_i="//trim(str(i))//".txt"
    call visualisation_qgis(filename, myrank, reshape(model3d(i:i, 1:ny, 1:nz), (/ny, nz/)), ny, nz, cell_size)
  enddo

  ! Write XZ profiles.
  do j = 1, ny, jstep
    ind = this%grid_full%ind(1, j, 1)
    y = int(this%grid_full%get_Y_cell_center(ind))

    filename = trim(name_prefix)//"model_XZ_y="//trim(str(y))//"_j="//trim(str(j))//".txt"
    call visualisation_qgis(filename, myrank, reshape(model3d(1:nx, j:j, 1:nz), (/nx, nz/)), nx, nz, cell_size)
  enddo

  ! Write XY profiles.
  do k = 1, nz, kstep
    ind = this%grid_full%ind(1, 1, k)
    z = int(this%grid_full%get_Z_cell_center(ind))

    filename = trim(name_prefix)//"model_XY_z="//trim(str(z))//"_k="//trim(str(k))//".txt"
    call visualisation_qgis(filename, myrank, reshape(model3d(1:nx, 1:ny, k:k), (/nx, ny/)), nx, ny, cell_size)
  enddo

  ! Deallocate local arrays.
  deallocate(model3d)

end subroutine model_write_qgis

end module model_IO
