
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

!================================================================================================
! This module contains only the I/O operations (reading/writing) with the model class.
!
! Vitaliy Ogarko, UWA, Australia.
!================================================================================================
module model_IO

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use paraview
  use string, only: str
  use parallel_tools
  use string
  use model
  use grid
  use parameters_inversion

  implicit none

  private

  public :: set_model
  private :: model_read
  public :: read_model_grid
  public :: model_write

  public :: set_model_bounds
  private :: read_bound_constraints
  public :: read_damping_gradient_weights

  private :: model_write_voxels_format
  private :: model_write_paraview

  logical, parameter, private :: WRITE_UNSTRUCTURED_GRID_PARAVIEW_MODEL = .false.

contains

!========================================================================================
! Sets the model values: via a constant from Parfile or via reading it from file.
!========================================================================================
subroutine set_model(model, model_type, model_val, model_file, myrank, nbproc)
  type(t_model), intent(inout) :: model
  integer, intent(in) :: model_type, myrank, nbproc
  real(kind=CUSTOM_REAL), intent(in) :: model_val
  character(len=256), intent(in) :: model_file

  if (myrank == 0) then
    if (model_type == 1) then
      ! Setting a constant value.
      model%val_full = model_val

    else if (model_type == 2) then
      ! Reading from file.
      call model_read(model, model_file, myrank)

    else
      call exit_MPI("Unknown model type in set_model!", myrank, model_type)
    endif

    ! Units conversion.
    model%val_full = model%val_full * model%units_mult

    if (model%ncomponents == 3  .and. model%grid_full%z_axis_dir /= 1) then
      ! Flip the Z-axis direction.
      model%val_full(:, 3) = -model%val_full(:, 3)
    endif
  endif

  ! Distribute the model values among CPUs.
  call model%distribute(myrank, nbproc)

end subroutine set_model

!================================================================================================
! Read the model from a file.
!================================================================================================
subroutine model_read(model, file_name, myrank)
  class(t_model), intent(inout) :: model
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank

  integer :: i, nelements_read
  integer :: ierr
  character(len=256) :: msg
  real(kind=CUSTOM_REAL) :: dummy(6), val(model%ncomponents)
  integer :: i_, j_, k_

  if (myrank == 0) then
  ! Reading the full model by the master CPU only.
    print *, 'Reading model from file ', trim(file_name)

    open(10, file=trim(file_name), status='old', action='read', iostat=ierr, iomsg=msg)

    if (ierr /= 0) call exit_MPI("Error in opening the model file! path=" &
                                 //file_name//" iomsg="//msg, myrank, ierr)

    read(10, *, iostat=ierr) nelements_read
    if (ierr /= 0) call exit_MPI("Problem while reading the model file!", myrank, ierr)

    ! Sanity check.
    if (model%nelements_total /= nelements_read) &
      call exit_MPI("The grid is not correctly defined!"//new_line('a') &
                    //"nelements="//str(model%nelements)//new_line('a') &
                    //"nelements_read="//str(nelements_read)//new_line('a') &
                    //"nelements_total="//str(model%nelements_total), myrank, 0)

    ! Reading the model only (without grid).
    do i = 1, model%nelements_total
      ! Note we read an array of val.
      read(10, *, iostat=ierr) dummy, val, i_, j_, k_

      ! Set the model value.
      model%val_full(i, :) = val

      if (ierr /= 0) call exit_MPI("Problem while reading the model file in model_read_voxels!", myrank, ierr)
    enddo

    close(10)
  endif

end subroutine model_read

!================================================================================================
! Read the model grid from a file.
!================================================================================================
subroutine read_model_grid(grid, nmodel_components, file_name, myrank)
  class(t_grid), intent(inout) :: grid
  integer, intent(in) :: nmodel_components
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank

  integer :: p, nelements_read
  integer :: i, j, k
  integer :: nelements_total
  integer :: ierr
  character(len=256) :: msg
  real(kind=CUSTOM_REAL) :: tmp
  real(kind=CUSTOM_REAL) :: val(nmodel_components)
  logical :: correct_order

  if (myrank == 0) then
  ! Reading the full grid by master CPU only.
    print *, 'Reading model grid from file ', trim(file_name)

    open(10, file=trim(file_name), status='old', action='read', iostat=ierr, iomsg=msg)

    if (ierr /= 0) call exit_MPI("Error in opening the model grid file! path=" &
                                 //file_name//" iomsg="//msg, myrank, ierr)

    read(10, *, iostat=ierr) nelements_read
    if (ierr /= 0) call exit_MPI("Problem while reading the model grid file!", myrank, ierr)

    nelements_total = grid%nx * grid%ny * grid%nz

    ! Sanity check.
    if (nelements_total /= nelements_read) &
      call exit_MPI("The grid is not correctly defined!"//new_line('a') &
                    //"nelements_read="//str(nelements_read)//new_line('a') &
                    //"nelements_total="//str(nelements_total), myrank, 0)

    correct_order = .true.

    ! Reading the full grid.
    do p = 1, nelements_total
      read(10, *, iostat=ierr) grid%X1(p), grid%X2(p), &
                               grid%Y1(p), grid%Y2(p), &
                               grid%Z1(p), grid%Z2(p), &
                               val, &
                               i, j, k

      if (ierr /= 0) call exit_MPI("Problem while reading the grid file in grid_read!", myrank, ierr)

      ! Sanity check.
      if (i < 1 .or. j < 1 .or. k < 1 .or. &
          i > grid%nx .or. j > grid%ny .or. k > grid%nz) then
        call exit_MPI("The model grid dimensions in the Parfile are inconsistent with the model 3D indexes!"//new_line('a') &
                  //"i ="//str(i)//new_line('a') &
                  //"j ="//str(j)//new_line('a') &
                  //"k ="//str(k), myrank, 0)
      endif

      ! Sanity check.
      if (grid%X1(p) >= grid%X2(p) .or. &
          grid%Y1(p) >= grid%Y2(p) .or. &
          grid%Z1(p) >= grid%Z2(p)) then
        call exit_MPI("The grid is not correctly defined (X1 >= X2 or Y1 >= Y2 or Z1 >= Z2)!", myrank, 0)
      endif

      ! Test that the grid cell order is i-j-k.
      if (p == 1 .and. (i /= 1 .or. j /= 1 .or. k /= 1)) correct_order = .false.
      if (p == 2 .and. (i /= 2 .or. j /= 1 .or. k /= 1)) correct_order = .false.
      if (p == grid%nx + 1 .and. (i /= 1 .or. j /= 2 .or. k /= 1)) correct_order = .false.
      if (p == grid%nx * grid%ny + 1 .and. (i /= 1 .or. j /= 1 .or. k /= 2)) correct_order = .false.

      if (.not. correct_order) then
        call exit_MPI("Wrong cell order in the model grid file! Use the i-j-k order (i is the fastest index)!", myrank, 0)
      endif

    enddo
    close(10)

    ! Flip the Z-axis direction.
    if (grid%z_axis_dir /= 1) then
      do p = 1, nelements_total
        tmp = grid%Z1(p)
        grid%Z1(p) = -grid%Z2(p)
        grid%Z2(p) = -tmp
      enddo
    endif

    print *, 'Xmin, Xmax, SizeX =', grid%get_Xmin(), grid%get_Xmax(), grid%get_Xmax()- grid%get_Xmin()
    print *, 'Ymin, Ymax, SizeY =', grid%get_Ymin(), grid%get_Ymax(), grid%get_Ymax()- grid%get_Ymin()
    print *, 'Zmin, Zmax, SizeZ =', grid%get_Zmin(), grid%get_Zmax(), grid%get_Zmax()- grid%get_Zmin()
  endif ! myrank == 0

  ! Broadcast full grid to all CPUs.
  call grid%broadcast(myrank)

end subroutine read_model_grid

!========================================================================================
! Sets the model bounds (for the ADMM).
!========================================================================================
subroutine set_model_bounds(ipar, model, problem_type, myrank, nbproc)
  type(t_parameters_inversion), intent(in) :: ipar
  integer, intent(in) :: problem_type
  integer, intent(in) :: myrank, nbproc
  type(t_model), intent(inout) :: model

  integer :: i

  ! Allocate bound arrays.
  call model%allocate_bound_arrays(ipar%nlithos, myrank)

  if (ipar%admm_bound_type == 1) then
    ! Global bounds - define from Parfile parameters.
    do i = 1, model%nelements
      model%min_bound(:, i) = ipar%admm_bounds(problem_type)%val(1::2)
      model%max_bound(:, i) = ipar%admm_bounds(problem_type)%val(2::2)
    enddo
    model%bound_weight(:) = 1.d0

    ! Sanity check.
    do i = 1, ipar%nlithos
      if (model%min_bound(i, 1) > model%max_bound(i, 1)) then
        call exit_MPI("Wrong admm bounds: define bounds as: min1 max1 ... minN maxN.", myrank, 0)
      endif
    enddo
  else
    ! Local bounds - read from file.
    call read_bound_constraints(model, ipar%bounds_ADMM_file(problem_type), myrank, nbproc)
  endif

  ! Units conversion.
  model%min_bound = model%min_bound * model%units_mult
  model%max_bound = model%max_bound * model%units_mult

end subroutine set_model_bounds

!==========================================================================================================
! Read the local bound constraints (for the ADMM).
!==========================================================================================================
subroutine read_bound_constraints(model, file_name, myrank, nbproc)
  class(t_model), intent(inout) :: model
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank, nbproc

  integer :: ierr, nsmaller, ind, i, j
  integer :: nelements_read, nlithos_read
  character(len=256) :: msg
  real(kind=CUSTOM_REAL) :: model_bounds(2 * model%nlithos)

  if (myrank == 0) print *, 'Reading local bound constraints from file ', trim(file_name)

  open(10, file=trim(file_name), status='old', action='read', iostat=ierr, iomsg=msg)
  if (ierr /= 0) call exit_MPI("Error in opening the bound constraints file! path=" &
                 //file_name//" iomsg="//msg, myrank, ierr)

  read(10, *, iostat=ierr) nelements_read, nlithos_read
  if (ierr /= 0) call exit_MPI("Problem while reading the bound constraints file!", myrank, ierr)

  ! Sanity check.
  if (model%nelements_total /= nelements_read) &
    call exit_MPI("The constraints are not correctly defined!"//new_line('a') &
          //"nelements="//str(model%nelements)//new_line('a') &
          //"nelements_read="//str(nelements_read)//new_line('a') &
          //"nelements_total="//str(model%nelements_total), myrank, 0)

  if (model%nlithos /= nlithos_read) &
    call exit_MPI("The constraints are not correctly defined!"//new_line('a') &
          //"nlithos="//str(model%nlithos)//new_line('a') &
          //"nlithos_read="//str(nlithos_read), myrank, 0)

  if (myrank == 0) print *, 'Read nelements, nlithos = ', nelements_read, nlithos_read

  ! The number of elements on CPUs with rank smaller than myrank.
  nsmaller = get_nsmaller(model%nelements, myrank, nbproc)

  ! Reading.
  do i = 1, model%nelements_total
    if (i > nsmaller .and. i <= nsmaller + model%nelements) then
      ind = i - nsmaller
      read(10, *, iostat=ierr) model_bounds, model%bound_weight(ind)

      if (ierr /= 0) then
        call exit_MPI("Problem with reading the bound constraints for pixel i = "//str(i), myrank, ierr)
      endif

      model%min_bound(:, ind) = model_bounds(1::2)
      model%max_bound(:, ind) = model_bounds(2::2)

      ! Sanity check.
      do j = 1, model%nlithos
        if (model%min_bound(j, ind) > model%max_bound(j, ind)) then
          call exit_MPI("Wrong admm bounds: define bounds as: min1 max1 ... minN maxN.", myrank, ind)
        endif
      enddo

    else
      ! Skip line.
      read(10, *)
    endif

    if (i > nsmaller + model%nelements) then
      exit
    endif
  enddo

  close(10)

end subroutine read_bound_constraints

!==========================================================================================================
! Read damping gradient local weights.
!==========================================================================================================
subroutine read_damping_gradient_weights(model, file_name, myrank)
  class(t_model), intent(inout) :: model
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank

  integer :: ierr, i, nelements_read
  character(len=256) :: msg

  if (myrank == 0) print *, 'Reading damping gradient weights from file ', trim(file_name)

  open(10, file=trim(file_name), status='old', action='read', iostat=ierr, iomsg=msg)
  if (ierr /= 0) call exit_MPI("Error in opening the damping gradient weights file! path=" &
                 //file_name//" iomsg="//msg, myrank, ierr)

  read(10, *, iostat=ierr) nelements_read
  if (ierr /= 0) call exit_MPI("Problem while reading the damping gradient weights file!", myrank, ierr)

  ! Sanity check.
  if (model%nelements_total /= nelements_read) &
    call exit_MPI("The damping gradient weights are not correctly defined!"//new_line('a') &
          //"nelements="//str(model%nelements)//new_line('a') &
          //"nelements_read="//str(nelements_read)//new_line('a') &
          //"nelements_total="//str(model%nelements_total), myrank, 0)

  ! Reading.
  do i = 1, model%nelements_total
      read(10, *, iostat=ierr) &
        model%damping_grad_weight(i, 1), model%damping_grad_weight(i, 2), model%damping_grad_weight(i, 3)

      if (ierr /= 0) &
        call exit_MPI("Problem with reading the damping gradient weights for pixel i = "//str(i), myrank, ierr)
  enddo

  close(10)

end subroutine read_damping_gradient_weights

!======================================================================================================
! Write the model snapshots for visualization.
!======================================================================================================
subroutine model_write(model, name_prefix, gather_full_model, write_voxet, myrank, nbproc)
  class(t_model), intent(inout) :: model
  character(len=*), intent(in) :: name_prefix
  logical, intent(in) :: gather_full_model, write_voxet
  integer, intent(in) :: myrank, nbproc

  if (gather_full_model) then
    call model%update_full(.false., myrank, nbproc)
  endif

  ! Write the model in vtk format.
  call model_write_paraview(model, name_prefix, myrank)

  if (write_voxet) then
    ! Write the full model in voxels format.
    call model_write_voxels_format(model, name_prefix//"voxet_full.txt", myrank)
  endif

end subroutine model_write

!================================================================================================
! Write the full model (in voxels format) to file.
!================================================================================================
subroutine model_write_voxels_format(model, file_name, myrank)
  class(t_model), intent(in) :: model
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank

  character(len=256) :: filename_full
  integer :: i, ierr
  ! Temporary array for writing the model to file.
  real(kind=CUSTOM_REAL), allocatable :: val_full(:, :)

  if (myrank == 0) then
    ! Create a directory.
    call execute_command_line('mkdir -p "'//trim(path_output)//'/Voxet"')

    filename_full = trim(path_output)//"/Voxet/"//file_name

    print *, 'Writing the full model to file ', trim(filename_full)

    allocate(val_full(model%nelements_total, model%ncomponents), source=0._CUSTOM_REAL, stat=ierr)

    ! Units conversion.
    val_full = model%val_full / model%units_mult

    open(27, file=trim(filename_full), access='stream', form='formatted', status='replace', action='write')

    write(27, *) model%nelements_total

    ! Write the full model array.
    ! Use the loop instead of writing the whole array to write model components in different columns.
    write(27, *) (val_full(i, :), new_line("A"), i = 1, model%nelements_total)

    close(27)
    deallocate(val_full)
  endif

end subroutine model_write_voxels_format

!======================================================================================================
! Write the model snapshots in Paraview format for visualization.
!======================================================================================================
subroutine model_write_paraview(model, name_prefix, myrank)
  class(t_model), intent(in) :: model
  character(len=*), intent(in) :: name_prefix
  integer, intent(in) :: myrank

  integer :: nx, ny, nz
  character(len=256) :: filename
  logical :: INVERT_Z_AXIS

  ! Write files my master CPU only.
  if (myrank /= 0) return

  print *, 'Writing models for Paraview visualization for ', trim(name_prefix)

  nx = model%grid_full%nx
  ny = model%grid_full%ny
  nz = model%grid_full%nz

  INVERT_Z_AXIS = .true.

  if (WRITE_UNSTRUCTURED_GRID_PARAVIEW_MODEL) then
    ! Write the full model (using unstructured grid format).
    ! In this format each cell has constant value, while in the structured grid case the values get smoothed between the cell centers.
    ! But this format needs to allocate ~6 times more memory, and the output file has much larger size.
    filename = trim(name_prefix)//"model3D_full_lego.vtk"
    call visualisation_paraview_legogrid(filename, myrank, model%nelements_total, model%ncomponents, model%val_full, &
                                         model%grid_full%X1, model%grid_full%Y1, model%grid_full%Z1, &
                                         model%grid_full%X2, model%grid_full%Y2, model%grid_full%Z2, &
                                         model%grid_full%nx, model%grid_full%ny, model%grid_full%nz, &
                                         1, nx, 1, ny, 1, nz, &
                                         INVERT_Z_AXIS, model%units_mult, model%vtk_label)
  endif

  ! Write the full model using structured grid.
  filename = trim(name_prefix)//"model3D_full.vtk"
  call visualisation_paraview_struct_grid(filename, myrank, model%nelements_total, model%ncomponents, model%val_full, &
                                       model%grid_full%X1, model%grid_full%Y1, model%grid_full%Z1, &
                                       model%grid_full%X2, model%grid_full%Y2, model%grid_full%Z2, &
                                       model%grid_full%nx, model%grid_full%ny, model%grid_full%nz, &
                                       1, nx, 1, ny, 1, nz, &
                                       INVERT_Z_AXIS, model%units_mult, model%vtk_label)

  ! Write the x-profile of the model.
  filename = trim(name_prefix)//"model3D_half_x.vtk"
  call visualisation_paraview_legogrid(filename, myrank, model%nelements_total, model%ncomponents, model%val_full, &
                                       model%grid_full%X1, model%grid_full%Y1, model%grid_full%Z1, &
                                       model%grid_full%X2, model%grid_full%Y2, model%grid_full%Z2, &
                                       model%grid_full%nx, model%grid_full%ny, model%grid_full%nz, &
                                       nx / 2 + 1, nx / 2 + 1, 1, ny, 1, nz, &
                                       INVERT_Z_AXIS, model%units_mult, model%vtk_label)

  ! Write the y-profile of the model.
  filename = trim(name_prefix)//"model3D_half_y.vtk"
  call visualisation_paraview_legogrid(filename, myrank, model%nelements_total, model%ncomponents, model%val_full, &
                                       model%grid_full%X1, model%grid_full%Y1, model%grid_full%Z1, &
                                       model%grid_full%X2, model%grid_full%Y2, model%grid_full%Z2, &
                                       model%grid_full%nx, model%grid_full%ny, model%grid_full%nz, &
                                       1, nx, ny / 2 + 1, ny / 2 + 1, 1, nz, &
                                       INVERT_Z_AXIS, model%units_mult, model%vtk_label)

  ! Write the z-profile of the model.
  filename = trim(name_prefix)//"model3D_half_z.vtk"
  call visualisation_paraview_legogrid(filename, myrank, model%nelements_total, model%ncomponents, model%val_full, &
                                       model%grid_full%X1, model%grid_full%Y1, model%grid_full%Z1, &
                                       model%grid_full%X2, model%grid_full%Y2, model%grid_full%Z2, &
                                       model%grid_full%nx, model%grid_full%ny, model%grid_full%nz, &
                                       1, nx, 1, ny, nz / 2 + 1, nz / 2 + 1, &
                                       INVERT_Z_AXIS, model%units_mult, model%vtk_label)
end subroutine model_write_paraview

end module model_IO
