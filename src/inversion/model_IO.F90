
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
  use model_base

  implicit none

  private

  public :: model_read
  public :: model_write
  public :: model_read_bound_constraints

  private :: model_read_voxels_format
  private :: model_write_voxels_format
  private :: model_write_paraview

contains

!==========================================================================================================
! Read the full model and grid, and then broadcast to all CPUs.
!==========================================================================================================
subroutine model_read(model, file_name, read_grid, myrank, nbproc)
  class(t_model_base), intent(inout) :: model
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank, nbproc
  logical, intent(in) :: read_grid

  integer :: ierr

  call model_read_voxels_format(model, file_name, read_grid, myrank, nbproc)

  if (read_grid) then
    ! Broadcast full grid to all CPUs.
    call model%grid_full%broadcast(myrank)
  endif

  ! Broadcast full model to all CPUs.
  call MPI_Bcast(model%val_full, model%nelements_total, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("Error in MPI_Bcast in model_read!", myrank, ierr)

  if (myrank == 0) then
    print *, 'Xmin, Xmax =', model%get_Xmin(), model%get_Xmax()
    print *, 'Ymin, Ymax =', model%get_Ymin(), model%get_Ymax()
  endif

end subroutine model_read

!==========================================================================================================
! Read the local bound constraints (for ADMM).
!==========================================================================================================
subroutine model_read_bound_constraints(model, file_name, myrank, nbproc)
  class(t_model_base), intent(inout) :: model
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
  nsmaller = pt%get_nsmaller(model%nelements, myrank, nbproc)

  ! Reading.
  do i = 1, model%nelements_total
    if (i > nsmaller .and. i <= nsmaller + model%nelements) then
      ind = i - nsmaller
      read(10, *, iostat=ierr) &
        model%min_local_bound(ind, :), model%max_local_bound(ind, :), model%local_bound_constraints_weight(ind)
    else
      read(10, '(A)', iostat=ierr) dummy_line
    endif

    if (i > nsmaller + model%nelements) then
      exit
    endif
  enddo

  close(10)

end subroutine model_read_bound_constraints

!================================================================================================
! Read the full model and grid in voxels format.
!================================================================================================
subroutine model_read_voxels_format(model, file_name, read_grid, myrank, nbproc)
  class(t_model_base), intent(inout) :: model
  character(len=*), intent(in) :: file_name
  logical, intent(in) :: read_grid
  integer, intent(in) :: myrank, nbproc

  integer :: i, nelements_read
  integer :: ierr
  character(len=256) :: msg
  real(kind=CUSTOM_REAL) :: dummy, val, cov
  integer :: i_, j_, k_

  ! Displacement for mpi_scatterv.
  integer :: displs(nbproc)
  ! The number of elements on every CPU for mpi_scatterv.
  integer :: nelements_at_cpu(nbproc)
  type(t_parallel_tools) :: pt
  real(kind=CUSTOM_REAL), allocatable :: cov_full(:)

  if (myrank == 0) then
  ! Reading the full model and grid by master CPU only.
    allocate(cov_full(model%nelements_total), stat=ierr)

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

    if (read_grid) then
    ! Reading the full grid and the model.
      do i = 1, model%nelements_total
        read(10, *, iostat=ierr) model%grid_full%X1(i), model%grid_full%X2(i), &
                                 model%grid_full%Y1(i), model%grid_full%Y2(i), &
                                 model%grid_full%Z1(i), model%grid_full%Z2(i), &
                                 model%val_full(i), &
                                 model%grid_full%i_(i), model%grid_full%j_(i), model%grid_full%k_(i), &
                                 cov_full(i)

        ! Sanity check.
        if (model%grid_full%i_(i) < 1 .or. &
            model%grid_full%j_(i) < 1 .or. &
            model%grid_full%k_(i) < 1 .or. &
            model%grid_full%i_(i) > model%grid_full%nx .or. &
            model%grid_full%j_(i) > model%grid_full%ny .or. &
            model%grid_full%k_(i) > model%grid_full%nz) then

          call exit_MPI("The model grid dimensions in the Parfile are inconsistent with the model 3D indexes!"//new_line('a') &
                    //"i="//str(model%grid_full%i_(i))//new_line('a') &
                    //"j="//str(model%grid_full%j_(i))//new_line('a') &
                    //"k="//str(model%grid_full%k_(i)), myrank, 0)
        endif

        ! Store 1D grid index of the model parameter.
        model%grid_full%ind(model%grid_full%i_(i), model%grid_full%j_(i), model%grid_full%k_(i)) = i

        ! Sanity check.
        if (model%grid_full%X1(i) > model%grid_full%X2(i) .or. &
            model%grid_full%Y1(i) > model%grid_full%Y2(i) .or. &
            model%grid_full%Z1(i) > model%grid_full%Z2(i)) then
          call exit_MPI("The grid is not correctly defined (X1>X2 or Y1>Y2 or Z1>Z2)!", myrank, 0)
        endif

        if (ierr /= 0) call exit_MPI("Problem while reading the model file in model_read_voxels!", myrank, ierr)
      enddo

    else
    ! Reading the model only (without grid).
      do i = 1, model%nelements_total
        read(10, *, iostat=ierr) dummy, dummy, dummy, dummy, dummy, dummy, val, i_, j_, k_, cov

        ! Set the model value.
        model%val_full(i) = val

        ! Set the covariance value.
        cov_full(i) = cov

        if (ierr /= 0) call exit_MPI("Problem while reading the model file in model_read_voxels!", myrank, ierr)
      enddo
    endif

    close(10)
  endif

  !------------------------------------------------------------------------------
  ! Distribute the covarianve values among CPUS.
  !------------------------------------------------------------------------------
  ! Partitioning for MPI_Scatterv.
  call pt%get_mpi_partitioning(model%nelements, displs, nelements_at_cpu, myrank, nbproc)

  call MPI_Scatterv(cov_full, nelements_at_cpu, displs, CUSTOM_MPI_TYPE, &
                    model%cov, model%nelements, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (allocated(cov_full)) deallocate(cov_full)

end subroutine model_read_voxels_format

!======================================================================================================
! Write the model snapshots for visualization.
!======================================================================================================
subroutine model_write(model, name_prefix, gather_full_model, myrank, nbproc)
  class(t_model_base), intent(inout) :: model
  character(len=*), intent(in) :: name_prefix
  logical, intent(in) :: gather_full_model
  integer, intent(in) :: myrank, nbproc

  type(t_parallel_tools) :: pt

  ! Note: model_write function uses values from val_full array.

  if (gather_full_model) then
    call pt%get_full_array(model%val, model%nelements, model%val_full, .false., myrank, nbproc)
  endif

  ! Write the model in vtk format.
  call model_write_paraview(model, name_prefix, myrank)

  ! Write the full model in voxels format.
  call model_write_voxels_format(model, name_prefix//"voxet_full.txt", myrank)

end subroutine model_write

!================================================================================================
! Write the full model and grid in voxels format.
! Using the same format as in model_read_voxels_format subroutine.
!================================================================================================
subroutine model_write_voxels_format(model, file_name, myrank)
  class(t_model_base), intent(in) :: model
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank

  integer :: i
  character(len=256) :: filename_full

  if (myrank == 0) then

    call system('mkdir -p '//trim(path_output)//"/Voxet/")

    filename_full = trim(path_output)//"/Voxet/"//file_name

    print *, 'Writing the full model to file ', trim(filename_full)

    open(27, file=trim(filename_full), access='stream', form='formatted', status='unknown', action='write')

    write (27, *) model%nelements_total

    ! Writing the full grid and the model.
    do i = 1, model%nelements_total
      write (27, *) model%grid_full%X1(i), model%grid_full%X2(i), &
                    model%grid_full%Y1(i), model%grid_full%Y2(i), &
                    model%grid_full%Z1(i), model%grid_full%Z2(i), &
                    model%val_full(i), &
                    model%grid_full%i_(i), model%grid_full%j_(i), model%grid_full%k_(i)
    enddo

    close(27)
  endif

end subroutine model_write_voxels_format

!======================================================================================================
! Write the model snapshots in Paraview format for visualization.
!======================================================================================================
subroutine model_write_paraview(model, name_prefix, myrank)
  class(t_model_base), intent(in) :: model
  character(len=*), intent(in) :: name_prefix
  integer, intent(in) :: myrank

  integer :: nx, ny, nz
  character(len=256) :: filename

  ! Write files my master CPU only.
  if (myrank /= 0) return

  print *, 'Writing models for Paraview visualization for ', trim(name_prefix)

  nx = model%grid_full%nx
  ny = model%grid_full%ny
  nz = model%grid_full%nz

  ! Write the full model.
  filename = trim(name_prefix)//"model3D_full.vtk"
  call visualisation_paraview_legogrid(filename, myrank, model%nelements_total, model%val_full, &
                                       model%grid_full%X1, model%grid_full%Y1, model%grid_full%Z1, &
                                       model%grid_full%X2, model%grid_full%Y2, model%grid_full%Z2, &
                                       model%grid_full%i_, model%grid_full%j_, model%grid_full%k_, &
                                       1, nx, 1, ny, 1, nz, &
                                       .true.)

  ! Write the x-profile of the model.
  filename = trim(name_prefix)//"model3D_half_x.vtk"
  call visualisation_paraview_legogrid(filename, myrank, model%nelements_total, model%val_full, &
                                       model%grid_full%X1, model%grid_full%Y1, model%grid_full%Z1, &
                                       model%grid_full%X2, model%grid_full%Y2, model%grid_full%Z2, &
                                       model%grid_full%i_, model%grid_full%j_, model%grid_full%k_, &
                                       nx / 2, nx / 2, 1, ny, 1, nz, &
                                       .true.)

  ! Write the y-profile of the model.
  filename = trim(name_prefix)//"model3D_half_y.vtk"
  call visualisation_paraview_legogrid(filename, myrank, model%nelements_total, model%val_full, &
                                       model%grid_full%X1, model%grid_full%Y1, model%grid_full%Z1, &
                                       model%grid_full%X2, model%grid_full%Y2, model%grid_full%Z2, &
                                       model%grid_full%i_, model%grid_full%j_, model%grid_full%k_, &
                                       1, nx, ny / 2, ny / 2, 1, nz, &
                                       .true.)

  ! Write the z-profile of the model.
  filename = trim(name_prefix)//"model3D_half_z.vtk"
  call visualisation_paraview_legogrid(filename, myrank, model%nelements_total, model%val_full, &
                                       model%grid_full%X1, model%grid_full%Y1, model%grid_full%Z1, &
                                       model%grid_full%X2, model%grid_full%Y2, model%grid_full%Z2, &
                                       model%grid_full%i_, model%grid_full%j_, model%grid_full%k_, &
                                       1, nx, 1, ny, nz / 2, nz / 2, &
                                       .true.)
end subroutine model_write_paraview

end module model_IO
