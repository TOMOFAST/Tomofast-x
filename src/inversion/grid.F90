
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
! A class to work with model grids for forward and inverse problems.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!================================================================================================
module grid

  use, intrinsic :: iso_c_binding

  use global_typedefs
  use mpi_tools, only: exit_MPI

  implicit none

  private

  type, public :: t_grid
    ! Grid dimensions.
    integer :: nx, ny, nz

    ! Direction of the Z-axis (1 = down, -1 = up).
    integer :: z_axis_dir

    ! Shared memory variables.
    integer :: comm_shm
    integer :: shm_rank
    integer :: win_X1, win_Y1, win_Z1
    integer :: win_X2, win_Y2, win_Z2
    type(c_ptr) :: baseptr_X1, baseptr_Y1, baseptr_Z1
    type(c_ptr) :: baseptr_X2, baseptr_Y2, baseptr_Z2

    integer :: is_first_group

    ! The starting and ending coordinates of the 3D prism for each model grid cell.
    real(kind=CUSTOM_REAL), pointer, dimension(:) :: X1, Y1, Z1
    real(kind=CUSTOM_REAL), pointer, dimension(:) :: X2, Y2, Z2

  contains
    private

    procedure, public, pass :: allocate => grid_allocate
    procedure, public, pass :: deallocate => grid_deallocate

    procedure, public, pass :: get_hx => grid_get_hx
    procedure, public, pass :: get_hy => grid_get_hy
    procedure, public, pass :: get_hz => grid_get_hz

    procedure, public, pass :: get_X_cell_center => grid_get_X_cell_center
    procedure, public, pass :: get_Y_cell_center => grid_get_Y_cell_center
    procedure, public, pass :: get_Z_cell_center => grid_get_Z_cell_center

    procedure, public, pass :: get_cell_volume => grid_get_cell_volume

    procedure, public, pass :: get_Xmin => grid_get_Xmin
    procedure, public, pass :: get_Xmax => grid_get_Xmax
    procedure, public, pass :: get_Ymin => grid_get_Ymin
    procedure, public, pass :: get_Ymax => grid_get_Ymax
    procedure, public, pass :: get_Zmin => grid_get_Zmin
    procedure, public, pass :: get_Zmax => grid_get_Zmax

  end type t_grid

  !-------------------------------------------------------------------------------------
  ! A memory efficient grid for calculating the gradients.
  ! It only stores cell sizes along each dimension, i.e., the total memory is O(nx + ny + nz).
  type, public :: t_grad_grid

    ! Cell sizes along each dimension.
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: dX, dY, dZ

    ! Full grid dimensions.
    integer :: nx, ny, nz
  contains
    private
    procedure, public, pass :: init => grad_grid_init
    procedure, public, pass :: get_ind => grad_grid_get_ind

  end type t_grad_grid

contains

!=======================================================================================
! Allocate grid arrays.
!=======================================================================================
subroutine grid_allocate(this, nx, ny, nz, z_axis_dir, myrank)
  class(t_grid), intent(inout) :: this
  integer, intent(in) :: nx, ny, nz, z_axis_dir, myrank

  integer :: nelements_total
  integer :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: size
  integer :: disp_unit
  integer :: shm_group_size
!  integer :: color, key
  integer :: shape_arr(1)

  this%nx = nx
  this%ny = ny
  this%nz = nz

  this%z_axis_dir = z_axis_dir

  nelements_total = nx * ny * nz

  !------------------------------------------------------------------------------
  ! Shared memory allocation.
  !------------------------------------------------------------------------------
!  ! Split MPI ranks into even and odd groups
!  if (mod(myrank, 2) == 0) then
!     color = 0  ! Even group
!  else
!     color = 1  ! Odd group
!  end if
!  key = myrank  ! Maintain original rank order in new communicators
!
!  call MPI_Comm_split(MPI_COMM_WORLD, color, key, this%comm_shm, ierr)

  ! Create a shared memory communicator.
  call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, this%comm_shm, ierr)

  ! Define the first group's flag.
  this%is_first_group = 0
  if (myrank == 0) then
     this%is_first_group = 1
  endif
  ! Broadcast the flag of the first group to all ranks
  call MPI_Bcast(this%is_first_group, 1, MPI_INT, 0, this%comm_shm, ierr)

  call MPI_Comm_rank(this%comm_shm, this%shm_rank, ierr)
  call MPI_Comm_size(this%comm_shm, shm_group_size, ierr)

  if (this%shm_rank == 0) then
    print *, 'myrank, shm_rank, shm_group_size =', myrank, this%shm_rank, shm_group_size
  endif

  ! Define the size of the shared memory array (allocate only on rank 0 of shared communicator).
  if (this%shm_rank == 0) then
     size = nelements_total * CUSTOM_REAL
  else
     size = 0  ! Other ranks request zero size (they will query instead).
  end if
  disp_unit = CUSTOM_REAL

  ! Allocate shared memory.
  call MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, this%comm_shm, this%baseptr_X1, this%win_X1, ierr)
  call MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, this%comm_shm, this%baseptr_Y1, this%win_Y1, ierr)
  call MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, this%comm_shm, this%baseptr_Z1, this%win_Z1, ierr)
  call MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, this%comm_shm, this%baseptr_X2, this%win_X2, ierr)
  call MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, this%comm_shm, this%baseptr_Y2, this%win_Y2, ierr)
  call MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, this%comm_shm, this%baseptr_Z2, this%win_Z2, ierr)

  if (ierr /= 0) call exit_MPI("Shared memory allocation error in grid_allocate!", myrank, ierr)

  ! Query the shared memory pointer on all ranks
  call MPI_Win_shared_query(this%win_X1, 0, size, disp_unit, this%baseptr_X1, ierr)
  call MPI_Win_shared_query(this%win_Y1, 0, size, disp_unit, this%baseptr_Y1, ierr)
  call MPI_Win_shared_query(this%win_Z1, 0, size, disp_unit, this%baseptr_Z1, ierr)
  call MPI_Win_shared_query(this%win_X2, 0, size, disp_unit, this%baseptr_X2, ierr)
  call MPI_Win_shared_query(this%win_Y2, 0, size, disp_unit, this%baseptr_Y2, ierr)
  call MPI_Win_shared_query(this%win_Z2, 0, size, disp_unit, this%baseptr_Z2, ierr)

  if (ierr /= 0) call exit_MPI("Shared memory query error in grid_allocate!", myrank, ierr)

  shape_arr(1) = nelements_total

  ! Convert C pointer to Fortran pointer.
  call c_f_pointer(this%baseptr_X1, this%X1, shape_arr)
  call c_f_pointer(this%baseptr_Y1, this%Y1, shape_arr)
  call c_f_pointer(this%baseptr_Z1, this%Z1, shape_arr)
  call c_f_pointer(this%baseptr_X2, this%X2, shape_arr)
  call c_f_pointer(this%baseptr_Y2, this%Y2, shape_arr)
  call c_f_pointer(this%baseptr_Z2, this%Z2, shape_arr)

end subroutine grid_allocate

!=======================================================================================
! Dealocate grid arrays.
!=======================================================================================
subroutine grid_deallocate(this)
  class(t_grid), intent(inout) :: this
  integer :: ierr

  ! Shared memory deallocation.
  call MPI_Win_free(this%win_X1, ierr)
  call MPI_Win_free(this%win_Y1, ierr)
  call MPI_Win_free(this%win_Z1, ierr)
  call MPI_Win_free(this%win_X2, ierr)
  call MPI_Win_free(this%win_Y2, ierr)
  call MPI_Win_free(this%win_Z2, ierr)

  if (ierr /= 0) call exit_MPI("Shared memory deallocation error in grid_deallocate!", 0, ierr)

end subroutine grid_deallocate

!===================================================================================
! Returns grid step hx.
!===================================================================================
pure function grid_get_hx(this, i) result(res)
  class(t_grid), intent(in) :: this
  integer, intent(in) :: i
  real(kind=CUSTOM_REAL) :: res

  res = abs(this%X2(i) - this%X1(i))

end function grid_get_hx

!===================================================================================
! Returns grid step hy.
!===================================================================================
pure function grid_get_hy(this, i) result(res)
  class(t_grid), intent(in) :: this
  integer, intent(in) :: i
  real(kind=CUSTOM_REAL) :: res

  res = abs(this%Y2(i) - this%Y1(i))

end function grid_get_hy

!===================================================================================
! Returns grid step hz.
!===================================================================================
pure function grid_get_hz(this, i) result(res)
  class(t_grid), intent(in) :: this
  integer, intent(in) :: i
  real(kind=CUSTOM_REAL) :: res

  res = abs(this%Z2(i) - this%Z1(i))

end function grid_get_hz

!===================================================================================
! Returns X-coordinate of the cell center.
!===================================================================================
pure function grid_get_X_cell_center(this, i) result(res)
  class(t_grid), intent(in) :: this
  integer, intent(in) :: i
  real(kind=CUSTOM_REAL) :: res

  res = 0.5d0 * (this%X1(i) + this%X2(i))

end function grid_get_X_cell_center

!===================================================================================
! Returns Y-coordinate of the cell center.
!===================================================================================
pure function grid_get_Y_cell_center(this, i) result(res)
  class(t_grid), intent(in) :: this
  integer, intent(in) :: i
  real(kind=CUSTOM_REAL) :: res

  res = 0.5d0 * (this%Y1(i) + this%Y2(i))

end function grid_get_Y_cell_center

!===================================================================================
! Returns Z-coordinate of the cell center.
!===================================================================================
pure function grid_get_Z_cell_center(this, i) result(res)
  class(t_grid), intent(in) :: this
  integer, intent(in) :: i
  real(kind=CUSTOM_REAL) :: res

  res = 0.5d0 * (this%Z1(i) + this%Z2(i))

end function grid_get_Z_cell_center

!===================================================================================
! Returns the cell volume.
!===================================================================================
pure function grid_get_cell_volume(this, i) result(res)
  class(t_grid), intent(in) :: this
  integer, intent(in) :: i
  real(kind=CUSTOM_REAL) :: res

  res = abs((this%X2(i) - this%X1(i)) * &
            (this%Y2(i) - this%Y1(i)) * &
            (this%Z2(i) - this%Z1(i)))

end function grid_get_cell_volume

!================================================================================================
! Get the minimum X-coordinate of the grid.
!================================================================================================
pure function grid_get_Xmin(this) result(res)
  class(t_grid), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = minval(this%X1)
end function grid_get_Xmin

!================================================================================================
! Get the maximum X-coordinate of the grid.
!================================================================================================
pure function grid_get_Xmax(this) result(res)
  class(t_grid), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = maxval(this%X2)
end function grid_get_Xmax

!================================================================================================
! Get the minimum Y-coordinate of the grid.
!================================================================================================
pure function grid_get_Ymin(this) result(res)
  class(t_grid), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = minval(this%Y1)
end function grid_get_Ymin

!================================================================================================
! Get the maximum Y-coordinate of the grid.
!================================================================================================
pure function grid_get_Ymax(this) result(res)
  class(t_grid), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = maxval(this%Y2)
end function grid_get_Ymax

!================================================================================================
! Get the minimum Z-coordinate of the grid.
!================================================================================================
pure function grid_get_Zmin(this) result(res)
  class(t_grid), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = minval(this%Z1)
end function grid_get_Zmin

!================================================================================================
! Get the maximum Z-coordinate of the grid.
!================================================================================================
pure function grid_get_Zmax(this) result(res)
  class(t_grid), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = maxval(this%Z2)
end function grid_get_Zmax

!============================================================================
! Initialize the gradient grid from a full model grid.
! The full model grid only needs to be allocated on master CPU.
!============================================================================
subroutine grad_grid_init(this, grid, myrank)
  class(t_grad_grid), intent(inout) :: this
  type(t_grid), intent(in) :: grid
  integer, intent(in) :: myrank
  integer :: i, j, k, ind, ierr

  this%nx = grid%nx
  this%ny = grid%ny
  this%nz = grid%nz

  ierr = 0

  allocate(this%dX(this%nx), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%dY(this%ny), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%dZ(this%nz), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in grad_grid_init!", myrank, ierr)

  ! Build the grid from a full model grid on master CPU.
  ! Note that we assume structured grid, i.e., dX depends only on i, dY depends only on j, and dZ depends only on k.
  if (myrank == 0) then
    do i = 1, this%nx
      ind = this%get_ind(i, 1, 1)
      this%dX(i) = grid%get_hx(ind)
    enddo
    do j = 1, this%ny
      ind = this%get_ind(1, j, 1)
      this%dY(j) = grid%get_hy(ind)
    enddo
    do k = 1, this%nz
      ind = this%get_ind(1, 1, k)
      this%dZ(k) = grid%get_hz(ind)
    enddo
  endif

  ierr = 0

  ! Broadcast the grid from the master to all CPUs.
  call MPI_Bcast(this%dX, this%nx, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%dY, this%ny, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%dZ, this%nz, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("Error in MPI_Bcast in grad_grid_init!", myrank, ierr)

end subroutine grad_grid_init

!============================================================================
! Returns 1D element index based on 3D one.
! If the 3D index in out of bounds, returns index = -1.
!============================================================================
pure function grad_grid_get_ind(this, i, j, k) result(ind)
  class(t_grad_grid), intent(in) :: this
  integer, intent(in) :: i, j, k
  integer :: ind

  if (i < 1 .or. i > this%nx .or. &
      j < 1 .or. j > this%ny .or. &
      k < 1 .or. k > this%nz) then

    ! Set invalid index.
    ind = -1
    return
  endif

  ! Note: we assume the i-j-k order of cells in the model grid.
  ind = i + (j - 1) * this%nx + (k - 1) * this%nx * this%ny

end function grad_grid_get_ind

end module grid
