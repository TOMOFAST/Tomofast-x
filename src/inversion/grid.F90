
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
! A class to work with grids for forward and inverse problems.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!================================================================================================
module grid

  use global_typedefs
  use mpi_tools, only: exit_MPI

  implicit none

  private

  type, public :: t_grid

    ! The beginning and the ending coordinates of a 3D prism for every model element.
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: X1, Y1, Z1
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: X2, Y2, Z2

    ! 3D index of the grid element.
    integer, dimension(:), allocatable :: i_, j_, k_

    ! Grid size (local).
    integer :: nelements

    ! 1D index of the grid element (makes sense only for the full grid).
    integer, allocatable :: ind(:, :, :)

    ! Full grid dimensions.
    integer :: nx, ny, nz

  contains
    private

    procedure, public, pass :: initialize => grid_initialize
    procedure, public, pass :: broadcast => grid_broadcast

    procedure, public, pass :: get_ind => grid_get_ind
    procedure, public, pass :: get_k_top => grid_get_k_top

    procedure, public, pass :: get_hx => grid_get_hx
    procedure, public, pass :: get_hy => grid_get_hy
    procedure, public, pass :: get_hz => grid_get_hz

    procedure, public, pass :: get_X_cell_center => grid_get_X_cell_center
    procedure, public, pass :: get_Y_cell_center => grid_get_Y_cell_center
    procedure, public, pass :: get_Z_cell_center => grid_get_Z_cell_center

    procedure, public, pass :: get_cell_volume => grid_get_cell_volume

  end type t_grid

contains

!=======================================================================================
! Allocate grid arrays.
!=======================================================================================
subroutine grid_initialize(this, nelements, nx, ny, nz, myrank)
  class(t_grid), intent(inout) :: this
  integer, intent(in) :: nelements, nx, ny, nz, myrank
  integer :: ierr

  this%nx = nx
  this%ny = ny
  this%nz = nz

  this%nelements  = nelements

  ierr = 0

  if (.not. allocated(this%X1)) allocate(this%X1(this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%X2)) allocate(this%X2(this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%Y1)) allocate(this%Y1(this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%Y2)) allocate(this%Y2(this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%Z1)) allocate(this%Z1(this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%Z2)) allocate(this%Z2(this%nelements), source=0._CUSTOM_REAL, stat=ierr)

  if (.not. allocated(this%i_)) allocate(this%i_(this%nelements), source=0, stat=ierr)
  if (.not. allocated(this%j_)) allocate(this%j_(this%nelements), source=0, stat=ierr)
  if (.not. allocated(this%k_)) allocate(this%k_(this%nelements), source=0, stat=ierr)

  if (.not. allocated(this%ind)) allocate(this%ind(nx, ny, nz), source=0, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in grid_initialize!", myrank, ierr)

end subroutine grid_initialize

!==============================================================================================
! Broadcasts grid arrays from master CPU to all.
!==============================================================================================
subroutine grid_broadcast(this, myrank)
  class(t_grid), intent(inout) :: this
  integer, intent(in) :: myrank
  integer :: ierr

  ierr = 0

  call MPI_Bcast(this%X1, this%nelements, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%X2, this%nelements, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%Y1, this%nelements, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%Y2, this%nelements, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%Z1, this%nelements, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%Z2, this%nelements, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%i_, this%nelements, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%j_, this%nelements, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%k_, this%nelements, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%ind, this%nx * this%ny * this%nz, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("Error in MPI_Bcast in grid_broadcast!", myrank, ierr)

end subroutine grid_broadcast

!============================================================================
! Returns 1D element index based on 3D one.
! If the 3D index in out of bounds, returns index = -1.
!============================================================================
pure function grid_get_ind(this, i, j, k) result(index)
  class(t_grid), intent(in) :: this
  integer, intent(in) :: i, j, k

  integer :: index

  if (i < 1 .or. i > this%nx .or. &
      j < 1 .or. j > this%ny .or. &
      k < 1 .or. k > this%nz) then

    ! Set invalid index.
    index = -1
    return
  endif

  index = this%ind(i, j, k)

end function grid_get_ind

!=======================================================================================
! Returns k-index that corresponds to the top layer.
! Assume z-axis directs from top to bottom, i.e., smallest z-coordinate at the top.
!=======================================================================================
function grid_get_k_top(this) result(res)
  class(t_grid), intent(in) :: this
  integer :: res

  integer, save :: k_calculated = - 1
  integer :: k, k_top, ind_top

  if (k_calculated == - 1) then
    ! Find a k-index corresponding to the top layer.
    k_top = 1
    ind_top = this%ind(1, 1, k_top)
    do k = 2, this%nz
      if (this%Z1(this%ind(1, 1, k)) < this%Z1(ind_top)) then
        k_top = k
        ind_top = this%ind(1, 1, k)
      endif
    enddo

    k_calculated = k_top
    res = k_top
  else
    res = k_calculated
  end if

end function grid_get_k_top

!===================================================================================
! Returns grid step hx.
!===================================================================================
pure function grid_get_hx(this) result(res)
  class(t_grid), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = abs(this%X2(1) - this%X1(1))

end function grid_get_hx

!===================================================================================
! Returns grid step hy.
!===================================================================================
pure function grid_get_hy(this) result(res)
  class(t_grid), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = abs(this%Y2(1) - this%Y1(1))

end function grid_get_hy

!===================================================================================
! Returns grid step hz.
!===================================================================================
pure function grid_get_hz(this) result(res)
  class(t_grid), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = abs(this%Z2(1) - this%Z1(1))

end function grid_get_hz

!===================================================================================
! Returns X-coordinate of the cell center.
!===================================================================================
pure function grid_get_X_cell_center(this, cell_index) result(res)
  class(t_grid), intent(in) :: this
  integer, intent(in) :: cell_index
  real(kind=CUSTOM_REAL) :: res

  res = 0.5d0 * (this%X1(cell_index) + this%X2(cell_index))

end function grid_get_X_cell_center

!===================================================================================
! Returns Y-coordinate of the cell center.
!===================================================================================
pure function grid_get_Y_cell_center(this, cell_index) result(res)
  class(t_grid), intent(in) :: this
  integer, intent(in) :: cell_index
  real(kind=CUSTOM_REAL) :: res

  res = 0.5d0 * (this%Y1(cell_index) + this%Y2(cell_index))

end function grid_get_Y_cell_center

!===================================================================================
! Returns Z-coordinate of the cell center.
!===================================================================================
pure function grid_get_Z_cell_center(this, cell_index) result(res)
  class(t_grid), intent(in) :: this
  integer, intent(in) :: cell_index
  real(kind=CUSTOM_REAL) :: res

  res = 0.5d0 * (this%Z1(cell_index) + this%Z2(cell_index))

end function grid_get_Z_cell_center

!===================================================================================
! Returns the cell volume.
!===================================================================================
pure function grid_get_cell_volume(this, cell_index) result(res)
  class(t_grid), intent(in) :: this
  integer, intent(in) :: cell_index
  real(kind=CUSTOM_REAL) :: res

  res = abs((this%X2(cell_index) - this%X1(cell_index)) * &
            (this%Y2(cell_index) - this%Y1(cell_index)) * &
            (this%Z2(cell_index) - this%Z1(cell_index)))

end function grid_get_cell_volume

end module grid
