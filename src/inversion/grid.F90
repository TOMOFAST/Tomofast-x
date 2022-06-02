
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
! Vitaliy Ogarko, UWA, CET, Australia.
!================================================================================================
module grid

  use global_typedefs
  use mpi_tools, only: exit_MPI

  implicit none

  private

  type, public :: t_grid

    ! The beginning and the ending coordinates of a 3D prism for every model element.
    real(kind=4), dimension(:), allocatable :: X1, Y1, Z1
    real(kind=4), dimension(:), allocatable :: X2, Y2, Z2

    ! 3D index of the grid element.
    integer, dimension(:), allocatable :: i_, j_, k_

    ! 1D index of the grid element (makes sense only for the full grid).
    integer, allocatable :: ind(:, :, :)

    ! Full grid dimensions.
    integer :: nx, ny, nz

  contains
    private

    procedure, public, pass :: allocate => grid_allocate
    procedure, public, pass :: deallocate => grid_deallocate
    procedure, public, pass :: broadcast => grid_broadcast

    procedure, public, pass :: get_ind => grid_get_ind

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
subroutine grid_allocate(this, nx, ny, nz, myrank)
  class(t_grid), intent(inout) :: this
  integer, intent(in) :: nx, ny, nz, myrank

  integer :: nelements_total
  integer :: ierr

  this%nx = nx
  this%ny = ny
  this%nz = nz

  nelements_total = nx * ny * nz

  ierr = 0

  allocate(this%X1(nelements_total), source=0., stat=ierr)
  allocate(this%X2(nelements_total), source=0., stat=ierr)
  allocate(this%Y1(nelements_total), source=0., stat=ierr)
  allocate(this%Y2(nelements_total), source=0., stat=ierr)
  allocate(this%Z1(nelements_total), source=0., stat=ierr)
  allocate(this%Z2(nelements_total), source=0., stat=ierr)

  allocate(this%i_(nelements_total), source=0, stat=ierr)
  allocate(this%j_(nelements_total), source=0, stat=ierr)
  allocate(this%k_(nelements_total), source=0, stat=ierr)

  allocate(this%ind(nx, ny, nz), source=0, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in grid_initialize!", myrank, ierr)

end subroutine grid_allocate

!=======================================================================================
! Dealocate grid arrays.
!=======================================================================================
subroutine grid_deallocate(this)
  class(t_grid), intent(inout) :: this

  deallocate(this%X1)
  deallocate(this%X2)
  deallocate(this%Y1)
  deallocate(this%Y2)
  deallocate(this%Z1)
  deallocate(this%Z2)

  deallocate(this%i_)
  deallocate(this%j_)
  deallocate(this%k_)

  deallocate(this%ind)

end subroutine grid_deallocate

!==============================================================================================
! Broadcasts grid arrays from master CPU to all.
!==============================================================================================
subroutine grid_broadcast(this, myrank)
  class(t_grid), intent(inout) :: this
  integer, intent(in) :: myrank

  integer :: nelements_total
  integer :: ierr

  nelements_total = this%nx * this%ny * this%nz

  ierr = 0

  call MPI_Bcast(this%X1, nelements_total, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%X2, nelements_total, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%Y1, nelements_total, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%Y2, nelements_total, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%Z1, nelements_total, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%Z2, nelements_total, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%i_, nelements_total, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%j_, nelements_total, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%k_, nelements_total, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%ind, nelements_total, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

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

end module grid
