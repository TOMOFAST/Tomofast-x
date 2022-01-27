
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

!==========================================================================
! A class to calculate parameter gradients.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2016.
!==========================================================================
module gradient

  use global_typedefs, only: CUSTOM_REAL
  use vector
  use grid

  implicit none

  private

  character(len=4), public, parameter :: BWD_TYPE = "BWD1"
  character(len=4), public, parameter :: FWD_TYPE = "FWD1"
  character(len=4), public, parameter :: CNT_TYPE = "CNT2"
  character(len=4), public, parameter :: FWD2_TYPE = "FWD2"
  character(len=4), public, parameter :: FWD3_TYPE = "FWD3"

  type, public :: t_gradient
    private

  contains
    private

    procedure, public, nopass :: get_grad => gradient_get_grad
    procedure, public, nopass :: get_par => gradient_get_par

    procedure, public, nopass :: get_der_type => gradient_get_der_type

  end type t_gradient

contains

!==============================================================================================
! Map derivative type to local types.
!==============================================================================================
pure function gradient_get_der_type(der_type) result(res)
  integer, intent(in) :: der_type

  character(len=4) :: res

  if (der_type == 0) then
    res = BWD_TYPE
  else if (der_type == 1) then
    res = FWD_TYPE
  else if (der_type == 2) then
    res = CNT_TYPE
  else if (der_type == 3) then
    res = FWD2_TYPE
  else if (der_type == 6) then
    res = FWD2_TYPE
  else if (der_type == 7) then
    res = FWD3_TYPE
  else
    res = FWD_TYPE
  endif

end function gradient_get_der_type

!==============================================================================================
! Returns the model gradient at pixel (i, j, k).
!==============================================================================================
recursive function gradient_get_grad(val, grid, i, j, k, der_type) result(res)
  real(kind=CUSTOM_REAL), intent(in) :: val(:)
  type(t_grid), intent(in) :: grid
  integer, intent(in) :: i, j, k
  character(len=4), intent(in) :: der_type

  type(t_vector) :: res, res_low
  integer :: p, t
  real(kind=CUSTOM_REAL) :: c(4)
  logical :: do_x, do_y, do_z
  integer :: ind
  type(t_gradient) :: grad

  ind = grid%get_ind(i, j, k)

  if (der_type == BWD_TYPE) then
  ! Backward difference scheme. O(h).
    res%x = (grad%get_par(val, grid, i, j, k) - grad%get_par(val, grid, i - 1, j, k)) / grid%get_hx(ind)
    res%y = (grad%get_par(val, grid, i, j, k) - grad%get_par(val, grid, i, j - 1, k)) / grid%get_hy(ind)
    res%z = (grad%get_par(val, grid, i, j, k) - grad%get_par(val, grid, i, j, k - 1)) / grid%get_hz(ind)

  else if (der_type == FWD_TYPE) then
  ! Forward difference scheme. O(h).
    res%x = (grad%get_par(val, grid, i + 1, j, k) - grad%get_par(val, grid, i, j, k)) / grid%get_hx(ind)
    res%y = (grad%get_par(val, grid, i, j + 1, k) - grad%get_par(val, grid, i, j, k)) / grid%get_hy(ind)
    res%z = (grad%get_par(val, grid, i, j, k + 1) - grad%get_par(val, grid, i, j, k)) / grid%get_hz(ind)

  else if (der_type == CNT_TYPE) then
  ! Central difference scheme. O(h^2).
    res%x = (grad%get_par(val, grid, i + 1, j, k) - grad%get_par(val, grid, i - 1, j, k)) / 2.d0 / grid%get_hx(ind)
    res%y = (grad%get_par(val, grid, i, j + 1, k) - grad%get_par(val, grid, i, j - 1, k)) / 2.d0 / grid%get_hy(ind)
    res%z = (grad%get_par(val, grid, i, j, k + 1) - grad%get_par(val, grid, i, j, k - 1)) / 2.d0 / grid%get_hz(ind)

  else if (der_type == FWD2_TYPE) then
  ! Forward difference scheme using three points. O(h^2).
    ! Use O(h) formula if there are less than two cells away from boundary.
    if (i >= grid%nx - 1) then
      res%x = (grad%get_par(val, grid, i + 1, j, k) - grad%get_par(val, grid, i, j, k)) / grid%get_hx(ind)
    else
      res%x = (- 1.d0 * grad%get_par(val, grid, i + 2, j, k) + 4.d0 * grad%get_par(val, grid, i + 1, j, k) &
               - 3.d0 * grad%get_par(val, grid, i, j, k)) / grid%get_hx(ind) / 2.d0
    endif

    if (j >= grid%ny - 1) then
      res%y = (grad%get_par(val, grid, i, j + 1, k) - grad%get_par(val, grid, i, j, k)) / grid%get_hy(ind)
    else
      res%y = (- 1.d0 * grad%get_par(val, grid, i, j + 2, k) + 4.d0 * grad%get_par(val, grid, i, j + 1, k) &
               - 3.d0 * grad%get_par(val, grid, i, j, k)) / grid%get_hy(ind) / 2.d0
    endif

    if (k >= grid%nz - 1) then
      res%z = (grad%get_par(val, grid, i, j, k + 1) - grad%get_par(val, grid, i, j, k)) / grid%get_hz(ind)
    else
      res%z = (- 1.d0 * grad%get_par(val, grid, i, j, k + 2) + 4.d0 * grad%get_par(val, grid, i, j, k + 1) &
               - 3.d0 * grad%get_par(val, grid, i, j, k)) / grid%get_hz(ind) / 2.d0
    endif

  else if (der_type == FWD3_TYPE) then
  ! Forward difference scheme using four points. (Derived from central difference).
  ! This is actually a central difference scheme at u[i] of fourth order,
  ! and it should be OK to use it at position u[i+1/2], i.e., shifted by 1/2 cell,
  ! considering as forward difference of third order at u[i]. See Eq. (27) in Ref. [1].
  ! [1] J. Virieux, V. Etienne, V. Cruz-Atienza, R. Brossier, E. Chaljub, et al.. Modelling Seismic
  !     Wave Propagation for Geophysical Imaging. Seismic Waves - Research and Analysis, Masaki
  !     Kanao, 253-304, Chap.13, 2012, 978-953-307-944-8.
  !
  ! u'[i + 1/2] = (- u[i + 2] + 27 u[i + 1] - 27 u[i] + u[i - 1]) / (24 dx).

    c(1) = - 1.d0
    c(2) = + 27.d0
    c(3) = - 27.d0
    c(4) = + 1.d0

    do_x = .false.
    do_y = .false.
    do_z = .false.

    ! Identify if we are close to boundary. If yes, then we will be using a lower order formula.
    if (i + 2 <= grid%nx .and. i - 1 >= 1) do_x = .true.
    if (j + 2 <= grid%ny .and. j - 1 >= 1) do_y = .true.
    if (k + 2 <= grid%nz .and. k - 1 >= 1) do_z = .true.

    res = 0.d0

    do p = 1, 4
      t = 3 - p ! = 2, 1, 0, -1

      if (do_x) res%x = res%x + c(p) * grad%get_par(val, grid, i + t, j, k)
      if (do_y) res%y = res%y + c(p) * grad%get_par(val, grid, i, j + t, k)
      if (do_z) res%z = res%z + c(p) * grad%get_par(val, grid, i, j, k + t)
    enddo

    if (do_x) res%x = res%x / 24.d0 / grid%get_hx(ind)
    if (do_y) res%y = res%y / 24.d0 / grid%get_hy(ind)
    if (do_z) res%z = res%z / 24.d0 / grid%get_hz(ind)

    ! When we are close to the boundary then use low order finite difference.
    res_low = grad%get_grad(val, grid, i, j, k, FWD_TYPE)

    if (.not. do_x) res%x = res_low%x
    if (.not. do_y) res%y = res_low%y
    if (.not. do_z) res%z = res_low%z

  else
    print *, "Unknown derivative type in gradient_get_grad!, der_type =", der_type
    stop
  endif

end function gradient_get_grad

!==================================================================================
! Returns a parameter value by its 3D index, applying the boundary conditions.
!==================================================================================
function gradient_get_par(val, grid, i, j, k) result(res)
  real(kind=CUSTOM_REAL), intent(in) :: val(:)
  type(t_grid), intent(in) :: grid
  integer, intent(in) :: i, j, k

  integer :: ib, jb, kb, index
  real(kind=CUSTOM_REAL) :: res

  if (i < 0 .or. i > grid%nx + 1 .or. &
      j < 0 .or. j > grid%ny + 1 .or. &
      k < 0 .or. k > grid%nz + 1) then
    print *, 'Wrong index in cross_gradient_get_par! nx, ny, nz, i, j, k =', &
             grid%nx, grid%ny, grid%nz, i, j, k
    stop
  endif

  ib = i
  jb = j
  kb = k

  ! Set right-side boundary conditions (for the central and forward difference scheme).
!  if (i == grid%nx + 1) ib = grid%nx
!  if (j == grid%ny + 1) jb = grid%ny
!  if (k == grid%nz + 1) kb = grid%nz
!
!  ! Set left-side boundary (for the central and backward difference scheme).
!  if (i == 0) ib = 1
!  if (j == 0) jb = 1
!  if (k == 0) kb = 1

  if (i == grid%nx + 1 .or. j == grid%ny + 1 .or. k == grid%nz + 1) then
    res = 0.d0
    return
  endif

  if (i == 0 .or. j == 0 .or. k == 0) then
    res = 0.d0
    return
  endif

  ! Calculate 1D (linear) index.
  index = grid%get_ind(ib, jb, kb)

  res = val(index)

end function gradient_get_par

end module gradient
