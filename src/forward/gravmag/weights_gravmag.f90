
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

!===============================================================================================
! A class to calculate weights for the sensitivity matrix.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!===============================================================================================
module weights_gravmag

  use global_typedefs
  use parameters_gravmag
  use parameters_inversion
  use inversion_arrays
  use grid
  use mpi_tools, only: exit_MPI
  use data_gravmag

  implicit none

  private

  type, public :: t_weights
    private

  contains
    private

    procedure, public, nopass :: calculate => weights_calculate
    procedure, public, nopass :: normalize_depth_weight

    procedure, private, nopass :: calculate_depth_weight

  end type t_weights

contains

!===================================================================================
! Calculates the weights for inversion.
!===================================================================================
subroutine weights_calculate(par, iarr, data, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  integer, intent(in) :: myrank, nbproc
  type(t_data), intent(in) :: data
  type(t_inversion_arrays), intent(inout) :: iarr

  integer :: i, j
  integer :: ii, jj, kk, ind
  real(kind=CUSTOM_REAL) :: distX_sq(2), distY_sq(2), distZ_sq(2)
  real(kind=CUSTOM_REAL) :: dhx, dhy, dhz, dfactor
  real(kind=CUSTOM_REAL) :: Rij(8), R0
  real(kind=CUSTOM_REAL) :: dVj
  real(kind=CUSTOM_REAL) :: integral, wr

  if (myrank == 0) print *, 'Calculating the depth weight, type = ', par%depth_weighting_type

  !--------------------------------------------------------------------------------
  ! Calculate the normalized depth weight.
  !--------------------------------------------------------------------------------
  if (par%depth_weighting_type == 1) then
  ! Depth weighting.

    ! Use empirical function 1/(z+z0)**(beta/2).
    do i = 1, par%nelements
      iarr%damping_weight(i) = calculate_depth_weight(iarr%model%grid, par%depth_weighting_power, par%Z0, i, myrank)
    enddo

  else if (par%depth_weighting_type == 2) then
  ! Distance weighting.
  ! The implementation is based on the Eq.(10) in https://www.eoas.ubc.ca/ubcgif/iag/sftwrdocs/grav3d/grav3d-manual.pdf

    ! A small constant for integral validity.
    R0 = 0.1d0 ! 0.1 meter

    ! A factor to move the cell corner points inside a cell.
    dfactor = 0.25d0

    do i = 1, par%nelements
      ! Cell colume.
      dVj = iarr%model%grid%get_cell_volume(i)

      ! Shifts to make the integral points to lie inside the cell volume.
      dhx = dfactor * abs(iarr%model%grid%X2(i) - iarr%model%grid%X1(i))
      dhy = dfactor * abs(iarr%model%grid%Y2(i) - iarr%model%grid%Y1(i))
      dhz = dfactor * abs(iarr%model%grid%Z2(i) - iarr%model%grid%Z1(i))

      wr = 0.d0

      do j = 1, par%ndata

        ! The squared 1D distances along each cell dimension.
        distX_sq(1) = (iarr%model%grid%X1(i) + dhx - data%X(j))**2.d0
        distY_sq(1) = (iarr%model%grid%Y1(i) + dhy - data%Y(j))**2.d0
        distZ_sq(1) = (iarr%model%grid%Z1(i) + dhz - data%Z(j))**2.d0

        distX_sq(2) = (iarr%model%grid%X2(i) - dhx - data%X(j))**2.d0
        distY_sq(2) = (iarr%model%grid%Y2(i) - dhy - data%Y(j))**2.d0
        distZ_sq(2) = (iarr%model%grid%Z2(i) - dhz - data%Z(j))**2.d0

        ! Calculate the distance from 8 points inside a cell to the data location.
        ind = 0
        do ii = 1, 2
          do jj = 1, 2
            do kk = 1, 2
              ind = ind + 1
              Rij(ind) = sqrt(distX_sq(ii) + distY_sq(jj) + distZ_sq(kk))
            enddo
          enddo
        enddo

        integral = 0.d0
        do ind = 1, 8
          integral = integral + 1.d0 / (Rij(ind) + R0)**par%depth_weighting_power
        enddo
        integral = integral * dVj / 8.d0

        wr = wr + integral**2.d0

      enddo ! data loop

      iarr%damping_weight(i) = (1.d0 / sqrt(dVj)) * wr**(1.d0 / 4.d0)
    enddo ! cells loop

  else
    call exit_MPI("Not known depth weight type!", myrank, par%depth_weighting_type)
  endif

  ! Normalize the depth weight.
  call normalize_depth_weight(iarr%damping_weight, myrank, nbproc)

  !--------------------------------------------------------------------------------
  ! Calculate the matrix column weight.
  !--------------------------------------------------------------------------------

  ! This condition essentially leads to the system:
  !
  ! | S W^{-1} | d(Wm)
  ! |  alpha I |
  !
  do i = 1, par%nelements
    if (iarr%damping_weight(i) /= 0.d0) then
      iarr%column_weight(i) = 1.d0 / iarr%damping_weight(i)
    else
      call exit_MPI("Zero damping weight! Exiting.", myrank, 0)
    endif
  enddo

  if (myrank == 0) print *, 'Finished calculating the depth weight.'

end subroutine weights_calculate

!===================================================================================
! Calculates the depth weight for a pixel using empirical function.
!===================================================================================
function calculate_depth_weight(grid, beta, Z0, i, myrank) result(weight)
  type(t_grid), intent(in) :: grid
  real(kind=CUSTOM_REAL), intent(in) :: beta
  integer, intent(in) :: i, myrank
  real(kind=CUSTOM_REAL), intent(in) :: Z0
  real(kind=CUSTOM_REAL) :: weight

  real(kind=CUSTOM_REAL) :: depth

  ! Depth to the middle of the voxel.
  depth = grid%get_Z_cell_center(i)

  if (depth + Z0 > 0.d0) then
    weight = (depth + Z0)**(- beta / 2.d0)
  else
    print *, depth
    print *, Z0
    call exit_MPI("Error: non-positive depth in calculate_depth_weight!", myrank, 0)
  endif

end function calculate_depth_weight

!===================================================================================
! Normalizes the depth weight.
!===================================================================================
subroutine normalize_depth_weight(damping_weight, myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), intent(inout) :: damping_weight(:)

  integer :: ierr
  real(kind=CUSTOM_REAL) :: norm, norm_glob

  ! Find the maximum value in the depth weight array.
  norm = maxval(damping_weight)
  if (nbproc > 1) then
    call mpi_allreduce(norm, norm_glob, 1, CUSTOM_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD, ierr)
    norm = norm_glob
  endif

  if (norm /= 0) then
    ! Normalize.
    damping_weight = damping_weight / norm
  else
    call exit_MPI("Zero damping weight norm! Exiting.", myrank, 0)
  endif

end subroutine normalize_depth_weight

end module weights_gravmag
