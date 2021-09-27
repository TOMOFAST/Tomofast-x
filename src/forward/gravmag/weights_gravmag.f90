
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
! A class to calculate weights for sensitivity matrix and damping.
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
subroutine weights_calculate(par, iarr, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  integer, intent(in) :: myrank, nbproc
  type(t_inversion_arrays), intent(inout) :: iarr
  integer :: i

  ! NOTE:
  ! We use the power-function depth weighting with wavelet compression to first estimate the nnz.
  ! The actual depth weight in this case will be computed in calculate_sensitivity().

  if (par%compression_type == 0 .and. par%depth_weighting_type == 3) then
    call exit_MPI("For the sensitivity-based depth weighting activate the wavelet compression!", myrank, 0)
  endif

  if (myrank == 0) print *, 'Calculating the depth weight.'

  !--------------------------------------------------------------------------------
  ! Calculate the normalized depth weight.
  !--------------------------------------------------------------------------------

  ! Use empirical function 1/(z+z0)**(beta/2).
  do i = 1, par%nelements
    iarr%damping_weight(i) = calculate_depth_weight(iarr%model%grid, par%beta, par%Z0, i, myrank)
  enddo

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
