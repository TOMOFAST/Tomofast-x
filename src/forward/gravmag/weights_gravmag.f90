
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

    procedure, private, nopass :: calculate_depth_weight
    procedure, private, nopass :: calculate_depth_weight_sensit
    procedure, private, nopass :: normalize_depth_weight

  end type t_weights

contains

!===================================================================================
! Calculates the weights for inversion.
!===================================================================================
subroutine weights_calculate(par, iarr, xdata, ydata, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  real(kind=CUSTOM_REAL), intent(in) :: xdata(:), ydata(:)
  integer, intent(in) :: myrank, nbproc
  type(t_inversion_arrays), intent(inout) :: iarr
  integer :: i, ierr

  if (par%depth_weighting_type /= 1) then
    print *, "Currently only the depth weight type = 1 is supported!"
    stop
  endif

  if (myrank == 0) print *, 'Calculating the depth weight...'

  !--------------------------------------------------------------------------------
  ! Calculate the damping weight as the normalized depth weight.
  !--------------------------------------------------------------------------------

  if (par%depth_weighting_type == 1) then

    ! Method I: use empirical function 1/(z+z0)**(beta/2).
    do i = 1, par%nelements
      iarr%damping_weight(i) = calculate_depth_weight(iarr%model%grid, par%beta, par%Z0, i, myrank)
    enddo

  else if (par%depth_weighting_type == 2) then

    if (myrank == 0) print *, 'Error: Not supported case!'
    stop

    ! Method II: use only sensitivity values directly below the data (i.e., z-column).
    ! Calculate damping weight using sensitivity kernel.
    call calculate_depth_weight_sensit(iarr%model%grid, xdata, ydata, iarr%sensitivity, iarr%damping_weight, &
                                       iarr%nelements, iarr%ndata, myrank)

  else if (par%depth_weighting_type == 3) then

    ! Temporary disable this case as wavelet compression requires to weight the sensitivity line before compression.
    if (myrank == 0) print *, 'Error: Not supported (due to wavelet compression)!'
    stop

    ! Method III: scale model by the integrated sensitivities, see
    ! [1] (!!) Yaoguo Li, Douglas W. Oldenburg., Joint inversion of surface and three-component borehole magnetic data, 2000.
    ! [2] Portniaguine and Zhdanov (2002).
    ! For discussion on different weightings see also:
    !   [1] M. Pilkington, Geophysics, vol. 74, no. 1, 2009.
    !   [2] F. Cella and M. Fedi, Geophys. Prospecting, 2012, 60, 313-336.

    ! Integrated sensitivity matrix (diagonal).
    call iarr%matrix_sensit%get_integrated_sensit(iarr%damping_weight)
    iarr%damping_weight = sqrt(iarr%damping_weight)

  else
    call exit_MPI("Unknown depth weighting type!", myrank, 0)
  endif

  ! Normalize the depth weight.
  call normalize_depth_weight(iarr%damping_weight, myrank, nbproc)

  !--------------------------------------------------------------------------------
  ! Calculate the matrix column weight.
  !--------------------------------------------------------------------------------

  ! This condition essentially leads to the system:
  !
  ! | S W^{-1} | d(Wm)
  ! |    I     |
  !
  do i = 1, par%nelements
    if (iarr%damping_weight(i) /= 0.d0) then
      iarr%column_weight(i) = 1.d0 / iarr%damping_weight(i)
    else
      !iarr%column_weight(i) = 1.d0
      call exit_MPI("Zero damping weight! Exiting.", myrank, 0)
    endif
  enddo

  ! This condition essentially leads to the system:
  !
  ! | S | dm
  ! | W |
  !iarr%column_weight = 1._CUSTOM_REAL


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

!========================================================================================
! Calculates the damping weight using sensitivity kernel below the data location.
!========================================================================================
subroutine calculate_depth_weight_sensit(grid, xdata, ydata, sensit, damping_weight, &
                                         nelements, ndata, myrank)
  type(t_grid), intent(in) :: grid
  real(kind=CUSTOM_REAL), intent(in) :: xdata(:), ydata(:)
  real(kind=CUSTOM_REAL), intent(in) :: sensit(:, :)
  integer, intent(in) :: nelements, ndata, myrank
  real(kind=CUSTOM_REAL), intent(out) :: damping_weight(:)

  integer :: p, i, idata

  ! Loop over local elements.
  do p = 1, nelements

    ! Search for data corresponding to the element:
    !   the data (X, Y) position is inside a grid-cell (X, Y) position.
    idata = 0
    do i = 1, ndata
      if (xdata(i) >= grid%X1(p) .and. xdata(i) <= grid%X2(p) .and. &
          ydata(i) >= grid%Y1(p) .and. ydata(i) <= grid%Y2(p)) then
        idata = i
        exit
      endif
    enddo

    if (idata == 0) then
      call exit_MPI("Error: Not found data corresponding to the pixel!", myrank, 0)
    endif

    if (sensit(p, idata) >= 0) then
      damping_weight(p) = sqrt(sensit(p, idata))
    else
      call exit_MPI("Error: Negative sensitivity, cannot calculate damping weight!", myrank, 0)
    endif
  enddo

end subroutine calculate_depth_weight_sensit

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
