
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
! A class to calculate depth weights for the sensitivity kernel.
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
  use parallel_tools

  implicit none

  private

  public :: calculate_depth_weight
  public :: apply_local_depth_weighting

  private :: normalize_depth_weight
  private :: calc_depth_weight_pixel

contains

!===================================================================================
! Calculates the depth weight for sensitivity kernel.
!===================================================================================
subroutine calculate_depth_weight(par, iarr, grid_full, data, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  type(t_grid), intent(in) :: grid_full
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
  real(kind=CUSTOM_REAL) :: dist, mindist
  integer :: nsmaller, p

  if (myrank == 0) print *, 'Calculating the depth weight, type = ', par%depth_weighting_type

  ! The number of elements on CPUs with rank smaller than myrank.
  nsmaller = get_nsmaller(par%nelements, myrank, nbproc)

  !--------------------------------------------------------------------------------
  ! Calculate the normalized depth weight.
  !--------------------------------------------------------------------------------
  if (par%depth_weighting_type == 1) then
  ! Depth weighting.

    ! Use empirical function 1/(z+z0)**(power/2).
    do i = 1, par%nelements
      ! Full grid index.
      p = nsmaller + i
      iarr%column_weight(i) = calc_depth_weight_pixel(grid_full, par%depth_weighting_power, par%Z0, p, myrank)
    enddo

  else if (par%depth_weighting_type == 2) then
  ! Distance weighting.
  ! Based on Eq.(19) in Li & Oldenburg, Joint inversion of surface and three‐component borehole magnetic data, GEOPHYSICS 65, 540–552 (2000).

    ! A small constant for integral validity.
    R0 = 0.1d0 ! 0.1 meter

    ! A factor to move the cell corner points inside a cell.
    dfactor = 0.25d0

    do i = 1, par%nelements
      ! Full grid index.
      p = nsmaller + i

      ! Cell colume.
      dVj = grid_full%get_cell_volume(p)

      ! Shifts to make the integral points to lie inside the cell volume.
      dhx = dfactor * abs(grid_full%X2(p) - grid_full%X1(p))
      dhy = dfactor * abs(grid_full%Y2(p) - grid_full%Y1(p))
      dhz = dfactor * abs(grid_full%Z2(p) - grid_full%Z1(p))

      wr = 0.d0

      do j = 1, par%ndata

        ! The squared 1D distances along each cell dimension.
        distX_sq(1) = (grid_full%X1(p) + dhx - data%X(j))**2.d0
        distY_sq(1) = (grid_full%Y1(p) + dhy - data%Y(j))**2.d0
        distZ_sq(1) = (grid_full%Z1(p) + dhz - data%Z(j))**2.d0

        distX_sq(2) = (grid_full%X2(p) - dhx - data%X(j))**2.d0
        distY_sq(2) = (grid_full%Y2(p) - dhy - data%Y(j))**2.d0
        distZ_sq(2) = (grid_full%Z2(p) - dhz - data%Z(j))**2.d0

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

      iarr%column_weight(i) = (1.d0 / sqrt(dVj)) * wr**(par%depth_weighting_beta / 4.d0)
    enddo ! cells loop

  else if (par%depth_weighting_type == 3) then
    ! A small constant for division validity.
    R0 = 0.01d0

    do i = 1, par%nelements
      ! Full grid index.
      p = nsmaller + i

      ! Calculate the minimum distance from a cell center to data.
      mindist = 1.d30
      do j = 1, par%ndata
        dist = sqrt((grid_full%get_X_cell_center(p) - data%X(j))**2.d0 + &
                    (grid_full%get_Y_cell_center(p) - data%Y(j))**2.d0 + &
                    (grid_full%get_Z_cell_center(p) - data%Z(j))**2.d0)

        if (dist < mindist) then
          mindist = dist
        endif
      enddo ! data loop

      iarr%column_weight(i) = sqrt(1.d0 / (mindist + R0)**par%depth_weighting_power)
    enddo ! cells loop

  else
    call exit_MPI("Not known depth weight type!", myrank, par%depth_weighting_type)
  endif

  !--------------------------------------------------------------------------------
  ! Scale the sensitivity kernel with the cell volume.
  !--------------------------------------------------------------------------------
  do i = 1, par%nelements
      ! Full grid index.
      p = nsmaller + i

      iarr%column_weight(i) = iarr%column_weight(i) * sqrt(grid_full%get_cell_volume(p))
  enddo

  ! Normalize the depth weight.
  call normalize_depth_weight(iarr%column_weight, myrank, nbproc)

  !--------------------------------------------------------------------------------
  ! Calculate the matrix column weight.
  !--------------------------------------------------------------------------------

  ! This condition essentially leads to the system:
  !
  ! | S W^{-1} | d(Wm)
  ! |  alpha I |
  !
  do i = 1, par%nelements
    if (iarr%column_weight(i) /= 0.d0) then
      iarr%column_weight(i) = 1.d0 / iarr%column_weight(i)
    else
      call exit_MPI("Zero damping weight! Exiting.", myrank, 0)
    endif
  enddo

  if (myrank == 0) print *, 'Finished calculating the depth weight.'

end subroutine calculate_depth_weight

!===================================================================================
! Calculates the depth weight for a pixel using empirical function.
!===================================================================================
function calc_depth_weight_pixel(grid, power, Z0, i, myrank) result(weight)
  type(t_grid), intent(in) :: grid
  real(kind=CUSTOM_REAL), intent(in) :: power
  integer, intent(in) :: i, myrank
  real(kind=CUSTOM_REAL), intent(in) :: Z0
  real(kind=CUSTOM_REAL) :: weight

  real(kind=CUSTOM_REAL) :: depth

  ! Depth to the middle of the voxel.
  depth = grid%get_Z_cell_center(i)

  if (depth + Z0 > 0.d0) then
    weight = (depth + Z0)**(- power / 2.d0)
  else
    print *, depth, Z0
    call exit_MPI("Error: non-positive depth in calc_depth_weight_pixel!", myrank, 0)
  endif

end function calc_depth_weight_pixel

!===================================================================================
! Normalizes the depth weight.
!===================================================================================
subroutine normalize_depth_weight(depth_weight, myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), intent(inout) :: depth_weight(:)

  integer :: ierr
  real(kind=CUSTOM_REAL) :: norm, norm_glob

  ! Find the maximum value in the depth weight array.
  norm = maxval(depth_weight)
  if (nbproc > 1) then
    call mpi_allreduce(norm, norm_glob, 1, CUSTOM_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD, ierr)
    norm = norm_glob
  endif

  if (norm /= 0) then
    ! Normalize.
    depth_weight = depth_weight / norm
  else
    call exit_MPI("Zero depth weight norm! Exiting.", myrank, 0)
  endif

end subroutine normalize_depth_weight

!===================================================================================
! Precondition the depth weight with local weights read from file.
!===================================================================================
subroutine apply_local_depth_weighting(par, column_weight, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(inout) :: column_weight(par%nelements)

  character(len=256) :: msg
  integer :: ierr
  integer :: i, ind
  integer :: nsmaller, nelements_total, nelements_read
  real(kind=CUSTOM_REAL) :: local_weight

  if (par%apply_local_weight > 0) then

    open(10, file=trim(par%local_weight_file), status='old', action='read', iostat=ierr, iomsg=msg)
    if (ierr /= 0) call exit_MPI("Error in opening the local weight file! path=" &
                   //par%local_weight_file//" iomsg="//msg, myrank, ierr)

    read(10, *, iostat=ierr) nelements_read

    nelements_total = par%nx * par%ny * par%nz

    ! Sanity check.
    if (nelements_total /= nelements_read) &
      call exit_MPI("The local weight is not correctly defined!", myrank, 0)

    ! The number of elements on CPUs with rank smaller than myrank.
    nsmaller = get_nsmaller(par%nelements, myrank, nbproc)

    ! Reading local weights from file.
    do i = 1, nelements_total
      if (i > nsmaller .and. i <= nsmaller + par%nelements) then
        read(10, *, iostat=ierr) local_weight

        if (ierr /= 0) &
          call exit_MPI("Problem with reading the local weight!", myrank, ierr)

        ind = i - nsmaller

        ! Apply local weight.
        if (local_weight /= 0.d0) then
          ! Divide the column weight which makes the depth weight multiplied.
          column_weight(ind) = column_weight(ind) / local_weight
        else
          column_weight(ind) = 0.d0
        endif
      else
        ! Skip the line.
        read(10, *)
      endif

      if (i > nsmaller + par%nelements) then
        exit
      endif
    enddo
  endif

end subroutine apply_local_depth_weighting

end module weights_gravmag
