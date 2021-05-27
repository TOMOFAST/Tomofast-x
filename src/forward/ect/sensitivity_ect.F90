
!========================================================================
!
!                    T O M O F A S T X  Version 1.0
!                  ----------------------------------
!
!              Main authors: Vitaliy Ogarko, Roland Martin,
!                   Jeremie Giraud, Dimitri Komatitsch.
! CNRS, France, and University of Western Australia.
! (c) CNRS, France, and University of Western Australia. January 2018
!
! This software is a computer program whose purpose is to perform
! capacitance, gravity, magnetic, or joint gravity and magnetic tomography.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

module sensitivity_ect

  use utils
  use global_typedefs
  use parameters_ect
  use data_ect
  use geometry, only: compute_shifted_grid_coordinates
  use paraview_ect

  implicit none

  private

  public :: calculate_sensitivity
  public :: visualise_sensitivity
  public :: calculate_weights
  public :: calculate_capa_sens

  private :: calculate_gradients

contains

!===============================================================================================
! Calculate potential gradients, needed for the sensitivity.
!===============================================================================================
pure subroutine calculate_gradients(phisensit, r_step_inv, theta_step_inv, z_step_inv, &
                                    i, j, k, ielectrode, &
                                    gradphi_r, gradphi_theta, gradphi_z)

  real(kind=CUSTOM_REAL), intent(in) :: phisensit(0:, 0:, 0:, :)
  ! Inverse gradient steps in each direction.
  real(kind=CUSTOM_REAL), intent(in) :: r_step_inv, theta_step_inv, z_step_inv
  integer, intent(in) :: i, j, k, ielectrode

  real(kind=CUSTOM_REAL), intent(out) :: gradphi_r, gradphi_theta, gradphi_z

  real(kind=CUSTOM_REAL) :: phi_iel_1, phi_iel_2, phi_iel_3, phi_iel_4
  real(kind=CUSTOM_REAL) :: phi_iel_5, phi_iel_6, phi_iel_7, phi_iel_8
  real(kind=CUSTOM_REAL) :: phi_T, phi_B, phi_N, phi_S, phi_E, phi_W

  phi_iel_1 = phisensit(i, j, k, ielectrode)
  phi_iel_2 = phisensit(i - 1, j, k, ielectrode)
  phi_iel_3 = phisensit(i - 1, j + 1, k, ielectrode)
  phi_iel_4 = phisensit(i, j + 1, k, ielectrode)
  phi_iel_5 = phisensit(i, j, k - 1, ielectrode)
  phi_iel_6 = phisensit(i - 1, j, k - 1, ielectrode)
  phi_iel_7 = phisensit(i - 1, j + 1, k - 1, ielectrode)
  phi_iel_8 = phisensit(i, j + 1, k - 1, ielectrode)

  phi_T = avg4(phi_iel_1, phi_iel_2, phi_iel_3, phi_iel_4)
  phi_B = avg4(phi_iel_5, phi_iel_6, phi_iel_7, phi_iel_8)

  gradphi_z = (phi_T - phi_B) * z_step_inv

  phi_N = avg4(phi_iel_1, phi_iel_4, phi_iel_8, phi_iel_5)
  phi_S = avg4(phi_iel_2, phi_iel_3, phi_iel_7, phi_iel_6)

  gradphi_r = (phi_N - phi_S) * r_step_inv

  phi_E = avg4(phi_iel_3, phi_iel_4, phi_iel_8, phi_iel_7)
  phi_W = avg4(phi_iel_1, phi_iel_2, phi_iel_6, phi_iel_5)

  gradphi_theta = (phi_E - phi_W) * theta_step_inv

end subroutine calculate_gradients

!===================================================================================================================
! Calculates capacitance values based on the sensitivity.
! Note: only those electrode pairs are calculated which are allowed in accept_electrode_pair(), i.e.,
! by default some pairs data are skipped for inversion.
!===================================================================================================================
subroutine calculate_capa_sens(par, sensit, capacity, permit, myrank, nbproc)
  type(t_parameters_ect), intent(in) :: par
  ! Sensitivity kernel.
  real(kind=CUSTOM_REAL), intent(in) :: sensit(:, :)
  ! Permittivity values.
  real(kind=CUSTOM_REAL), intent(in) :: permit(0:, 0:, 0:)
  integer, intent(in) :: myrank, nbproc

  ! Capacitance data.
  real(kind=CUSTOM_REAL), intent(out) :: capacity(:, :)

  real(kind=CUSTOM_REAL) :: capacity_loc, capacity_glob
  integer :: i, j, k, l, l1, el
  integer :: ielectrode, jelectrode
  integer :: ierr

  ! Note: the same loop as in calculate_sensitivity.
  l = 0
  l1 = 0
  do jelectrode = 1, par%nel
    do ielectrode = 1, par%nel
      l = l + 1

      ! Not all data are used for inversion.
      if (par%accept_electrode_pair(ielectrode, jelectrode)) l1 = l1 + 1

      capacity_loc = 0._CUSTOM_REAL
      el = 0

      do k = 1, par%dims%nzlocal + (myrank + 1) / nbproc
        do j = 1, par%dims%ntheta
          do i = 2, par%dims%nr + 1

            if (par%accept_electrode_pair(ielectrode, jelectrode)) then
              el = el + 1
              capacity_loc = capacity_loc + permit(i, j, k) * sensit(el, l1)
            endif

          enddo ! i
        enddo ! j
      enddo ! k

      ! Compute global capacitance.
      if (nbproc > 1) then
        call mpi_allreduce(capacity_loc, capacity_glob, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      else
        capacity_glob = capacity_loc
      endif
      capacity(ielectrode, jelectrode) = capacity_glob

      !if (myrank == 0) print *, l, ielectrode, jelectrode, capacity(ielectrode, jelectrode)

    enddo ! ielectrode
  enddo ! jelectrode

  ! Writing to file the capacitance data.
  call write_data_to_file(par, capacity, 'capa-with-sens2', myrank)

end subroutine

!===========================================================================================
! Calculates sensitivity values (preliminary, not scaled).
!===========================================================================================
pure subroutine calculate_sensitivity(par, sensit, phisensit, cell_volume, &
                                      r, theta, z, myrank, nbproc)
  type(t_parameters_ect), intent(in) :: par
  real(kind=CUSTOM_REAL), intent(in) :: phisensit(0:, 0:, 0:, :)
  real(kind=CUSTOM_REAL), intent(in) :: cell_volume(:, :, :)
  ! Steps along r, theta and z.
  real(kind=CUSTOM_REAL), intent(in) :: r(0:)
  real(kind=CUSTOM_REAL), intent(in) :: theta(0:)
  real(kind=CUSTOM_REAL), intent(in) :: z(0:)
  integer, intent(in) :: myrank, nbproc

  ! Sensitivity kernel.
  real(kind=CUSTOM_REAL), intent(out) :: sensit(:, :)

  ! Local variables.
  real(kind=CUSTOM_REAL) :: dtheta, gphi_ij, sensit_pixel
  real(kind=CUSTOM_REAL) :: gradphi_theta_i, gradphi_r_i, gradphi_z_i
  real(kind=CUSTOM_REAL) :: gradphi_theta_j, gradphi_r_j, gradphi_z_j
  real(kind=CUSTOM_REAL) :: r_step_inv, theta_step_inv, z_step_inv
  integer :: i, j, k, k1, l1, el
  integer :: ielectrode, jelectrode

  sensit = 0._CUSTOM_REAL

  l1 = 0
  do jelectrode = 1, par%nel
    do ielectrode = 1, par%nel

      ! Not all data are used for inversion.
      if (.not. par%accept_electrode_pair(ielectrode, jelectrode)) then
        cycle
      endif

      l1 = l1 + 1
      el = 0

      do k = 1, par%dims%nzlocal + (myrank + 1) / nbproc
        k1 = myrank * par%dims%nzlocal + k
        z_step_inv = 1._CUSTOM_REAL / (z(k1) - z(k1 - 1))

        do j = 1, par%dims%ntheta
          dtheta = abs(theta(j + 1) - theta(j))

          do i = 2, par%dims%nr + 1
            theta_step_inv = 1._CUSTOM_REAL / (r(2 * i - 2) * dtheta)
            r_step_inv = 1._CUSTOM_REAL / (r(2 * i - 1) - r(2 * i - 3))

            call calculate_gradients(phisensit, r_step_inv, theta_step_inv, z_step_inv, &
                                     i, j, k, ielectrode, &
                                     gradphi_r_i, gradphi_theta_i, gradphi_z_i)

            call calculate_gradients(phisensit, r_step_inv, theta_step_inv, z_step_inv, &
                                     i, j, k, jelectrode, &
                                     gradphi_r_j, gradphi_theta_j, gradphi_z_j)

            gphi_ij = gradphi_r_i * gradphi_r_j + gradphi_theta_i * gradphi_theta_j + gradphi_z_i * gradphi_z_j

            sensit_pixel = gphi_ij * cell_volume(i, j, k)

            ! Store sensitivity kernel and corresponding scaling factor (cell volume).
            ! Note the i-index starts at i=2, since no sensitivity data is calculated near the axis,
            ! because the gradient cannot be calculated (no potential values at the axis i=0).
            el = el + 1

            sensit(el, l1) = par%permit0 * sensit_pixel
          enddo ! i
        enddo ! j
      enddo ! k

    enddo ! ielectrode
  enddo ! jelectrode

end subroutine calculate_sensitivity

!===================================================================================================
! Calculates the column and damping weights, needed for inversion.
!===================================================================================================
pure subroutine calculate_weights(par, cell_volume, column_weight, damping_weight, myrank, nbproc)
  type(t_parameters_ect), intent(in) :: par
  real(kind=CUSTOM_REAL), intent(in) :: cell_volume(:, :, :)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), intent(out) :: damping_weight(:)
  real(kind=CUSTOM_REAL), intent(out) :: column_weight(:)

  integer :: i, j, k, k1, el

  el = 0

  ! Note: use the same loop dimensions as in calculate_sensitivity().

  do k = 1, par%dims%nzlocal + (myrank + 1) / nbproc
    k1 = myrank * par%dims%nzlocal + k
    do j = 1, par%dims%ntheta
      do i = 2, par%dims%nr + 1
        el = el + 1

        ! Scale sensitivity matrix with the cell volume.
        column_weight(el) = 1._CUSTOM_REAL / cell_volume(i, j, k)

        ! Needed to have a proper norm of the model: Integral{|m|^p dv}.
        damping_weight(el) = sqrt(cell_volume(i, j, k))

        ! Vitaliy tried to scale with cell_volume_avg to keep sens values of same order but results are the same.
        ! When doing it adjust cell_volume_avg calculation for parallel version in calculate_cell_volumes().
      enddo
    enddo
  enddo

end subroutine calculate_weights

!=============================================================================================
! Write Paraview profile for the sensitivity.
!=============================================================================================
subroutine visualise_sensitivity(par, sensit, column_weight, r, theta, z, myrank, nbproc)
  type(t_parameters_ect), intent(in) :: par
  real(kind=CUSTOM_REAL), intent(in) :: sensit(:, :)
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  real(kind=CUSTOM_REAL), intent(in) :: r(0:)
  real(kind=CUSTOM_REAL), intent(in) :: theta(0:)
  real(kind=CUSTOM_REAL), intent(in) :: z(0:)
  integer, intent(in) :: myrank,nbproc

  type(t_dimensions) :: dims
  real(kind=CUSTOM_REAL), allocatable :: sensit_pair(:, :, :)
  ! Coordinates of points in the grid for Paraview visualization.
  real(kind=CUSTOM_REAL), allocatable :: xgrid(:, :, :)
  real(kind=CUSTOM_REAL), allocatable :: ygrid(:, :, :)
  real(kind=CUSTOM_REAL), allocatable :: zgrid(:, :, :)

  logical, parameter :: use_logscale = .true.
  ! TODO: This constant is used for the grid size 36^3.
  ! TODO: For larger sized, the cell size is smaller which will affect the magnitude of scaled sensitivity values.
  ! TODO: We should divide scaling_factor (column_weight) by average cell volume to keep the same magnitude.
  real(kind=CUSTOM_REAL), parameter :: min_val_logscale = 1e-1_CUSTOM_REAL
  integer :: ielectrode,jelectrode,l
  integer :: i,j,k,el

  dims = par%dims

  allocate(sensit_pair(0:dims%nr+1, 0:dims%ntheta+1, 0:dims%nzlocal+1))
  allocate(xgrid(0:dims%nr+1, 0:dims%ntheta+1, 0:dims%nz+1))
  allocate(ygrid(0:dims%nr+1, 0:dims%ntheta+1, 0:dims%nz+1))
  allocate(zgrid(0:dims%nr+1, 0:dims%ntheta+1, 0:dims%nz+1))

  call compute_shifted_grid_coordinates(dims%nr, dims%ntheta, dims%nz, r, theta, z, xgrid, ygrid, zgrid)

  ! Initialize non-calculated sensitivity values, e.g. for i=1.
  if (use_logscale) then
    sensit_pair = log10(min_val_logscale)
  else
    sensit_pair = 0.0_CUSTOM_REAL
  endif

  l = 0
  do jelectrode = 1, par%nel
    do ielectrode = 1, par%nel

      if (par%accept_electrode_pair(ielectrode, jelectrode)) then
        l = l + 1

        ! Visualize only selected electrode pairs.
        if (ielectrode == 7 .and. jelectrode == 1) then

          el = 0
          ! Use same loop as in calculate_sensitivity().
          do k=1,dims%nzlocal + (myrank+1)/nbproc
            do j=1,dims%ntheta
              do i=2,dims%nr+1
                el = el + 1

                ! Use same index mapping as in calculate_sensitivity().
                !sensit_pair(i,j,k) = sensit(i-1,j,k,l) * column_weight(el)
                sensit_pair(i,j,k) = sensit(el, l) * column_weight(el)

                if (use_logscale) then
                  ! There are positive and negative values.
                  sensit_pair(i,j,k) = abs(sensit_pair(i,j,k))

                  ! Modify values for plotting in log-log scale, since real values start from zero.
                  if (sensit_pair(i,j,k) < min_val_logscale) sensit_pair(i,j,k) = min_val_logscale

                  sensit_pair(i,j,k) = log10(sensit_pair(i,j,k))
                endif
              enddo
            enddo
          enddo

          ! Apply periodic BCs.
          call enforce_pb(dims%nr,dims%ntheta,dims%nzlocal,sensit_pair)

          ! Write Paraview data files for visualization.
          call paraview_write_2d_profiles(myrank,"sensit2d_",ielectrode,&
                                          dims%nr,dims%ntheta,dims%nz,dims%nzlocal,sensit_pair,xgrid,ygrid,zgrid)

          call paraview_write_3d_profiles(myrank,"sensit3d_",ielectrode,&
                                          dims%nr,dims%ntheta,dims%nzlocal,sensit_pair,xgrid,ygrid,zgrid)
        endif
      endif
    enddo
  enddo

  deallocate(sensit_pair)
  deallocate(xgrid)
  deallocate(ygrid)
  deallocate(zgrid)

end subroutine visualise_sensitivity

end module sensitivity_ect

