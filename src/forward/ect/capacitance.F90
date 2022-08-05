
!========================================================================
!
!                          T o m o f a s t - x
!                        -----------------------
!
!              Main authors: Vitaliy Ogarko, Roland Martin,
!                   Jeremie Giraud, Dimitri Komatitsch.
!
! (c) CNRS, France, and University of Western Australia. January 2018
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

!=========================================================================================
! Capacitance (array capacity) data are computed directly at each electrode
! by integrating the normal gradient of the potential on each electrode.
!=========================================================================================
module capacitance

  use global_typedefs
  use parameters_ect
  use utils

  implicit none

  private

  public :: capacitance_computed_directly

  private :: validate_capacitance
  private :: capacitance_2d_surf
  !private :: capacitance_3d_surf
  private :: capacitance_plane

contains

!=========================================================================================
! Calculate capacitance on the vertical 2D surface,
! e.g., electrodes (two-sides) or sidewalls (one-side).
!=========================================================================================
subroutine capacitance_2d_surf(nr,ntheta,nz,nzlocal,phi,theta,r,z,&
                               permit0,permit_air,permit_isolated_tube,ir,&
                               itheta_start,itheta_end,iz_start,iz_end,myrank,capacity)

  integer, intent(in) :: nr,ntheta,nz,nzlocal
  real(kind=CUSTOM_REAL), intent(in) :: phi(0:nr+1,0:ntheta+1,0:nzlocal+1)
  real(kind=CUSTOM_REAL), intent(in) :: theta(0:ntheta+1)
  real(kind=CUSTOM_REAL), intent(in) :: r(0:2*nr+1)
  real(kind=CUSTOM_REAL), intent(in) :: z(0:nz+1)

  real(kind=CUSTOM_REAL), intent(in) :: permit0,permit_air,permit_isolated_tube
  integer, intent(in) :: ir
  integer, intent(in) :: itheta_start,itheta_end
  integer, intent(in) :: iz_start,iz_end
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: capacity

  real(kind=CUSTOM_REAL) :: dtheta,capac_loc,capac,dz
  integer :: i,j,k,k1

  capacity = 0._CUSTOM_REAL

  i = ir
  do k=1,nzlocal
    k1=myrank*nzlocal+k

    if (k1 >= iz_start .and. k1 <= iz_end) then

      capac = 0._CUSTOM_REAL

      ! Use trapezoidal integration rule (with varying step).
      dz = getstep(k1,1,iz_start,iz_end,nz+1,z)

      do j=itheta_start,itheta_end

        if (i == nr+1) then
        ! Capacitance on the R3 (only inner side).
          capac_loc = &
            permit0*permit_air*&
            der(phi(i,j,k),phi(i-1,j,k),phi(i-2,j,k),phi(i-3,j,k),-(r(2*i-1)-r(2*i-3)))
        else
        ! Capacitance on the R2 (both sides).
          capac_loc = &
            permit0*permit_isolated_tube*&
            der(phi(i,j,k),phi(i-1,j,k),phi(i-2,j,k),phi(i-3,j,k),-(r(2*i-1)-r(2*i-3)))&
            -permit0*permit_air*&
            der(phi(i,j,k),phi(i+1,j,k),phi(i+2,j,k),phi(i+3,j,k),(r(2*i+1)-r(2*i-1)))
        endif

        ! Use trapezoidal integration rule (with varying step).
        dtheta = getstep(j,1,itheta_start,itheta_end,ntheta+1,theta)

        !print *, 'dtheta=', k, j, dtheta*180._CUSTOM_REAL/PI

        capac_loc = capac_loc*dtheta*r(2*i-1)*dz

        capac = capac + capac_loc
      enddo ! loop on j

      capacity = capacity + capac

    endif

  enddo ! loop on k

end subroutine capacitance_2d_surf

!===========================================================================================
! Calculate capacitance on the top/bottom walls.
! plane_type=1 for the bottom wall, plane_type=2 for the top wall.
!===========================================================================================
subroutine capacitance_plane(nr, ntheta, nz, nzlocal, phi, theta, r, z, permit, &
                             permit0, permit_air, permit_isolated_tube, i1, i2, &
                             myrank, capacity, plane_type)

  integer, intent(in) :: nr, ntheta, nz, nzlocal
  real(kind=CUSTOM_REAL), intent(in) :: phi(0:nr+1, 0:ntheta+1, 0:nzlocal+1)
  real(kind=CUSTOM_REAL), intent(in) :: theta(0:ntheta+1)
  real(kind=CUSTOM_REAL), intent(in) :: r(0:2*nr+1)
  real(kind=CUSTOM_REAL), intent(in) :: z(0:nz+1)
  real(kind=CUSTOM_REAL), intent(in) :: permit(0:, 0:, 0:)

  real(kind=CUSTOM_REAL), intent(in) :: permit0, permit_air, permit_isolated_tube
  integer, intent(in) :: i1, i2, plane_type
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: capacity

  real(kind=CUSTOM_REAL) :: dtheta, perm, capac_loc, normal, dr
  integer :: i_start, i_end, j_start, j_end
  integer :: i, j, k, k1

  capacity = 0._CUSTOM_REAL

  i_start = 1
  i_end = nr

  j_start = 1
  j_end = ntheta+1

  do k=0,nzlocal+1,nzlocal+1
    k1=myrank*nzlocal+k

    if ((plane_type == 1 .and. k1 == 0) .or. &
        (plane_type == 2 .and. k1 == nz+1)) then

      do i=i_start,i_end

        ! Use trapezoidal integration rule (with varying step).
        dr = getstep(2*i-1,2,2*i_start-1,2*i_end-1,2*nr+1,r)

        do j=j_start,j_end

          ! Use trapezoidal integration rule (with varying step).
          dtheta = getstep(j,1,j_start,j_end,ntheta+1,theta)

          if (i < i1) then
            perm = permit0 * permit(i, j, k)
          else if (i == i1) then
            perm = permit0 * avg2(permit(i, j, k), permit_isolated_tube)
          else if (i > i1 .and. i < i2) then
            perm = permit0 * permit_isolated_tube
          else if (i == i2) then
            perm = permit0 * avg2(permit_isolated_tube, permit_air)
          else if (i > i2) then
            perm = permit0 * permit_air
          endif

          if (plane_type == 1) then
            normal = 1.0_CUSTOM_REAL
            capac_loc = -perm*der(phi(i,j,k),phi(i,j,k+1),phi(i,j,k+2),phi(i,j,k+3),z(k+1)-z(k))
          else if (plane_type == 2) then
            normal = -1.0_CUSTOM_REAL
            capac_loc = -perm*der(phi(i,j,k),phi(i,j,k-1),phi(i,j,k-2),phi(i,j,k-3),z(k-1)-z(k))
          endif

          ! TODO: Confirm area calculation.
          if (i < i_end) then
            capac_loc = capac_loc*normal*0.5_CUSTOM_REAL*dtheta*((r(2*i-1)+dr)**2-r(2*i-1)**2)
          else if (i == i_end) then
            capac_loc = capac_loc*normal*0.5_CUSTOM_REAL*dtheta*(r(2*i-1)**2-(r(2*i-1)-dr)**2)
          endif

          capacity = capacity + capac_loc

        enddo ! loop on j
      enddo ! loop on i

    endif

  enddo ! loop on k

end subroutine capacitance_plane

!===================================================================================================
! Compute capacitance by direct integration over the electrode/walls.
!===================================================================================================
subroutine capacitance_computed_directly(phi, theta,r, z, permit, i1, i2, &
                                         par, capacity, ielectrode, elecs, guards, myrank, nbproc)
  type(t_parameters_ect), intent(in) :: par
  type(t_electrode), intent(in) :: elecs(:)
  type(t_guards), intent(in) :: guards

  real(kind=CUSTOM_REAL), intent(in) :: phi(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: theta(0:)
  real(kind=CUSTOM_REAL), intent(in) :: r(0:)
  real(kind=CUSTOM_REAL), intent(in) :: z(0:)
  real(kind=CUSTOM_REAL), intent(in) :: permit(0:, 0:, 0:)
  integer, intent(in) :: i1, i2
  integer, intent(in) :: myrank, nbproc
  integer, intent(in) :: ielectrode

  real(kind=CUSTOM_REAL), intent(out) :: capacity(:, :)

  ! Local variables.
  integer :: ierr, jelectrode
  real(kind=CUSTOM_REAL) :: capacity_loc, capacity_glob
  real(kind=CUSTOM_REAL) :: capacity_low, capacity_up, capacity_walls
  real(kind=CUSTOM_REAL) :: capacity_bottom, capacity_top
  real(kind=CUSTOM_REAL) :: capacity_sum

  ! Low guard capacitance calculation.
  call capacitance_2d_surf(par%dims%nr, par%dims%ntheta, par%dims%nz, par%dims%nzlocal, &
            phi, theta, r, z, &
            par%permit0, par%permit_air, par%permit_isolated_tube, i2, &
            1, par%dims%ntheta + 1, 1, guards%lower_izmax, myrank, capacity_low)

!  call capacitance_3d_surf(par%dims%nr,par%dims%ntheta,par%dims%nz,par%dims%nzlocal,&
!                           phi,theta,r,z,&
!                           par%permit0,par%permit_air,par%permit_isolated_tube,i2,&
!                           1,par%dims%ntheta+1,1,guards%lower_izmax,&
!                           myrank,capacity_low,2)

  capacity_loc = capacity_low

  if (nbproc > 1) then
    call mpi_allreduce(capacity_loc,capacity_glob,1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD,ierr)
  else
    capacity_glob = capacity_loc
  endif
  capacity_low = capacity_glob

  ! Upper guard capacitance calculation.
  call capacitance_2d_surf(par%dims%nr, par%dims%ntheta, par%dims%nz, par%dims%nzlocal, &
                           phi, theta, r, z, &
                           par%permit0, par%permit_air, par%permit_isolated_tube, i2, &
                           1, par%dims%ntheta + 1, guards%upper_izmin, par%dims%nz, myrank, capacity_up)

!  call capacitance_3d_surf(par%dims%nr,par%dims%ntheta,par%dims%nz,par%dims%nzlocal,&
!                           phi,theta,r,z,&
!                           par%permit0,par%permsit_air,par%permit_isolated_tube,i2,&
!                           1,par%dims%ntheta+1,guards%upper_izmin,par%dims%nz,&
!                           myrank,capacity_up,3)

  capacity_loc = capacity_up

  if (nbproc > 1) then
    call mpi_allreduce(capacity_loc,capacity_glob,1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  else
    capacity_glob = capacity_loc
  endif
  capacity_up = capacity_glob

  ! Bottom wall capacitance calculation.
  call capacitance_plane(par%dims%nr, par%dims%ntheta, par%dims%nz, par%dims%nzlocal, &
                         phi, theta, r, z, permit, &
                         par%permit0, par%permit_air, par%permit_isolated_tube, i1, i2, &
                         myrank, capacity_bottom, 1)

  capacity_loc = capacity_bottom

  if (nbproc > 1) then
    call mpi_allreduce(capacity_loc,capacity_glob,1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD,ierr)
  else
    capacity_glob = capacity_loc
  endif
  capacity_bottom = capacity_glob

  ! Top wall capacitance calculation.
  call capacitance_plane(par%dims%nr, par%dims%ntheta, par%dims%nz, par%dims%nzlocal, &
                         phi, theta, r, z, permit, &
                         par%permit0, par%permit_air, par%permit_isolated_tube, i1, i2, &
                         myrank, capacity_top, 2)

  capacity_loc = capacity_top

  if (nbproc > 1) then
    call mpi_allreduce(capacity_loc,capacity_glob,1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD,ierr)
  else
    capacity_glob = capacity_loc
  endif
  capacity_top = capacity_glob

  ! Side wall capacitance calculation.
  call capacitance_2d_surf(par%dims%nr, par%dims%ntheta, par%dims%nz, par%dims%nzlocal, &
                           phi, theta, r, z, &
                           par%permit0, par%permit_air, par%permit_isolated_tube, par%dims%nr + 1, &
                           1, par%dims%ntheta + 1, 1, par%dims%nz, myrank, capacity_walls)

  capacity_loc = capacity_walls

  if (nbproc > 1) then
    call mpi_allreduce(capacity_loc,capacity_glob,1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  else
    capacity_glob = capacity_loc
  endif
  capacity_walls = capacity_glob

  ! Electrodes elementary capacitances calculation.
  do jelectrode = 1, par%nel

!    if (nint(par%sens%space_electrodes) > 0) then
!      call capacitance_3d_surf(par%dims%nr,par%dims%ntheta,par%dims%nz,par%dims%nzlocal,&
!                               phi,theta,r,z,&
!                               par%permit0,par%permit_air,par%permit_isolated_tube,i2,&
!                               elecs(jelectrode)%ithetamin,elecs(jelectrode)%ithetamax,&
!                               elecs(jelectrode)%izmin,elecs(jelectrode)%izmax,&
!                               myrank,capacity(ielectrode,jelectrode),1)
!    else
      call capacitance_2d_surf(par%dims%nr, par%dims%ntheta, par%dims%nz, par%dims%nzlocal, &
                               phi, theta, r, z, &
                               par%permit0, par%permit_air, par%permit_isolated_tube, i2, &
                               elecs(jelectrode)%ithetamin, elecs(jelectrode)%ithetamax, &
                               elecs(jelectrode)%izmin, elecs(jelectrode)%izmax, &
                               myrank, capacity(ielectrode, jelectrode))
!    endif
  enddo

  do jelectrode = 1, par%nel
    capacity_loc = capacity(ielectrode, jelectrode)
    if (nbproc > 1) then
      call mpi_allreduce(capacity_loc, capacity_glob, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    else
      capacity_glob = capacity_loc
    endif
    capacity(ielectrode, jelectrode) = capacity_glob

    ! Sanity check.
    call validate_capacitance(capacity, ielectrode, jelectrode)
  enddo

! END CAPACITANCE CALCULATIONS

#ifndef SUPPRESS_OUTPUT
  if (myrank == 0) then
    print *, 'CAPACITANCE, electrode = ', ielectrode

    capacity_sum = 0._CUSTOM_REAL
    do jelectrode = 1, par%nel
      print *, ielectrode, jelectrode, capacity(ielectrode, jelectrode)
      capacity_sum = capacity_sum + capacity(ielectrode, jelectrode)
    enddo

    print *, 'Lower guard capacitance =', capacity_low
    print *, 'Upper guard capacitance =', capacity_up
    print *, '  Side wall capacitance =', capacity_walls
    print *, 'Bottom wall capacitance =', capacity_bottom
    print *, '   Top wall capacitance =', capacity_top
    print *, 'Capacitance sum (only electrodes) = ', capacity_sum, &
            '  Scaled by self-capa = ', capacity_sum / abs(capacity(ielectrode, ielectrode))

    capacity_sum = capacity_sum + capacity_low + capacity_up + capacity_walls + capacity_bottom + capacity_top
    print *, 'Capacitance sum TOTAL = ', capacity_sum, &
            '  scaled by self-capa = ', capacity_sum / abs(capacity(ielectrode, ielectrode))

#ifdef USE_FLUSH6
    flush(6)
#endif
  endif
#endif

end subroutine capacitance_computed_directly

!===================================================================================================
! Validate calculated capacitance values: Cij < 0 for i /= j and Cii > 0.
!===================================================================================================
subroutine validate_capacitance(capacity, ielectrode, jelectrode)
  real(kind=CUSTOM_REAL), intent(in) :: capacity(:, :)
  integer, intent(in) :: ielectrode, jelectrode

  if (ielectrode /= jelectrode) then
    if (capacity(ielectrode, jelectrode) >= 0.d0) then
      print *, 'Sanity check failed: non-negative capacitance Cij (for i /= j)! i, j, Cij =', &
        ielectrode, jelectrode, capacity(ielectrode, jelectrode)
      stop
    endif
  else
    if (capacity(ielectrode, jelectrode) <= 0.d0) then
      print *, 'Sanity check failed: non-positive self capacitance Cii! i, Cij =', &
        ielectrode, capacity(ielectrode, jelectrode)
      stop
    endif
  endif

end subroutine validate_capacitance

!=========================================================================================
! NOT USED (Sanity check validate_capacitance (Cij < 0 for i /= j) is not passing).
!
! Calculate integral of potential gradient on the 3D surface around the electrode/guard.
! surf_type=1 for electrodes,
! surf_type=2 for the lower guard,
! surf_type=3 for the upper guard.
!=========================================================================================
!subroutine capacitance_3d_surf(nr,ntheta,nz,nzlocal,phi,theta,r,z,&
!                               permit0,permit_air,permit_isolated_tube,i2,&
!                               itheta_start,itheta_end,iz_start,iz_end,myrank,capacity,surf_type)
!  integer, intent(in) :: nr,ntheta,nz,nzlocal
!  real(kind=CUSTOM_REAL), intent(in) :: phi(0:nr+1,0:ntheta+1,0:nzlocal+1)
!  real(kind=CUSTOM_REAL), intent(in) :: theta(0:ntheta+1)
!  real(kind=CUSTOM_REAL), intent(in) :: r(0:2*nr+1)
!  real(kind=CUSTOM_REAL), intent(in) :: z(0:nz+1)
!
!  real(kind=CUSTOM_REAL), intent(in) :: permit0,permit_air,permit_isolated_tube
!  integer, intent(in) :: i2,itheta_start,itheta_end,iz_start,iz_end
!  integer, intent(in) :: myrank
!  integer, intent(in) :: surf_type
!
!  real(kind=CUSTOM_REAL), intent(out) :: capacity
!
!  real(kind=CUSTOM_REAL) :: dtheta,capac_loc,dz,dr
!  real(kind=CUSTOM_REAL) :: perm,normal
!  integer :: i_start,i_end,j_start,j_end,k_start,k_end
!  integer :: i,j,k,k1
!
!  ! DISABLED: Sanity check validate_capacitance (Cij < 0 for i /= j) is not passing.
!  print *, 'Subroutine capacitance_3d_surf() should not be used before it is fixed!'
!  stop
!
!  i_start = i2 - 1
!  i_end   = i2 + 1
!
!  if (surf_type == 1) then
!    j_start = itheta_start - 1
!    j_end   = itheta_end + 1
!  else
!    ! For guards do exactly full circle 360 degrees, from 1 to ntheta+1.
!    j_start = itheta_start
!    j_end   = itheta_end
!  endif
!
!  k_start = iz_start - 1
!  k_end   = iz_end + 1
!
!  capacity = 0._CUSTOM_REAL
!
!  ! Vertical surfaces (parallel to the electrode/guard).
!  do k=1,nzlocal
!    k1=myrank*nzlocal+k
!
!    if (k1 >= k_start .and. k1 <= k_end) then
!      ! Use trapezoidal integration rule (with varying step).
!      dz = getstep(k1,1,k_start,k_end,nz+1,z)
!
!      do j=j_start,j_end
!        ! Use trapezoidal integration rule (with varying step).
!        dtheta = getstep(j,1,j_start,j_end,ntheta+1,theta)
!
!        ! Inner (i_start) and outer surfaces (i_end).
!        do i=i_start,i_end,i_end-i_start
!
!          if (i < i2) then
!            normal = -1._CUSTOM_REAL
!            perm = permit0*permit_isolated_tube
!          else if (i > i2) then
!            normal = 1._CUSTOM_REAL
!            perm = permit0*permit_air
!          endif
!
!          capac_loc = -perm*(phi(i+1,j,k)-phi(i-1,j,k))/(r(2*i+1)-r(2*i-3))
!
!          capac_loc = capac_loc*normal*dtheta*r(2*i-1)*dz
!
!          capacity = capacity + capac_loc
!        enddo
!      enddo
!    endif
!  enddo
!
!  ! Top & bottom surfaces.
!  do k=1,nzlocal
!    k1=myrank*nzlocal+k
!
!    if ((surf_type == 1 .and. (k1 == k_start .or. k1 == k_end)) .or. &
!        (surf_type == 2 .and. k1 == k_end) .or. &
!        (surf_type == 3 .and. k1 == k_start)) then
!    ! Exclude bottom surface for the lower guard, and top surface for the upper guard.
!
!      do j=j_start,j_end
!        ! Use trapezoidal integration rule (with varying step).
!        dtheta = getstep(j,1,j_start,j_end,ntheta+1,theta)
!
!        do i=i_start,i_end
!
!          ! Use trapezoidal integration rule (with varying step).
!          dr = getstep(2*i-1,2,2*i_start-1,2*i_end-1,2*nr+1,r)
!
!          if (i < i2) then
!            perm = permit0*permit_isolated_tube
!          else if (i > i2) then
!            perm = permit0*permit_air
!          else
!            perm = permit0*avg2(permit_isolated_tube,permit_air)
!          endif
!
!          if (k1 == k_start) then
!            normal = -1._CUSTOM_REAL
!          else if (k1 == k_end) then
!            normal = 1._CUSTOM_REAL
!          endif
!
!          capac_loc = -perm*(phi(i,j,k+1)-phi(i,j,k-1))/(z(k1+1)-z(k1-1))
!
!          ! TODO: Confirm area calculation.
!          if (i < i_end) then
!            capac_loc = capac_loc*normal*0.5_CUSTOM_REAL*dtheta*((r(2*i-1)+dr)**2-r(2*i-1)**2)
!          else if (i == i_end) then
!            capac_loc = capac_loc*normal*0.5_CUSTOM_REAL*dtheta*(r(2*i-1)**2-(r(2*i-1)-dr)**2)
!          endif
!
!!          if (i == i_start) then
!!            capac_loc = capac_loc*normal*0.5_CUSTOM_REAL*dtheta*(r(2*i+1)**2-r(2*i-1)**2)*0.5_CUSTOM_REAL
!!          else if (i == i_end) then
!!            capac_loc = capac_loc*normal*0.5_CUSTOM_REAL*dtheta*(r(2*i-1)**2-r(2*i-3)**2)*0.5_CUSTOM_REAL
!!          else
!!            capac_loc = capac_loc*normal*0.5_CUSTOM_REAL*dtheta*(r(2*i+1)**2-r(2*i-3)**2)*0.5_CUSTOM_REAL
!!          endif
!
!          capacity = capacity + capac_loc
!        enddo
!      enddo
!    endif
!  enddo
!
!
!  if (surf_type == 1) then
!  ! Only for electrodes. Not for guards, because they are periodic.
!
!    ! Side surfaces (perpendicular to electrode).
!    do k=1,nzlocal
!      k1=myrank*nzlocal+k
!
!      if (k1 >= k_start .and. k1 <= k_end) then
!        ! Use trapezoidal integration rule (with varying step).
!        dz = getstep(k1,1,k_start,k_end,nz+1,z)
!
!        do j=j_start,j_end,j_end-j_start
!
!          if (j == j_start) then
!            normal = -1._CUSTOM_REAL
!          else if (j == j_end) then
!            normal = 1._CUSTOM_REAL
!          endif
!
!          do i=i_start,i_end
!            ! Use trapezoidal integration rule (with varying step).
!            dr = getstep(2*i-1,2,2*i_start-1,2*i_end-1,2*nr+1,r)
!
!            if (i < i2) then
!              perm = permit0*permit_isolated_tube
!            else if (i > i2) then
!              perm = permit0*permit_air
!            else
!              perm = permit0*avg2(permit_isolated_tube,permit_air)
!            endif
!
!            if (j-1 < 1) then
!              capac_loc = perm*(phi(i,j+1,k)-phi(i,j-1+ntheta,k))/(theta(j+1)-theta(j-1+ntheta)-2._CUSTOM_REAL*PI)/r(2*i-1)
!            else
!              capac_loc = perm*(phi(i,j+1,k)-phi(i,j-1,k))/(theta(j+1)-theta(j-1))/r(2*i-1)
!            endif
!
!            capac_loc = capac_loc*normal*dr*dz
!
!            capacity = capacity + capac_loc
!          enddo
!        enddo ! j
!      endif
!    enddo
!  endif
!
!end subroutine capacitance_3d_surf

end module capacitance
