
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

!======================================================================================================
! In the Finite Volume formulation one has the nodes of the unknowns located at the center of a cell,
! and the nodes of the fluxes located in the middle of each edge of the cell.
! Permittivities are defined at the same location at a corner of a cell.
! It is necessary to do that because one can not impose two different permittivities
! at a sharp interface of the tube (at the electrodes and at i = i1).
!======================================================================================================
module geometry

  use global_typedefs
  use parameters_ect

  implicit none

  private

  public :: get_refined_ntheta
  public :: set_geometry
  public :: compute_grid_coordinates
  public :: compute_shifted_grid_coordinates
  public :: refine_geometry
  public :: calculate_cell_volumes

  private :: calc_electrode_coordinates

contains

!========================================================
! Returns the number of points per electrode.
!========================================================
function get_npelec(ntheta, nel)
  integer, intent(in) :: ntheta, nel
  integer :: get_npelec

  get_npelec = 2*ntheta/nel/3
  return
end function get_npelec

!========================================================
! Returns the number of points per gaps.
!========================================================
function get_npgaps(ntheta, nel)
  integer, intent(in) :: ntheta, nel
  integer :: get_npgaps

  get_npgaps = ntheta/nel - get_npelec(ntheta, nel)
  return
end function get_npgaps

!========================================================
! Returns the ntheta for the refined mesh.
!========================================================
function get_refined_ntheta(ntheta, nel)
  integer, intent(in) :: ntheta, nel
  integer :: get_refined_ntheta

  get_refined_ntheta = ntheta + get_npelec(ntheta,nel)-1 + 2*(get_npgaps(ntheta,nel)+1)
  return
end function get_refined_ntheta

!=============================================================================================
! Calculate the bottom and top (z-direction) grid coordinates of the electrodes.
!=============================================================================================
subroutine calc_electrode_coordinates(nrings, nel, kguards, elecs)
  integer, intent(in) :: nrings, nel, kguards
  type(t_electrode), intent(inout) :: elecs(:)

  integer :: ielectrode
  integer :: height, iring

  ! TODO: Add (vertical) gaps between the rings.

  if (mod(2 * kguards, nrings) /= 0) then
    print *, "2 * kguards is not a multiple of nrings; exiting..."
    stop
  endif

  if (mod(nel, nrings) /= 0) then
    print *, "nel is not a multiple of nrings; exiting..."
    stop
  endif

  height = 2 * kguards / nrings

  do ielectrode = 1, nel
    ! The ring number to which this electrode belongs, e.g.,
    ! 1-12 electrodes belong to the 1st rings, 13-24 electrodes to the 2nd ring etc.
    iring = (ielectrode - 1) / (nel / nrings) + 1

    elecs(ielectrode)%izmin = kguards + 1 + (iring - 1) * height
    elecs(ielectrode)%izmax = elecs(ielectrode)%izmin + height - 1

    !print *, '---->', ielectrode, iring, elecs(ielectrode)%izmin, elecs(ielectrode)%izmax
  enddo

end subroutine calc_electrode_coordinates

!=============================================================================================
! Sets the system geometry.
!=============================================================================================
subroutine set_geometry(ilevel, par, nr, ntheta, nz, r, dr, theta, dz, z, izspace, kguards, &
                        i1, i2, elecs, guards)
  ! Parameters.
  type(t_parameters_ect), intent(in) :: par
  ! Dimensions (locally).
  integer, intent(in) :: nr, ntheta, nz, kguards
  integer, intent(in) :: ilevel
  ! theta, r and z angle, radius and size of each point, step along r and z
  real(kind=CUSTOM_REAL), intent(out) :: r(0:)
  real(kind=CUSTOM_REAL), intent(out) :: dr(0:)
  real(kind=CUSTOM_REAL), intent(out) :: theta(0:)
  real(kind=CUSTOM_REAL), intent(out) :: dz(0:)
  real(kind=CUSTOM_REAL), intent(out) :: z(0:)
  integer, intent(out) :: izspace
  integer, intent(out) :: i1, i2
  type(t_electrode), intent(out) :: elecs(:)
  type(t_guards), intent(out) :: guards

  ! Local variables --------------------------------------
  real(kind=CUSTOM_REAL) :: dtheta_elec, dtheta_gap, heiguard
  real(kind=CUSTOM_REAL) :: dr1, dr2, dr3, dz1, dz2, dz3, dz4, dz5
  integer :: i, j, k, i1_read, i2_read
  integer :: npelec, npgaps, ielec, ielectrode, jelectrode
  integer :: nr1, nr2
  integer :: igaps_on_theta,igaps_on_z
  integer :: jleft
  integer :: nelr

  !=======================================================
  ! Compute radii.
  !=======================================================

  if (par%ifixed_elecgeo == 1) then
    nr1 = int(par%sens%radiusin * par%dims%nr / par%sens%radiusoutout)
    nr2 = int((par%sens%radiusout - par%sens%radiusin) * par%dims%nr / par%sens%radiusoutout)
    i1_read = nr1-mod(nr1,8)
    i2_read = i1_read+nr2-mod(nr2,8)
    nr1 = i1_read
    nr2 = i2_read

    i1 = i1_read/2**(ilevel-1)
    i2 = i2_read/2**(ilevel-1)
    nr1 = i1
    nr2 = i2
  else
    i2 = 8*nr/9
    i1 = 3*i2/4
    nr1 = i1
    nr2 = i2
  endif

! VO VO RM RM 6 Jan 2015 fixed the mapping according to Mapping #3 in the document Mappings.pdf
! VO VO RM RM cell index i --> radius index 2i-1,
! VO VO RM RM with ghost points i=0 and i=nr+1 separated by 2dr from regular points i=1 and i=nr.
! VO VO RM RM To be consistent with matrix_init and r(2i1-1)=R1, r(2i2-1)=R2 and r(2(nr+1)-1)=R3
  dr1 = (par%sens%radiusin) / dble(2*nr1)
  dr2 = (par%sens%radiusout - par%sens%radiusin) / dble(2*nr2-2*nr1)
  dr3 = (par%sens%radiusoutout - par%sens%radiusout) / dble(2*nr-2*nr2+2)

  do i=0,2*i1-1
    r(i) = dr1+dble(i)*dr1
  enddo

  do i=2*i1,2*i2-1
    r(i) = par%sens%radiusin+dble(i-2*nr1+1)*dr2
  enddo

  do i=2*i2,2*nr+1
    r(i) = par%sens%radiusout+dble(i-2*nr2+1)*dr3
  enddo
! VO VO RM RM 6 Jan 2015 End of the fix.

   !Set equal dr everywhere.
!  do i=0,2*nr+1
!    r(i) = dble(i)*sens%radiusoutout/dble(2*nr+1)
!  enddo

! VO VO 9 Jan fixed dr(0)=r(0)-r(-1) according to radius index Mapping #3.
  dr(0) = r(0)
  do i=1,2*nr+1
    dr(i) = r(i)-r(i-1)
  enddo

  !=======================================================
  ! Compute theta angles.
  !=======================================================

  ! Number of electrodes per ring.
  nelr = par%nel / par%nrings

  igaps_on_theta = nint(par%sens%space_electrodes)
  !if (myrank == 0) print *,'ilevel, igaps_on_theta=', ilevel, igaps_on_theta

! VO VO 8 Jan 2015 Fixed to have the same angle at the tips of electrodes at different level.
! VO VO 8 Jan 2015 Added unit tests for this check marked with (***)
! VO VO 8 Jan 2015 Made the calculation robust, i.e., avoiding iterative summation in the loop.
! VO VO 14 Jan 2015 Made sure that there is at least one point on the gap, i.e., a point without Dirichlet BCs.
! VO VO So for ntheta=36 electrode points will be 1-2, 4-5, 7-8, etc.
  if (igaps_on_theta /= 0) then
    ! The number of points per electrode.
    npelec = get_npelec(ntheta, nelr)
    ! The number of points per gap, i.e., points without Dirichlet BCs.
    npgaps = get_npgaps(ntheta, nelr)
    ! Compute angles.
    dtheta_elec = (2._CUSTOM_REAL*PI/nelr - par%sens%space_electrodes*PI/180._CUSTOM_REAL) / dble(npelec-1)
    dtheta_gap = par%sens%space_electrodes*PI/180._CUSTOM_REAL / dble(npgaps+1)

    theta(0) = dtheta_gap
    do j=1,ntheta+1
      ielec = (j-1) / (npelec+npgaps)
      jleft = (j-1) - ielec*(npelec+npgaps)

      if (jleft <= npelec-1) then
        theta(j) = -2._CUSTOM_REAL*PI/nelr*ielec - dtheta_elec*jleft
      else
        theta(j) = -2._CUSTOM_REAL*PI/nelr*ielec - dtheta_elec*(npelec-1) - dtheta_gap*(jleft-npelec+1)
      endif
    enddo

    ! Calculate the first and the last indices in the theta-array for the electrodes.
    do ielectrode=1,nelr
      elecs(ielectrode)%ithetamin = (ielectrode-1)*(npelec+npgaps)+1
      elecs(ielectrode)%ithetamax = (ielectrode-1)*(npelec+npgaps)+npelec
    enddo
! VO VO 14 Jan 2015 End of the fix.

  else
    ! The number of points per electrode.
    npelec = ntheta/nelr
    ! Compute angle.
    dtheta_elec = (2._CUSTOM_REAL*PI/nelr)/dble(npelec)

    do j=0,ntheta+1
      theta(j) = dtheta_elec*dble(1-j)
    enddo

    ! Calculate the first and the last indices in the theta-array for the electrodes.
    do ielectrode=1,nelr
      elecs(ielectrode)%ithetamin = (ielectrode-1)*(npelec)+1
      elecs(ielectrode)%ithetamax = (ielectrode-1)*(npelec)+npelec
    enddo
  endif
! VO VO 8 Jan 2015 End of the fix.

  if (par%nrings > 1) then
  ! Copy the theta angles from the first to the other electrode rings.
    do ielectrode = nelr + 1, par%nel
      jelectrode = mod(ielectrode, nelr)
      if (jelectrode == 0) jelectrode = nelr

      elecs(ielectrode)%ithetamin = elecs(jelectrode)%ithetamin
      elecs(ielectrode)%ithetamax = elecs(jelectrode)%ithetamax
    enddo
  endif

  !-------------------------------------------------------------------
  ! Calculate top and bottom z-coordinates of the electrodes.
  !-------------------------------------------------------------------
  call calc_electrode_coordinates(par%nrings, par%nel, kguards, elecs)

  !=======================================================
  ! Compute height interval along z-axis, and spacing.
  !=======================================================

! VO VO 13 Jan 2015 Fixed problems of different electrode top and bottom position at different levels.
  igaps_on_z = nint(par%sens%space_elec_guards / (par%sens%heicyl/dble(nz+1)))

  if (igaps_on_z /= 0) then

    ! izspace = the number of points on the gap, i.e., those which will not get any Dirichlet BCs.
    if (par%ifixed_elecgeo == 1) then
      izspace = par%dims%nz/36 / 2**(ilevel-1)
    else
      izspace = igaps_on_z
    endif

    dz1 = (0.25_CUSTOM_REAL*par%sens%heicyl - par%sens%space_elec_guards) / dble(kguards+1-izspace-1)
    dz2 = par%sens%space_elec_guards / dble(izspace+1)
    dz3 = 0.5_CUSTOM_REAL*par%sens%heicyl / dble(3*kguards-kguards-1)
    dz4 = dz2
    dz5 = dz1

    ! lower guard.
    do k=0,kguards+1-izspace-1
      dz(k) = dz1
      z(k)  = dble(k)*dz1
    enddo

    ! lower gap.
    do k=kguards+1-izspace,kguards+1
      dz(k) = dz2
      z(k)  = (0.25_CUSTOM_REAL*par%sens%heicyl - par%sens%space_elec_guards) + dble(k-kguards+izspace)*dz2
    enddo

    ! electrode.
    do k=kguards+2,3*kguards
      dz(k) = dz3
      z(k)  = 0.25_CUSTOM_REAL*par%sens%heicyl + dble(k-kguards-1)*dz3
    enddo

    ! upper gap.
    do k=3*kguards+1,3*kguards+izspace+1
      dz(k) = dz4
      z(k)  = 0.75_CUSTOM_REAL*par%sens%heicyl + dble(k-3*kguards)*dz4
    enddo

    ! upper guard.
    do k=3*kguards+izspace+2,nz+1
      dz(k) = dz5
      z(k)  = (0.75_CUSTOM_REAL*par%sens%heicyl + par%sens%space_elec_guards) + dble(k-3*kguards-izspace-1)*dz5
    enddo
! VO VO 13 Jan 2015 End of the fix.

  else
! VO VO 9 Jan 2015 Fixed problems of different electrode top and bottom position at different levels.
    izspace = 0

    ! The guard height.
    heiguard = 0.25_CUSTOM_REAL*par%sens%heicyl

    dz1 = heiguard/dble(kguards+1)
    dz2 = (par%sens%heicyl-2._CUSTOM_REAL*heiguard) / dble(3*kguards-kguards-1)
    dz3 = dz1

    do k=0,kguards+1
      dz(k) = dz1
      z(k)  = dble(k)*dz1
    enddo

    do k=kguards+2,3*kguards
      dz(k) = dz2
      z(k)  = heiguard + dble(k-kguards-1)*dz2
    enddo

    do k=3*kguards+1,nz+1
      dz(k) = dz3
      z(k)  = (par%sens%heicyl - heiguard) + dble(k-3*kguards)*dz3
    enddo
! VO VO 9 Jan 2015 End of the fix.
  endif

  !------------------------------------------------------------------
  ! Calculate the z grid coordinates for the guards.
  guards%lower_izmax = kguards+1-izspace-1
  guards%upper_izmin = 3*kguards+izspace+1

end subroutine set_geometry

!===========================================================================================
! Compute grid coordinates.
! This grid is used for model initialization and visualization.
!===========================================================================================
subroutine compute_grid_coordinates(nr, ntheta, nz, r, theta, z, xgrid, ygrid, zgrid)
  integer, intent(in) :: nr, ntheta, nz
  real(kind=CUSTOM_REAL), intent(in) :: r(0:)
  real(kind=CUSTOM_REAL), intent(in) :: theta(0:)
  real(kind=CUSTOM_REAL), intent(in) :: z(0:)

  real(kind=CUSTOM_REAL), intent(out) :: xgrid(0:,0:,0:)
  real(kind=CUSTOM_REAL), intent(out) :: ygrid(0:,0:,0:)
  real(kind=CUSTOM_REAL), intent(out) :: zgrid(0:,0:,0:)

  ! Local variables.
  integer :: i,j,k
  real(kind=CUSTOM_REAL), allocatable :: x(:,:)
  real(kind=CUSTOM_REAL), allocatable :: y(:,:)

  allocate(x(0:2*nr+1,0:ntheta+1))
  allocate(y(0:2*nr+1,0:ntheta+1))

  ! x-y-z coordinates of r-theta points
  do j=0,ntheta+1
    do i=0,2*nr+1
      x(i,j) = r(i)*cos(theta(j))
      y(i,j) = r(i)*sin(theta(j))
    enddo
  enddo

  ! x-y-z coordinates of r-theta points 1 over 2 points
  ! DG DG and DK DK: create a cutplane in (r,theta) direction, and extend
  ! it in the z direction
  ! DG DG and DK DK: staggered grid, thus only use corner points
  ! VO VO 6 Jan 2015 fixed according to the new, correct radii indices mapping i-->2*i-1
  do k=0,nz+1
    do j=0,ntheta+1
      xgrid(0,j,k) = 0._CUSTOM_REAL
      ygrid(0,j,k) = 0._CUSTOM_REAL
      zgrid(0,j,k) = z(k)

      do i=1,nr+1
        xgrid(i,j,k) = x(2*i-1,j)
        ygrid(i,j,k) = y(2*i-1,j)
        zgrid(i,j,k) = z(k)
      enddo
    enddo
  enddo

  deallocate(x)
  deallocate(y)

end subroutine compute_grid_coordinates

!===========================================================================================
! Calculate the volume of the cells which is used to calculate the sensitivity or
! to scale the columns of the sensitivity matrix and the model.
!===========================================================================================
subroutine calculate_cell_volumes(dims, volume, volume_avg, r, theta, z, myrank, nbproc)
  type(t_dimensions), intent(in) :: dims
  real(kind=CUSTOM_REAL), intent(in) :: r(0:)
  real(kind=CUSTOM_REAL), intent(in) :: theta(0:)
  real(kind=CUSTOM_REAL), intent(in) :: z(0:)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), intent(out) :: volume(:,:,:)
  real(kind=CUSTOM_REAL), intent(out) :: volume_avg

  integer :: i,j,k,k1
  real(kind=CUSTOM_REAL) :: dz,dtheta
  integer :: n_avg

  volume_avg = 0._CUSTOM_REAL
  n_avg = 0

  do k=1,dims%nzlocal + (myrank+1)/nbproc
    k1 = myrank*dims%nzlocal+k
    dz = z(k1)-z(k1-1)

    do j=1,dims%ntheta
      dtheta = abs(theta(j+1)-theta(j))

      do i=1,dims%nr+1
        if (i == 1) then
          ! Initialize also these elements (not used). Note that r(-1) = 0 by definition.
          volume(i,j,k) = 0.5_CUSTOM_REAL * r(2*i-1)**2 * dtheta * dz
        else
          volume(i,j,k) = 0.5_CUSTOM_REAL * (r(2*i-1)**2 - r(2*i-3)**2) * dtheta * dz
        endif

        volume_avg = volume_avg + volume(i,j,k)
        n_avg = n_avg + 1
      enddo
    enddo
  enddo

  ! TODO: Gather volume_avg for parallel version.
  volume_avg = volume_avg / dble(n_avg)

end subroutine calculate_cell_volumes

!===========================================================================================
! Shift grid coordinates to the nodes where the permittivity is defined.
! This is needed to for example visualize the sensitivity data which is calculated on those nodes.
!===========================================================================================
subroutine compute_shifted_grid_coordinates(nr,ntheta,nz,r,theta,z,xgrid,ygrid,zgrid)
  integer, intent(in) :: nr, ntheta, nz
  real(kind=CUSTOM_REAL), intent(in) :: r(0:)
  real(kind=CUSTOM_REAL), intent(in) :: theta(0:)
  real(kind=CUSTOM_REAL), intent(in) :: z(0:)

  real(kind=CUSTOM_REAL), intent(out) :: xgrid(0:,0:,0:)
  real(kind=CUSTOM_REAL), intent(out) :: ygrid(0:,0:,0:)
  real(kind=CUSTOM_REAL), intent(out) :: zgrid(0:,0:,0:)

  real(kind=CUSTOM_REAL), allocatable :: r_shifted(:)
  real(kind=CUSTOM_REAL), allocatable :: theta_shifted(:)
  real(kind=CUSTOM_REAL), allocatable :: z_shifted(:)
  integer :: i,j,k

  allocate(r_shifted(0:2*nr+1))
  allocate(theta_shifted(0:ntheta+1))
  allocate(z_shifted(0:nz+1))

  ! Shift the radius r.
  ! Permittivity permit(i,:,:) is located between phi(i-1,:,:) and phi(i,:,:).
  do i=1,2*nr+1
    r_shifted(i) = r(i-1)
  enddo
  ! Since we have a mapping i-->2*i-1 and therefore r(-1)=0.
  r_shifted(0) = 0._CUSTOM_REAL

  ! Shift the angle theta.
  ! Permittivity permit(:,j,:) is located between phi(:,j,:) and phi(:,j+1,:).
  do j=1,ntheta
    theta_shifted(j) = theta(j) + 0.5_CUSTOM_REAL*(theta(j+1)-theta(j))
  enddo

  ! Apply periodic BCs.
  theta_shifted(0) = theta_shifted(ntheta) + 2._CUSTOM_REAL*PI
  theta_shifted(ntheta+1) = theta_shifted(1) - 2._CUSTOM_REAL*PI

  ! Shift the height z.
  ! Permittivity permit(:,:,k) is located between phi(:,:,k-1) and phi(:,:,k).
  do k=1,nz+1
    z_shifted(k) = z(k) - 0.5_CUSTOM_REAL*(z(k)-z(k-1))
  enddo
  ! We do not have permittivity defined between z(0)=0 and z(-1).
  ! So just define z_shifted(0) to do not have undefined elements.
  z_shifted(0) = 0._CUSTOM_REAL

  call compute_grid_coordinates(nr,ntheta,nz,r_shifted,theta_shifted,z_shifted,xgrid,ygrid,zgrid)

  deallocate(z_shifted)
  deallocate(theta_shifted)
  deallocate(r_shifted)

end subroutine compute_shifted_grid_coordinates

!===========================================================================================
! Increase twice the theta-dimension for the electrode ielectrode and its gaps.
subroutine refine_geometry(ielectrode,ntheta,theta,nel,elecs,&
                           ntheta0,theta0,elecs0,ithetamap)

  integer, intent(in) :: ielectrode,ntheta,nel,ntheta0
  real(kind=CUSTOM_REAL), intent(inout) :: theta0(0:), theta(0:)
  type(t_electrode), intent(inout) :: elecs(:)
  type(t_electrode), intent(inout) :: elecs0(:)
  integer, intent(inout) :: ithetamap(0:)

  integer, save :: originalValuesStored = 0
  integer, save :: npelec,npgaps

  integer :: j,l,el
  logical :: toInsert

  if (originalValuesStored == 0) then
    theta0(:) = theta(:)
    elecs0(:) = elecs(:)

    npelec = get_npelec(ntheta0, nel)
    npgaps = get_npgaps(ntheta0, nel)

    originalValuesStored = 1
  endif

  l=0
  do j=1,ntheta0+1
    toInsert = .false.

    if (ielectrode > 1 .and. ielectrode < 12) then
      if (j > elecs0(ielectrode-1)%ithetamax .and. j <= elecs0(ielectrode+1)%ithetamin) toInsert = .true.
    else if (ielectrode == 1) then
      if ((j > elecs0(ielectrode)%ithetamin .and. j <= elecs0(ielectrode+1)%ithetamin) &
          .or. (j > elecs0(nel)%ithetamax)) toInsert = .true.
    else if (ielectrode == 12) then
      if (j > elecs0(ielectrode-1)%ithetamax) toInsert = .true.
    endif

    if (toInsert) then
    ! Insert an extra grid dimension between theta(j-1) and theta(j).
      l=l+1
      theta(l) = 0.5_CUSTOM_REAL*(theta0(j-1)+theta0(j))
    endif

    l=l+1
    theta(l) = theta0(j)

    ithetamap(j) = l

    ! Sanity check.
    if (j == ntheta0+1) then
      if (l /= ntheta+1) then
        print *, "Error in mesh refinement! j=", j, ", l=", l
        stop
      endif
    endif
  enddo

  do el=1,nel
    if (el < ielectrode) then
        elecs(el)%ithetamin = elecs0(el)%ithetamin
        elecs(el)%ithetamax = elecs0(el)%ithetamax
    else if (el == ielectrode) then
      if (ielectrode == 1) then
        elecs(el)%ithetamax = elecs0(el)%ithetamax + (npelec-1)
      else  if (ielectrode > 1) then
        elecs(el)%ithetamin = elecs0(el)%ithetamin + npgaps+1
        elecs(el)%ithetamax = elecs0(el)%ithetamax + npgaps + npelec
      endif
    else if (el > ielectrode) then
      if (ielectrode == 1) then
        elecs(el)%ithetamin = elecs0(el)%ithetamin + npgaps + npelec
        elecs(el)%ithetamax = elecs0(el)%ithetamax + npgaps + npelec
      else
        elecs(el)%ithetamin = elecs0(el)%ithetamin + 2*(npgaps+1) + (npelec-1)
        elecs(el)%ithetamax = elecs0(el)%ithetamax + 2*(npgaps+1) + (npelec-1)
      endif
    endif
  enddo

  ! Periodic boundaries.
  theta(0) = theta(ntheta) + 2.0_CUSTOM_REAL*PI
  ithetamap(0) = ithetamap(ntheta0)

end subroutine refine_geometry

end module geometry
