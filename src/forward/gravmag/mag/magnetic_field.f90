
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

!========================================================================
! Custom Tomofast-x module for computation of total magnetic field.
! This module computes the magnetic kernel which can be used
! to calculate the total field at any point outside the magnetizing volume.
!========================================================================
module magnetic_field
  use global_typedefs
  use grid

  implicit none
  private

  ! Degree to radian.
  double precision, parameter :: d2rad = PI / 180.d0

  ! The vacuum magnetic permeability.
  double precision, parameter :: mu0 = 4.d0 * PI * 1.d-7

  ! Conversion from tesla to nanotesla.
  double precision, parameter :: T2nT = 1.d+9

  ! Main class.
  type, public :: t_magnetic_field
    private

    ! Ebxternal field intensity in nT.
    double precision :: intensity

    ! Magnetic field vectors.
    double precision, dimension(3) :: magv

    contains
      private

      procedure, public, pass     :: initialize => magnetic_field_initialize
      procedure, public, pass     :: magprism => magnetic_field_magprism

      procedure, private, nopass  :: sharmbox
      procedure, private, nopass  :: dircos

    end type t_magnetic_field

contains

!==============================================================================
! Initializes the magnetic field by calculating and storing their direction
! cosines.
! Currently the ambient magnetic field components are not used
!==============================================================================
subroutine magnetic_field_initialize(this, mi, md, theta, intensity)
  class(t_magnetic_field), intent(inout) :: this
  double precision, intent(in) :: mi, md, theta, intensity

  double precision :: ma, mb, mc

  this%intensity = intensity
  call this%dircos(mi, md, theta, ma, mb, mc)

  this%magv(1) = ma
  this%magv(2) = mb
  this%magv(3) = mc

end subroutine magnetic_field_initialize

!==============================================================================
!  Subroutine DIRCOS computes direction cosines from inclination
!  and declination.
!
!  Input parameters:
!    incl:  inclination in degrees positive below horizontal.
!    decl:  declination in degrees positive east of true north.
!    azim:  azimuth of x axis in degrees positive east of north.
!
!  Output parameters:
!    a, b, c:  the three direction cosines.
!==============================================================================
subroutine dircos(incl, decl, azim, a, b, c)
  double precision, intent(in)    :: incl, decl, azim
  double precision, intent(out)   :: a, b, c

  double precision :: xincl, xdecl, xazim
  double precision :: incl2, decl2

  ! Convert North to cartesian X-axis.
  decl2 = mod(450.d0 - decl,  360.d0)
  incl2 = incl

  xincl = incl2 * d2rad
  xdecl = decl2 * d2rad
  xazim = azim * d2rad

  a = cos(xincl) * cos(xdecl - xazim)
  b = cos(xincl) * sin(xdecl - xazim)
  c = sin(xincl)

end subroutine dircos

!===================================================================================================
! Calculates the magnetic sensitivity kernel.
! The model is either susceptibility (scalar) or magnetisation vector.
! The data is either TMI or three-component data.
! Uses X for East, Y for North, and Z points downwards.
!===================================================================================================
subroutine magnetic_field_magprism(this, nelements, nmodel_components, ndata_components, grid, Xdata, Ydata, Zdata, sensit_line)
  class(t_magnetic_field), intent(in)     :: this
  integer, intent(in)                     :: nelements, nmodel_components, ndata_components
  type(t_grid), intent(in)                :: grid
  real(kind=CUSTOM_REAL), intent(in)      :: Xdata, Ydata, Zdata

  real(kind=CUSTOM_REAL), intent(out)     :: sensit_line(nelements, nmodel_components, ndata_components)

  integer :: i, k
  real(kind=SENSIT_REAL) :: tx(3), ty(3), tz(3)
  double precision :: mx, my, mz

  integer :: j
  real(kind=SENSIT_REAL) :: temp_tx(3), temp_ty(3), temp_tz(3)
  real(kind=SENSIT_REAL) :: temp_x1(6), temp_x2(6), temp_y1(6), temp_y2(6), temp_z1(6), temp_z2(6)
  real(kind=SENSIT_REAL) :: width, min_clr

  do i = 1, nelements
    ! Calculate the magnetic tensor.

    ! Check if the point is inside the model grid.
    if ((grid%X1(i) < Xdata) .and. (grid%X2(i) > Xdata) .and. &
        (grid%Y1(i) < Ydata) .and. (grid%Y2(i) > Ydata) .and. &
        (grid%Z1(i) < Zdata) .and. (grid%Z2(i) > Zdata)) then

        ! Default void width.
        width = 0.1

        ! Drillhole observation point is not guaranteed to be at the enter of the voxel so its face clearance needs to be checked.
        ! Check if width actually exceeds the observation point's cleareance to each face of the voxel.
        ! If so, set void width to 50% of the minimum clearance, otherwise the default width is set.
        min_clr = min(abs(Xdata - grid%X1(i)), abs(Xdata - grid%X2(i)), &
                      abs(Ydata - grid%Y1(i)), abs(Ydata - grid%Y2(i)), &
                      abs(Zdata - grid%Z1(i)), abs(Zdata - grid%Z2(i)))

        if (width > min_clr) width = 0.5 * min_clr

        ! Calculate the 6 subvoxels coordinates.

        ! Top.
        temp_x1(1) = grid%X1(i)
        temp_x2(1) = grid%X2(i)
        temp_y1(1) = grid%Y1(i)
        temp_y2(1) = grid%Y2(i)
        temp_z1(1) = grid%Z1(i)
        temp_z2(1) = Zdata - width

        ! Bottom.
        temp_x1(2) = grid%X1(i)
        temp_x2(2) = grid%X2(i)
        temp_y1(2) = grid%Y1(i)
        temp_y2(2) = grid%Y2(i)
        temp_z1(2) = Zdata + width
        temp_z2(2) = grid%Z2(i)

        ! West.
        temp_x1(3) = grid%X1(i)
        temp_x2(3) = Xdata - width
        temp_y1(3) = grid%Y1(i)
        temp_y2(3) = grid%Y2(i)
        temp_z1(3) = Zdata - width
        temp_z2(3) = Zdata + width

        ! East.
        temp_x1(4) = Xdata + width
        temp_x2(4) = grid%X2(i)
        temp_y1(4) = grid%Y1(i)
        temp_y2(4) = grid%Y2(i)
        temp_z1(4) = Zdata - width
        temp_z2(4) = Zdata + width

        ! South.
        temp_x1(5) = Xdata - width
        temp_x2(5) = Xdata + width
        temp_y1(5) = grid%Y1(i)
        temp_y2(5) = Ydata - width
        temp_z1(5) = Zdata - width
        temp_z2(5) = Zdata + width

        ! North.
        temp_x1(6) = Xdata - width
        temp_x2(6) = Xdata + width
        temp_y1(6) = Ydata + width
        temp_y2(6) = grid%Y2(i)
        temp_z1(6) = Zdata - width
        temp_z2(6) = Zdata + width

        tx = 0.d0
        ty = 0.d0
        tz = 0.d0

        do j = 1, 6
            call this%sharmbox(real(Xdata, SENSIT_REAL), &
                               real(Ydata, SENSIT_REAL), &
                               real(Zdata, SENSIT_REAL), &
                               real(temp_x1(j), SENSIT_REAL), &
                               real(temp_y1(j), SENSIT_REAL), &
                               real(temp_z1(j), SENSIT_REAL), &
                               real(temp_x2(j), SENSIT_REAL), &
                               real(temp_y2(j), SENSIT_REAL), &
                               real(temp_z2(j), SENSIT_REAL), &
                               temp_tx, temp_ty, temp_tz)

            tx = tx + temp_tx
            ty = ty + temp_ty
            tz = tz + temp_tz
        enddo
    else
    ! Point is outside the model grid.

        call this%sharmbox(real(Xdata, SENSIT_REAL), &
                           real(Ydata, SENSIT_REAL), &
                           real(Zdata, SENSIT_REAL), &
                           real(grid%X1(i), SENSIT_REAL), &
                           real(grid%Y1(i), SENSIT_REAL), &
                           real(grid%Z1(i), SENSIT_REAL), &
                           real(grid%X2(i), SENSIT_REAL), &
                           real(grid%Y2(i), SENSIT_REAL), &
                           real(grid%Z2(i), SENSIT_REAL), &
                           tx, ty, tz)
    endif

    if (nmodel_components == 1) then
    ! Susceptibility model.

      mx = sum(tx * this%magv)
      my = sum(ty * this%magv)
      mz = sum(tz * this%magv)

      if (ndata_components == 1) then
        sensit_line(i, 1, 1) = mx * this%magv(1) + my * this%magv(2) + mz * this%magv(3)

      else if (ndata_components == 3) then
        sensit_line(i, 1, 1) = mx
        sensit_line(i, 1, 2) = my
        sensit_line(i, 1, 3) = mz

      else
        print *, "Wrong number of data components in magnetic_field_magprism!"
        stop
      endif

    else if (nmodel_components == 3) then
    ! Magnetisation model (Mx, My, Mz).

      if (ndata_components == 1) then
        do k = 1, 3
          sensit_line(i, k, 1) = tx(k) * this%magv(1) + ty(k) * this%magv(2) + tz(k) * this%magv(3)
        enddo

      else if (ndata_components == 3) then
        do k = 1, 3
          sensit_line(i, k, 1) = tx(k)
          sensit_line(i, k, 2) = ty(k)
          sensit_line(i, k, 3) = tz(k)
        enddo

      else
        print *, "Wrong number of data components in magnetic_field_magprism!"
        stop
      endif

    else
      print *, "Wrong number of model components in magnetic_field_magprism!"
      stop
    endif
  enddo

  if (nmodel_components == 1) then
    sensit_line = this%intensity * sensit_line

  else if (nmodel_components == 3) then
    ! The input magnetisation vector is in A/m and the output data is in nanotesla.
    sensit_line = (mu0 * T2nT) * sensit_line
  endif

  ! Convert to SI.
  sensit_line = sensit_line / (4.d0 * PI)

end subroutine magnetic_field_magprism

!===================================================================================
! Calculates the magnetic tensor.
! Uses the algorithm proposed by P. Vallabh Sharma in his 1966 paper:
! [Rapid Computation of Magnetic Anomalies and Demagnetization Effects Caused by Arbitrary Shape]
!
! Units:
!   coordinates:        m
!   field intensity:    nT
!   incl/decl/azi:      deg
!   mag suscept.:       cgs
!
! Inputs:
!   x0, y0, z0      coordinates of the observation point
!   x1, y1, z1      coordinates of one of the corners on the top face, where z1 is the depth
!   x2, y2, z2      coordinates of the opposite corner on the bottom face, where z2 is the depth
!
! Outputs:
!   tx = [txx txy txz]
!   ty = [tyx tyy tyz]
!   tz = [tzx tzy tzz]
!   components of the magnetic tensor
!===================================================================================
subroutine sharmbox(x0, y0, z0, x1, y1, z1, x2, y2, z2, ts_x, ts_y, ts_z)
  real(kind=SENSIT_REAL), intent(in) :: x0, y0, z0, x1, y1, z1, x2, y2, z2

  real(kind=SENSIT_REAL), intent(out) :: ts_x(3), ts_y(3), ts_z(3)

  real(kind=SENSIT_REAL) :: rx1, rx2, ry1, ry2, rz1, rz2
  real(kind=SENSIT_REAL) :: rx1sq, rx2sq, ry1sq, ry2sq, rz1sq, rz2sq
  real(kind=SENSIT_REAL) :: arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8
  real(kind=SENSIT_REAL) :: R1, R2, R3, R4
  real(kind=SENSIT_REAL) :: eps

  eps = 0.

  ! Relative coordinates to obs.
  ! Voxel runs from x1 to x2, y1 to y2, z1 to z2.
  ! Observation points at obs_x obs_y obs_z.
  rx1 = x1 - x0 + eps; ! rx1 = -u2
  rx2 = x2 - x0 + eps; ! rx2 = -u1
  ry1 = y1 - y0 + eps; ! ry1 = -v2
  ry2 = y2 - y0 + eps; ! ry2 = -v1
  rz1 = z1 - z0 + eps; ! rz1 = -w2
  rz2 = z2 - z0 + eps; ! rz2 = -w1

  ! Verify that the data does not lie on the boundary.
  ! We don't check for dZ since having dZ = 0 is fine as long as dX and dY are not zero.
  if (rx1 == 0. .or. rx2 == 0.) then
    print *, "The model grid X-boundary coincides with the data position, for data at x0, y0, z0 = ", x0, y0, z0
    stop
  endif

  if (ry1 == 0. .or. ry2 == 0.) then
    print *, "The model grid Y-boundary coincides with the data position, for data at x0, y0, z0 = ", x0, y0, z0
    stop
  endif

  ! Squares.
  rx1sq = rx1 ** 2; rx2sq = rx2 ** 2;
  ry1sq = ry1 ** 2; ry2sq = ry2 ** 2;
  rz1sq = rz1 ** 2; rz2sq = rz2 ** 2;

  R1 = ry2sq + rx2sq ! -v1**2 + -u1**2 -> R1
  R2 = ry2sq + rx1sq ! -v1**2 + -u2**2 -> R3
  R3 = ry1sq + rx2sq ! -v2**2 + -u1**2 -> R2
  R4 = ry1sq + rx1sq ! -v2**2 + -u2**2 -> R4

  arg1 = sqrt(rz2sq + R2)
  arg2 = sqrt(rz2sq + R1)
  arg3 = sqrt(rz1sq + R1)
  arg4 = sqrt(rz1sq + R2)
  arg5 = sqrt(rz2sq + R3)
  arg6 = sqrt(rz2sq + R4)
  arg7 = sqrt(rz1sq + R4)
  arg8 = sqrt(rz1sq + R3)

  ! ts_xx
  ts_x(1) = atan2(ry1 * rz2, (rx2 * arg5 + eps)) - &
            atan2(ry2 * rz2, (rx2 * arg2 + eps)) + &
            atan2(ry2 * rz1, (rx2 * arg3 + eps)) - &
            atan2(ry1 * rz1, (rx2 * arg8 + eps)) + &
            atan2(ry2 * rz2, (rx1 * arg1 + eps)) - &
            atan2(ry1 * rz2, (rx1 * arg6 + eps)) + &
            atan2(ry1 * rz1, (rx1 * arg7 + eps)) - &
            atan2(ry2 * rz1, (rx1 * arg4 + eps))

  ! ts_yx
  ts_y(1) = log((rz2 + arg2 + eps) / (rz1 + arg3 + eps)) - &
            log((rz2 + arg1 + eps) / (rz1 + arg4 + eps)) + &
            log((rz2 + arg6 + eps) / (rz1 + arg7 + eps)) - &
            log((rz2 + arg5 + eps) / (rz1 + arg8 + eps))

  ! ts_yy
  ts_y(2) = atan2(rx1 * rz2, (ry2 * arg1 + eps)) - &
            atan2(rx2 * rz2, (ry2 * arg2 + eps)) + &
            atan2(rx2 * rz1, (ry2 * arg3 + eps)) - &
            atan2(rx1 * rz1, (ry2 * arg4 + eps)) + &
            atan2(rx2 * rz2, (ry1 * arg5 + eps)) - &
            atan2(rx1 * rz2, (ry1 * arg6 + eps)) + &
            atan2(rx1 * rz1, (ry1 * arg7 + eps)) - &
            atan2(rx2 * rz1, (ry1 * arg8 + eps))

  ! Following computations do not reuse variables so it may be
  ! faster to just compute them directly instead of storing them.
  ! It does help legibility, however.
  R1 = ry2sq + rz1sq
  R2 = ry2sq + rz2sq
  R3 = ry1sq + rz1sq
  R4 = ry1sq + rz2sq

  arg1 = sqrt(rx1sq + R1)
  arg2 = sqrt(rx2sq + R1)
  arg3 = sqrt(rx1sq + R2)
  arg4 = sqrt(rx2sq + R2)
  arg5 = sqrt(rx1sq + R3)
  arg6 = sqrt(rx2sq + R3)
  arg7 = sqrt(rx1sq + R4)
  arg8 = sqrt(rx2sq + R4)

  ! ts_yz
  ts_y(3) = log((rx1 + arg1 + eps) / (rx2 + arg2 + eps)) - &
            log((rx1 + arg3 + eps) / (rx2 + arg4 + eps)) + &
            log((rx1 + arg7 + eps) / (rx2 + arg8 + eps)) - &
            log((rx1 + arg5 + eps) / (rx2 + arg6 + eps))

  R1 = rx2sq + rz1sq
  R2 = rx2sq + rz2sq
  R3 = rx1sq + rz1sq
  R4 = rx1sq + rz2sq

  arg1 = sqrt(ry1sq + R1)
  arg2 = sqrt(ry2sq + R1)
  arg3 = sqrt(ry1sq + R2)
  arg4 = sqrt(ry2sq + R2)
  arg5 = sqrt(ry1sq + R3)
  arg6 = sqrt(ry2sq + R3)
  arg7 = sqrt(ry1sq + R4)
  arg8 = sqrt(ry2sq + R4)

  ! ts_xz
  ts_x(3) = log((ry1 + arg1 + eps) / (ry2 + arg2 + eps)) - &
            log((ry1 + arg3 + eps) / (ry2 + arg4 + eps)) + &
            log((ry1 + arg7 + eps) / (ry2 + arg8 + eps)) - &
            log((ry1 + arg5 + eps) / (ry2 + arg6 + eps))

  ! Filling the rest of the tensor.
  ! ts_zz
  ts_z(3) = -1 * (ts_x(1) + ts_y(2)) ! Gauss

  ! ts_zy
  ts_z(2) = ts_y(3)

  ! ts_xy
  ts_x(2) = ts_y(1)

  ! ts_zx
  ts_z(1) = ts_x(3)

end subroutine sharmbox

end module magnetic_field
