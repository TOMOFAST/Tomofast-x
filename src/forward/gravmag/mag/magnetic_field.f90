
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
! This module only computes the magnetic tensor which
! can be used to calculate the total field at any point
! outside the magnetizing volume.
!========================================================================
module magnetic_field
  use global_typedefs
  use grid

  implicit none
  private

  ! Degree to radian.
  double precision, parameter :: d2rad = PI / 180.d0

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

      procedure, private, nopass  :: guobox
      procedure, private, nopass  :: mbox2

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

    ! Convert North to cartesian-X.
    decl2 = mod(450.d0 - decl,  360.d0)
    incl2 = -incl

    xincl = incl2 * d2rad
    xdecl = decl2 * d2rad
    xazim = azim * d2rad

    a = dcos(xincl) * dcos(xdecl - xazim)
    b = dcos(xincl) * dsin(xdecl - xazim)
    c = dsin(xincl)

end subroutine dircos

!===================================================================================================
! Calculates the magnetic tensor for each point and is returned to the calling process
! This subroutine is meant to perform mbox on a set of voxels before returning their
! respective magnetic tensor flattened in vector form.
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

    do i = 1, nelements
        ! Calculate the magnetic tensor.
        call this%mbox2(real(Xdata, SENSIT_REAL), &
                           real(Ydata, SENSIT_REAL), &
                           real(Zdata, SENSIT_REAL), &
                           real(grid%X1(i), SENSIT_REAL), &
                           real(grid%Y1(i), SENSIT_REAL), &
                           real(grid%Z1(i), SENSIT_REAL), &
                           real(grid%X2(i), SENSIT_REAL), &
                           real(grid%Y2(i), SENSIT_REAL), &
                           real(grid%Z2(i), SENSIT_REAL), &
                           tx, ty, tz)

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
    endif

    ! Convert to SI.
    sensit_line = sensit_line / (4.d0 * PI)

end subroutine magnetic_field_magprism

!===================================================================================
! Rewritten MBOX to use the algorithm proposed by P. Vallabh Sharma in his 1966 paper
! [Rapid Computation of Magnetic Anomalies and Demagnetization Effects Caused by Arbitrary Shape]
! Requires the DIRCOS subroutine.
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
    real(kind=SENSIT_REAL) :: four_pi
    logical :: l_inside

    eps = 1.e-8
    four_pi = 4 * 3.1415926535897932385_SENSIT_REAL

    l_inside = .false.

    ! Relative coordinates to obs.
    ! Voxel runs from x1 to x2, y1 to y2, z1 to z2.
    ! Observation points at obs_x obs_y obs_z.
    rx1 = x1 - x0 + eps; ! rx1 = -u2
    rx2 = x2 - x0 + eps; ! rx2 = -u1
    ry1 = y1 - y0 + eps; ! ry1 = -v2
    ry2 = y2 - y0 + eps; ! ry2 = -v1
    rz1 = z1 - z0 + eps; ! rz1 = -w2
    rz2 = z2 - z0 + eps; ! rz2 = -w1

    ! Squares.
    rx1sq = rx1 ** 2; rx2sq = rx2 ** 2;
    ry1sq = ry1 ** 2; ry2sq = ry2 ** 2;
    rz1sq = rz1 ** 2; rz2sq = rz2 ** 2;

    R1 = ry2sq + rx2sq ! -v1**2 + -u1**2 -> R1
    R2 = ry2sq + rx1sq ! -v1**2 + -u2**2 -> R3
    R3 = ry1sq + rx2sq ! -v2**2 + -u1**2 -> R2
    R4 = ry1sq + rx1sq ! -v2**2 + -u2**2 -> R4
    arg1 = SQRT(rz2sq + R2)
    arg2 = SQRT(rz2sq + R1)
    arg3 = SQRT(rz1sq + R1)
    arg4 = SQRT(rz1sq + R2)
    arg5 = SQRT(rz2sq + R3)
    arg6 = SQRT(rz2sq + R4)
    arg7 = SQRT(rz1sq + R4)
    arg8 = SQRT(rz1sq + R3)

    ! ts_xx
    ts_x(1) =   ATAN2(ry1 * rz2, (rx2 * arg5 + eps)) - &
                ATAN2(ry2 * rz2, (rx2 * arg2 + eps)) + &
                ATAN2(ry2 * rz1, (rx2 * arg3 + eps)) - &
                ATAN2(ry1 * rz1, (rx2 * arg8 + eps)) + &
                ATAN2(ry2 * rz2, (rx1 * arg1 + eps)) - &
                ATAN2(ry1 * rz2, (rx1 * arg6 + eps)) + &
                ATAN2(ry1 * rz1, (rx1 * arg7 + eps)) - &
                ATAN2(ry2 * rz1, (rx1 * arg4 + eps))

    ! ts_yx
    ts_y(1) =   LOG((rz2 + arg2 + eps) / (rz1 + arg3 + eps)) - &
                LOG((rz2 + arg1 + eps) / (rz1 + arg4 + eps)) + &
                LOG((rz2 + arg6 + eps) / (rz1 + arg7 + eps)) - &
                LOG((rz2 + arg5 + eps) / (rz1 + arg8 + eps))

    ! ts_yy
    ts_y(2) =   ATAN2(rx1 * rz2, (ry2 * arg1 + eps)) - &
                ATAN2(rx2 * rz2, (ry2 * arg2 + eps)) + &
                ATAN2(rx2 * rz1, (ry2 * arg3 + eps)) - &
                ATAN2(rx1 * rz1, (ry2 * arg4 + eps)) + &
                ATAN2(rx2 * rz2, (ry1 * arg5 + eps)) - &
                ATAN2(rx1 * rz2, (ry1 * arg6 + eps)) + &
                ATAN2(rx1 * rz1, (ry1 * arg7 + eps)) - &
                ATAN2(rx2 * rz1, (ry1 * arg8 + eps))

    ! Following computations do not reuse variables so it may be
    ! faster to just compute them directly instead of storing them.
    ! It does help legibility, however.
    R1 = ry2sq + rz1sq
    R2 = ry2sq + rz2sq
    R3 = ry1sq + rz1sq
    R4 = ry1sq + rz2sq
    arg1 = SQRT(rx1sq + R1)
    arg2 = SQRT(rx2sq + R1)
    arg3 = SQRT(rx1sq + R2)
    arg4 = SQRT(rx2sq + R2)
    arg5 = SQRT(rx1sq + R3)
    arg6 = SQRT(rx2sq + R3)
    arg7 = SQRT(rx1sq + R4)
    arg8 = SQRT(rx2sq + R4)

    ! ts_yz
    ts_y(3) =   LOG((rx1 + arg1 + eps) / (rx2 + arg2 + eps)) - &
                LOG((rx1 + arg3 + eps) / (rx2 + arg4 + eps)) + &
                LOG((rx1 + arg7 + eps) / (rx2 + arg8 + eps)) - &
                LOG((rx1 + arg5 + eps) / (rx2 + arg6 + eps))

    R1 = rx2sq + rz1sq
    R2 = rx2sq + rz2sq
    R3 = rx1sq + rz1sq
    R4 = rx1sq + rz2sq
    arg1 = SQRT(ry1sq + R1)
    arg2 = SQRT(ry2sq + R1)
    arg3 = SQRT(ry1sq + R2)
    arg4 = SQRT(ry2sq + R2)
    arg5 = SQRT(ry1sq + R3)
    arg6 = SQRT(ry2sq + R3)
    arg7 = SQRT(ry1sq + R4)
    arg8 = SQRT(ry2sq + R4)

    ! ts_xz
    ts_x(3) =   LOG((ry1 + arg1 + eps) / (ry2 + arg2 + eps)) - &
                LOG((ry1 + arg3 + eps) / (ry2 + arg4 + eps)) + &
                LOG((ry1 + arg7 + eps) / (ry2 + arg8 + eps)) - &
                LOG((ry1 + arg5 + eps) / (ry2 + arg6 + eps))

    ! Checking if point is inside the voxel.
    ! If so, use poisson's relation.
    if (x0 >= x1 .and. x0 <= x2) then
        if (y0 >= y1 .and. y0 <= y2) then
            if (z0 >= z1 .and. z0 <= z2) then
                l_inside = .true.
            endif
        endif
    endif

    ! Filling the rest of the tensor.
    ! ts_zz
    if (l_inside) then
        print *, "Observation point inside target voxel!"
        print *, "Obs: ", x0, y0, z0
        print *, "Voxel: ", x1, x2, y1, y2, z1, z2

        ts_z(3) = -1 * (ts_x(1) + ts_y(2) + four_pi) ! poisson
    else
        ts_z(3) = -1 * (ts_x(1) + ts_y(2)) ! gauss
    endif

    ! ts_zy
    ts_z(2) = ts_y(3)

    ! ts_xy
    ts_x(2) = ts_y(1)

    ! ts_zx
    ts_z(1) = ts_x(3)

end subroutine sharmbox

!================================================================================
subroutine guobox(x, y, z, x1, y1, z1, x2, y2, z2, ts_x, ts_y, ts_z)
  real(kind=SENSIT_REAL), intent(in) :: x, y, z, x1, y1, z1, x2, y2, z2

  real(kind=SENSIT_REAL), intent(out) :: ts_x(3), ts_y(3), ts_z(3)

  real(kind=SENSIT_REAL) :: sx(2), sy(2), sz(2)
  real(kind=SENSIT_REAL) :: R
  real(kind=SENSIT_REAL) :: ss(2), s
  real(kind=SENSIT_REAL) :: four_pi
  integer :: i, j, k

  ! Ranges of the source point in cuboid.
  sx(1) = x1
  sx(2) = x2

  sy(1) = y1
  sy(2) = y2

  sz(1) = z1
  sz(2) = z2

  ts_x = 0._SENSIT_REAL
  ts_y = 0._SENSIT_REAL
  ts_z = 0._SENSIT_REAL

  ss(1) = -1._SENSIT_REAL
  ss(2) = 1._SENSIT_REAL

  ! Calculate the integral.
  do i = 1, 2
    do j = 1, 2
      do k = 1, 2
        ! Sign of the integral component.
        s = ss(i) * ss(j) * ss(k)

        R = sqrt((x - sx(i))**2 + (y - sy(j))**2 + (z - sz(k))**2)

        ts_x(1) = ts_x(1) - s * atan2((sx(i) - x) * (sy(j) - y), (sx(i) - x)**2 + R * (sz(k) - z) + (sz(k) - z)**2)
        ts_x(2) = ts_x(2) + s * log(sz(k) - z + R)
        ts_x(3) = ts_x(3) + s * log(sy(j) - y + R)

        ts_y(2) = ts_y(2) - s * atan2((sx(i) - x) * (sy(j) - y), (sy(j) - y)**2 + R * (sz(k) - z) + (sz(k) - z)**2)
        ts_y(3) = ts_y(3) + s * log(sx(i) - x + R)

        ts_z(3) = ts_z(3) - s * atan2((sx(i) - x) * (sy(j) - y), R * (sz(k) - z))
      enddo
    enddo
  enddo

  ! Fixing the total trace to make it zero, and to align the expression with Guo et al (2004).
  ts_x(1) = - ts_x(1)
  ts_y(2) = - ts_y(2)

  ! Filling the symmetric components.
  ts_y(1) = ts_x(2)
  ts_z(1) = ts_x(3)
  ts_z(2) = ts_y(3)

end subroutine guobox

!================================================================================
subroutine mbox2(x, y, z, x1, y1, z1, x2, y2, z2, ts_x, ts_y, ts_z)
  real(kind=SENSIT_REAL), intent(in) :: x, y, z, x1, y1, z1, x2, y2, z2

  real(kind=SENSIT_REAL), intent(out) :: ts_x(3), ts_y(3), ts_z(3)

  real(kind=SENSIT_REAL) :: dx1, dx2, dy1, dy2, dz1, dz2
  real(kind=SENSIT_REAL) :: dx1dx1, dx2dx2, dy1dy1, dy2dy2, dz1dz1, dz2dz2
  real(kind=SENSIT_REAL) :: R1_, R2_, R3_, R4_
  real(kind=SENSIT_REAL) :: r1, r2, r3, r4, r5, r6, r7, r8
  real(kind=SENSIT_REAL) :: arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10
  real(kind=SENSIT_REAL) :: arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20
  real(kind=SENSIT_REAL) :: arg21, arg22, arg23, arg24, arg25, arg26, arg27, arg28, arg29, arg30
  real(kind=SENSIT_REAL) :: arg31, arg32, arg33, arg34, arg35, arg36, arg37, arg38, arg39, arg40
  real(kind=SENSIT_REAL) :: arg1_, arg6_, arg11_, arg16_, arg21_, arg26_, arg31_, arg36_
  real(kind=SENSIT_REAL) :: tol1
  real(kind=SENSIT_REAL) :: four_pi
  logical :: l_inside

  four_pi = 4 * 3.1415926535897932385_SENSIT_REAL

  ! Tolerance for numerical stability over nodes and edges.
  tol1 = 1.e-10

  ! Relative coordinates to obs.
  dx1 = x1 - x
  dx2 = x2 - x

  dy1 = y1 - y
  dy2 = y2 - y

  dz1 = z1 - z
  dz2 = z2 - z

  ! Squared difference.
  dx2dx2 = dx2**2
  dx1dx1 = dx1**2

  dy2dy2 = dy2**2
  dy1dy1 = dy1**2

  dz2dz2 = dz2**2
  dz1dz1 = dz1**2

  ! The 2D radius component squared of corner nodes.
  R1_ = dy2dy2 + dx2dx2
  R2_ = dy2dy2 + dx1dx1
  R3_ = dy1dy1 + dx2dx2
  R4_ = dy1dy1 + dx1dx1

  ! Radius to each cell node
  r1 = sqrt(dz2dz2 + R2_)
  r2 = sqrt(dz2dz2 + R1_)
  r3 = sqrt(dz1dz1 + R1_)
  r4 = sqrt(dz1dz1 + R2_)
  r5 = sqrt(dz2dz2 + R3_)
  r6 = sqrt(dz2dz2 + R4_)
  r7 = sqrt(dz1dz1 + R4_)
  r8 = sqrt(dz1dz1 + R3_)

  ! Argument calculations.
  arg1_ = dx1 + dy2 + r1
  arg1 = dy2 + dz2 + r1
  arg2 = dx1 + dz2 + r1
  arg3 = dx1 + r1
  arg4 = dy2 + r1
  arg5 = dz2 + r1

  arg6_ = dx2 + dy2 + r2
  arg6 = dy2 + dz2 + r2
  arg7 = dx2 + dz2 + r2
  arg8 = dx2 + r2
  arg9 = dy2 + r2
  arg10 = dz2 + r2

  arg11_ = dx2 + dy2 + r3
  arg11 = dy2 + dz1 + r3
  arg12 = dx2 + dz1 + r3
  arg13 = dx2 + r3
  arg14 = dy2 + r3
  arg15 = dz1 + r3

  arg16_ = dx1 + dy2 + r4
  arg16 = dy2 + dz1 + r4
  arg17 = dx1 + dz1 + r4
  arg18 = dx1 + r4
  arg19 = dy2 + r4
  arg20 = dz1 + r4

  arg21_ = dx2 + dy1 + r5
  arg21 = dy1 + dz2 + r5
  arg22 = dx2 + dz2 + r5
  arg23 = dx2 + r5
  arg24 = dy1 + r5
  arg25 = dz2 + r5

  arg26_ = dx1 + dy1 + r6
  arg26 = dy1 + dz2 + r6
  arg27 = dx1 + dz2 + r6
  arg28 = dx1 + r6
  arg29 = dy1 + r6
  arg30 = dz2 + r6

  arg31_ = dx1 + dy1 + r7
  arg31 = dy1 + dz1 + r7
  arg32 = dx1 + dz1 + r7
  arg33 = dx1 + r7
  arg34 = dy1 + r7
  arg35 = dz1 + r7

  arg36_ = dx2 + dy1 + r8
  arg36 = dy1 + dz1 + r8
  arg37 = dx2 + dz1 + r8
  arg38 = dx2 + r8
  arg39 = dy1 + r8
  arg40 = dz1 + r8

  ! Calculate tensor components.
  ts_x(1) = &
        (-2 * atan2(dx1, arg1 + tol1)) &
      - (-2 * atan2(dx2, arg6 + tol1)) &
      + (-2 * atan2(dx2, arg11 + tol1)) &
      - (-2 * atan2(dx1, arg16 + tol1)) &
      + (-2 * atan2(dx2, arg21 + tol1)) &
      - (-2 * atan2(dx1, arg26 + tol1)) &
      + (-2 * atan2(dx1, arg31 + tol1)) &
      - (-2 * atan2(dx2, arg36 + tol1))

  ts_x(2) = &
        log(arg5) - log(arg10) &
      + log(arg15) - log(arg20) &
      + log(arg25) - log(arg30) &
      + log(arg35) - log(arg40)

  ts_x(3) = &
        log(arg4) - log(arg9) &
      + log(arg14) - log(arg19) &
      + log(arg24) - log(arg29) &
      + log(arg34) - log(arg39)

  ts_y(2) = &
        (-2 * atan2(dy2, arg2 + tol1)) &
      - (-2 * atan2(dy2, arg7 + tol1)) &
      + (-2 * atan2(dy2, arg12 + tol1)) &
      - (-2 * atan2(dy2, arg17 + tol1)) &
      + (-2 * atan2(dy1, arg22 + tol1)) &
      - (-2 * atan2(dy1, arg27 + tol1)) &
      + (-2 * atan2(dy1, arg32 + tol1)) &
      - (-2 * atan2(dy1, arg37 + tol1))

  ts_y(3) = &
        log(arg3) - log(arg8) &
      + log(arg13) - log(arg18) &
      + log(arg23) - log(arg28) &
      + log(arg33) - log(arg38)

  ! Filling the symmetric components.
  ts_y(1) = ts_x(2)
  ts_z(1) = ts_x(3)
  ts_z(2) = ts_y(3)

  l_inside = .false.

  ! Checking if point is inside the voxel.
  ! If so, use Poisson's relation.
  if (x >= x1 .and. x <= x2) then
    if (y >= y1 .and. y <= y2) then
      if (z >= z1 .and. z <= z2) then
        l_inside = .true.
      endif
    endif
  endif

  ! Filling the rest of the diagonal.
  if (l_inside) then
    print *, "Observation point inside target voxel!"
    print *, "Obs: ", x, y, z
    print *, "Voxel: ", x1, x2, y1, y2, z1, z2

    ts_z(3) = -1 * (ts_x(1) + ts_y(2) + four_pi) ! poisson
  else
    ts_z(3) = -1 * (ts_x(1) + ts_y(2)) ! gauss
  endif

  ! Flip the sign.
  ts_x = - ts_x
  ts_y = - ts_y
  ts_z = - ts_z

end subroutine mbox2

end module magnetic_field
