
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
      procedure, private, nopass  :: sharmbox2
      procedure, private, nopass  :: dircos

      procedure, private, nopass  :: tensorbox
      procedure, private, nopass  :: calc_tensor_iii
      procedure, private, nopass  :: calc_tensor_iij

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

  real(kind=SENSIT_REAL) :: mtensor(10)!mtensor(3, 3, 3)
  real(kind=SENSIT_REAL) :: internal_mtensor(10)!internal_mtensor(3, 3, 3)

  integer :: j
  real(kind=SENSIT_REAL) :: temp_tx(3), temp_ty(3), temp_tz(3)
  real(kind=SENSIT_REAL) :: temp_x1(6), temp_x2(6), temp_y1(6), temp_y2(6), temp_z1(6), temp_z2(6)
  real(kind=SENSIT_REAL) :: width, min_clr
  logical :: inside_x, inside_y, inside_z

  !real(kind=SENSIT_REAL) :: znodes(grid%nx * grid%ny * 12), dummy_znodes(12), tensorZnodes(grid%nx * grid%ny * 40)
  real(kind=SENSIT_REAL), allocatable :: znodes(:), dummy_znodes(:)
  logical :: l_calcznodes
  integer :: nele_xylayer, xy_ind, len_znode

  ! allocate and initialize znode array depending on components (sharmbox or tensorbox)
  len_znode = MERGE(12, 40, ndata_components <= 3)
  nele_xylayer = grid%nx * grid%ny

  allocate(znodes(nele_xylayer * len_znode))
  allocate(dummy_znodes(nele_xylayer * len_znode))
  znodes = -9999.9
  dummy_znodes = -9999.9

  do i = 1, nelements
    ! Clear mtensor for each data observation point
    mtensor = 0.d0

    ! Find the cell's index on the X-Y plane
    xy_ind = MERGE(MOD(i, nele_xylayer), nele_xylayer, (MOD(i, nele_xylayer) > 0))

    ! Check if the znode values need to be calculated
    ! They need to be calculated if:
    ! 1. It's the first layer or
    ! 2. The znodes above the current voxel are invalid
    !l_calcznodes = .true.
    l_calcznodes = (i <= nele_xylayer) .or. (znodes((xy_ind - 1) * len_znode + 1) < -9999.0)

    ! Calculate the magnetic tensor.

    ! Check if the point is inside the model grid.
    inside_x = (grid%X1(i) < Xdata) .and. (grid%X2(i) > Xdata)
    inside_y = (grid%Y1(i) < Ydata) .and. (grid%Y2(i) > Ydata)
    inside_z = (grid%Z1(i) < Zdata) .and. (grid%Z2(i) > Zdata)

    if (inside_x .and. inside_y .and. inside_z) then

        ! update znode array with invalids
        znodes(xy_ind:xy_ind + len_znode) = -9999.9

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
        internal_mtensor = 0.d0

        do j = 1, 6
          if (ndata_components <= 3) then
            call this%sharmbox2(real(Xdata, SENSIT_REAL), &
                               real(Ydata, SENSIT_REAL), &
                               real(Zdata, SENSIT_REAL), &
                               real(temp_x1(j), SENSIT_REAL), &
                               real(temp_y1(j), SENSIT_REAL), &
                               real(temp_z1(j), SENSIT_REAL), &
                               real(temp_x2(j), SENSIT_REAL), &
                               real(temp_y2(j), SENSIT_REAL), &
                               real(temp_z2(j), SENSIT_REAL), &
                               temp_tx, temp_ty, temp_tz)!, dummy_znodes, .true., 1)!, 1)

            tx = tx + temp_tx
            ty = ty + temp_ty
            tz = tz + temp_tz

          else
            call this%tensorbox(real(Xdata, SENSIT_REAL), &
                               real(Ydata, SENSIT_REAL), &
                               real(Zdata, SENSIT_REAL), &
                               real(temp_x1(j), SENSIT_REAL), &
                               real(temp_y1(j), SENSIT_REAL), &
                               real(temp_z1(j), SENSIT_REAL), &
                               real(temp_x2(j), SENSIT_REAL), &
                               real(temp_y2(j), SENSIT_REAL), &
                               real(temp_z2(j), SENSIT_REAL), &
                               internal_mtensor, dummy_znodes, .true., 1)

            mtensor = mtensor + internal_mtensor

          endif
        enddo

    else
    ! Point is outside the voxel.
        if (ndata_components <= 3) then
          call this%sharmbox2(real(Xdata, SENSIT_REAL), &
                            real(Ydata, SENSIT_REAL), &
                            real(Zdata, SENSIT_REAL), &
                            real(grid%X1(i), SENSIT_REAL), &
                            real(grid%Y1(i), SENSIT_REAL), &
                            real(grid%Z1(i), SENSIT_REAL), &
                            real(grid%X2(i), SENSIT_REAL), &
                            real(grid%Y2(i), SENSIT_REAL), &
                            real(grid%Z2(i), SENSIT_REAL), &
                            tx, ty, tz)!, znodes, l_calcznodes, xy_ind)!, nele_xylayer)

        else
          call this%tensorbox(real(Xdata, SENSIT_REAL), &
                            real(Ydata, SENSIT_REAL), &
                            real(Zdata, SENSIT_REAL), &
                            real(grid%X1(i), SENSIT_REAL), &
                            real(grid%Y1(i), SENSIT_REAL), &
                            real(grid%Z1(i), SENSIT_REAL), &
                            real(grid%X2(i), SENSIT_REAL), &
                            real(grid%Y2(i), SENSIT_REAL), &
                            real(grid%Z2(i), SENSIT_REAL), &
                            mtensor, znodes, l_calcznodes, xy_ind)

        endif
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

      else if (ndata_components == 5) then
        ! Note: Mz kernel term multiplied by -1 to switch to elevation space as required by the equations
        ! bxx kernal parts - bxx dropped for 5c tests
        !sensit_line(i, 1, 1) = mtensor(1, 1, 1)
        !sensit_line(i, 2, 1) = mtensor(1, 1, 2)
        !sensit_line(i, 3, 1) = mtensor(1, 1, 3) * -1.0

        ! byy kernal parts
        sensit_line(i, 1, 1) = mtensor(4)!mtensor(1, 2, 2)
        sensit_line(i, 2, 1) = mtensor(7)!mtensor(2, 2, 2)
        sensit_line(i, 3, 1) = mtensor(8) * -1.0!mtensor(2, 2, 3) * -1.0

        ! bzz kernal parts
        sensit_line(i, 1, 2) = mtensor(6)!mtensor(1, 3, 3)
        sensit_line(i, 2, 2) = mtensor(9)!mtensor(2, 3, 3)
        sensit_line(i, 3, 2) = mtensor(10) * -1.0!mtensor(3, 3, 3) * -1.0

        ! bxy kernal parts
        sensit_line(i, 1, 3) = mtensor(2)!mtensor(1, 1, 2)
        sensit_line(i, 2, 3) = mtensor(4)!mtensor(1, 2, 2)
        sensit_line(i, 3, 3) = mtensor(5) * -1.0!mtensor(1, 2, 3) * -1.0

        ! byz kernal parts - flipped to correct inverted output
        sensit_line(i, 1, 4) = mtensor(5) * -1.0!mtensor(1, 2, 3) * -1.0
        sensit_line(i, 2, 4) = mtensor(8) * -1.0!mtensor(2, 2, 3) * -1.0
        sensit_line(i, 3, 4) = mtensor(9)!mtensor(2, 3, 3)

        ! bxz kernal parts - flipped to correct inverted output
        sensit_line(i, 1, 5) = mtensor(3) * -1.0!mtensor(1, 1, 3) * -1.0
        sensit_line(i, 2, 5) = mtensor(5) * -1.0!mtensor(1, 2, 3) * -1.0
        sensit_line(i, 3, 5) = mtensor(6)!mtensor(1, 3, 3)

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

  ! Cleanup
  deallocate(dummy_znodes)
  deallocate(znodes)

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
subroutine sharmbox(x0, y0, z0, x1, y1, z1, x2, y2, z2, ts_x, ts_y, ts_z, znodes, l_calcznodes, xy_ind)!, nele_xylayer)
  real(kind=SENSIT_REAL), intent(in) :: x0, y0, z0, x1, y1, z1, x2, y2, z2

  real(kind=SENSIT_REAL), intent(out) :: ts_x(3), ts_y(3), ts_z(3)

  real(kind=SENSIT_REAL) :: rx1, rx2, ry1, ry2, rz1, rz2
  real(kind=SENSIT_REAL) :: rx1sq, rx2sq, ry1sq, ry2sq, rz1sq, rz2sq
  real(kind=SENSIT_REAL) :: arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8
  real(kind=SENSIT_REAL) :: R1, R2, R3, R4
  real(kind=SENSIT_REAL) :: eps
  real(kind=SENSIT_REAL) :: at1, at2, at3, at4, at5, at6, at7, at8, lg1, lg2, lg3, lg4

  real(kind=SENSIT_REAL), intent(inout) :: znodes(:)
  logical, intent(in) :: l_calcznodes
  integer, intent(in) :: xy_ind!, nele_xylayer
  integer :: znodes_i!, prev_znodes_i

  !print *, x1, x2
  !print *, y1, y2
  !print *, z1, z2

  eps = 1.e-8
  znodes_i = (xy_ind - 1) * 12 + 1
  !prev_znodes_i = (i_ele - nele_xylayer - 1) * 12 + 1

  !print *, xy_ind, znodes_i

  ! Check for bad znode values, force recalc if so
  !if ((.not. l_firstlayer) .and. (znodes(znodes_i) < -9999.0)) then
  !  l_firstlayer = .true.
  !endif

  !l_firstlayer = ()

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

  ! Lets assume:
  ! C1 = SW corner of voxel top face (+Z axis)
  ! C2 = SE corner of voxel top face
  ! C3 = NW corner of voxel top face
  ! C4 = NE corner of voxel top face
  ! C5 = SW corner of voxel bottom face
  ! C6 = SE corner of voxel bottom face
  ! C7 = NW corner of voxel bottom face
  ! C8 = NE corner of voxel bottom face

  arg1 = sqrt(rz2sq + R2) ! C7
  arg2 = sqrt(rz2sq + R1) ! C8
  arg3 = sqrt(rz1sq + R1) ! C4
  arg4 = sqrt(rz1sq + R2) ! C3
  arg5 = sqrt(rz2sq + R3) ! C6
  arg6 = sqrt(rz2sq + R4) ! C5
  arg7 = sqrt(rz1sq + R4) ! C1
  arg8 = sqrt(rz1sq + R3) ! C2

  at1 = atan2(ry1 * rz2, (rx2 * arg5 + eps))
  at2 = atan2(ry2 * rz2, (rx2 * arg2 + eps))
  !at3 = atan2(ry2 * rz1, (rx2 * arg3 + eps))
  !at4 = atan2(ry1 * rz1, (rx2 * arg8 + eps))
  at5 = atan2(ry2 * rz2, (rx1 * arg1 + eps))
  at6 = atan2(ry1 * rz2, (rx1 * arg6 + eps))
  !at7 = atan2(ry1 * rz1, (rx1 * arg7 + eps))
  !at8 = atan2(ry2 * rz1, (rx1 * arg4 + eps))

  if (l_calcznodes) then
    !arg3 = sqrt(rz1sq + R1) ! C4 ! Commented out because these are still needed for ts_yx
    !arg4 = sqrt(rz1sq + R2) ! C3
    !arg7 = sqrt(rz1sq + R4) ! C1
    !arg8 = sqrt(rz1sq + R3) ! C2

    at3 = atan2(ry2 * rz1, (rx2 * arg3 + eps))
    at4 = atan2(ry1 * rz1, (rx2 * arg8 + eps))
    at7 = atan2(ry1 * rz1, (rx1 * arg7 + eps))
    at8 = atan2(ry2 * rz1, (rx1 * arg4 + eps))

  else
    at4 = znodes(znodes_i)
    at3 = znodes(znodes_i + 1)
    at8 = znodes(znodes_i + 2)
    at7 = znodes(znodes_i + 3)

  endif

  znodes(znodes_i) = at1
  znodes(znodes_i + 1) = at2
  znodes(znodes_i + 2) = at5
  znodes(znodes_i + 3) = at6


  !print *, 'xx at1', at1
  !print *, 'xx at2', at2
  !print *, 'xx at3', at3
  !print *, 'xx at4', at4
  !print *, 'xx at5', at5
  !print *, 'xx at6', at6
  !print *, 'xx at7', at7
  !print *, 'xx at8', at8

  ! ts_xx
  ! C6 - C8 + C4 - C2 + C7 - C5 + C1 - C3
  !ts_x(1) = atan2(ry1 * rz2, (rx2 * arg5 + eps)) - & ! atan2(-inf/inf)
  !          atan2(ry2 * rz2, (rx2 * arg2 + eps)) + & ! -atan2(-inf/inf)
  !          atan2(ry2 * rz1, (rx2 * arg3 + eps)) - &
  !          atan2(ry1 * rz1, (rx2 * arg8 + eps)) + &
  !          atan2(ry2 * rz2, (rx1 * arg1 + eps)) - & ! atan2(-inf/inf)
  !          atan2(ry1 * rz2, (rx1 * arg6 + eps)) + & ! -atan2(-inf/inf)
  !          atan2(ry1 * rz1, (rx1 * arg7 + eps)) - &
  !          atan2(ry2 * rz1, (rx1 * arg4 + eps))

  ts_x(1) = at1 - at2 + at3 - at4 + at5 - at6 + at7 - at8

  ! mapping
  ! at1 -> at4
  ! at2 -> at3
  ! at5 -> at8
  ! at6 -> at7

  lg1 = log((rz2 + arg2 + eps) / (rz1 + arg3 + eps))
  lg2 = log((rz2 + arg1 + eps) / (rz1 + arg4 + eps))
  lg3 = log((rz2 + arg6 + eps) / (rz1 + arg7 + eps))
  lg4 = log((rz2 + arg5 + eps) / (rz1 + arg8 + eps))

  !print *, 'yx lg1', lg1
  !print *, 'yx lg2', lg2
  !print *, 'yx lg3', lg3
  !print *, 'yx lg4', lg4

  ! ts_yx
  !ts_y(1) = log((rz2 + arg2 + eps) / (rz1 + arg3 + eps)) - &
  !          log((rz2 + arg1 + eps) / (rz1 + arg4 + eps)) + &
  !          log((rz2 + arg6 + eps) / (rz1 + arg7 + eps)) - &
  !          log((rz2 + arg5 + eps) / (rz1 + arg8 + eps))

  ts_y(1) = lg1 - lg2 + lg3 - lg4

  ! mapping
  ! no matches

  at1 = atan2(rx1 * rz2, (ry2 * arg1 + eps))
  at2 = atan2(rx2 * rz2, (ry2 * arg2 + eps))
  !at3 = atan2(rx2 * rz1, (ry2 * arg3 + eps))
  !at4 = atan2(rx1 * rz1, (ry2 * arg4 + eps))
  at5 = atan2(rx2 * rz2, (ry1 * arg5 + eps))
  at6 = atan2(rx1 * rz2, (ry1 * arg6 + eps))
  !at7 = atan2(rx1 * rz1, (ry1 * arg7 + eps))
  !at8 = atan2(rx2 * rz1, (ry1 * arg8 + eps))

  if (l_calcznodes) then
    at3 = atan2(rx2 * rz1, (ry2 * arg3 + eps))
    at4 = atan2(rx1 * rz1, (ry2 * arg4 + eps))
    at7 = atan2(rx1 * rz1, (ry1 * arg7 + eps))
    at8 = atan2(rx2 * rz1, (ry1 * arg8 + eps))

  else
    at4 = znodes(znodes_i + 4)
    at3 = znodes(znodes_i + 5)
    at8 = znodes(znodes_i + 6)
    at7 = znodes(znodes_i + 7)

  endif

  znodes(znodes_i + 4) = at1
  znodes(znodes_i + 5) = at2
  znodes(znodes_i + 6) = at5
  znodes(znodes_i + 7) = at6

  !print *, 'yy at1', at1
  !print *, 'yy at2', at2
  !print *, 'yy at3', at3
  !print *, 'yy at4', at4
  !print *, 'yy at5', at5
  !print *, 'yy at6', at6
  !print *, 'yy at7', at7
  !print *, 'yy at8', at8

  ! ts_yy
  !ts_y(2) = atan2(rx1 * rz2, (ry2 * arg1 + eps)) - &
  !          atan2(rx2 * rz2, (ry2 * arg2 + eps)) + &
  !          atan2(rx2 * rz1, (ry2 * arg3 + eps)) - &
  !          atan2(rx1 * rz1, (ry2 * arg4 + eps)) + &
  !          atan2(rx2 * rz2, (ry1 * arg5 + eps)) - &
  !          atan2(rx1 * rz2, (ry1 * arg6 + eps)) + &
  !          atan2(rx1 * rz1, (ry1 * arg7 + eps)) - &
  !          atan2(rx2 * rz1, (ry1 * arg8 + eps))

  ts_y(2) = at1 - at2 + at3 - at4 + at5 - at6 + at7 - at8

  ! mapping
  ! at1 -> at4
  ! at2 -> at3
  ! at5 -> at8
  ! at6 -> at7

  ! Following computations do not reuse variables so it may be
  ! faster to just compute them directly instead of storing them.
  ! It does help legibility, however.
  R1 = ry2sq + rz1sq
  R2 = ry2sq + rz2sq
  R3 = ry1sq + rz1sq
  R4 = ry1sq + rz2sq

  !arg1 = sqrt(rx1sq + R1)
  !arg2 = sqrt(rx2sq + R1)
  arg3 = sqrt(rx1sq + R2)
  arg4 = sqrt(rx2sq + R2)
  !arg5 = sqrt(rx1sq + R3)
  !arg6 = sqrt(rx2sq + R3)
  arg7 = sqrt(rx1sq + R4)
  arg8 = sqrt(rx2sq + R4)

  !lg1 = log((rx1 + arg1 + eps) / (rx2 + arg2 + eps))
  lg2 = log((rx1 + arg3 + eps) / (rx2 + arg4 + eps))
  lg3 = log((rx1 + arg7 + eps) / (rx2 + arg8 + eps))
  !lg4 = log((rx1 + arg5 + eps) / (rx2 + arg6 + eps))

  if (l_calcznodes) then
    arg1 = sqrt(rx1sq + R1)
    arg2 = sqrt(rx2sq + R1)
    arg5 = sqrt(rx1sq + R3)
    arg6 = sqrt(rx2sq + R3)

    lg1 = log((rx1 + arg1 + eps) / (rx2 + arg2 + eps))
    lg4 = log((rx1 + arg5 + eps) / (rx2 + arg6 + eps))

  else
    lg1 = znodes(znodes_i + 8)
    lg4 = znodes(znodes_i + 9)

  endif

  znodes(znodes_i + 8) = lg2
  znodes(znodes_i + 9) = lg3

  !print *, 'yz lg1', lg1
  !print *, 'yz lg2', lg2
  !print *, 'yz lg3', lg3
  !print *, 'yz lg4', lg4

  ! ts_yz
  !ts_y(3) = log((rx1 + arg1 + eps) / (rx2 + arg2 + eps)) - &
  !          log((rx1 + arg3 + eps) / (rx2 + arg4 + eps)) + &
  !          log((rx1 + arg7 + eps) / (rx2 + arg8 + eps)) - &
  !          log((rx1 + arg5 + eps) / (rx2 + arg6 + eps))

  ts_y(3) = lg1 - lg2 + lg3 - lg4

  ! mapping
  ! lg2 -> lg1
  ! lg3 -> lg4

  R1 = rx2sq + rz1sq
  R2 = rx2sq + rz2sq
  R3 = rx1sq + rz1sq
  R4 = rx1sq + rz2sq

  !arg1 = sqrt(ry1sq + R1)
  !arg2 = sqrt(ry2sq + R1)
  arg3 = sqrt(ry1sq + R2)
  arg4 = sqrt(ry2sq + R2)
  !arg5 = sqrt(ry1sq + R3)
  !arg6 = sqrt(ry2sq + R3)
  arg7 = sqrt(ry1sq + R4)
  arg8 = sqrt(ry2sq + R4)

  !lg1 = log((ry1 + arg1 + eps) / (ry2 + arg2 + eps))
  lg2 = log((ry1 + arg3 + eps) / (ry2 + arg4 + eps))
  lg3 = log((ry1 + arg7 + eps) / (ry2 + arg8 + eps))
  !lg4 = log((ry1 + arg5 + eps) / (ry2 + arg6 + eps))

  if (l_calcznodes) then
    arg1 = sqrt(ry1sq + R1)
    arg2 = sqrt(ry2sq + R1)
    arg5 = sqrt(ry1sq + R3)
    arg6 = sqrt(ry2sq + R3)

    lg1 = log((ry1 + arg1 + eps) / (ry2 + arg2 + eps))
    lg4 = log((ry1 + arg5 + eps) / (ry2 + arg6 + eps))

  else
    lg1 = znodes(znodes_i + 10)
    lg4 = znodes(znodes_i + 11)

  endif

  znodes(znodes_i + 10) = lg2
  znodes(znodes_i + 11) = lg3

  !print *, 'xz lg1', lg1
  !print *, 'xz lg2', lg2
  !print *, 'xz lg3', lg3
  !print *, 'xz lg4', lg4

  ! ts_xz
  !ts_x(3) = log((ry1 + arg1 + eps) / (ry2 + arg2 + eps)) - &
  !          log((ry1 + arg3 + eps) / (ry2 + arg4 + eps)) + &
  !          log((ry1 + arg7 + eps) / (ry2 + arg8 + eps)) - &
  !          log((ry1 + arg5 + eps) / (ry2 + arg6 + eps))

  ts_x(3) = lg1 - lg2 + lg3 - lg4

  ! mapping
  ! lg2 -> lg1
  ! lg3 -> lg4

  ! Filling the rest of the tensor.
  ! ts_zz
  ts_z(3) = -1 * (ts_x(1) + ts_y(2)) ! Gauss

  ! ts_zy
  ts_z(2) = ts_y(3)

  ! ts_xy
  ts_x(2) = ts_y(1)

  ! ts_zx
  ts_z(1) = ts_x(3)

  !print *, 'tsx', ts_x
  !print *, 'tsy', ts_y
  !print *, 'tsz', ts_z

end subroutine sharmbox

subroutine sharmbox2(x0, y0, z0, x1, y1, z1, x2, y2, z2, ts_x, ts_y, ts_z)!, znodes, l_calcznodes, xy_ind)!, nele_xylayer)
  real(kind=SENSIT_REAL), intent(in) :: x0, y0, z0, x1, y1, z1, x2, y2, z2
  real(kind=SENSIT_REAL), intent(out) :: ts_x(3), ts_y(3), ts_z(3)
  !logical, intent(in) :: l_calcznodes
  !integer, intent(in) :: xy_ind

  integer :: i, j, k, ind, sig
  real(kind=SENSIT_REAL) :: rx, ry, rz, rxsq, rysq, rzsq, r, eps
  real(kind=SENSIT_REAL) :: rx_arr(2), ry_arr(2), rz_arr(2)
  !real(kind=SENSIT_REAL) :: rxsq_arr(2), rysq_arr(2), rzsq_arr(2)
  real(kind=SENSIT_REAL) :: cur_node(5,8), tmp_sum(5)
  !real(kind=SENSIT_REAL), intent(inout) :: znodes(:)

  ! do something with unsued variables so the compiler is happy
  !if (l_calcznodes) then
  !  znodes = xy_ind
  !  znodes = -9999.9
  !endif

  eps = 1e-12

  rx_arr = (/ x1 - x0, x2 - x0 /)
  ry_arr = (/ y1 - y0, y2 - y0 /)
  rz_arr = (/ z1 - z0, z2 - z0 /)

  do i = 1, 2
    rx = rx_arr(i)
    rxsq = rx * rx

    do j = 1, 2
      ry = ry_arr(j)
      rysq = ry * ry

      do k = 1, 2
        ! i j k   ind corner  xx    yx    yy    yz    xz     sum sig
        ! ==========================================================
        ! 1 1 1   1   C1       at7  -lg32  at7  -lg41 -lg41  3   -
        ! 1 1 2   2   C5      -at6   lg31 -at6   lg31  lg31  4   +
        ! 1 2 1   3   C3      -at8   lg22 -at4   lg11  lg42  4   +
        ! 1 2 2   4   C7       at5  -lg21  at1  -lg21 -lg32  5   -
        ! 2 1 1   5   C2      -at4   lg42 -at8   lg42  lg11  4   +
        ! 2 1 2   6   C6       at1  -lg41  at5  -lg32 -lg21  5   -
        ! 2 2 1   7   C4       at3  -lg12  at3  -lg12 -lg12  5   -
        ! 2 2 2   8   C8      -at2   lg11 -at2   lg22  lg22  6   +

        ! atan2s: at1 - at2 + at3 - at4 + at5 - at6 + at7 - at8

        ! logs:   (lg11 - lg12) - (lg21 - lg22) + (lg31 - lg32) - (lg41 - lg42)

        ind = (i - 1) * 4 + (k - 1) * 2 + (j - 1) * 1 + 1
        sig = (-1) ** (i + j + k)

        rz = rz_arr(k)
        rzsq = rz * rz

        r = sqrt(rxsq + rysq + rzsq)

        cur_node(1, ind) = -sig * atan2(ry * rz, rx * r + eps) ! xx
        cur_node(2, ind) =  sig * log(rz + r + eps)            ! yx
        cur_node(3, ind) = -sig * atan2(rx * rz, ry * r + eps) ! yy
        cur_node(4, ind) =  sig * log(rx + r + eps)            ! yz
        cur_node(5, ind) =  sig * log(ry + r + eps)            ! xz
      enddo
    enddo
  enddo

  ! finalize tensor
  tmp_sum = SUM(cur_node, DIM=2)
  ts_x(1) = tmp_sum(1) ! xx
  ts_y(1) = tmp_sum(2) ! yx
  ts_y(2) = tmp_sum(3) ! yy
  ts_y(3) = tmp_sum(4) ! yz
  ts_x(3) = tmp_sum(5) ! xz

  ts_z(3) = -(ts_x(1) + ts_y(2)) ! zz = -(xx + yy)
  ts_z(2) = ts_y(3) ! zy = yz
  ts_x(2) = ts_y(1) ! xy = yx
  ts_z(1) = ts_x(3) ! zx = xz


end subroutine sharmbox2

! ======================================================================
! Calculates the kernels required for magnetic tensor mag
! When multiplied by the magnetisation vector, gives the derivative of the mag components in nT/m
!
! Units:
!   coordinates:        m
!   field intensity:    nT
!   incl/decl/azi:      deg
!   magnetisation:      Amp/m
!
! Inputs
!   x0, y0, z0    observation point coordinates
!   x1, x2        prism west-east coordinates respectively
!   y1, y2        prism south-north coordinates respectively
!   z1, z2        prism top-bottom coordinates respectively
!
! returns
!   mtensor_ijk   where i,j,k, are of the order {x,y,z} => {1,2,3}
!
! ======================================================================
subroutine tensorbox(x0, y0, z0, x1, y1, z1, x2, y2, z2, mtensor, tensorZnodes, l_calcznodes, xy_ind)
  real(kind=SENSIT_REAL), intent(in) :: x0, y0, z0, x1, y1, z1, x2, y2, z2
  real(kind=SENSIT_REAL), intent(out) :: mtensor(10)!mtensor(3, 3, 3)

  real(kind=SENSIT_REAL) :: temp_mtensor(10)!temp_mtensor(3, 3, 3)
  real(kind=SENSIT_REAL) :: rx, ry, rz
  real(kind=SENSIT_REAL) :: rx_sq, ry_sq, rz_sq, r
  real(kind=SENSIT_REAL) :: rx_arr(2), ry_arr(2), rz_arr(2)
  real(kind=SENSIT_REAL) :: rx_sq_arr(2), ry_sq_arr(2), rz_sq_arr(2)
  integer :: i, j, k

  real(kind=SENSIT_REAL), intent(inout) :: tensorZnodes(:)
  logical, intent(in) :: l_calcznodes
  integer, intent(in) :: xy_ind
  integer :: znodes_i, offset

  !print *, x1, x2
  !print *, y1, y2
  !print *, z1, z2

  print *, 'why are you in here'

  temp_mtensor = 0.0
  znodes_i = (xy_ind - 1) * 40 + 1

  ! Relative easting
  rx_arr(1) = (x2 - x0) ! East face
  rx_arr(2) = (x1 - x0) ! West face
  rx_sq_arr = rx_arr**2

  ! Relative northing
  ry_arr(1) = (y2 - y0) ! North face
  ry_arr(2) = (y1 - y0) ! South face
  ry_sq_arr = ry_arr**2

  ! Relative z
  ! Order of operations flipped to account for formula working in elevation space
  rz_arr(1) = (z0 - z1) ! Top face
  rz_arr(2) = (z0 - z2) ! Bottom face
  rz_sq_arr = rz_arr**2

  do i = 0, 1
    rx = rx_arr(i+1)
    rx_sq = rx_sq_arr(i+1)

    do j = 0, 1
      ry = ry_arr(j+1)
      ry_sq = ry_sq_arr(j+1)

      do k = 0, 1
        !print *, i, j, k

        ! i j
        ! 0 0 -> 0
        ! 0 1 -> 1
        ! 1 0 -> 2
        ! 1 1 -> 3
        offset = ((i * 2) + (j * 1)) * 10 + znodes_i

        ! control flags
        ! calc = T, k = 0 -> calc
        ! calc = T, k = 1 -> calc, save
        ! calc = F, k = 0 -> load
        ! calc = F, k = 1 -> calc, save

        if (.not. ((l_calcznodes == .false.) .and. (k == 0))) then

          rz = rz_arr(k+1)
          rz_sq = rz_sq_arr(k+1)

          ! Dist
          r = sqrt(rx_sq + ry_sq + rz_sq)

          ! Calculate the 10 required kernels
          ! order
          ! 1   2   3   4   5   6   7   8   9   10
          ! xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz

          ! u_xxx
          !temp_mtensor(1, 1, 1) = calc_tensor_iii(rx, ry, rz, r)
          temp_mtensor(1) = calc_tensor_iii(rx, ry, rz, r)

          ! u_xxy
          !temp_mtensor(1, 1, 2) = calc_tensor_iij(rx, ry, rz, r)
          temp_mtensor(2) = calc_tensor_iij(rx, ry, rz, r)

          ! u_xxz
          !temp_mtensor(1, 1, 3) = calc_tensor_iij(rx, rz, ry, r)
          temp_mtensor(3) = calc_tensor_iij(rx, rz, ry, r)

          ! u_xyy
          !temp_mtensor(1, 2, 2) = calc_tensor_iij(ry, rx, rz, r)
          temp_mtensor(4) = calc_tensor_iij(ry, rx, rz, r)

          ! u_xyz
          !temp_mtensor(1, 2, 3) = -1.0/r
          temp_mtensor(5) = -1.0/r

          ! u_xzz
          !temp_mtensor(1, 3, 3) = calc_tensor_iij(rz, rx, ry, r)
          temp_mtensor(6) = calc_tensor_iij(rz, rx, ry, r)

          ! u_yyy
          !temp_mtensor(2, 2, 2) = calc_tensor_iii(ry, rz, rx, r)
          temp_mtensor(7) = calc_tensor_iii(ry, rz, rx, r)

          ! u_yyz
          !temp_mtensor(2, 2, 3) = calc_tensor_iij(ry, rz, rx, r)
          temp_mtensor(8) = calc_tensor_iij(ry, rz, rx, r)

          ! u_yzz
          !temp_mtensor(2, 3, 3) = calc_tensor_iij(rz, ry, rx, r)
          temp_mtensor(9) = calc_tensor_iij(rz, ry, rx, r)

          ! u_zzz
          !temp_mtensor(3, 3, 3) = calc_tensor_iii(rz, ry, rx, r)
          temp_mtensor(10) = calc_tensor_iii(rz, ry, rx, r)

          if (k == 1) then
            !print *, 'save'
            !tensorZnodes(znodes_i + offset)     = temp_mtensor(1, 1, 1)
            !tensorZnodes(znodes_i + offset + 1) = temp_mtensor(1, 1, 2)
            !tensorZnodes(znodes_i + offset + 2) = temp_mtensor(1, 1, 3)
            !tensorZnodes(znodes_i + offset + 3) = temp_mtensor(1, 2, 2)
            !tensorZnodes(znodes_i + offset + 4) = temp_mtensor(1, 2, 3)
            !tensorZnodes(znodes_i + offset + 5) = temp_mtensor(1, 3, 3)
            !tensorZnodes(znodes_i + offset + 6) = temp_mtensor(2, 2, 2)
            !tensorZnodes(znodes_i + offset + 7) = temp_mtensor(2, 2, 3)
            !tensorZnodes(znodes_i + offset + 8) = temp_mtensor(2, 3, 3)
            !tensorZnodes(znodes_i + offset + 9) = temp_mtensor(3, 3, 3)
            tensorZnodes(offset : offset + 9) = temp_mtensor

          endif

        else
          !print *, 'load'
          !temp_mtensor(1, 1, 1) = tensorZnodes(znodes_i + offset)
          !temp_mtensor(1, 1, 2) = tensorZnodes(znodes_i + offset + 1)
          !temp_mtensor(1, 1, 3) = tensorZnodes(znodes_i + offset + 2)
          !temp_mtensor(1, 2, 2) = tensorZnodes(znodes_i + offset + 3)
          !temp_mtensor(1, 2, 3) = tensorZnodes(znodes_i + offset + 4)
          !temp_mtensor(1, 3, 3) = tensorZnodes(znodes_i + offset + 5)
          !temp_mtensor(2, 2, 2) = tensorZnodes(znodes_i + offset + 6)
          !temp_mtensor(2, 2, 3) = tensorZnodes(znodes_i + offset + 7)
          !temp_mtensor(2, 3, 3) = tensorZnodes(znodes_i + offset + 8)
          !temp_mtensor(3, 3, 3) = tensorZnodes(znodes_i + offset + 9)
          temp_mtensor = tensorZnodes(offset : offset + 9)

        endif
        !print *, 'test'
        !print *, 'xxx', temp_mtensor(1, 1, 1)
        !print *, 'xxy', temp_mtensor(1, 1, 2)
        !print *, 'xxz', temp_mtensor(1, 1, 3)
        !print *, 'xyy', temp_mtensor(1, 2, 2)
        !print *, 'xyz', temp_mtensor(1, 2, 3)
        !print *, 'xzz', temp_mtensor(1, 3, 3)
        !print *, 'yyy', temp_mtensor(2, 2, 2)
        !print *, 'yyz', temp_mtensor(2, 2, 3)
        !print *, 'yzz', temp_mtensor(2, 3, 3)
        !print *, 'zzz', temp_mtensor(3, 3, 3)

        ! Aggregate the tensor over all prisms for a specific data observation point
        mtensor = mtensor + (-1.0)**(i + j + k) * temp_mtensor

      enddo
    enddo
  enddo

end subroutine tensorbox

! ======================================================================
! evaluates the diagonal 3rd order kernels of a prism
! ======================================================================
pure function calc_tensor_iii(r_i, r_j, r_k, r) result(val)
  real(kind=SENSIT_REAL), intent(in) :: r_i, r_j, r_k, r
  real(kind=SENSIT_REAL) :: i_sq, j_sq, k_sq, num, den
  real(kind=SENSIT_REAL) :: val
  logical :: bool1, bool2

  bool1 = (r_i == 0.d0) .and. (r_j == 0.d0)
  bool2 = (r_i == 0.d0) .and. (r_k == 0.d0)

  if (bool1 .or. bool2) then
    val = 0.d0

  else
    i_sq = r_i * r_i
    j_sq = r_j * r_j
    k_sq = r_k * r_k

    num = -r_j * r_k * (2.0 * i_sq + j_sq + k_sq)
    den = (i_sq + j_sq) * (i_sq + k_sq) * r
    val = num / den

  endif

end function calc_tensor_iii

! ======================================================================
! evaluates the non-diagonal 3rd order kernels of a prism
! ======================================================================
pure function calc_tensor_iij(r_i, r_j, r_k, r) result(val)
  real(kind=SENSIT_REAL), intent(in) :: r_i, r_j, r_k, r
  real(kind=SENSIT_REAL) :: val

  if ((r_i == 0.0) .and. (r_j == 0.0)) then
    val = 0.0
  else
    val = -r_i / ((r_k + r) * r)
  endif

end function calc_tensor_iij

end module magnetic_field
