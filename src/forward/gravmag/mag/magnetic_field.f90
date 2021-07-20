
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

!========================================================================================
! A class to calculate sensitivity values for magnetic field.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!========================================================================================
module magnetic_field

  use global_typedefs
  use grid

  implicit none

  private

  ! The vacuum magnetic permeability (SI).
  double precision, parameter :: mu0 = PI * 4.d-7

  ! Degrees to radians conversion factor.
  double precision, parameter :: d2rad = PI / 180.d0

  ! Main class.
  type, public :: t_magnetic_field
    private

    ! For optimization of mbox() function calculation precompute these parameters.
    double precision :: fm1, fm2, fm3, fm4, fm5, fm6

    ! External/background (average regional) field intensity (nT).
    double precision :: intensity

  contains
    private

    procedure, public, pass :: initialize => magnetic_field_initialize
    procedure, public, pass :: magprism => magnetic_field_magprism

    procedure, private, pass :: mbox
    procedure, private, nopass :: dircos

  end type t_magnetic_field

contains

!==============================================================================
! Initialization.
! See parameters description in mbox() function.
!==============================================================================
subroutine magnetic_field_initialize(this, mi, md, fi, fd, theta, intensity)
  class(t_magnetic_field), intent(inout) :: this
  double precision, intent(in) :: mi, md, fi, fd, theta, intensity

  double precision :: ma, mb, mc
  double precision :: fa, fb, fc

  this%intensity = intensity

  call this%dircos(mi, md, theta, ma, mb, mc)
  call this%dircos(fi, fd, theta, fa, fb, fc)

  this%fm1 = ma * fb + mb * fa
  this%fm2 = ma * fc + mc * fa
  this%fm3 = mb * fc + mc * fb
  this%fm4 = ma * fa
  this%fm5 = mb * fb
  this%fm6 = mc * fc

end subroutine magnetic_field_initialize

!====================================================================================================
! Calculates the sensitivity kernel.
!====================================================================================================
subroutine magnetic_field_magprism(this, nelements, data_j, grid, Xdata, Ydata, Zdata, sensit_line)
  class(t_magnetic_field), intent(in) :: this
  integer, intent(in) :: nelements, data_j
  type(t_grid), intent(in) :: grid
  real(kind=CUSTOM_REAL), intent(in) :: Xdata(:), Ydata(:), Zdata(:)

  real(kind=CUSTOM_REAL), intent(out) :: sensit_line(:)

  integer :: i
  real(kind=CUSTOM_REAL) :: t1, t2

  do i = 1, nelements

    call this%mbox(Xdata(data_j), Ydata(data_j), Zdata(data_j), &
                   grid%X1(i), grid%Y1(i), grid%Z1(i), grid%X2(i), grid%Y2(i), t1)

    call this%mbox(Xdata(data_j), Ydata(data_j), Zdata(data_j), &
                   grid%X1(i), grid%Y1(i), grid%Z2(i), grid%X2(i), grid%Y2(i), t2)

    sensit_line(i) = (t1 - t2) / (4.d0 * PI)

  enddo

end subroutine magnetic_field_magprism

!===================================================================================
!  Subroutine MBOX computes the total field anomaly of an infinitely
!  extended rectangular prism. Sides of prism are parallel to x,y,z
!  axes, and z is vertical down. Bottom of prism extends to infinity.
!  Two calls to mbox can provide the anomaly of a prism with finite
!  thickness; e.g.,
!
!     call mbox(x0,y0,z0,x1,y1,z1,x2,y2,mi,md,fi,fd,m,theta,t1)
!     call mbox(x0,y0,z0,x1,y1,z2,x2,y2,mi,md,fi,fd,m,theta,t2)
!     t=t1-t2
!
!  Requires subroutine DIRCOS.
!  Method from Ref. [1]. For more details see also Refs. [2-3].
!  [1] B. K. Bhattacharyya, Magnetic anomalies due to prism-shaped bodies with arbitrary polarization, Geophysics, 29, 1964.
!  [2] B. K. Bhattacharyya, Computer modeling in gravity and magnetic interpretation, Geophysics, 43, 1978.
!  [3] B. K. Bhattacharyya, A generalized multibody model for inversion of magnetic anomalies, Geophysics, 45, 1980.
!
!  Input parameters:
!    Observation point is (x0,y0,z0). Prism extends from x1 to
!    x2, y1 to y2, and z1 to infinity in x, y, and z directions,
!    respectively.  Magnetization defined by inclination mi,
!    declination md, intensity m.  Ambient field defined by
!    inclination fi and declination fd.  X axis has declination
!    theta. Distance units are irrelevant but must be consistent.
!    Angles are in degrees, with inclinations positive below
!    horizontal and declinations positive east of true north.
!    Magnetization in A/m.
!
!  Output parameters:
!    Total field anomaly t, in nT.
!===================================================================================
subroutine mbox(this, x0, y0, z0, x1, y1, z1, x2, y2, t)
  class(t_magnetic_field), intent(in) :: this
  double precision, intent(in) :: x0, y0, z0, x1, y1, z1, x2, y2
  double precision, intent(out) :: t

  ! Local variables.
  double precision :: alpha(2), beta(2)
  double precision :: alphabeta, alphasq, arg1, arg2, arg3, arg4
  double precision :: h, hsq, r0, r0h, r0sq, sign, tatan, tlog
  integer :: i, j

  alpha(1) = x1 - x0
  alpha(2) = x2 - x0

  beta(1) = y1 - y0
  beta(2) = y2 - y0

  h = z1 - z0

  ! Sanity check.
  if (h < 0.d0) then
    print *, "Error in mbox: negative depth h =", h
    stop
  endif

  t = 0.
  hsq = h**2

  do i = 1, 2
    alphasq = alpha(i)**2
    do j = 1, 2
      sign = 1.
      if (i /= j) sign = -1.

      r0sq = alphasq + beta(j)**2 + hsq
      r0 = dsqrt(r0sq)
      r0h = r0 * h
      alphabeta = alpha(i) * beta(j)

      arg1 = (r0 - alpha(i)) / (r0 + alpha(i))
      arg2 = (r0 - beta(j)) / (r0 + beta(j))
      arg3 = alphasq + r0h + hsq
      arg4 = r0sq + r0h - alphasq

      tlog = this%fm3 * dlog(arg1) / 2. + this%fm2 * dlog(arg2) / 2. - this%fm1 * dlog(r0 + h)
      tatan = - this%fm4 * datan2(alphabeta, arg3) - this%fm5 * datan2(alphabeta, arg4) + this%fm6 * datan2(alphabeta, r0h)

      t = t + sign * (tlog + tatan)
    enddo
  enddo

  !t = t * mu0 * this%intensity
  t = t * this%intensity

end subroutine mbox


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
!    a,b,c:  the three direction cosines.
!==============================================================================
subroutine dircos(incl, decl, azim, a, b, c)
  double precision, intent(in) :: incl, decl, azim
  double precision, intent(out) :: a, b, c

  double precision :: xincl, xdecl, xazim

  xincl = incl * d2rad
  xdecl = decl * d2rad
  xazim = azim * d2rad

  a = dcos(xincl) * dcos(xdecl - xazim)
  b = dcos(xincl) * dsin(xdecl - xazim)
  c = dsin(xincl)

end subroutine dircos

end module magnetic_field
