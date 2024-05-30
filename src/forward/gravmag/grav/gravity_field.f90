
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

module gravity_field

  use global_typedefs, only: CUSTOM_REAL, SIZE_REAL, PI
  use mpi_tools, only: exit_MPI
  use grid

  implicit none

  private

  ! Gravitational constant (SI units).
  real(kind=CUSTOM_REAL), parameter :: G_grav = 6.674e-11

  public :: graviprism_full
  public :: graviprism_z
  public :: gradiprism_full
  public :: gradiprism_zz

contains

!==========================================================================================
! Compute all components of gravity field.
! Surface topography is included and taken into account in the calculations.
! IMPORTANT: We purposely do all the calculations for "integral_computed" in double precision below,
! for accuracy reasons, and then convert the result back to single precision at the endif needed.
!==========================================================================================
subroutine graviprism_full(nelements, grid, Xdata, Ydata, Zdata, LineX, LineY, LineZ, myrank)
  integer, intent(in) :: nelements
  type(t_grid), intent(in) :: grid
  real(kind=CUSTOM_REAL), intent(in) :: Xdata, Ydata, Zdata
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: LineX(nelements)
  real(kind=CUSTOM_REAL), intent(out) :: LineY(nelements)
  real(kind=CUSTOM_REAL), intent(out) :: LineZ(nelements)

  ! Local variables.
  double precision, parameter :: twopi = 2.d0 * pi
  double precision, parameter, dimension(2) :: signo = (/-1.0d0, 1.0d0/)

  double precision :: Xdatad, Ydatad, Zdatad
  double precision :: lamb_over_sigmad
  double precision :: XX(2), YY(2), ZZ(2)
  double precision :: arg1, arg2, arg3, arg4, arg5, arg6
  double precision :: gx, gy, gz
  double precision :: dmu, Rs, sigmad
  integer :: i, k, l, m

  if (CUSTOM_REAL == SIZE_REAL) then
    Xdatad = dble(Xdata)
    Ydatad = dble(Ydata)
    Zdatad = dble(Zdata)
  else
    Xdatad = Xdata
    Ydatad = Ydata
    Zdatad = Zdata
  endif

  sigmad = 1.d0
  lamb_over_sigmad = G_grav / sigmad

  do i = 1, nelements

    XX(1) = Xdatad - dble(grid%X1(i))
    XX(2) = Xdatad - dble(grid%X2(i))
    YY(1) = Ydatad - dble(grid%Y1(i))
    YY(2) = Ydatad - dble(grid%Y2(i))
    ZZ(1) = Zdatad - dble(grid%Z1(i))
    ZZ(2) = Zdatad - dble(grid%Z2(i))

    gx = 0.d0
    gy = 0.d0
    gz = 0.d0

    do K=1,2
      do L=1,2
        do M=1,2
          dmu = signo(K) * signo(L) * signo(M)

          Rs = sqrt(XX(K)**2 + YY(L)**2 + ZZ(M)**2)

          arg1 = atan2(YY(L)*ZZ(M), (XX(K)*Rs))
          arg2 = atan2(XX(K)*ZZ(M), (YY(L)*Rs))
          arg3 = atan2(XX(K)*YY(L), (ZZ(M)*Rs))

          if (arg1 < 0) then
            arg1 = arg1 + twopi
          endif
          if (arg2 < 0) then
            arg2 = arg2 + twopi
          endif
          if (arg3 < 0) then
            arg3 = arg3 + twopi
          endif

          arg4 = Rs + XX(K)
          arg5 = Rs + YY(L)
          arg6 = Rs + ZZ(M)

          if (arg4 <= 0.) then
            call exit_MPI("Data coordinate coincides with model grid boundary (YZ). Adjust the model grid!", myrank, 0)
          endif
          if (arg5 <= 0.) then
            call exit_MPI("Data coordinate coincides with model grid boundary (XZ). Adjust the model grid!", myrank, 0)
          endif
          if (arg6 <= 0.) then
            call exit_MPI("Data coordinate coincides with model grid boundary (XY). Adjust the model grid!", myrank, 0)
          endif

          arg4 = log(arg4)
          arg5 = log(arg5)
          arg6 = log(arg6)

          gx = gx + dmu * (XX(K)*arg1 - YY(L)*arg6 - ZZ(M)*arg5)
          gy = gy + dmu * (YY(L)*arg2 - ZZ(M)*arg4 - XX(K)*arg6)
          gz = gz + dmu * (ZZ(M)*arg3 - XX(K)*arg5 - YY(L)*arg4)
        enddo
      enddo
    enddo

    ! Store element in full matrix line
    if (CUSTOM_REAL == SIZE_REAL) then
      LineX(i) = sngl(gx*lamb_over_sigmad)
      LineY(i) = sngl(gy*lamb_over_sigmad)
      LineZ(i) = sngl(gz*lamb_over_sigmad)
    else
      LineX(i) = gx*lamb_over_sigmad
      LineY(i) = gy*lamb_over_sigmad
      LineZ(i) = gz*lamb_over_sigmad
    endif

  enddo

end subroutine graviprism_full

!==========================================================================================
! Compute Z-components of gravity field (a reduced version of graviprism_full).
!==========================================================================================
subroutine graviprism_z(nelements, grid, Xdata, Ydata, Zdata, LineZ, myrank)
  integer, intent(in) :: nelements
  type(t_grid), intent(in) :: grid
  real(kind=CUSTOM_REAL), intent(in) :: Xdata, Ydata, Zdata
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: LineZ(nelements)

  ! Local variables.
  double precision, parameter :: twopi = 2.d0 * pi
  double precision, parameter, dimension(2) :: signo = (/-1.0d0, 1.0d0/)

  double precision :: Xdatad, Ydatad, Zdatad
  double precision :: lamb_over_sigmad
  double precision :: XX(2), YY(2), ZZ(2)
  double precision :: arg3, arg4, arg5
  double precision :: gz
  double precision :: dmu, Rs, sigmad
  integer :: i, k, l, m

  if (CUSTOM_REAL == SIZE_REAL) then
    Xdatad = dble(Xdata)
    Ydatad = dble(Ydata)
    Zdatad = dble(Zdata)
  else
    Xdatad = Xdata
    Ydatad = Ydata
    Zdatad = Zdata
  endif

  sigmad = 1.d0
  lamb_over_sigmad = G_grav / sigmad

  do i = 1, nelements

    XX(1) = Xdatad - dble(grid%X1(i))
    XX(2) = Xdatad - dble(grid%X2(i))
    YY(1) = Ydatad - dble(grid%Y1(i))
    YY(2) = Ydatad - dble(grid%Y2(i))
    ZZ(1) = Zdatad - dble(grid%Z1(i))
    ZZ(2) = Zdatad - dble(grid%Z2(i))

    gz = 0.d0

    do K=1,2
      do L=1,2
        do M=1,2
          dmu = signo(K) * signo(L) * signo(M)

          Rs = sqrt(XX(K)**2 + YY(L)**2 + ZZ(M)**2)

          arg3 = atan2(XX(K)*YY(L), (ZZ(M)*Rs))

          if (arg3 < 0) then
            arg3 = arg3 + twopi
          endif

          arg4 = Rs + XX(K)
          arg5 = Rs + YY(L)

          if (arg4 <= 0.) then
            call exit_MPI("Data coordinate coincides with model grid boundary (YZ). Adjust the model grid!", myrank, 0)
          endif
          if (arg5 <= 0.) then
            call exit_MPI("Data coordinate coincides with model grid boundary (XZ). Adjust the model grid!", myrank, 0)
          endif

          arg4 = log(arg4)
          arg5 = log(arg5)

          gz = gz + dmu * (ZZ(M)*arg3 - XX(K)*arg5 - YY(L)*arg4)
        enddo
      enddo
    enddo

    ! Store element in full matrix line
    if (CUSTOM_REAL == SIZE_REAL) then
      LineZ(i) = sngl(gz*lamb_over_sigmad)
    else
      LineZ(i) = gz*lamb_over_sigmad
    endif

  enddo

end subroutine graviprism_z

!==========================================================================================
! Compute the full gravity tensor required for gravity gradiometry
! Code adapted from matlab code provided in the paper Computation of the gravity field and
! its gradient: Some applications by Dubey and Tiwari in 2015
! expected units/conventions:
! distance unit: meters
! Z axis convention: positive down
! density: kg/m3
! output unit: m s-2 m-1
!==========================================================================================
subroutine gradiprism_full(nelements, grid, Xdata, Ydata, Zdata, LineXX, LineXY, LineYY, LineZX, LineYZ, LineZZ, myrank)
  integer, intent(in) :: nelements
  type(t_grid), intent(in) :: grid
  real(kind=CUSTOM_REAL), intent(in) :: Xdata, Ydata, Zdata
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: LineXX(nelements)
  real(kind=CUSTOM_REAL), intent(out) :: LineXY(nelements)
  real(kind=CUSTOM_REAL), intent(out) :: LineYY(nelements)
  real(kind=CUSTOM_REAL), intent(out) :: LineZX(nelements)
  real(kind=CUSTOM_REAL), intent(out) :: LineYZ(nelements)
  real(kind=CUSTOM_REAL), intent(out) :: LineZZ(nelements)

  ! Local variables.
  double precision, parameter :: twopi = 2.d0 * pi
  double precision, parameter, dimension(2) :: signo = (/-1.0d0, 1.0d0/)

  double precision :: Xdatad, Ydatad, Zdatad
  double precision :: lamb_over_sigmad
  double precision :: XX(2), YY(2), ZZ(2)
  double precision :: vxx, vxy, vyy, vzx, vyz, vzz
  double precision :: gxx, gxy, gyy, gzx, gyz, gzz
  double precision :: dmu, Rs, sigmad
  integer :: i, k, l, m

  if (CUSTOM_REAL == SIZE_REAL) then
    Xdatad = dble(Xdata)
    Ydatad = dble(Ydata)
    Zdatad = dble(Zdata)
  else
    Xdatad = Xdata
    Ydatad = Ydata
    Zdatad = Zdata
  endif

  sigmad = 1.d0
  lamb_over_sigmad = G_grav / sigmad

  do i = 1, nelements

    XX(1) = Xdatad - dble(grid%X1(i))
    XX(2) = Xdatad - dble(grid%X2(i))
    YY(1) = Ydatad - dble(grid%Y1(i))
    YY(2) = Ydatad - dble(grid%Y2(i))
    ZZ(1) = -( Zdatad - dble(grid%Z1(i)) )
    ZZ(2) = -( Zdatad - dble(grid%Z2(i)) )

    gxx = 0.d0
    gxy = 0.d0
    gyy = 0.d0
    gzx = 0.d0
    gyz = 0.d0
    gzz = 0.d0

    do K=1,2
      do L=1,2
        do M=1,2
          dmu = signo(K) * signo(L) * signo(M)

          Rs = sqrt(XX(K)**2 + YY(L)**2 + ZZ(M)**2)

          vxx = atan((XX(K)*YY(L))/((XX(K)**2)+Rs*ZZ(M)+ZZ(M)**2))
          vxy = log(Rs+ZZ(M))
          vyy = atan((XX(K)*YY(L))/(Rs**2+Rs*ZZ(M)-XX(K)**2))
          vzx = 0.5*log((Rs-YY(L)/(Rs+YY(L))))
          vyz = 0.5*log((Rs-XX(K)/(Rs+XX(K))))
          vzz = -1.*atan(XX(K)*YY(L)/(Rs*ZZ(M)))

          gxx = gxx + dmu*vxx
          gxy = gxy + dmu*vxy
          gyy = gyy + dmu*vyy
          gzx = gzx + dmu*vzx
          gyz = gyz + dmu*vyz
          gzz = gzz + dmu*vzz

        enddo
      enddo
    enddo

    ! Store element in full matrix line
    if (CUSTOM_REAL == SIZE_REAL) then
      LineXX(i) = sngl(gxx*lamb_over_sigmad)
      LineXY(i) = sngl(gxy*lamb_over_sigmad)
      LineYY(i) = sngl(gyy*lamb_over_sigmad)
      LineZX(i) = sngl(gzx*lamb_over_sigmad)
      LineYZ(i) = sngl(gyz*lamb_over_sigmad)
      LineZZ(i) = sngl(gzz*lamb_over_sigmad)
    else
      LineXX(i) = gxx*lamb_over_sigmad
      LineXY(i) = gxy*lamb_over_sigmad
      LineYY(i) = gyy*lamb_over_sigmad
      LineZX(i) = gzx*lamb_over_sigmad
      LineYZ(i) = gyz*lamb_over_sigmad
      LineZZ(i) = gzz*lamb_over_sigmad
    endif

  enddo

end subroutine gradiprism_full

!==========================================================================================
! Reduced version of gradiprism_full
!==========================================================================================
subroutine gradiprism_zz(nelements, grid, Xdata, Ydata, Zdata, LineZZ, myrank)
  integer, intent(in) :: nelements
  type(t_grid), intent(in) :: grid
  real(kind=CUSTOM_REAL), intent(in) :: Xdata, Ydata, Zdata
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: LineZZ(nelements)

  ! Local variables.
  double precision, parameter :: twopi = 2.d0 * pi
  double precision, parameter, dimension(2) :: signo = (/-1.0d0, 1.0d0/)

  double precision :: Xdatad, Ydatad, Zdatad
  double precision :: lamb_over_sigmad
  double precision :: XX(2), YY(2), ZZ(2)
  double precision :: vzz
  double precision :: gzz
  double precision :: dmu, Rs, sigmad
  integer :: i, k, l, m

  if (CUSTOM_REAL == SIZE_REAL) then
    Xdatad = dble(Xdata)
    Ydatad = dble(Ydata)
    Zdatad = dble(Zdata)
  else
    Xdatad = Xdata
    Ydatad = Ydata
    Zdatad = Zdata
  endif

  sigmad = 1.d0
  lamb_over_sigmad = G_grav / sigmad

  do i = 1, nelements

    XX(1) = Xdatad - dble(grid%X1(i))
    XX(2) = Xdatad - dble(grid%X2(i))
    YY(1) = Ydatad - dble(grid%Y1(i))
    YY(2) = Ydatad - dble(grid%Y2(i))
    ZZ(1) = -( Zdatad - dble(grid%Z1(i)) )
    ZZ(2) = -( Zdatad - dble(grid%Z2(i)) )

    gzz = 0.d0

    do K=1,2
      do L=1,2
        do M=1,2
          dmu = signo(K) * signo(L) * signo(M)

          Rs = sqrt(XX(K)**2 + YY(L)**2 + ZZ(M)**2)

          vzz = -1.*atan2(XX(K)*YY(L), (Rs*ZZ(M)))

          gzz = gzz + dmu*vzz

        enddo
      enddo
    enddo

    ! Store element in full matrix line
    if (CUSTOM_REAL == SIZE_REAL) then
      LineZZ(i) = sngl(gzz*lamb_over_sigmad)
    else
      LineZZ(i) = gzz*lamb_over_sigmad
    endif

  enddo

end subroutine gradiprism_zz

end module gravity_field
