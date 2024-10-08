
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

  double precision :: XX(2), YY(2), ZZ(2)
  double precision :: arg1, arg2, arg3, arg4, arg5, arg6
  double precision :: gx, gy, gz
  double precision :: dmu, Rs
  integer :: i, k, l, m

  do i = 1, nelements

    XX(1) = dble(Xdata - grid%X1(i))
    XX(2) = dble(Xdata - grid%X2(i))
    YY(1) = dble(Ydata - grid%Y1(i))
    YY(2) = dble(Ydata - grid%Y2(i))
    ZZ(1) = dble(Zdata - grid%Z1(i))
    ZZ(2) = dble(Zdata - grid%Z2(i))

    gx = 0.d0
    gy = 0.d0
    gz = 0.d0

    do K = 1, 2
      do L = 1, 2
        do M = 1, 2
          dmu = signo(K) * signo(L) * signo(M)

          Rs = sqrt(XX(K)**2 + YY(L)**2 + ZZ(M)**2)

          arg1 = atan2(YY(L) * ZZ(M), XX(K) * Rs)
          arg2 = atan2(XX(K) * ZZ(M), YY(L) * Rs)
          arg3 = atan2(XX(K) * YY(L), ZZ(M) * Rs)

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

          gx = gx + dmu * (XX(K) * arg1 - YY(L) * arg6 - ZZ(M) * arg5)
          gy = gy + dmu * (YY(L) * arg2 - ZZ(M) * arg4 - XX(K) * arg6)
          gz = gz + dmu * (ZZ(M) * arg3 - XX(K) * arg5 - YY(L) * arg4)
        enddo
      enddo
    enddo

    ! Store element in full matrix line
    LineX(i) = real(G_grav * gx, CUSTOM_REAL)
    LineY(i) = real(G_grav * gy, CUSTOM_REAL)
    LineZ(i) = real(G_grav * gz, CUSTOM_REAL)
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

  double precision :: XX(2), YY(2), ZZ(2)
  double precision :: arg3, arg4, arg5
  double precision :: gz
  double precision :: dmu, Rs
  integer :: i, k, l, m

  do i = 1, nelements

    XX(1) = dble(Xdata - grid%X1(i))
    XX(2) = dble(Xdata - grid%X2(i))
    YY(1) = dble(Ydata - grid%Y1(i))
    YY(2) = dble(Ydata - grid%Y2(i))
    ZZ(1) = dble(Zdata - grid%Z1(i))
    ZZ(2) = dble(Zdata - grid%Z2(i))

    gz = 0.d0

    do K = 1, 2
      do L = 1, 2
        do M = 1, 2
          dmu = signo(K) * signo(L) * signo(M)

          Rs = sqrt(XX(K)**2 + YY(L)**2 + ZZ(M)**2)

          arg3 = atan2(XX(K) * YY(L), ZZ(M) * Rs)

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

          gz = gz + dmu * (ZZ(M) * arg3 - XX(K) * arg5 - YY(L) * arg4)
        enddo
      enddo
    enddo

    ! Store element in full matrix line
    LineZ(i) = real(G_grav * gz, CUSTOM_REAL)
  enddo

end subroutine graviprism_z

!==============================================================================================
! Compute the full gravity tensor required for gravity gradiometry
! Code adapted from matlab code provided in the paper Computation of the gravity field and
! its gradient: Some applications by Dubey and Tiwari in 2015.
! Expected units/conventions:
!   distance unit: meters
!   Z axis convention: positive down
!   density: kg/m3
!   output unit: m s-2 m-1
!==============================================================================================
subroutine gradiprism_full(nelements, grid, Xdata, Ydata, Zdata, LineXX, LineYY, LineZZ, LineXY, LineYZ, LineZX, myrank)
  integer, intent(in) :: nelements
  type(t_grid), intent(in) :: grid
  real(kind=CUSTOM_REAL), intent(in) :: Xdata, Ydata, Zdata
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: LineXX(nelements)
  real(kind=CUSTOM_REAL), intent(out) :: LineYY(nelements)
  real(kind=CUSTOM_REAL), intent(out) :: LineZZ(nelements)
  real(kind=CUSTOM_REAL), intent(out) :: LineXY(nelements)
  real(kind=CUSTOM_REAL), intent(out) :: LineYZ(nelements)
  real(kind=CUSTOM_REAL), intent(out) :: LineZX(nelements)

  ! Local variables.
  double precision, parameter :: twopi = 2.d0 * pi
  double precision, parameter, dimension(2) :: signo = (/-1.0d0, 1.0d0/)

  double precision :: XX(2), YY(2), ZZ(2)
  double precision :: vxx, vxy, vyy, vzx, vyz, vzz
  double precision :: gxx, gxy, gyy, gzx, gyz, gzz
  double precision :: arg1, arg2, arg3
  double precision :: arg21, arg22, arg31, arg32
  double precision :: dmu, Rs
  integer :: i, k, l, m

  do i = 1, nelements

    XX(1) = dble(Xdata - grid%X1(i))
    XX(2) = dble(Xdata - grid%X2(i))
    YY(1) = dble(Ydata - grid%Y1(i))
    YY(2) = dble(Ydata - grid%Y2(i))
    ZZ(1) = -dble(Zdata - grid%Z1(i))
    ZZ(2) = -dble(Zdata - grid%Z2(i))

    gxx = 0.d0
    gxy = 0.d0
    gyy = 0.d0
    gzx = 0.d0
    gyz = 0.d0
    gzz = 0.d0

    do K = 1, 2
      do L = 1, 2
        do M = 1, 2
          dmu = signo(K) * signo(L) * signo(M)

          Rs = sqrt(XX(K)**2 + YY(L)**2 + ZZ(M)**2)

          vxx = atan2(XX(K) * YY(L), XX(K)**2 + Rs * ZZ(M) + ZZ(M)**2)
          vyy = atan2(XX(K) * YY(L), Rs**2 + Rs * ZZ(M) - XX(K)**2)
          vzz = -atan2(XX(K) * YY(L), Rs * ZZ(M))

          if (vxx < 0) then
            vxx = vxx + twopi
          endif
          if (vyy < 0) then
            vyy = vyy + twopi
          endif
          if (vzz < 0) then
            vzz = vzz + twopi
          endif

          arg1 = Rs + ZZ(M)
          arg21 = Rs - YY(L)
          arg22 = Rs + YY(L)
          arg31 = Rs - XX(K)
          arg32 = Rs + XX(K)

          if (arg22 == 0. .or. arg32 == 0.) then
            call exit_MPI("Zero denominator in gradiprism_full! Adjust the model grid.", myrank, 0)
          endif

          arg2 = arg21 / arg22
          arg3 = arg31 / arg32

          if (arg1 <= 0. .or. arg2 <= 0. .or. arg3 <= 0.) then
            call exit_MPI("Bad log argument in gradiprism_full! Adjust the model grid.", myrank, 0)
          endif

          vxy = log(arg1)
          vzx = 0.5d0 * log(arg2)
          vyz = 0.5d0 * log(arg3)

          gxx = gxx + dmu * vxx
          gyy = gyy + dmu * vyy
          gzz = gzz + dmu * vzz
          gxy = gxy + dmu * vxy
          gyz = gyz + dmu * vyz
          gzx = gzx + dmu * vzx
        enddo
      enddo
    enddo

    ! Store element in full matrix line
    LineXX(i) = real(G_grav * gxx, CUSTOM_REAL)
    LineXY(i) = real(G_grav * gxy, CUSTOM_REAL)
    LineYY(i) = real(G_grav * gyy, CUSTOM_REAL)
    LineZX(i) = real(G_grav * gzx, CUSTOM_REAL)
    LineYZ(i) = real(G_grav * gyz, CUSTOM_REAL)
    LineZZ(i) = real(G_grav * gzz, CUSTOM_REAL)
  enddo

end subroutine gradiprism_full

!==========================================================================================
! Reduced version of gradiprism_full
!==========================================================================================
subroutine gradiprism_zz(nelements, grid, Xdata, Ydata, Zdata, LineZZ)
  integer, intent(in) :: nelements
  type(t_grid), intent(in) :: grid
  real(kind=CUSTOM_REAL), intent(in) :: Xdata, Ydata, Zdata

  real(kind=CUSTOM_REAL), intent(out) :: LineZZ(nelements)

  ! Local variables.
  double precision, parameter :: twopi = 2.d0 * pi
  double precision, parameter, dimension(2) :: signo = (/-1.0d0, 1.0d0/)

  double precision :: XX(2), YY(2), ZZ(2)
  double precision :: vzz
  double precision :: gzz
  double precision :: dmu, Rs
  integer :: i, k, l, m

  do i = 1, nelements

    XX(1) = dble(Xdata - grid%X1(i))
    XX(2) = dble(Xdata - grid%X2(i))
    YY(1) = dble(Ydata - grid%Y1(i))
    YY(2) = dble(Ydata - grid%Y2(i))
    ZZ(1) = -dble(Zdata - grid%Z1(i))
    ZZ(2) = -dble(Zdata - grid%Z2(i))

    gzz = 0.d0

    do K = 1, 2
      do L = 1, 2
        do M = 1, 2
          dmu = signo(K) * signo(L) * signo(M)

          Rs = sqrt(XX(K)**2 + YY(L)**2 + ZZ(M)**2)

          vzz = -atan2(XX(K) * YY(L), Rs * ZZ(M))

          if (vzz < 0) then
            vzz = vzz + twopi
          endif

          gzz = gzz + dmu * vzz
        enddo
      enddo
    enddo

    ! Store element in full matrix line
    LineZZ(i) = real(G_grav * gzz, CUSTOM_REAL)
  enddo

end subroutine gradiprism_zz

end module gravity_field
