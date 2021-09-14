
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

contains

!==========================================================================================
! Compute all components of gravity field.
! Surface topography is included and taken into account in the calculations.
! IMPORTANT: We purposely do all the calculations for "integral_computed" in double precision below,
! for accuracy reasons, and then convert the result back to single precision at the endif needed.
!==========================================================================================
subroutine graviprism_full(nelements, ncomponents, grid, Xdata, Ydata, Zdata, LineX, LineY, LineZ, myrank)
  integer, intent(in) :: nelements, ncomponents
  type(t_grid), intent(in) :: grid
  real(kind=CUSTOM_REAL), intent(in) :: Xdata, Ydata, Zdata
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: LineX(:)
  real(kind=CUSTOM_REAL), intent(out) :: LineY(:)
  real(kind=CUSTOM_REAL), intent(out) :: LineZ(:)

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
            call exit_MPI("Bad point field (arg4).", myrank, 0)
          endif
          if (arg5 <= 0.) then
            call exit_MPI("Bad point field (arg5).", myrank, 0)
          endif
          if (arg6 <= 0.) then
            call exit_MPI("Bad point field (arg6).", myrank, 0)
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
      LineZ(i) = sngl(gz*lamb_over_sigmad)

      if (ncomponents == 3) then
        LineX(i) = sngl(gx*lamb_over_sigmad)
        LineY(i) = sngl(gy*lamb_over_sigmad)
      endif
    else
      LineZ(i) = gz*lamb_over_sigmad

      if (ncomponents == 3) then
        LineX(i) = gx*lamb_over_sigmad
        LineY(i) = gy*lamb_over_sigmad
      endif
    endif

  enddo

end subroutine graviprism_full

end module gravity_field

