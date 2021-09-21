
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

!===============================================================================================
! Unit tests for ECT (forward problem) part.
!
! Author: Vitaliy Ogarko, UWA, CET, Australia, 2015.
!===============================================================================================
module tests_ect

  use global_typedefs
  use ftnunit

  use parameters_ect
  use geometry
  use flag
  use matrix
  use solution
  use conjugate_gradient
  use laplace
  use utils
  use paraview

  implicit none

  private

  public :: test_geometry_all
  public :: test_boundary_conditions_all
  public :: test_analytical_comparison_all

  private :: test_geometry
  private :: test_boundary_conditions
  private :: test_analytical_comparison

  private :: init_params
  private :: set_dimensions_multigrid

  ! TODO: make them local.
  integer :: myrank, nbproc

contains

subroutine test_geometry_all(myrank__,nbproc__)
  integer, intent(in) :: myrank__,nbproc__

  ! Geometry of the cylinder.
  real(kind=CUSTOM_REAL), parameter :: radiusin     = 0.045
  real(kind=CUSTOM_REAL), parameter :: radiusout    = 0.06
  real(kind=CUSTOM_REAL), parameter :: radiusoutout = 0.07
  real(kind=CUSTOM_REAL), parameter :: heicyl       = 0.2

  myrank = myrank__
  nbproc = nbproc__

  ! Minimum problem size without multigrid.
  call test_geometry(12, 1, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(12, 1, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))

  ! Minimum problem size with multigrid 2 levels.
  call test_geometry(24, 2, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(24, 2, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))

  call test_geometry(36, 1, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(36, 1, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.))
  call test_geometry(36, 1, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))
  call test_geometry(36, 1, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.01))

  ! Minimum problem size with multigird 3 levels.
  call test_geometry(48, 3, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(48, 3, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(48, 3, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))

  call test_geometry(72, 2, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(72, 2, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(72, 2, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.))
  call test_geometry(72, 2, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.))
  call test_geometry(72, 2, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))
  call test_geometry(72, 2, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))
  call test_geometry(72, 2, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.01))
  call test_geometry(72, 2, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.01))

  ! Minimum problem size with multigird 4 levels.
  call test_geometry(96, 4, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(96, 4, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(96, 4, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))
  call test_geometry(96, 4, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.02))

  call test_geometry(144, 3, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(144, 3, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(144, 3, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.))
  call test_geometry(144, 3, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.))
  call test_geometry(144, 3, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))
  call test_geometry(144, 3, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))
  call test_geometry(144, 3, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.01))
  call test_geometry(144, 3, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.01))

  call test_geometry(288, 3, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(288, 3, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(288, 3, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.))
  call test_geometry(288, 3, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.))
  call test_geometry(288, 3, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))
  call test_geometry(288, 3, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))
  call test_geometry(288, 3, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.01))
  call test_geometry(288, 3, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.01))

  call test_geometry(288, 4, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(288, 4, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_geometry(288, 4, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.))
  call test_geometry(288, 4, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.))
  call test_geometry(288, 4, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))
  call test_geometry(288, 4, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))
  call test_geometry(288, 4, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.01))
  call test_geometry(288, 4, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.01))

end subroutine test_geometry_all

!==============================================================================================
subroutine test_boundary_conditions_all(myrank__,nbproc__)
  integer, intent(in) :: myrank__,nbproc__

  ! Geometry of the cylinder.
  real(kind=CUSTOM_REAL), parameter :: radiusin     = 0.045
  real(kind=CUSTOM_REAL), parameter :: radiusout    = 0.06
  real(kind=CUSTOM_REAL), parameter :: radiusoutout = 0.07
  real(kind=CUSTOM_REAL), parameter :: heicyl       = 0.2

  myrank = myrank__
  nbproc = nbproc__

  call test_boundary_conditions(12, 1, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_boundary_conditions(12, 1, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))

  call test_boundary_conditions(24, 2, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_boundary_conditions(24, 2, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))

  call test_boundary_conditions(72, 2, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_boundary_conditions(72, 2, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.))
  call test_boundary_conditions(72, 2, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.))
  call test_boundary_conditions(72, 2, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.))
  call test_boundary_conditions(72, 2, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))
  call test_boundary_conditions(72, 2, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.01))
  call test_boundary_conditions(72, 2, 0, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.01))
  call test_boundary_conditions(72, 2, 1, t_sensor(radiusin,radiusout,radiusoutout,heicyl,3.,0.01))

end subroutine test_boundary_conditions_all

!==============================================================================================
subroutine test_analytical_comparison_all(myrank__,nbproc__)
  integer, intent(in) :: myrank__,nbproc__

  ! Geometry of the cylinder.
  real(kind=CUSTOM_REAL), parameter :: radiusin     = 0.045
  real(kind=CUSTOM_REAL), parameter :: radiusout    = 0.06
  real(kind=CUSTOM_REAL), parameter :: radiusoutout = 0.07
  real(kind=CUSTOM_REAL), parameter :: heicyl       = 0.2

  ! Analytical solutions are known only without gaps.
  type (t_sensor), parameter :: sensor = t_sensor(radiusin,radiusout,radiusoutout,heicyl,0.,0.)

  myrank = myrank__
  nbproc = nbproc__

  !call test_analytical_comparison(1, 24, 1, 1000, sensor)
  call test_analytical_comparison(1, 36, 1, 1000, sensor)
  !call test_analytical_comparison(1, 72, 1, 1000, sensor)
  !call test_analytical_comparison(1, 144, 1, 3000, sensor)

  !call test_analytical_comparison(2, 24, 1, 1000, sensor)
  !call test_analytical_comparison(2, 36, 1, 1000, sensor)
  !call test_analytical_comparison(2, 72, 1, 1000, sensor)
  !call test_analytical_comparison(2, 144, 1, 3000, sensor)

  !call test_analytical_comparison(3, 24, 1, 1000, sensor)
  !call test_analytical_comparison(3, 72, 1, 1000, sensor)
  !call test_analytical_comparison(3, 144, 1, 3000, sensor)

  !call test_analytical_comparison(4, 24, 1, 1000, sensor)
  !call test_analytical_comparison(4, 36, 1, 1000, sensor)
  !call test_analytical_comparison(4, 72, 1, 1000, sensor)
  !call test_analytical_comparison(4, 144, 1, 3000, sensor)
  !call test_analytical_comparison(4, 288, 1, 10000, sensor)

end subroutine test_analytical_comparison_all

!==============================================================================================
! Initialize the system parameters for testing functions.
!==============================================================================================
subroutine init_params(size, nbproc, par, nr_read, nz_read, nzlocal_read, ntheta_read, kguards_read, ifixed_elecgeo)

  integer, intent(in) :: size, nbproc, ifixed_elecgeo

  type(t_parameters_ect), intent(out) :: par
  integer, intent(out) :: nr_read, nz_read, nzlocal_read, ntheta_read, kguards_read

  nr_read                   = size
  ntheta_read               = size
  nz_read                   = size
  kguards_read              = size/4
  nzlocal_read              = nz_read/nbproc

  par%dims%nr = nr_read
  par%dims%ntheta = ntheta_read
  par%dims%nz = nz_read
  par%dims%kguards = kguards_read
  par%dims%nzlocal = nzlocal_read

  par%nel = 12
  par%nrings = 1

  par%ifixed_elecgeo = ifixed_elecgeo

end subroutine init_params

!==============================================================================================
! Used to initialize the multigrid dimensions for testing functions.
!==============================================================================================
subroutine set_dimensions_multigrid(ilevel_coarse, nr_read, nz_read, nzlocal_read, ntheta_read, kguards_read, &
                                nr_at_level, ntheta_at_level, nz_at_level, nzlocal_at_level, kguards_at_level)

  integer, intent(in) :: ilevel_coarse

  integer, intent(out) :: nr_at_level(ilevel_coarse)
  integer, intent(out) :: ntheta_at_level(ilevel_coarse)
  integer, intent(out) :: nz_at_level(ilevel_coarse)
  integer, intent(out) :: nzlocal_at_level(ilevel_coarse)
  integer, intent(out) :: kguards_at_level(ilevel_coarse)

  integer :: nr_read, nz_read, nzlocal_read, ntheta_read, kguards_read
  integer :: ilevel
  ! For the multigrid method start from level=1.
  integer, parameter :: ilevel_fine = 1

  do ilevel = ilevel_fine, ilevel_coarse
    nr_at_level(ilevel)      = nr_read / 2**(ilevel-1)
    ntheta_at_level(ilevel)  = ntheta_read / 2**(ilevel-1)
    nz_at_level(ilevel)      = nz_read / 2**(ilevel-1)
    nzlocal_at_level(ilevel) = nzlocal_read / 2**(ilevel-1)
    kguards_at_level(ilevel) = kguards_read / 2**(ilevel-1)
  enddo

end subroutine set_dimensions_multigrid


!==============================================================================================
! Tests of the system geometry.
!==============================================================================================
subroutine test_geometry(size,ilevel_coarse,ifixed_elecgeo,sens)

  integer, intent(in) :: size
  integer, intent(in) :: ilevel_coarse, ifixed_elecgeo
  type(t_sensor), intent(in) :: sens

  type(t_parameters_ect) :: par
  integer :: ilevel
  integer :: nr_read, nz_read, nzlocal_read, ntheta_read, kguards_read
  integer, dimension(:), allocatable :: i1, i2
  type(t_electrode_p), dimension(:), allocatable :: elecs
  type(t_guards), dimension(:), allocatable :: guards
  type(real_p), dimension(:), allocatable :: r, z, dr, dz, theta
  integer, dimension(:), allocatable :: nr_at_level, ntheta_at_level, nz_at_level, nzlocal_at_level, kguards_at_level
  integer, dimension(:), allocatable :: izspace
  integer :: ielectrode

  integer :: nr, ntheta, nz, kguards
  integer :: i, j, k
  real(kind=CUSTOM_REAL) :: z_sum
  integer :: gap_size

  ! For the multigrid method start from level=1.
  integer, parameter :: ilevel_fine = 1
  !--------------------------------------------------------------------------

  !print *, ' size, ilevel_coarse, ifixed_elecgeo, space_electrodes, space_elec_guards = '
  !print *, '     ', size, ilevel_coarse, ifixed_elecgeo, sens%space_electrodes, sens%space_elec_guards

  ! Use nz instead of nzlocal which is not needed for geometry.
  call init_params(size, 1, par, nr_read, nz_read, nzlocal_read, ntheta_read, kguards_read, ifixed_elecgeo)
  par%sens = sens

  !--------------------------------------------------------------------------
  allocate(i1(ilevel_coarse))
  allocate(i2(ilevel_coarse))
  allocate(nr_at_level(ilevel_coarse))
  allocate(ntheta_at_level(ilevel_coarse))
  allocate(nz_at_level(ilevel_coarse))
  allocate(nzlocal_at_level(ilevel_coarse))
  allocate(kguards_at_level(ilevel_coarse))
  allocate(izspace(ilevel_coarse))

  allocate(theta(ilevel_coarse))
  allocate(dr(ilevel_coarse))
  allocate(dz(ilevel_coarse))
  allocate(r(ilevel_coarse))
  allocate(z(ilevel_coarse))
  allocate(elecs(ilevel_coarse))
  allocate(guards(ilevel_coarse))

  ! Use nz instead of nzlocal which is not needed for geometry.
  call set_dimensions_multigrid(ilevel_coarse, nr_read, nz_read, nz_read, ntheta_read, kguards_read, &
                              nr_at_level, ntheta_at_level, nz_at_level, nzlocal_at_level, kguards_at_level)

  do ilevel = ilevel_fine, ilevel_coarse
    nr      = nr_at_level(ilevel)
    ntheta  = ntheta_at_level(ilevel)
    nz      = nz_at_level(ilevel)

    allocate(theta(ilevel)%p(0:ntheta+1))
    allocate(dr(ilevel)%p(0:2*nr+1))
    allocate(dz(ilevel)%p(0:nz+1))
    allocate(r(ilevel)%p(0:2*nr+1))
    allocate(z(ilevel)%p(0:nz+1))
    allocate(elecs(ilevel)%p(par%nel))
  enddo

  do ilevel = ilevel_fine, ilevel_coarse
    ! (!!!) The function being tested.
    call set_geometry(ilevel,par,nr_at_level(ilevel),ntheta_at_level(ilevel),nz_at_level(ilevel),&
                      r(ilevel)%p,dr(ilevel)%p,theta(ilevel)%p,dz(ilevel)%p,&
                      z(ilevel)%p,izspace(ilevel),kguards_at_level(ilevel),&
                      i1(ilevel),i2(ilevel),elecs(ilevel)%p,guards(ilevel))
  enddo

  ! Tests.
  do ilevel = ilevel_fine, ilevel_coarse
    nr      = nr_at_level(ilevel)
    ntheta  = ntheta_at_level(ilevel)
    nz      = nz_at_level(ilevel)
    kguards = kguards_at_level(ilevel)

    !------------------------------------------------------------------------------------------
    ! Tests of heights.
    !------------------------------------------------------------------------------------------

    ! The height is monotonically increasing.
    do k=1,nz+1
      call assert_true(z(ilevel)%p(k-1) < z(ilevel)%p(k), "z(k) value is increasing with k.")
    enddo

    ! Assume that the electrode bottom is at the kguards+1, and the electrode top is at the 3*kguards.

    ! z(k) at the electrode bottom is equal to the 1/4 cylinder height.
    call assert_comparable_real(z(ilevel)%p(kguards+1), 0.25_CUSTOM_REAL*sens%heicyl, tol, &
                                "z(k) at the electrode bottom is equal to the 1/4 cylinder height.")

    ! z(k) at the electrode top is equal to the 3/4 cylinder height.
    call assert_comparable_real(z(ilevel)%p(3*kguards), 0.75_CUSTOM_REAL*sens%heicyl, tol, &
                                "z(k) at the electrode top is equal to the 3/4 cylinder height.")

    if (izspace(ilevel) /= 0) then
      ! Correct z(k) at the top of the lower guard.
      call assert_comparable_real(z(ilevel)%p(kguards+1-izspace(ilevel)-1), &
                                  0.25_CUSTOM_REAL*sens%heicyl-sens%space_elec_guards, tol, &
                                  "Correct z(k) at the top of the lower guard.")

      ! Correct z(k) at the bottom of the upper guard.
      call assert_comparable_real(z(ilevel)%p(3*kguards+izspace(ilevel)+1), &
                                  0.75_CUSTOM_REAL*sens%heicyl+sens%space_elec_guards, tol, &
                                  "Correct z(k) at the bottom of the upper guard.")
    endif

    ! TODO: Add the test that if we have z-gaps on 1st level than we should have them on all levels???

    ! z(k) at k=0 is equal to zero.
    call assert_comparable_real(z(ilevel)%p(0), 0._CUSTOM_REAL, tol, "z(k) at k=0 is equal to zero.")

    ! z(k) at k=nz+1 is equal to the cylinder height.
    call assert_comparable_real(z(ilevel)%p(nz+1), sens%heicyl, tol, "z(k) at k=nz+1 is equal to the cylinder height.")

    ! Integral of dz is equal to the cylinder height.
    z_sum = 0._CUSTOM_REAL
    do k=1,nz+1
      z_sum = z_sum + dz(ilevel)%p(k)

      ! Correct dz(k) = z(k) - z(k-1).
      call assert_comparable_real(z(ilevel)%p(k), z(ilevel)%p(k-1)+dz(ilevel)%p(k), tol, "Correct dz(k).")
    enddo
    call assert_comparable_real(z_sum, sens%heicyl, tol, "Integral of dz is equal to the cylinder height.")

    !------------------------------------------------------------------------------------------
    ! Tests of radii.
    !------------------------------------------------------------------------------------------

    ! The radius is monotonically increasing.
    do i=1,2*nr+1
      call assert_true(r(ilevel)%p(i-1) < r(ilevel)%p(i), "r(i) value is increasing with i.")
    enddo

    ! Radii at i1, i2, and 2*nr+1 correspond to radiusin, radiusout, and radiusoutout.
    call assert_comparable_real(r(ilevel)%p(2*i1(ilevel)-1), sens%radiusin, tol, "Radius at i=i1 is radiusin.")
    call assert_comparable_real(r(ilevel)%p(2*i2(ilevel)-1), sens%radiusout, tol, "Radius at i=i2 is radiusout.")
    call assert_comparable_real(r(ilevel)%p(2*(nr+1)-1), sens%radiusoutout, tol, "Radius at i=nr+1 is radiusoutout.")

    ! Tests for a new radii Mapping #3: cell index i --> radius index 2i-1,
    ! with ghost points i=0 and i=nr+1 separated by 2dr from regular points i=1 and i=nr.
    call assert_comparable_real(2._CUSTOM_REAL*r(ilevel)%p(0), r(ilevel)%p(1), tol, &
                                "Correct radius index mapping near the axis (0 and 1).")

    call assert_comparable_real(3._CUSTOM_REAL*r(ilevel)%p(0), r(ilevel)%p(2), tol, &
                                "Correct radius index mapping near the axis (0 and 2).")

    call assert_comparable_real(2._CUSTOM_REAL*r(ilevel)%p(1), r(ilevel)%p(3), tol, &
                                "Correct radius index mapping near the axis (1 and 3).")

    call assert_comparable_real(dr(ilevel)%p(0), r(ilevel)%p(0), tol, "Correct dr(0).")

    !------------------------------------------------------------------------------------------
    ! Tests of angles.
    !------------------------------------------------------------------------------------------

    ! The angle is monotonically decreasing.
    do j=1,ntheta+1
      call assert_true(theta(ilevel)%p(j-1) > theta(ilevel)%p(j), "theta(j) value is decreasing with j.")
    enddo

    ! Correct theta values at the angle = 0 and -2Pi.
    call assert_comparable_real(theta(ilevel)%p(1), 0._CUSTOM_REAL, tol, "Correct theta(1) value.")
    call assert_comparable_real(theta(ilevel)%p(ntheta+1), -2._CUSTOM_REAL*PI, tol, "Correct theta(ntheta+1) value.")

    ! Correct periodic boundary conditions in theta.
    call assert_comparable_real(theta(ilevel)%p(0), theta(ilevel)%p(ntheta)+2._CUSTOM_REAL*PI, tol, &
                                "Correct periodic boundary conditions in theta 0.")

    call assert_comparable_real(theta(ilevel)%p(1), theta(ilevel)%p(ntheta+1)+2._CUSTOM_REAL*PI, tol, &
                                "Correct periodic boundary conditions in theta 1.")

    ! (***)
    ! Note that we can only satisfy the below tests for a system with gaps,
    ! so technically multigrid cannot be applied to systems without gaps.
    if (ilevel > 1 .and. nint(sens%space_electrodes) > 0) then
      do ielectrode=1,par%nel
        ! The min angle at the tips of the electrodes is the same at every level.
        call assert_comparable_real(theta(ilevel)%p(elecs(ilevel)%p(ielectrode)%ithetamin), &
                                    theta(ilevel-1)%p(elecs(ilevel-1)%p(ielectrode)%ithetamin), tol, &
                                    "The min angle at the tips of the electrodes is the same at every level.")

        ! The max angle at the tips of the electrodes is the same at every level.
        call assert_comparable_real(theta(ilevel)%p(elecs(ilevel)%p(ielectrode)%ithetamax), &
                                    theta(ilevel-1)%p(elecs(ilevel-1)%p(ielectrode)%ithetamax), tol, &
                                    "The max angle at the tips of the electrodes is the same at every level.")
      enddo
    endif

    ! I. The checks between electrodes, except a pair 12 and 1:

    do ielectrode=2,par%nel
      ! The same number of points at every electrode.
      call assert_equal_int(elecs(ilevel)%p(ielectrode)%ithetamax-elecs(ilevel)%p(ielectrode)%ithetamin, &
                            elecs(ilevel)%p(ielectrode-1)%ithetamax-elecs(ilevel)%p(ielectrode-1)%ithetamin, &
                            "The same number of points at every electrode.")

      ! Correct angle of gaps between electrodes.
      if (nint(sens%space_electrodes) > 0) then
        call assert_comparable_real(theta(ilevel)%p(elecs(ilevel)%p(ielectrode)%ithetamin)- &
                                    theta(ilevel)%p(elecs(ilevel)%p(ielectrode-1)%ithetamax), &
                                    -sens%space_electrodes*PI/180._CUSTOM_REAL, &
                                    tol, "Correct angle of gaps between electrodes.")
      endif

      gap_size = elecs(ilevel)%p(ielectrode)%ithetamin-elecs(ilevel)%p(ielectrode-1)%ithetamax
      if (nint(sens%space_electrodes) > 0) then
        ! The gaps between electrodes are present.
        call assert_true(gap_size > 1, "The gaps between electrodes: at least one point on the gap.")
      else
        ! The gaps between electrodes are absent.
        ! Commented this test since we actually have gaps in this case, i.e.,
        ! when electrode grid points are 1-3, 4-6, 7-9 (for grid) 36^3, there is a gap between points 3 and 4.
        !call assert_equal_int(gap_size, 0, "The gaps between electrodes are absent.")
      endif
    enddo

    ! II. The rest checks between electrode pair 12 and 1:

    if (nint(sens%space_electrodes) > 0) then
      ! Correct angle of gaps between electrodes.
      call assert_comparable_real(theta(ilevel)%p(elecs(ilevel)%p(1)%ithetamin)- &
                                  theta(ilevel)%p(elecs(ilevel)%p(12)%ithetamax), &
                                  -sens%space_electrodes*PI/180._CUSTOM_REAL+2._CUSTOM_REAL*PI, &
                                  tol, "Correct angle of gaps between electrodes (1 and 12).")
    endif

  enddo ! End Tests.

  !----------------------------------------------------------------------------------
  ! Deallocate all the arrays
  !----------------------------------------------------------------------------------

  do ilevel=ilevel_fine, ilevel_coarse
    deallocate(theta(ilevel)%p)
    deallocate(dr(ilevel)%p)
    deallocate(dz(ilevel)%p)
    deallocate(r(ilevel)%p)
    deallocate(z(ilevel)%p)
    deallocate(elecs(ilevel)%p)
  enddo

  ! deallocate all pointer arrays
  deallocate(i1)
  deallocate(i2)
  deallocate(nr_at_level)
  deallocate(ntheta_at_level)
  deallocate(nz_at_level)
  deallocate(nzlocal_at_level)
  deallocate(kguards_at_level)

  deallocate(izspace)
  deallocate(theta)
  deallocate(dr)
  deallocate(dz)
  deallocate(r)
  deallocate(z)
  deallocate(elecs)
  deallocate(guards)

end subroutine test_geometry


!==============================================================================================
! Tests of the boundary conditions and right hand side.
!==============================================================================================
subroutine test_boundary_conditions(size, ilevel_coarse, ifixed_elecgeo, sens)

  integer, intent(in) :: size
  integer, intent(in) :: ilevel_coarse, ifixed_elecgeo
  type(t_sensor), intent(in) :: sens

  type(t_parameters_ect) :: par
  integer :: ilevel
  integer :: nr_read, nz_read, nzlocal_read, ntheta_read, kguards_read

  integer, dimension(:), allocatable :: i1, i2
  type(t_electrode_p), dimension(:), allocatable :: elecs
  type(t_guards), dimension(:), allocatable :: guards
  type(real_p), dimension(:), allocatable :: r, z, dr, dz, theta
  integer, dimension(:), allocatable :: nr_at_level, ntheta_at_level, nz_at_level, kguards_at_level
  integer, dimension(:), allocatable :: nzlocal_at_level
  type(real3_p), dimension(:), allocatable :: val, phi, x2
  type(int3_p), dimension(:), allocatable :: flag
  integer, dimension(:), allocatable :: izspace
  integer :: ielectrode

  integer :: nr, ntheta, nz, nzlocal
  integer :: j, k, k1

  ! For the multigrid method start from level=1.
  integer, parameter :: ilevel_fine = 1
  !--------------------------------------------------------------------------

  !print *, ' size, ilevel_coarse, ifixed_elecgeo, space_electrodes, space_elec_guards = '
  !print *, '     ', size, ilevel_coarse, ifixed_elecgeo, sens%space_electrodes, sens%space_elec_guards

  call init_params(size, nbproc, par, nr_read, nz_read, nzlocal_read, ntheta_read, kguards_read, ifixed_elecgeo)
  par%sens = sens

  !--------------------------------------------------------------------------
  allocate(i1(ilevel_coarse))
  allocate(i2(ilevel_coarse))
  allocate(nr_at_level(ilevel_coarse))
  allocate(ntheta_at_level(ilevel_coarse))
  allocate(nz_at_level(ilevel_coarse))
  allocate(nzlocal_at_level(ilevel_coarse))
  allocate(kguards_at_level(ilevel_coarse))
  allocate(izspace(ilevel_coarse))

  allocate(theta(ilevel_coarse))
  allocate(dr(ilevel_coarse))
  allocate(dz(ilevel_coarse))
  allocate(r(ilevel_coarse))
  allocate(z(ilevel_coarse))
  allocate(elecs(ilevel_coarse))
  allocate(guards(ilevel_coarse))
  allocate(flag(ilevel_coarse))
  allocate(val(ilevel_coarse))
  allocate(phi(ilevel_coarse))
  allocate(x2(ilevel_coarse))

  call set_dimensions_multigrid(ilevel_coarse, nr_read, nz_read, nzlocal_read, ntheta_read, kguards_read, &
                    nr_at_level, ntheta_at_level, nz_at_level, nzlocal_at_level, kguards_at_level)

  do ilevel = ilevel_fine, ilevel_coarse
    nr      = nr_at_level(ilevel)
    ntheta  = ntheta_at_level(ilevel)
    nz      = nz_at_level(ilevel)
    nzlocal = nzlocal_at_level(ilevel)

    allocate(theta(ilevel)%p(0:ntheta+1))
    allocate(dr(ilevel)%p(0:2*nr+1))
    allocate(dz(ilevel)%p(0:nz+1))
    allocate(r(ilevel)%p(0:2*nr+1))
    allocate(z(ilevel)%p(0:nz+1))
    allocate(elecs(ilevel)%p(par%nel))
    allocate(flag(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(val(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(phi(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(x2(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
  enddo

  ! The tests depend on the geometry() function.
  do ilevel = ilevel_fine, ilevel_coarse
    call set_geometry(ilevel,par,nr_at_level(ilevel),ntheta_at_level(ilevel),nz_at_level(ilevel),&
                      r(ilevel)%p,dr(ilevel)%p,theta(ilevel)%p,dz(ilevel)%p,&
                      z(ilevel)%p,izspace(ilevel),kguards_at_level(ilevel),&
                      i1(ilevel),i2(ilevel),elecs(ilevel)%p,guards(ilevel))
  enddo

  !----------------------------------------------------------------------------------
  ! TESTS of flag_init_outer_bc() function: depends on geometry()
  !----------------------------------------------------------------------------------

  do ilevel = ilevel_fine, ilevel_coarse
    call flag_init_bc(flag(ilevel)%p,nr_at_level(ilevel),ntheta_at_level(ilevel), &
                      nzlocal_at_level(ilevel),i1(ilevel),i2(ilevel),par%nel, &
                      elecs(ilevel)%p,guards(ilevel),myrank,nbproc)
  enddo

  ! Tests.
  do ilevel = ilevel_fine, ilevel_coarse
    nr      = nr_at_level(ilevel)
    ntheta  = ntheta_at_level(ilevel)
    nzlocal = nzlocal_at_level(ilevel)

    ! Check the BC's at electrodes.
    do k=0,nzlocal+1
      k1=myrank*nzlocal+k

      do ielectrode=1,par%nel
        do j=elecs(ilevel)%p(ielectrode)%ithetamin,elecs(ilevel)%p(ielectrode)%ithetamax
          ! Skip the guards.
          if (k1 >= elecs(ilevel)%p(ielectrode)%izmin .and. k1 <= elecs(ilevel)%p(ielectrode)%izmax) then
            call assert_equal_int(flag(ilevel)%p(i2(ilevel),j,k), 1, "Correct BC's flag(i,j,k) at electrodes.")
          endif
        enddo
      enddo
    enddo

  enddo ! End Tests.

  !----------------------------------------------------------------------------------
  ! TESTS of solution_init() function: depends on geometry() and flag_init_outer_bc()
  !----------------------------------------------------------------------------------

  do ielectrode=1,par%nel
    nr      = nr_at_level(ilevel_fine)
    ntheta  = ntheta_at_level(ilevel_fine)
    nzlocal = nzlocal_at_level(ilevel_fine)

    ! (!!!) The function being tested.
    ! It is called only on the fine level!
    call solution_init(phi(ilevel_fine)%p,val(ilevel_fine)%p,flag(ilevel_fine)%p,&
                       x2(ilevel_fine)%p,par%dims,&
                       myrank,i2(ilevel_fine),elecs(ilevel_fine)%p(ielectrode))

    ! The solution is set to zero everywhere at the outermost cylinder i=nr+1.
    do k=0,nzlocal+1
      do j=0,ntheta+1
        call assert_comparable_real(val(ilevel_fine)%p(nr+1,j,k), 0._CUSTOM_REAL, tol, &
                                    "The solution is set to zero everywhere at the outermost cylinder.")
      enddo
    enddo

    ! The solution is set to unity everywhere at the electrode.
    do k=0,nzlocal+1
      k1=myrank*nzlocal+k

      do j=elecs(ilevel_fine)%p(ielectrode)%ithetamin,elecs(ilevel_fine)%p(ielectrode)%ithetamax
        ! Skip the guards.
        if (k1 >= elecs(ilevel_fine)%p(ielectrode)%izmin .and. k1 <= elecs(ilevel_fine)%p(ielectrode)%izmax) then
          call assert_comparable_real(val(ilevel_fine)%p(i2(ilevel_fine),j,k), 1._CUSTOM_REAL, tol, &
                                      "The solution is set to unity everywhere at the electrode.")
        endif
      enddo
    enddo

    ! The solution is set to zero everywhere at the guards.
    do k=0,nzlocal+1
      k1=myrank*nzlocal+k

      if (k1 <= guards(ilevel_fine)%lower_izmax .and. k1 >= guards(ilevel_fine)%upper_izmin) then
        do j=0,ntheta+1
          call assert_comparable_real(val(ilevel_fine)%p(i2(ilevel_fine),j,k), 0._CUSTOM_REAL, tol, &
                                      "The solution is set to zero everywhere at the guards.")
        enddo
      endif
    enddo

  enddo


  !----------------------------------------------------------------------------------
  ! Deallocate all the arrays
  !----------------------------------------------------------------------------------

  do ilevel=ilevel_fine, ilevel_coarse
    deallocate(theta(ilevel)%p)
    deallocate(dr(ilevel)%p)
    deallocate(dz(ilevel)%p)
    deallocate(r(ilevel)%p)
    deallocate(z(ilevel)%p)
    deallocate(elecs(ilevel)%p)
    deallocate(flag(ilevel)%p)
    deallocate(val(ilevel)%p)
    deallocate(phi(ilevel)%p)
    deallocate(x2(ilevel)%p)
  enddo

  ! deallocate all pointer arrays
  deallocate(i1)
  deallocate(i2)
  deallocate(nr_at_level)
  deallocate(ntheta_at_level)
  deallocate(nz_at_level)
  deallocate(nzlocal_at_level)
  deallocate(kguards_at_level)

  deallocate(izspace)
  deallocate(theta)
  deallocate(dr)
  deallocate(dz)
  deallocate(r)
  deallocate(z)
  deallocate(elecs)
  deallocate(guards)
  deallocate(flag)
  deallocate(val)
  deallocate(phi)
  deallocate(x2)

end subroutine test_boundary_conditions

!==============================================================================================
! Compare the solution to analytical ones for the Laplace equation in cylinder with BC's:
! (1) u=1 on the side wall, and u=0 on the top and bottom walls [M. A. Pinsky, Example 3.5.1],
! (2) u=1 on the top wall, and u=0 on the bottom and side walls [M. A. Pinsky, Example 3.5.2].
!==============================================================================================
subroutine test_analytical_comparison(bc_type, size, ilevel_coarse, itmax, sens)

  ! The type boundary conditions -- to test different analytical solutions.
  integer, intent(in) :: bc_type

  integer, intent(in) :: size
  integer, intent(in) :: ilevel_coarse, itmax
  type(t_sensor), intent(in) :: sens

  type(t_parameters_ect) :: par
  type(t_dimensions), allocatable :: dims_at_level(:)
  integer :: ilevel
  integer :: nr_read, nz_read, nzlocal_read, ntheta_read, kguards_read

  integer, dimension(:), allocatable :: i1, i2
  type(t_electrode_p), dimension(:), allocatable :: elecs
  type(t_guards), dimension(:), allocatable :: guards
  type(real_p), dimension(:), allocatable :: r, z, dr, dz, theta
  integer, dimension(:), allocatable :: nr_at_level, ntheta_at_level, nz_at_level, kguards_at_level
  integer, dimension(:), allocatable :: nzlocal_at_level
  type(real3_p), dimension(:), allocatable :: xgrid, ygrid, zgrid
  type(real3_p), dimension(:), allocatable :: val, phi, x2, permit
  type(int3_p), dimension(:), allocatable :: flag
  integer, dimension(:), allocatable :: izspace
  real(kind=CUSTOM_REAL) permit_isolated_tube, permit_air
  type(real4_p), dimension(:), allocatable :: a
  type(real3_p), dimension(:), allocatable :: b_RHS
  integer :: iter_pcg, ifixed_elecgeo
  type (t_pcg_auxarrays) :: pcg_auxarrays
  ! Analytical solution.
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: asol

  integer :: nr, ntheta, nz, nzlocal
  integer :: i, j, k, imax, kmax
  integer :: ierr
  !character(len=40) :: name

  ! For the multigrid method start from level=1.
  integer, parameter :: ilevel_fine = 1

  ! Frequency of printing information during iterations.
  integer, parameter :: output_frequency = 20
  !--------------------------------------------------------------------------

  !print *, ' bc_type, size, ilevel_coarse, itmax = '
  !print *, '     ', bc_type, size, ilevel_coarse, itmax

  ! No electrodes in analytical solutions.
  ifixed_elecgeo = 0

  ! RM RM added this 2 March 2015
  par%iprecond = 1
  par%omega1 = 0.8_CUSTOM_REAL

  call init_params(size, nbproc, par, nr_read, nz_read, nzlocal_read, ntheta_read, kguards_read, ifixed_elecgeo)
  par%sens = sens

  !--------------------------------------------------------------------------
  allocate(i1(ilevel_coarse))
  allocate(i2(ilevel_coarse))
  allocate(nr_at_level(ilevel_coarse))
  allocate(ntheta_at_level(ilevel_coarse))
  allocate(nz_at_level(ilevel_coarse))
  allocate(nzlocal_at_level(ilevel_coarse))
  allocate(dims_at_level(ilevel_coarse))
  allocate(kguards_at_level(ilevel_coarse))
  allocate(izspace(ilevel_coarse))

  allocate(theta(ilevel_coarse))
  allocate(dr(ilevel_coarse))
  allocate(dz(ilevel_coarse))
  allocate(r(ilevel_coarse))
  allocate(z(ilevel_coarse))
  allocate(xgrid(ilevel_coarse))
  allocate(ygrid(ilevel_coarse))
  allocate(zgrid(ilevel_coarse))
  allocate(elecs(ilevel_coarse))
  allocate(guards(ilevel_coarse))
  allocate(flag(ilevel_coarse))
  allocate(val(ilevel_coarse))
  allocate(phi(ilevel_coarse))
  allocate(x2(ilevel_coarse))
  allocate(permit(ilevel_coarse))
  allocate(a(ilevel_coarse))
  allocate(b_RHS(ilevel_coarse))

  call set_dimensions_multigrid(ilevel_coarse, nr_read, nz_read, nzlocal_read, ntheta_read, kguards_read, &
                                nr_at_level, ntheta_at_level, nz_at_level, nzlocal_at_level, kguards_at_level)

  do ilevel = ilevel_fine, ilevel_coarse
    nr      = nr_at_level(ilevel)
    ntheta  = ntheta_at_level(ilevel)
    nz      = nz_at_level(ilevel)
    nzlocal = nzlocal_at_level(ilevel)

    ! TODO: keep only dims_at_level (remove nr_at_level, ...).
    dims_at_level(ilevel)%nr = nr
    dims_at_level(ilevel)%ntheta = ntheta
    dims_at_level(ilevel)%nz = nz
    dims_at_level(ilevel)%nzlocal = nzlocal

    allocate(theta(ilevel)%p(0:ntheta+1))
    allocate(dr(ilevel)%p(0:2*nr+1))
    allocate(dz(ilevel)%p(0:nz+1))
    allocate(r(ilevel)%p(0:2*nr+1))
    allocate(z(ilevel)%p(0:nz+1))
    allocate(xgrid(ilevel)%p(0:nr+1,0:ntheta+1,0:nz+1))
    allocate(ygrid(ilevel)%p(0:nr+1,0:ntheta+1,0:nz+1))
    allocate(zgrid(ilevel)%p(0:nr+1,0:ntheta+1,0:nz+1))
    allocate(elecs(ilevel)%p(par%nel))
    allocate(flag(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(val(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(phi(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(x2(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(permit(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(a(ilevel)%p(1:nr,1:ntheta,1:nzlocal,7))
    allocate(b_RHS(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
  enddo

  allocate(pcg_auxarrays%p(0:nr_read+1,0:ntheta_read+1,0:nzlocal_read+1))
  allocate(pcg_auxarrays%r(0:nr_read+1,0:ntheta_read+1,0:nzlocal_read+1))
  allocate(pcg_auxarrays%inv_diagA(1:nr_read,1:ntheta_read,1:nzlocal_read))
  allocate(pcg_auxarrays%z(0:nr_read+1,0:ntheta_read+1,0:nzlocal_read+1))
  allocate(pcg_auxarrays%tmp(0:nr_read+1,0:ntheta_read+1,0:nzlocal_read+1))

  allocate(asol(0:nr_read+1,0:ntheta_read+1,0:nz_read+1))
  !----------------------------------------------------------------------------------

  do ilevel = ilevel_fine, ilevel_coarse
    call set_geometry(ilevel,par,nr_at_level(ilevel),ntheta_at_level(ilevel),nz_at_level(ilevel),&
                      r(ilevel)%p,dr(ilevel)%p,theta(ilevel)%p,dz(ilevel)%p,&
                      z(ilevel)%p,izspace(ilevel),kguards_at_level(ilevel),&
                      i1(ilevel),i2(ilevel),elecs(ilevel)%p,guards(ilevel))

    ! TODO: initialize xgrid,ygrid,zgrid
  enddo

  permit_air = 1._CUSTOM_REAL
  permit_isolated_tube = 1._CUSTOM_REAL

  do ilevel=ilevel_fine, ilevel_coarse
    nr      = nr_at_level(ilevel)
    ntheta  = ntheta_at_level(ilevel)
    nz      = nz_at_level(ilevel)

    ! Set permittivity everywhere to constant.
    permit(ilevel)%p(:,:,:) = permit_air

    ! The position of the side wall.
    imax = nr+1

    ! The position of the top wall.
    kmax = nz+1

    ! Set Dirichlet BCs flag.
    flag(ilevel)%p(:,:,:) = 0

    do k=0,nz+1
      do j=0,ntheta+1
        do i=0,nr+1
          ! Side wall.
          if (i == imax) flag(ilevel)%p(i,j,k) = 1
          ! Bottom wall.
          if (k == 0) flag(ilevel)%p(i,j,k) = 1
          ! Top wall.
          if (k == kmax) flag(ilevel)%p(i,j,k) = 1
        enddo
      enddo
    enddo
  enddo

  ! (+) Set boundary conditions (through the right-hand-side of Ax=b_RHS) --------
  b_RHS(ilevel_fine)%p(:,:,:) = 0._CUSTOM_REAL

  if (bc_type == 1) then
    do k=1,kmax-1
      do j=0,ntheta_read+1
        ! Side wall excluding first and last points (bottom and top walls).
        b_RHS(ilevel_fine)%p(imax,j,k) = 1._CUSTOM_REAL
      enddo
    enddo
  else if (bc_type == 2) then
    do j=0,ntheta_read+1
      do i=0,imax-1
        ! Top wall excluding last point (side wall).
        b_RHS(ilevel_fine)%p(i,j,kmax) = 1._CUSTOM_REAL
      enddo
    enddo
  else if (bc_type == 3) then
    do k=(kmax-1)/4+1,3*(kmax-1)/4
      do j=0,ntheta_read+1
        ! Side wall excluding guards.
        b_RHS(ilevel_fine)%p(imax,j,k) = 1._CUSTOM_REAL
      enddo
    enddo
  else if (bc_type == 4) then
    do k=1,kmax-1
      do j=0,ntheta_read+1
        ! Gradually side walls.
        ! f(x)=4x(1-x), f(0)=f(1)=0, f(1/2)=1, f'(1/2)=0.
        b_RHS(ilevel_fine)%p(imax,j,k) = &
                  4._CUSTOM_REAL*dble(k)/dble(kmax)*(1._CUSTOM_REAL - dble(k)/dble(kmax))
      enddo
    enddo
  else
    print *, "Unknown boundary conditions type!"
    stop
  endif
  ! (-) Set boundary conditions ------------------------------------------------------------

  ! Initialize matrix.
  do ilevel = ilevel_fine, ilevel_coarse
    call matrix_init(a(ilevel)%p, flag(ilevel)%p, permit(ilevel)%p, permit_isolated_tube, permit_air, r(ilevel)%p, &
                     theta(ilevel)%p, dz(ilevel)%p, dims_at_level(ilevel), i1(ilevel), i2(ilevel), myrank)
  enddo

  ! Smoother.
  do k=1,nzlocal_read
    do j=1,ntheta_read
      do i=1,nr_read
        pcg_auxarrays%inv_diagA(i,j,k)=1._CUSTOM_REAL/a(ilevel_fine)%p(i,j,k,3)
      enddo
    enddo
  enddo

  ! Initial guess.
  x2(ilevel_fine)%p(:,:,:) = b_RHS(ilevel_fine)%p(:,:,:)

  ! Solve the system numerically.
  call solver_pcg(a(ilevel_fine)%p,b_RHS(ilevel_fine)%p,x2(ilevel_fine)%p,NORM_L2,par%iprecond,par%omega1,&
                  1.e-12_CUSTOM_REAL,itmax,iter_pcg,output_frequency,.false.,&
                  nr_at_level(ilevel_fine),ntheta_at_level(ilevel_fine),nzlocal_at_level(ilevel_fine),&
                  ierr,myrank,nbproc,pcg_auxarrays)

  ! Enforce periodic boundary conditions.
  call enforce_pb(nr_read, ntheta_read, nzlocal_read, x2(ilevel_fine)%p)

  ! Find analytical solution.
  if (bc_type == 1) then
    call solution1(nr_read,ntheta_read,nz_read,r(ilevel_fine)%p,z(ilevel_fine)%p,asol,imax,kmax,1)
  else if (bc_type == 2) then
    call solution2(nr_read,ntheta_read,nz_read,r(ilevel_fine)%p,z(ilevel_fine)%p,asol,imax,kmax)
  else if (bc_type == 3) then
    call solution1(nr_read,ntheta_read,nz_read,r(ilevel_fine)%p,z(ilevel_fine)%p,asol,imax,kmax,2)
  else if (bc_type == 4) then
    call solution1(nr_read,ntheta_read,nz_read,r(ilevel_fine)%p,z(ilevel_fine)%p,asol,imax,kmax,3)
  endif

  !------------------------------------------------------------------------------------------------
  nr      = nr_at_level(ilevel_fine)
  ntheta  = ntheta_at_level(ilevel_fine)
  nz      = nz_at_level(ilevel_fine)

!  path_output = "output_tests/"
!
!  write (name, '("phi3d_half_nz",i3.3,"_numerical.vtk")') nz
!  call visualisation_paraview(name, myrank, nr, ntheta, nz, &
!        x2(ilevel_fine)%p, xgrid(ilevel_fine)%p, ygrid(ilevel_fine)%p, zgrid(ilevel_fine)%p, &
!        0, imax, ntheta/2+1, ntheta+1, 0, kmax, 1, ntheta/4, 1, 'POINT_DATA')
!
!  write (name, '("phi3d_half_nz",i3.3,"_analytical.vtk")') nz
!  call visualisation_paraview(name, myrank, nr, ntheta, nz, &
!        asol, xgrid(ilevel_fine)%p, ygrid(ilevel_fine)%p, zgrid(ilevel_fine)%p, &
!        0, imax, ntheta/2+1, ntheta+1, 0, kmax, 1, ntheta/4, 1, 'POINT_DATA')

  ! Calculate the error.
  do k=0,kmax
    do j=0,ntheta_read+1
      do i=0,imax
        if (i == 0) then
          asol(i,j,k) = 0._CUSTOM_REAL
        endif

        ! !! TEMPORARY
!        if (flag(ilevel_fine)%p(i,j,k) == 1) then
!          asol(i,j,k) = 0._CUSTOM_REAL
!        endif

        ! The relative error.
        if (asol(i,j,k) /= 0._CUSTOM_REAL) then
          asol(i,j,k) = abs((asol(i,j,k) - x2(ilevel_fine)%p(i,j,k)) / asol(i,j,k))

          ! Assert the error is small.
          call assert_true(asol(i,j,k) < 5.d-2, "Large error in test_analytical_comparison.")
        endif
      enddo
    enddo
  enddo

!  write (name, '("phi3d_half_nz",i3.3,"_error.vtk")') nz
!  call visualisation_paraview(name, myrank, nr, ntheta, nz, &
!        asol, xgrid(ilevel_fine)%p, ygrid(ilevel_fine)%p, zgrid(ilevel_fine)%p, &
!        0, imax, ntheta/2+1, ntheta+1, 0, kmax, 1, ntheta/4, 1, 'POINT_DATA')

  !----------------------------------------------------------------------------------
  ! Deallocate all the arrays
  !----------------------------------------------------------------------------------

  do ilevel=ilevel_fine, ilevel_coarse
    deallocate(theta(ilevel)%p)
    deallocate(dr(ilevel)%p)
    deallocate(dz(ilevel)%p)
    deallocate(r(ilevel)%p)
    deallocate(z(ilevel)%p)
    deallocate(xgrid(ilevel)%p)
    deallocate(ygrid(ilevel)%p)
    deallocate(zgrid(ilevel)%p)
    deallocate(elecs(ilevel)%p)
    deallocate(flag(ilevel)%p)
    deallocate(val(ilevel)%p)
    deallocate(phi(ilevel)%p)
    deallocate(x2(ilevel)%p)
    deallocate(permit(ilevel)%p)
    deallocate(a(ilevel)%p)
    deallocate(b_RHS(ilevel)%p)
  enddo

  ! deallocate all pointer arrays
  deallocate(i1)
  deallocate(i2)
  deallocate(nr_at_level)
  deallocate(ntheta_at_level)
  deallocate(nz_at_level)
  deallocate(nzlocal_at_level)
  deallocate(kguards_at_level)

  deallocate(izspace)
  deallocate(theta)
  deallocate(dr)
  deallocate(dz)
  deallocate(r)
  deallocate(z)
  deallocate(xgrid)
  deallocate(ygrid)
  deallocate(zgrid)
  deallocate(elecs)
  deallocate(guards)
  deallocate(flag)
  deallocate(val)
  deallocate(phi)
  deallocate(x2)
  deallocate(permit)
  deallocate(a)
  deallocate(b_RHS)

  deallocate(asol)

end subroutine test_analytical_comparison

end module tests_ect

