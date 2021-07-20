
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

module solution

  use global_typedefs, only: CUSTOM_REAL
  use utils, only: enforce_pb
  use parameters_ect, only: t_dimensions, t_electrode
  use paraview_ect

  implicit none

  private

  public :: solution_init
  public :: visualize_solution

contains

!===========================================================================================
! Initialize initial solution x2 and phi, and val used to set the right hand side.
!===========================================================================================
subroutine solution_init(phi, val, flag, x2, dims, myrank, i2, elec)

  type(t_dimensions), intent(in) :: dims
  integer, intent(in) :: myrank
  integer, intent(in) :: i2
  type(t_electrode), intent(in) :: elec
  integer, intent(in) :: flag(0:, 0:, 0:)

  real(kind=CUSTOM_REAL), intent(out) :: phi(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(out) :: val(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(out) :: x2(0:, 0:, 0:)

  integer i, j, k, k1

  ! Phase 1: fill val array (needed to assemble the RHS in source_RHS().
  ! Solution phi is imposed to val=1 for ielectrode and 0 for all other electrodes.

  ! Initialize val to zero.
  val = 0._CUSTOM_REAL

  do k=0,dims%nzlocal+1
    k1 = myrank*dims%nzlocal+k
    ! Skip the guards.
    if (k1 >= elec%izmin .and. k1 <= elec%izmax) then
      do j = elec%ithetamin, elec%ithetamax
        ! VO VO Here flag(i,j,k) must be unity by construction.

        ! DG DG this is technically speaking not correct for finite volumes: We have unknowns in
        ! DG DG volume centers, thus no solution values on the boundary (which are on volume edges.
        ! DG DG no idea why it works nontheless
        val(i2,j,k) = 1._CUSTOM_REAL
      enddo

      ! VO VO Fix the electrode #12 last point for the case without gaps between electrode.
      ! VO VO Because the point j=1 is the same as j=ntheta+1 (which is the last point of the electrode),
      ! VO VO it must have phi=1.
      ! VO VO This condition should be true only for the case electrode=12 and no gaps between the electrodes.
      if (elec%ithetamax == dims%ntheta+1) val(i2,1,k) = 1._CUSTOM_REAL
    endif
  enddo

  ! Periodic conditions of potentials.
  call enforce_pb(dims%nr, dims%ntheta, dims%nzlocal, val)

  ! Zero out phi.
  phi = 0._CUSTOM_REAL

!  ! Read initial guess from file if requested (this overwrites phi, except at the ghost cells).
!  if (read_guess_from_file) then
!    open(unit=991,file='gradient_output_true_exact_001',status='old',action='read',iostat=ierr)
!    if (ierr /= 0) call exit_MPI("something went wrong reading initial guess",myrank,19)
!    do k=1,nzlocal
!      do j=1,ntheta
!        do i=1,nr
!          read(991,*) phi(i,j,k)
!        enddo
!      enddo
!    enddo
!    close(991)
!  endif

  ! Initial guess for the linear solver is also just a copy of phi, with known Dirichlet values
  ! directly inserted, so that Dirichlet values are no longer treated as unknowns by the linear solver.
  do k=0,dims%nzlocal+1
    do j=0,dims%ntheta+1
      do i=0,dims%nr+1
        x2(i,j,k) = phi(i,j,k)
        ! Dirichlet values.
        if (flag(i,j,k) == 1) x2(i,j,k) = val(i,j,k)
      enddo
    enddo
  enddo

end subroutine solution_init

!===========================================================================================
! Create Paraview files for the potential phi visualization.
!===========================================================================================
subroutine visualize_solution(dims, phi, xgrid, ygrid, zgrid, ielectrode, myrank)
  type(t_dimensions), intent(in) :: dims
  real(kind=CUSTOM_REAL), intent(in) :: phi(0:, 0:, 0:)
  ! Coordinates of points in the grid.
  real(kind=CUSTOM_REAL), intent(in) :: xgrid(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: ygrid(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: zgrid(0:, 0:, 0:)
  integer, intent(in) :: ielectrode, myrank

  call paraview_write_2d_profiles(myrank, "phi2d_", ielectrode, &
                                  dims%nr,dims%ntheta,dims%nz,dims%nzlocal, &
                                  phi, xgrid, ygrid, zgrid)

  call paraview_write_3d_profiles(myrank, "phi3d_", ielectrode, &
                                  dims%nr,dims%ntheta,dims%nzlocal, &
                                  phi, xgrid, ygrid, zgrid)

end subroutine visualize_solution

end module solution
