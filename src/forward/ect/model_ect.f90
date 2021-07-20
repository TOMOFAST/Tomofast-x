
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

module model_ect

  use global_typedefs
  use parameters_ect
  use utils
  use geometry, only: compute_shifted_grid_coordinates, compute_grid_coordinates
  use paraview_ect
  use mpi_tools

  implicit none

  private

  public :: model_init
  public :: write_model_to_1D_array
  public :: init_model_from_inversion
  public :: visualize_model

contains

!====================================================================================================
! This routine sets the permittivities in the three parts of the cylinder, and also inside bubbles.
! i2: location of receivers and senders corresponding to radius=radiusout (R2).
! i1: location of the frontier between tube and inner core corresponding to radius=radiusin (R1).
!====================================================================================================
subroutine model_init(par, permit, dims, i1, i2, permit_bubble, xb, yb, zb, rb, r, theta, z, myrank)

  type(t_parameters_ect), intent(in) :: par
  type(t_dimensions), intent(in) :: dims
  integer, intent(in) :: i1,i2
  real(kind=CUSTOM_REAL), intent(in) :: permit_bubble(:)
  real(kind=CUSTOM_REAL), intent(in) :: xb(:)
  real(kind=CUSTOM_REAL), intent(in) :: yb(:)
  real(kind=CUSTOM_REAL), intent(in) :: zb(:)
  real(kind=CUSTOM_REAL), intent(in) :: rb(:)
  real(kind=CUSTOM_REAL), intent(in) :: r(0:)
  real(kind=CUSTOM_REAL), intent(in) :: theta(0:)
  real(kind=CUSTOM_REAL), intent(in) :: z(0:)
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: permit(0:, 0:, 0:)

  integer :: ib, i, j, k, k1, ierr

  real(kind=CUSTOM_REAL), allocatable :: xgrid(:, :, :)
  real(kind=CUSTOM_REAL), allocatable :: ygrid(:, :, :)
  real(kind=CUSTOM_REAL), allocatable :: zgrid(:, :, :)

  allocate(xgrid(0:dims%nr+1, 0:dims%ntheta+1, 0:dims%nz+1), source=0._CUSTOM_REAL, stat=ierr)
  allocate(ygrid(0:dims%nr+1, 0:dims%ntheta+1, 0:dims%nz+1), source=0._CUSTOM_REAL, stat=ierr)
  allocate(zgrid(0:dims%nr+1, 0:dims%ntheta+1, 0:dims%nz+1), source=0._CUSTOM_REAL, stat=ierr)

  ! VO VO 16.07.2015: Fixed a bug. Grid for setting the permittivities is shifted from the potential phi grid.
  call compute_shifted_grid_coordinates(dims%nr, dims%ntheta, dims%nz, r, theta, z, xgrid, ygrid, zgrid)

  if (par%num_bubbles <= 0) then
    ! If we do not have bubbles, set permit_matrix everywhere up to radiusin, i.e., up to i=i1.
    do k=0,dims%nzlocal+1
      do j=0,dims%ntheta+1
        do i=0,i1
          permit(i,j,k) = par%permit_matrix
        enddo
      enddo
    enddo
  else
    ! If we have bubbles, they get a different permittivity.
    do k=0,dims%nzlocal+1
      k1=myrank*dims%nzlocal+k
      do j=0,dims%ntheta+1
        do i=0,i1
          permit(i,j,k) = par%permit_oil

          do ib=1,par%num_bubbles
            ! Check if we are inside a bubble.
            if ((xgrid(i,j,k1)-xb(ib))**2 + (ygrid(i,j,k1)-yb(ib))**2 + (zgrid(i,j,k1)-zb(ib))**2 <= rb(ib)**2) &
              permit(i,j,k) = permit_bubble(ib)
          enddo
        enddo
      enddo
    enddo
  endif

  ! Liquid (often oil) between radii radiusin and radiusout (i.e. the area of the tube of the real instrument under study).
  do k=0,dims%nzlocal+1
    do j=0,dims%ntheta+1
      do i=i1+1,i2
        permit(i,j,k) = par%permit_isolated_tube
      enddo
    enddo
  enddo

  ! Air between radiusout and radiusoutout (between electrodes and screen).
  do k=0,dims%nzlocal+1
    do j=0,dims%ntheta+1
      do i=i2+1,dims%nr+1
        permit(i,j,k) = par%permit_air
      enddo
    enddo
  enddo

  ! Apply periodic boundary conditions.
  call enforce_pb(dims%nr, dims%ntheta, dims%nzlocal, permit)

  deallocate(zgrid)
  deallocate(ygrid)
  deallocate(xgrid)

end subroutine model_init

!====================================================================================================
! Write the (prior) model for the inversion.
!====================================================================================================
subroutine write_model_to_1D_array(permit, dims, model, myrank, nbproc)
  type(t_dimensions), intent(in) :: dims
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(in) :: permit(0:, 0:, 0:)

  real(kind=CUSTOM_REAL), intent(out) :: model(:)

  integer :: i, j, k, l

  ! NOTE: Use same loop as in calculate_sensitivity(), to match the dimensions!

  l = 0
  do k = 1, dims%nzlocal + (myrank + 1) / nbproc
    do j = 1, dims%ntheta
      do i = 2, dims%nr + 1
        l = l + 1
        model(l) = permit(i ,j, k)
      enddo
    enddo
  enddo

end subroutine write_model_to_1D_array

!====================================================================================================
! Initialize the model from inversion:
! 1) Convert from 1D array (model) to 3D array (permit)
! 2) Initialize i=0, i=1 values.
! 3) Initialize k1=0, k1=nz+1 values
! 4) Initialize j=0, j=ntheta+1 values (periodic BC's).
! 5) MPI send/receive the k=0, k=nzlocal+1 values.
!====================================================================================================
subroutine init_model_from_inversion(par,dims,permit,model,i1,i2,myrank,nbproc)

  type(t_parameters_ect), intent(in) :: par
  type(t_dimensions), intent(in) :: dims
  integer, intent(in) :: myrank,nbproc
  integer, intent(in) :: i1,i2
  ! Permittivity model from inversion.
  real(kind=CUSTOM_REAL), intent(in) :: model(:)

  real(kind=CUSTOM_REAL), intent(out) :: permit(0:,0:,0:)

  integer :: i,j,k,l,ierr
  real(kind=CUSTOM_REAL) :: permit_avg

  ! NOTE: Use same loop as in calculate_sensitivity(), to match the dimensions!

  l = 0
  do k=1,dims%nzlocal + (myrank+1)/nbproc
    do j=1,dims%ntheta
      do i=2,dims%nr+1
        l = l+1
        permit(i,j,k) = model(l)

        ! NOTE, manually change the permittivity!
        ! VO VO commented: this leads to a bad misfit behavior - it increases after 1st iteration,
        ! VO VO: but the images look better, so maybe a problem in the forward problem.
        !if (permit(i,j,k) > par%permit_matrix) permit(i,j,k) = par%permit_matrix

        ! Lower bound modification makes sense physically.
        ! At least we must not have negative permittivity!
        ! And for the 1-bubble case it goes to negative values in the center.
        if (permit(i,j,k) < par%permit0) permit(i,j,k) = par%permit0
      enddo
    enddo
  enddo

  ! Initialize i=0 & i=1 values using average permittivity at i=2 ring.
  do k=1,dims%nzlocal + (myrank+1)/nbproc
    permit_avg = 0._CUSTOM_REAL
    do j=1,dims%ntheta
      permit_avg = permit_avg + permit(2,j,k)
    enddo
    permit_avg = permit_avg / dble(dims%ntheta)

    permit(0,:,k) = permit_avg
    permit(1,:,k) = permit_avg
  enddo

  if (myrank == 0) then
    ! Initialize k=0 values.
    permit(0:i1,:,0) = par%permit_matrix
    permit(i1+1:i2,:,0) = par%permit_isolated_tube
    permit(i2+1:dims%nr+1,:,0) = par%permit_air
  endif

  if (myrank == nbproc-1) then
    ! Initialize k=nz+1 values.
    permit(0:i1,:,dims%nzlocal+1) = par%permit_matrix
    permit(i1+1:i2,:,dims%nzlocal+1) = par%permit_isolated_tube
    permit(i2+1:dims%nr+1,:,dims%nzlocal+1) = par%permit_air
  endif

  ! Initialize known sensor permittivity values between R1 and R3.
  permit(i1+1:i2,:,:) = par%permit_isolated_tube
  permit(i2+1:dims%nr+1,:,:) = par%permit_air

  ! Initialize j=0, j=ntheta+1 values (apply periodic boundary conditions).
  call enforce_pb(dims%nr, dims%ntheta, dims%nzlocal, permit)

  ! Initialize k=0, k=nzlocal+1 values.
  if (nbproc > 1) call mpisendrecv(permit, dims%nr, dims%ntheta, dims%nzlocal, myrank, nbproc, ierr)

end subroutine init_model_from_inversion

!====================================================================================================
! Create Paraview files for permittivity model visualization.
! For this, first create the shifted grid where the permittivity is defined.
!====================================================================================================
subroutine visualize_model(permit, r, theta, z, dims, ielectrode, myrank)
  integer, intent(in) :: ielectrode, myrank
  type(t_dimensions), intent(in) :: dims

  real(kind=CUSTOM_REAL), intent(in) :: permit(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: r(0:)
  real(kind=CUSTOM_REAL), intent(in) :: theta(0:)
  real(kind=CUSTOM_REAL), intent(in) :: z(0:)

  integer :: ierr

  real(kind=CUSTOM_REAL), allocatable :: xgrid(:, :, :)
  real(kind=CUSTOM_REAL), allocatable :: ygrid(:, :, :)
  real(kind=CUSTOM_REAL), allocatable :: zgrid(:, :, :)

  allocate(xgrid(0:dims%nr+1, 0:dims%ntheta+1, 0:dims%nz+1), source=0._CUSTOM_REAL, stat=ierr)
  allocate(ygrid(0:dims%nr+1, 0:dims%ntheta+1, 0:dims%nz+1), source=0._CUSTOM_REAL, stat=ierr)
  allocate(zgrid(0:dims%nr+1, 0:dims%ntheta+1, 0:dims%nz+1), source=0._CUSTOM_REAL, stat=ierr)

  call compute_shifted_grid_coordinates(dims%nr, dims%ntheta, dims%nz, r, theta, z, xgrid, ygrid, zgrid)
  !call compute_grid_coordinates(dims%nr, dims%ntheta, dims%nz, r, theta, z, xgrid, ygrid, zgrid)

  call paraview_write_2d_profiles(myrank, "model2d_", ielectrode, &
                                  dims%nr, dims%ntheta, dims%nz, dims%nzlocal, &
                                  permit, xgrid, ygrid, zgrid)

  call paraview_write_3d_profiles(myrank, "model3d_", ielectrode, &
                                  dims%nr, dims%ntheta, dims%nzlocal, &
                                  permit, xgrid, ygrid, zgrid)

  deallocate(zgrid)
  deallocate(ygrid)
  deallocate(xgrid)

end subroutine visualize_model

end module model_ect
