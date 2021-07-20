
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

module paraview_ect

  use global_typedefs, only: CUSTOM_REAL, path_output
  use mpi_tools, only: exit_MPI
  use paraview

  implicit none

  private

  public :: paraview_write_2d_profiles
  public :: paraview_write_3d_profiles

contains

!============================================================================================================
! This subroutine writes the 2d horizontal profiles (cuts at different hight Z) for Paraview visualization.
!============================================================================================================
subroutine paraview_write_2d_profiles(myrank, name_prefix, ielectrode, nr, ntheta, nz, nzlocal, &
                                      sol, xgrid, ygrid, zgrid)
  ! MPI rank of this process.
  integer, intent(in) :: myrank
  ! Prefix for the file name.
  character(len=*), intent(in) :: name_prefix
  ! Current electrode.
  integer, intent(in) :: ielectrode
  ! Dimensions of the problem (r, theta and z direction, the model is split in the z direction for MPI).
  integer, intent(in) :: nr, ntheta, nz, nzlocal
  ! Solution (e.g. electric potential or permittivity).
  real(kind=CUSTOM_REAL), intent(in) :: sol(0:, 0:, 0:)
  ! Coordinates of points in the grid.
  real(kind=CUSTOM_REAL), intent(in) :: xgrid(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: ygrid(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: zgrid(0:, 0:, 0:)

  character(len=60) :: name
  integer :: k, k1

  do k=1,nzlocal
    k1 = myrank*nzlocal+k

    !if (k1 == 0 .or. k1 == nz/4 .or. k1 == nz/2 .or. k1 == 3*nz/4 .or. k1 == nz) then
    if (k1 == nz/4 .or. k1 == nz/2 .or. k1 == 3*nz/4 .or. &
        k1 == 3*nz/8 .or. k1 == 5*nz/8) then ! Added these heights for two electrode rings.

      ! Make a filename.
      write(name,'("nz",i3.3,"_k",i3.3,"_e",i2.2,"_r",i2.2,".vtk")') nz, k1, ielectrode, myrank

      name = name_prefix//name

      call visualisation_paraview(name, myrank, nr, ntheta, nzlocal, sol, xgrid, ygrid, zgrid, &
                                  0, nr+1, 0, ntheta, k, k+1, 1, 1, 1, 'POINT_DATA')
    endif
  enddo

end subroutine paraview_write_2d_profiles

!================================================================================================
! This subroutine writes the 3d horizontal profiles for Paraview visualization.
!================================================================================================
subroutine paraview_write_3d_profiles(myrank, name_prefix, ielectrode, nr, ntheta, nzlocal, &
                                      sol, xgrid, ygrid, zgrid)
  ! MPI rank of this process.
  integer, intent(in) :: myrank
  ! Prefix for the file name.
  character(len=*), intent(in) :: name_prefix
  ! Current electrode.
  integer, intent(in) :: ielectrode
  ! Dimensions of the problem (r, theta and z direction, the model is split in the z direction for MPI).
  integer, intent(in) :: nr, ntheta, nzlocal
  ! Solution (e.g. electric potential or permittivity).
  real(kind=CUSTOM_REAL), intent(in) :: sol(0:, 0:, 0:)
  ! Coordinates of points in the grid.
  real(kind=CUSTOM_REAL), intent(in) :: xgrid(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: ygrid(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: zgrid(0:, 0:, 0:)

  character(len=60) :: name
  integer :: kmin, kmax

  ! To visualize the model.
  kmin = 1
  kmax = nzlocal + 1

  ! Visualize a vertical profile of the cylinder.
  write(name, '("half_nzlocal",i3.3,"_e",i2.2,"_r",i2.2,".vtk")') nzlocal, ielectrode, myrank

  name = name_prefix//name

  call visualisation_paraview(name, myrank, nr, ntheta, nzlocal, &
                              sol, xgrid, ygrid, zgrid(:, :, (myrank * nzlocal + kmin - 1):(myrank * nzlocal + kmax)), &
                              0, nr + 1, ntheta / 2 + 1, ntheta + 1, kmin, kmax, 1, ntheta / 4, 1, 'POINT_DATA')

end subroutine paraview_write_3d_profiles

end module paraview_ect
