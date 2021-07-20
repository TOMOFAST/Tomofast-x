
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

!
! Assembles the matrix, cell-centered finite volumes.
!
! a note on zero padding of the matrix bands:
! We only have A(1,j,k,2)=0 for all j and k because it ocrresponds to zero flux
! at the axis of the cylinder. This is not the case of A(1,j,k,4).
! And we have also A(i,j,1,6)=0 and A(i,j,nz,7)=0 for all i,j because they
! correspond respecrtivley to zero flux at the bottom and top surfaces

!RM RM added this 30 dec 2014 to define the conventions once for all
!The convention we must take once for all is :

!1)flag= 1 and also value =1 or 0 (and thus also phi) on the electrodes at i2
!except in the gaps spaces if they exist
!2)flag=1 and value =0 (and then phi=0) at nr+1 and on the guards
!3) null flux (grad phi . n =0 ) for radius = 0 (we must know exactly if it is
!r(0) or r(1)) ... For radius = 0 (axis) phi is not defined at the centre but the
!flux yes (dphi/dr *S.n =0 because S=0 at radius=0)
!4) grad. phi (top) =0 at k=nz for i<i1 and flux = 0 (bottom) also at k= 1 for i<i1

!This is why we compute solutions for 1<=i<=nr,  1<=j<=ntheta and 1<=k <= nz,
!because we define top or bottom surface fluxes = 0 for k=1 or nz ( 1<=i<i1) or
!for i=1. And solutions are assigned to 0 on guards and electrodes or outer
!screen (i=nr+1), and also at k=0 or nz+1 for i>=i1.

!Be aware that we make a sum of 6 fluxes (2 along r, 2 along theta and 2 along z)
!across the six faces. If we are on top (k=nz, and i<i1) then the top flux is
!zero (element a(i,j,k,7)=0) and phi(nz+1)=0 for i>=i1. Same thing if we are at
!the bottom (k=1, a(i,j,k,6)=0, i<i1) and phi(0)=0 for i>=i1.

!This is the convention. Maybe there are some mistakes but we must be sure that
!we are solving the whole stuff from 1<=i<=nr,  1<=j<=ntheta and 1<=k <= nz and
!that we impose outer Dirichlet BC phi=0 at i=nr+1 and also  k=0 and k=nz+1 for
!i>=i1. The null fluxes are applied at the the null surface located at the axis
!that belongs to cell i=1 to calculate phi(i=1,j,k).
!And null top surface fluxes are applied at cell k=nz to compute phi at k=nz
!(i<i1) and also bottom surface fluxes at cell k=1 to compute phi at k=1 (i<i1).

module matrix

  use global_typedefs, only: CUSTOM_REAL
  use parameters_ect
  use utils, only: avg2, avg4

  implicit none

  private

  public :: matrix_init

contains

!=================================================================================================================
! Building 7-diagonal matrix A.
!=================================================================================================================
pure subroutine matrix_init(a, flag, permit, permit_isolated_tube, permit_air, r, theta, dz, dims, i1, i2, myrank)
  type(t_dimensions), intent(in) :: dims
  integer, intent(in) :: flag(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: permit(0:, 0:, 0:)
  real(kind=CUSTOM_REAL), intent(in) :: r(0:)
  real(kind=CUSTOM_REAL), intent(in) :: theta(0:)
  real(kind=CUSTOM_REAL), intent(in) :: dz(0:)
  real(kind=CUSTOM_REAL), intent(in) :: permit_air,permit_isolated_tube
  integer, intent(in) :: i1, i2
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: a(:, :, :, :)

  ! Local variables.

  ! The alphas are shortcuts for the matrix entries in the seven bands.
  real(kind=CUSTOM_REAL) :: alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7
  real(kind=CUSTOM_REAL) :: per1,per2,per4,per5,per6,per7
  real(kind=CUSTOM_REAL) :: dtheta,dthetae_inv,dthetaw_inv
  real(kind=CUSTOM_REAL) :: dr_n,dr_s
  real(kind=CUSTOM_REAL) :: area_th_e,area_th_w
  real(kind=CUSTOM_REAL) :: area_r_n,area_r_s
  real(kind=CUSTOM_REAL) :: area_z_t,area_z_b

  integer :: i, j, k, k1
  real(kind=CUSTOM_REAL) :: dzavg

  ! Zero out every matrix band.
  a = 0._CUSTOM_REAL

  ! Construction of the matrix.
  do k=1,dims%nzlocal
    k1=myrank*dims%nzlocal+k

    dzavg = 0.5_CUSTOM_REAL*(dz(k1+1)+dz(k1))

    do j=1,dims%ntheta
      dthetae_inv = 1.0_CUSTOM_REAL/abs((theta(j+1)-theta(j)))
      dthetaw_inv = 1.0_CUSTOM_REAL/abs((theta(j)-theta(j-1)))
      dtheta = abs(theta(j+1)-theta(j-1))*0.5_CUSTOM_REAL

      do i=1,dims%nr

        if (flag(i,j,k) /= 1) then

          ! Surface area along angle.
          area_th_e = (r(2*i)-r(2*i-2))*dzavg
          area_th_w = area_th_e

          ! Surface area along radius rdtheta*dz.
          area_r_n = r(2*i)*dtheta*dzavg
          area_r_s = r(2*i-2)*dtheta*dzavg

          ! Surface area along z on top and bottom interfaces rdtheta*dr.
          area_z_t = (r(2*i)**2 -r(2*i-2)**2)*0.5_CUSTOM_REAL*dtheta
          area_z_b = area_z_t

          dr_n = r(2*i+1)-r(2*i-1)
          if (i > 1) dr_s = r(2*i-1)-r(2*i-3)
          if (i == 1) dr_s = dr_n

          ! Matrix elements inside the tube.
          if (i < i1) then
!            per5 = (avg2(dr(2*i)*permit(i+1,j,k),dr(2*i-1)*permit(i,j,k)) + &
!                   avg2(dr(2*i)*permit(i+1,j,k+1),dr(2*i-1)*permit(i,j,k+1))) / (dr(2*i)+dr(2*i-1))
!
!            per1 = (avg2(dr(2*i)*permit(i+1,j-1,k),dr(2*i-1)*permit(i,j-1,k)) + &
!                    avg2(dr(2*i)*permit(i+1,j-1,k+1),dr(2*i-1)*permit(i,j-1,k+1))) / (dr(2*i)+dr(2*i-1))

!            per5=(0.5_CUSTOM_REAL*(dr(2*i)*permit(i+1,j,k)+dr(2*i-1)*permit(i,j,k)) &
!                 + 0.5_CUSTOM_REAL*(dr(2*i)*permit(i+1,j,k+1)+dr(2*i-1)*permit(i,j,k+1))) &
!                 / (dr(2*i)+dr(2*i-1))
!             per1=(0.5_CUSTOM_REAL*(dr(2*i)*permit(i+1,j-1,k)+dr(2*i-1)*permit(i,j-1,k)) &
!                 + 0.5_CUSTOM_REAL*(dr(2*i)*permit(i+1,j-1,k+1)+dr(2*i-1)*permit(i,j-1,k+1))) &
!                 / (dr(2*i)+dr(2*i-1))

            per5 = avg4(permit(i+1,j,k),permit(i,j,k),permit(i+1,j,k+1),permit(i,j,k+1))
            per1 = avg4(permit(i+1,j-1,k),permit(i,j-1,k),permit(i+1,j-1,k+1),permit(i,j-1,k+1))

            per4 = avg4(permit(i+1,j,k),permit(i+1,j-1,k),permit(i+1,j,k+1),permit(i+1,j-1,k+1))
            per2 = avg4(permit(i,j,k),permit(i,j-1,k),permit(i,j,k+1),permit(i,j-1,k+1))

            per7 = avg4(permit(i,j,k+1),permit(i+1,j,k+1),permit(i,j-1,k+1),permit(i+1,j-1,k+1))
            per6 = avg4(permit(i,j,k),permit(i+1,j,k),permit(i,j-1,k),permit(i+1,j-1,k))

          else if (i == i1) then
            per2 = avg4(permit(i,j,k),permit(i,j-1,k),permit(i,j,k+1),permit(i,j-1,k+1))
            per4 = permit_isolated_tube
            per7 = avg2(per2,per4)
            per6 = per7
            per5 = per7
            per1 = per7

          ! Radial matrix elements in the insulated tube and air screen.
          else if (i > i1 .and. i < i2) then
            per2 = permit_isolated_tube
            per4 = permit_isolated_tube
            per7 = permit_isolated_tube
            per6 = permit_isolated_tube
            per5 = permit_isolated_tube
            per1 = permit_isolated_tube

          else if (i == i2) then
            per2 = permit_isolated_tube
            per4 = permit_air
            per7 = avg2(permit_isolated_tube,permit_air)
            per6 = per7
            per5 = per7
            per1 = per7

          else if (i > i2) then
            per2 = permit_air
            per4 = permit_air
            per7 = permit_air
            per6 = permit_air
            per5 = permit_air
            per1 = permit_air
          endif

          ! radius
          alpha2 = -per2/dr_s*area_r_s
          alpha4 = -per4/dr_n*area_r_n
          ! height
          alpha7 = -per7/dz(k1+1)*area_z_t
          alpha6 = -per6/dz(k1)*area_z_b
          ! angle
          alpha5 = -per5/r(2*i-1)*dthetae_inv*area_th_e
          alpha1 = -per1/r(2*i-1)*dthetaw_inv*area_th_w

          ! The flux at the axis is zero (analytically).
          if (i == 1) alpha2 = 0._CUSTOM_REAL

          ! Neumann boundary conditions.
          if (flag(i,j,k-1) == 2) alpha6 = 0._CUSTOM_REAL
          if (flag(i,j,k+1) == 2) alpha7 = 0._CUSTOM_REAL
          if (flag(i+1,j,k) == 2) alpha4 = 0._CUSTOM_REAL

          !! DK DK Roland told us that this is a standard thing in finite volumes:
          !! DK DK sum of all the contributions in the stencil for the diagonal
          alpha3 = -(alpha1+alpha2+alpha4+alpha5+alpha6+alpha7)

          a(i,j,k,1) = alpha1
          a(i,j,k,2) = alpha2
          a(i,j,k,3) = alpha3
          a(i,j,k,4) = alpha4
          a(i,j,k,5) = alpha5
          a(i,j,k,6) = alpha6
          a(i,j,k,7) = alpha7

        else !if (flag(i,j,k) == 1) then
        ! Dirichlet boundary conditions are realized as unit rows in the matrix
        ! and exact values in the RHS, so that the values are inserted into
        ! the solution vector during the solution process of the linear system
        ! (but not during a single matrix-vector multiplication, which would be
        ! the case if we would also insert Dirichlet values directly into the
        ! initial guess for the solution vector)
          a(i,j,k,3) = 1._CUSTOM_REAL
          a(i,j,k,1) = 0._CUSTOM_REAL
          a(i,j,k,2) = 0._CUSTOM_REAL
          a(i,j,k,4) = 0._CUSTOM_REAL
          a(i,j,k,5) = 0._CUSTOM_REAL
          a(i,j,k,6) = 0._CUSTOM_REAL
          a(i,j,k,7) = 0._CUSTOM_REAL
        endif

      enddo
    enddo
  enddo

end subroutine matrix_init

end module matrix

