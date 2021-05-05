
!========================================================================
!
!                    T O M O F A S T X  Version 1.0
!                  ----------------------------------
!
!              Main authors: Vitaliy Ogarko, Roland Martin,
!                   Jeremie Giraud, Dimitri Komatitsch.
! CNRS, France, and University of Western Australia.
! (c) CNRS, France, and University of Western Australia. January 2018
!
! This software is a computer program whose purpose is to perform
! capacitance, gravity, magnetic, or joint gravity and magnetic tomography.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

!===============================================================================
! Contains analytical solutions to 3D Laplace's equation in cylinder.
!
! Author: Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!===============================================================================
module laplace

  use global_typedefs
  use trootj

  implicit none

  private

  public :: solution1
  public :: solution2

contains

!===============================================================================
! Reference: M. A. Pinsky, Partial differential equations, Example 3.5.1.
! TODO: add the reference for type=3 general solution.
!===============================================================================
subroutine solution1(nr,ntheta,nz,r,z,u,imax,kmax,type)
  integer, intent(in) :: nr,ntheta,nz,imax,kmax,type
  real(kind=CUSTOM_REAL), intent(in) :: r(0:2*nr+1)
  real(kind=CUSTOM_REAL), intent(in) :: z(0:nz+1)

  real(kind=CUSTOM_REAL), intent(out) :: u(0:nr+1,0:ntheta+1,0:nz+1)

  integer :: i,j,k,n
  real(kind=CUSTOM_REAL) :: radius, Rmax, L, coef, An

  u(:,:,:) = 0.d0

  ! Radius at which we have set the BC's u=1.
  Rmax = r(2*imax-1)

  ! Height at which we have set the BC's u=0.
  L = z(kmax)

  do k=0,kmax
    do j=0,ntheta+1
      do i=0,imax
        if (i == 0) then
          radius = 0.d0
        else
          radius = r(2*i-1)
        endif

        ! Without guards.
        if (type == 1) then
          ! Only odd n-values for series calculation are used.
          do n=1,100,2
            coef = dble(n)*PI/L
            u(i,j,k) = u(i,j,k) + sin(coef*z(k)) * BESSI0(coef*radius) / (dble(n)*BESSI0(coef*Rmax))
          enddo
          u(i,j,k) = u(i,j,k)*4.d0/PI

        ! With guards.
        else if (type == 2) then
          ! Only odd n-values for series calculation are used.
          do n=1,100,2
            coef = dble(n)*PI/L
            u(i,j,k) = u(i,j,k) + sin(coef*z(k)) * BESSI0(coef*radius) / (dble(n)*BESSI0(coef*Rmax)) * &
                                  (cos(dble(n)*PI/4.d0) - cos(3.d0*dble(n)*PI/4.d0))
          enddo
          u(i,j,k) = u(i,j,k)*2.d0/PI

        ! With gradual increase and decrease: Ta(z) = 4z/L(1-z/L)
        else if (type == 3) then
          ! Only odd n-values for series calculation are used.
          do n=1,100,2
            ! Partial An coefficient, excluding a constant factor (4L/PI^3)
            An = - (PI*dble(n)*sin(PI*dble(n)) + 2.d0*cos(PI*dble(n)) - 2.d0) / dble(n**3)

            coef = dble(n)*PI/L
            u(i,j,k) = u(i,j,k) + An * sin(coef*z(k)) * BESSI0(coef*radius) / BESSI0(coef*Rmax)
          enddo
          u(i,j,k) = u(i,j,k)*8.d0/PI**3
        endif

        ! This is analytically correct, since sin(n*PI)=0, but numerically it is not exactly zero.
        if (k == kmax) u(i,j,k) = 0.d0
      enddo
    enddo
  enddo

end subroutine solution1

!======================================================================================
! Reference: M. A. Pinsky, Partial differential equations, Example 3.5.2.
! Note there is a mistake in the book in the last equations, a term with z is missing.
!======================================================================================
subroutine solution2(nr,ntheta,nz,r,z,u,imax,kmax)
  integer, intent(in) :: nr,ntheta,nz,imax,kmax
  real(kind=CUSTOM_REAL), intent(in) :: r(0:2*nr+1)
  real(kind=CUSTOM_REAL), intent(in) :: z(0:nz+1)

  real(kind=CUSTOM_REAL), intent(out) :: u(0:nr+1,0:ntheta+1,0:nz+1)

  integer :: i,j,k,n
  real(kind=CUSTOM_REAL) :: radius, Rmax, L

  integer, parameter :: NK = 50
  INTEGER IER(NK)
  REAL(KIND=8) JZERO(NK)

  u(:,:,:) = 0.d0

  ! Radius at i=imax (where we have set the BC's u=0).
  Rmax = r(2*imax-1)

  ! Height at which we have set the BC's u=1.
  L = z(kmax)

  ! Calculate the roots of Bessel's function.
  call ROOTJ(0,NK,JZERO,IER)

  do n=1,NK
    if (IER(n) /= 0) then
      print *, "ERROR with Bessel function roots calculation!"
      stop
    endif
  enddo

  do k=0,kmax
    do j=0,ntheta+1
      do i=0,imax
        if (i == 0) then
          radius = 0.d0
        else
          radius = r(2*i-1)
        endif

        do n=1,NK
          u(i,j,k) = u(i,j,k) + BESSJ0(radius*JZERO(n)/Rmax) * sinh(JZERO(n)*z(k)/Rmax) / &
                                (JZERO(n) * BESSJ1(JZERO(n)) * sinh(JZERO(n)*L/Rmax))
        enddo

        u(i,j,k) = u(i,j,k)*2.d0

        ! This is analytically correct, since BESSJ0(JZERO(n))=0, but numerically it is not exactly zero.
        if (i == imax) u(i,j,k) = 0.d0
      enddo
    enddo
  enddo

end subroutine solution2

end module laplace
