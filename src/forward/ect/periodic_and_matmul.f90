
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

module m_matmul

  use global_typedefs, only: CUSTOM_REAL
  use utils, only: enforce_pb

  implicit none

  private

  public :: periodic_and_matmul
  private :: matmult

contains

! DG DG and DK DK to RM RM: it would be nice to have a similar routine in this file called defect()
! DG DG and DK DK to RM RM: that computes r=b-Ax. But unfortunately, this is not possible,
! DG DG and DK DK to RM RM: without having this routine perform the MPI_sendreceive() internally,
! DG DG and DK DK to RM RM: which means code restructuring, and then the benefit of not moving
! DG DG and DK DK to RM RM: an array twice through memory would be lost anyway???

!======================================================================================
! Computes only the local (at current CPU) part of r = Ax.
!
! NOTE 1: Arrays A and x have both intent 'in'. Therefore we split this function from
! periodic_and_matmul() where x has intent 'inout'. This potentially allows compiler for
! better optimization/vectorization of the loop.
!
! NOTE 2: The outer values (i,j,k=0, i=nr+1 etc) of r array should be set (zero out)
! outside this subroutine once and for all. Since having r = 0 inside the subroutine is
! too expensive (takes ~5% of the total CPU time).
!======================================================================================
pure subroutine matmult(A, x, r, nr, ntheta, nzlocal)
  integer, intent(in) :: nr, ntheta, nzlocal
  ! Vitaliy checked with callgrind and Allinea profilers that if we use assumed-shape arrays here,
  ! then it slows down the performance.
  real(kind=CUSTOM_REAL), intent(in) :: A(nr, ntheta, nzlocal, 7)
  real(kind=CUSTOM_REAL), intent(in) :: x(0:nr+1, 0:ntheta+1, 0:nzlocal+1)
  real(kind=CUSTOM_REAL), intent(out) :: r(0:nr+1, 0:ntheta+1, 0:nzlocal+1)

  ! Loop variables.
  integer i, j, k

  ! General multiplication: all BC's are included in a(i,j,k,:) coefficients.
  ! Vitaliy checked that this loop is vectorized by GNU compiler (gcc 4.9.3),
  ! by looking at the output with the flag -fopt-info-vec-optimized=vec.info.
  ! Also checked that the inner loop is vectorized by Intel compiler.
  do k = 1, nzlocal
    do j = 1, ntheta
      do i = 1, nr
        r(i, j, k) = a(i, j, k, 3) * x(i, j, k) + &
                     a(i, j, k, 2) * x(i - 1, j, k) + a(i, j, k, 4) * x(i + 1, j, k) + &
                     a(i, j, k, 1) * x(i, j - 1, k) + a(i, j, k, 5) * x(i, j + 1, k) + &
                     a(i, j, k, 6) * x(i, j, k - 1) + a(i, j, k, 7) * x(i, j, k + 1)
      enddo
    enddo
  enddo

end subroutine matmult

!======================================================================================
! Enforces periodic condition in theta on the x vector,
! then computes only the local part of r = Ax.
!======================================================================================
pure subroutine periodic_and_matmul(A, x, r, nr, ntheta, nzlocal)
  integer, intent(in) :: nr, ntheta, nzlocal
  ! Vitaliy checked with callgrind and Allinea profilers that if we use assumed-shape arrays here,
  ! then it slows down the performance.
  real(kind=CUSTOM_REAL), intent(in) :: A(nr, ntheta, nzlocal, 7)
  real(kind=CUSTOM_REAL), intent(inout) :: x(0:nr+1, 0:ntheta+1, 0:nzlocal+1)
  real(kind=CUSTOM_REAL), intent(out) :: r(0:nr+1, 0:ntheta+1, 0:nzlocal+1)

  ! Enforce periodic condition in theta on the x vector.
  call enforce_pb(nr, ntheta, nzlocal, x)

  ! Compute the local part of r = Ax.
  call matmult(A, x, r, nr, ntheta, nzlocal)

end subroutine periodic_and_matmul

end module m_matmul

