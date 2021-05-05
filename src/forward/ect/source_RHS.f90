
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

module source

  use global_typedefs
  use parameters_ect

  implicit none

  private

  public :: source_RHS

contains

! Prepare the source vector b for the linear system a*x=b to be solved by the linear solver.
subroutine source_RHS(b,flag,dims,val)
  type(t_dimensions), intent(in) :: dims
  integer, intent(in) :: flag(0:,0:,0:)
  real(kind=CUSTOM_REAL), intent(in) :: val(0:,0:,0:)

  real(kind=CUSTOM_REAL), intent(out) :: b(0:,0:,0:)

  integer :: i,j,k

  ! Initialize b to zero.
  b = 0._CUSTOM_REAL

  ! Compute the source vector b of the system a*x=b.
  do k=1,dims%nzlocal
    do j=1,dims%ntheta
      do i=1,dims%nr
        !! DK DK again, this memory copy is useless (if we fuse b and val everywhere)
        if (flag(i,j,k) == 1) b(i,j,k) = val(i,j,k)
      enddo
    enddo
  enddo

end subroutine source_RHS

end module source

