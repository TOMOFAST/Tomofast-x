
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

module flag

  use global_typedefs
  use parameters_ect

  implicit none

  private

  public :: flag_init_bc

contains

! flag(i,j,k) = 0 for all inner points and some ghost points.
! flag(i,j,k) = 1 if potentials val and phi have to be IMPOSED (to 0 or 1); Dirichlet
! flag(i,j,k) = 2 if gradient potentials val and phi have to be IMPOSED at the outer boundaries(to 0 or 1); Neumann
subroutine flag_init_bc(flag,nr,ntheta,nzlocal,&
                        i1,i2,nel,elecs,guards,myrank,nbproc)

  ! Problem size in each dimension.
  integer, intent(in) :: nr, ntheta, nzlocal, nel
  integer, intent(in) :: i1, i2
  integer, intent(in) :: myrank, nbproc
  type(t_electrode), intent(in) :: elecs(:)
  type(t_guards), intent(in) :: guards

  ! Flag array to be filled, including ghost cells.
  integer, intent(out) :: flag(0:,0:,0:)

  ! Local variables.
  integer :: i, j, k, k1, jelectrode

  ! debug comments from RM on why flags are weird for multigrid:
  !
  ! What I do think is that we must think all the code in terms of even dimensions or not.
  ! i1 fron r=0 to R1 should be even
  ! i2 in R2 too
  ! nr too
  ! the number of points on electrodes, gaps and guards  too
  ! And then we must adapt the r, theta and z accordingly. We must find a way to do that
  ! by subdomains, even if it is not exactly the same as in finite elements.
  !
  ! The thing might come from the fact that I must introduce gaps between electrodes
  ! and between electrodes and guards.
  !
  ! Between electrodes: In flag_init_BC.f90 we detect the points
  ! on the electrodes and I withdraw 2 points left and 2 points at right tip of each
  ! electrode to make the gap in a way or another. We should adapt the theta I guess
  ! in order to have exactly the right length of the electrodes (the physical gap
  ! spacing is given in the figure of the proguide).
  !
  ! Between guards and electrodes: same thing but we start guards at a izspace
  ! number of points above or below the electrodes. Maybe a bad shift of +1 or -1
  ! is hidden there and then again the exact heigt of the electrodes and guards
  ! should given according to the right spacing.



  ! note on the dimensions:
  ! k=0, k=nzlocal+1    ghost cells for MPI in z direction
  !                     (except for k=0 and myrank=0: BCs at the very bottom of the device;
  !                     and k=lzlocal+1 and myrank=nbproc-1: BCs at the very top of the device)
  ! j=0, j=ntheta+1     used for periodicity
  ! i=0, i=nr+1         i=0 is the singularity at the cyclinder axis
  !                     i=nr+1 corresponds to outermost radius and Dirichlet BCs
  !
  ! Only Neumann conditions and Dirichlet conditions at the electrodes are actually
  ! inserted in the matrix and RHS, all other Dirichlet conditions are treated as
  ! known values and are excluded a priori (they are used during flux computations
  ! in matrix_init() and source_RHS() to get values at the immediate neighbours correct,
  ! but they are never actually included in the linear system)

  ! Boundary conditions coming from the cylinder itself.
  do k=0,nzlocal+1
    do j=0,ntheta+1
      do i=0,nr+1

        ! Mark everything as "inner point".
        flag(i,j,k) = 0

        ! Ghost cells on the very exterior (r direction) boundary of the cylinder: Dirichlet.
        ! will be removed later
        if (i == nr+1) flag(i,j,k) = 1

        ! At the bottom end of the device: Neumann.
        ! Exception: Dirichlet for all i>=i1 (outward from i1).
        if (myrank == 0 .and. k == 0) then
          flag(i,j,k) = 2
          if (i >= i1) flag(i,j,k) = 1
        endif

        ! Same for the top end of the device.
        if (myrank == nbproc-1 .and. k == nzlocal+1) then
          flag(i,j,k) = 2
          if (i >= i1) flag(i,j,k) = 1
        endif

      enddo
    enddo
  enddo


  ! Boundary conditions coming from the electrodes.
  ! Note that actual values (if an electrode fires or receives) are set later.
  do jelectrode=1,nel
    do k=0,nzlocal+1
      k1 = myrank*nzlocal+k
      do j=0,ntheta+1
        if (j >= elecs(jelectrode)%ithetamin .and. j <= elecs(jelectrode)%ithetamax .and. &
            k1 >= elecs(jelectrode)%izmin .and. k1 <= elecs(jelectrode)%izmax) &
            flag(i2,j,k) = 1
      enddo
    enddo
  enddo

  ! Boundary conditions coming from the guards.
  do k=0,nzlocal+1
    k1 = myrank*nzlocal+k
    do j=0,ntheta+1
      if (k1 <= guards%lower_izmax) flag(i2,j,k) = 1
      if (k1 >= guards%upper_izmin) flag(i2,j,k) = 1
    enddo
  enddo

  ! Ensure periodicity of the flag array in the theta direction, including ghost cells.
  do k=0,nzlocal+1
    do i=0,nr+1
      flag(i,ntheta+1,k) = flag(i,1,k)
      flag(i,0,k) = flag(i,ntheta,k)
    enddo
  enddo

end subroutine flag_init_bc

end module flag

