
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


module global_typedefs

  !use mpi
  ! VO VO: Comment from Dimitri from partII revision 405:
  ! replaced "use mpi" back with "include 'mpif.h'" because of problems on two systems
  ! (that have no 'module load' system, which is common);
  ! please do NOT change that back in the future.

  implicit none

  include 'mpif.h'

  !-----------------------------------------------------------------------
  ! Choose between single and double precision for the whole code;
  ! solver in single or double precision depending on the machine (4 or 8 bytes)

  ! set to MPI_REAL to run in single precision
  ! set to MPI_DOUBLE_PRECISION to run in double precision
  ! integer, parameter :: CUSTOM_MPI_TYPE = MPI_REAL
  integer, parameter :: CUSTOM_MPI_TYPE = MPI_DOUBLE_PRECISION

  integer, parameter :: SIZE_REAL = 4, SIZE_DOUBLE = 8
  integer, parameter :: SIZE_INT  = 4

  ! set to SIZE_REAL to run in single precision
  ! set to SIZE_DOUBLE to run in double precision (increases memory size by 2 and makes code 30% slower or so)
  ! integer, parameter :: CUSTOM_REAL = SIZE_REAL
  integer, parameter :: CUSTOM_REAL = SIZE_DOUBLE
  !-----------------------------------------------------------------------

  ! Path to the output folder, set in the Parfile.
  character(len=256) :: path_output

  ! PI
  real(kind=CUSTOM_REAL), parameter :: PI = 3.1415926535897932385_CUSTOM_REAL

  ! Tolerance for comparing real numbers in unit tests.
  real(kind=CUSTOM_REAL), parameter :: tol = 1.e-12_CUSTOM_REAL

  ! Type of norm used for convergence control.
  integer, parameter :: NORM_L2 = 1
  integer, parameter :: NORM_MAX = 2

  ! To have arrays of pointers (which are not directly supported in Fortran),
  ! we need arrays of types which only contain pointers.
  type real_p
    real(kind=CUSTOM_REAL), dimension(:), pointer :: p
  end type real_p

  type real2_p
    real(kind=CUSTOM_REAL), dimension(:,:), pointer :: p
  end type real2_p

  type real3_p
    real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: p
  end type real3_p

  type real4_p
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), pointer :: p
  end type real4_p

  type int_p
    integer, dimension(:), pointer :: p
  end type int_p

  type int3_p
    integer, dimension(:,:,:), pointer :: p
  end type int3_p

end module global_typedefs

