
!========================================================================
!
!                          T o m o f a s t - x
!                        -----------------------
!
!           Authors: Vitaliy Ogarko, Jeremie Giraud, Roland Martin.
!
!               (c) 2021 The University of Western Australia.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

module global_typedefs

  use mpi

  implicit none

  ! Try this if 'use mpi' does not work, e.g., there is no 'module load' system.
  !include 'mpif.h'

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

  ! The precision of the stored sensitivity kernel values (single / double).
  integer, parameter :: MATRIX_PRECISION = SIZE_REAL

  ! The precision used for calculating the sensitivity kernel (mag only).
  integer, parameter :: SENSIT_REAL = SIZE_DOUBLE
  !-----------------------------------------------------------------------

  ! Path to the output folder, set in the Parfile.
  character(len=256) :: path_output

  ! A factor to multiply data by on input and output, e.g 1.0E-5 if input is in mGal or 1.0E-6 if input is in um/s2, default = 1.0 (m/s2).
  real(kind=CUSTOM_REAL) :: data_units_mult

  ! A factor to multiply model by on input and output, e.g 1.0E3 if input model is T/m3, default = 1.0 (kg/m3).
  real(kind=CUSTOM_REAL) :: model_units_mult

  ! PI
  real(kind=CUSTOM_REAL), parameter :: PI = 3.1415926535897932385_CUSTOM_REAL

  ! Tolerance for comparing real numbers in unit tests.
  real(kind=CUSTOM_REAL), parameter :: tol = merge(1.e-12_CUSTOM_REAL, 1.e-6_CUSTOM_REAL, MATRIX_PRECISION == SIZE_DOUBLE)

  ! An auxiliarily type that stores a 1D array.
  type t_real1d
    real(kind=CUSTOM_REAL), allocatable :: val(:)
  end type t_real1d

  ! An auxiliarily type that stores a 2D array.
  type t_real2d
    real(kind=CUSTOM_REAL), allocatable :: val(:, :)
  end type t_real2d

end module global_typedefs

