
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
! If you use this code for your own research, please cite these articles:
!
! [1] V. Ogarko, J. Giraud, R. Martin, and M. Jessell (2021),
!     "Disjoint interval bound constraints using the alternating direction method of multipliers
!     for geologically constrained inversion: Application to gravity data,"
!     GEOPHYSICS 86: G1-G11. https://doi.org/10.1190/geo2019-0633.1
!
! [2] J. Giraud, V. Ogarko, R. Martin, M. Lindsay, M. Jessell (2021),
!     "Structural, petrophysical and geological constraints in potential field inversion 
!     using the Tomofast-x open-source code", Geoscientific Model Development Discussions, 
!     https://doi.org/10.5194/gmd-2021-14
!
!========================================================================

!===============================================================================================
! A class for solving Loop3d problem.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!===============================================================================================
module problem_loop3d

  use global_typedefs
  use parameters_grav
  use parameters_inversion
  use inversion_arrays
  use sensitivity_gravmag
  use sparse_matrix

  implicit none

  private

  public :: solve_problem_loop3d

contains

!===================================================================================
! Solves loop3d problem.
!===================================================================================
subroutine solve_problem_loop3d(par, ipar, myrank, nbproc)
  type(t_parameters_grav), intent(inout) :: par
  type(t_parameters_inversion), intent(inout) :: ipar
  integer, intent(in) :: myrank, nbproc

  integer(kind=8) :: nnz
  integer :: nelements

  type(t_sparse_matrix) :: matrix
  integer :: nl, nl_empty

  type(t_inversion_arrays) :: iarr

  if (myrank == 0) print *, "Solving loop3d problem."

  ! Read the sensitivity metadata file to define the nnz.
  call read_sensitivity_metadata(par, nnz, nelements, 1, myrank, nbproc)

  ! Update the nelements.
  par%nelements = nelements
  ipar%nelements = nelements

  if (myrank == 0) print *, "Read nnz, nelements =", nnz, nelements

  ! Memory allocation for auxiliarily inversion arrays.
  call iarr%allocate_aux(par%nelements, par%ndata, par%ndata_components, myrank)

  nl = par%ndata
  nl_empty = 0
  call matrix%initialize(nl, par%nelements, nnz, myrank, nl_empty)

  ! Reading the sensitivity kernel and depth weight from files.
  call read_sensitivity_kernel(par, matrix, iarr%column_weight, ipar%problem_weight(1), 1, myrank, nbproc)

end subroutine solve_problem_loop3d

end module problem_loop3d
