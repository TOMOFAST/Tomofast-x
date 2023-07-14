
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
  use lsqr_solver
  use string, only: str

  implicit none

  private

  public :: solve_problem_loop3d

  private :: read_b_RHS

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
  integer :: ierr

  type(t_sparse_matrix) :: matrix
  integer :: nl, nl_empty

  type(t_inversion_arrays) :: iarr
  real(kind=CUSTOM_REAL), allocatable :: b_RHS(:)
  real(kind=CUSTOM_REAL), allocatable :: delta_model(:)

  if (myrank == 0) print *, "Solving loop3d problem."

  ! Read the sensitivity metadata file to define the nnz.
  call read_sensitivity_metadata(par, nnz, nelements, 1, myrank, nbproc)

  ! Update the nelements.
  par%nelements = nelements
  ipar%nelements = nelements

  ! Allocate memory.
  allocate(b_RHS(par%ndata), source=0._CUSTOM_REAL, stat=ierr)
  allocate(delta_model(par%nelements), source=0._CUSTOM_REAL, stat=ierr)

  if (myrank == 0) print *, "Read nnz, nelements =", nnz, nelements

  ! Memory allocation for auxiliarily inversion arrays.
  call iarr%allocate_aux(par%nelements, par%ndata, par%ndata_components, myrank)

  nl = par%ndata
  nl_empty = 0
  call matrix%initialize(nl, par%nelements, nnz, myrank, nl_empty)

  ! Reading the sensitivity kernel and depth weight from files.
  call read_sensitivity_kernel(par, matrix, iarr%column_weight, ipar%problem_weight(1), 1, myrank, nbproc)

  ! Reading the right-hand side vector b.
  call read_b_RHS(par, b_RHS, myrank, nbproc)

  !-------------------------------------------------------------------------------------
  ! Parallel sparse inversion.
  !-------------------------------------------------------------------------------------
  call matrix%finalize(myrank)

  delta_model = 0._CUSTOM_REAL
  call lsqr_solve(size(b_RHS), size(delta_model), ipar%niter, ipar%rmin, ipar%gamma, &
                  matrix, b_RHS, delta_model, myrank)
  !-------------------------------------------------------------------------------------

end subroutine solve_problem_loop3d

!===================================================================================
! Reads the right-hand side vector b.
!===================================================================================
subroutine read_b_RHS(par, b, myrank, nbproc)
  type(t_parameters_grav), intent(inout) :: par
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), intent(out) :: b(:)

  integer :: Nb, ierr
  character(len=256) :: filename, filename_full
  character(len=256) :: msg

  ! Define the file name.
  filename = "sensit_grav_"//trim(str(nbproc))//"_b"

  filename_full = trim(par%sensit_path)//filename

  if (myrank == 0) print *, "Reading the right-hand side file ", trim(filename_full)

  ! Open the file.
  open(78, file=trim(filename_full), form='unformatted', status='old', action='read', access='stream', &
       iostat=ierr, iomsg=msg)

  if (ierr /= 0) call exit_MPI("Error in opening the right-hand side file! path=" &
                                //trim(filename_full)//", iomsg="//msg, myrank, ierr)

  read(78) Nb
  read(78) b

  close(78)

  if (myrank == 0) print *, "Read Nb =", Nb
end subroutine read_b_RHS

end module problem_loop3d
