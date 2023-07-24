
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
  use model
  use model_IO
  use admm_method
  use string, only: str

  implicit none

  private

  public :: solve_problem_loop3d

  private :: read_b_RHS
  private :: write_final_model

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
  integer :: Na

  type(t_sparse_matrix) :: matrix
  integer :: nl, nl_empty

  type(t_model) :: Qx
  type(t_inversion_arrays) :: iarr
  type(t_admm_method) :: admm_method
  real(kind=CUSTOM_REAL), allocatable :: x0_ADMM(:)
  real(kind=CUSTOM_REAL), allocatable :: b0(:)
  real(kind=CUSTOM_REAL), allocatable :: b(:)
  real(kind=CUSTOM_REAL), allocatable :: Mx(:)
  real(kind=CUSTOM_REAL), allocatable :: model(:)
  real(kind=CUSTOM_REAL), allocatable :: delta_model(:)

  if (myrank == 0) print *, "Solving loop3d problem."

  ! Read the sensitivity metadata file to define the nnz.
  call read_sensitivity_metadata(par, nnz, nelements, 1, myrank, nbproc)

  ! Update the nelements.
  par%nelements = nelements
  ipar%nelements = nelements

  ! Some extra initialisation.
  ipar%ndata(1) = par%ndata
  ipar%nx = par%nx
  ipar%ny = par%ny
  ipar%nz = par%nz
  ipar%nelements_total = ipar%nx * ipar%ny * ipar%nz
  ipar%nmodel_components = 1
  ipar%ndata_components = 1
  ipar%compression_type = par%compression_type

  if (myrank == 0) print *, "Read nnz, nelements =", nnz, nelements

  ! Allocate memory.
  allocate(b0(par%ndata), source=0._CUSTOM_REAL, stat=ierr)
  allocate(b(par%ndata), source=0._CUSTOM_REAL, stat=ierr)
  allocate(Mx(par%ndata), source=0._CUSTOM_REAL, stat=ierr)
  allocate(model(par%nelements), source=0._CUSTOM_REAL, stat=ierr)
  allocate(delta_model(par%nelements), source=0._CUSTOM_REAL, stat=ierr)

  ! Memory allocation for auxiliarily inversion arrays.
  call iarr%allocate_aux(par%nelements, par%ndata, par%ndata_components, myrank)

  ! Reading the ADMM bounds.
  if (ipar%admm_type > 0) then
    call Qx%initialize(ipar%nelements, ipar%nmodel_components, myrank, nbproc)
    call Qx%allocate_bound_arrays(ipar%nlithos, myrank)
    call read_bound_constraints(Qx, ipar%bounds_ADMM_file(1), myrank, nbproc)
  endif

  ! Init the ADMM arrays.
  if (ipar%admm_type > 0) then
    call admm_method%initialize(par%nelements, myrank)
    allocate(x0_ADMM(par%nelements), source=0._CUSTOM_REAL, stat=ierr)
  endif

  !-------------------------------------------------------------------------------------
  ! Reading the matrix.
  !-------------------------------------------------------------------------------------
  nl = par%ndata
  nl_empty = 0
  call matrix%initialize(nl, par%nelements, nnz, myrank, nl_empty)

  ! Reading the sensitivity kernel and depth weight from files.
  call read_sensitivity_kernel(par, matrix, iarr%column_weight, ipar%problem_weight(1), 1, myrank, nbproc)
  call matrix%finalize(myrank)

  ! Reading the right-hand side vector b.
  call read_b_RHS(par, b0, myrank, nbproc)

  ! Size of the A-matrix (in a big matrix with A and Q vertically stacked).
  Na = par%ndata - ipar%nelements_total

  ! Starting model.
  model = 0.d0

  ! Calculate the forward problem.
  ! TODO: Can reuse the b-array for Mx.
  call matrix%mult_vector(model, Mx)

  !-------------------------------------------------------------------------------------
  ! Calculate the ADMM constraints.
  !-------------------------------------------------------------------------------------
  Qx%val(:, 1) = Mx(Na + 1 : par%ndata)

  call admm_method%iterate_admm_arrays(Qx%nlithos, &
                                       Qx%min_local_bound, Qx%max_local_bound, &
                                       Qx%val(:, 1), x0_ADMM, myrank)

  !-------------------------------------------------------------------------------------
  ! Build the right-hand side part.
  !-------------------------------------------------------------------------------------
  ! Data residual part (including the smoothing term).
  b(1:Na) = b0(1:Na) - Mx(1:Na)

  ! ADMM constraints.
  b(Na + 1 : par%ndata) = - (Qx%val(:, 1) - x0_ADMM)

  !-------------------------------------------------------------------------------------
  ! Parallel sparse inversion.
  !-------------------------------------------------------------------------------------
  delta_model = 0._CUSTOM_REAL
  call lsqr_solve(size(b), size(delta_model), ipar%niter, ipar%rmin, ipar%gamma, &
                  matrix, b, delta_model, myrank)
  !-------------------------------------------------------------------------------------

  ! Write result to a file.
  call write_final_model(delta_model, myrank)

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
  if (myrank == 0) print *, "Read Nb =", Nb

  read(78) b(1:Nb)

  close(78)

end subroutine read_b_RHS

!===================================================================================
! Writes the final model.
! TODO: Use this temporariry - to replace with the model class writers.
!===================================================================================
subroutine write_final_model(model, myrank)
  integer, intent(in) :: myrank
  real(kind=CUSTOM_REAL), intent(in) :: model(:)

  integer :: i, nelements
  character(len=256) :: filename, filename_full

  if (myrank == 0) then
    filename = "model_final.txt"
    filename_full = trim(path_output)//filename

    print *, 'Writing the full model to file ', trim(filename_full)

    open(27, file=trim(filename_full), access='stream', form='formatted', status='replace', action='write')

    ! Write the full model array.
    nelements = size(model)
    write(27, *) (model(i), new_line("A"), i = 1, nelements)

    close(27)
  endif
end subroutine write_final_model

end module problem_loop3d
