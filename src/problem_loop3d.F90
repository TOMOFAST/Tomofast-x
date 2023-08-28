
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
  private :: read_Q_size
  private :: write_paraview_model

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
  integer :: A_size, Q_size, it
  real(kind=CUSTOM_REAL) :: cost_data1, cost_data2, cost_data
  real(kind=CUSTOM_REAL) :: model_norm, cost_data_model
  real(kind=CUSTOM_REAL) :: cost_admm1, cost_admm2, cost_admm

  type(t_sparse_matrix) :: matrix
  integer :: nl, nl_empty

  type(t_model) :: Qx
  type(t_admm_method) :: admm_method
  real(kind=CUSTOM_REAL), allocatable :: x0_ADMM(:)
  real(kind=CUSTOM_REAL), allocatable :: b0(:)
  real(kind=CUSTOM_REAL), allocatable :: b(:)
  real(kind=CUSTOM_REAL), allocatable :: Mx(:)
  real(kind=CUSTOM_REAL), allocatable :: model(:)
  real(kind=CUSTOM_REAL), allocatable :: delta_model(:)
  real(kind=CUSTOM_REAL), allocatable :: column_weight(:)

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
  allocate(column_weight(par%nelements), source=0._CUSTOM_REAL, stat=ierr)

  ! Reading the ADMM bounds.
  if (ipar%admm_type > 0) then
    ! Read the size of the Q matrix (number of rows).
    Q_size = read_Q_size(ipar%bounds_ADMM_file(1), myrank)
    if (myrank == 0) print *, 'Q_size =', Q_size

    call Qx%initialize(Q_size, 1, myrank, nbproc)
    call Qx%allocate_bound_arrays(ipar%nlithos, myrank)
    call read_bound_constraints(Qx, ipar%bounds_ADMM_file(1), myrank, nbproc)
  endif

  ! Init the ADMM arrays.
  if (ipar%admm_type > 0) then
    call admm_method%initialize(Q_size, myrank)
    allocate(x0_ADMM(Q_size), source=0._CUSTOM_REAL, stat=ierr)
  endif

  !-------------------------------------------------------------------------------------
  ! Reading the matrix.
  !-------------------------------------------------------------------------------------
  nl = par%ndata
  nl_empty = 0
  call matrix%initialize(nl, par%nelements, nnz, myrank, nl_empty)

  ! Reading the sensitivity kernel and depth weight from files.
  call read_sensitivity_kernel(par, matrix, column_weight, ipar%problem_weight(1), 1, myrank, nbproc)
  call matrix%finalize(myrank)

  ! Reading the right-hand side vector b.
  call read_b_RHS(par, b0, myrank, nbproc)

  ! Size of the A-matrix (in a big matrix with A and Q vertically stacked).
  A_size = par%ndata - Q_size

  if (myrank == 0) print *, 'A_size =', A_size

  ! Apply the ADMM weight to the Q-matrix.
  call matrix%mult_rows(A_size + 1, par%ndata, ipar%rho_ADMM(1))

  ! Starting model.
  model = 0.d0

  ! Create the costs file.
  if (myrank == 0) then
    open(10, file=trim(path_output)//'/costs.txt', access='stream', form='formatted', status='replace', action='write')
  endif

  ! Major inversion loop.
  do it = 1, ipar%ninversions

    if (myrank == 0) print *, 'it =', it

    ! Calculate the forward problem.
    ! Note: Can reuse the b-array for Mx.
    call matrix%mult_vector(model, Mx)

    !-------------------------------------------------------------------------------------
    ! Calculate the ADMM constraints.
    !-------------------------------------------------------------------------------------
    ! Note we scale back the ADMM weight as we want to apply constraints on the original Qx.
    Qx%val(:, 1) = Mx(A_size + 1 : par%ndata) / ipar%rho_ADMM(1)

    call admm_method%iterate_admm_arrays(Qx%nlithos, &
                                         Qx%min_local_bound, Qx%max_local_bound, &
                                         Qx%val(:, 1), x0_ADMM, myrank)

    !-------------------------------------------------------------------------------------
    ! Build the right-hand side part.
    !-------------------------------------------------------------------------------------
    ! Data residual part (including the smoothing term).
    b(1:A_size) = b0(1:A_size) - Mx(1:A_size)

    ! ADMM constraints.
    b(A_size + 1 : par%ndata) = - ipar%rho_ADMM(1) * (Qx%val(:, 1) - x0_ADMM)

    !-------------------------------------------------------------------------------------
    ! Calculate the costs.
    !-------------------------------------------------------------------------------------
    cost_data1 = norm2(b(1:A_size))
    cost_data2 = norm2(b0(1:A_size))
    model_norm = norm2(model)

    ! Calculate the relative cost of the data+reg term.
    cost_data = -1.d0
    if (cost_data2 > 0.d0) cost_data = cost_data1 / cost_data2

    ! Calculate the relative cost of the data+reg term - scaled by the model norm.
    if (model_norm > 0.d0) then
      cost_data_model = cost_data1 / model_norm
    else
      cost_data_model = 0.d0
    endif

    cost_admm1 = norm2(Qx%val(:, 1) - admm_method%z)
    cost_admm2 = norm2(admm_method%z)

    ! Calculate the relative cost of the ADMM term.
    cost_admm = -1.d0
    if (cost_admm2 > 0.d0) cost_admm = cost_admm1 / cost_admm2

    if (myrank == 0) then
      print *, 'cost (data+reg) =', cost_data
      print *, 'cost2 (data+reg) =', cost_data_model
      print *, 'cost (ADMM) =', cost_admm

      ! Write costs to file.
      write(10, *) it, cost_data, cost_data_model, cost_admm
    endif

    !-------------------------------------------------------------------------------------
    ! Parallel sparse inversion.
    !-------------------------------------------------------------------------------------
    delta_model = 0._CUSTOM_REAL
    call lsqr_solve(size(b), size(delta_model), ipar%niter, ipar%rmin, ipar%gamma, &
                    matrix, b, delta_model, myrank)
    !-------------------------------------------------------------------------------------

    ! Update the model.
    model = model + delta_model
  enddo

  if (myrank == 0) close(10)

  ! Write solution to a file.
  call write_final_model(model, "model_final.txt", myrank)

  ! Write constraints Qx to a file.
  call write_final_model(Qx%val(:, 1), "Qx_final.txt", myrank)

  ! Write models for paraview.
  call write_paraview_model(ipar, model, myrank, nbproc)

end subroutine solve_problem_loop3d

!===================================================================================
! Reads the right-hand side vector b.
!===================================================================================
subroutine write_paraview_model(ipar, val, myrank, nbproc)
  type(t_parameters_inversion), intent(in) :: ipar
  real(kind=CUSTOM_REAL), intent(in) :: val(:)
  integer, intent(in) :: myrank, nbproc

  type(t_model) :: vmodel
  integer :: i, j, k, p

  call vmodel%initialize(ipar%nelements, 1, myrank, nbproc)
  call vmodel%grid_full%allocate(ipar%nx, ipar%ny, ipar%nz, myrank)

  ! Init a dummy grid.
  p = 1
  do k = 1, ipar%nz
    do j = 1, ipar%ny
      do i = 1, ipar%nx
        vmodel%grid_full%X1(p) = dble(i)
        vmodel%grid_full%X2(p) = dble(i + 1)
        vmodel%grid_full%Y1(p) = dble(j)
        vmodel%grid_full%Y2(p) = dble(j + 1)
        vmodel%grid_full%Z1(p) = dble(k)
        vmodel%grid_full%Z2(p) = dble(k + 1)

        vmodel%grid_full%i_(p) = i
        vmodel%grid_full%j_(p) = j
        vmodel%grid_full%k_(p) = k
        p = p + 1
      enddo
    enddo
  enddo

  vmodel%val(:, 1) = val
  call model_write(vmodel, "loop", .true., .false., myrank, nbproc)

end subroutine write_paraview_model

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
!===================================================================================
subroutine write_final_model(model, filename, myrank)
  integer, intent(in) :: myrank
  character(len=*), intent(in) :: filename
  real(kind=CUSTOM_REAL), intent(in) :: model(:)

  integer :: i, nelements
  character(len=256) :: filename_full

  if (myrank == 0) then
    filename_full = trim(path_output)//trim(filename)

    print *, 'Writing data to file ', trim(filename_full)

    open(27, file=trim(filename_full), access='stream', form='formatted', status='replace', action='write')

    ! Write the full model array.
    nelements = size(model)
    write(27, *) (model(i), new_line("A"), i = 1, nelements)

    close(27)
  endif
end subroutine write_final_model

!===================================================================================
! Reads the size of Q matrix (the number of rows).
!===================================================================================
function read_Q_size(file_name, myrank) result(Q_size)
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank
  integer :: Q_size, nlithos, ierr
  character(len=256) :: msg

  open(10, file=trim(file_name), status='old', action='read', iostat=ierr, iomsg=msg)
  if (ierr /= 0) call exit_MPI("Error in opening the bound constraints file! path=" &
                 //file_name//" iomsg="//msg, myrank, ierr)

  read(10, *, iostat=ierr) Q_size, nlithos
  close(10)

end function read_Q_size

end module problem_loop3d
