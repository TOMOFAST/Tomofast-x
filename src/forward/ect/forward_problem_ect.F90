
!========================================================================
!
!                          T o m o f a s t - x
!                        -----------------------
!
!              Main authors: Vitaliy Ogarko, Roland Martin,
!                   Jeremie Giraud, Dimitri Komatitsch.
!
! (c) CNRS, France, and University of Western Australia. January 2018
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

module forward_problem_ect

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use parameters_ect
  use geometry
  use flag
  use matrix
  use source
  use solution
  use conjugate_gradient
  use sanity_check
  use data_ect
  use capacitance
  use sensitivity_ect
  use model_ect

  implicit none

  private

  public :: solve_forward_problem_ect

contains

!=================================================================================================
! Main function to solve forward ECT problem.
!=================================================================================================
subroutine solve_forward_problem_ect(par, sensit, cdata, model, column_weight, damping_weight, &
                                    myrank, nbproc)
  ! MPI variables.
  integer, intent(in) :: myrank, nbproc
  ! Input parameters
  type(t_parameters_ect), intent(in) :: par

  ! Permittivity model.
  real(kind=CUSTOM_REAL), intent(inout) :: model(:)
  ! Sensitivity matrix.
  real(kind=CUSTOM_REAL), intent(out) :: sensit(:, :)
  ! Capacitance data.
  real(kind=CUSTOM_REAL), intent(out) :: cdata(:)
  ! Inversion weights to scale sensitivity matrix columns (the cell volume).
  real(kind=CUSTOM_REAL), intent(out) :: column_weight(:)
  ! Inversion weights to scale damping.
  real(kind=CUSTOM_REAL), intent(out) :: damping_weight(:)

  !----------------------------------------------------------------------
  ! Local variables.

  ! various loop variables and counters
  integer :: i, j, k, iter_pcg
  integer :: ierr

  !
  ! Multigrid variables, loop counters and arrays, containing per-level dimensions
  !
  ! TODO: use t_dimensions in the parameters of all functions which are using at_level variables.
  type(t_dimensions), dimension(:), allocatable :: dims_at_level
  integer :: ilevel_fine, ilevel
  integer :: nr, ntheta, nz, nzlocal

  !
  ! GEOMETRY
  !
  integer, dimension(:), allocatable :: izspace

  ! theta, r and z angle, radius and side for each point
  type(real_p), dimension(:), allocatable :: theta
  ! Step along r et z.
  type(real_p), dimension(:), allocatable :: dr, dz, r, z
  ! Coordinates of points in the grid.
  type(real3_p), dimension(:), allocatable :: xgrid, ygrid, zgrid
  ! Indices in r-direction of the various radii of the device.
  integer, dimension(:), allocatable :: i1, i2
  ! Grid parameters of the electrodes at different multigrid levels.
  type(t_electrode_p), dimension(:), allocatable :: elecs
  ! Grid parameters of the guards at different multigrid levels.
  type(t_guards), dimension(:), allocatable :: guards

  ! Original theta before mesh refinement.
  type(real_p), dimension(:), allocatable :: theta0
  ! Original electrode grid parameters.
  type(t_electrode_p), dimension(:), allocatable :: elecs0

  ! Store mapped indexes of original and refined theta arrays.
  type(int_p), dimension(:), allocatable :: ithetamap

  !
  ! SOLUTION
  !
  ! counters
  integer :: ielectrode
  ! flags/markers for boundary conditions
  type(int3_p), dimension(:), allocatable :: flag
  ! value of the potential at the electrodes
  type(real3_p), dimension(:), allocatable :: val
  ! solution phi
  type(real3_p), dimension(:), allocatable :: phi

  !
  ! MATRIX and vector arrays
  !
  ! matrix
  type(real4_p), dimension(:), allocatable :: a
  ! source right and side b_RHS of Ax = b equation
  type(real3_p), dimension(:), allocatable :: b_RHS !! DK DK fuse with val? or maybe not? will see ! DG DG no, makes MG much more complicated
  ! temporary solution
  type(real3_p), dimension(:), allocatable :: x2 !! DK DK fuse with phi???

  !
  ! Configuration of the linear solver
  !
  ! auxiliary arrays for conjugate gradient solver
  type (t_pcg_auxarrays) :: pcg_auxarrays
  ! auxiliary arrays for multigrid solver
  !type (t_mg_auxarrays), dimension(:), allocatable :: mg_auxarrays
  logical :: suppress_output

  !
  !  MODEL
  !
  ! position, radius and permittivity of bubbles
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xb, yb, zb, rb, permit_bubble

  ! model of permittivity on the mesh
  type(real3_p), dimension(:), allocatable :: permit

  !
  ! Variables needed for sensitivity kernel calculation
  !
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: phisensit
  ! Calculated capacitance data set.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: capacity
  ! Cell volume.
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: cell_volume
  ! Average volume of all cells.
  real(kind=CUSTOM_REAL) :: cell_volume_avg

#ifdef SUPPRESS_OUTPUT
  suppress_output = .true.
#else
  suppress_output = .false.
#endif

!=============================================================
!==== BUBBLES ================================================
!=============================================================

  ! If needed, allocate memory for the bubble model and read it;
  ! again, only the master does file IO and then broadcasts.
  if (par%num_bubbles > 0) then
    ! TODO: create a type t_bubble
    allocate(xb(par%num_bubbles))
    allocate(yb(par%num_bubbles))
    allocate(zb(par%num_bubbles))
    allocate(rb(par%num_bubbles))
    allocate(permit_bubble(par%num_bubbles))

    ! TODO: Move to model.f90
    if (myrank == 0) then
      open(unit=11, file=par%filename_bubbles, status='old', action='read', iostat=ierr)
      if (ierr /= 0) call exit_MPI("Error in reading the bubbles file!", myrank, ierr)

      do j = 1, par%num_bubbles
        read(11,*) xb(j), yb(j), zb(j), rb(j), permit_bubble(j)
        if (.not. suppress_output .and. myrank == 0) &
          print *,'bubble ', j, xb(j), yb(j), zb(j), rb(j), permit_bubble(j)
      enddo
      close(11)
    endif

    call MPI_Bcast(xb,par%num_bubbles,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(yb,par%num_bubbles,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(zb,par%num_bubbles,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(rb,par%num_bubbles,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(permit_bubble,par%num_bubbles,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  endif

!=============================================================
!==== PARTITIONING FOR GEOMETRIC MULTRIGRID ==================
!=============================================================

  ! ilevel_fine is always the finest grid, ilevel_coarse is determined by
  ! subsequently halving all dimensions, and stopping once one of them becomes
  ! too small and/or no longer divisible by two
  ilevel_fine = 1
  ! DG DG current implementation sets the number of multigrid levels via parfile
  ! DG DG current implementation does not support semi-coarsening yet

  ! sanity checks
  ! DG DG will move some of these to the new sanity checker once semi-coarsening
  !       is implemented
  if (par%ilevel_coarse < ilevel_fine) &
    call exit_MPI("coarser levels cannot be larger than finer levels; exiting...",myrank,04)
  if (par%linear_solver/=LINSOLV_MG .and. par%ilevel_coarse/=ilevel_fine) &
    call exit_MPI("different levels are only possible for multigrid solvers; exiting...",myrank,05)
  if (par%linear_solver==LINSOLV_MG .and. &
     (par%dims%nr/=par%dims%ntheta .or. par%dims%nr/=par%dims%nz)) &
    call exit_MPI("for multigrid,nr ntheta and nz must be the same; exiting...",myrank,06)

!=============================================================
!==== DATA ALLOCATION ========================================
!=============================================================

  ! allocate arrays with per-level dimensions
  ! Fortran knowledge exploitation: in the singlegrid case, ilevel_coarse is
  ! always ==1, so this does not allocate 2D arrays in that case. Not an
  ! efficiency deficit, since the compiler is smart enough (and the Fortran
  ! standard is smart enough) to get rid of (*,1) indirections at compile time
  allocate(dims_at_level(par%ilevel_coarse))

  ! TODO: Overload the / operator for t_dimensions to perform this in one operation.
  ! Determine dimensions and perform sanity checks.
  do ilevel=ilevel_fine, par%ilevel_coarse
    dims_at_level(ilevel)%nr      = par%dims%nr / 2**(ilevel-1)
    dims_at_level(ilevel)%ntheta  = par%dims%ntheta / 2**(ilevel-1)
    if(mod(dims_at_level(ilevel)%ntheta, par%nel/par%nrings) /= 0) then
      call exit_MPI("ntheta at a coarser level is not a multiple of nel; exiting...",myrank,07)
    endif
    dims_at_level(ilevel)%nz      = par%dims%nz / 2**(ilevel-1)
    dims_at_level(ilevel)%nzlocal = par%dims%nzlocal / 2**(ilevel-1)
    dims_at_level(ilevel)%kguards = par%dims%kguards / 2**(ilevel-1)

    ! Print out problem size information.
    if (.not. suppress_output .and. myrank == 0) print *, 'lvl,nr,ntheta,nz,nzloc,kguards:',&
      ilevel,dims_at_level(ilevel)%nr,dims_at_level(ilevel)%ntheta,dims_at_level(ilevel)%nz,&
      dims_at_level(ilevel)%nzlocal,dims_at_level(ilevel)%kguards
  enddo
  ! DG DG TODO note to self sanity checks if smallest problem n*_at_*(ilevel_coarse) <= some_threshold


  ! Allocate pointer arrays for everything that is level-dependent.
  allocate(izspace(par%ilevel_coarse))
  allocate(theta(par%ilevel_coarse))
  allocate(dr(par%ilevel_coarse))
  allocate(dz(par%ilevel_coarse))
  allocate(r(par%ilevel_coarse))
  allocate(z(par%ilevel_coarse))
  allocate(xgrid(par%ilevel_coarse))
  allocate(ygrid(par%ilevel_coarse))
  allocate(zgrid(par%ilevel_coarse))
  allocate(flag(par%ilevel_coarse))
  allocate(i1(par%ilevel_coarse))
  allocate(i2(par%ilevel_coarse))
  allocate(val(par%ilevel_coarse))
  allocate(phi(par%ilevel_coarse))
  allocate(a(par%ilevel_coarse))
  allocate(b_RHS(par%ilevel_coarse))
  allocate(x2(par%ilevel_coarse))
  allocate(permit(par%ilevel_coarse))
  allocate(elecs(par%ilevel_coarse))
  allocate(guards(par%ilevel_coarse))

  if (par%irefine == 1) then
    allocate(theta0(par%ilevel_coarse))
    allocate(ithetamap(par%ilevel_coarse))
    allocate(elecs0(par%ilevel_coarse))
  endif

  ! Multigrid needs more aux arrays.
  if (par%linear_solver == LINSOLV_MG) then
    !allocate(mg_auxarrays(par%ilevel_coarse))
  endif

  ! Allocate data arrays for each level.
  ! Note that only in the multigrid case, ilevel_fine /= ilevel_coarse.
  do ilevel = ilevel_fine, par%ilevel_coarse
    nr = dims_at_level(ilevel)%nr
    ntheta = dims_at_level(ilevel)%ntheta
    nz = dims_at_level(ilevel)%nz
    nzlocal = dims_at_level(ilevel)%nzlocal

    allocate(theta(ilevel)%p(0:ntheta+1))
    allocate(dr(ilevel)%p(0:2*nr+1))
    allocate(dz(ilevel)%p(0:nz+1))
    allocate(r(ilevel)%p(0:2*nr+1))
    allocate(z(ilevel)%p(0:nz+1))
    allocate(xgrid(ilevel)%p(0:nr+1,0:ntheta+1,0:nz+1))
    allocate(ygrid(ilevel)%p(0:nr+1,0:ntheta+1,0:nz+1))
    allocate(zgrid(ilevel)%p(0:nr+1,0:ntheta+1,0:nz+1))
    allocate(flag(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(val(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(phi(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(a(ilevel)%p(1:nr,1:ntheta,1:nzlocal,7))
    allocate(b_RHS(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(x2(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(permit(ilevel)%p(0:nr+1,0:ntheta+1,0:nzlocal+1))
    allocate(elecs(ilevel)%p(par%nel))

    if (par%irefine == 1) then
      allocate(theta0(ilevel)%p(0:ntheta+1))
      allocate(ithetamap(ilevel)%p(0:ntheta+1))
      allocate(elecs0(ilevel)%p(par%nel))
    endif

!    ! Auxiliary arrays only needed for multigrid.
!    if (par%linear_solver == LINSOLV_MG) then
!      allocate(mg_auxarrays(ilevel)%c(0:nr+1,0:ntheta+1,0:nzlocal+1))
!      allocate(mg_auxarrays(ilevel)%d(0:nr+1,0:ntheta+1,0:nzlocal+1))
!      allocate(mg_auxarrays(ilevel)%inv_diagA(1:nr,1:ntheta,1:nzlocal))
!      allocate(mg_auxarrays(ilevel)%cprime_r(nr))
!      allocate(mg_auxarrays(ilevel)%dprime_r(nr))
!      allocate(mg_auxarrays(ilevel)%xprime_r(nr))
!      allocate(mg_auxarrays(ilevel)%cprime_theta(ntheta))
!      allocate(mg_auxarrays(ilevel)%dprime_theta(ntheta))
!      allocate(mg_auxarrays(ilevel)%xprime_theta(ntheta))
!      allocate(mg_auxarrays(ilevel)%cprime_z(nzlocal))
!      allocate(mg_auxarrays(ilevel)%dprime_z(nzlocal))
!      allocate(mg_auxarrays(ilevel)%xprime_z(nzlocal))
!    endif
  enddo

  ! Arrays in case multigrid solver using PCG as coarse grid solver.
  if (par%linear_solver == LINSOLV_MG .and. par%coarse_solver == LINSOLV_PCG) then
    nr = dims_at_level(par%ilevel_coarse)%nr
    ntheta = dims_at_level(par%ilevel_coarse)%ntheta
    nzlocal = dims_at_level(par%ilevel_coarse)%nzlocal

    allocate(pcg_auxarrays%p(0:nr+1, 0:ntheta+1, 0:nzlocal+1))
    allocate(pcg_auxarrays%r(0:nr+1, 0:ntheta+1, 0:nzlocal+1))
    allocate(pcg_auxarrays%z(0:nr+1, 0:ntheta+1, 0:nzlocal+1))
    allocate(pcg_auxarrays%tmp(0:nr+1, 0:ntheta+1, 0:nzlocal+1))
    allocate(pcg_auxarrays%inv_diagA(1:nr, 1:ntheta, 1:nzlocal))
  endif

  ! Arrays that are only needed for conjugate gradient solver.
  if (par%linear_solver == LINSOLV_PCG) then
    allocate(pcg_auxarrays%p(0:par%dims%nr+1, 0:par%dims%ntheta+1, 0:par%dims%nzlocal+1))
    allocate(pcg_auxarrays%r(0:par%dims%nr+1, 0:par%dims%ntheta+1, 0:par%dims%nzlocal+1))
    allocate(pcg_auxarrays%z(0:par%dims%nr+1, 0:par%dims%ntheta+1, 0:par%dims%nzlocal+1))
    allocate(pcg_auxarrays%tmp(0:par%dims%nr+1, 0:par%dims%ntheta+1, 0:par%dims%nzlocal+1))
    allocate(pcg_auxarrays%inv_diagA(1:par%dims%nr, 1:par%dims%ntheta, 1:par%dims%nzlocal))
  endif

  ! Arrays for capacitance and sensitivity, only needed on the finest grid.
  allocate(capacity(par%nel,par%nel),source=0.0_CUSTOM_REAL)
  allocate(cell_volume(par%dims%nr+1,par%dims%ntheta0,par%dims%nzlocal + (myrank+1)/nbproc),source=0.0_CUSTOM_REAL)
  allocate(phisensit(0:par%dims%nr+1,0:par%dims%ntheta0+1,0:par%dims%nzlocal+1,par%nel),source=0.0_CUSTOM_REAL,stat=ierr)

  ! Check that the allocation of the last large array went well,
  ! to see if the machine has enough memory per proc to run the problem under study.
  if(ierr /= 0) call exit_MPI("Dynamic memory allocation error, you probably ran out of memory; exiting...",myrank,09)

!=============================================================
!==== Done with all allocations. =============================
!=============================================================

  ! Define geometry, needed on every level.
  do ilevel=ilevel_fine, par%ilevel_coarse
    call set_geometry(ilevel,par,dims_at_level(ilevel)%nr,dims_at_level(ilevel)%ntheta,dims_at_level(ilevel)%nz,&
                      r(ilevel)%p,dr(ilevel)%p,theta(ilevel)%p,dz(ilevel)%p,&
                      z(ilevel)%p,izspace(ilevel),dims_at_level(ilevel)%kguards,&
                      i1(ilevel),i2(ilevel),elecs(ilevel)%p,guards(ilevel))

    if (.not. suppress_output .and. myrank == 0) print *,'GEOMETRY_INIT, level, index of R1 and R2, and izspace',ilevel,&
                            'i1=',i1(ilevel),'i2=',i2(ilevel),'izspace=',izspace(ilevel)
#ifdef USE_FLUSH6
    call flush(6)
#endif
  enddo

  ! Calculate cell volumes that also act as scaling factor for sensitivity kernel and model.
  call calculate_cell_volumes(par%dims,cell_volume,cell_volume_avg,&
                              r(ilevel_fine)%p,theta(ilevel_fine)%p,z(ilevel_fine)%p,myrank,nbproc)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Loop over nel electrodes = nel to solve forward problems.
  !! They all have the same matrix, if no mesh refinement used.
  !! But all have different RHS and initial guess for the solution.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (myrank == 0) print *, 'BEGIN ELECTRODE LOOP'

  do ielectrode = 1, par%nel

    if ((par%irefine == 1) .or. (ielectrode == 1)) then

      ! Perform geometry refinement near the firing electrode.
      if (par%irefine == 1) then
        call refine_geometry(ielectrode,dims_at_level(ilevel_fine)%ntheta,theta(ilevel_fine)%p,&
                             par%nel,elecs(ilevel_fine)%p,&
                             par%dims%ntheta0,theta0(ilevel_fine)%p,elecs0(ilevel_fine)%p,&
                             ithetamap(ilevel_fine)%p)
      endif

      ! Compute grid coordinates used for model initialization or visualization.
      do ilevel=ilevel_fine, par%ilevel_coarse
        call compute_grid_coordinates(dims_at_level(ilevel)%nr,dims_at_level(ilevel)%ntheta,dims_at_level(ilevel)%nz,&
                                      r(ilevel)%p,theta(ilevel)%p,z(ilevel)%p,&
                                      xgrid(ilevel)%p,ygrid(ilevel)%p,zgrid(ilevel)%p)
      enddo

      ! Compute boundary sites by initializing a flag array for all the FV cell centers, including ghost cells.
      ! This is needed on every level, and will be used to assemble the matrix, the RHS and the initial guess later.
      do ilevel=ilevel_fine, par%ilevel_coarse
        call flag_init_bc(flag(ilevel)%p,dims_at_level(ilevel)%nr,dims_at_level(ilevel)%ntheta,dims_at_level(ilevel)%nzlocal,&
                          i1(ilevel),i2(ilevel),par%nel,elecs(ilevel)%p,guards(ilevel),myrank,nbproc)
        if (.not. suppress_output .and. myrank == 0) print *, 'FLAG_INIT, level=', ilevel
      enddo

      if (par%read_model_from_inversion == 1) then
        ! NOTE: This is implemented only for the single-grid case (the full model from inversion used).
        ilevel = ilevel_fine

        ! Initialize the model from inversion.
        call init_model_from_inversion(par, dims_at_level(ilevel), permit(ilevel)%p, model, &
                                       i1(ilevel), i2(ilevel), myrank, nbproc)

        if (.not. suppress_output .and. myrank == 0) print *, 'END MODEL READING'
      else
        ! Generate a synthetic model.
        do ilevel = ilevel_fine, par%ilevel_coarse
          call model_init(par,permit(ilevel)%p,dims_at_level(ilevel),i1(ilevel),i2(ilevel),&
                          permit_bubble,xb,yb,zb,rb,r(ilevel_fine)%p,theta(ilevel_fine)%p,z(ilevel_fine)%p,myrank)

          if (.not. suppress_output .and. myrank == 0) print *, 'MODEL_INIT, level=', ilevel
        enddo
      endif

      ! Write the model to 1D array for output (for inversion).
      call write_model_to_1D_array(permit(ilevel_fine)%p, par%dims, model, myrank, nbproc)

#ifndef SUPPRESS_OUTPUT
      ! Write Paraview images of the model.
      call visualize_model(permit(ilevel_fine)%p, r(ilevel_fine)%p, theta(ilevel_fine)%p, z(ilevel_fine)%p, &
                           par%dims, ielectrode, myrank)
#endif

      ! Matrix construction (assembly), needed on every level.
      do ilevel = ilevel_fine, par%ilevel_coarse
        call matrix_init(a(ilevel)%p, flag(ilevel)%p, permit(ilevel)%p, &
                         par%permit_isolated_tube, par%permit_air, r(ilevel)%p, &
                         theta(ilevel)%p, dz(ilevel)%p, dims_at_level(ilevel), i1(ilevel), i2(ilevel), myrank)

        !call test_matrix_symmetry(a(ilevel)%p,nr_at_level(ilevel),ntheta_at_level(ilevel),&
        !                          nz_at_level(ilevel),nzlocal_at_level(ilevel),flag(ilevel)%p,myrank)
      enddo
      if (.not. suppress_output .and. myrank == 0) print *, 'END MATRIX INIT'

      ! Preconditioner / smoother construction, different for PCG and MG solver.
      if (par%linear_solver == LINSOLV_PCG) then
        ! inverse of diagonal of A
        do k=1,par%dims%nzlocal
          do j=1,par%dims%ntheta
            do i=1,par%dims%nr
              pcg_auxarrays%inv_diagA(i,j,k) = 1._CUSTOM_REAL / a(ilevel_fine)%p(i,j,k,3)
            enddo
          enddo
        enddo
      else if (par%linear_solver == LINSOLV_MG) then
        ! inverse of diagonal of A
        do ilevel=ilevel_fine, par%ilevel_coarse
          do k=1,dims_at_level(ilevel)%nzlocal
            do j=1,dims_at_level(ilevel)%ntheta
              do i=1,dims_at_level(ilevel)%nr
                !mg_auxarrays(ilevel)%inv_diagA(i,j,k) = 1._CUSTOM_REAL / a(ilevel)%p(i,j,k,3)
              enddo
            enddo
          enddo
        enddo
        ! coarse grid solver conjugate gradient uses Jacobi preconditioning as well
        if (par%coarse_solver == LINSOLV_PCG) then
          do k=1,dims_at_level(par%ilevel_coarse)%nzlocal
            do j=1,dims_at_level(par%ilevel_coarse)%ntheta
              do i=1,dims_at_level(par%ilevel_coarse)%nr
                pcg_auxarrays%inv_diagA(i,j,k) = 1._CUSTOM_REAL / a(par%ilevel_coarse)%p(i,j,k,3)
              enddo
            enddo
          enddo
        endif
      endif
    endif
    !-----------------------------------------------------------------------------------------

    ! Assign potentials val at electrode locations for ielectrode problem.
    ! Only needed on the finest grid.
    call solution_init(phi(ilevel_fine)%p,val(ilevel_fine)%p,flag(ilevel_fine)%p,x2(ilevel_fine)%p,par%dims,&
                       myrank,i2(ilevel_fine),elecs(ilevel_fine)%p(ielectrode))

    if (.not. suppress_output .and. myrank == 0) print *, 'END SOLUTION INIT'

    ! Source = right-hand-side term b of Ax = b equation, only needed on the finest grid.
    call source_RHS(b_RHS(ilevel_fine)%p,flag(ilevel_fine)%p,par%dims,val(ilevel_fine)%p)

    if (.not. suppress_output .and. myrank == 0) print *, 'END SOURCE_RHS'

    ! Solve linear system with matrix a, source b_RHS and potential solution x2.
    if (par%linear_solver == LINSOLV_PCG) then
      call solver_pcg(a(ilevel_fine)%p, b_RHS(ilevel_fine)%p, x2(ilevel_fine)%p, &
                      par%itypenorm, par%iprecond, par%omega1, par%tol, par%itmax, iter_pcg, &
                      par%output_frequency, .true., &
                      dims_at_level(ilevel_fine)%nr, dims_at_level(ilevel_fine)%ntheta, dims_at_level(ilevel_fine)%nzlocal, &
                      ierr, myrank, nbproc, pcg_auxarrays)
    else
      ! Copy Dirichlet values into initial guess before solving.
      do k=0,dims_at_level(ilevel_fine)%nzlocal+1
        do j=0,dims_at_level(ilevel_fine)%ntheta+1
          do i=0,dims_at_level(ilevel_fine)%nr+1
            if (flag(ilevel_fine)%p(i,j,k) == 1) x2(ilevel_fine)%p(i,j,k) = b_RHS(ilevel_fine)%p(i,j,k)
          enddo
        enddo
      enddo

!      call solver_mg(ilevel_fine,par%ilevel_coarse,a,b_RHS,x2,mg_auxarrays,pcg_auxarrays,flag,&
!                     dims_at_level%nr,dims_at_level%ntheta,dims_at_level%nzlocal,&
!                     par%tol,par%itmax,par%smoother,par%npresmooth,par%npostsmooth,par%omega,&
!                     par%coarse_solver,par%itmax_coarse,par%iprecond,par%omega1,par%tol_coarse,&
!                     par%output_frequency,ierr,myrank,nbproc)
    endif

    !call dump_vector_ascii(x2(ilevel_fine),nr_read,ntheta_read,nzlocal_read,myrank,"x2.txt")
    !call dump_vector_gmv(x2(ilevel_fine),nr_read,ntheta_read,nzlocal_read,myrank,linear_solver)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! From here on, there is no more level-dependent computation,
    !! everything is performed on the initial, finest grid.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Convert solution of the linear system back to potential phi.
    do k=0,par%dims%nzlocal+1
      do j=1,par%dims%ntheta
        do i=1,par%dims%nr
          phi(ilevel_fine)%p(i,j,k) = x2(ilevel_fine)%p(i,j,k)
        enddo
      enddo
    enddo

    ! Enforce periodic boundary conditions for the potential solution.
    ! DG DG this might be redundant, depends on the last operations the
    ! DG DG linear solver performed, e.g. if it finished with a matrix-vector multiplication
    do k=0,par%dims%nzlocal+1
      do i=0,par%dims%nr+1
        phi(ilevel_fine)%p(i,0,k) = phi(ilevel_fine)%p(i,par%dims%ntheta,k)
        phi(ilevel_fine)%p(i,par%dims%ntheta+1,k) = phi(ilevel_fine)%p(i,1,k)
      enddo
    enddo

#ifndef SUPPRESS_OUTPUT
    if (ielectrode == 1) then
      ! Output the electric potential data for visualization in Paraview.
      call visualize_solution(par%dims, phi(ilevel_fine)%p, &
                              xgrid(ilevel_fine)%p, ygrid(ilevel_fine)%p, zgrid(ilevel_fine)%p, &
                              ielectrode, myrank)
    endif
#endif

    ! Capacitance calculations for current electrode ielectrode, using surface integral approximation.
    call capacitance_computed_directly(phi(ilevel_fine)%p, theta(ilevel_fine)%p, r(ilevel_fine)%p, &
                                       z(ilevel_fine)%p, permit(ilevel_fine)%p, i1(ilevel_fine), i2(ilevel_fine), &
                                       par, capacity, ielectrode, elecs(ilevel_fine)%p, guards(ilevel_fine), myrank, nbproc)

    ! Potential calculations local to each process for current electrode problem ielectrode.
    ! Needed to compute the sensitivity kernel at the end of the routine.
    if (par%irefine == 1) then
      ! Map potential back on the original unrefined mesh -- use subset of nodes defined by ithetamap
      do k=0,par%dims%nzlocal+1
        do j=0,par%dims%ntheta0+1
          do i=0,par%dims%nr+1
            phisensit(i,j,k,ielectrode) = phi(ilevel_fine)%p(i,ithetamap(ilevel_fine)%p(j),k)
          enddo
        enddo
      enddo

#ifndef SUPPRESS_OUTPUT
      ! (+) For visualization on the original mesh ------------------------------------------------------------
      ! TODO: Move to a separate function.
!      do k=0,nzlocal_read+1
!        do j=0,ntheta_read0+1
!          do i=0,nr_read+1
!            xgrid(ilevel_fine)%p(i,j,k) = xgrid(ilevel_fine)%p(i,ithetamap(ilevel_fine)%p(j),k)
!            ygrid(ilevel_fine)%p(i,j,k) = ygrid(ilevel_fine)%p(i,ithetamap(ilevel_fine)%p(j),k)
!            zgrid(ilevel_fine)%p(i,j,k) = zgrid(ilevel_fine)%p(i,ithetamap(ilevel_fine)%p(j),k)
!          enddo
!        enddo
!      enddo
!
!      ! Back to original theta.
!      theta(ilevel_fine)%p(:) = theta0(ilevel_fine)%p(:)
!
!      call paraview_write_2d_profiles(myrank, "phi2d_", ielectrode, &
!          nr_at_level(ilevel_fine), ntheta_read0, nz_at_level(ilevel_fine), nzlocal_at_level(ilevel_fine), &
!          phisensit(:,:,:,ielectrode), xgrid(ilevel_fine)%p, ygrid(ilevel_fine)%p, zgrid(ilevel_fine)%p)
      ! (-) For visualization on the original mesh  ------------------------------------------------------------
#endif

    else
    ! No mesh refinement.
      do k=0,par%dims%nzlocal+1
        do j=0,par%dims%ntheta+1
          do i=0,par%dims%nr+1
            phisensit(i,j,k,ielectrode) = phi(ilevel_fine)%p(i,j,k)
          enddo
        enddo
      enddo
    endif

  enddo ! end of external loop on the nel electrodes, i.e., end of solving the nel forward problems

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Done with solving all forward problems.
  !! Compute sensitivity kernel now.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (.not. suppress_output .and. myrank == 0) call write_data_to_file(par, capacity, 'capa-at-elec', myrank)

  ! Write 1D data array (Cij, i>j) for inversion.
  call write_data_for_inversion(par, capacity, cdata)

  !-----------------------------------------------------------------------------
  if (par%irefine == 1) then
    ! Back to the original unrefined mesh.
    theta(ilevel_fine)%p = theta0(ilevel_fine)%p
  endif

  ! Calculations of sensitivity kernel,
  ! with approximation of pixels in each elementary volume.
  ! Also calculates the scaling factor for inversion.
  if (myrank == 0) print *, 'Calculating ECT sensitivity...'
  call calculate_sensitivity(par, sensit, phisensit, cell_volume, &
                             r(ilevel_fine)%p, theta(ilevel_fine)%p, z(ilevel_fine)%p, myrank, nbproc)

  ! Calculates weights needed for inversion.
  call calculate_weights(par, cell_volume, column_weight, damping_weight, myrank, nbproc)

#ifndef SUPPRESS_OUTPUT
  ! Calculate capacitance based on the sensitivity.
  ! NOTE: We do not need these values for calculations, but to compare with capacitance from electrodes, to show accuracy.
  call calculate_capa_sens(par, sensit, capacity, permit(ilevel_fine)%p, myrank, nbproc)
#endif

#ifndef SUPPRESS_OUTPUT
  ! Write Paraview profiles for the sensitivity.
  call visualise_sensitivity(par, sensit, column_weight, &
                             r(ilevel_fine)%p, theta(ilevel_fine)%p, z(ilevel_fine)%p, myrank, nbproc)
#endif

  !-----------------------------------------------------------------------------
  ! Deallocate all the arrays.
  do ilevel = ilevel_fine, par%ilevel_coarse

    deallocate(theta(ilevel)%p)
    deallocate(dr(ilevel)%p)
    deallocate(dz(ilevel)%p)
    deallocate(r(ilevel)%p)
    deallocate(z(ilevel)%p)
    deallocate(xgrid(ilevel)%p)
    deallocate(ygrid(ilevel)%p)
    deallocate(zgrid(ilevel)%p)
    deallocate(flag(ilevel)%p)
    deallocate(val(ilevel)%p)
    deallocate(phi(ilevel)%p)
    deallocate(a(ilevel)%p)
    deallocate(b_RHS(ilevel)%p)
    deallocate(x2(ilevel)%p)
    deallocate(permit(ilevel)%p)
    deallocate(elecs(ilevel)%p)

    if (par%irefine == 1) then
      deallocate(theta0(ilevel)%p)
      deallocate(ithetamap(ilevel)%p)
      deallocate(elecs0(ilevel)%p)
    endif

!    if (par%linear_solver == LINSOLV_MG) then
!      deallocate(mg_auxarrays(ilevel)%c)
!      deallocate(mg_auxarrays(ilevel)%d)
!      deallocate(mg_auxarrays(ilevel)%inv_diagA)
!      deallocate(mg_auxarrays(ilevel)%cprime_r)
!      deallocate(mg_auxarrays(ilevel)%dprime_r)
!      deallocate(mg_auxarrays(ilevel)%xprime_r)
!      deallocate(mg_auxarrays(ilevel)%cprime_theta)
!      deallocate(mg_auxarrays(ilevel)%dprime_theta)
!      deallocate(mg_auxarrays(ilevel)%xprime_theta)
!      deallocate(mg_auxarrays(ilevel)%cprime_z)
!      deallocate(mg_auxarrays(ilevel)%dprime_z)
!      deallocate(mg_auxarrays(ilevel)%xprime_z)
!    endif
  enddo

  ! Deallocate all pointer arrays.
  deallocate(izspace)
  deallocate(theta)
  deallocate(dr)
  deallocate(dz)
  deallocate(r)
  deallocate(z)
  deallocate(xgrid)
  deallocate(ygrid)
  deallocate(zgrid)
  deallocate(flag)
  deallocate(i1)
  deallocate(i2)
  deallocate(val)
  deallocate(phi)
  deallocate(a)
  deallocate(b_RHS)
  deallocate(x2)
  deallocate(permit)
  deallocate(elecs)
  deallocate(guards)

  if (par%irefine == 1) then
    deallocate(theta0)
    deallocate(ithetamap)
    deallocate(elecs0)
  endif

  if (par%linear_solver == LINSOLV_MG) then
    !deallocate(mg_auxarrays)
  endif
  deallocate(dims_at_level)

  ! Deallocate all non-level-dependent data.
  deallocate(capacity)
  deallocate(phisensit)
  deallocate(cell_volume)

  if (par%num_bubbles > 0) then
    deallocate(xb)
    deallocate(yb)
    deallocate(zb)
    deallocate(rb)
    deallocate(permit_bubble)
  endif

  if (par%linear_solver == LINSOLV_PCG) then
    deallocate(pcg_auxarrays%p)
    deallocate(pcg_auxarrays%r)
    deallocate(pcg_auxarrays%inv_diagA)
    deallocate(pcg_auxarrays%z)
  endif

  if (par%linear_solver == LINSOLV_MG .and. par%coarse_solver == LINSOLV_PCG) then
    deallocate(pcg_auxarrays%p)
    deallocate(pcg_auxarrays%r)
    deallocate(pcg_auxarrays%inv_diagA)
    deallocate(pcg_auxarrays%z)
  endif

end subroutine solve_forward_problem_ect

end module forward_problem_ect

