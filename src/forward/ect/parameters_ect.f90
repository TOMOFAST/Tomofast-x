
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

module parameters_ect

  use global_typedefs

  implicit none

  ! Type of linear solver.
  integer, parameter :: LINSOLV_PCG = 1
  integer, parameter :: LINSOLV_MG = 2

  !-------------------------------------------------------------------
  ! Dimensions of the problem.
  type, public :: t_dimensions
    ! r,theta and z direction
    integer :: nr
    integer :: ntheta
    integer :: nz
    ! The model is split in the z direction for MPI.
    integer :: nzlocal
    ! Original ntheta before mesh refinement.
    integer :: ntheta0
    ! TODO: get rid of it -- set it via specifying the physical length of the guards.
    ! The number of points of the electrical guards on which potential is set to zero.
    ! Two rings of guards are located one at the bottom and another one on top,
    ! they are just a little bit less than 5 cm each and separated
    ! from the electrodes by space_elec_guards = 0.002 cm
    integer :: kguards
  end type t_dimensions

  !-------------------------------------------------------------------
  ! Geometry parameters of the ECT sensor.
  type, public :: t_sensor
    ! R1, radius of inner domain
    real(kind=CUSTOM_REAL) :: radiusin
    ! R2, radius of electrodes ring
    real(kind=CUSTOM_REAL) :: radiusout
    ! R3, outer radius
    real(kind=CUSTOM_REAL) :: radiusoutout
    ! height of the cylindrical sensor
    real(kind=CUSTOM_REAL) :: heicyl
    ! theta-space between electrodes
    real(kind=CUSTOM_REAL) :: space_electrodes
    ! z-space between guards and electrodes
    real(kind=CUSTOM_REAL) :: space_elec_guards
  end type t_sensor

  !-------------------------------------------------------------------
  ! Grid parameters of the electrode.
  type, public :: t_electrode
    ! Minimal and maximal grid nodes in theta-direction.
    integer :: ithetamin, ithetamax

    ! Minimal and maximal grid nodes in z-direction.
    integer :: izmin, izmax
  end type t_electrode

  ! Array of pointers to represent parameters of electrodes at different multigrid level.
  ! We need this because we have an array of electrodes (e.g. 12) at different multigrid level.
  type, public :: t_electrode_p
    type(t_electrode), dimension(:), pointer :: p
  end type t_electrode_p

  !-------------------------------------------------------------------
  ! Grid parameters of the guards.
  type, public :: t_guards
    ! The maximal z-coordinate of the lower guard.
    integer :: lower_izmax
    ! The minimal z-coordinate of the upper guard.
    integer :: upper_izmin
  end type t_guards

  !-------------------------------------------------------------------
  ! Main ECT input parameters read from the Parfile.
  type, public :: t_parameters_ect
    ! Dimensions of the problem.
    type(t_dimensions) :: dims

    ! Number of electrodes.
    integer :: nel
    ! Number of electrode rings.
    integer :: nrings
    ! Flag: fixed electrodes by geometry.
    integer :: ifixed_elecgeo
    ! Flag: mesh refinement (NO = 0, YES = 1).
    integer :: irefine
    ! Geometry parameters of the ECT sensor.
    type (t_sensor) :: sens

    ! Number of bubbles to create in the physical model
    ! num_bubbles == 0 means no bubbles at all and then it takes permit_matrix as
    ! default permittivities inside the sensor (r < R1).
    integer :: num_bubbles
    ! File name to read the bubble model from (in case num_bubbles > 0)
    character(len=256) :: filename_bubbles
    ! absolute permittivity
    real(kind=CUSTOM_REAL) :: permit0
    ! permittivity in the inner tube when no bubbles.
    real(kind=CUSTOM_REAL) :: permit_matrix
    ! permittivity of oil
    real(kind=CUSTOM_REAL) :: permit_oil
    ! permittivity of air for r>R2  and permittivity of the isolated tube for R1<r<R2
    real(kind=CUSTOM_REAL) :: permit_isolated_tube, permit_air

    ! Flag: if we use the model from input array (from inversion) or generate synthetic model.
    integer :: read_model_from_inversion

    ! type of the linear solver (LINSOLVE_PCG or LINSOLVE_MG)
    integer :: linear_solver
    ! damping parameter for PCG preconditioning
    integer :: iprecond
    ! damping parameter for PCG preconditioning
    real(kind=CUSTOM_REAL) :: omega1
    ! max number of iterations, type of norm to use, frequency at which
    ! convergence is checked and printed out (the latter only for PCG)
    integer itmax,itypenorm,output_frequency
    ! tolerance error on the residual for the linear solver
    real(kind=CUSTOM_REAL) :: tol

    ! number of multigrid levels
    integer :: ilevel_coarse
    integer :: coarse_solver

    ! read initial guess phi0 from file?
    logical :: read_guess_from_file

  contains
    procedure, pass :: broadcast => t_parameters_ect_broadcast
    procedure, pass :: accept_electrode_pair => t_parameters_ect_accept_electrode_pair
    procedure, pass :: get_ndata => t_parameters_ect_get_ndata
    procedure, pass :: get_electrode_index => t_parameters_ect_get_electrode_index
  end type t_parameters_ect

contains

!====================================================================================
! MPI broadcast parameters that are read from a Parfile.
!====================================================================================
subroutine t_parameters_ect_broadcast(this)
  class(t_parameters_ect), intent(in) :: this
  integer :: ierr

  call MPI_Bcast(this%read_guess_from_file,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%sens%radiusin,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%sens%radiusout,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%sens%radiusoutout,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%sens%heicyl,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%sens%space_elec_guards,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%sens%space_electrodes,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%ifixed_elecgeo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%iprecond,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%omega1,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%irefine,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%permit0,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%permit_air,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%permit_isolated_tube,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%permit_matrix,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%permit_oil,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%linear_solver,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%itypenorm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%tol,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%itmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%output_frequency,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%dims%nr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%dims%ntheta,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%dims%nz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%dims%kguards,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%nel,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%nrings,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(this%num_bubbles,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  call MPI_Bcast(this%ilevel_coarse,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  !call MPI_Bcast(this%smoother,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  !call MPI_Bcast(this%npresmooth,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  !call MPI_Bcast(this%npostsmooth,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  !call MPI_Bcast(this%omega,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)
  !call MPI_Bcast(this%coarse_solver,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  !call MPI_Bcast(this%itmax_coarse,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  !call MPI_Bcast(this%tol_coarse,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ierr)

end subroutine t_parameters_ect_broadcast

!==================================================================================================
! Determines if the pair of electrode is used for inversion.
!==================================================================================================
pure function t_parameters_ect_accept_electrode_pair(this, ielectrode, jelectrode) result(accept)
  class(t_parameters_ect), intent(in) :: this
  integer, intent(in) :: ielectrode, jelectrode
  logical :: accept

  ! Number of electrodes per ring.
  integer :: nelr, ielec1, jelec1

  nelr = this%nel / this%nrings

  accept = .true.

  ! NOTE: Uncomment this to accept all electrode pairs.
  !return

  ! Exclude symmetric pairs (e.g. 1-4 and 4-1) and self-capacitances (e.g. 1-1, 2-2 etc).
  if (ielectrode > jelectrode) then

    ! Exclude adjacent electrodes.
!    if (ielectrode - jelectrode == 1) accept = .false.        ! e.g. 2-1
!    if (ielectrode - jelectrode == nelr - 1) accept = .false. ! e.g. 12-1
!    if (ielectrode - jelectrode == nelr) accept = .false.     ! e.g. 13-1
!    if (ielectrode - jelectrode == nelr + 1) accept = .false. ! e.g. 14-1
!    if (ielectrode - jelectrode == 2 * nelr - 1) accept = .false. ! e.g. 24-1

    ! The code below excludes 'mapped' adjacent electrodes, e.g.,
    ! electrodes 1 and 25 (from 1st and 3rd rings) are considered adjacent,
    ! since the are located at the same theta-positions.
    ! Note, for 1- and 2-ring cases it is the same as the code above.

    ! Map electrodes to the first ring (e.g., 13 --> 1, 24 --> 12 etc).
    ielec1 = mod(ielectrode, nelr)
    if (ielec1 == 0) ielec1 = nelr

    jelec1 = mod(jelectrode, nelr)
    if (jelec1 == 0) jelec1 = nelr

    if (abs(ielec1 - jelec1) <= 1 &
        .or. abs(ielec1 - jelec1) >= nelr - 1) accept = .false.

  else
    accept = .false.
  endif

end function t_parameters_ect_accept_electrode_pair

!======================================================================================
! Returns the number of data using the rule of electrode pair acceptance,
! defined in accept_electrode_pair().
!======================================================================================
pure function t_parameters_ect_get_ndata(this) result(ndata)
  class(t_parameters_ect), intent(in) :: this
  integer :: ndata
  integer :: ielectrode, jelectrode

  ndata = 0

  do jelectrode = 1, this%nel
    do ielectrode = 1, this%nel
      if (this%accept_electrode_pair(ielectrode, jelectrode)) then
        ndata = ndata + 1
      endif
    enddo
  enddo

end function t_parameters_ect_get_ndata

!==================================================================================================
! Returns electrode 2D index (pair(i,j)) for a given 1D linear data_index.
! if elec_type = 1 then returns index of 1st electrode,
! if elec_type = 2 then returns index of 2nd electrode.
!==================================================================================================
pure function t_parameters_ect_get_electrode_index(this, data_index, elec_type) result(elec_index)
  class(t_parameters_ect), intent(in) :: this
  integer, intent(in) :: data_index, elec_type
  integer :: elec_index
  integer :: ielectrode, jelectrode, l

  elec_index = 0
  l = 0

  do jelectrode=1,this%nel
    do ielectrode=1,this%nel
      if (this%accept_electrode_pair(ielectrode,jelectrode)) then
        l = l+1

        if (l == data_index) then
          if (elec_type == 1) then
            elec_index = ielectrode
          else if (elec_type == 2) then
            elec_index = jelectrode
          endif
          return
        endif

      endif
    enddo
  enddo

end function t_parameters_ect_get_electrode_index

end module parameters_ect
