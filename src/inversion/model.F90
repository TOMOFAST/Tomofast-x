
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

!================================================================================================
! A class that extends t_model_IO to work with parallel models (split between CPUs) for inversion.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!================================================================================================
module model

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use paraview
  use grid
  use string, only: str
  use parallel_tools
  use string
  use model_base
  use model_IO
  use sparse_matrix
  use wavelet_transform

  implicit none

  private

  type, extends(t_model_IO), public :: t_model

  contains
    private

    procedure, public, pass :: update => model_update
    procedure, public, pass :: calculate_data => model_calculate_data

  end type t_model

  public :: rescale_model

contains

!======================================================================================================
! Update model after inversion.
!======================================================================================================
subroutine model_update(this, delta_model)
  class(t_model), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: delta_model(:)

  integer :: i

  do i = 1, this%nelements
    this%val(i) = this%val(i) + delta_model(i)
  enddo

end subroutine model_update

!======================================================================================================
! Calculate the (linear) data using original (not scaled)
! sensitivity kernel (S) and model (m) as d = S * m.
!======================================================================================================
subroutine model_calculate_data(this, ndata, matrix_sensit, column_weight, data, compression_type, myrank, nbproc)
  class(t_model), intent(in) :: this
  integer, intent(in) :: ndata, compression_type, myrank, nbproc
  type(t_sparse_matrix), intent(in) :: matrix_sensit
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)

  real(kind=CUSTOM_REAL), intent(out) :: data(:)

  real(kind=CUSTOM_REAL), allocatable :: model_scaled(:)
  real(kind=CUSTOM_REAL), allocatable :: model_scaled_full(:)
  integer :: i, ierr
  integer :: nsmaller
  type(t_parallel_tools) :: pt

  allocate(model_scaled(this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_calculate_data!", myrank, ierr)

  ! Rescale the model to calculate data, as we store the depth-weighted sensitivity kernel.
  do i = 1, this%nelements
    model_scaled(i) = this%val(i) / column_weight(i)
  enddo

  if (compression_type > 0) then
  ! Apply wavelet transform to the model to calculate data using compressed sensitivity kernel.

    if (nbproc > 1) then
    ! Parallel version.

      ! Allocate memory for the full model.
      allocate(model_scaled_full(this%nelements_total), source=0._CUSTOM_REAL, stat=ierr)
      if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_calculate_data!", myrank, ierr)

      ! Gather the full model from all processors.
      call pt%get_full_array(model_scaled, this%nelements, model_scaled_full, .true., myrank, nbproc)

      ! Compress the full model.
      call Haar3D(model_scaled_full, this%grid%nx, this%grid%ny, this%grid%nz)

      ! Extract the local model part.
      nsmaller = pt%get_nsmaller(this%nelements, myrank, nbproc)
      model_scaled = model_scaled_full(nsmaller + 1 : nsmaller + this%nelements)

      deallocate(model_scaled_full)

    else
    ! Serial version.
      call Haar3D(model_scaled, this%grid%nx, this%grid%ny, this%grid%nz)
    endif
  endif

  ! Calculate data: d = S * m
  call matrix_sensit%mult_vector(model_scaled, data)

  deallocate(model_scaled)

  ! NOTE: Not sure this is correct way of calling the function (manual says MPI_IN_PLACE should be used on the root only).
  !       But this works somehow.
  call MPI_Allreduce(MPI_IN_PLACE, data, ndata, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(data, ndata, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI error in model_calculate_data!", myrank, ierr)

end subroutine model_calculate_data

!================================================================================================
! Weights the model parameters.
!================================================================================================
subroutine rescale_model(model, weight, nc)
  real(kind=CUSTOM_REAL), intent(inout) :: model(:)
  integer, intent(in) :: nc
  real(kind=CUSTOM_REAL), intent(in) :: weight(:)

  integer :: i

  do i = 1, nc
    model(i) = model(i) * weight(i)
  enddo
end subroutine rescale_model

end module model

