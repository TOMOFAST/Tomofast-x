
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

!==========================================================================
! A class to compare the models.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2017.
!==========================================================================
module compare_models

  use global_typedefs, only: CUSTOM_REAL
  use mpi_tools, only: exit_MPI
  use vector
  use model
  use parallel_tools
  use gradient

  implicit none

  private

  type, public :: t_compare_models
    private


  contains
    private

    procedure, public, nopass :: compare => compare_models_compare

  end type t_compare_models

contains

!=====================================================================================================
! Compare two models (final0 and final) by calculating the sum of cross gradient norms.
!=====================================================================================================
subroutine compare_models_compare(model, der_type, res, myrank, nbproc)
  class(t_model), intent(in) :: model
  integer, intent(in) :: der_type
  integer, intent(in) :: myrank, nbproc

  type(t_parallel_tools) :: pt
  type(t_vector) :: m1_grad, m2_grad, cross_grad
  type(t_gradient) :: grad
  integer :: i, j, k, p
  integer :: ierr

  real(kind=CUSTOM_REAL), intent(out) :: res

  real(kind=CUSTOM_REAL), allocatable :: model_val_full0(:)

  ierr = 0

  ! Allocate here an array for the full model (final0).
  ! TODO: For memory efficiency it is better to work with local models and pass border elements to calculate the gradients.
  allocate(model_val_full0(model%nelements_total), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in cross_gradient_compare_models!", myrank, ierr)

  ! MPI gather the full model.
  call pt%get_full_array(model%val_final0, model%nelements, model_val_full0, .true., myrank, nbproc)

  res = 0.d0

  do p = 1, model%nelements_total
    ! Grid (i,j,k)-index in the full grid.
    i = model%grid_full%i_(p)
    j = model%grid_full%j_(p)
    k = model%grid_full%k_(p)

    ! Calculate model gradients.
    m1_grad = grad%get_grad(model_val_full0, model%grid_full, i, j, k, grad%get_der_type(der_type))
    m2_grad = grad%get_grad(model%val_full, model%grid_full, i, j, k, grad%get_der_type(der_type))

    ! Calculate the cross-product between model gradients.
    cross_grad = m1_grad%cross_product(m2_grad)

    res = res + cross_grad%get_norm()
  enddo

  res = res / dble(model%nelements_total)

  deallocate(model_val_full0)

end subroutine compare_models_compare

end module compare_models
