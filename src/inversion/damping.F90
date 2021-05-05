
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

!===========================================================================================
! A class to add (inversion) damping contribution to the matrix and right hand side.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!===========================================================================================
module damping

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use sparse_matrix
  use parallel_tools
  use model
  use grid

  implicit none

  private

  type, public :: t_damping
    private

    ! Number of model parameters.
    integer :: nelements
    ! Weight of the damping.
    real(kind=CUSTOM_REAL) :: alpha
    ! Weight of the whole problem (damping + misfit) in joint inversion.
    real(kind=CUSTOM_REAL) :: problem_weight
    ! Power p of Lp norm (for LSQR method). Use p=2. for pure LSQR.
    real(kind=CUSTOM_REAL) :: norm_power
    ! Cost of the model objective function phi_m.
    real(kind=CUSTOM_REAL) :: cost

  contains
    private

    procedure, public, pass :: initialize => damping_initialize
    procedure, public, pass :: add => damping_add
    procedure, public, pass :: get_cost => damping_get_cost

    procedure, private, pass :: add_RHS => damping_add_RHS
    procedure, private, pass :: get_norm_multiplier => damping_get_norm_multiplier

  end type t_damping

contains

!===========================================================================================
! Initialization.
!===========================================================================================
subroutine damping_initialize(this, nelements, alpha, problem_weight, norm_power, myrank)
  class(t_damping), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: alpha, problem_weight, norm_power
  integer, intent(in) :: nelements, myrank

  this%nelements = nelements
  this%alpha = alpha
  this%problem_weight = problem_weight
  this%norm_power = norm_power
end subroutine damping_initialize

!===========================================================================================
! 1) Adds damping (below the sensitivity kernel) which is identity matrix times alpha:
!     S_new = (  S )
!             ( wD ),
!    where D is damping matrix, w is damping parameter,
!    applying the same scaling as applied to the sensitivity matrix.
! 2) Adds the corresponding contribution to the right hand side 'b_RHS'.
!
! Tested in unit_tests.f90 in test_damping_identity_matrix().
!===========================================================================================
subroutine damping_add(this, matrix, b_RHS, column_weight, damping_weight, &
                       model, model_ref, param_shift, myrank, nbproc)
  class(t_damping), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  real(kind=CUSTOM_REAL), intent(in) :: damping_weight(:)
  type(t_model), intent(in) :: model
  real(kind=CUSTOM_REAL), intent(in) :: model_ref(:)
  integer, intent(in) :: param_shift
  integer, intent(in) :: myrank, nbproc

  type(t_sparse_matrix), intent(inout) :: matrix
  real(kind=CUSTOM_REAL), intent(inout) :: b_RHS(:)

  integer :: i, nsmaller, nelements_total
  integer :: row_beg, row_end
  real(kind=CUSTOM_REAL) :: value, alpha
  type(t_parallel_tools) :: pt

  ! The number of elements on CPUs with rank smaller than myrank.
  nsmaller = pt%get_nsmaller(this%nelements, myrank, nbproc)

  ! The total number of elements.
  nelements_total = pt%get_total_number_elements(this%nelements, myrank, nbproc)

  ! First matrix row (in the big matrix) of the damping matrix that will be added.
  row_beg = matrix%get_current_row_number() + 1

  ! Add empty lines.
  call matrix%add_empty_rows(nsmaller, myrank)

  ! Add lines with damping.
  do i = 1, this%nelements
    call matrix%new_row(myrank)

    ! Map only the model perturbation variable.
    value = this%alpha * this%problem_weight * damping_weight(i) * &
            column_weight(i) * &
            this%get_norm_multiplier(model%val(i), model_ref(i))

    ! Apply model covariance (diagonal), which is equivalent of having local alpha.
    value = value * model%cov(i)

    call matrix%add(value, param_shift + i, myrank)
  enddo

  ! Add empty lines.
  call matrix%add_empty_rows(nelements_total - this%nelements - nsmaller, myrank)

  ! Last matrix row (in the big matrix) of the added damping matrix.
  row_end = matrix%get_current_row_number()

  ! Sanity check.
  if (row_end - row_beg + 1 /= nelements_total) &
    call exit_MPI("Sanity check failed in damping_add!", myrank, 0)

  ! Add the damping contribution to the right hand side.
  call this%add_RHS(b_RHS(row_beg:row_end), model, model_ref, damping_weight, myrank, nbproc)

  ! Calculate the damping cost.
  this%cost = sum(b_RHS(row_beg:row_end)**2)

end subroutine damping_add

!=============================================================================================
! Adds damping contribution in the right hand side.
! model_ref - reference model.
!=============================================================================================
subroutine damping_add_RHS(this, b_RHS, model, model_ref, damping_weight, myrank, nbproc)
  class(t_damping), intent(in) :: this
  type(t_model), intent(in) :: model
  real(kind=CUSTOM_REAL), intent(in) :: model_ref(:)
  real(kind=CUSTOM_REAL), intent(in) :: damping_weight(:)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), intent(inout) :: b_RHS(:)
  type(t_parallel_tools) :: pt
  real(kind=CUSTOM_REAL) :: alpha
  integer :: i

  do i = 1, this%nelements
    b_RHS(i) = - this%alpha * this%problem_weight * damping_weight(i) * &
               (model%val(i) - model_ref(i)) * &
               this%get_norm_multiplier(model%val(i), model_ref(i))

    ! Apply model covariance (diagonal), which is equivalent of having local alpha.
    b_RHS(i) = b_RHS(i) * model%cov(i)
  enddo

  ! Gather full right hand side.
  call pt%get_full_array_in_place(this%nelements, b_RHS, .true., myrank, nbproc)

end subroutine damping_add_RHS

!===========================================================================================
! Returns model objective function cost (norm).
!===========================================================================================
pure function damping_get_cost(this) result(res)
  class(t_damping), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = this%cost

end function damping_get_cost

!===========================================================================================
! Returns a multiplier (for one pixel) to change L2 norm to Lp, in the LSQR method.
!===========================================================================================
pure function damping_get_norm_multiplier(this, model, model_ref) result(res)
  class(t_damping), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: model, model_ref
  real(kind=CUSTOM_REAL) :: res

  if (model /= model_ref) then
    res = (abs(model - model_ref))**(this%norm_power / 2.d0 - 1.d0)
  else
    res = 1._CUSTOM_REAL
  endif

end function damping_get_norm_multiplier

end module damping
