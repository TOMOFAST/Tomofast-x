
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

module forward_problem_gravmag

  use global_typedefs, only: CUSTOM_REAL
  use data_gravmag
  use sensitivity_gravmag
  use parameters_gravmag
  use inversion_arrays

  implicit none

  private

  public :: solve_forward_problem
  public :: calculate_kernel_size_and_weight

contains

!=========================================================================================
! Solves forward problem for gravity or magnetism.
!=========================================================================================
subroutine solve_forward_problem(par, iarr, data, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  type(t_data), intent(inout) :: data
  type(t_inversion_arrays), intent(inout) :: iarr
  integer, intent(in) :: myrank, nbproc

  type(t_sensitivity_gravmag) :: sens

  ! Calculate sensitivity kernel (analytically).
  call sens%calculate_sensitivity_kernel(par, iarr%model%grid, data, iarr%column_weight, iarr%matrix_sensit, myrank, nbproc)

  ! Calculate the data using sensitivity (S) and prior model (m) as d = S * m.
  call iarr%model%calculate_data(par%ndata, iarr%matrix_sensit, iarr%column_weight, data%val_calc, par%compression_type, &
                                 myrank, nbproc)

end subroutine solve_forward_problem

!=========================================================================================
! Calcualtes the size of the sensitivity kernel.
!=========================================================================================
subroutine calculate_kernel_size_and_weight(par, iarr, data, nnz, weight_multiplier, myrank, nbproc)
  class(t_parameters_base), intent(in) :: par
  type(t_data), intent(inout) :: data
  type(t_inversion_arrays), intent(inout) :: iarr
  real(kind=CUSTOM_REAL) :: weight_multiplier
  integer, intent(in) :: myrank, nbproc
  integer, intent(out) :: nnz

  type(t_sensitivity_gravmag) :: sens

  nnz = sens%predict_sensit_kernel_size(par, iarr%model%grid, data, iarr%column_weight, myrank, nbproc)

  if (par%depth_weighting_type == 3) then
    ! Precondition the column weights.
    iarr%column_weight = weight_multiplier * iarr%column_weight
  endif

end subroutine calculate_kernel_size_and_weight

end module forward_problem_gravmag
