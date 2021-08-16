
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

contains

!=========================================================================================
! Solves forward problem for gravity or magnetism.
!=========================================================================================
subroutine solve_forward_problem(par, iarr, data, myrank)
  class(t_parameters_base), intent(in) :: par
  type(t_data), intent(inout) :: data
  type(t_inversion_arrays), intent(inout) :: iarr
  integer, intent(in) :: myrank

  type(t_sensitivity_gravmag) :: sens

  if (par%calc_data_directly == 1) then
    ! Calculate data directly without computing sensitivity matrix.
    call sens%calc_data_directly(par, iarr%model, data, myrank)

    return
  endif

  ! Calculate sensitivity kernel (analytically).
  call sens%calculate_sensitivity(par, iarr%model%grid, data, iarr%matrix_sensit, myrank)

  ! Calculate the data using sensitivity (S) and prior model (m) as d = S * m.
  call iarr%model%calculate_data(par%ndata, iarr%matrix_sensit, data%val_calc, myrank)

end subroutine solve_forward_problem

end module forward_problem_gravmag
