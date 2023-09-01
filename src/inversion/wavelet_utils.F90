
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

!=============================================================================
! Utils to use wavelet compression with the parallel model vector.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!=============================================================================
module wavelet_utils

  use global_typedefs
  use parallel_tools
  use wavelet_transform

  implicit none

  private

  public :: apply_wavelet_transform

contains

!======================================================================================
! Apply forward/inverse transfrom to the parallel model vector.
!======================================================================================
subroutine apply_wavelet_transform(nelements, nx, ny, nz, ncomponents, v, model_full, FWD, &
                                   compression_type, nproblems, SOLVE_PROBLEM, myrank, nbproc)
  integer, intent(in) :: nelements, nx, ny, nz, ncomponents
  logical, intent(in) :: FWD
  integer, intent(in) :: compression_type
  integer, intent(in) :: nproblems
  logical, intent(in) :: SOLVE_PROBLEM(nproblems)
  integer, intent(in) :: myrank, nbproc

  ! Buffer for wavelet transform.
  real(kind=CUSTOM_REAL), intent(inout) :: model_full(nx * ny * nz)
  ! Converted to/from wavelet domain result.
  real(kind=CUSTOM_REAL), intent(inout) :: v(nelements, ncomponents, nproblems)

  integer :: i, k

  do i = 1, nproblems
    if (SOLVE_PROBLEM(i)) then
      ! Loop over the model components.
      do k = 1, ncomponents
        call get_full_array(v(:, k, i), nelements, model_full, .false., myrank, nbproc)

        if (myrank == 0) then
          if (FWD) then
            call forward_wavelet(model_full, nx, ny, nz, compression_type)
          else
            call inverse_wavelet(model_full, nx, ny, nz, compression_type)
          endif
        endif

        call scatter_full_array(nelements, model_full, v(:, k, i), myrank, nbproc)
      enddo
    endif
  enddo

end subroutine apply_wavelet_transform

end module wavelet_utils
