
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
! A class to smooth the initial model.
! This is needed e.g. for cross-gradient joint inversion to have non-zero gradients.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!================================================================================================
module filter

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use model

  implicit none

  private

  type, public :: t_filter
    private

  contains
    private

    procedure, public, nopass :: apply_Gaussian_filter => apply_Gaussian_filter

  end type t_filter

contains

!======================================================================================================
! Applies Gaussian filter to the model.
! sigma - standard deviation.
! lfgauss - size of the Gaussian window.
! type: 1 - arithmetic average, 2 - log-average (same as geometric mean).
!======================================================================================================
subroutine apply_Gaussian_filter(model, lfgauss, sigma, type, myrank)
  integer, intent(in) :: lfgauss
  real(kind=CUSTOM_REAL), intent(in) :: sigma
  integer, intent(in) :: type, myrank
  type(t_model), intent(inout) :: model

  integer :: i_line
  integer :: i_line_1
  integer :: i, j, k
  integer :: i1, j1, k1
  integer :: nbx, nby, nbz
  integer :: ierr
  real(kind=CUSTOM_REAL) :: lx, ly, lz, dx, dy, dz, avg_cell_size
  real(kind=CUSTOM_REAL) :: dx1_square, dy1_square, dz1_square
  real(kind=CUSTOM_REAL) :: gauss, gauss_const, scale

  real(kind=CUSTOM_REAL), allocatable :: model_temp(:)

  allocate(model_temp(model%nelements_total), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in apply_Gaussian_filter!", myrank, ierr)

  gauss_const = 1.d0 / (2.d0 * PI * sigma**2)

  dx = model%grid_full%get_hx()
  dy = model%grid_full%get_hy()
  dz = model%grid_full%get_hz()

  avg_cell_size = (dx + dy + dz) / 3.d0

  ! The distance is measured in the units of (average) cell size.
  ! Note: we assume same cell sizes along one dimension.
  lx = dx / avg_cell_size
  ly = dy / avg_cell_size
  lz = dz / avg_cell_size

  nbx = model%grid_full%nx
  nby = model%grid_full%ny
  nbz = model%grid_full%nz

  model_temp = 0.d0

  do k = 1, nbz
    do j = 1, nby
      do i = 1, nbx
        i_line = model%grid_full%ind(i, j, k)

        scale = 0.d0

        do k1 = max(k - lfgauss, 1), min(k + lfgauss, nbz)
          dz1_square = (dble(k1 - k) * lz)**2

          do j1 = max(j - lfgauss, 1), min(j + lfgauss, nby)
            dy1_square = (dble(j1 - j) * ly)**2

            do i1 = max(i - lfgauss, 1), min(i + lfgauss, nbx)
              dx1_square = (dble(i1 - i) * lx)**2

              gauss = gauss_const * exp(- (dx1_square + dy1_square + dz1_square) / (2.d0 * sigma**2))

              i_line_1 = model%grid_full%ind(i1, j1, k1)

              if (type == 1) then
              ! Arithmetic average.
                model_temp(i_line) = model_temp(i_line) + gauss * model%val_full(i_line_1)

              else if (type == 2) then
              ! Log-average.
                if (model%val_full(i_line_1) > 0.d0) then
                  model_temp(i_line) = model_temp(i_line) + gauss * log(model%val_full(i_line_1))
                else
                  print *, 'Bad log argument in apply_Gaussian_filter!'
                  stop
                endif
              endif

              scale = scale + gauss
            enddo
          enddo
        enddo

        model_temp(i_line) = model_temp(i_line) / scale

      enddo
    enddo
  enddo

  if (type == 1) then
    model%val_full = model_temp
  else if (type == 2) then
    model%val_full = exp(model_temp)
  endif

  deallocate(model_temp)

end subroutine apply_Gaussian_filter

end module filter
