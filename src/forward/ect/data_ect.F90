
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

module data_ect

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use parameters_ect

  implicit none

  private

  public :: check_data
  public :: ect_normalization
  public :: write_data_for_inversion
  public :: write_data_to_file

contains

!================================================================================
! Check that data_high > data_low.
!================================================================================
subroutine check_data(par, ndata, data_low, data_high, myrank)
  type(t_parameters_ect), intent(in) :: par
  integer, intent(in) :: ndata
  real(kind=CUSTOM_REAL), intent(in) :: data_low(:)
  real(kind=CUSTOM_REAL), intent(in) :: data_high(:)
  integer, intent(in) :: myrank

  integer :: i, failed

  failed = 0

  ! Check data.
  do i = 1, ndata
    if (abs(data_high(i)) <= abs(data_low(i))) then
      print *, "i, data_low, data_high =", i, data_low(i), data_high(i), &
              " Test failed at pair", par%get_electrode_index(i, 1), par%get_electrode_index(i, 2)
      failed = failed + 1
    endif
  enddo

  if (failed /= 0) then
    print *, "Data check (data_high <= data_low) failed n=", failed, " times!"
    call exit_MPI("Data check failed: data_high <= data_low!", myrank, 0)
  else
    print *, "Data check (data_high <= data_low) passed!"
  endif
end subroutine

!============================================================================================
! Apply normalization typical in ECT: data_new = (data - data_low) / (data_high - data_low),
! where data_low and data_high correspond to data obtained with an empty and full sensors.
! For the sensitivity kernel,  (dCij/dEpsilon)_new = dCij/dEpsilon/(data_high - data_low),
!============================================================================================
subroutine ect_normalization(ndata, sensit, residuals, data_high, data_low)

  integer, intent(in) :: ndata
  real(kind=CUSTOM_REAL), intent(in) :: data_high(:)
  real(kind=CUSTOM_REAL), intent(in) :: data_low(:)

  ! Sensitivity matrix.
  real(kind=CUSTOM_REAL), intent(inout) :: sensit(:, :)
  real(kind=CUSTOM_REAL), intent(inout) :: residuals(:)

  integer :: i

  ! Normalize data and sensitivity kernel.
  do i = 1, ndata
    residuals(i) = residuals(i) / (data_high(i) - data_low(i))

    sensit(:, i) = sensit(:, i) / (data_high(i) - data_low(i))
  enddo

end subroutine ect_normalization

!==========================================================================================
! Write the capacitance data to file.
!==========================================================================================
subroutine write_data_to_file(par, capacity, file_name, myrank)
  ! Parameters.
  type(t_parameters_ect), intent(in) :: par
  ! Capacitance data.
  real(kind=CUSTOM_REAL), intent(in) :: capacity(:, :)
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: myrank

  integer :: ielectrode, jelectrode, l, ierr
  character(len=160) :: file_name_full

  if (myrank == 0) then
    file_name_full = trim(path_output)//file_name
    open(16,file=trim(file_name_full), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) call exit_MPI("Error in writing "//trim(file_name_full), myrank, ierr)

    l = 0
    do jelectrode = 1, par%nel
      do ielectrode = 1, par%nel
        l = l + 1
        write(16,*) l, ielectrode, jelectrode, capacity(ielectrode, jelectrode)
      enddo
    enddo
    close(16)
  endif

end subroutine write_data_to_file

!==========================================================================================
! Write 1D data array with non-repeated data, excluding self-capacitances (Cij, i>j),
! to be used for inversion.
!==========================================================================================
subroutine write_data_for_inversion(par, capacity, cdata)
  ! Parameters.
  type(t_parameters_ect), intent(in) :: par
  ! Capacitance data.
  real(kind=CUSTOM_REAL), intent(in) :: capacity(:, :)
  ! Capacitance data 1D.
  real(kind=CUSTOM_REAL), intent(out) :: cdata(:)

  integer :: ielectrode, jelectrode, l

  l = 0
  do jelectrode = 1, par%nel
    do ielectrode = 1, par%nel

      if (par%accept_electrode_pair(ielectrode, jelectrode)) then
        l = l + 1
        cdata(l) = capacity(ielectrode, jelectrode)
      endif

    enddo
  enddo

end subroutine write_data_for_inversion

end module data_ect
