
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

module source

  use global_typedefs
  use parameters_ect

  implicit none

  private

  public :: source_RHS

contains

! Prepare the source vector b for the linear system a*x=b to be solved by the linear solver.
subroutine source_RHS(b,flag,dims,val)
  type(t_dimensions), intent(in) :: dims
  integer, intent(in) :: flag(0:,0:,0:)
  real(kind=CUSTOM_REAL), intent(in) :: val(0:,0:,0:)

  real(kind=CUSTOM_REAL), intent(out) :: b(0:,0:,0:)

  integer :: i,j,k

  ! Initialize b to zero.
  b = 0._CUSTOM_REAL

  ! Compute the source vector b of the system a*x=b.
  do k=1,dims%nzlocal
    do j=1,dims%ntheta
      do i=1,dims%nr
        !! DK DK again, this memory copy is useless (if we fuse b and val everywhere)
        if (flag(i,j,k) == 1) b(i,j,k) = val(i,j,k)
      enddo
    enddo
  enddo

end subroutine source_RHS

end module source

