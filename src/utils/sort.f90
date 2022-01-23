
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

!==============================================================================================================
! Functions for sorting an array of numbers.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!==============================================================================================================
module sort

  use global_typedefs, only: CUSTOM_REAL

  implicit none

  private

  public :: quicksort

contains

!================================================================
! Quicksort.
! VO: Took this from https://gist.github.com/t-nissie/479f0f16966925fa29ea
!================================================================
recursive subroutine quicksort(a, first, last)
  integer, intent(in) :: first, last
  real(kind=CUSTOM_REAL), intent(inout) :: a(:)

  real(kind=CUSTOM_REAL) :: x, t
  integer :: i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort

end module sort
