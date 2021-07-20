
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

!==============================================================================================================
! Contains various helper functions used in the code, which are not dependent on other modules,
! do not logically belong to a particular module, and can be easily moved here.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015.
!==============================================================================================================
module utils

  use global_typedefs, only: CUSTOM_REAL, NORM_L2

  implicit none

  private

  public :: compute_snrm
  public :: enforce_pb
  public :: der
  public :: avg2
  public :: avg4
  public :: getstep

contains

!==============================================================================================================
!! DK DK and DG DG: this is *not* BLAS1 SNRM2() because it has no square root at the end
!! DK DK and DG DG: and there is also a if statement to choose between L1 and L2 norm

!! DK DK and DG DG: it seems to us that the L1 norm does not work anyway because in the calling program
!! DK DK and DG DG: we then *always* take the square root of the result, which makes no sense for L1.

!! DK DK and DG DG: also, the mpi_reduce which is called in the main program after calls to this routine
!! DK DK and DG DG: always uses MPI_SUM and not MPI_MAX, thus again it would not work for L1.

!! DK DK and DG DG: thus, there is room for improvement here because if we only consider L2 then we can probably
!! DK DK and DG DG: switch to a BLAS1 call (?)
!==============================================================================================================
subroutine compute_snrm(sx,itypenorm,nr,ntheta,nz,i_begin,j_begin,k_begin,i_end,j_end,k_end,snrm)

  integer, intent(in) :: nr,ntheta,nz
  real(kind=CUSTOM_REAL), intent(in) :: sx(0:nr+1,0:ntheta+1,0:nz+1)
  !! select the norm to use depending on this parameter "itypenorm"
  integer, intent(in) :: itypenorm
  integer, intent(in) :: i_begin,j_begin,k_begin,i_end,j_end,k_end
  real(kind=CUSTOM_REAL), intent(out) :: snrm

  integer :: i,j,k,isamax,jsamax,ksamax

  if (itypenorm == NORM_L2) then
    !! L2 norm here, without square root calculation at the end
    snrm=0._CUSTOM_REAL
    do k=k_begin,k_end
      do j=j_begin,j_end
        do i=i_begin,i_end
          snrm=snrm+sx(i,j,k)**2
        enddo
      enddo
    enddo

  else
    ! DK DK and DG DG: L1 norm (max norm) here
    isamax=1
    jsamax=1
    ksamax=1
    do k=k_begin,k_end
      do j=j_begin,j_end
        do i=i_begin,i_end
          if(abs(sx(i,j,k))>abs(sx(isamax,jsamax,ksamax)))then
            isamax=i
            jsamax=j
            ksamax=k
          endif
        enddo
      enddo
    enddo
    snrm=abs(sx(isamax,jsamax,ksamax))

  endif

end subroutine compute_snrm

!==========================================================================
! Enforce periodic boundary conditions on x (in cylindrical coordinates).
!==========================================================================
pure subroutine enforce_pb(nr, ntheta, nz, x)
  integer, intent(in) :: nr, ntheta, nz
  real(kind=CUSTOM_REAL), intent(inout) :: x(0:nr+1, 0:ntheta+1, 0:nz+1)

  integer :: i, k

  do k = 0, nz + 1
    do i = 0, nr + 1
      x(i, 0, k) = x(i, ntheta, k)
      x(i, ntheta + 1, k) = x(i, 1, k)
    enddo
  enddo
end subroutine enforce_pb

!====================================================================
! Four-points stencil first derivative.
!====================================================================
function der(f0,f1,f2,f3,h)
  real(kind=CUSTOM_REAL), intent(in) :: f0, f1, f2, f3
  real(kind=CUSTOM_REAL), intent(in) :: h
  real(kind=CUSTOM_REAL) :: der

  ! 2nd order interpolation.
  der = (2._CUSTOM_REAL*f3 - 9._CUSTOM_REAL*f2 + 18._CUSTOM_REAL*f1 - &
         11._CUSTOM_REAL*f0) / (6._CUSTOM_REAL*h)

  ! 1st order interpolation.
  !der = (f1 - f0) / h

  return
end function der

!====================================================================
! Physically accurate average of 2 values.
!====================================================================
pure function avg2(a1, a2) result(res)
  real(kind=CUSTOM_REAL), intent(in) :: a1, a2
  real(kind=CUSTOM_REAL) :: res

  ! Option 1: linear.
!  res = 0.5_CUSTOM_REAL*(a1+a2)

  ! Option 2: quadratic - more accurate.
  if (a1 + a2 == 0._CUSTOM_REAL) then
    res = 0._CUSTOM_REAL
  else
    res = 2._CUSTOM_REAL * (a1 * a2) / (a1 + a2)
  endif
end function avg2

!========================================================================
! Physically accurate average of 4 values.
!========================================================================
pure function avg4(a1, a2, a3, a4) result(res)
  real(kind=CUSTOM_REAL), intent(in) :: a1, a2, a3, a4
  real(kind=CUSTOM_REAL) :: denom
  real(kind=CUSTOM_REAL) :: res

  ! Option 1: linear.
!  res = 0.25_CUSTOM_REAL * (a1 + a2 + a3 + a4)

  ! Option 2: quadratic - more accurate.
  !denom = a2*a3*a4 + a1*a3*a4 + a1*a2*a4 + a1*a2*a3
  denom = a3 * a4 * (a1 + a2) + a1 * a2 * (a3 + a4)
  if (denom == 0._CUSTOM_REAL) then
    res = 0._CUSTOM_REAL
  else
    res = 4._CUSTOM_REAL * (a1 * a2 * a3 * a4) / denom
  endif
end function avg4

!=================================================================================
! Trapezoidal integration rule (with varying step).
! Return integration step.
! di -- change in index i of data array a, di=1 for all except radius, where di=2.
!=================================================================================
function getstep(i,di,istart,iend,n,a)
  integer, intent(in) :: i,di,istart,iend,n
  real(kind=CUSTOM_REAL), intent(in) :: a(0:n)
  real(kind=CUSTOM_REAL) :: getstep

  ! Sanity check.
  if (i < 0 .or. i > n) then
    print *, 'Error in getstep(). Exiting.'
    stop
  endif

  if (i == istart) then
    getstep = 0.5_CUSTOM_REAL*abs(a(i+di) - a(i))
  else if (i == iend) then
    getstep = 0.5_CUSTOM_REAL*abs(a(i) - a(i-di))
  else
    getstep = 0.5_CUSTOM_REAL*abs(a(i+di) - a(i-di))
  endif

  return
end function getstep

end module utils
