
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

module wavelet_transform

  use global_typedefs, only: CUSTOM_REAL

  implicit none

  private

  public :: forward_wavelet
  public :: inverse_wavelet

  public :: Haar3D
  public :: iHaar3D

  public :: DaubD43D
  public :: iDaubD43D

contains

!=====================================================================================================
! Forward wavelet transform.
!=====================================================================================================
subroutine forward_wavelet(s, n1, n2, n3, wavelet_type)
  integer, intent(in) :: n1, n2, n3, wavelet_type
  real(kind=CUSTOM_REAL), intent(inout) :: s(n1, n2, n3)

  select case(wavelet_type)
    case(1)
      call Haar3D(s, n1, n2, n3)
    case(2)
      call DaubD43D(s, n1, n2, n3)
    case default
      print *, 'Unknown wavelet type!'
      stop
  end select

end subroutine forward_wavelet

!=====================================================================================================
! Inverse wavelet transform.
!=====================================================================================================
subroutine inverse_wavelet(s, n1, n2, n3, wavelet_type)
  integer, intent(in) :: n1, n2, n3, wavelet_type
  real(kind=CUSTOM_REAL), intent(inout) :: s(n1, n2, n3)

  select case(wavelet_type)
    case(1)
      call iHaar3D(s, n1, n2, n3)
    case(2)
      call iDaubD43D(s, n1, n2, n3)
    case default
      print *, 'Unknown wavelet type!'
      stop
  end select

end subroutine inverse_wavelet

!=====================================================================================================
! Haar wavelet transform (adapted from code by Sebastien Chevrot).
!=====================================================================================================
subroutine Haar3D(s, n1, n2, n3)
  integer, intent(in) :: n1, n2, n3
  real(kind=CUSTOM_REAL), intent(inout) :: s(n1, n2, n3)

  integer :: i, ic, L, il, ig, ngmin, ngmax
  integer :: istep, step_incr, step2, nscale, ng

  ! Loop over the 3 dimensions.
  do ic = 1, 3
    if (ic == 1) then
      nscale = int(log(real(n1, CUSTOM_REAL)) / log(2._CUSTOM_REAL))
      L = n1
    else if (ic == 2) then
      nscale = int(log(real(n2, CUSTOM_REAL)) / log(2._CUSTOM_REAL))
      L = n2
    else
      nscale = int(log(real(n3, CUSTOM_REAL)) / log(2._CUSTOM_REAL))
      L = n3
    endif

    ! Loop over the scales.
    do istep = 1, nscale
      step_incr = 2**istep
      ngmin = step_incr / 2 + 1
      ngmax = ngmin + int((L - ngmin) / step_incr) * step_incr
      ng = (ngmax - ngmin) / step_incr + 1
      step2 = step_incr

      ! Predict.
      ig = ngmin
      il = 1
      do i = 1, ng
        if (ic == 1) then
          s(ig, :, :) = s(ig, :, :) - s(il, :, :)
        else if (ic == 2) then
          s(:, ig, :) = s(:, ig, :) - s(:, il, :)
        else
          s(:, :, ig) = s(:, :, ig) - s(:, :, il)
        endif
        il = il + step2
        ig = ig + step2
      enddo

      ! Update.
      ig = ngmin
      il = 1
      do i = 1, ng
        if (ic == 1) then
          s(il, :, :) = s(il, :, :) + s(ig, :, :) / 2._CUSTOM_REAL
        else if (ic == 2) then
          s(:, il, :) = s(:, il, :) + s(:, ig, :) / 2._CUSTOM_REAL
        else
          s(:, :, il) = s(:, :, il) + s(:, :, ig) / 2._CUSTOM_REAL
        endif
        il = il + step2
        ig = ig + step2
      enddo

      ! Normalization.
      ig = ngmin
      il = 1
      do i = 1, ng
        if (ic == 1) then
          s(il, :, :) = s(il, :, :) * sqrt(2._CUSTOM_REAL)
          s(ig, :, :) = s(ig, :, :) / sqrt(2._CUSTOM_REAL)
        else if (ic == 2) then
          s(:, il, :) = s(:, il, :) * sqrt(2._CUSTOM_REAL)
          s(:, ig, :) = s(:, ig, :) / sqrt(2._CUSTOM_REAL)
        else
          s(:, :, il) = s(:, :, il) * sqrt(2._CUSTOM_REAL)
          s(:, :, ig) = s(:, :, ig) / sqrt(2._CUSTOM_REAL)
        endif
        il = il + step2
        ig = ig + step2
      enddo
    enddo
  enddo

end subroutine Haar3D

!=====================================================================================================
! Inverse Haar transform (adapted from code by Sebastien Chevrot).
!=====================================================================================================
subroutine iHaar3D(s, n1, n2, n3)
  integer, intent(in) :: n1, n2, n3
  real(kind=CUSTOM_REAL), intent(inout) :: s(n1, n2, n3)

  integer :: i, ic, L, il, ig, ngmin, ngmax
  integer :: istep, step_incr, step2, nscale, ng

  ! Loop over the 3 dimensions.
  do ic = 1, 3
    if (ic == 1) then
      nscale = int(log(real(n1, CUSTOM_REAL)) / log(2._CUSTOM_REAL))
      L = n1
    else if (ic == 2) then
      nscale = int(log(real(n2, CUSTOM_REAL)) / log(2._CUSTOM_REAL))
      L = n2
    else
      nscale = int(log(real(n3, CUSTOM_REAL)) / log(2._CUSTOM_REAL))
      L = n3
    endif

    ! Loop over the scales.
    do istep = nscale, 1, -1
      step_incr = 2**istep
      ngmin = step_incr / 2 + 1
      ngmax = ngmin + int((L - ngmin) / step_incr) * step_incr
      ng = (ngmax - ngmin) / step_incr + 1
      step2 = step_incr

      ! Normalization.
      ig = ngmin
      il = 1
      do i = 1, ng
        if (ic == 1) then
          s(il, :, :) = s(il, :, :) / sqrt(2._CUSTOM_REAL)
          s(ig, :, :) = s(ig, :, :) * sqrt(2._CUSTOM_REAL)
        else if (ic == 2) then
          s(:, il, :) = s(:, il, :) / sqrt(2._CUSTOM_REAL)
          s(:, ig, :) = s(:, ig, :) * sqrt(2._CUSTOM_REAL)
        else
          s(:, :, il) = s(:, :, il) / sqrt(2._CUSTOM_REAL)
          s(:, :, ig) = s(:, :, ig) * sqrt(2._CUSTOM_REAL)
        endif
        il = il + step2
        ig = ig + step2
      enddo

      ! Update.
      ig = ngmin
      il = 1
      do i = 1, ng
        if (ic == 1) then
          s(il, :, :) = s(il, :, :) - s(ig, :, :) / 2._CUSTOM_REAL
        else if (ic == 2) then
          s(:, il, :) = s(:, il, :) - s(:, ig, :) / 2._CUSTOM_REAL
        else
          s(:, :, il) = s(:, :, il) - s(:, :, ig) / 2._CUSTOM_REAL
        endif
        il = il + step2
        ig = ig + step2
      enddo

      ! Predict.
      ig = ngmin
      il = 1
      do i = 1, ng
        if (ic == 1) then
          s(ig, :, :) = s(ig, :, :) + s(il, :, :)
        else if (ic == 2) then
          s(:, ig, :) = s(:, ig, :) + s(:, il, :)
        else
          s(:, :, ig) = s(:, :, ig) + s(:, :, il)
        endif
        il = il + step2
        ig = ig + step2
      enddo
    enddo
  enddo

end subroutine iHaar3D

!=====================================================================================================
! Daubechies D4 transform.
! Boundary conditions (the edge problem) are adapted from the algorithm by Ian Kaplan (2001) found at:
! http://www.bearcave.com/misl/misl_tech/wavelets/daubechies/index.html
!=====================================================================================================
subroutine DaubD43D(s, n1, n2, n3)
  integer, intent(in) :: n1, n2, n3
  real(kind=CUSTOM_REAL), intent(inout) :: s(n1, n2, n3)

  integer :: i, ic, L, il, ig, ngmin, ngmax, ilmax
  integer :: istep, step_incr, step2, nscale, ng
  real(kind=CUSTOM_REAL) :: c0, c1, c2, c3, c4

  c0 = sqrt(3._CUSTOM_REAL)
  c1 = sqrt(3._CUSTOM_REAL) / 4._CUSTOM_REAL
  c2 = (sqrt(3._CUSTOM_REAL) - 2._CUSTOM_REAL) / 4._CUSTOM_REAL
  c3 = (sqrt(3._CUSTOM_REAL) - 1._CUSTOM_REAL) / sqrt(2._CUSTOM_REAL)
  c4 = (sqrt(3._CUSTOM_REAL) + 1._CUSTOM_REAL) / sqrt(2._CUSTOM_REAL)

  ! Loop over the 3 dimensions.
  do ic = 1, 3
     if (ic == 1) then
        nscale = int(log(real(n1, CUSTOM_REAL)) / log(2._CUSTOM_REAL))
        L = n1
     else if (ic == 2) then
        nscale = int(log(real(n2, CUSTOM_REAL)) / log(2._CUSTOM_REAL))
        L = n2
     else
        nscale = int(log(real(n3, CUSTOM_REAL)) / log(2._CUSTOM_REAL))
        L = n3
     endif

     ! Loop over the scales.
     do istep = 1, nscale
        step_incr = 2**istep
        ngmin = step_incr / 2 + 1
        ngmax = ngmin + int((L - ngmin) / step_incr) * step_incr
        ng = (ngmax - ngmin) / step_incr + 1
        step2 = step_incr
        ! The last index of the odd indexes (corresponding to the first half of an array). Corresponds to S[half-1].
        ilmax = 1 + (ng - 1) * step2

        ! Update 1.
        ig = ngmin
        il = 1
        do i = 1, ng
           if (ic == 1) then
              s(il, :, :) = s(il, :, :) + s(ig, :, :) * c0
           else if (ic == 2) then
              s(:, il, :) = s(:, il, :) + s(:, ig, :) * c0
           else
              s(:, :, il) = s(:, :, il) + s(:, :, ig) * c0
           endif
           il = il + step2
           ig = ig + step2
        enddo

        ! Predict.
        ! Applying boundary conditions.
        il = 1
        ig = ngmin
        if (ic == 1) then
           s(ig, :, :) = s(ig, :, :) - s(il, :, :) * c1 - s(ilmax, :, :) * c2
        else if (ic == 2) then
           s(:, ig, :) = s(:, ig, :) - s(:, il, :) * c1 - s(:, ilmax, :) * c2
        else
           s(:, :, ig) = s(:, :, ig) - s(:, :, il) * c1 - s(:, :, ilmax) * c2
        endif

        ig = ngmin + step2
        il = 1 + step2
        do i = 1, ng - 1
           if (ic == 1) then
              s(ig, :, :) = s(ig, :, :) - s(il, :, :) * c1 - s(il - step2, :, :) * c2
           else if (ic == 2) then
              s(:, ig, :) = s(:, ig, :) - s(:, il, :) * c1 - s(:, il - step2, :) * c2
           else
              s(:, :, ig) = s(:, :, ig) - s(:, :, il) * c1 - s(:, :, il - step2) * c2
           endif
           il = il + step2
           ig = ig + step2
        enddo

        ! Update 2.
        ig = ngmin
        il = 1
        do i = 1, ng - 1
           if (ic == 1) then
              s(il, :, :) = s(il, :, :) - s(ig + step2, :, :)
           else if (ic == 2) then
              s(:, il, :) = s(:, il, :) - s(:, ig + step2, :)
           else
              s(:, :, il) = s(:, :, il) - s(:, :, ig + step2)
           endif
           il = il + step2
           ig = ig + step2
        enddo

        ! Applying boundary conditions (S[half-1] = S[half-1] - S[half]).
        ig = ngmin
        il = ilmax
        if (ic == 1) then
          s(il, :, :) = s(il, :, :) - s(ig, :, :)
        else if (ic == 2) then
          s(:, il, :) = s(:, il, :) - s(:, ig, :)
        else
          s(:, :, il) = s(:, :, il) - s(:, :, ig)
        endif

        ! Normalization.
        ig = ngmin
        il = 1
        do i = 1, ng
           if (ic == 1) then
              s(il, :, :) = s(il, :, :) * c3
              s(ig, :, :) = s(ig, :, :) * c4
           else if (ic == 2) then
              s(:, il, :) = s(:, il, :) * c3
              s(:, ig, :) = s(:, ig, :) * c4
           else
              s(:, :, il) = s(:, :, il) * c3
              s(:, :, ig) = s(:, :, ig) * c4
           endif
           il = il + step2
           ig = ig + step2
        enddo
     enddo
  enddo

end subroutine DaubD43D

!=====================================================================================================
! Inverse Daubechies D4 transform.
! Boundary conditions (the edge problem) are adapted from the algorithm by Ian Kaplan (2001) found at
! http://www.bearcave.com/misl/misl_tech/wavelets/daubechies/index.html
!=====================================================================================================
subroutine iDaubD43D(s, n1, n2, n3)
  integer, intent(in) :: n1, n2, n3
  real(kind=CUSTOM_REAL), intent(inout) :: s(n1, n2, n3)

  integer :: i, ic, L, il, ig, ngmin, ngmax, ilmax
  integer :: istep, step_incr, step2, nscale, ng
  real(kind=CUSTOM_REAL) :: c0, c1, c2, c3, c4

  c0 = sqrt(3._CUSTOM_REAL)
  c1 = sqrt(3._CUSTOM_REAL) / 4._CUSTOM_REAL
  c2 = (sqrt(3._CUSTOM_REAL) - 2._CUSTOM_REAL) / 4._CUSTOM_REAL
  c3 = (sqrt(3._CUSTOM_REAL) - 1._CUSTOM_REAL) / sqrt(2._CUSTOM_REAL)
  c4 = (sqrt(3._CUSTOM_REAL) + 1._CUSTOM_REAL) / sqrt(2._CUSTOM_REAL)

  ! Loop over the 3 dimensions.
  do ic = 1, 3
     if (ic == 1) then
        nscale = int(log(real(n1, CUSTOM_REAL)) / log(2._CUSTOM_REAL))
        L = n1
     else if (ic == 2) then
        nscale = int(log(real(n2, CUSTOM_REAL)) / log(2._CUSTOM_REAL))
        L = n2
     else
        nscale = int(log(real(n3, CUSTOM_REAL)) / log(2._CUSTOM_REAL))
        L = n3
     endif

     ! Loop over the scales.
     do istep = nscale, 1, -1
        step_incr = 2**istep
        ngmin = step_incr / 2 + 1
        ngmax = ngmin + int((L - ngmin) / step_incr) * step_incr
        ng = (ngmax - ngmin) / step_incr + 1
        step2 = step_incr
        ! The last index of the odd indexes (corresponding to the first half of an array). Corresponds to S[half-1].
        ilmax = 1 + (ng - 1) * step2

        ! Normalization.
        ig = ngmin
        il = 1
        do i = 1, ng
           if (ic == 1) then
              s(il, :, :) = s(il, :, :) * c4
              s(ig, :, :) = s(ig, :, :) * c3
           else if (ic == 2) then
              s(:, il, :) = s(:, il, :) * c4
              s(:, ig, :) = s(:, ig, :) * c3
           else
              s(:, :, il) = s(:, :, il) * c4
              s(:, :, ig) = s(:, :, ig) * c3
           endif
           il = il + step2
           ig = ig + step2
        enddo

        ! Update 2.
        ig = ngmin + (ng - 2) * step2
        il = 1 + (ng - 2) * step2
        do i = 1, ng - 1
           if (ic == 1) then
              s(il, :, :) = s(il, :, :) + s(ig + step2, :, :)
           else if (ic == 2) then
              s(:, il, :) = s(:, il, :) + s(:, ig + step2, :)
           else
              s(:, :, il) = s(:, :, il) + s(:, :, ig + step2)
           endif
           il = il - step2
           ig = ig - step2
        enddo

        ! Applying boundary conditions (S[half-1] = S[half-1] - S[half]).
        ig = ngmin
        il = ilmax
        if (ic == 1) then
          s(il, :, :) = s(il, :, :) + s(ig, :, :)
        else if (ic == 2) then
          s(:, il, :) = s(:, il, :) + s(:, ig, :)
        else
          s(:, :, il) = s(:, :, il) + s(:, :, ig)
        endif

        ! Predict.
        ig = ngmin + (ng) * step2
        il = 1 + (ng) * step2
        do i = 1, ng - 1
           il = il - step2
           ig = ig - step2
           if (ic == 1) then
              s(ig, :, :) = s(ig, :, :) + s(il, :, :) * c1 + s(il - step2, :, :) * c2
           else if (ic == 2) then
              s(:, ig, :) = s(:, ig, :) + s(:, il, :) * c1 + s(:, il - step2, :) * c2
           else
              s(:, :, ig) = s(:, :, ig) + s(:, :, il) * c1 + s(:, :, il - step2) * c2
           endif
        enddo

        ! Applying boundary conditions.
        il = 1
        ig = ngmin
        if (ic == 1) then
           s(ig, :, :) = s(ig, :, :) + s(il, :, :) * c1 + s(ilmax, :, :) * c2
        else if (ic == 2) then
           s(:, ig, :) = s(:, ig, :) + s(:, il, :) * c1 + s(:, ilmax, :) * c2
        else
           s(:, :, ig) = s(:, :, ig) + s(:, :, il) * c1 + s(:, :, ilmax) * c2
        endif

        ! Update 1.
        ig = ngmin
        il = 1
        do i = 1, ng
           if (ic == 1) then
              s(il, :, :) = s(il, :, :) - s(ig, :, :) * c0
           else if (ic == 2) then
              s(:, il, :) = s(:, il, :) - s(:, ig, :) * c0
           else
              s(:, :, il) = s(:, :, il) - s(:, :, ig) * c0
           endif
           il = il + step2
           ig = ig + step2
        enddo
     enddo
  enddo

end subroutine iDaubD43D

end module wavelet_transform
