
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
  use mpi_tools, only: exit_MPI
  use parallel_tools

  implicit none

  private

  public :: Haar3D
  public :: Haar3D_serial
  public :: iHaar3D
  public :: DaubD43D
  public :: iDaubD43D

  private :: Haar3D_all
  private :: calculate_parallel_array_bounds

contains

!=====================================================================================================
! Subroutine for computing bounds for array slices for parallel wavelets.
!=====================================================================================================
subroutine calculate_parallel_array_bounds(n, i1, i2, n_loc, myrank, nbproc)
  integer, intent(in) :: n, myrank, nbproc
  integer, intent(out) :: i1, i2, n_loc

  type(t_parallel_tools) :: pt

  ! Sanity check.
  if (nbproc > n) then
    ! TODO: Need to adjust the subroutine for this case.
    call exit_MPI("Number of CPUs is greater than array size in calculate_parallel_array_bounds!", myrank, 0)
  endif

  ! Define MPI partitioning.
  n_loc = pt%calculate_nelements_at_cpu(n, myrank, nbproc)

  if (myrank < nbproc - 1) then
    i1 = myrank * n_loc + 1
    i2 = myrank * n_loc + n_loc
  else
  ! Last rank.
    i1 = n - n_loc + 1
    i2 = n
  endif
end subroutine calculate_parallel_array_bounds

!=====================================================================================================
! Routine needed for Haar wavelet transform parallelisation.
! By calling it we convert 3D array slices to 1D arrays, needed for parallelisation.
! Author: Vitaliy Ogarko.
!=====================================================================================================
subroutine Haar3D_all(n, a, b, myrank, nbproc)
  integer, intent(in) :: n, myrank, nbproc
  real(kind=CUSTOM_REAL), intent(inout) :: a(n)
  real(kind=CUSTOM_REAL), intent(inout) :: b(n)

  type(t_parallel_tools) :: pt
  integer :: i1, i2, n_loc

  ! Calculate array bounds to process at this CPU.
  call calculate_parallel_array_bounds(n, i1, i2, n_loc, myrank, nbproc)

  ! Predict.
  a(i1:i2) = a(i1:i2) - b(i1:i2)

  ! Update.
  b(i1:i2) = b(i1:i2) + a(i1:i2) / 2._CUSTOM_REAL

  ! Normalize.
  a(i1:i2) = a(i1:i2) / sqrt(2._CUSTOM_REAL)
  b(i1:i2) = b(i1:i2) * sqrt(2._CUSTOM_REAL)

  ! Update the full array on all CPUs.
  call pt%get_full_array_in_place2(n_loc, a, myrank, nbproc)
  call pt%get_full_array_in_place2(n_loc, b, myrank, nbproc)

end subroutine Haar3D_all

!=====================================================================================================
! Sanity check for the min/max bounds of the 3D array slices.
!=====================================================================================================
pure function check_array_bounds(n, ig, igmax, il, ilmax) result (passed)
  integer, intent(in) :: n, ig, igmax, il, ilmax
  logical :: passed

  passed = .true.
  if (ig < 1 .or. ig > n) passed = .false.
  if (il < 1 .or. il > n) passed = .false.
  if (igmax < 1 .or. igmax > n) passed = .false.
  if (ilmax < 1 .or. ilmax > n) passed = .false.
end function check_array_bounds

!=====================================================================================================
! Parallel Haar wavelet transform.
! Author: Vitaliy Ogarko
!=====================================================================================================
subroutine Haar3D(s, n1, n2, n3, myrank, nbproc)
  integer, intent(in) :: n1, n2, n3
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(inout) :: s(n1, n2, n3)

  integer :: ic, L, il, ig, ngmin, ngmax
  integer :: istep, step_incr, step2, nscale, ng
  integer :: igmax, ilmax

  ! Auxiliary arrays.
  real(kind=CUSTOM_REAL), allocatable :: a(:, :, :)
  real(kind=CUSTOM_REAL), allocatable :: b(:, :, :)

  ! Loop over the 3 dimensions.
  do ic = 1, 3

    ! Loop over the scales.
    if (ic == 1) then
      nscale = int(log(real(n1,CUSTOM_REAL))/log(2._CUSTOM_REAL))
      L = n1
    else if (ic == 2) then
      nscale = int(log(real(n2,CUSTOM_REAL))/log(2._CUSTOM_REAL))
      L = n2
    else
      nscale = int(log(real(n3,CUSTOM_REAL))/log(2._CUSTOM_REAL))
      L = n3
    endif

    do istep = 1, nscale
      step_incr = 2**istep
      ngmin = step_incr/2+1
      ngmax = ngmin+int((L-ngmin)/step_incr)*step_incr
      ng = (ngmax-ngmin)/step_incr+1
      step2 = step_incr

      ig = ngmin
      il = 1
      igmax = ig + step2 * (ng - 1)
      ilmax = il + step2 * (ng - 1)

      ! Sanity check.
      if (.not. check_array_bounds(L, ig, igmax, il, ilmax)) then
        print *, ic, istep, L, ig, igmax, il, ilmax, myrank
        call exit_MPI("Bad array slices bounds in Haar3D!", myrank, 0)
      endif

      if (ic == 1) then
        allocate(a(ng, n2, n3))
        allocate(b(ng, n2, n3))

        a = s(ig:igmax:step2, 1:n2, 1:n3)
        b = s(il:ilmax:step2, 1:n2, 1:n3)

        call Haar3D_all(ng * n2 * n3, a, b, myrank, nbproc)

        s(ig:igmax:step2, 1:n2, 1:n3) = a
        s(il:ilmax:step2, 1:n2, 1:n3) = b

      else if (ic == 2) then
        allocate(a(n1, ng, n3))
        allocate(b(n1, ng, n3))

        a = s(1:n1, ig:igmax:step2, 1:n3)
        b = s(1:n1, il:ilmax:step2, 1:n3)

        call Haar3D_all(n1 * ng * n3, a, b, myrank, nbproc)

        s(1:n1, ig:igmax:step2, 1:n3) = a
        s(1:n1, il:ilmax:step2, 1:n3) = b

      else
        allocate(a(n1, n2, ng))
        allocate(b(n1, n2, ng))

        a = s(1:n1, 1:n2, ig:igmax:step2)
        b = s(1:n1, 1:n2, il:ilmax:step2)

        call Haar3D_all(n1 * n2 * ng, a, b, myrank, nbproc)

        s(1:n1, 1:n2, ig:igmax:step2) = a
        s(1:n1, 1:n2, il:ilmax:step2) = b
      endif

      deallocate(a)
      deallocate(b)

    enddo
  enddo

end subroutine Haar3D

!=====================================================================================================
! Haar wavelet transform adapted from code by Sebastien Chevrot.
!=====================================================================================================
subroutine Haar3D_serial(s,n1,n2,n3)

  integer, intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL), intent(inout) :: s(n1,n2,n3)

  integer :: i,i1,i2,i3,ic,L,il,ig,ngmin,ngmax
  integer :: istep,step_incr,step2,nscale,ng

! Loop over the 3 dimensions
  do ic = 1,3

! Loop over the scales
    if (ic==1) then
      nscale = int(log(real(n1,CUSTOM_REAL))/log(2._CUSTOM_REAL))
      L = n1
    else if (ic==2) then
      nscale = int(log(real(n2,CUSTOM_REAL))/log(2._CUSTOM_REAL))
      L = n2
    else
      nscale = int(log(real(n3,CUSTOM_REAL))/log(2._CUSTOM_REAL))
      L = n3
    endif
    do istep = 1,nscale
      step_incr = 2**istep
      ngmin = step_incr/2+1
      ngmax = ngmin+int((L-ngmin)/step_incr)*step_incr
      ng = (ngmax-ngmin)/step_incr+1
      step2 = step_incr
!-------------- Predict
      ig = ngmin
      il = 1
      do i = 1,ng
        if (ic==1) then
          forall(i2 = 1:n2, i3 = 1:n3) s(ig,i2,i3) = s(ig,i2,i3)-s(il,i2,i3)
        else if (ic==2) then
          forall(i1 = 1:n1, i3 = 1:n3) s(i1,ig,i3) = s(i1,ig,i3)-s(i1,il,i3)
        else
          forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,ig) = s(i1,i2,ig)-s(i1,i2,il)
        endif
        il = il+step2
        ig = ig+step2
      enddo
!------------- Update
      ig = ngmin
      il = 1
      do i = 1,ng
        if (ic==1) then
          forall(i2 = 1:n2, i3 = 1:n3) s(il,i2,i3) = s(il,i2,i3)+s(ig,i2,i3)/2._CUSTOM_REAL
        else if (ic==2) then
          forall(i1 = 1:n1, i3 = 1:n3) s(i1,il,i3) = s(i1,il,i3)+s(i1,ig,i3)/2._CUSTOM_REAL
        else
          forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,il) = s(i1,i2,il)+s(i1,i2,ig)/2._CUSTOM_REAL
        endif
        il = il+step2
        ig = ig+step2
      enddo
! Normalization
      ig = ngmin
      il = 1
      do i = 1,ng
        if (ic==1) then
          forall(i2 = 1:n2, i3 = 1:n3) s(il,i2,i3) = s(il,i2,i3)*sqrt(2._CUSTOM_REAL)
          forall(i2 = 1:n2, i3 = 1:n3) s(ig,i2,i3) = s(ig,i2,i3)/sqrt(2._CUSTOM_REAL)
        else if (ic==2) then
          forall(i1 = 1:n1, i3 = 1:n3) s(i1,il,i3) = s(i1,il,i3)*sqrt(2._CUSTOM_REAL)
          forall(i1 = 1:n1, i3 = 1:n3) s(i1,ig,i3) = s(i1,ig,i3)/sqrt(2._CUSTOM_REAL)
        else
          forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,il) = s(i1,i2,il)*sqrt(2._CUSTOM_REAL)
          forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,ig) = s(i1,i2,ig)/sqrt(2._CUSTOM_REAL)
        endif
        il = il+step2
        ig = ig+step2
      enddo
    enddo
  enddo

end subroutine Haar3D_serial

!=====================================================================================================
! Inverse Haar transform adapted from code by Sebastien Chevrot.
!=====================================================================================================
subroutine iHaar3D(s,n1,n2,n3)

  integer, intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL), intent(inout) :: s(n1,n2,n3)

  integer :: i,i1,i2,i3,ic,L,il,ig,ngmin,ngmax
  integer :: istep,step_incr,step2,nscale,ng

! Loop over the 3 dimensions
  do ic = 1,3

! Loop over the scales
    if (ic==1) then
      nscale = int(log(real(n1,CUSTOM_REAL))/log(2._CUSTOM_REAL))
      L = n1
    else if (ic==2) then
      nscale = int(log(real(n2,CUSTOM_REAL))/log(2._CUSTOM_REAL))
      L = n2
    else
      nscale = int(log(real(n3,CUSTOM_REAL))/log(2._CUSTOM_REAL))
      L = n3
    endif
    do istep = nscale,1,-1
      step_incr = 2**istep
      ngmin = step_incr/2+1
      ngmax = ngmin+int((L-ngmin)/step_incr)*step_incr
      ng = (ngmax-ngmin)/step_incr+1
      step2 = step_incr
!-------------- Normalization
      ig = ngmin
      il = 1
      do i = 1,ng
        if (ic==1) then
          forall(i2 = 1:n2, i3 = 1:n3) s(il,i2,i3) = s(il,i2,i3)/sqrt(2._CUSTOM_REAL)
          forall(i2 = 1:n2, i3 = 1:n3) s(ig,i2,i3) = s(ig,i2,i3)*sqrt(2._CUSTOM_REAL)
        else if (ic==2) then
          forall(i1 = 1:n1, i3 = 1:n3) s(i1,il,i3) = s(i1,il,i3)/sqrt(2._CUSTOM_REAL)
          forall(i1 = 1:n1, i3 = 1:n3) s(i1,ig,i3) = s(i1,ig,i3)*sqrt(2._CUSTOM_REAL)
        else
          forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,il) = s(i1,i2,il)/sqrt(2._CUSTOM_REAL)
          forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,ig) = s(i1,i2,ig)*sqrt(2._CUSTOM_REAL)
        endif
        il = il+step2
        ig = ig+step2
      enddo
!------------- Update
      ig = ngmin
      il = 1
      do i = 1,ng
        if (ic==1) then
          forall(i2 = 1:n2, i3 = 1:n3) s(il,i2,i3) = s(il,i2,i3)-s(ig,i2,i3)/2._CUSTOM_REAL
        else if (ic==2) then
          forall(i1 = 1:n1, i3 = 1:n3) s(i1,il,i3) = s(i1,il,i3)-s(i1,ig,i3)/2._CUSTOM_REAL
        else
          forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,il) = s(i1,i2,il)-s(i1,i2,ig)/2._CUSTOM_REAL
        endif
        il = il+step2
        ig = ig+step2
      enddo
!-------------- Predict
      ig = ngmin
      il = 1
      do i = 1,ng
        if (ic==1) then
          forall(i2 = 1:n2, i3 = 1:n3) s(ig,i2,i3) = s(ig,i2,i3)+s(il,i2,i3)
        else if (ic==2) then
          forall(i1 = 1:n1, i3 = 1:n3) s(i1,ig,i3) = s(i1,ig,i3)+s(i1,il,i3)
        else
          forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,ig) = s(i1,i2,ig)+s(i1,i2,il)
        endif
        il = il+step2
        ig = ig+step2
      enddo
    enddo
  enddo

end subroutine iHaar3D

!=====================================================================================================
! Daubechies D4 transform from algorithm found at
! http://www.bearcave.com/misl/misl_tech/wavelets/daubechies/index.html
!=====================================================================================================
subroutine DaubD43D(s,n1,n2,n3)

  integer, intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL), intent(inout) :: s(n1,n2,n3)

  integer :: i,i1,i2,i3,ic,L,il,ig,ngmin,ngmax
  integer :: istep,step_incr,step2,nscale,ng

! Loop over the 3 dimensions
  do ic = 1,3

! Loop over the scales
     if (ic==1) then
        nscale = int(log(real(n1,CUSTOM_REAL))/log(2._CUSTOM_REAL))
        L = n1
     else if (ic==2) then
        nscale = int(log(real(n2,CUSTOM_REAL))/log(2._CUSTOM_REAL))
        L = n2
     else
        nscale = int(log(real(n3,CUSTOM_REAL))/log(2._CUSTOM_REAL))
        L = n3
     endif
     do istep = 1,nscale
        step_incr = 2**istep
        ngmin = step_incr/2+1
        ngmax = ngmin+int((L-ngmin)/step_incr)*step_incr
        ng = (ngmax-ngmin)/step_incr+1
        step2 = step_incr

        !--------------Update 1
        ig = ngmin
        il = 1
        do i = 1,ng
           if (ic==1) then
              forall(i2 = 1:n2, i3 = 1:n3) s(il,i2,i3) = s(il,i2,i3)+s(ig,i2,i3)*sqrt(3._CUSTOM_REAL)
           else if (ic==2) then
              forall(i1 = 1:n1, i3 = 1:n3) s(i1,il,i3) = s(i1,il,i3)+s(i1,ig,i3)*sqrt(3._CUSTOM_REAL)
           else
              forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,il) = s(i1,i2,il)+s(i1,i2,ig)*sqrt(3._CUSTOM_REAL)
           endif
           il = il+step2
           ig = ig+step2
        enddo

        !-------------- Predict
        il=1
        ig=ngmin
        if (ic==1) then
           forall(i2 = 1:n2, i3 = 1:n3) s(ig,i2,i3) = s(ig,i2,i3) - s(il,i2,i3)*sqrt(3._CUSTOM_REAL)/4._CUSTOM_REAL
        else if (ic==2) then
           forall(i1 = 1:n1, i3 = 1:n3) s(i1,ig,i3) = s(i1,ig,i3) - s(i1,il,i3) *sqrt(3._CUSTOM_REAL)/4._CUSTOM_REAL
        else
           forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,ig) = s(i1,i2,ig) - s(i1,i2,il) *sqrt(3._CUSTOM_REAL)/4._CUSTOM_REAL
        endif

        ig = ngmin + step2
        il = 1 + step2
        do i = 1,ng-1
           if (ic==1) then
              forall(i2 = 1:n2, i3 = 1:n3) s(ig,i2,i3) = s(ig,i2,i3) - s(il,i2,i3)*sqrt(3._CUSTOM_REAL)/4._CUSTOM_REAL &
              - s(il-step2,i2,i3)*(sqrt(3._CUSTOM_REAL)-2._CUSTOM_REAL)/4._CUSTOM_REAL
           else if (ic==2) then
              forall(i1 = 1:n1, i3 = 1:n3) s(i1,ig,i3) = s(i1,ig,i3)-s(i1,il,i3) *sqrt(3._CUSTOM_REAL)/4._CUSTOM_REAL &
              - s(i1,il-step2,i3)*(sqrt(3._CUSTOM_REAL)-2._CUSTOM_REAL)/4._CUSTOM_REAL
           else
              forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,ig) = s(i1,i2,ig)-s(i1,i2,il) *sqrt(3._CUSTOM_REAL)/4._CUSTOM_REAL &
              - s(i1,i2,il-step2)*(sqrt(3._CUSTOM_REAL)-2._CUSTOM_REAL)/4._CUSTOM_REAL
           endif
           il = il+step2
           ig = ig+step2
        enddo


        !------------- Update 2

        ig = ngmin
        il = 1
        do i = 1,ng-1
           if (ic==1) then
              forall(i2 = 1:n2, i3 = 1:n3) s(il,i2,i3) = s(il,i2,i3)-s(ig+step2,i2,i3)
           else if (ic==2) then
              forall(i1 = 1:n1, i3 = 1:n3) s(i1,il,i3) = s(i1,il,i3)-s(i1,ig+step2,i3)
           else
              forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,il) = s(i1,i2,il)-s(i1,i2,ig+step2)
           endif
           il = il+step2
           ig = ig+step2
        enddo

        !--------------- Normalization
        ig = ngmin
        il = 1
        do i = 1,ng
           if (ic==1) then
              forall(i2 = 1:n2, i3 = 1:n3) s(il,i2,i3) = s(il,i2,i3) * (sqrt(3._CUSTOM_REAL) - 1._CUSTOM_REAL) /sqrt(2._CUSTOM_REAL)
              forall(i2 = 1:n2, i3 = 1:n3) s(ig,i2,i3) = s(ig,i2,i3) * (sqrt(3._CUSTOM_REAL) + 1._CUSTOM_REAL) /sqrt(2._CUSTOM_REAL)
           else if (ic==2) then
              forall(i1 = 1:n1, i3 = 1:n3) s(i1,il,i3) = s(i1,il,i3) * (sqrt(3._CUSTOM_REAL) - 1._CUSTOM_REAL) /sqrt(2._CUSTOM_REAL)
              forall(i1 = 1:n1, i3 = 1:n3) s(i1,ig,i3) = s(i1,ig,i3) * (sqrt(3._CUSTOM_REAL) + 1._CUSTOM_REAL) /sqrt(2._CUSTOM_REAL)
           else
              forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,il) = s(i1,i2,il) * (sqrt(3._CUSTOM_REAL) - 1._CUSTOM_REAL) /sqrt(2._CUSTOM_REAL)
              forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,ig) = s(i1,i2,ig) * (sqrt(3._CUSTOM_REAL) + 1._CUSTOM_REAL) /sqrt(2._CUSTOM_REAL)
           endif
           il = il+step2
           ig = ig+step2
        enddo
     enddo
  enddo

end subroutine DaubD43D

!=====================================================================================================
! Inverse Daubechies D4 transform from algorithm found at
! http://www.bearcave.com/misl/misl_tech/wavelets/daubechies/index.html
!=====================================================================================================
subroutine iDaubD43D(s,n1,n2,n3)

  integer, intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL), intent(inout) :: s(n1,n2,n3)

  integer :: i,i1,i2,i3,ic,L,il,ig,ngmin,ngmax
  integer :: istep,step_incr,step2,nscale,ng

! Loop over the 3 dimensions
  do ic = 1,3

     ! Loop over the scales
     if (ic==1) then
        nscale = int(log(real(n1,CUSTOM_REAL))/log(2._CUSTOM_REAL))
        L = n1
     else if (ic==2) then
        nscale = int(log(real(n2,CUSTOM_REAL))/log(2._CUSTOM_REAL))
        L = n2
     else
        nscale = int(log(real(n3,CUSTOM_REAL))/log(2._CUSTOM_REAL))
        L = n3
     endif
     do istep = nscale,1,-1
        step_incr = 2**istep
        ngmin = step_incr/2+1
        ngmax = ngmin+int((L-ngmin)/step_incr)*step_incr
        ng = (ngmax-ngmin)/step_incr+1
        step2 = step_incr

        !--------------- Normalization
        ig = ngmin
        il = 1
        do i = 1,ng
           if (ic==1) then
              forall(i2 = 1:n2, i3 = 1:n3) s(il,i2,i3) = s(il,i2,i3) * (sqrt(3._CUSTOM_REAL) + 1._CUSTOM_REAL) /sqrt(2._CUSTOM_REAL)
              forall(i2 = 1:n2, i3 = 1:n3) s(ig,i2,i3) = s(ig,i2,i3) * (sqrt(3._CUSTOM_REAL) - 1._CUSTOM_REAL) /sqrt(2._CUSTOM_REAL)
           else if (ic==2) then
              forall(i1 = 1:n1, i3 = 1:n3) s(i1,il,i3) = s(i1,il,i3) * (sqrt(3._CUSTOM_REAL) + 1._CUSTOM_REAL) /sqrt(2._CUSTOM_REAL)
              forall(i1 = 1:n1, i3 = 1:n3) s(i1,ig,i3) = s(i1,ig,i3) * (sqrt(3._CUSTOM_REAL) - 1._CUSTOM_REAL) /sqrt(2._CUSTOM_REAL)
           else
              forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,il) = s(i1,i2,il) * (sqrt(3._CUSTOM_REAL) + 1._CUSTOM_REAL) /sqrt(2._CUSTOM_REAL)
              forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,ig) = s(i1,i2,ig) * (sqrt(3._CUSTOM_REAL) - 1._CUSTOM_REAL) /sqrt(2._CUSTOM_REAL)
           endif
           il = il+step2
           ig = ig+step2
        enddo


        !------------- Update 2

        ig = ngmin+(ng-2)*step2
        il = 1+(ng-2)*step2
        do i = 1,ng-1
           if (ic==1) then
              forall(i2 = 1:n2, i3 = 1:n3) s(il,i2,i3) = s(il,i2,i3)+s(ig+step2,i2,i3)
           else if (ic==2) then
              forall(i1 = 1:n1, i3 = 1:n3) s(i1,il,i3) = s(i1,il,i3)+s(i1,ig+step2,i3)
           else
              forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,il) = s(i1,i2,il)+s(i1,i2,ig+step2)
           endif
           il = il-step2
           ig = ig-step2
        enddo

        !-------------- Predict

        ig = ngmin + (ng)*step2
        il = 1 + (ng)*step2
        do i = 1,ng-1

           il = il-step2
           ig = ig-step2
           if (ic==1) then
              forall(i2 = 1:n2, i3 = 1:n3) s(ig,i2,i3) = s(ig,i2,i3) + s(il,i2,i3)*sqrt(3._CUSTOM_REAL)/4._CUSTOM_REAL &
              + s(il-step2,i2,i3)*(sqrt(3._CUSTOM_REAL)-2._CUSTOM_REAL)/4._CUSTOM_REAL
           else if (ic==2) then
              forall(i1 = 1:n1, i3 = 1:n3) s(i1,ig,i3) = s(i1,ig,i3) + s(i1,il,i3) *sqrt(3._CUSTOM_REAL)/4._CUSTOM_REAL &
              + s(i1,il-step2,i3)*(sqrt(3._CUSTOM_REAL)-2._CUSTOM_REAL)/4._CUSTOM_REAL
           else
              forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,ig) = s(i1,i2,ig) + s(i1,i2,il) *sqrt(3._CUSTOM_REAL)/4._CUSTOM_REAL &
              + s(i1,i2,il-step2)*(sqrt(3._CUSTOM_REAL)-2._CUSTOM_REAL)/4._CUSTOM_REAL
           endif

        enddo

        il=1
        ig=ngmin
        if (ic==1) then
           forall(i2 = 1:n2, i3 = 1:n3) s(ig,i2,i3) = s(ig,i2,i3) + s(il,i2,i3)*sqrt(3._CUSTOM_REAL)/4._CUSTOM_REAL
        else if (ic==2) then
           forall(i1 = 1:n1, i3 = 1:n3) s(i1,ig,i3) = s(i1,ig,i3) + s(i1,il,i3) *sqrt(3._CUSTOM_REAL)/4._CUSTOM_REAL
        else
           forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,ig) = s(i1,i2,ig) + s(i1,i2,il) *sqrt(3._CUSTOM_REAL)/4._CUSTOM_REAL
        endif

        !--------------Update 1
        ig = ngmin
        il = 1
        do i = 1,ng
           if (ic==1) then
              forall(i2 = 1:n2, i3 = 1:n3) s(il,i2,i3) = s(il,i2,i3)-s(ig,i2,i3)*sqrt(3._CUSTOM_REAL)
           else if (ic==2) then
              forall(i1 = 1:n1, i3 = 1:n3) s(i1,il,i3) = s(i1,il,i3)-s(i1,ig,i3)*sqrt(3._CUSTOM_REAL)
           else
              forall(i1 = 1:n1, i2 = 1:n2) s(i1,i2,il) = s(i1,i2,il)-s(i1,i2,ig)*sqrt(3._CUSTOM_REAL)
           endif
           il = il+step2
           ig = ig+step2
        enddo

     enddo
  enddo

end subroutine iDaubD43D

end module wavelet_transform
