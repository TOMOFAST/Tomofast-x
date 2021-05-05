
!========================================================================
!
!                    T O M O F A S T X  Version 1.0
!                  ----------------------------------
!
!              Main authors: Vitaliy Ogarko, Roland Martin,
!                   Jeremie Giraud, Dimitri Komatitsch.
! CNRS, France, and University of Western Australia.
! (c) CNRS, France, and University of Western Australia. January 2018
!
! This software is a computer program whose purpose is to perform
! capacitance, gravity, magnetic, or joint gravity and magnetic tomography.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

!================================================================================================
! Functions to generate noise.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!================================================================================================
module noise

  use global_typedefs

  implicit none

  private

  public :: init_random_seed
  public :: generate_Gaussian_noise
  public :: generate_uniform_noise

contains

!==================================================================================
! Initializes a random seed from the system clock at every run (fortran 95 code).
! http://stackoverflow.com/questions/18754438/generating-random-numbers-in-a-fortran-module
!==================================================================================
subroutine init_random_seed()
  integer :: i, n, clock
  integer, allocatable :: seed(:)

  call RANDOM_SEED(size = n)
  allocate(seed(n))

  call SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call RANDOM_SEED(PUT = seed)

  deallocate(seed)

end subroutine init_random_seed

!============================================================================
! Generates Gaussian noise.
!============================================================================
function generate_Gaussian_noise(sigma, mu) result(noise)
  double precision, intent(in) :: sigma, mu

  double precision :: rand1, rand2
  double precision :: gaussian1, gaussian2
  double precision :: noise

  call random_number(rand1)
  call random_number(rand2)

  ! Uncorrelated Gaussian (normal) random numbers (by Box & Muller method).
  gaussian1 = sqrt(- 2.d0 * log(rand1)) * cos(2.d0 * PI * rand2)
  gaussian2 = sqrt(- 2.d0 * log(rand1)) * sin(2.d0 * PI * rand2)

  noise = gaussian1 * sigma + mu

end function generate_Gaussian_noise

!============================================================================
! Generates uniform noise.
!============================================================================
function generate_uniform_noise(scale) result(noise)
  double precision, intent(in) :: scale

  double precision :: rand
  double precision :: noise

  call random_number(rand)

  noise = (rand - 0.5d0) * scale

end function generate_uniform_noise

end module noise
