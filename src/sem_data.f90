!
! sem_data.f90 - Shared Data Module for SEM_08 
! 
! Replaces COMMON blocks with modern Fortran 2008 module
! Created: August 16, 2025
!
! Replaces:
!   common /index/ iu,iv,ip,iom        -> Field indices
!   common /jacpar/alp,bet,rv          -> Jacobi polynomial parameters
!
module sem_data
  implicit none
  
  ! Field indices (replaces common /index/)
  ! These define the position of each variable in the solution vector
  integer, parameter :: iu = 1     ! U-velocity component index
  integer, parameter :: iv = 2     ! V-velocity component index  
  integer, parameter :: ip = 3     ! Pressure component index
  integer, parameter :: iom = 4    ! Vorticity component index
  
  ! Jacobi polynomial parameters (replaces common /jacpar/)
  ! These are working variables used in spectral basis function computations
  real(8) :: alp     ! Alpha parameter for Jacobi polynomials
  real(8) :: bet     ! Beta parameter for Jacobi polynomials
  real(8) :: rv      ! Working variable for polynomial evaluation
  
end module sem_data
