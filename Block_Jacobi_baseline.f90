! =============================================================================
! main_timing_demo_fixed.f90
! Compares the performance of PCG with no preconditioner, a direct inverse
! Block Jacobi preconditioner, and a CGS-based Block Jacobi preconditioner.
!
! To compile: gfortran -fdefault-real-8 -O3 -mtune=native main_timing_demo_fixed.f90 -o solver -llapack -lblas
! To run:     ./solver
! =============================================================================
program timing_demo
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  integer, parameter :: n = 8
  integer, parameter :: block_size = 4
  real(real64) :: A(n, n), b(n), x(n)
  real(real64) :: time_start, time_end
  integer :: i

  ! Abstract interface for the outer preconditioner
  abstract interface
    subroutine preconditioner_interface(A_matrix, r, z)
      import :: real64, n, block_size
      implicit none
      real(real64), intent(in)  :: A_matrix(n, n)
      real(real64), intent(in)  :: r(n)
      real(real64), intent(out) :: z(n)
    end subroutine preconditioner_interface
  end interface

  ! --- Define an 8x8 ill-conditioned system ---
  A = 0.0_real64
  ! Block 1
  do i = 1, block_size
    A(i, i) = 2.0_real64
    if (i > 1) A(i, i-1) = -1.0_real64
    if (i < block_size) A(i, i+1) = -1.0_real64
  end do
  ! Block 2 (scaled to create ill-conditioning between blocks)
  do i = block_size + 1, n
    A(i, i) = 200.0_real64
    if (i > block_size + 1) A(i, i-1) = -100.0_real64
    if (i < n) A(i, i+1) = -100.0_real64
  end do

  ! Define b such that the exact solution is x = (1, 1, ..., 1)
  x = 1.0_real64
  b = matmul(A, x)

  write(*, '(A)') "Comparing Preconditioner Performance on an 8x8 System"
  write(*, '(80("="))')

  ! ===========================================================================
  ! Scenario 1: PCG with No Preconditioning
  ! ===========================================================================
  write(*, '(A)') "SCENARIO 1: PCG with No Preconditioner"
  x = 0.0_real64
  call cpu_time(time_start)
  call pcg_solver(A, b, x, apply_identity_preconditioner)
  call cpu_time(time_end)
  write(*, '(A, F10.6, A)') "Elapsed Time: ", time_end - time_start, " seconds"
  write(*, '(80("="))')

  ! ===========================================================================
  ! Scenario 2: PCG with Block Jacobi (Direct 4x4 Inverse)
  ! ===========================================================================
  write(*, '(A)') "SCENARIO 2: PCG with Block Jacobi (Direct Inverse)"
  x = 0.0_real64
  call cpu_time(time_start)
  call pcg_solver(A, b, x, apply_block_jacobi_direct)
  call cpu_time(time_end)
  write(*, '(A, F10.6, A)') "Elapsed Time: ", time_end - time_start, " seconds"
  write(*, '(80("="))')

  ! ===========================================================================
  ! Scenario 3: PCG with Block Jacobi (Inner CGS Solver)
  ! ===========================================================================
  write(*, '(A)') "SCENARIO 3: PCG with Block Jacobi (Inner CGS Solver)"
  x = 0.0_real64
  call cpu_time(time_start)
  call pcg_solver(A, b, x, apply_block_jacobi_cgs)
  call cpu_time(time_end)
  write(*, '(A, F10.6, A)') "Elapsed Time: ", time_end - time_start, " seconds"
  write(*, '(80("="))')


contains

  ! ===========================================================================
  ! PRECONDITIONER 1: IDENTITY (simulates NO preconditioning)
  ! ===========================================================================
  subroutine apply_identity_preconditioner(A_matrix, r, z)
    implicit none
    real(real64), intent(in)  :: A_matrix(n, n)
    real(real64), intent(in)  :: r(n)
    real(real64), intent(out) :: z(n)
    z = r
  end subroutine apply_identity_preconditioner

  ! ===========================================================================
  ! PRECONDITIONER 2: DIRECT INVERSE
  ! ===========================================================================
  subroutine apply_block_jacobi_direct(A_matrix, r, z)
    implicit none
    real(real64), intent(in)  :: A_matrix(n, n)
    real(real64), intent(in)  :: r(n)
    real(real64), intent(out) :: z(n)
    
    real(real64) :: block_A(block_size, block_size)
    integer :: i, offset
    
    do i = 1, n / block_size
      offset = (i - 1) * block_size
      block_A = A_matrix(offset+1 : offset+block_size, offset+1 : offset+block_size)
      z(offset+1 : offset+block_size) = matmul(invert_4x4(block_A), r(offset+1 : offset+block_size))
    end do
  end subroutine apply_block_jacobi_direct

  ! ===========================================================================
  ! PRECONDITIONER 3: INNER CGS SOLVER
  ! ===========================================================================
  subroutine apply_block_jacobi_cgs(A_matrix, r, z)
    implicit none
    real(real64), intent(in)  :: A_matrix(n, n)
    real(real64), intent(in)  :: r(n)
    real(real64), intent(out) :: z(n)
    
    real(real64) :: block_A(block_size, block_size)
    real(real64) :: block_r(block_size), block_z(block_size)
    integer :: i, offset
    
    do i = 1, n / block_size
      offset = (i - 1) * block_size
      block_A = A_matrix(offset+1 : offset+block_size, offset+1 : offset+block_size)
      block_r = r(offset+1 : offset+block_size)
      call cgs_solver(block_A, block_r, block_z)
      z(offset+1 : offset+block_size) = block_z
    end do
  end subroutine apply_block_jacobi_cgs

  ! ===========================================================================
  ! HELPER: Computes the inverse of a 4x4 matrix
  ! ===========================================================================
  function invert_4x4(A_in) result(A_inv)
      implicit none
      real(real64), intent(in) :: A_in(4,4)
      real(real64) :: A_inv(4,4)
      real(real64) :: det
      ! Implementation using adjugate method
      A_inv(1, 1) = A_in(2, 2) * (A_in(3, 3) * A_in(4, 4) - A_in(3, 4) * A_in(4, 3)) - A_in(2, 3) * (A_in(3, 2) * A_in(4, 4) - A_in(3, 4) * A_in(4, 2)) + A_in(2, 4) * (A_in(3, 2) * A_in(4, 3) - A_in(3, 3) * A_in(4, 2))
      A_inv(2, 1) = -A_in(2, 1) * (A_in(3, 3) * A_in(4, 4) - A_in(3, 4) * A_in(4, 3)) + A_in(2, 3) * (A_in(3, 1) * A_in(4, 4) - A_in(3, 4) * A_in(4, 1)) - A_in(2, 4) * (A_in(3, 1) * A_in(4, 3) - A_in(3, 3) * A_in(4, 1))
      A_inv(3, 1) = A_in(2, 1) * (A_in(3, 2) * A_in(4, 4) - A_in(3, 4) * A_in(4, 2)) - A_in(2, 2) * (A_in(3, 1) * A_in(4, 4) - A_in(3, 4) * A_in(4, 1)) + A_in(2, 4) * (A_in(3, 1) * A_in(4, 2) - A_in(3, 2) * A_in(4, 1))
      A_inv(4, 1) = -A_in(2, 1) * (A_in(3, 2) * A_in(4, 3) - A_in(3, 3) * A_in(4, 2)) + A_in(2, 2) * (A_in(3, 1) * A_in(4, 3) - A_in(3, 3) * A_in(4, 1)) - A_in(2, 3) * (A_in(3, 1) * A_in(4, 2) - A_in(3, 2) * A_in(4, 1))
      A_inv(1, 2) = -A_in(1, 2) * (A_in(3, 3) * A_in(4, 4) - A_in(3, 4) * A_in(4, 3)) + A_in(1, 3) * (A_in(3, 2) * A_in(4, 4) - A_in(3, 4) * A_in(4, 2)) - A_in(1, 4) * (A_in(3, 2) * A_in(4, 3) - A_in(3, 3) * A_in(4, 2))
      A_inv(2, 2) = A_in(1, 1) * (A_in(3, 3) * A_in(4, 4) - A_in(3, 4) * A_in(4, 3)) - A_in(1, 3) * (A_in(3, 1) * A_in(4, 4) - A_in(3, 4) * A_in(4, 1)) + A_in(1, 4) * (A_in(3, 1) * A_in(4, 3) - A_in(3, 3) * A_in(4, 1))
      A_inv(3, 2) = -A_in(1, 1) * (A_in(3, 2) * A_in(4, 4) - A_in(3, 4) * A_in(4, 2)) + A_in(1, 2) * (A_in(3, 1) * A_in(4, 4) - A_in(3, 4) * A_in(4, 1)) - A_in(1, 4) * (A_in(3, 1) * A_in(4, 2) - A_in(3, 2) * A_in(4, 1))
      A_inv(4, 2) = A_in(1, 1) * (A_in(3, 2) * A_in(4, 3) - A_in(3, 3) * A_in(4, 2)) - A_in(1, 2) * (A_in(3, 1) * A_in(4, 3) - A_in(3, 3) * A_in(4, 1)) + A_in(1, 3) * (A_in(3, 1) * A_in(4, 2) - A_in(3, 2) * A_in(4, 1))
      A_inv(1, 3) = A_in(1, 2) * (A_in(2, 3) * A_in(4, 4) - A_in(2, 4) * A_in(4, 3)) - A_in(1, 3) * (A_in(2, 2) * A_in(4, 4) - A_in(2, 4) * A_in(4, 2)) + A_in(1, 4) * (A_in(2, 2) * A_in(4, 3) - A_in(2, 3) * A_in(4, 2))
      A_inv(2, 3) = -A_in(1, 1) * (A_in(2, 3) * A_in(4, 4) - A_in(2, 4) * A_in(4, 3)) + A_in(1, 3) * (A_in(2, 1) * A_in(4, 4) - A_in(2, 4) * A_in(4, 1)) - A_in(1, 4) * (A_in(2, 1) * A_in(4, 3) - A_in(2, 3) * A_in(4, 1))
      A_inv(3, 3) = A_in(1, 1) * (A_in(2, 2) * A_in(4, 4) - A_in(2, 4) * A_in(4, 2)) - A_in(1, 2) * (A_in(2, 1) * A_in(4, 4) - A_in(2, 4) * A_in(4, 1)) + A_in(1, 4) * (A_in(2, 1) * A_in(4, 2) - A_in(2, 2) * A_in(4, 1))
      A_inv(4, 3) = -A_in(1, 1) * (A_in(2, 2) * A_in(4, 3) - A_in(2, 3) * A_in(4, 2)) + A_in(1, 2) * (A_in(2, 1) * A_in(4, 3) - A_in(2, 3) * A_in(4, 1)) - A_in(1, 3) * (A_in(2, 1) * A_in(4, 2) - A_in(2, 2) * A_in(4, 1))
      A_inv(1, 4) = -A_in(1, 2) * (A_in(2, 3) * A_in(3, 4) - A_in(2, 4) * A_in(3, 3)) + A_in(1, 3) * (A_in(2, 2) * A_in(3, 4) - A_in(2, 4) * A_in(3, 2)) - A_in(1, 4) * (A_in(2, 2) * A_in(3, 3) - A_in(2, 3) * A_in(3, 2))
      A_inv(2, 4) = A_in(1, 1) * (A_in(2, 3) * A_in(3, 4) - A_in(2, 4) * A_in(3, 3)) - A_in(1, 3) * (A_in(2, 1) * A_in(3, 4) - A_in(2, 4) * A_in(3, 1)) + A_in(1, 4) * (A_in(2, 1) * A_in(3, 3) - A_in(2, 3) * A_in(3, 1))
      A_inv(3, 4) = -A_in(1, 1) * (A_in(2, 2) * A_in(3, 4) - A_in(2, 4) * A_in(3, 2)) + A_in(1, 2) * (A_in(2, 1) * A_in(3, 4) - A_in(2, 4) * A_in(3, 1)) - A_in(1, 4) * (A_in(2, 1) * A_in(3, 2) - A_in(2, 2) * A_in(3, 1))
      A_inv(4, 4) = A_in(1, 1) * (A_in(2, 2) * A_in(3, 3) - A_in(2, 3) * A_in(3, 2)) - A_in(1, 2) * (A_in(2, 1) * A_in(3, 3) - A_in(2, 3) * A_in(3, 1)) + A_in(1, 3) * (A_in(2, 1) * A_in(3, 2) - A_in(2, 2) * A_in(3, 1))
      det = A_in(1, 1) * A_inv(1, 1) + A_in(1, 2) * A_inv(2, 1) + A_in(1, 3) * A_inv(3, 1) + A_in(1, 4) * A_inv(4, 1)
      A_inv = A_inv / det
  end function invert_4x4

  ! ===========================================================================
  ! The Inner Solver: CGS
  ! ===========================================================================
  subroutine cgs_solver(A_block, b_block, x_block)
    implicit none
    real(real64), intent(in)    :: A_block(block_size, block_size)
    real(real64), intent(in)    :: b_block(block_size)
    real(real64), intent(inout) :: x_block(block_size)

    real(real64) :: r(block_size), r_hat(block_size), p(block_size)
    real(real64) :: u(block_size), q(block_size), v_hat(block_size)
    real(real64) :: rho, rho_old, alpha, beta
    real(real64), parameter :: inner_tol = 1.0e-6
    integer, parameter :: max_inner_iter = 15
    integer :: iter
    
    x_block = 0.0_real64
    r = b_block - matmul(A_block, x_block)
    r_hat = r
    p = 0.0_real64; u = 0.0_real64; rho_old = 1.0_real64
    
    cgs_loop: do iter = 1, max_inner_iter
        rho = dot_product(r_hat, r)
        if (abs(rho) < 1.0e-20) exit cgs_loop
        
        if (iter == 1) then
            u = r; p = u
        else
            beta = rho / rho_old
            u = r + beta * q
            p = u + beta * (q + beta * p)
        end if
        
        v_hat = matmul(A_block, p)
        alpha = rho / dot_product(r_hat, v_hat)
        q = u - alpha * v_hat
        x_block = x_block + alpha * (u + q)
        r = r - alpha * matmul(A_block, u + q)
        
        if (norm2(r) < inner_tol) exit cgs_loop
        rho_old = rho
    end do cgs_loop
  end subroutine cgs_solver

  ! ===========================================================================
  ! The Outer Solver: PCG
  ! ===========================================================================
  subroutine pcg_solver(A_matrix, b_vec, x_vec, preconditioner)
    implicit none
    real(real64), intent(in)    :: A_matrix(n, n)
    real(real64), intent(in)    :: b_vec(n)
    real(real64), intent(inout) :: x_vec(n)
    procedure(preconditioner_interface) :: preconditioner
    real(real64) :: r(n), p(n), z(n), Ap(n)
    real(real64) :: r_dot_z, r_dot_z_new, alpha, beta
    real(real64), parameter :: outer_tol = 1.0e-12
    integer, parameter :: max_outer_iter = 100
    integer :: iter
    
    r = b_vec - matmul(A_matrix, x_vec)
    call preconditioner(A_matrix, r, z)
    p = z
    r_dot_z = dot_product(r, z)

    write(*, '(A7, A18)') "Iter", "Residual Norm"
    write(*, '(I5, E18.6)') 0, norm2(r)
    
    pcg_loop: do iter = 1, max_outer_iter
      Ap = matmul(A_matrix, p)
      alpha = r_dot_z / dot_product(p, Ap)
      x_vec = x_vec + alpha * p
      r     = r     - alpha * Ap
      if (norm2(r) < outer_tol) then
        write(*, '(I5, E18.6)') iter, norm2(r)
        write(*,'(A)') "CONVERGED"
        return
      end if
      call preconditioner(A_matrix, r, z)
      r_dot_z_new = dot_product(r, z)
      beta = r_dot_z_new / r_dot_z
      p = z + beta * p
      r_dot_z = r_dot_z_new
      write(*, '(I5, E18.6)') iter, norm2(r)
    end do pcg_loop
    write(*, '(A, I5, A)') "FAILED TO CONVERGE in ", max_outer_iter, " iterations."
  end subroutine pcg_solver

end program timing_demo
