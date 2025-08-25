! =============================================================================
! main_timing_demo_fixed.f90
! Compares the performance of PCG with no preconditioner, a direct inverse
! Block Jacobi preconditioner, and a CGS-based Block Jacobi preconditioner.
!
! To compile: gfortran -fdefault-real-8 -O3 -mtune=native Block_Jacobi.f90 -o solver -llapack -lblas
! To run:     ./solver
! =============================================================================
program timing_demo
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  integer, parameter :: block_size = 16
  integer            :: n
  real(real64), allocatable :: A(:, :), b(:), x(:)
  real(real64) :: time_start, time_end
  integer :: i, j, offset

  ! Variables for condition number calculation
  real(real64) :: norm_A, norm_A_inv, cond_num
  real(real64), allocatable :: A_inv(:,:), work(:)
  integer, allocatable :: ipiv(:)
  integer :: lwork, info
  real(real64) :: dlange
  external :: dlange, dgetrf, dgetri

  ! Abstract interface for the outer preconditioner
  abstract interface
    subroutine preconditioner_interface(A_matrix, r, z)
      import :: real64
      real(real64), intent(in)  :: A_matrix(:, :)
      real(real64), intent(in)  :: r(:)
      real(real64), intent(out) :: z(:)
    end subroutine preconditioner_interface
  end interface

  ! Set the size of the system. Must be a multiple of block_size.
  n = 256
  if (mod(n, block_size) /= 0) then
    write(*,*) "Error: n must be a multiple of block_size (4)."
    stop
  end if

  allocate(A(n, n), b(n), x(n))

  ! --- Define a 64x64 matrix with a target condition number of ~10^15 ---
  A = 0.0_real64
  do i = 1, n
    A(i, i) = 10.0_real64**(15.0_real64 * (dble(i) - 1.0_real64) / (dble(n) - 1.0_real64))
  end do
  ! Add small off-diagonal elements to maintain symmetry
  do i = 1, n
    do j = i + 1, n
      A(i, j) = 1.0e-4_real64
      A(j, i) = 1.0e-4_real64
    end do
  end do

  ! --- Calculate Condition Number ---
  A_inv = A
  allocate(ipiv(n), work(n))
  norm_A = dlange('I', n, n, A, n, work)
  call dgetrf(n, n, A_inv, n, ipiv, info)
  if (info /= 0) then
    write(*,*) "Error: LU factorization failed."
    stop
  end if
  lwork = n
  call dgetri(n, A_inv, n, ipiv, work, lwork, info)
  if (info /= 0) then
    write(*,*) "Error: Matrix inversion failed."
    stop
  end if
  norm_A_inv = dlange('I', n, n, A_inv, n, work)
  cond_num = norm_A * norm_A_inv
  write(*, '(A, E12.5)') "Condition Number (inf-norm): ", cond_num
  deallocate(ipiv, work)

  ! Define b such that the exact solution is x = (1, 1, ..., 1)
  x = 1.0_real64
  b = matmul(A, x)

  write(*, '(A, I0, A, I0, A)') "Comparing Preconditioner Performance on an ", n, "x", n, " System"
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
  ! Scenario 2: PCG with Block Jacobi (Direct Inverse)
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
  ! The Outer Solver: PCG
  ! ===========================================================================
  subroutine pcg_solver(A_matrix, b_vec, x_vec, preconditioner)
    real(real64), intent(in)    :: A_matrix(:, :)
    real(real64), intent(in)    :: b_vec(:)
    real(real64), intent(inout) :: x_vec(:)
    procedure(preconditioner_interface) :: preconditioner
    
    integer :: n_outer
    real(real64), allocatable :: r(:), p(:), z(:), Ap(:)
    real(real64) :: r_dot_z, r_dot_z_new, alpha, beta
    real(real64), parameter :: outer_tol = 1.0e-12, epsilon = 1.0e-30
    integer, parameter :: max_outer_iter = 1000
    integer :: iter
    
    n_outer = size(b_vec)
    allocate(r(n_outer), p(n_outer), z(n_outer), Ap(n_outer))
    
    r = b_vec - matmul(A_matrix, x_vec)
    call preconditioner(A_matrix, r, z)
    p = z
    r_dot_z = dot_product(r, z)

    write(*, '(A7, A18)') "Iter", "Residual Norm"
    write(*, '(I5, E18.6)') 0, norm2(r)
    
    pcg_loop: do iter = 1, max_outer_iter
      Ap = matmul(A_matrix, p)
      alpha = r_dot_z / (dot_product(p, Ap) + epsilon)
      x_vec = x_vec + alpha * p
      r     = r     - alpha * Ap
      if (norm2(r) < outer_tol) then
        write(*, '(I5, E18.6)') iter, norm2(r)
        write(*,'(A)') "CONVERGED"
        deallocate(r, p, z, Ap)
        return
      end if
      call preconditioner(A_matrix, r, z)
      r_dot_z_new = dot_product(r, z)
      beta = r_dot_z_new / (r_dot_z + epsilon)
      p = z + beta * p
      r_dot_z = r_dot_z_new
      write(*, '(I5, E18.6)') iter, norm2(r)
    end do pcg_loop
    write(*, '(A, I5, A)') "FAILED TO CONVERGE in ", max_outer_iter, " iterations."
    deallocate(r, p, z, Ap)
  end subroutine pcg_solver

  ! ===========================================================================
  ! PRECONDITIONER 1: IDENTITY (simulates NO preconditioning)
  ! ===========================================================================
  subroutine apply_identity_preconditioner(A_matrix, r, z)
    real(real64), intent(in)  :: A_matrix(:, :)
    real(real64), intent(in)  :: r(:)
    real(real64), intent(out) :: z(:)
    z = r
  end subroutine apply_identity_preconditioner

  ! ===========================================================================
  ! PRECONDITIONER 2: DIRECT INVERSE
  ! ===========================================================================
  subroutine apply_block_jacobi_direct(A_matrix, r, z)
    real(real64), intent(in)  :: A_matrix(:, :)
    real(real64), intent(in)  :: r(:)
    real(real64), intent(out) :: z(:)
    
    integer :: n_jacobi, i, offset, info
    real(real64) :: block_A(block_size, block_size)
    real(real64) :: block_z(block_size)
    integer      :: ipiv(block_size)
    
    n_jacobi = size(r)
    
    do i = 1, n_jacobi / block_size
      offset = (i - 1) * block_size
      block_A = A_matrix(offset+1 : offset+block_size, offset+1 : offset+block_size)
      
      ! Add a small regularization term to the diagonal for stability
      do j = 1, block_size
          block_A(j, j) = block_A(j, j) + 1.0e-9_real64
      end do
      
      block_z = r(offset+1 : offset+block_size)
      
      ! Solve the block system A_ii * z_i = r_i using LAPACK
      call dgesv(block_size, 1, block_A, block_size, ipiv, block_z, block_size, info)
      if (info /= 0) then
        write(*,*) "WARNING: dgesv failed for block ", i, " with info code ", info
        ! If solver fails, return the original residual block (no preconditioning for this block)
        z(offset+1 : offset+block_size) = r(offset+1 : offset+block_size)
      else
        z(offset+1 : offset+block_size) = block_z
      end if
    end do
  end subroutine apply_block_jacobi_direct

  ! ===========================================================================
  ! PRECONDITIONER 3: INNER CGS SOLVER
  ! ===========================================================================
  subroutine apply_block_jacobi_cgs(A_matrix, r, z)
    real(real64), intent(in)  :: A_matrix(:, :)
    real(real64), intent(in)  :: r(:)
    real(real64), intent(out) :: z(:)
    
    integer :: n_cgs, i, j, offset
    real(real64) :: block_A(block_size, block_size)
    real(real64) :: block_r(block_size), block_z(block_size)
    real(real64) :: diag_inv(block_size)
    real(real64), parameter :: epsilon = 1.0e-30
    
    n_cgs = size(r)
    
    do i = 1, n_cgs / block_size
      offset = (i - 1) * block_size
      block_A = A_matrix(offset+1 : offset+block_size, offset+1 : offset+block_size)
      block_r = r(offset+1 : offset+block_size)
      
      ! Add a small regularization term to the diagonal for stability
      do j = 1, block_size
          block_A(j, j) = block_A(j, j) + 1.0e-9_real64
      end do
      
      ! --- Diagonal Preconditioning ---
      ! 1. Extract and invert the diagonal of the block
      do j = 1, block_size
        diag_inv(j) = 1.0_real64 / (block_A(j, j) + epsilon)
      end do
      
      ! 2. Precondition the block matrix and vector
      do j = 1, block_size
        block_A(j, :) = block_A(j, :) * diag_inv(j)
      end do
      block_r = block_r * diag_inv
      ! ------------------------------------

      call cgs_solver(block_A, block_r, block_z)
      z(offset+1 : offset+block_size) = block_z
    end do
  end subroutine apply_block_jacobi_cgs

  ! ===========================================================================
  ! The Inner Solver: CGS
  ! ===========================================================================
  subroutine cgs_solver(A_block, b_block, x_block)
    real(real64), intent(in)    :: A_block(:, :)
    real(real64), intent(in)    :: b_block(:)
    real(real64), intent(inout) :: x_block(:)

    integer :: n_inner
    real(real64), allocatable :: r(:), r_hat(:), p(:)
    real(real64), allocatable :: u(:), q(:), v_hat(:)
    real(real64) :: rho, rho_old, alpha, beta
    real(real64), parameter :: inner_tol = 1.0e-10, epsilon = 1.0e-30
    integer, parameter :: max_inner_iter = 50
    integer :: iter
    
    n_inner = size(b_block)
    allocate(r(n_inner), r_hat(n_inner), p(n_inner), u(n_inner), q(n_inner), v_hat(n_inner))

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
            beta = rho / (rho_old + epsilon)
            u = r + beta * q
            p = u + beta * (q + beta * p)
        end if
        
        v_hat = matmul(A_block, p)
        alpha = rho / (dot_product(r_hat, v_hat) + epsilon)
        q = u - alpha * v_hat
        x_block = x_block + alpha * (u + q)
        r = r - alpha * matmul(A_block, u + q)
        
        if (norm2(r) < inner_tol) exit cgs_loop
        rho_old = rho
    end do cgs_loop
    deallocate(r, r_hat, p, u, q, v_hat)
  end subroutine cgs_solver

end program timing_demo