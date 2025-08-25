!
! test_derv_2D_enhanced.f90 - Enhanced test for 2D derivative matrices
!
! Tests multiple analytical functions and provides convergence analysis
!
program test_derv_2D_enhanced
    implicit none
    
    ! Parameters
    integer, parameter :: max_n = 16
    integer, parameter :: max_dim = max_n + 1
    integer, parameter :: max_nodes_2d = max_dim * max_dim
    integer, parameter :: num_tests = 3
    
    ! Variables
    integer :: n, test_case
    integer :: nnodes_1d, nnodes_2d
    real(8) :: x_1d(0:max_dim), w_1d(0:max_dim)
    real(8) :: d_1d(0:max_dim, 0:max_dim)
    real(8) :: dx_2d(max_nodes_2d, max_nodes_2d)
    real(8) :: dy_2d(max_nodes_2d, max_nodes_2d)
    real(8) :: x_2d(max_nodes_2d), y_2d(max_nodes_2d)
    real(8) :: f_vals(max_nodes_2d), fx_exact(max_nodes_2d), fy_exact(max_nodes_2d)
    real(8) :: fx_computed(max_nodes_2d), fy_computed(max_nodes_2d)
    real(8) :: max_error_x, max_error_y
    character(30) :: test_names(num_tests)
    
    ! Test function names
    test_names(1) = 'sin(πx)cos(πy)'
    test_names(2) = 'polynomial: x³y² + xy³'
    test_names(3) = 'exp(x+y)'
    
    write(*,*) '========================================'
    write(*,*) 'Enhanced 2D Derivative Matrix Testing'
    write(*,*) '========================================'
    write(*,*)
    
    ! Test each function
    do test_case = 1, num_tests
        write(*,*) 'Test Function: ', trim(test_names(test_case))
        write(*,'(A8,A15,A15,A12)') 'Order N', 'Max Error dx', 'Max Error dy', 'DOF (N+1)²'
        write(*,'(A8,A15,A15,A12)') '-------', '------------', '------------', '-----------'
        
        do n = 2, 12, 2
            nnodes_1d = n + 1
            nnodes_2d = nnodes_1d * nnodes_1d
            
            ! Setup 1D grid and derivative matrix
            call setup_1d_grid(n, x_1d, w_1d, d_1d, max_dim)
            
            ! Create 2D grid and derivative matrices
            call create_2d_grid(n, x_1d, x_2d, y_2d, max_dim, max_nodes_2d)
            call create_2d_derivative_matrices(n, d_1d, dx_2d, dy_2d, max_dim, max_nodes_2d)
            
            ! Test specific function
            select case (test_case)
            case (1)
                call test_function_sincos(nnodes_2d, x_2d, y_2d, f_vals, fx_exact, fy_exact, max_nodes_2d)
            case (2)
                call test_function_polynomial(nnodes_2d, x_2d, y_2d, f_vals, fx_exact, fy_exact, max_nodes_2d)
            case (3)
                call test_function_exponential(nnodes_2d, x_2d, y_2d, f_vals, fx_exact, fy_exact, max_nodes_2d)
            end select
            
            ! Compute derivatives
            call apply_derivative_matrix(nnodes_2d, dx_2d, f_vals, fx_computed, max_nodes_2d)
            call apply_derivative_matrix(nnodes_2d, dy_2d, f_vals, fy_computed, max_nodes_2d)
            
            ! Compute errors
            call compute_errors(nnodes_2d, fx_exact, fy_exact, fx_computed, fy_computed, &
                               max_error_x, max_error_y, max_nodes_2d)
            
            write(*,'(I8,2E15.6,I12)') n, max_error_x, max_error_y, nnodes_2d
        enddo
        write(*,*)
    enddo
    
    ! Test convergence rate
    write(*,*) 'Convergence Rate Analysis for sin(πx)cos(πy):'
    write(*,'(A8,A15,A15,A12,A12)') 'Order N', 'Error dx', 'Rate dx', 'Error dy', 'Rate dy'
    write(*,'(A8,A15,A15,A12,A12)') '-------', '--------', '-------', '--------', '-------'
    
    call convergence_analysis()
    
end program test_derv_2D_enhanced

!**********************************************************************
subroutine test_function_polynomial(nnodes, x, y, f, fx, fy, ndim)
!**********************************************************************
    implicit none
    integer, intent(in) :: nnodes, ndim
    real(8), intent(in) :: x(ndim), y(ndim)
    real(8), intent(out) :: f(ndim), fx(ndim), fy(ndim)
    integer :: i
    
    ! Test function: f(x,y) = x³y² + xy³
    ! df/dx = 3x²y² + y³
    ! df/dy = 2x³y + 3xy²
    
    do i = 1, nnodes
        f(i) = x(i)**3 * y(i)**2 + x(i) * y(i)**3
        fx(i) = 3.0d0 * x(i)**2 * y(i)**2 + y(i)**3
        fy(i) = 2.0d0 * x(i)**3 * y(i) + 3.0d0 * x(i) * y(i)**2
    enddo
    
end subroutine test_function_polynomial

!**********************************************************************
subroutine test_function_exponential(nnodes, x, y, f, fx, fy, ndim)
!**********************************************************************
    implicit none
    integer, intent(in) :: nnodes, ndim
    real(8), intent(in) :: x(ndim), y(ndim)
    real(8), intent(out) :: f(ndim), fx(ndim), fy(ndim)
    integer :: i
    
    ! Test function: f(x,y) = exp(x + y)
    ! df/dx = exp(x + y)
    ! df/dy = exp(x + y)
    
    do i = 1, nnodes
        f(i) = exp(x(i) + y(i))
        fx(i) = exp(x(i) + y(i))
        fy(i) = exp(x(i) + y(i))
    enddo
    
end subroutine test_function_exponential

!**********************************************************************
subroutine convergence_analysis()
!**********************************************************************
    implicit none
    integer, parameter :: max_n = 16, max_dim = max_n + 1, max_nodes_2d = max_dim * max_dim
    integer :: n, nnodes_1d, nnodes_2d
    real(8) :: x_1d(0:max_dim), w_1d(0:max_dim), d_1d(0:max_dim, 0:max_dim)
    real(8) :: dx_2d(max_nodes_2d, max_nodes_2d), dy_2d(max_nodes_2d, max_nodes_2d)
    real(8) :: x_2d(max_nodes_2d), y_2d(max_nodes_2d)
    real(8) :: f_vals(max_nodes_2d), fx_exact(max_nodes_2d), fy_exact(max_nodes_2d)
    real(8) :: fx_computed(max_nodes_2d), fy_computed(max_nodes_2d)
    real(8) :: errors_x(10), errors_y(10), rate_x, rate_y
    integer :: i, n_vals(10) = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
    
    ! Compute errors for different orders
    do i = 1, 6  ! Test up to N=12
        n = n_vals(i)
        nnodes_1d = n + 1
        nnodes_2d = nnodes_1d * nnodes_1d
        
        call setup_1d_grid(n, x_1d, w_1d, d_1d, max_dim)
        call create_2d_grid(n, x_1d, x_2d, y_2d, max_dim, max_nodes_2d)
        call create_2d_derivative_matrices(n, d_1d, dx_2d, dy_2d, max_dim, max_nodes_2d)
        call test_function_sincos(nnodes_2d, x_2d, y_2d, f_vals, fx_exact, fy_exact, max_nodes_2d)
        call apply_derivative_matrix(nnodes_2d, dx_2d, f_vals, fx_computed, max_nodes_2d)
        call apply_derivative_matrix(nnodes_2d, dy_2d, f_vals, fy_computed, max_nodes_2d)
        call compute_errors(nnodes_2d, fx_exact, fy_exact, fx_computed, fy_computed, &
                           errors_x(i), errors_y(i), max_nodes_2d)
        
        if (i > 1) then
            rate_x = log(errors_x(i-1)/errors_x(i)) / log(real(n_vals(i))/real(n_vals(i-1)))
            rate_y = log(errors_y(i-1)/errors_y(i)) / log(real(n_vals(i))/real(n_vals(i-1)))
            write(*,'(I8,2E15.6,2F12.2)') n, errors_x(i), rate_x, errors_y(i), rate_y
        else
            write(*,'(I8,2E15.6,2A12)') n, errors_x(i), errors_y(i), '     ---', '     ---'
        endif
    enddo
    
end subroutine convergence_analysis

! Include all the subroutines from the original test_derv_2D.f90
! (I'll reference them to avoid duplication in this context)

!**********************************************************************
subroutine setup_1d_grid(n, x, w, d, ndim)
!**********************************************************************
    implicit none
    integer, intent(in) :: n, ndim
    real(8), intent(out) :: x(0:ndim), w(0:ndim), d(0:ndim, 0:ndim)
    call jacobl(n, 0.0d0, 0.0d0, x, ndim)
    call quad(n, x, w, ndim)
    call derv(n, x, d, ndim)
end subroutine setup_1d_grid

!**********************************************************************
subroutine create_2d_grid(n, x_1d, x_2d, y_2d, ndim_1d, ndim_2d)
!**********************************************************************
    implicit none
    integer, intent(in) :: n, ndim_1d, ndim_2d
    real(8), intent(in) :: x_1d(0:ndim_1d)
    real(8), intent(out) :: x_2d(ndim_2d), y_2d(ndim_2d)
    integer :: i, j, node_idx
    
    node_idx = 1
    do j = 0, n
        do i = 0, n
            x_2d(node_idx) = x_1d(i)
            y_2d(node_idx) = x_1d(j)
            node_idx = node_idx + 1
        enddo
    enddo
end subroutine create_2d_grid

!**********************************************************************
subroutine create_2d_derivative_matrices(n, d_1d, dx_2d, dy_2d, ndim_1d, ndim_2d)
!**********************************************************************
!
! Purpose:
!   Constructs 2D derivative matrices in x and y directions using tensor 
!   products of 1D differentiation matrices for spectral element methods
!
! Mathematical Background:
!   For a 2D spectral element with polynomial order N in both directions:
!   - Total DOF = (N+1)^2 nodes arranged in tensor product grid
!   - Node indexing: global_node = i + j*(N+1) where i,j ∈ [0,N]
!   - 2D derivatives computed via tensor products:
!     * ∂/∂x: D_x = D_1D ⊗ I_y (differentiate in x, identity in y)
!     * ∂/∂y: D_y = I_x ⊗ D_1D (identity in x, differentiate in y)
!
! Implementation:
!   - For x-derivative: Apply 1D derivative matrix along x-direction 
!     while keeping y-coordinate constant
!   - For y-derivative: Apply 1D derivative matrix along y-direction 
!     while keeping x-coordinate constant
!
! Matrix Structure:
!   - Size: (N+1)^2 × (N+1)^2 for each direction
!   - Sparsity: Each row has exactly (N+1) non-zero entries
!   - Storage: Dense format used here; sparse storage possible for large N
!
! Arguments:
!   n        (in)  : Polynomial order (number of intervals = N+1 points)
!   d_1d     (in)  : 1D differentiation matrix from derv() subroutine
!   dx_2d    (out) : 2D derivative matrix for ∂/∂x direction
!   dy_2d    (out) : 2D derivative matrix for ∂/∂y direction  
!   ndim_1d  (in)  : Dimension of 1D arrays (≥ N+1)
!   ndim_2d  (in)  : Dimension of 2D arrays (≥ (N+1)^2)
!
! Example Usage:
!   To compute derivatives: df_dx = matmul(dx_2d, f_2d)
!                          df_dy = matmul(dy_2d, f_2d)
!   where f_2d contains function values at 2D tensor product grid points
!
! Notes:
!   - Grid points must be ordered as: node(i,j) = i + j*(N+1)
!   - Compatible with Legendre-Gauss-Lobatto (LGL) collocation points
!   - Achieves spectral accuracy for smooth functions
!
! Author: Generated for LSSEM_F90 project
! Date: August 2025
!
!**********************************************************************
    implicit none
    integer, intent(in) :: n, ndim_1d, ndim_2d
    real(8), intent(in) :: d_1d(0:ndim_1d, 0:ndim_1d)
    real(8), intent(out) :: dx_2d(ndim_2d, ndim_2d), dy_2d(ndim_2d, ndim_2d)
    integer :: i, j, k, l, row_idx, col_idx
    
    dx_2d = 0.0d0
    dy_2d = 0.0d0
    
    do j = 0, n
        do i = 0, n
            row_idx = i + 1 + j * (n + 1)
            do k = 0, n
                col_idx = k + 1 + j * (n + 1)
                dx_2d(row_idx, col_idx) = d_1d(i, k)
            enddo
        enddo
    enddo
    
    do j = 0, n
        do i = 0, n
            row_idx = i + 1 + j * (n + 1)
            do l = 0, n
                col_idx = i + 1 + l * (n + 1)
                dy_2d(row_idx, col_idx) = d_1d(j, l)
            enddo
        enddo
    enddo
end subroutine create_2d_derivative_matrices

!**********************************************************************
subroutine test_function_sincos(nnodes, x, y, f, fx, fy, ndim)
!**********************************************************************
    implicit none
    integer, intent(in) :: nnodes, ndim
    real(8), intent(in) :: x(ndim), y(ndim)
    real(8), intent(out) :: f(ndim), fx(ndim), fy(ndim)
    real(8), parameter :: pi = 3.14159265358979323846d0
    integer :: i
    
    do i = 1, nnodes
        f(i) = sin(pi * x(i)) * cos(pi * y(i))
        fx(i) = pi * cos(pi * x(i)) * cos(pi * y(i))
        fy(i) = -pi * sin(pi * x(i)) * sin(pi * y(i))
    enddo
end subroutine test_function_sincos

!**********************************************************************
subroutine apply_derivative_matrix(nnodes, d_matrix, f, df, ndim)
!**********************************************************************
    implicit none
    integer, intent(in) :: nnodes, ndim
    real(8), intent(in) :: d_matrix(ndim, ndim), f(ndim)
    real(8), intent(out) :: df(ndim)
    integer :: i, j
    
    do i = 1, nnodes
        df(i) = 0.0d0
        do j = 1, nnodes
            df(i) = df(i) + d_matrix(i, j) * f(j)
        enddo
    enddo
end subroutine apply_derivative_matrix

!**********************************************************************
subroutine compute_errors(nnodes, fx_exact, fy_exact, fx_computed, fy_computed, &
                         max_error_x, max_error_y, ndim)
!**********************************************************************
    implicit none
    integer, intent(in) :: nnodes, ndim
    real(8), intent(in) :: fx_exact(ndim), fy_exact(ndim)
    real(8), intent(in) :: fx_computed(ndim), fy_computed(ndim)
    real(8), intent(out) :: max_error_x, max_error_y
    integer :: i
    real(8) :: error
    
    max_error_x = 0.0d0
    max_error_y = 0.0d0
    
    do i = 1, nnodes
        error = abs(fx_exact(i) - fx_computed(i))
        if (error > max_error_x) max_error_x = error
        
        error = abs(fy_exact(i) - fy_computed(i))
        if (error > max_error_y) max_error_y = error
    enddo
end subroutine compute_errors
