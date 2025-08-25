!
! compare_even_odd_orders.f90 - Compare accuracy between even and odd polynomial orders
!
program compare_even_odd_orders
    implicit none
    
    ! Parameters
    integer, parameter :: max_n = 22
    integer, parameter :: max_dim = max_n + 1
    integer, parameter :: max_nodes_2d = max_dim * max_dim
    
    ! Variables
    integer :: n
    integer :: nnodes_1d, nnodes_2d
    real(8) :: x_1d(0:max_dim), w_1d(0:max_dim)
    real(8) :: d_1d(0:max_dim, 0:max_dim)
    real(8) :: dx_2d(max_nodes_2d, max_nodes_2d)
    real(8) :: dy_2d(max_nodes_2d, max_nodes_2d)
    real(8) :: x_2d(max_nodes_2d), y_2d(max_nodes_2d)
    real(8) :: f_vals(max_nodes_2d), fx_exact(max_nodes_2d), fy_exact(max_nodes_2d)
    real(8) :: fx_computed(max_nodes_2d), fy_computed(max_nodes_2d)
    real(8) :: max_error_x, max_error_y
    logical :: has_zero_point
    integer :: i
    character(5) :: parity_str, zero_str
    
    write(*,*) '====================================================='
    write(*,*) 'Comparison: Even vs Odd Polynomial Orders'
    write(*,*) 'Test Function: sin(πx)cos(πy)'
    write(*,*) '====================================================='
    write(*,*)
    write(*,'(A8,A8,A15,A15,A12,A10)') 'Order N', 'Parity', 'Max Error dx', 'Max Error dy', 'DOF (N+1)²', 'x=0 Point'
    write(*,'(A8,A8,A15,A15,A12,A10)') '-------', '------', '------------', '------------', '-----------', '---------'
    
    ! Test both even and odd orders from 2 to 21
    do n = 2, 21
        nnodes_1d = n + 1
        nnodes_2d = nnodes_1d * nnodes_1d
        
        ! Setup 1D grid and derivative matrix
        call setup_1d_grid(n, x_1d, w_1d, d_1d, max_dim)
        
        ! Check if there's a point at x=0
        has_zero_point = .false.
        do i = 0, n
            if (abs(x_1d(i)) < 1.0d-12) then
                has_zero_point = .true.
                exit
            endif
        enddo
        
        ! Create 2D grid and derivative matrices
        call create_2d_grid(n, x_1d, x_2d, y_2d, max_dim, max_nodes_2d)
        call create_2d_derivative_matrices(n, d_1d, dx_2d, dy_2d, max_dim, max_nodes_2d)
        
        ! Test trigonometric function
        call test_function_sincos(nnodes_2d, x_2d, y_2d, f_vals, fx_exact, fy_exact, max_nodes_2d)
        
        ! Compute derivatives
        call apply_derivative_matrix(nnodes_2d, dx_2d, f_vals, fx_computed, max_nodes_2d)
        call apply_derivative_matrix(nnodes_2d, dy_2d, f_vals, fy_computed, max_nodes_2d)
        
        ! Compute errors
        call compute_errors(nnodes_2d, fx_exact, fy_exact, fx_computed, fy_computed, &
                           max_error_x, max_error_y, max_nodes_2d)
        
        ! Set output strings
        if (mod(n, 2) == 0) then
            parity_str = 'Even'
        else
            parity_str = 'Odd '
        endif
        
        if (has_zero_point) then
            zero_str = 'Yes'
        else
            zero_str = 'No '
        endif
        
        write(*,'(I8,A8,2E15.6,I12,A10)') n, parity_str, max_error_x, max_error_y, nnodes_2d, zero_str
    enddo
    
    write(*,*)
    write(*,*) 'Key Observations:'
    write(*,*) '1. Even orders (N=2,4,6,...) include a collocation point at x=0'
    write(*,*) '2. Odd orders (N=3,5,7,...) do NOT include x=0 as a collocation point'
    write(*,*) '3. Both achieve spectral convergence rates'
    write(*,*) '4. Even orders may have slight advantages for symmetric functions'
    write(*,*)
    
    ! Compare specific pairs
    write(*,*) 'Direct Comparison - Adjacent Orders:'
    write(*,'(A12,A12,A15,A15)') 'Even/Odd', 'DOF Ratio', 'Error Ratio dx', 'Error Ratio dy'
    write(*,'(A12,A12,A15,A15)') '--------', '---------', '--------------', '--------------'
    
    call compare_adjacent_orders()
    
end program compare_even_odd_orders

!**********************************************************************
subroutine compare_adjacent_orders()
!**********************************************************************
    implicit none
    integer, parameter :: max_n = 22, max_dim = max_n + 1, max_nodes_2d = max_dim * max_dim
    integer :: n_even, n_odd, dof_even, dof_odd
    real(8) :: x_1d(0:max_dim), w_1d(0:max_dim), d_1d(0:max_dim, 0:max_dim)
    real(8) :: dx_2d(max_nodes_2d, max_nodes_2d), dy_2d(max_nodes_2d, max_nodes_2d)
    real(8) :: x_2d(max_nodes_2d), y_2d(max_nodes_2d)
    real(8) :: f_vals(max_nodes_2d), fx_exact(max_nodes_2d), fy_exact(max_nodes_2d)
    real(8) :: fx_computed(max_nodes_2d), fy_computed(max_nodes_2d)
    real(8) :: error_x_even, error_y_even, error_x_odd, error_y_odd
    real(8) :: dof_ratio, error_ratio_x, error_ratio_y
    character(8) :: pair_str
    
    ! Compare pairs: (4,5), (6,7), (8,9), etc.
    do n_even = 4, 16, 2
        n_odd = n_even + 1
        
        ! Test even order
        call setup_1d_grid(n_even, x_1d, w_1d, d_1d, max_dim)
        call create_2d_grid(n_even, x_1d, x_2d, y_2d, max_dim, max_nodes_2d)
        call create_2d_derivative_matrices(n_even, d_1d, dx_2d, dy_2d, max_dim, max_nodes_2d)
        dof_even = (n_even + 1) * (n_even + 1)
        call test_function_sincos(dof_even, x_2d, y_2d, f_vals, fx_exact, fy_exact, max_nodes_2d)
        call apply_derivative_matrix(dof_even, dx_2d, f_vals, fx_computed, max_nodes_2d)
        call apply_derivative_matrix(dof_even, dy_2d, f_vals, fy_computed, max_nodes_2d)
        call compute_errors(dof_even, fx_exact, fy_exact, fx_computed, fy_computed, &
                           error_x_even, error_y_even, max_nodes_2d)
        
        ! Test odd order
        call setup_1d_grid(n_odd, x_1d, w_1d, d_1d, max_dim)
        call create_2d_grid(n_odd, x_1d, x_2d, y_2d, max_dim, max_nodes_2d)
        call create_2d_derivative_matrices(n_odd, d_1d, dx_2d, dy_2d, max_dim, max_nodes_2d)
        dof_odd = (n_odd + 1) * (n_odd + 1)
        call test_function_sincos(dof_odd, x_2d, y_2d, f_vals, fx_exact, fy_exact, max_nodes_2d)
        call apply_derivative_matrix(dof_odd, dx_2d, f_vals, fx_computed, max_nodes_2d)
        call apply_derivative_matrix(dof_odd, dy_2d, f_vals, fy_computed, max_nodes_2d)
        call compute_errors(dof_odd, fx_exact, fy_exact, fx_computed, fy_computed, &
                           error_x_odd, error_y_odd, max_nodes_2d)
        
        ! Compute ratios
        dof_ratio = real(dof_odd) / real(dof_even)
        error_ratio_x = error_x_even / error_x_odd
        error_ratio_y = error_y_even / error_y_odd
        
        write(pair_str, '(I1,A1,I2)') n_even, '/', n_odd
        write(*,'(A12,F12.2,2E15.6)') pair_str, dof_ratio, error_ratio_x, error_ratio_y
    enddo
    
end subroutine compare_adjacent_orders

! Include common subroutines (same as before)
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
