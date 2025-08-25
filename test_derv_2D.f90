!
! test_derv_2D.f90 - Test 2D derivative matrices using tensor products
!
! This program tests the construction of 2D derivative matrices from 1D
! differentiation matrices and validates spectral accuracy
!
program test_derv_2D
    implicit none
    
    ! Parameters
    integer, parameter :: max_n = 16
    integer, parameter :: max_dim = max_n + 1
    integer, parameter :: max_nodes_2d = max_dim * max_dim
    
    ! Variables
    integer :: n, i, j, k
    integer :: nnodes_1d, nnodes_2d
    real(8) :: x_1d(0:max_dim), w_1d(0:max_dim)
    real(8) :: d_1d(0:max_dim, 0:max_dim)
    real(8) :: dx_2d(max_nodes_2d, max_nodes_2d)
    real(8) :: dy_2d(max_nodes_2d, max_nodes_2d)
    real(8) :: x_2d(max_nodes_2d), y_2d(max_nodes_2d)
    real(8) :: f_vals(max_nodes_2d), fx_exact(max_nodes_2d), fy_exact(max_nodes_2d)
    real(8) :: fx_computed(max_nodes_2d), fy_computed(max_nodes_2d)
    real(8) :: error_x, error_y, max_error_x, max_error_y
    
    write(*,*) '================================'
    write(*,*) 'Testing 2D Derivative Matrices'
    write(*,*) '================================'
    write(*,*)
    
    ! Test for different polynomial orders
    write(*,'(A8,A15,A15)') 'Order N', 'Max Error dx', 'Max Error dy'
    write(*,'(A8,A15,A15)') '-------', '------------', '------------'
    
    do n = 2, 12, 2  ! Test even orders from 2 to 12
        nnodes_1d = n + 1
        nnodes_2d = nnodes_1d * nnodes_1d
        
        ! Step 1: Get 1D collocation points and derivative matrix
        call setup_1d_grid(n, x_1d, w_1d, d_1d, max_dim)
        
        ! Step 2: Create 2D grid and derivative matrices
        call create_2d_grid(n, x_1d, x_2d, y_2d, max_dim, max_nodes_2d)
        call create_2d_derivative_matrices(n, d_1d, dx_2d, dy_2d, max_dim, max_nodes_2d)
        
        ! Step 3: Test with known function f(x,y) = sin(π*x)*cos(π*y)
        call test_function_sincos(nnodes_2d, x_2d, y_2d, f_vals, fx_exact, fy_exact, max_nodes_2d)
        
        ! Step 4: Compute derivatives using 2D matrices
        call apply_derivative_matrix(nnodes_2d, dx_2d, f_vals, fx_computed, max_nodes_2d)
        call apply_derivative_matrix(nnodes_2d, dy_2d, f_vals, fy_computed, max_nodes_2d)
        
        ! Step 5: Compute errors
        call compute_errors(nnodes_2d, fx_exact, fy_exact, fx_computed, fy_computed, &
                           max_error_x, max_error_y, max_nodes_2d)
        
        write(*,'(I8,2E15.6)') n, max_error_x, max_error_y
    enddo
    
    write(*,*)
    write(*,*) 'Test completed successfully!'
    
end program test_derv_2D

!**********************************************************************
subroutine setup_1d_grid(n, x, w, d, ndim)
!**********************************************************************
    implicit none
    integer, intent(in) :: n, ndim
    real(8), intent(out) :: x(0:ndim), w(0:ndim), d(0:ndim, 0:ndim)
    
    ! Get Legendre-Gauss-Lobatto points
    call jacobl(n, 0.0d0, 0.0d0, x, ndim)
    
    ! Get quadrature weights
    call quad(n, x, w, ndim)
    
    ! Get 1D derivative matrix
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
    
    ! Create 2D tensor product grid
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
    
    ! Initialize matrices
    dx_2d = 0.0d0
    dy_2d = 0.0d0
    
    ! Create x-derivative matrix (differentiate in x, identity in y)
    do j = 0, n          ! y-direction index
        do i = 0, n      ! x-direction index
            row_idx = i + 1 + j * (n + 1)
            do k = 0, n  ! x-direction basis function
                col_idx = k + 1 + j * (n + 1)
                dx_2d(row_idx, col_idx) = d_1d(i, k)
            enddo
        enddo
    enddo
    
    ! Create y-derivative matrix (identity in x, differentiate in y)
    do j = 0, n          ! y-direction index
        do i = 0, n      ! x-direction index
            row_idx = i + 1 + j * (n + 1)
            do l = 0, n  ! y-direction basis function
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
    
    ! Test function: f(x,y) = sin(π*x) * cos(π*y)
    ! df/dx = π * cos(π*x) * cos(π*y)
    ! df/dy = -π * sin(π*x) * sin(π*y)
    
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
    
    ! Matrix-vector multiplication: df = D * f
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
