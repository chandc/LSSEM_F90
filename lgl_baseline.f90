!
!  lgl_baseline.f90 - Legendre-Gauss-Lobatto Utilities
!  Extracted from SEM_base_2D.f and enhanced with 2D capabilities
!  
!  Contains:
!  - legen:  Legendre polynomial evaluation
!  - quad:   Quadrature weight calculation
!  - derv:   1D differentiation matrix construction  
!  - jacobl: Gauss-Lobatto-Legendre quadrature points
!  - jacobf: Jacobi polynomial evaluation
!  - create_2d_derivative_matrices: 2D derivative matrices via tensor products
!  - create_2d_grid: 2D tensor product grid construction
!
!  Updated: August 2025 - Added 2D spectral element utilities
!
!**********************************************************************
      subroutine legen(al,alp,n,xc,ndim)
!**********************************************************************
!
! Purpose:
!   Evaluates Legendre polynomials and their derivatives at a given point
!   using recurrence relations for efficient computation
!
! Mathematical Background:
!   Legendre polynomials P_k(x) satisfy the recurrence relation:
!   (k+1)*P_{k+1}(x) = (2k+1)*x*P_k(x) - k*P_{k-1}(x)
!   with initial conditions: P_0(x) = 1, P_1(x) = x
!   
!   Their derivatives satisfy:
!   (k+1)*P'_{k+1}(x) = (2k+1)*(P_k(x) + x*P'_k(x)) - k*P'_{k-1}(x)
!   with initial conditions: P'_0(x) = 0, P'_1(x) = 1
!
! Arguments:
!   al(0:ndim) (out) : Legendre polynomial values P_0(xc), ..., P_n(xc)
!   alp(0:ndim)(out) : Legendre polynomial derivatives P'_0(xc), ..., P'_n(xc)
!   n          (in)  : Maximum polynomial degree to evaluate
!   xc         (in)  : Point at which to evaluate polynomials
!   ndim       (in)  : Dimension of arrays al and alp (must be ≥ n+1)
!
! Usage:
!   Used internally by quad() and derv() subroutines for computing
!   quadrature weights and differentiation matrices
!
! Notes:
!   - Recurrence relations provide numerical stability
!   - Avoids factorial computations that would cause overflow
!   - Efficient O(n) algorithm
!
! Author: Original SEM_base_2D.f, enhanced documentation August 2025
!
!**********************************************************************
      implicit none
      integer, intent(in) :: n, ndim
      real(8), intent(in) :: xc
      real(8), intent(out) :: al(0:ndim), alp(0:ndim)
      integer :: k, kp, km
!
      al(0) = 1.
      al(1) = xc
      alp(0) = 0
      alp(1) = 1.
!
      do k=1,n
        kp = k + 1
        km = k - 1
        al(kp) = (2.*k+1.)*xc*al(k)/kp - k*al(km)/kp
      enddo
!
      do k=1,n
        kp = k + 1
        km = k - 1
        alp(kp) = (2.*k+1.)*(al(k)+xc*alp(k))/kp - &
                  k*alp(km)/kp
      enddo
!
      return
      end
!
!**********************************************************************
      subroutine quad(n,x,w,ndim)
!**********************************************************************
!
! Purpose:
!   Computes Legendre-Gauss-Lobatto quadrature weights for given
!   collocation points
!
! Mathematical Background:
!   For Legendre-Gauss-Lobatto quadrature, the weights are given by:
!   w_i = 2 / (n*(n+1)*[P_n(x_i)]²)
!   where P_n(x_i) is the n-th Legendre polynomial evaluated at x_i
!
!   These weights satisfy the property:
!   ∫₋₁¹ f(x) dx ≈ Σᵢ₌₀ⁿ w_i * f(x_i)
!   which is exact for polynomials of degree ≤ 2n-1
!
! Arguments:
!   n         (in)  : Polynomial degree (number of intervals)
!   x(0:ndim) (in)  : LGL collocation points from jacobl() subroutine
!   w(0:ndim) (out) : Corresponding quadrature weights
!   ndim      (in)  : Dimension of arrays x and w (must be ≥ n+1)
!
! Algorithm:
!   1. For each collocation point x_i
!   2. Evaluate Legendre polynomial P_n(x_i) using legen()
!   3. Compute weight using formula w_i = 2/(n*(n+1)*P_n(x_i)²)
!   4. Add small value to prevent division by zero
!
! Usage:
!   Typically called after jacobl() to get complete quadrature rule:
!   call jacobl(n, 0.0d0, 0.0d0, x, ndim)
!   call quad(n, x, w, ndim)
!
! Notes:
!   - Weights are positive and sum to 2 (length of interval [-1,1])
!   - End point weights are typically smaller than interior weights
!   - Essential for spectral element numerical integration
!
! Author: Original SEM_base_2D.f, enhanced documentation August 2025
!
!**********************************************************************
      implicit none
      integer, parameter :: nn = 100
      integer, intent(in) :: n, ndim
      real(8), intent(in) :: x(0:ndim)
      real(8), intent(out) :: w(0:ndim)
      real(8) :: alp1(0:nn), al1(0:nn)
      real(8) :: small, xc
      integer :: k
!
!  determine the Gauss Quadrature weighting factors
!
      small = 1.0e-30
      do k=0,n
       xc = x(k)
       call legen(al1,alp1,n,xc,nn)  
       w(k) = 2. / &
              ( n*(n+1)*al1(n)*al1(n) + small )
      enddo
      return
      end
!
!**********************************************************************
      subroutine derv(nterm,x,d,ndim)
!**********************************************************************
!
! Purpose:
!   Constructs the spectral differentiation matrix for Legendre-Gauss-Lobatto
!   collocation points using Lagrange polynomial derivatives
!
! Mathematical Background:
!   For spectral collocation methods, the derivative of a function f at
!   collocation point x_i is approximated as:
!   f'(x_i) ≈ Σⱼ₌₀ⁿ D_ij * f(x_j)
!   
!   The differentiation matrix elements are:
!   D_ij = P_n(x_i) / (P_n(x_j) * (x_i - x_j))  for i ≠ j
!   D_ii = 0  for interior points
!   D_00 = -n*(n+1)/4  (left boundary)
!   D_nn = +n*(n+1)/4  (right boundary)
!
! Arguments:
!   nterm       (in)  : Polynomial degree (same as n in other subroutines)
!   x(0:ndim)   (in)  : LGL collocation points from jacobl()
!   d(0:ndim,0:ndim)(out): Differentiation matrix D_ij
!   ndim        (in)  : Dimension of arrays (must be ≥ nterm+1)
!
! Algorithm:
!   1. For each pair of collocation points (x_i, x_j)
!   2. If i ≠ j: Compute D_ij using Legendre polynomial ratio
!   3. If i = j: Set D_ii = 0 for interior points
!   4. Apply special boundary conditions for endpoints
!
! Usage:
!   Essential for spectral element methods:
!   call jacobl(n, 0.0d0, 0.0d0, x, ndim)
!   call derv(n, x, d, ndim)
!   ! Then: df_dx = matmul(d, f) gives derivatives
!
! Notes:
!   - Matrix has exact sum-to-zero property for interior rows
!   - Achieves spectral accuracy for smooth functions
!   - Forms the foundation for 2D derivative matrices via tensor products
!   - Boundary point derivatives have special treatment
!
! Author: Original SEM_base_2D.f, enhanced documentation August 2025
!
!**********************************************************************
      implicit none
      integer, parameter :: nn = 100
      integer, intent(in) :: nterm, ndim
      real(8), intent(in) :: x(0:ndim)
      real(8), intent(out) :: d(0:ndim,0:ndim)
      real(8) :: al1(0:nn), alp1(0:nn), al2(0:nn), alp2(0:nn)
      real(8) :: xi, xj, ann
      integer :: i, j
!
!  determine the derivative at the collocation points
!
      do i=0,nterm
        xi = x(i)
        call legen(al1,alp1,nterm,xi,nn)  
       do j=0,nterm
        xj = x(j)
        call legen(al2,alp2,nterm,xj,nn)  
        if(i == j) then
         d(i,j) = 0
        else
         d(i,j) = al1(nterm)/(al2(nterm)*(xi-xj))
        endif
       enddo
      enddo
!
      ann = 0.25*nterm*(nterm+1)
      d(0,0) = -ann
      d(nterm,nterm) =  ann
      return
      end
!
!**********************************************************************
      subroutine jacobl(n,alpha,beta,xcol,ndim)
!**********************************************************************
!
! Purpose:
!   Computes Gauss-Lobatto collocation points for Jacobi polynomials
!   using Newton-Raphson iteration with high precision
!
! Mathematical Background:
!   Gauss-Lobatto points are the roots of (1-x²)P'_n(x) = 0 where
!   P_n(x) is the n-th Jacobi polynomial with parameters α and β.
!   
!   Special cases:
!   - α = β = 0:    Legendre-Gauss-Lobatto points (most common in SEM)
!   - α = β = -0.5: Chebyshev-Gauss-Lobatto points
!
!   Properties:
!   - Always includes boundary points x = ±1
!   - Points cluster near boundaries for better resolution
!   - Symmetric about x = 0 for α = β
!   - Even n: includes point at x = 0
!   - Odd n: no point exactly at x = 0
!
! Arguments:
!   n            (in)  : Polynomial degree (number of intervals)
!   alpha        (in)  : Jacobi polynomial parameter α
!   beta         (in)  : Jacobi polynomial parameter β  
!   xcol(0:ndim) (out) : Computed collocation points, ordered from -1 to +1
!   ndim         (in)  : Dimension of xcol array (must be ≥ n+1)
!
! Algorithm:
!   1. Set boundary points: x_0 = -1, x_n = +1
!   2. Use trigonometric initial guess for interior points
!   3. Apply Newton-Raphson iteration to find roots
!   4. Use symmetry to reduce computational cost
!   5. Reorder points from -1 to +1
!
! Convergence:
!   - Tolerance: eps = 1.0e-12 (near machine precision)
!   - Maximum iterations: kstop = 10 (typically converges in 2-3 iterations)
!   - Uses deflation to avoid previously found roots
!
! Usage:
!   For Legendre-Gauss-Lobatto points (most common):
!   call jacobl(n, 0.0d0, 0.0d0, xcol, ndim)
!   
!   For Chebyshev-Gauss-Lobatto points:
!   call jacobl(n, -0.5d0, -0.5d0, xcol, ndim)
!
! Notes:
!   - Uses sem_data module for global variables alp, bet, rv
!   - Foundation for all spectral element computations
!   - Points provide optimal accuracy for polynomial interpolation
!   - Achieves exponential convergence for smooth functions
!
! Author: Original SEM_base_2D.f, enhanced documentation August 2025
!
!**********************************************************************
      use sem_data
      implicit none
      
      ! Arguments
      integer, intent(in) :: n, ndim
      real(8), intent(in) :: alpha, beta
      real(8), intent(out) :: xcol(0:ndim)
      
      ! Local variables
      real(8) :: xjac(200)
      integer, parameter :: kstop = 10
      real(8), parameter :: eps = 1.0e-12
      integer :: np, nh, npp, j, k, kk, i, jm
      real(8) :: pnp1p, pdnp1p, pnp, pdnp, pnm1p, pdnm1
      real(8) :: pnp1m, pdnp1m, pnm, pdnm, pnm1m
      real(8) :: det, rp, rm, a, b
      real(8) :: dth, cd, sd, cs, ss, cssave
      real(8) :: x, pnp1, pdnp1, pn, pdn, pnm1
      real(8) :: poly, pder, recsum, delx
!
      alp = alpha
      bet = beta
      rv = 1. + alp
      np = n + 1
!
      call jacobf(np,pnp1p,pdnp1p,pnp,pdnp,pnm1p,pdnm1,1.0)
      call jacobf(np,pnp1m,pdnp1m,pnm,pdnm,pnm1m,pdnm1,-1.0)
      det = pnp*pnm1m - pnm*pnm1p
      rp = -pnp1p
      rm = -pnp1m
      a = (rp*pnm1m - rm*pnm1p)/det
      b = (rm*pnp   - rp*pnm)/det
      xjac(1) = 1.0
      nh = (n+1)/2
      dth = 3.14159265/(2*n+1)
      cd = cos(2.*dth)
      sd = sin(2.*dth)
      cs = cos(dth)
      ss = sin(dth)
!
      do j=2,nh
       x = cs
       do k=1,kstop
        call jacobf(np,pnp1,pdnp1,pn,pdn,pnm1,pdnm1,x)
        poly = pnp1 + a*pn + b*pnm1
        pder = pdnp1 + a*pdn + b*pdnm1
        recsum = 0.0
        jm = j - 1
        do i=1,jm
          recsum = recsum + 1.0/(x-xjac(i))
        enddo
        delx = -poly/(pder-recsum*poly)
        x = x +delx
        if(abs(delx) < eps) exit
       enddo
        xjac(j) = x
        cssave = cs*cd - ss*sd
        ss = cs*sd + ss*cd
        cs = cssave
      enddo
        xjac(np) = -1.0
        npp = n + 2
        do i=2,nh
          xjac(npp-i) = -xjac(i)
        enddo
        if(n /= 2*(n/2)) then 
        do k=0,n
         kk = n - k + 1
         xcol(k) = xjac(kk)
         enddo
        else
        xjac(nh+1) = 0.0
        do k=0,n
         kk = n - k + 1
         xcol(k) = xjac(kk)
         enddo
        endif
!
      return
        end
!
!**********************************************************************
        subroutine jacobf(n,poly,pder,polym1,pderm1,polym2,pderm2,x)
!**********************************************************************
!
! Purpose:
!   Evaluates Jacobi polynomial P_n^(α,β)(x) and its derivative at a given point,
!   along with the previous two polynomials P_{n-1} and P_{n-2}
!
! Mathematical Background:
!   Jacobi polynomials are orthogonal polynomials defined on [-1,1] with
!   weight function w(x) = (1-x)^α (1+x)^β for α,β > -1.
!   
!   Three-term recurrence relation:
!   2n(n+α+β)(2n+α+β-2) P_n = [(2n+α+β-1)(α²-β²) + (2n+α+β-2)(2n+α+β-1)(2n+α+β)x] P_{n-1}
!                              - 2(n+α-1)(n+β-1)(2n+α+β) P_{n-2}
!   
!   Special cases:
!   - α = β = 0:    Legendre polynomials (most common in SEM)
!   - α = β = -0.5: Chebyshev polynomials of first kind
!   - α = β = 0.5:  Chebyshev polynomials of second kind
!
!   Starting values:
!   - P_0(x) = 1
!   - P_1(x) = (α+β+2)/2 * x + (α-β)/2
!
! Arguments:
!   n        (in)  : Degree of polynomial to evaluate
!   x        (in)  : Point at which to evaluate polynomial [-1 ≤ x ≤ 1]
!   poly     (out) : Value of P_n^(α,β)(x)
!   pder     (out) : Value of derivative P'_n^(α,β)(x)
!   polym1   (out) : Value of P_{n-1}^(α,β)(x)
!   pderm1   (out) : Value of P'_{n-1}^(α,β)(x)
!   polym2   (out) : Value of P_{n-2}^(α,β)(x)
!   pderm2   (out) : Value of P'_{n-2}^(α,β)(x)
!
! Algorithm:
!   1. Initialize P_0(x) = 1, P'_0(x) = 0
!   2. Compute P_1(x) = rv*x where rv = (α+β+2)/2
!   3. Apply three-term recurrence relation for k = 2,3,...,n
!   4. Compute derivatives using recurrence for derivatives
!   5. Track previous two polynomials and their derivatives
!
! Global Variables Used:
!   alp : α parameter from sem_data module
!   bet : β parameter from sem_data module  
!   rv  : (α+β+2)/2 scaling factor from sem_data module
!
! Usage:
!   For Legendre polynomial evaluation (α = β = 0):
!   call jacobf(5, poly, pder, polym1, pderm1, polym2, pderm2, 0.5d0)
!   
!   For Chebyshev polynomial evaluation (α = β = -0.5):
!   call jacobf(3, poly, pder, polym1, pderm1, polym2, pderm2, 0.0d0)
!
! Notes:
!   - Essential for Newton-Raphson iteration in root-finding algorithms
!   - Used extensively in Gauss-Lobatto point computation
!   - Provides both function values and derivatives in single call
!   - Maintains numerical stability through careful coefficient computation
!   - Returns previous polynomials needed for deflation in root-finding
!
! Numerical Properties:
!   - Stable evaluation using forward recurrence
!   - O(n) computational complexity
!   - Maintains orthogonality properties to machine precision
!
! Author: Original SEM_base_2D.f, enhanced documentation August 2025
!
!**********************************************************************
        use sem_data
        implicit none
        
        ! Arguments  
        integer, intent(in) :: n
        real(8), intent(in) :: x
        real(8), intent(out) :: poly, pder, polym1, pderm1, polym2, pderm2
        
        ! Local variables
        real(8) :: apb, psave, pdsave, polylst, pderlst
        real(8) :: a1, a2, a3, a4, b3, polyn, pdern
        integer :: k
        
        apb = alp + bet
        poly = 1.0
        pder = 0.0
        psave = 0.0
        pdsave = 0.0
        if(n == 0) return
        polylst = poly
        pderlst = pder
        poly = rv*x
        pder = rv
        if(n == 1) return
        do k=2,n
          a1 = 2.*k*(k+apb)*(2.*k+apb-2.)
          a2 = (2.*k+apb-1.)*(alp**2-bet**2)
          b3 = (2.*k+apb-2.)
          a3 = b3*(b3+1.)*(b3+2.)
          a4 = 2.*(k+alp-1.)*(k+bet-1.)*(2.*k+apb)
          polyn = ((a2+a3*x)*poly - a4*polylst) / a1
          pdern = ((a2+a3*x)*pder - a4*pderlst + a3*poly) / a1
          psave = polylst
          pdsave = pderlst
          polylst = poly
          poly = polyn
          pderlst = pder
          pder = pdern
        enddo
        polym1 = polylst
        pderm1 = pderlst
        polym2 = psave
        pderm2 = pdsave
        return
        end

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
! Tensor Product Implementation:
!   For node (i,j) → global index k = i + j*(N+1):
!   - X-derivative: affects nodes with same j, different i
!   - Y-derivative: affects nodes with same i, different j
!   
!   Matrix elements:
!   - dx_2d(k,l): contribution from node l to x-derivative at node k
!   - dy_2d(k,l): contribution from node l to y-derivative at node k
!
! Arguments:
!   n         (in)  : Polynomial degree (order = n+1)
!   d_1d      (in)  : 1D differentiation matrix [(n+1) × (n+1)]
!   dx_2d     (out) : 2D x-derivative matrix [(n+1)² × (n+1)²]
!   dy_2d     (out) : 2D y-derivative matrix [(n+1)² × (n+1)²]
!   ndim_1d   (in)  : Leading dimension of d_1d array (≥ n+1)
!   ndim_2d   (in)  : Leading dimension of 2D matrices (≥ (n+1)²)
!
! Algorithm:
!   1. Initialize matrices to zero
!   2. For each target node (i,j):
!      - X-derivative: copy d_1d row i to appropriate positions
!      - Y-derivative: copy d_1d row j to appropriate positions
!   3. Apply tensor product structure systematically
!
! Storage Requirements:
!   - Input:  1D matrix: (n+1)² elements
!   - Output: 2D matrices: 2 × (n+1)⁴ elements total
!   - Memory scales as O(n⁴) for 2D problems
!
! Usage:
!   ! First compute 1D derivative matrix
!   call derv(n, d_1d, xcol, ndim_1d)
!   
!   ! Then create 2D matrices
!   call create_2d_derivative_matrices(n, d_1d, dx_2d, dy_2d, ndim_1d, ndim_2d)
!   
!   ! Apply to solution vector
!   dudx_2d = matmul(dx_2d, u_2d)  ! ∂u/∂x
!   dudy_2d = matmul(dy_2d, u_2d)  ! ∂u/∂y
!
! Applications:
!   - Navier-Stokes equations: velocity gradients, vorticity
!   - Heat equation: temperature gradients, heat flux
!   - Wave equations: spatial derivatives for time stepping
!   - Poisson equations: Laplacian operator construction
!
! Performance Notes:
!   - Matrices are sparse with (n+1)² non-zeros per row
!   - Consider sparse storage for large n
!   - Assembly cost: O(n⁴), usage cost: O(n⁴) per application
!
! Accuracy:
!   - Inherits spectral accuracy from 1D differentiation
!   - Exponential convergence for smooth functions
!   - Machine precision for polynomial functions up to degree n
!
! Author: Enhanced SEM library August 2025
!
!**********************************************************************
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
!   - This subroutine complements the 1D utilities: legen, quad, derv, jacobl
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
!
!     Initialize matrices to zero
!
      dx_2d = 0.0d0
      dy_2d = 0.0d0
!
!     Create x-derivative matrix (differentiate in x, identity in y)
!     For each y-line (j constant), apply 1D derivative matrix in x-direction
!
      do j = 0, n          ! y-direction index (constant for each row)
        do i = 0, n        ! x-direction index (row in 1D derivative)
          row_idx = i + 1 + j * (n + 1)
          do k = 0, n      ! x-direction basis function (column in 1D derivative)
            col_idx = k + 1 + j * (n + 1)
            dx_2d(row_idx, col_idx) = d_1d(i, k)
          enddo
        enddo
      enddo
!
!     Create y-derivative matrix (identity in x, differentiate in y)
!     For each x-line (i constant), apply 1D derivative matrix in y-direction
!
      do j = 0, n          ! y-direction index (row in 1D derivative)
        do i = 0, n        ! x-direction index (constant for each row)
          row_idx = i + 1 + j * (n + 1)
          do l = 0, n      ! y-direction basis function (column in 1D derivative)
            col_idx = i + 1 + l * (n + 1)
            dy_2d(row_idx, col_idx) = d_1d(j, l)
          enddo
        enddo
      enddo
!
      return
      end

!**********************************************************************
      subroutine create_2d_grid(n, x_1d, x_2d, y_2d, ndim_1d, ndim_2d)
!**********************************************************************
!
! Purpose:
!   Creates 2D tensor product grid from 1D Legendre-Gauss-Lobatto collocation 
!   points for spectral element methods
!
! Mathematical Background:
!   For spectral element methods, 2D computational domains use tensor product
!   grids constructed from 1D collocation points. Given 1D LGL points 
!   ξ_i ∈ [-1,1] for i = 0,1,...,N, the 2D grid consists of all combinations:
!   
!   Grid points: (ξ_i, ξ_j) for i,j ∈ {0,1,...,N}
!   Total points: (N+1)² nodes per element
!   
!   Properties:
!   - Clustered near boundaries for better resolution
!   - Optimal for high-order polynomial interpolation
!   - Compatible with Gauss-Lobatto quadrature
!   - Enables spectral accuracy for smooth functions
!
! Grid Ordering Convention:
!   Linear indexing: node_index = i + j*(N+1) + 1 (Fortran 1-based)
!   This creates row-major ordering:
!   - Row 0: (0,0), (1,0), (2,0), ..., (N,0)  [y = ξ_0 = -1]
!   - Row 1: (0,1), (1,1), (2,1), ..., (N,1)  [y = ξ_1]
!   - ...
!   - Row N: (0,N), (1,N), (2,N), ..., (N,N)  [y = ξ_N = +1]
!
! Arguments:
!   n        (in)  : Polynomial degree (number of intervals)
!   x_1d     (in)  : 1D LGL collocation points [ξ_0, ξ_1, ..., ξ_N]
!   x_2d     (out) : x-coordinates of 2D grid points [(N+1)² total]
!   y_2d     (out) : y-coordinates of 2D grid points [(N+1)² total]
!   ndim_1d  (in)  : Leading dimension of x_1d array (≥ N+1)
!   ndim_2d  (in)  : Leading dimension of 2D arrays (≥ (N+1)²)
!
! Algorithm:
!   1. Initialize node counter to 1 (Fortran indexing)
!   2. Outer loop: j = 0 to N (y-direction, rows)
!   3. Inner loop: i = 0 to N (x-direction, columns)  
!   4. Assign: x_2d(node) = ξ_i, y_2d(node) = ξ_j
!   5. Increment node counter
!
! Grid Visualization (N=2 example):
!   Node layout with (x,y) coordinates:
!   7:(ξ₀,ξ₂) 8:(ξ₁,ξ₂) 9:(ξ₂,ξ₂)    [y = +1]
!   4:(ξ₀,ξ₁) 5:(ξ₁,ξ₁) 6:(ξ₂,ξ₁)    [y = ξ₁]  
!   1:(ξ₀,ξ₀) 2:(ξ₁,ξ₀) 3:(ξ₂,ξ₀)    [y = -1]
!   [x = -1]  [x = ξ₁]  [x = +1]
!
! Usage:
!   ! Generate 1D LGL points
!   call jacobl(n, 0.0d0, 0.0d0, x_1d, ndim_1d)
!   
!   ! Create 2D tensor product grid
!   call create_2d_grid(n, x_1d, x_2d, y_2d, ndim_1d, ndim_2d)
!   
!   ! Use for function evaluation
!   do k = 1, (n+1)**2
!     u_2d(k) = func(x_2d(k), y_2d(k))
!   enddo
!
! Applications:
!   - Spectral element mesh generation
!   - Function interpolation and evaluation
!   - Numerical integration setup
!   - Boundary condition enforcement
!   - Visualization and post-processing
!
! Compatibility:
!   - Works with create_2d_derivative_matrices() subroutine
!   - Compatible with standard SEM data structures
!   - Ordering matches matrix-vector operations
!   - Suitable for parallel decomposition
!
! Performance:
!   - O((N+1)²) complexity (optimal)
!   - Minimal memory overhead
!   - Cache-friendly access pattern
!
! Author: Enhanced SEM library August 2025
!
!**********************************************************************
      implicit none
      integer, intent(in) :: n, ndim_1d, ndim_2d
      real(8), intent(in) :: x_1d(0:ndim_1d)
      real(8), intent(out) :: x_2d(ndim_2d), y_2d(ndim_2d)
      integer :: i, j, node_idx
!
!     Create 2D tensor product grid with proper ordering
!
      node_idx = 1
      do j = 0, n         ! y-direction (outer loop)
        do i = 0, n       ! x-direction (inner loop)
          x_2d(node_idx) = x_1d(i)
          y_2d(node_idx) = x_1d(j)
          node_idx = node_idx + 1
        enddo
      enddo
!
      return
      end
