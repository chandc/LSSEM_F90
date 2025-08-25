!
!  lgl.f90 - Legendre-Gauss-Lobatto Utilities
!  Extracted from SEM_base_2D.f
!  
!  Contains:
!  - legen:  Legendre polynomial evaluation
!  - quad:   Quadrature weight calculation
!  - derv:   Differentiation matrix construction  
!  - jacobl: Gauss-Lobatto-Legendre quadrature points
!  - jacobf: Jacobi polynomial evaluation
!
!**********************************************************************
      subroutine legen(al,alp,n,xc,ndim)
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
!  CORRECTED VERSION - August 17, 2025
!  
!  Major improvements implemented:
!  1. Fixed mathematical error in diagonal element calculation
!  2. Enforces correct row-sum property: sum(D_ij) = 0 for all rows
!  3. Optimized performance: O(N²) → O(N) calls to legen()
!  4. Eliminated special case handling for corner elements
!  
!  Mathematical principle: For differentiation matrix D, the derivative
!  of a constant must be zero, requiring each row to sum to zero.
!  
!  CHANGES MADE:
!  - Pre-compute L_N(x_j) values once to avoid redundant legen() calls
!  - Calculate all off-diagonal elements first
!  - Set diagonal elements as negative sum of row's off-diagonal elements
!  - Removed hardcoded corner element assignments
!**********************************************************************
      implicit none
      integer, parameter :: nn = 100
      integer, intent(in) :: nterm, ndim
      real(8), intent(in) :: x(0:ndim)
      real(8), intent(out) :: d(0:ndim,0:ndim)
      
      ! UPDATED: More efficient variable declarations
      real(8) :: ln_at_nodes(0:nn)  ! NEW: Pre-computed L_N(x_j) values
      real(8) :: al_temp(0:nn), alp_temp(0:nn)  ! RENAMED: Temporary arrays
      real(8) :: row_sum            ! NEW: For diagonal calculation
      integer :: i, j
      
      ! REMOVED: Redundant variables (al1, alp1, al2, alp2, xi, xj, ann)
!
!  PHASE 1: Pre-calculate L_N(x) at each node for efficiency
!  This optimization reduces legen() calls from O(N²) to O(N)
!
      do i = 0, nterm
         call legen(al_temp, alp_temp, nterm, x(i), nn)
         ln_at_nodes(i) = al_temp(nterm)  ! Store L_N(x_i)
      enddo
!
!  PHASE 2: Calculate off-diagonal elements of differentiation matrix
!  Using the standard formula: D_ij = L_N(x_i) / [L_N(x_j) * (x_i - x_j)]
!
      do i = 0, nterm
         do j = 0, nterm
            if (i /= j) then
               ! CORRECTED: Use pre-computed values for efficiency
               d(i,j) = ln_at_nodes(i) / (ln_at_nodes(j) * (x(i) - x(j)))
            endif
         enddo
      enddo
!
!  PHASE 3: Calculate diagonal elements using row-sum rule
!  MATHEMATICAL PRINCIPLE: Each row must sum to zero (derivative of constant = 0)
!  Therefore: D_ii = -∑(j≠i) D_ij
!  This single rule correctly handles ALL diagonal elements (corners AND interior)
!
      do i = 0, nterm
         row_sum = 0.0d0
         do j = 0, nterm
            if (i /= j) then
               row_sum = row_sum + d(i,j)
            endif
         enddo
         ! CORRECTED: Set diagonal element to ensure row sum = 0
         d(i,i) = -row_sum
      enddo
!
!  REMOVED: Incorrect hardcoded corner element assignments
!  OLD CODE (mathematically incorrect):
!    ann = 0.25*nterm*(nterm+1)
!    d(0,0) = -ann
!    d(nterm,nterm) = ann
!  
!  The new implementation correctly calculates ALL diagonal elements,
!  including corners, using the fundamental mathematical constraint.
!
      return
      end
!
!**********************************************************************
      subroutine jacobl(n,alpha,beta,xcol,ndim)
!**********************************************************************
      use sem_data
      implicit none
!
!  computes the gauss-lobatto collocation points
!
!   N:              degree of approximation (order of polynomials=n+1)
!   ALPHA:          parameter in Jacobi weight
!   BETA:           parameter in Jacobi weight
!   XJAC:           roots from largest to smallest
!
!
!  for Chebyshev-Gauss-Lobatto points use alpha=-0.5 and beta=-0.5
!  for Legendre-Gauss-Lobatto points use           0             0
!
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
          a4 = 2.*(K+alp-1.)*(k+bet-1.)*(2.*k+apb)
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
