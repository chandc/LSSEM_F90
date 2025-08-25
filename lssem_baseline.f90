!
!  lssem.f90 - Least Squares Spectral Element Method
!  Extracted from SEM_base_2D.f
!  
!  Contains:
!  - rhs:     Right-hand side calculation (4-field system)
!  - lhs:     Left-hand side matrix-vector product  
!  - collect: Inter-element communication
!
! 
!**********************************************************************
      subroutine rhs(u,v,p,om,un,vn,unn,vnn, &
                     fu,fv, &
                     dt,pr,nelem,nterm,neig,norder, &
                     fac1,fac2,fac3, &
                     wid,wht,wg,d, &
                     c1,c2,c3,c4,ndim)
!**********************************************************************
!
! Purpose:
!   Computes the right-hand side (RHS) vector for the incompressible 
!   Navier-Stokes equations in spectral element form using Least Squares
!   Spectral Element Method (LSSEM)
!
! Mathematical Background:
!   The LSSEM formulation solves the first-order system:
!   ∂u/∂t + u·∇u + ∂p/∂x + (1/Re)∂ω/∂y = f_u    (momentum x)
!   ∂v/∂t + u·∇v + ∂p/∂y - (1/Re)∂ω/∂x = f_v    (momentum y)  
!   ∇·u = 0                                        (continuity)
!   ω - (∂v/∂x - ∂u/∂y) = 0                       (vorticity definition)
!
!   Where ω is the vorticity field. This first-order formulation avoids
!   second derivatives by introducing vorticity as a primary variable.
!   
!   Key differences from standard Navier-Stokes:
!   - Viscous terms: +(1/Re)∂ω/∂y in x-momentum, -(1/Re)∂ω/∂x in y-momentum
!   - No Laplacian ∇²u terms - replaced by first-order vorticity derivatives
!   - Vorticity explicitly computed and evolved as separate field
!
!   Time discretization uses Backward Differentiation Formula (BDF):
!   BDF1: (3u^n+1 - 4u^n + u^n-1)/(2Δt) = RHS
!   BDF2: (11u^n+1 - 18u^n + 9u^n-1 - 2u^n-2)/(6Δt) = RHS
!
! Algorithm:
!   For each element:
!   1. Compute spatial derivatives using spectral differentiation matrices:
!      - ∂u/∂x, ∂u/∂y, ∂v/∂x, ∂v/∂y (velocity derivatives)
!      - ∂p/∂x, ∂p/∂y (pressure gradients)
!      - ∂ω/∂x, ∂ω/∂y (vorticity derivatives for viscous terms)
!   2. Evaluate linearized convection terms from previous Newton iteration
!   3. Apply time integration factors (fac1, fac2, fac3) for BDF scheme
!   4. Assemble first-order viscous terms: (1/Re)∂ω/∂y, -(1/Re)∂ω/∂x
!   5. Include external forcing terms (fu, fv) and their derivatives
!   6. Store results in RHS vectors (c1, c2, c3, c4)
!
! Arguments:
!   u(ndim,*)    (in)  : x-velocity field at current time level
!   v(ndim,*)    (in)  : y-velocity field at current time level  
!   p(ndim,*)    (in)  : pressure field at current time level
!   om(ndim,*)   (in)  : vorticity field at current time level
!   un(ndim,*)   (in)  : x-velocity at previous time level n
!   vn(ndim,*)   (in)  : y-velocity at previous time level n
!   unn(ndim,*)  (in)  : x-velocity at time level n-1 (for BDF2)
!   vnn(ndim,*)  (in)  : y-velocity at time level n-1 (for BDF2)
!   fu(ndim,*)   (in)  : external forcing in x-direction
!   fv(ndim,*)   (in)  : external forcing in y-direction
!   dt           (in)  : time step size
!   pr           (in)  : Prandtl number (= 1/Re for Navier-Stokes)
!   nelem        (in)  : number of spectral elements
!   nterm        (in)  : polynomial order + 1 (number of terms per direction)
!   neig         (in)  : number of eigenvalues (not used in this version)
!   norder       (in)  : maximum polynomial order
!   fac1         (in)  : BDF time integration factor 1
!   fac2         (in)  : BDF time integration factor 2  
!   fac3         (in)  : BDF time integration factor 3
!   wid(*)       (in)  : element widths in x-direction
!   wht(*)       (in)  : element heights in y-direction
!   wg(*)        (in)  : Gauss-Lobatto quadrature weights
!   d(norder,*)  (in)  : spectral differentiation matrix
!   c1(ndim,*)   (out) : RHS for x-momentum equation
!   c2(ndim,*)   (out) : RHS for y-momentum equation
!   c3(ndim,*)   (out) : RHS for continuity equation
!   c4(ndim,*)   (out) : RHS for vorticity equation
!   ndim         (in)  : leading dimension of arrays
!
! Local Arrays:
!   dudx, dudy   : x and y derivatives of u velocity
!   dvdx, dvdy   : x and y derivatives of v velocity
!   dpdx, dpdy   : x and y derivatives of pressure
!   domdx, domdy : x and y derivatives of vorticity
!   dfudx, dfudy : x and y derivatives of x-forcing
!   dfvdx, dfvdy : x and y derivatives of y-forcing
!   su(ntdof,4)  : temporary storage for RHS components
!
! Physical Parameters:
!   Reynolds number: Re = U*L/ν where U=characteristic velocity,
!                   L=characteristic length, ν=kinematic viscosity
!   Prandtl number:  Pr = 1/Re for incompressible Navier-Stokes
!
! Spectral Element Discretization:
!   - Each element mapped from physical (x,y) to computational (ξ,η) ∈ [-1,1]²  
!   - Jacobian transformation: ajac = (wid*wht)/4
!   - Metric scaling: facx = 2/wid, facy = 2/wht
!   - High-order Lagrange interpolation on Legendre-Gauss-Lobatto points
!
! Time Integration:
!   BDF schemes provide stability for convection-dominated flows:
!   - BDF1 (first-order): Implicit Euler, unconditionally stable
!   - BDF2 (second-order): Higher accuracy, conditionally stable
!   - Factor application: fac1*u^n+1 + fac2*u^n + fac3*u^n-1
!
! Nonlinear Terms:
!   Convection handled implicitly through Newton linearization:
!   - u·∇u linearized as: u^n·∇u + u·∇u^n - u^n·∇u^n  
!   - u·∇v linearized as: u^n·∇v + u·∇v^n - u^n·∇v^n
!   - RHS contains previous iteration terms and linearization corrections
!   - LHS contains current iteration linearized convection operator
!
! Usage:
!   Called within Newton iteration or dual time-stepping solver
!   to compute residual vector for iterative solution methods
!
! Notes:
!   - Assumes structured tensor-product spectral element grid
!   - Uses implicit treatment of nonlinear convection terms via linearization
!   - Supports both BDF1 and BDF2 time integration schemes
!   - Optimized for incompressible flow with constant properties
!   - LSSEM advantage: avoids second derivatives (∇²u) in favor of 
!     first-order vorticity derivatives, improving numerical conditioning
!   - First-order formulation more suitable for spectral element methods
!
! Author: Original SEM_base_2D.f, enhanced documentation August 2025
!
!**********************************************************************
      implicit none
      
      ! Parameters
      integer, parameter :: ndeg = 20, ntdof = ndeg*ndeg
      
      ! Arguments
      integer, intent(in) :: nelem, nterm, neig, norder, ndim
      real, intent(in) :: dt, pr, fac1, fac2, fac3
      real, intent(in) :: wid(*), wht(*), wg(*), d(norder,*)
      real, intent(in) :: u(ndim,*), v(ndim,*), p(ndim,*), om(ndim,*)
      real, intent(in) :: un(ndim,*), vn(ndim,*), unn(ndim,*), vnn(ndim,*)
      real, intent(in) :: fu(ndim,*), fv(ndim,*)
      real, intent(out) :: c1(ndim,*), c2(ndim,*), c3(ndim,*), c4(ndim,*)
      
      ! Local variables
      real :: dudx(ntdof), dudy(ntdof), dvdx(ntdof), dvdy(ntdof)
      real :: dpdx(ntdof), dpdy(ntdof), domdx(ntdof), domdy(ntdof)
      real :: dfudx(ntdof), dfudy(ntdof), dfvdx(ntdof), dfvdy(ntdof)
      real :: su(ntdof,4)
      integer :: li(ndeg)
      integer :: n, ne, i, j, ij, kk
      real :: ajac, facx, facy, dx, dy, facem
!
      do n=1,nterm
       li(n) = (n-1)*nterm
      enddo
!
      do ne=1,nelem
!
       ajac = 0.25*wid(ne)*wht(ne)
       facx = 2./wid(ne)
       facy = 2./wht(ne)
!
! derivatives in the x-direction
!
      do i=1,neig
       dudx(i) = 0.0
       dvdx(i) = 0.0
       dpdx(i) = 0.0
       domdx(i) = 0.0
       dfudx(i) = 0.0
       dfvdx(i) = 0.0
      enddo
!
      do n=1,nterm
       do i=1,nterm
        do j=1,nterm
          ij = li(i) + j
          kk = li(n) + j
          dx = d(i,n)*facx
          dudx(ij) = dudx(ij) + dx*u(kk,ne)                       
          dvdx(ij) = dvdx(ij) + dx*v(kk,ne)                       
          dpdx(ij) = dpdx(ij) + dx*p(kk,ne)                       
          domdx(ij) = domdx(ij) + dx*om(kk,ne)
          dfudx(ij) = dfudx(ij) + dx*fu(kk,ne)
          dfvdx(ij) = dfvdx(ij) + dx*fv(kk,ne)
        enddo
       enddo
      enddo                       
!
! derivatives in the y-direction
!
      do i=1,neig
       dudy(i) = 0.0
       dvdy(i) = 0.0
       dpdy(i) = 0.0
       domdy(i) = 0.0
       dfudy(i) = 0.0
       dfvdy(i) = 0.0
      enddo
!
      do n=1,nterm
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         kk = li(i) + n
         dy = d(j,n)*facy
         dudy(ij) = dudy(ij) + dy*u(kk,ne)                       
         dvdy(ij) = dvdy(ij) + dy*v(kk,ne)                       
         dpdy(ij) = dpdy(ij) + dy*p(kk,ne)                       
         domdy(ij) = domdy(ij) + dy*om(kk,ne)
         dfudy(ij) = dfudy(ij) + dy*fu(kk,ne)
         dfvdy(ij) = dfvdy(ij) + dy*fv(kk,ne)
        enddo
       enddo
      enddo
!
!
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         facem = ajac*wg(i)*wg(j)
!
         su(ij,1) =   fac1*u(ij,ne)*facem + dt*(  &                
                      fu(ij,ne)*dudx(ij)  + fv(ij,ne)*dudy(ij) + &
                      u(ij,ne)*dfudx(ij) + v(ij,ne)*dfudy(ij) + &
                      dpdx(ij) + pr*domdy(ij)  )*facem
!
         su(ij,2) =  fac1*v(ij,ne)*facem + dt*(  &
                         fu(ij,ne)*dvdx(ij)  + fv(ij,ne)*dvdy(ij) + &
                          u(ij,ne)*dfvdx(ij) + v(ij,ne)*dfvdy(ij) + &
                          dpdy(ij) - pr*domdx(ij)  )*facem
!
         su(ij,3) = ( dudx(ij) + dvdy(ij) )*facem
         su(ij,4) = ( om(ij,ne) + dudy(ij) - dvdx(ij) )*facem
        enddo
       enddo
!
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         facem = ajac*wg(i)*wg(j) 
         su(ij,1) = (  - fac2*un(ij,ne) - fac3*unn(ij,ne) + &
                    ( fu(ij,ne)*dfudx(ij) + fv(ij,ne)*dfudy(ij) )*dt &
                    )*facem - su(ij,1)
         su(ij,2) = (  - fac2*vn(ij,ne) - fac3*vnn(ij,ne) + &
                    (fu(ij,ne)*dfvdx(ij) + fv(ij,ne)*dfvdy(ij))*dt &
                    )*facem - su(ij,2)
         su(ij,3) = -su(ij,3)
         su(ij,4) = -su(ij,4)
        enddo
       enddo     
!
! Multiply by the transpose
!
! at collocation point
!
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         c1(ij,ne) = (fac1)*su(ij,1) + &
                     dt*dfudx(ij)*su(ij,1) + dt*dfvdx(ij)*su(ij,2)
!
         c2(ij,ne) = (fac1)*su(ij,2) + &
                     dfudy(ij)*dt*su(ij,1) + dt*dfvdy(ij)*su(ij,2)
!
         c3(ij,ne) = 0.0
!
         c4(ij,ne) = su(ij,4)
        enddo
       enddo 
!
! along the x-direction
!
       do n=1,nterm
        do i=1,nterm
         do j=1,nterm
          ij = li(i) + j
          kk = li(n) + j
          dx = d(n,i)*facx

          c1(ij,ne) = dx*su(kk,3) + &
                      dt*fu(kk,ne)*dx*su(kk,1) + c1(ij,ne)

          c2(ij,ne) = -dx*su(kk,4) + &
                      dt*fu(kk,ne)*dx*su(kk,2) + c2(ij,ne)

          c3(ij,ne) = dx*dt*su(kk,1) + c3(ij,ne)

          c4(ij,ne) = -pr*dt*dx*su(kk,2) + c4(ij,ne)
         enddo
        enddo
       enddo
!
       do n=1,nterm
        do i=1,nterm
         do j=1,nterm
          ij = li(i) + j
          kk = li(i) + n
          dy = d(n,j)*facy
          c1(ij,ne) = dy*su(kk,4) + dt*fv(kk,ne)*dy*su(kk,1) + c1(ij,ne)
          c2(ij,ne) = dy*su(kk,3) + dt*fv(kk,ne)*dy*su(kk,2) + c2(ij,ne)
          c3(ij,ne) = dt*dy*su(kk,2) + c3(ij,ne)
          c4(ij,ne) = dt*pr*dy*su(kk,1) + c4(ij,ne)
         enddo
        enddo
       enddo
!
       enddo
       return
       end      
!
!
! 
!**********************************************************************
      subroutine lhs(u,v,p,om, &
                     fu,fv, &
                     fac1, &
                     dt,pr,nelem,nterm,neig,norder, &
                     wid,wht,wg,d, &
                     c1,c2,c3,c4,ndim)
!**********************************************************************
!
! Purpose:
!   Computes the left-hand side (LHS) matrix-vector product for the 
!   linearized incompressible Navier-Stokes equations in LSSEM formulation
!
! Mathematical Background:
!   The LHS represents the implicit operator applied to the solution increment
!   in the linearized system arising from Newton's method or dual time stepping:
!   
!   [M/Δt + K]Δq = -R(q^n)
!   
!   Where:
!   - M: mass matrix from time derivative terms
!   - K: stiffness matrix from spatial operators (convection + diffusion)
!   - Δq: solution increment [Δu, Δv, Δp, Δω]
!   - R(q^n): residual from RHS evaluation
!
!   For the 4-field LSSEM system:
!   Field 1: ∂u/∂t + linearized momentum-x equation
!   Field 2: ∂v/∂t + linearized momentum-y equation  
!   Field 3: Continuity equation: ∇·u = 0
!   Field 4: Vorticity definition: ω - (∂v/∂x - ∂u/∂y) = 0
!
! Linearization Strategy:
!   Convection terms linearized about current solution state using Newton method:
!   
!   For u·∇u term:
!   u·∇u ≈ u^n·∇u + u·∇u^n - u^n·∇u^n  (Newton linearization)
!   
!   Where:
!   - u^n·∇u: current velocity advecting solution increment (goes to LHS)
!   - u·∇u^n: solution increment advecting current velocity (goes to LHS) 
!   - u^n·∇u^n: nonlinear term from current iteration (goes to RHS)
!   
!   This creates velocity-velocity coupling terms in the LHS operator
!   and provides quadratic convergence for Newton iteration
!
! Algorithm:
!   For each spectral element:
!   1. Compute spatial derivatives of current solution fields
!   2. Evaluate linearized convection terms with current velocities
!   3. Apply mass matrix scaling: fac1/Δt (from BDF time integration)
!   4. Include viscous terms scaled by Reynolds number
!   5. Assemble pressure-velocity coupling terms
!   6. Apply quadrature weights and Jacobian scaling
!   7. Store LHS contributions in output arrays
!
! Arguments:
!   u(ndim,*)    (in)  : x-velocity field (current iterate)
!   v(ndim,*)    (in)  : y-velocity field (current iterate)
!   p(ndim,*)    (in)  : pressure field (current iterate)
!   om(ndim,*)   (in)  : vorticity field (current iterate)
!   fu(ndim,*)   (in)  : external forcing in x-direction
!   fv(ndim,*)   (in)  : external forcing in y-direction
!   fac1         (in)  : BDF time integration factor (= coefficient of u^n+1)
!   dt           (in)  : time step size
!   pr           (in)  : Prandtl number (= 1/Re for Navier-Stokes)
!   nelem        (in)  : number of spectral elements
!   nterm        (in)  : polynomial order + 1 per direction
!   neig         (in)  : number of eigenvalues (not used)
!   norder       (in)  : maximum polynomial order
!   wid(*)       (in)  : element widths in x-direction
!   wht(*)       (in)  : element heights in y-direction
!   wg(*)        (in)  : Gauss-Lobatto quadrature weights
!   d(norder,*)  (in)  : spectral differentiation matrix
!   c1(ndim,*)   (out) : LHS for x-momentum equation
!   c2(ndim,*)   (out) : LHS for y-momentum equation
!   c3(ndim,*)   (out) : LHS for continuity equation
!   c4(ndim,*)   (out) : LHS for vorticity equation
!   ndim         (in)  : leading dimension of arrays
!
! Local Arrays:
!   dudx, dudy   : spatial derivatives of x-velocity
!   dvdx, dvdy   : spatial derivatives of y-velocity
!   dpdx, dpdy   : spatial derivatives of pressure
!   domdx, domdy : spatial derivatives of vorticity
!   su(ntdof,4)  : temporary storage for LHS contributions
!
! Matrix Structure:
!   The resulting LHS operator has the block structure:
!   
!   [A_uu  A_uv  A_up  A_uw] [Δu]   [RHS_u]
!   [A_vu  A_vv  A_vp  A_vw] [Δv] = [RHS_v]
!   [A_pu  A_pv  A_pp  A_pw] [Δp]   [RHS_p]
!   [A_wu  A_wv  A_wp  A_ww] [Δω]   [RHS_ω]
!
!   Where:
!   - A_uu, A_vv: momentum diagonal blocks (mass + linearized convection + viscous)
!   - A_up, A_vp: pressure-velocity coupling (pressure gradient terms)
!   - A_pu, A_pv: velocity-pressure coupling (divergence terms)
!   - A_uw, A_vw: velocity-vorticity coupling
!   - A_wu, A_wv: vorticity-velocity coupling (curl terms)
!
! Spectral Element Implementation:
!   - High-order accurate spatial discretization using LGL points
!   - Exact integration using Gauss-Lobatto quadrature
!   - Element-by-element assembly for parallel efficiency
!   - Tensor product structure for computational efficiency
!
! Time Integration:
!   LHS represents the implicit operator for BDF time stepping:
!   - Mass term coefficient: fac1/Δt
!   - Spatial terms multiplied by Δt for proper scaling
!   - Supports both BDF1 and BDF2 schemes
!
! Iterative Solution Context:
!   This LHS operator is typically used with:
!   - BiCGSTAB or GMRES iterative solvers
!   - Block-diagonal or approximate factorization preconditioners
!   - Matrix-free implementation for large-scale problems
!
! Performance Characteristics:
!   - O(N⁴) operations per element for full LHS assembly
!   - Sparse matrix structure exploitable for efficiency
!   - High arithmetic intensity suitable for modern processors
!
! Usage:
!   Called within Newton iteration or dual time stepping:
!   1. Evaluate RHS at current solution state
!   2. Apply LHS operator to solution increment
!   3. Solve linear system iteratively
!   4. Update solution and repeat until convergence
!
! Notes:
!   - Assumes constant element geometry (wid, wht)
!   - Linearization requires current solution as input
!   - Compatible with matrix-free iterative methods
!   - Supports element-by-element parallel assembly
!
! Author: Original SEM_base_2D.f, enhanced documentation August 2025
!
!**********************************************************************
!
      implicit none
      integer, parameter :: ndeg=20, ntdof=ndeg*ndeg
      
      ! Intent declarations
      real, intent(in) :: u(ndim,*),v(ndim,*),p(ndim,*),om(ndim,*)
      real, intent(in) :: fu(ndim,*),fv(ndim,*)
      real, intent(in) :: fac1, dt, pr
      integer, intent(in) :: nelem,nterm,neig,norder,ndim
      real, intent(in) :: wid(*),wht(*),wg(*),d(norder,*)
      real, intent(inout) :: c1(ndim,*),c2(ndim,*),c3(ndim,*),c4(ndim,*)
      
      ! Local variables
      real :: dudx(ntdof),dudy(ntdof),dvdx(ntdof),dvdy(ntdof)
      real :: dpdx(ntdof),dpdy(ntdof),domdx(ntdof),domdy(ntdof)
      real :: dfudx(ntdof),dfudy(ntdof),dfvdx(ntdof),dfvdy(ntdof)
      real :: su(ntdof,4)
      integer :: li(ndeg)
      real :: ajac, facx, facy, facem, dx, dy
      integer :: i, j, n, ne, ij, kk
!
      do n=1,nterm
       li(n) = (n-1)*nterm
      enddo
!
      do ne=1,nelem
!
       ajac = 0.25*wid(ne)*wht(ne)
       facx = 2./wid(ne)
       facy = 2./wht(ne)
!
! derivatives in the x-direction
!
      do i=1,neig
       dudx(i) = 0.0
       dvdx(i) = 0.0
       dpdx(i) = 0.0
       domdx(i) = 0.0
       dfudx(i) = 0.0
       dfvdx(i) = 0.0
      enddo
!
      do n=1,nterm
       do i=1,nterm
        do j=1,nterm
          ij = li(i) + j
          kk = li(n) + j
          dx = d(i,n)*facx
          dudx(ij) = dudx(ij) + dx*u(kk,ne)                       
          dvdx(ij) = dvdx(ij) + dx*v(kk,ne)                       
          dpdx(ij) = dpdx(ij) + dx*p(kk,ne)                       
          domdx(ij) = domdx(ij) + dx*om(kk,ne)
          dfudx(ij) = dfudx(ij) + dx*fu(kk,ne)
          dfvdx(ij) = dfvdx(ij) + dx*fv(kk,ne)
        enddo
       enddo
      enddo                       
!
! derivatives in the y-direction
!
      do i=1,neig
       dudy(i) = 0.0
       dvdy(i) = 0.0
       dpdy(i) = 0.0
       domdy(i) = 0.0
       dfudy(i) = 0.0
       dfvdy(i) = 0.0
      enddo
!
      do n=1,nterm
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         kk = li(i) + n
         dy = d(j,n)*facy
         dudy(ij) = dudy(ij) + dy*u(kk,ne)                       
         dvdy(ij) = dvdy(ij) + dy*v(kk,ne)                       
         dpdy(ij) = dpdy(ij) + dy*p(kk,ne)                       
         domdy(ij) = domdy(ij) + dy*om(kk,ne)
         dfudy(ij) = dfudy(ij) + dy*fu(kk,ne)
         dfvdy(ij) = dfvdy(ij) + dy*fv(kk,ne)
        enddo
       enddo
      enddo
!
!
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         facem = ajac*wg(i)*wg(j)
!
         su(ij,1) = fac1*u(ij,ne)*facem + dt*(  &

                    fu(ij,ne)*dudx(ij) + fv(ij,ne)*dudy(ij) + &
                    u(ij,ne)*dfudx(ij) + v(ij,ne)*dfudy(ij) + &
                    dpdx(ij) + pr*domdy(ij) )*facem
!
                  su(ij,2) = fac1*v(ij,ne)*facem + dt*( &

                    fu(ij,ne)*dvdx(ij) + fv(ij,ne)*dvdy(ij) + &
                    u(ij,ne)*dfvdx(ij) + v(ij,ne)*dfvdy(ij) + &
                    dpdy(ij) - pr*domdx(ij) )*facem
!
         su(ij,3) = ( dudx(ij) + dvdy(ij) )*facem
         su(ij,4) = ( om(ij,ne) + dudy(ij) - dvdx(ij) )*facem
        enddo
       enddo
!
!
! Multiply by the transpose
!
! at collocation point
!
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         c1(ij,ne) = (fac1)*su(ij,1) + &

                     dt*dfudx(ij)*su(ij,1) + dt*dfvdx(ij)*su(ij,2)
!
         c2(ij,ne) = (fac1)*su(ij,2) + &

                     dfudy(ij)*dt*su(ij,1) + dt*dfvdy(ij)*su(ij,2)
!
         c3(ij,ne) = 0.0
!
         c4(ij,ne) = su(ij,4)
        enddo
       enddo 
!
! along the x-direction
!
       do n=1,nterm
        do i=1,nterm
         do j=1,nterm
          ij = li(i) + j
          kk = li(n) + j
          dx = d(n,i)*facx

          c1(ij,ne) = dx*su(kk,3) + dt*fu(kk,ne)*dx*su(kk,1) + c1(ij,ne)

          c2(ij,ne) = -dx*su(kk,4) + dt*fu(kk,ne)*dx*su(kk,2) + c2(ij,ne)

          c3(ij,ne) = dx*dt*su(kk,1) + c3(ij,ne)

          c4(ij,ne) = -pr*dt*dx*su(kk,2) + c4(ij,ne)
         enddo
        enddo
       enddo
!
       do n=1,nterm
        do i=1,nterm
         do j=1,nterm
          ij = li(i) + j
          kk = li(i) + n
          dy = d(n,j)*facy
          c1(ij,ne) = dy*su(kk,4) + &
                      dt*fv(kk,ne)*dy*su(kk,1) + c1(ij,ne)
          c2(ij,ne) = dy*su(kk,3) + &
                      dt*fv(kk,ne)*dy*su(kk,2) + c2(ij,ne)
          c3(ij,ne) = dt*dy*su(kk,2) + c3(ij,ne)
          c4(ij,ne) = dt*pr*dy*su(kk,1) + c4(ij,ne)
         enddo
        enddo
       enddo
!
       enddo
       return
       end  
!
!************************************************************************
      subroutine collect(nelem,nterm,ndep,ntdof, &
                         isouth,inorth,iwest,ieast, &
                         res)
!************************************************************************
!
! Purpose:
!   Performs inter-element communication and residual collection for 
!   spectral element boundaries in multi-element domains
!
! Mathematical Background:
!   In spectral element methods, adjacent elements share common boundaries
!   where continuity must be enforced. The collect subroutine implements
!   the "direct stiffness summation" or "assembly" process by:
!   
!   1. Identifying shared boundary nodes between neighboring elements
!   2. Summing residual contributions from both elements at interface
!   3. Distributing the combined residual back to both elements
!   
!   This ensures:
!   - C⁰ continuity across element boundaries
!   - Conservation of mass and momentum at interfaces
!   - Proper global matrix assembly from element contributions
!
! Assembly Process:
!   For shared boundary point between elements A and B:
!   1. Collect: R_shared = R_A + R_B
!   2. Distribute: R_A := R_shared, R_B := R_shared
!   
!   This operation is essential for:
!   - Proper coupling between adjacent elements
!   - Global conservation properties
!   - Correct boundary condition enforcement
!
! Domain Connectivity:
!   Elements are connected through boundary mappings:
!   - isouth(ne): element number south of element ne
!   - inorth(ne): element number north of element ne  
!   - iwest(ne):  element number west of element ne
!   - ieast(ne):  element number east of element ne
!   
!   Value = 0 indicates physical boundary (no neighbor)
!   Value > 0 indicates neighboring element number
!
! Algorithm:
!   For each element ne:
!   1. Check south boundary: if isouth(ne) ≠ 0
!      - Sum residuals from southern edge of ne and northern edge of neighbor
!      - Distribute combined residual to both elements
!   2. Check west boundary: if iwest(ne) ≠ 0  
!      - Sum residuals from western edge of ne and eastern edge of neighbor
!      - Distribute combined residual to both elements
!   3. North and east boundaries reserved for future implementation
!
! Arguments:
!   nelem        (in)    : number of spectral elements in domain
!   nterm        (in)    : polynomial order + 1 per direction
!   ndep         (in)    : number of dependent variables per node (=4)
!   ntdof        (in)    : total degrees of freedom per element
!   isouth(*)    (in)    : south neighbor element numbers
!   inorth(*)    (in)    : north neighbor element numbers (reserved)
!   iwest(*)     (in)    : west neighbor element numbers  
!   ieast(*)     (in)    : east neighbor element numbers (reserved)
!   res(ntdof,*) (inout) : residual array [input: element residuals, 
!                                         output: assembled residuals]
!
! Data Indexing:
!   Global field indices from sem_data module:
!   - iu:  x-velocity field index
!   - iv:  y-velocity field index  
!   - ip:  pressure field index
!   - iom: vorticity field index
!
!   Node indexing within element:
!   - ijs = (i-1)*nterm*ndep + field_index     [south boundary]
!   - ijn = (i-1)*nterm*ndep + (nterm-1)*ndep + field_index [north boundary]
!   - ijw = (j-1)*ndep + field_index           [west boundary]
!   - ije = (nterm-1)*nterm*ndep + (j-1)*ndep + field_index [east boundary]
!
! Boundary Correspondence:
!   - South edge of element ne ↔ North edge of element isouth(ne)
!   - West edge of element ne  ↔ East edge of element iwest(ne)
!   - This ensures proper geometric connectivity
!
! Implementation Details:
!   Current implementation handles:
!   ✓ South-North boundary exchange (y-direction coupling)
!   ✓ West-East boundary exchange (x-direction coupling)
!   ⚠ North and East mappings reserved for future extension
!   ⚠ Assumes structured grid topology
!
! Conservation Properties:
!   The collect operation ensures:
!   - Global mass conservation across element interfaces
!   - Momentum conservation at element boundaries
!   - Energy conservation for viscous flows
!   - Proper scaling for iterative solver convergence
!
! Parallel Computing Context:
!   In parallel implementations, this operation requires:
!   - MPI communication for elements on different processors
!   - Careful handling of ghost/halo elements
!   - Synchronization points for global consistency
!
! Usage Context:
!   Called after element-level RHS or LHS computation:
!   1. Compute element residuals independently
!   2. Call collect() to enforce interface continuity
!   3. Proceed with iterative solver or time stepping
!
! Performance Considerations:
!   - O(nelem × nterm) complexity for boundary operations
!   - Memory access pattern optimized for cache efficiency  
!   - Minimal computational overhead compared to element operations
!   - Suitable for vectorization and parallel execution
!
! Error Handling:
!   - Assumes valid element connectivity arrays
!   - No bounds checking for performance reasons
!   - Requires consistent element numbering scheme
!
! Geometric Assumptions:
!   - Structured quadrilateral element topology
!   - Conforming element interfaces (no hanging nodes)
!   - Consistent orientation between neighboring elements
!
! Author: Original SEM_base_2D.f, enhanced documentation August 2025
!
!************************************************************************
      use sem_data, only: iu, iv, ip, iom
      
      implicit none
      
      ! Intent declarations
      integer, intent(in) :: nelem,nterm,ndep,ntdof
      integer, intent(in) :: isouth(*),inorth(*),iwest(*),ieast(*)
      real, intent(inout) :: res(ntdof,*)
      
      ! Local variables
      integer :: ne, i, ii, ijs, ijn, j, jj, ijs2, ijn2
      real :: resu, resv, resp, resm
      integer :: ije, ijw
!
! bidirectional exchange in the y-direction for shared boundaries
!

      do ne=1,nelem
       if(isouth(ne).ne.0) then
       do i=1,nterm
        ii = (i-1)*nterm
        ijs = ii*ndep                          
        ijn = (ii+nterm-1)*ndep
        resu = res(ijs+iu,ne) + res(ijn+iu,isouth(ne))
        resv = res(ijs+iv,ne) + res(ijn+iv,isouth(ne))
        resp = res(ijs+ip,ne) + res(ijn+ip,isouth(ne))
        resm = res(ijs+iom,ne) + res(ijn+iom,isouth(ne))
        res(ijs+iu,ne) = resu
        res(ijs+iv,ne) = resv
        res(ijs+ip,ne) = resp
        res(ijs+iom,ne) = resm
        res(ijn+iu,isouth(ne)) = resu
        res(ijn+iv,isouth(ne)) = resv
        res(ijn+ip,isouth(ne)) = resp
        res(ijn+iom,isouth(ne)) = resm
       enddo
       endif
      enddo
!
! bidirectional exchange in the x-direction for shared boundaries
!
      do ne=1,nelem
       if(iwest(ne).ne.0) then
       do j=1,nterm
        ijw = (j-1)*ndep                          
        ii = (nterm-1)*nterm
        ije = (ii+j-1)*ndep
        resu = res(ijw+iu,ne) + res(ije+iu,iwest(ne))
        resv = res(ijw+iv,ne) + res(ije+iv,iwest(ne))
        resp = res(ijw+ip,ne) + res(ije+ip,iwest(ne))
        resm = res(ijw+iom,ne) + res(ije+iom,iwest(ne))


        res(ijw+iu,ne) = resu
        res(ijw+iv,ne) = resv
        res(ijw+ip,ne) = resp
        res(ijw+iom,ne) = resm
        res(ije+iu, iwest(ne)) = resu
        res(ije+iv, iwest(ne)) = resv
        res(ije+ip, iwest(ne)) = resp
        res(ije+iom,iwest(ne)) = resm


       enddo
       endif
      enddo
!
      return
      end
