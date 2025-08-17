!
!**********************************************************
      subroutine dge(nelem,nterm,ndep,ntdof,norder, &
                     fac1,dt,pr, &
                     diag,wid,wht,wg, &
                     f,d)
!**********************************************************
      use sem_data, only: iu, iv, ip, iom
      implicit none
! --- Arguments ---
      integer :: nelem, nterm, ndep, ntdof, norder
      real :: fac1, dt, pr
      real :: diag(ntdof,*), wid(*), wht(*), wg(*)
      real :: f(ntdof,*), d(norder,*)
! --- Local Variables ---
      real :: aa(4,4)
      integer :: neig, nee, i, j, ne, ii, ij, k1, k2, kk1, lu, lv
      integer :: m, n, mm, mn, ip1
      real :: facx, facy, ajac, facem, uo, vo, dudx, dudy, dvdx, dvdy
      real :: dx, dy, shapx, shapy, shape
!
       neig = nterm*nterm
       nee = neig*ndep
!
      do i=1,ndep
       do j=1,ndep
        aa(i,j) = 0.0
       enddo
      enddo
!
      do 1000 ne=1,nelem
!
      do i=1,nee
       diag(i,ne) = 0.0
      enddo
!
      facx = 2./wid(ne)
      facy = 2./wht(ne)
      ajac = 0.25*wid(ne)*wht(ne)
!
      do i=1,nterm
       ii = (i-1)*nterm
       do j=1,nterm    
       ij = ii + j
!
      do k1=1,nterm  
      do k2=1,nterm
!
      facem = ajac*wg(k1)*wg(k2)
!
       kk1 = (k1-1)*nterm + k2
       kk1 = (kk1-1)*ndep
       lu = kk1 + iu
       lv = kk1 + iv
       uo = f(lu,ne)
       vo = f(lv,ne)
!
      dudx = 0.0
      dudy = 0.0
      dvdx = 0.0
      dvdy = 0.0
!
      do m=1,nterm
       mm = (m-1)*nterm
       do n=1,nterm    
       mn = (mm+n-1)*ndep
       lu = mn + iu
       lv = mn + iv
       shapx = 0.0
       shapy = 0.0
       if(k1.eq.m) shapx = 1.0
       if(k2.eq.n) shapy = 1.0
       dx = d(k1,m)*facx*shapy         
       dy = d(k2,n)*facy*shapx 
       dudx = dudx + dx*f(lu,ne)
       dudy = dudy + dy*f(lu,ne)
       dvdx = dvdx + dx*f(lv,ne)
       dvdy = dvdy + dy*f(lv,ne)
       enddo
      enddo
!
          shapx = 0.0
          shapy = 0.0
          if(k1.eq.i) shapx = 1.0
          if(k2.eq.j) shapy = 1.0
          dx = d(k1,i)*facx*shapy         
          dy = d(k2,j)*facy*shapx 
          shape = shapx*shapy
!
          aa(1,1) =  fac1*shape + (shape &
                    + uo*dx + vo*dy + dudx*shape)*dt
          aa(1,2) =  shape*dt*dudy
          aa(1,3) =  dx*dt
          aa(1,4) =  pr*dy*dt

          aa(2,1) =  shape*dvdx*dt
          aa(2,2) =  fac1*shape + (shape &
                    + uo*dx + vo*dy + dvdy*shape)*dt
          aa(2,3) =  dy*dt  
          aa(2,4) = -pr*dx*dt

          aa(3,1) = dx
          aa(3,2) = dy

          aa(4,1) =   dy
          aa(4,2) =  -dx
          aa(4,4) =  shape  

!
        ip1 = (ij-1)*ndep
        diag(ip1+1,ne) = diag(ip1+1,ne) + ( &
                          aa(1,1)*aa(1,1) + aa(2,1)*aa(2,1) + &
                          aa(3,1)*aa(3,1) + aa(4,1)*aa(4,1) )*facem
        diag(ip1+2,ne) = diag(ip1+2,ne) + ( &
                          aa(1,2)*aa(1,2) + aa(2,2)*aa(2,2) + &
                          aa(3,2)*aa(3,2) + aa(4,2)*aa(4,2) )*facem
        diag(ip1+3,ne) = diag(ip1+3,ne) + ( &
                          aa(1,3)*aa(1,3) + aa(2,3)*aa(2,3) + &
                          aa(3,3)*aa(3,3) + aa(4,3)*aa(4,3) )*facem
        diag(ip1+4,ne) = diag(ip1+4,ne) + ( &
                          aa(1,4)*aa(1,4) + aa(2,4)*aa(2,4) + &
                          aa(3,4)*aa(3,4) + aa(4,4)*aa(4,4) )*facem

!
      enddo
      enddo
      enddo
      enddo
!
1000  continue
!
      return
      end

!
!**********************************************************************
      subroutine precon(nelem, nee, ntdof, nem, &
                        r_in, r_out, diag, mask, small)
!**********************************************************************
! Apply the diagonal preconditioner M^-1 * r_in = r_out
! M is the diagonal matrix 'diag'.
!**********************************************************************
      implicit none
! --- Arguments ---
      integer :: ntdof, nem, nelem, nee
      double precision :: small
      real :: r_in(ntdof,nem), r_out(ntdof,nem)
      real :: diag(ntdof,nem)
      integer :: mask(ntdof,nem)
! --- Local Variables ---
      integer :: i, ne
!234567
       do ne=1,nelem
        do i=1,nee
          r_out(i,ne) = r_in(i,ne) / (diag(i,ne) + small)* mask(i,ne)
        enddo
      enddo

      return
      end

!
!**********************************************************************
      subroutine bicgstab(f, res, diag, mask, &
                          p, v_b, s, t, r_hat, phat, shat, &
                          nelem, nterm, ndep, ntdof, nem, norder, ndim, &
                          pr, dt, fac1, nitmax, tol, iprt, &
                          wid, wht, wg, d, small, &
                          u, v_vel, pp, om, fu, fv, &
                          u_res, v_res, p_res, om_res, &
                          isouth, inorth, iwest, ieast)
!**********************************************************************
! Solves A*x = b using Preconditioned BiCGSTAB (F77 Compatible).
!   f    = Solution vector (Input: guess, Output: result)
!   res  = Residual vector (Input: initial, Output: final)
!   diag = Preconditioner diagonal
!   mask = Boundary condition mask
!   p, v_b, s, t, r_hat, phat, shat = Workspace arrays
!   v_vel= Physical velocity 'v' for lhs
!   ...  = Other parameters & work arrays
!**********************************************************************
      implicit none

! --- Arguments ---
      integer :: ntdof, nem, nelem, nterm, ndep, norder, ndim
      integer :: nitmax, iprt
      double precision :: pr, dt, fac1, tol, small
      real :: f(ntdof,nem), res(ntdof,nem), diag(ntdof,nem)
      integer :: mask(ntdof,nem)
      real :: p(ntdof,nem), v_b(ntdof,nem), s(ntdof,nem)
      real :: t(ntdof,nem), r_hat(ntdof,nem)
      real :: phat(ntdof,nem), shat(ntdof,nem)
      real :: wid(*), wht(*), wg(*), d(norder,*)
      real :: u(ndim,*), v_vel(ndim,*), pp(ndim,*), om(ndim,*)
      real :: fu(ndim,*), fv(ndim,*)
      real :: u_res(ndim,*), v_res(ndim,*), p_res(ndim,*)
      real :: om_res(ndim,*)
      integer :: isouth(nem), inorth(nem), iwest(nem), ieast(nem)

! --- Local Variables ---
      integer :: nee, neig, i, ne, ij, kk, iter
      double precision :: rho, rho_old, alpha, beta, omega
      double precision :: res_norm, res0, dot_prod

! --- External Subroutines ---
      external lhs, precon

      neig = nterm * nterm
      nee = neig * ndep

! --- BiCGSTAB Initialization ---

! r_0 = b - A*x0 (This is 'res' on entry)
! Choose r_hat_0 = r_0
      do ne = 1, nelem
          do i = 1, nee
              r_hat(i,ne) = res(i,ne) * mask(i,ne)
              p(i,ne) = 0.0
              v_b(i,ne) = 0.0
          enddo
      enddo

      rho = 1.0
      alpha = 1.0
      omega = 1.0

! Calculate initial residual norm
      res0 = 0.0
      do ne = 1, nelem
          do i = 1, nee
              res0 = res0 + res(i,ne) * res(i,ne) * mask(i,ne)
          enddo
      enddo
      res0 = sqrt(res0)
      res_norm = res0

      if (iprt .eq. 1) print 105, 0, res_norm

! --- BiCGSTAB Iteration Loop ---

      do iter = 1, nitmax

          rho_old = rho

! rho_i = (r_hat_0)^T * r_{i-1}
          dot_prod = 0.0
          do ne = 1, nelem
              do i = 1, nee
                  dot_prod = dot_prod + r_hat(i,ne) * res(i,ne) &
                                      * mask(i,ne)
              enddo
          enddo
          rho = dot_prod

          if (abs(rho) .lt. small) then
              print *, 'BiCGSTAB Breakdown: rho = 0, iter = ', iter
              goto 999
          endif

          beta = (rho / rho_old) * (alpha / omega)

! p_i = r_{i-1} + beta * (p_{i-1} - omega * v_b_{i-1})
          do ne = 1, nelem
              do i = 1, nee
                  p(i,ne) = (res(i,ne) + beta * (p(i,ne) - &
                            omega * v_b(i,ne))) * mask(i,ne)
              enddo
          enddo

! Solve M * phat = p_i (Preconditioning Step 1)
! --- FIX 1a: Broke call precon into two lines ---
          call precon(nelem, nee, ntdof, nem, p, phat, diag, &
                      mask, small)

! v_b_i = A * phat (Matrix-Vector Product 1)
          do ne=1,nelem
             do ij=1,neig
                kk = (ij-1)*ndep
                u(ij,ne) = phat(kk+1,ne)
                v_vel(ij,ne) = phat(kk+2,ne)
                pp(ij,ne) = phat(kk+3,ne)
                om(ij,ne) = phat(kk+4,ne)
             enddo
          enddo
          call lhs(u,v_vel,pp,om,fu,fv, fac1, &
                   dt,pr,nelem,nterm,neig,norder, &
                   wid,wht,wg,d, &
                   u_res,v_res,p_res,om_res,ndim)
          do ne=1,nelem
             do i=1,neig
                ij=(i-1)*ndep
                v_b(ij+1,ne) = u_res(i,ne) * mask(ij+1,ne)
                v_b(ij+2,ne) = v_res(i,ne) * mask(ij+2,ne)
                v_b(ij+3,ne) = p_res(i,ne) * mask(ij+3,ne)
                v_b(ij+4,ne) = om_res(i,ne)* mask(ij+4,ne)
             enddo
          enddo

! --- ADD THIS CALL ---
          call collect(nelem,nterm,ndep,ntdof, &
                       isouth,inorth,iwest,ieast, &
                       v_b)
! --- END ADDITION ---          

! alpha = rho_i / (r_hat_0^T * v_b_i)
          dot_prod = 0.0
          do ne = 1, nelem
              do i = 1, nee
                  dot_prod = dot_prod + r_hat(i,ne) * v_b(i,ne) &
                                      * mask(i,ne)
              enddo
          enddo

          if (abs(dot_prod) .lt. small) then
              print *, 'BiCGSTAB Breakdown: r_hat_0_T * v = 0, iter=', &
                       iter
              goto 999
          endif
          alpha = rho / dot_prod

! s = r_{i-1} - alpha * v_b_i
          do ne = 1, nelem
              do i = 1, nee
                  s(i,ne) = (res(i,ne) - alpha * v_b(i,ne)) &
                           * mask(i,ne)
              enddo
          enddo

! Check convergence for 's' (early exit)
          res_norm = 0.0
          do ne = 1, nelem
              do i = 1, nee
                  res_norm = res_norm + s(i,ne) * s(i,ne) * mask(i,ne)
              enddo
          enddo
          res_norm = sqrt(res_norm)

          if (res_norm .lt. tol) then
              do ne = 1, nelem
                  do i = 1, nee
                      f(i,ne) = f(i,ne) + alpha * phat(i,ne) &
                               * mask(i,ne)
                  enddo
              enddo
              do ne = 1, nelem
                  do i = 1, nee
                      res(i,ne) = s(i,ne)
                  enddo
              enddo
              if (iprt .eq. 1) print 105, iter, res_norm
              print *, 'BiCGSTAB Converged (early)!'
              goto 999
          endif

! Solve M * shat = s (Preconditioning Step 2)
! --- FIX 1b: Broke call precon into two lines ---
          call precon(nelem, nee, ntdof, nem, s, shat, diag, &
                      mask, small)

! t = A * shat (Matrix-Vector Product 2)
          do ne=1,nelem
             do ij=1,neig
                kk = (ij-1)*ndep
                u(ij,ne) = shat(kk+1,ne)
                v_vel(ij,ne) = shat(kk+2,ne)
                pp(ij,ne) = shat(kk+3,ne)
                om(ij,ne) = shat(kk+4,ne)
             enddo
          enddo
          call lhs(u,v_vel,pp,om,fu,fv, fac1, &
                   dt,pr,nelem,nterm,neig,norder, &
                   wid,wht,wg,d, &
                   u_res,v_res,p_res,om_res,ndim)
          do ne=1,nelem
             do i=1,neig
                ij=(i-1)*ndep
                t(ij+1,ne) = u_res(i,ne) * mask(ij+1,ne)
                t(ij+2,ne) = v_res(i,ne) * mask(ij+2,ne)
                t(ij+3,ne) = p_res(i,ne) * mask(ij+3,ne)
                t(ij+4,ne) = om_res(i,ne)* mask(ij+4,ne)
             enddo
          enddo

! --- ADD THIS CALL ---
          call collect(nelem,nterm,ndep,ntdof, &
                       isouth,inorth,iwest,ieast, &
                       t)
! --- END ADDITION ---


! omega = (t^T * s) / (t^T * t)
          dot_prod = 0.0
          do ne = 1, nelem
              do i = 1, nee
                  dot_prod = dot_prod + t(i,ne) * t(i,ne) * mask(i,ne)
              enddo
          enddo

          if (abs(dot_prod) .lt. small) then
              print *, 'BiCGSTAB Breakdown: t_T * t = 0, iter=', iter
              goto 999
          endif

          omega = 0.0
          do ne = 1, nelem
              do i = 1, nee
                  omega = omega + t(i,ne) * s(i,ne) * mask(i,ne)
              enddo
          enddo
          omega = omega / dot_prod

          if (abs(omega) .lt. small) then
! --- FIX 2: Changed print statement to use FORMAT ---
              print 106, iter
              goto 999
          endif

! x_i = x_{i-1} + alpha * phat + omega * shat
          do ne = 1, nelem
              do i = 1, nee
! changed by DCC, May 26, 2025
!
                  f(i,ne) = f(i,ne) + (alpha * phat(i,ne) &
                           + omega * shat(i,ne)) &
                           * mask(i,ne)
              enddo
          enddo

! r_i = s - omega * t
          do ne = 1, nelem
              do i = 1, nee
                  res(i,ne) = (s(i,ne) - omega * t(i,ne)) * mask(i,ne)
              enddo
          enddo

! Check convergence
          res_norm = 0.0
          do ne = 1, nelem
              do i = 1, nee
                  res_norm = res_norm + res(i,ne) * res(i,ne) &
                           * mask(i,ne)
              enddo
          enddo
          res_norm = sqrt(res_norm)

          if (iprt .eq. 1) print 105, iter, res_norm

          if (res_norm .lt. tol) then
              print *, 'BiCGSTAB Converged!'
              goto 999
          endif

      enddo

      print *, 'BiCGSTAB failed to converge within ', nitmax, ' iters'
      print *, 'Final Residual Norm = ', res_norm

 999  continue

 105  format('BiCGSTAB iter = ',i5,5x,'residual= ',e14.6)
 106  format('BiCGSTAB Stagnation/Breakdown: omega=0, iter=', I5)

      return
      end