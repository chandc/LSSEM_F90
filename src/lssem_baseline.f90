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
         su(ij,1) =   fac1*u(ij,ne)*facem + dt*( u(ij,ne) + &                
                      fu(ij,ne)*dudx(ij)  + fv(ij,ne)*dudy(ij) + &
                      u(ij,ne)*dfudx(ij) + v(ij,ne)*dfudy(ij) + &
                      dpdx(ij) + pr*domdy(ij)  )*facem
!
         su(ij,2) =  fac1*v(ij,ne)*facem + dt*( v(ij,ne) + &
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
         su(ij,1) = ( dt*fu(ij,ne) - fac2*un(ij,ne) - fac3*unn(ij,ne) + &
                    ( fu(ij,ne)*dfudx(ij) + fv(ij,ne)*dfudy(ij) )*dt &
                    )*facem - su(ij,1)
         su(ij,2) = ( dt*fv(ij,ne) - fac2*vn(ij,ne) - fac3*vnn(ij,ne) + &
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
         c1(ij,ne) = (dt+fac1)*su(ij,1) + &
                     dt*dfudx(ij)*su(ij,1) + dt*dfvdx(ij)*su(ij,2)
!
         c2(ij,ne) = (dt+fac1)*su(ij,2) + &
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
         su(ij,1) = fac1*u(ij,ne)*facem + dt*( u(ij,ne) + &
                    fu(ij,ne)*dudx(ij) + fv(ij,ne)*dudy(ij) + &
                    u(ij,ne)*dfudx(ij) + v(ij,ne)*dfudy(ij) + &
                    dpdx(ij) + pr*domdy(ij) )*facem
!
         su(ij,2) = fac1*v(ij,ne)*facem + dt*( v(ij,ne) + &
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
         c1(ij,ne) = (dt+fac1)*su(ij,1) + &
                     dt*dfudx(ij)*su(ij,1) + dt*dfvdx(ij)*su(ij,2)
!
         c2(ij,ne) = (dt+fac1)*su(ij,2) + &
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
! Note: inorth and ieast arguments are reserved for future use
! Currently only isouth and iwest boundary exchanges are implemented
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
