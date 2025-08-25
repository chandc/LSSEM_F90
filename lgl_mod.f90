! lgl_mod.f90 - Legendre-Gauss-Lobatto Utilities (1-based indexing)
module lgl_mod
    use sem_data
    implicit none
contains
    subroutine legen(al, alp, n, xc)
        !**********************************************************************
        implicit none
        integer, intent(in) :: n
        real(rprec), intent(in) :: xc
        real(rprec), intent(out) :: al(n+1), alp(n+1)
        integer :: k
        !
        al(1) = 1._rprec
        alp(1) = 0._rprec
        if (n > 0) then
            al(2) = xc
            alp(2) = 1._rprec
        endif
        !
        do k=2,n
            al(k+1) = ( (2*k-1)*xc*al(k) - (k-1)*al(k-1) ) / k
            alp(k+1) = ( (2*k-1)*(al(k) + xc*alp(k)) - (k-1)*alp(k-1) ) / k
        enddo
        !
    end subroutine legen

    subroutine quad(n, x, w)
        !**********************************************************************
        implicit none
        integer, intent(in) :: n
        real(rprec), intent(in) :: x(n+1)
        real(rprec), intent(out) :: w(n+1)
        real(rprec) :: al(n+1), alp(n+1)
        real(rprec) :: small
        integer :: k
        !
        small = 1.0e-30_rprec
        do k=1,n+1
            call legen(al, alp, n, x(k))
            w(k) = 2._rprec / ( n*(n+1)*al(n+1)*al(n+1) + small )
        enddo
    end subroutine quad

    subroutine derv(n, x, d)
        !**********************************************************************
        implicit none
        integer, intent(in) :: n
        real(rprec), intent(in) :: x(n+1)
        real(rprec), intent(out) :: d(n+1, n+1)
        real(rprec) :: al_i(n+1), alp_i(n+1)
        real(rprec) :: al_j(n+1), alp_j(n+1)
        real(rprec) :: ann
        integer :: i, j
        !
        do i=1,n+1
            call legen(al_i, alp_i, n, x(i))
            do j=1,n+1
                if (i == j) then
                    d(i,j) = 0._rprec
                else
                    call legen(al_j, alp_j, n, x(j))
                    d(i,j) = al_i(n+1) / (al_j(n+1) * (x(i)-x(j)))
                endif
            enddo
        enddo
        !
        ann = 0.25_rprec*n*(n+1)
        d(1,1) = -ann
        d(n+1,n+1) = ann
    end subroutine derv

    subroutine jacobl(n, alpha, beta, xcol)
        !**********************************************************************
        implicit none
        integer, intent(in) :: n
        real(rprec), intent(in) :: alpha, beta
        real(rprec), intent(out) :: xcol(n+1)
        !
        real(rprec) :: xjac(n+1)
        integer, parameter :: kstop = 20
        real(rprec), parameter :: eps = 1.0e-12_rprec
        integer :: np, nh, j, k, i
        real(rprec) :: pnp1p, pdnp1p, pnp, pdnp, pnm1p, pdnm1
        real(rprec) :: pnp1m, pdnp1m, pnm, pdnm, pnm1m
        real(rprec) :: det, rp, rm, a, b
        real(rprec) :: dth, cd, sd, cs, ss, cssave
        real(rprec) :: x, pnp1, pdnp1, pn, pdn, pnm1
        real(rprec) :: poly, pder, recsum, delx
        !
        alp = alpha
        bet = beta
        np = n + 1
        !
        call jacobf(np, pnp1p, pdnp1p, pnp, pdnp, pnm1p, pdnm1, 1.0d0)
        call jacobf(np, pnp1m, pdnp1m, pnm, pdnm, pnm1m, pdnm1, -1.0d0)
        !
        det = pnp*pnm1m - pnm*pnm1p
        rp = -pnp1p
        rm = -pnp1m
        a = (rp*pnm1m - rm*pnm1p)/det
        b = (rm*pnp - rp*pnm)/det
        !
        xjac(1) = 1.0_rprec
        nh = (n+1)/2
        dth = acos(-1.0_rprec)/(2.0_rprec*n)
        !
        do j=2,nh+1
            x = cos((2*j-3)*dth)
            do k=1,kstop
                call jacobf(n, pnp1, pdnp1, pn, pdn, pnm1, pdnm1, x)
                poly = pnp1 + a*pn + b*pnm1
                pder = pdnp1 + a*pdn + b*pdnm1
                recsum = 0.0_rprec
                do i=1,j-1
                    recsum = recsum + 1.0_rprec/(x-xjac(i))
                enddo
                delx = -poly/(pder-recsum*poly)
                x = x + delx
                if(abs(delx) < eps) exit
            enddo
            xjac(j) = x
        enddo
        !
        xjac(n+1) = -1.0_rprec
        do i=2,nh+1
            xjac(n+2-i) = -xjac(i)
        enddo
        if (mod(n,2) /= 0) xjac(nh+1) = 0.0_rprec
        !
        do i=1,n+1
            xcol(i) = xjac(n+2-i)
        enddo
        !
    end subroutine jacobl

    subroutine jacobf(n, poly, pder, polym1, pderm1, polym2, pderm2, x)
        !**********************************************************************
        implicit none
        integer, intent(in) :: n
        real(rprec), intent(in) :: x
        real(rprec), intent(out) :: poly, pder, polym1, pderm1, polym2, pderm2
        !
        real(rprec) :: apb, psave, pdsave, polylst, pderlst
        real(rprec) :: a1, a2, a3, a4, b3, polyn, pdern
        integer :: k
        !
        apb = alp + bet
        poly = 1.0_rprec
        pder = 0.0_rprec
        psave = 0.0_rprec
        pdsave = 0.0_rprec
        if (n == 0) return
        !
        polylst = poly
        pderlst = pder
        poly = (alp - bet + (apb + 2)*x) / 2.0_rprec
        pder = (apb + 2) / 2.0_rprec
        if (n == 1) return
        !
        do k=2,n
            a1 = 2*k*(k+apb)*(2*k+apb-2)
            a2 = (2*k+apb-1)*(alp**2-bet**2)
            b3 = (2*k+apb-2)
            a3 = b3*(b3+1)*(b3+2)
            a4 = 2*(k+alp-1)*(k+bet-1)*(2*k+apb)
            polyn = ((a2+a3*x)*poly - a4*polylst) / a1
            pdern = ((a2+a3*x)*pder - a4*pderlst + a3*poly) / a1
            psave = polylst
            pdsave = pderlst
            polylst = poly
            poly = polyn
            pderlst = pder
            pder = pdern
        enddo
        !
        polym1 = polylst
        pderm1 = pderlst
        polym2 = psave
        pderm2 = pdsave
    end subroutine jacobf
end module lgl_mod
