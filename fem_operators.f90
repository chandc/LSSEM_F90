
! ==============================================================================
! Module: fem_operators
!
! Description:
!   This module provides subroutines for applying fundamental finite element
!   matrix operations (Mass, Stiffness/Gradient) to data fields within a
!   spectral element method framework. These operators act on a field on an
!   element-by-element basis.
!
! Contains:
!   - apply_mass_matrix: Applies the mass matrix (weighted identity).
!   - apply_gradient: Computes the gradient of a scalar field.
!   - apply_divergence: Computes the divergence of a vector field.
!
! ==============================================================================
MODULE fem_operators

  IMPLICIT NONE

CONTAINS

  ! ============================================================================
  ! apply_mass_matrix
  !
  ! Applies the elemental mass matrix to a field u.
  ! Operation: res = M * u
  ! where M is the diagonal mass matrix with entries (ajac * wg(i) * wg(j)).
  ! ============================================================================
  SUBROUTINE apply_mass_matrix(res, u, nelem, nterm, neig, &
                               wid, wht, wg, ndim)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nelem, nterm, neig, ndim
    REAL, INTENT(IN) :: u(ndim, nelem), wid(nelem), wht(nelem), wg(nterm)
    REAL, INTENT(OUT) :: res(ndim, nelem)

    INTEGER :: ne, i, j, ij
    REAL :: ajac, facem
    INTEGER, DIMENSION(nterm) :: li

    DO i = 1, nterm
      li(i) = (i - 1) * nterm
    END DO

    DO ne = 1, nelem
      ajac = 0.25 * wid(ne) * wht(ne)
      DO i = 1, nterm
        DO j = 1, nterm
          ij = li(i) + j
          facem = ajac * wg(i) * wg(j)
          res(ij, ne) = u(ij, ne) * facem
        END DO
      END DO
    END DO
  END SUBROUTINE apply_mass_matrix

  ! ============================================================================
  ! apply_gradient
  !
  ! Computes the gradient of a scalar field p.
  ! Operation: (res_x, res_y) = G * p
  ! where G is the discrete gradient operator.
  ! ============================================================================
  SUBROUTINE apply_gradient(res_x, res_y, p, nelem, nterm, neig, &
                            wid, wht, d, norder, ndim)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nelem, nterm, neig, norder, ndim
    REAL, INTENT(IN) :: p(ndim, nelem), wid(nelem), wht(nelem), d(norder, norder)
    REAL, INTENT(OUT) :: res_x(ndim, nelem), res_y(ndim, nelem)

    INTEGER :: ne, i, j, n, ij, kk
    REAL :: facx, facy, dx, dy
    INTEGER, DIMENSION(nterm) :: li

    DO n = 1, nterm
      li(n) = (n - 1) * nterm
    END DO

    res_x = 0.0
    res_y = 0.0

    DO ne = 1, nelem
      facx = 2.0 / wid(ne)
      facy = 2.0 / wht(ne)
      ! Derivatives in x-direction
      DO n = 1, nterm
        DO i = 1, nterm
          DO j = 1, nterm
            ij = li(i) + j
            kk = li(n) + j
            dx = d(i, n) * facx
            res_x(ij, ne) = res_x(ij, ne) + dx * p(kk, ne)
          END DO
        END DO
      END DO
      ! Derivatives in y-direction
      DO n = 1, nterm
        DO i = 1, nterm
          DO j = 1, nterm
            ij = li(i) + j
            kk = li(i) + n
            dy = d(j, n) * facy
            res_y(ij, ne) = res_y(ij, ne) + dy * p(kk, ne)
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE apply_gradient

  ! ============================================================================
  ! apply_divergence
  !
  ! Computes the divergence of a vector field (u, v).
  ! Operation: res = D * [u; v]
  ! where D is the discrete divergence operator (transpose of gradient).
  ! ============================================================================
  SUBROUTINE apply_divergence(res, u, v, nelem, nterm, neig, &
                              wid, wht, d, norder, ndim)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nelem, nterm, neig, norder, ndim
    REAL, INTENT(IN) :: u(ndim, nelem), v(ndim, nelem), wid(nelem), wht(nelem), d(norder, norder)
    REAL, INTENT(OUT) :: res(ndim, nelem)

    REAL, DIMENSION(ndim, nelem) :: dudx, dvdy
    CALL apply_gradient(dudx, res, u, nelem, nterm, neig, wid, wht, d, norder, ndim) ! res is dummy
    CALL apply_gradient(res, dvdy, v, nelem, nterm, neig, wid, wht, d, norder, ndim) ! res is dummy
    
    res = dudx + dvdy
  END SUBROUTINE apply_divergence

END MODULE fem_operators
