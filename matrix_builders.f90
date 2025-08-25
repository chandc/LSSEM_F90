!**********************************************************************
!> @module matrix_builders
!> @brief Contains subroutines to build elemental matrices for FEM.
!>
!> @details This module provides the logic for constructing the fundamental
!>          elemental matrices (Mass, Stiffness/Gradient) used in the
!>          spectral element method. These matrices are typically built
!>          once at the beginning of a simulation.
!**********************************************************************
MODULE matrix_builders
  IMPLICIT NONE
CONTAINS

!**********************************************************************
!> @brief Builds the diagonal elemental mass matrix for each element.
!>
!> @details In spectral methods with GLL quadrature, the mass matrix
!>          is diagonal ("lumped"), which is computationally efficient.
!>          This subroutine computes the diagonal entries, which are the
!>          product of the Jacobian and the quadrature weights.
!>
!> @param[out] M_e      The diagonal of the mass matrix for each element.
!> @param[in]  nelem    Number of elements.
!> @param[in]  nterm    Number of terms (polynomial order + 1).
!> @param[in]  ndim     Total degrees of freedom per element.
!> @param[in]  wid,wht  Width and height of each physical element.
!> @param[in]  wg       Gauss-Lobatto-Legendre quadrature weights.
!**********************************************************************
      SUBROUTINE build_mass_matrix(M_e, nelem, nterm, ndim, wid, wht, wg)
      IMPLICIT NONE

      ! Argument declarations
      INTEGER, INTENT(IN) :: nelem, nterm, ndim
      REAL, INTENT(IN) :: wid(nelem), wht(nelem), wg(nterm)
      REAL, INTENT(OUT) :: M_e(ndim, nelem)

      ! Local variables
      INTEGER :: ne, i, j, ij
      REAL :: ajac, facem
      INTEGER :: li(nterm)

      DO i = 1, nterm
        li(i) = (i - 1) * nterm
      END DO

      DO ne = 1, nelem
        ajac = 0.25 * wid(ne) * wht(ne)
        DO i = 1, nterm
          DO j = 1, nterm
            ij = li(i) + j
            facem = ajac * wg(i) * wg(j)
            M_e(ij, ne) = facem
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE build_mass_matrix

!**********************************************************************
!> @brief Builds the elemental stiffness matrices (discrete gradient).
!>
!> @details This subroutine constructs the dense elemental stiffness
!>          matrices Gx_e and Gy_e, which represent the partial
!>          derivatives with respect to x and y. These matrices are
!>          built from the 1D differentiation matrix 'd'.
!>
!> @param[out] Gx_e,Gy_e The elemental stiffness matrices for each element.
!> @param[in]  nelem     Number of elements.
!> @param[in]  nterm     Number of terms (polynomial order + 1).
!> @param[in]  ndim      Total degrees of freedom per element.
!> @param[in]  norder    Order of the differentiation matrix.
!> @param[in]  wid,wht   Width and height of each physical element.
!> @param[in]  d         The 1D differentiation matrix.
!**********************************************************************
      SUBROUTINE build_stiffness_matrix(Gx_e, Gy_e, nelem, nterm, ndim, &
                                        norder, wid, wht, d)
      IMPLICIT NONE

      ! Argument declarations
      INTEGER, INTENT(IN) :: nelem, nterm, ndim, norder
      REAL, INTENT(IN) :: wid(nelem), wht(nelem), d(norder, norder)
      REAL, INTENT(OUT) :: Gx_e(ndim, ndim, nelem)
      REAL, INTENT(OUT) :: Gy_e(ndim, ndim, nelem)

      ! Local variables
      INTEGER :: ne, i, j, n, row, col
      REAL :: facx, facy
      INTEGER :: li(nterm)

      DO i = 1, nterm
        li(i) = (i - 1) * nterm
      END DO

      DO ne = 1, nelem
        facx = 2.0 / wid(ne)
        facy = 2.0 / wht(ne)
        
        ! DEBUG: Check for NaN in scaling factors
        IF (facx /= facx .OR. facy /= facy) THEN
          PRINT *, 'WARNING: NaN in scaling factors - facx=', facx, ' facy=', facy
          PRINT *, 'WARNING: wid(', ne, ')=', wid(ne), ' wht(', ne, ')=', wht(ne)
        END IF
        
        Gx_e(:,:,ne) = 0.0
        Gy_e(:,:,ne) = 0.0

        DO i = 1, nterm
          DO j = 1, nterm
            row = li(i) + j

            ! Build Gx_e (x-derivative)
            DO n = 1, nterm
              col = li(n) + j
              Gx_e(row, col, ne) = d(i, n) * facx
              
              ! ROBUST: Clip extreme values to prevent Infinity/NaN propagation
              IF (ABS(Gx_e(row, col, ne)) > 1.0e15) THEN
                Gx_e(row, col, ne) = SIGN(1.0e15, Gx_e(row, col, ne))
                IF (ne == 1) PRINT *, 'WARNING: Clipped Gx_e(', row, ',', col, ',1) from d(', i, ',', n, ')=', d(i, n)
              END IF
            END DO

            ! Build Gy_e (y-derivative)
            DO n = 1, nterm
              col = li(i) + n
              Gy_e(row, col, ne) = d(j, n) * facy
              
              ! ROBUST: Clip extreme values to prevent Infinity/NaN propagation
              IF (ABS(Gy_e(row, col, ne)) > 1.0e15) THEN
                Gy_e(row, col, ne) = SIGN(1.0e15, Gy_e(row, col, ne))
                IF (ne == 1) PRINT *, 'WARNING: Clipped Gy_e(', row, ',', col, ',1) from d(', j, ',', n, ')=', d(j, n)
              END IF
            END DO
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE build_stiffness_matrix

END MODULE matrix_builders
