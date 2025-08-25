!**********************************************************************
!> @module elemental_operators
!> @brief Contains subroutines to apply elemental matrices to vectors.
!>
!> @details This module provides the computational kernels for applying
!>          the pre-built elemental matrices (Mass, Stiffness) to
!>          elemental data vectors. These operations form the core of
!>          the FEM calculation within the main solver loop.
!**********************************************************************
MODULE elemental_operators
  IMPLICIT NONE
CONTAINS

!**********************************************************************
!> @brief Applies the diagonal mass matrix to an elemental vector.
!**********************************************************************
      SUBROUTINE apply_mass_operator(res_e, M_e, u_e, ndim)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ndim
      REAL, INTENT(IN) :: M_e(ndim), u_e(ndim)
      REAL, INTENT(OUT) :: res_e(ndim)
      res_e = M_e * u_e
      RETURN
      END SUBROUTINE apply_mass_operator

!**********************************************************************
!> @brief Applies the stiffness matrices (gradient) to an elemental vector.
!**********************************************************************
      SUBROUTINE apply_gradient_operator(grad_x_e, grad_y_e, Gx_e, Gy_e, p_e, ndim)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ndim
      REAL, INTENT(IN) :: Gx_e(ndim, ndim), Gy_e(ndim, ndim), p_e(ndim)
      REAL, INTENT(OUT) :: grad_x_e(ndim), grad_y_e(ndim)
      ! For performance, replace with calls to BLAS DGEMV
      grad_x_e = MATMUL(Gx_e, p_e)
      grad_y_e = MATMUL(Gy_e, p_e)
      RETURN
      END SUBROUTINE apply_gradient_operator

!**********************************************************************
!> @brief Applies the transpose of the discrete differential operator.
!**********************************************************************
      SUBROUTINE apply_transpose_operator(c1_e, c2_e, c3_e, c4_e, &
                     su, fu_e, fv_e, dfudx, dfudy, dfvdx, dfvdy, &
                     fac1, dt, pr, ndim, Gx_e, Gy_e)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ndim
      REAL, INTENT(IN) :: fac1, dt, pr
      REAL, INTENT(IN) :: su(ndim, 4)
      REAL, INTENT(IN) :: fu_e(ndim), fv_e(ndim)
      REAL, INTENT(IN) :: dfudx(ndim), dfudy(ndim), dfvdx(ndim), dfvdy(ndim)
      REAL, INTENT(IN) :: Gx_e(ndim, ndim), Gy_e(ndim, ndim)
      REAL, INTENT(OUT) :: c1_e(ndim), c2_e(ndim), c3_e(ndim), c4_e(ndim)

      INTEGER :: i
      REAL, DIMENSION(ndim) :: grad_t_c1, grad_t_c2, grad_t_c3, grad_t_c4

      ! Local contributions
      c1_e = (dt + fac1) * su(:,1) + dt * dfudx * su(:,1) + dt * dfvdx * su(:,2)
      c2_e = (dt + fac1) * su(:,2) + dt * dfudy * su(:,1) + dt * dfvdy * su(:,2)
      c3_e = 0.0
      c4_e = su(:,4)

      ! Transpose-gradient part (using transposed matrices)
      ! grad_t_c1 = Gx_e' * su(:,3) + Gy_e' * su(:,4)
      grad_t_c1 = MATMUL(TRANSPOSE(Gx_e), su(:,3)) + MATMUL(TRANSPOSE(Gy_e), su(:,4))
      c1_e = c1_e + grad_t_c1

      ! grad_t_c2 = Gx_e' * (-su(:,4)) + Gy_e' * su(:,3)
      grad_t_c2 = MATMUL(TRANSPOSE(Gx_e), -su(:,4)) + MATMUL(TRANSPOSE(Gy_e), su(:,3))
      c2_e = c2_e + grad_t_c2

      ! grad_t_c3 = Gx_e' * su(:,1) + Gy_e' * su(:,2)
      grad_t_c3 = MATMUL(TRANSPOSE(Gx_e), su(:,1)) + MATMUL(TRANSPOSE(Gy_e), su(:,2))
      c3_e = c3_e + dt * grad_t_c3

      ! grad_t_c4 = Gx_e' * (-pr*su(:,2)) + Gy_e' * (pr*su(:,1))
      grad_t_c4 = MATMUL(TRANSPOSE(Gx_e), -pr*su(:,2)) + MATMUL(TRANSPOSE(Gy_e), pr*su(:,1))
      c4_e = c4_e + dt * grad_t_c4

      END SUBROUTINE apply_transpose_operator

END MODULE elemental_operators
