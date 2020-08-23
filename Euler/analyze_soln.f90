! analyze_soln.f90
! author: sunder

!-----------------------------------------------------------------------
! Analyzes the solution and prints L_inf error compared to the initial
! condition. Useful for smooth periodic test cases to check convergence
! rate
!-----------------------------------------------------------------------

SUBROUTINE analyze_soln
    USE ader_dg
    IMPLICIT NONE

    DOUBLE PRECISION, POINTER      :: uh_exact(:,:,:,:,:), error_density(:, :)
    INTEGER :: i, j, k, l
    DOUBLE PRECISION :: u0(nVar), xGP(nDim), corner(nDim), rho_ave, rho_ave_exact

    ALLOCATE(  uh_exact(nVar, nDOF(1), nDOF(2), IMAX, JMAX) )
    ALLOCATE(  error_density(IMAX, JMAX) )

    DO j = 1, JMAX
        DO i = 1, IMAX
            corner(:) = x(:,i,j) - 0.5d0*dx(:) ! coordinate of the lower left corner of the cell
            DO l = 1, nDOF(2)
                DO k = 1, nDOF(1)
                    xGP = corner + (/ xiGPN(k), xiGPN(l) /)*dx(:)
                    CALL InitialField(xGP, u0)
                    uh_exact(:, k, l, i, j) = u0(:)
                END DO
            END DO
        END DO
    END DO

    DO j = 1, JMAX
        DO i = 1, IMAX

            rho_ave= 0.0d0
            rho_ave_exact = 0.0d0

            DO l = 1, nDOF(2)
                DO k = 1, nDOF(1)
                    rho_ave_exact = rho_ave_exact + wGPN(l)*wGPN(k)*uh_exact(1, k, l, i, j)
                    rho_ave = rho_ave + wGPN(l)*wGPN(k)*uh(1, k, l, i, j)
                END DO
            END DO

            error_density(i,j) = ABS(rho_ave_exact - rho_ave)

        END DO
    END DO

    PRINT *, MAXVAL(error_density)

    DEALLOCATE( uh_exact, error_density)

END SUBROUTINE analyze_soln

