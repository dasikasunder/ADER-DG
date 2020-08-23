! surface_integral.f90
! author: sunder

!-----------------------------------------------------------------------
! Find the surface integral
!-----------------------------------------------------------------------

SUBROUTINE surface_integral(lduh, lFl, lFr, lGb, lGt)

    USE ader_dg
    IMPLICIT NONE
    ! Argument list
    DOUBLE PRECISION, INTENT(IN)    :: lFl(nVar,nDOF(2)), lFr(nVar,nDOF(2)), lGb(nVar,nDOF(1)), lGt(nVar,nDOF(1))
    DOUBLE PRECISION, INTENT(INOUT) :: lduh(nVar, nDOF(1), nDOF(2))       ! Spatial degrees of freedom

    ! Local variables
    INTEGER           :: i,j,iVar

    ! Now multiply the numerical fluxes on the surfaces with the test functions and compute the surface integrals

    ! x faces

    DO j = 1, nDOF(2)
        DO iVar = 1, nVar
            lduh(iVar,:,j) = lduh(iVar,:,j) - &
                wGPN(j)/dx(1)*( lFr(iVar,j)*FRCoeff - lFl(iVar,j)*FLCoeff )      ! left flux minus right flux
        ENDDO
    ENDDO

    ! y faces

    DO i = 1, nDOF(1)
        DO iVar = 1, nVar
            lduh(iVar,i,:) = lduh(iVar,i,:) - &
            wGPN(i)/dx(2)*( lGt(iVar,i)*FRCoeff - lGb(iVar,i)*FLCoeff )         ! bottom flux minus top flux
        ENDDO
    ENDDO

END SUBROUTINE surface_integral
