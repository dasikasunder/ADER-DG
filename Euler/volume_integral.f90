! volume_integral.f90
! author: sunder

!-----------------------------------------------------------------------
! Find the element volume integral
!-----------------------------------------------------------------------

SUBROUTINE volume_integral(lduh, lqhi, lFhi)
        USE ader_dg
        IMPLICIT NONE
        ! Argument list
        DOUBLE PRECISION, INTENT(IN)  :: lqhi(nVar,nDOF(1),nDOF(2))      ! space-time degrees of freedom
        DOUBLE PRECISION, INTENT(IN)  :: lFhi(nVar,nDim,nDOF(1),nDOF(2)) ! nonlinear flux tensor in each space-time DOF
        DOUBLE PRECISION, INTENT(OUT) :: lduh(nVar,nDOF(1),nDOF(2))      ! spatial degrees of freedom

        ! Local variables
        INTEGER           :: i,j

        ! Initialize the update DOF

        lduh = 0.0d0

        ! x - direction

        DO j = 1, nDOF(2)
            lduh(:,:,j) = lduh(:,:,j) + MATMUL( lFhi(:,1,:,j), TRANSPOSE(Kxi) )*wGPN(j)/dx(1)
        ENDDO

        ! y - direction

        DO i = 1, nDOF(1)
            lduh(:,i,:) = lduh(:,i,:) + MATMUL( lFhi(:,2,i,:), TRANSPOSE(Kxi) )*wGPN(i)/dx(2)
        ENDDO

        CONTINUE

END SUBROUTINE volume_integral
