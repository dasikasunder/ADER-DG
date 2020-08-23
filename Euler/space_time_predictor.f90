! space_time_predictor.f90
! author: sunder

!-----------------------------------------------------------------------
! Space time predictor function
!-----------------------------------------------------------------------

SUBROUTINE ader_space_time_predictor(lqhi, lFhi, lQbnd, lFbnd, luh)

    USE ader_dg
    IMPLICIT NONE
    ! Argument list
    DOUBLE PRECISION, INTENT(IN)  :: luh(nVar,nDOF(1),nDOF(2))              ! spatial degrees of freedom
    DOUBLE PRECISION, INTENT(OUT) :: lqhi(nVar,nDOF(1),nDOF(2))             ! time-averaged space-time degrees of freedom
    DOUBLE PRECISION, INTENT(OUT) :: lFhi(nVar,nDim,nDOF(1),nDOF(2))        ! time-averaged nonlinear flux tensor in each space-time DOF
    DOUBLE PRECISION, INTENT(OUT) :: lqbnd(nVar,4,nDOF(2))                  ! time-averaged space-time degrees of freedom
    DOUBLE PRECISION, INTENT(OUT) :: lFbnd(nVar,4,nDOF(2))                  ! time-averaged nonlinear flux tensor in each space-time DOF

    ! Local variables
    INTEGER  :: i,j,l,iVar,iDim, iter
    DOUBLE PRECISION :: rhs0(nVar,nDOF(1),nDOF(2),nDOF(0))      ! contribution of the initial condition to the known right hand side
    DOUBLE PRECISION :: rhs(nVar,nDOF(1),nDOF(2),nDOF(0))       ! known right hand side
    DOUBLE PRECISION :: lqh(nVar,nDOF(1),nDOF(2),nDOF(0))       ! space-time degrees of freedom
    DOUBLE PRECISION :: lFh(nVar,nDim,nDOF(1),nDOF(2),nDOF(0))  ! nonlinear flux tensor in each space-time DOF
    DOUBLE PRECISION :: lqhold(nVar,nDOF(1),nDOF(2),nDOF(0))    ! old space-time degrees of freedom

    ! Form an initial guess

    DO j = 1, nDOF(2)
        DO i = 1, nDOF(1)

            ! Trivial initial guess

            DO iVar = 1, nVar
                lqh(iVar,i,j,:) = luh(iVar,i,j)
            ENDDO

            ! Compute the contribution of the initial condition uh to the time update.

            DO iVar = 1, nVar
                rhs0(iVar,i,j,:) = wGPN(i)*wGPN(j)*F0(:)*luh(iVar,i,j)
            ENDDO

        ENDDO
        ENDDO

    ! Discrete Picard iterations.

    DO iter = 1, N+1
        ! save old space-time DOF
        lqhold = lqh
        DO l = 1, nDOF(0) ! loop over DOF in time

            ! Compute the fluxes (once these fluxes are available, the subsequent operations are independent from each other)
            DO j = 1, nDOF(2)
                DO i = 1, nDOF(1)
                    CALL PDEFlux(lqh(:,i,j,l), lFh(:,:,i,j,l))
                ENDDO
            ENDDO

            ! Compute the "derivatives" (contributions of the stiffness matrix)

            ! x direction (independent from the y derivatives)

            DO j = 1, nDOF(2)
                rhs(:,:,j,l) = rhs0(:,:,j,l) - wGPN(l)*wGPN(j)*dt/dx(1)*MATMUL( lFh(:,1,:,j,l), Kxi )
            ENDDO

            ! y direction (independent from the x derivatives)

            DO i = 1, nDOF(1)
                rhs(:,i,:,l) = rhs(:,i,:,l) - wGPN(l)*wGPN(i)*dt/dx(2)*MATMUL( lFh(:,2,i,:,l), Kxi )
            ENDDO

        ENDDO ! end loop over time DOF

        ! Multiply with (K1)^(-1) to get the discrete time integral of the discrete Picard iteration

        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                lqh(:,i,j,:) = 1.0d0/(wGPN(i)*wGPN(j))*MATMUL( rhs(:,i,j,:), TRANSPOSE(iK1) )
            ENDDO
        ENDDO
    ENDDO

    ! Immediately compute the time-averaged space-time polynomials

    DO j = 1, nDOF(2)
        DO i = 1, nDOF(1)
            lqhi(:,i,j) = MATMUL( lqh(:,i,j,:), wGPN )
            DO iDim = 1, nDim
                lFhi(:,iDim,i,j) = MATMUL( lFh(:,iDim,i,j,:), wGPN )
            ENDDO
        ENDDO
    ENDDO


    ! Compute the bounday-extrapolated values for Q and F*n

    lQbnd = 0.0d0
    lFbnd = 0.0d0

    ! x-direction: face 1 (left) and face 2 (right)

    DO j = 1, nDOF(2)
        lQbnd(:,1,j) = MATMUL( lqhi(:,:,j),   FLCoeff )   ! left
        lQbnd(:,2,j) = MATMUL( lqhi(:,:,j),   FRCoeff )   ! right
        lFbnd(:,1,j) = MATMUL( lFhi(:,1,:,j), FLCoeff )   ! left
        lFbnd(:,2,j) = MATMUL( lFhi(:,1,:,j), FRCoeff )   ! right
    ENDDO

    ! y-direction: face 3 (left) and face 4 (right)

    DO i = 1, nDOF(1)
        lQbnd(:,3,i) = MATMUL( lqhi(:,i,:),   FLCoeff )   ! left
        lQbnd(:,4,i) = MATMUL( lqhi(:,i,:),   FRCoeff )   ! right
        lFbnd(:,3,i) = MATMUL( lFhi(:,2,i,:), FLCoeff )   ! left
        lFbnd(:,4,i) = MATMUL( lFhi(:,2,i,:), FRCoeff )   ! right
        ENDDO

    CONTINUE

END SUBROUTINE ader_space_time_predictor
