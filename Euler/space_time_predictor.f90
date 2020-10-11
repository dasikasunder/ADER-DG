! space_time_predictor.f90
! author: sunder

!-----------------------------------------------------------------------
! Space time predictor function
!-----------------------------------------------------------------------

subroutine ader_space_time_predictor(lqhi, lFhi, lQbnd, lFbnd, luh)

    use ader_dg
    implicit none
    ! Argument list
    real, intent(in)  :: luh(nVar,nDOF(1),nDOF(2))              ! spatial degrees of freedom
    real, intent(out) :: lqhi(nVar,nDOF(1),nDOF(2))             ! time-averaged space-time degrees of freedom
    real, intent(out) :: lFhi(nVar,nDim,nDOF(1),nDOF(2))        ! time-averaged nonlinear flux tensor in each space-time DOF
    real, intent(out) :: lqbnd(nVar,4,nDOF(2))                  ! time-averaged space-time degrees of freedom
    real, intent(out) :: lFbnd(nVar,4,nDOF(2))                  ! time-averaged nonlinear flux tensor in each space-time DOF

    ! Local variables
    integer  :: i,j,l,iVar,iDim, iter
    real :: rhs0(nVar,nDOF(1),nDOF(2),nDOF(0))      ! contribution of the initial condition to the known right hand side
    real :: rhs(nVar,nDOF(1),nDOF(2),nDOF(0))       ! known right hand side
    real :: lqh(nVar,nDOF(1),nDOF(2),nDOF(0))       ! space-time degrees of freedom
    real :: lFh(nVar,nDim,nDOF(1),nDOF(2),nDOF(0))  ! nonlinear flux tensor in each space-time DOF
    real :: lqhold(nVar,nDOF(1),nDOF(2),nDOF(0))    ! old space-time degrees of freedom

    ! Form an initial guess

    do j = 1, nDOF(2)
        do i = 1, nDOF(1)

            ! Trivial initial guess

            do iVar = 1, nVar
                lqh(iVar,i,j,:) = luh(iVar,i,j)
            end do

            ! Compute the contribution of the initial condition uh to the time update.

            do iVar = 1, nVar
                rhs0(iVar,i,j,:) = wGPN(i)*wGPN(j)*F0(:)*luh(iVar,i,j)
            end do

        end do
    end do

    ! Discrete Picard iterations.

    do iter = 1, N+1
        ! save old space-time DOF
        lqhold = lqh
        do l = 1, nDOF(0) ! loop over DOF in time

            ! Compute the fluxes (once these fluxes are available, the subsequent operations are independent from each other)
            do j = 1, nDOF(2)
                do i = 1, nDOF(1)
                    call PDEFlux(lqh(:,i,j,l), lFh(:,:,i,j,l))
                end do
            end do

            ! Compute the "derivatives" (contributions of the stiffness matrix)

            ! x direction (independent from the y derivatives)

            do j = 1, nDOF(2)
                rhs(:,:,j,l) = rhs0(:,:,j,l) - wGPN(l)*wGPN(j)*dt/dx(1)*matmul( lFh(:,1,:,j,l), Kxi )
            end do

            ! y direction (independent from the x derivatives)

            do i = 1, nDOF(1)
                rhs(:,i,:,l) = rhs(:,i,:,l) - wGPN(l)*wGPN(i)*dt/dx(2)*matmul( lFh(:,2,i,:,l), Kxi )
            end do

        end do ! end loop over time DOF

        ! Multiply with (K1)^(-1) to get the discrete time integral of the discrete Picard iteration

        do j = 1, nDOF(2)
            do i = 1, nDOF(1)
                lqh(:,i,j,:) = 1.0/(wGPN(i)*wGPN(j))*matmul( rhs(:,i,j,:), transpose(iK1) )
            end do
        end do
    end do

    ! Immediately compute the time-averaged space-time polynomials

    do j = 1, nDOF(2)
        do i = 1, nDOF(1)
            lqhi(:,i,j) = matmul( lqh(:,i,j,:), wGPN )
            do iDim = 1, nDim
                lFhi(:,iDim,i,j) = matmul( lFh(:,iDim,i,j,:), wGPN )
            end do
        end do
    end do


    ! Compute the bounday-extrapolated values for Q and F*n

    lQbnd = 0.0
    lFbnd = 0.0

    ! x-direction: face 1 (left) and face 2 (right)

    do j = 1, nDOF(2)
        lQbnd(:,1,j) = matmul( lqhi(:,:,j),   FLCoeff )   ! left
        lQbnd(:,2,j) = matmul( lqhi(:,:,j),   FRCoeff )   ! right
        lFbnd(:,1,j) = matmul( lFhi(:,1,:,j), FLCoeff )   ! left
        lFbnd(:,2,j) = matmul( lFhi(:,1,:,j), FRCoeff )   ! right
    end do

    ! y-direction: face 3 (left) and face 4 (right)

    do i = 1, nDOF(1)
        lQbnd(:,3,i) = matmul( lqhi(:,i,:),   FLCoeff )   ! left
        lQbnd(:,4,i) = matmul( lqhi(:,i,:),   FRCoeff )   ! right
        lFbnd(:,3,i) = matmul( lFhi(:,2,i,:), FLCoeff )   ! left
        lFbnd(:,4,i) = matmul( lFhi(:,2,i,:), FRCoeff )   ! right
    end do

    continue

end subroutine ader_space_time_predictor
