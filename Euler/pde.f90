! pde.f90
! author: sunder

!-----------------------------------------------------------------------
! Routines related to the given PDE
!-----------------------------------------------------------------------

subroutine PDECons2Prim(Q, V)
    use ader_dg, only : nVar, nDim, EQN
    implicit none
    ! Argument list
    double precision, intent(IN)     :: Q(nVar)     ! vector of conserved quantities
    double precision, intent(OUT)    :: V(nVar)     ! primitive variables

    ! Local variables
    double precision                 :: p

    p = (EQN%gamma-1.0d0)*( Q(4) - 0.5d0*SUM(Q(2:3)**2)/Q(1) )    ! fluid pressure
    V(1) = Q(1)             ! fluid density
    V(2:3) = Q(2:3)/Q(1)    ! fluid velocity
    V(4)   = p              ! fluid pressure

end subroutine PDECons2Prim

!-----------------------------------------------------------------------
! Convert primitive vector to conserved vector
!-----------------------------------------------------------------------

subroutine PDEPrim2Cons(V, Q)
    use ader_dg, ONLY : nVar, nDim, EQN
    implicit none
    ! Argument list
    double precision, intent(in)      :: V(nVar)     ! primitive variables
    double precision, intent(out)     :: Q(nVar)     ! vector of conserved quantities

    Q(1)   = V(1)           ! fluid density
    Q(2:3) = V(1)*V(2:3)    ! momentum
    Q(4)   = V(4)/(EQN%GAMMA-1) + 0.5d0*V(1)*SUM(V(2:3)**2)   ! total energy = internal energy + kinetic energy

END subroutine PDEPrim2Cons

!-----------------------------------------------------------------------
! PDE flux tensor F(Q)
!-----------------------------------------------------------------------

subroutine PDEFlux(Q, F)
    use ader_dg, ONLY : nVar, nDim, EQN
    implicit none
    ! Argument list
    double precision, intent(in)  :: Q(nVar)
    double precision, intent(out) :: F(nVar,nDim)

    ! Local variables
    double precision :: p, irho

    irho = 1.0d0/Q(1)
    p = (EQN%GAMMA-1.0d0)*( Q(4) - 0.5d0*SUM(Q(2:3)**2)*irho )

    !IF(Q(1) .LT. 1.0d0-14) THEN
        !PRINT *, 'Negative Density', Q(1)
        !STOP 1
    !ENDIF

    !IF(p .LT. 1.0d0-12) THEN
        !PRINT *, 'Negative Pressure', p
        !STOP 1
    !ENDIF

    F(1,1) = Q(2)
    F(2,1) = irho*Q(2)*Q(2) + p
    F(3,1) = irho*Q(2)*Q(3)
    F(4,1) = irho*Q(2)*(Q(4) + p)

    F(1,2) = Q(3)
    F(2,2) = irho*Q(3)*Q(2)
    F(3,2) = irho*Q(3)*Q(3) + p
    F(4,2) = irho*Q(3)*(Q(4) + p)

end subroutine PDEFlux

!-----------------------------------------------------------------------
! Eigenvalues of the PDE in the normal direction nv
!-----------------------------------------------------------------------

subroutine PDEEigenvalues(Q, nv, Lambda)
    use ader_dg, ONLY : nVar, nDim, EQN
    implicit none
    ! Argument list
    double precision, intent(in)  :: Q(nVar), nv(nDim)
    double precision, intent(out) :: Lambda(nVar)

    ! Local variables
    double precision :: p, u, c

    u = ( Q(2)*nv(1) + Q(3)*nv(2) )/Q(1)                      ! normal velocity
    p = (EQN%gamma-1.0d0)*( Q(4) - 0.5d0*SUM(Q(2:3)**2)/Q(1) )    ! fluid pressure
    c = SQRT(EQN%gamma*p/Q(1))                                ! sound speed

    Lambda = (/ u-c, u, u, u+c /)                          ! The eigenvalues of the Euler equations

end subroutine PDEEigenvalues

!-----------------------------------------------------------------------
! LLF (Rusanov) Riemann solver
!-----------------------------------------------------------------------

subroutine RusanovFlux(QL, FL, QR, FR, nv, Flux)
    use ader_dg
    implicit none
    ! Argument list
    double precision, intent(in)  :: QL(nVar), FL(nVar), QR(nVar), FR(nVar), nv(nDim)
    double precision, intent(out) :: Flux(nVar)

    ! Local variables
    double precision :: smax, LL(nVar), LR(nVar)

    call PDEEigenvalues(QL,nv,LL)
    call PDEEigenvalues(QR,nv,LR)

    smax = max( maxval(abs(LL)), maxval(abs(LR)) )

    Flux = 0.5d0*(FL + FR) - 0.5d0*smax*(QR - QL)

end subroutine

!-----------------------------------------------------------------------
! Make sure the input state is physically valid
!-----------------------------------------------------------------------

subroutine PDEAssurePositivity(iErr,Q, Qout)
    use ader_dg, only : nVar, nDim, EQN
    implicit none
    ! Argument list
    double precision, intent(IN)     :: Q(nVar)
    double precision, intent(OUT)    :: Qout(nVar)
    integer, intent(OUT) :: iErr
    ! Local variables
    double precision :: p, irho

    iErr = 0
    Qout = Q ! we do not apply dirty tricks here at this point
    irho = 1.0d0/Q(1)

    ! 2D compressible Euler equations

    p = (EQN%GAMMA-1)*( Q(4) - 0.50d0*SUM(Q(2:3)**2)*irho )

    if ( Q(1) .le. 1.0d0-12) then
        iErr = -1
    end if

    if (p .le. 1.0d0-12) then
        iErr = -2
    end if
end subroutine PDEAssurePositivity
