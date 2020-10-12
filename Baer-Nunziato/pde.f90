! pde.f90
! author: sunder

!-----------------------------------------------------------------------
! Routines related to the given PDE
!-----------------------------------------------------------------------

subroutine PDECons2Prim(Q, V)
    use ader_dg, only : nVar, EQN
    implicit none
    ! Argument list
    real, intent(IN)     :: Q(nVar)     ! vector of conserved quantities
    real, intent(OUT)    :: V(nVar)     ! primitive variables

    ! Local variables
    real                 :: phi_s, phi_g, rho_s, rho_g, u_s, u_g, v_s, v_g, E_s, E_g, p_s, p_g

    phi_s = Q(9)
    phi_g = 1.0 - Q(9)

    rho_s = Q(1)/phi_s
    rho_g = Q(5)/phi_g

    u_s = Q(2)/Q(1); v_s = Q(3)/Q(1);
    u_g = Q(6)/Q(5); v_g = Q(7)/Q(5);

    E_s = Q(4)/phi_s
    E_g = Q(8)/phi_g

    p_s = (EQN%GAMMA_S - 1.0)*( E_s - 0.5*rho_s*(u_s*u_s + v_s*v_s) ) - EQN%GAMMA_S*EQN%PI_S
    p_g = (EQN%GAMMA_G - 1.0)*( E_g - 0.5*rho_g*(u_g*u_g + v_g*v_g) ) - EQN%GAMMA_G*EQN%PI_G

    V(1) = rho_s
    V(2) = u_s
    V(3) = v_s
    V(4) = p_s
    V(5) = rho_g
    V(6) = u_g
    V(7) = v_g
    V(8) = p_g
    V(9) = phi_s

end subroutine PDECons2Prim

!-----------------------------------------------------------------------
! Convert primitive vector to conserved vector
!-----------------------------------------------------------------------

subroutine PDEPrim2Cons(V, Q)
    use ader_dg, ONLY : nVar, EQN
    implicit none
    ! Argument list
    real, intent(in)      :: V(nVar)     ! primitive variables
    real, intent(out)     :: Q(nVar)     ! vector of conserved quantities

    Q(1) = V(9)*V(1)     ! phi_s*rho_s
    Q(2) = Q(1)*V(2)     ! phi_s*rho_s*vx_s
    Q(3) = Q(1)*V(3)     ! phi_s*rho_s*vy_s
    Q(4) = Q(1)*(0.5*(V(2)*V(2) + V(3)*V(3)) + (V(4) + EQN%GAMMA_S*EQN%PI_S)/(V(1)*(EQN%GAMMA_S-1.0)) ) ! phi_s*E_s

    Q(5) = (1.0-V(9))*V(5) ! phi_g*rho_g
    Q(6) = Q(5)*V(6)       ! phi_g*rho_g*vx_g
    Q(7) = Q(5)*V(7)       ! phi_g*rho_g*vy_g
    Q(8) = Q(5)*( 0.5*(V(6)*V(6) + V(7)*V(7)) + (V(8) + EQN%GAMMA_G*EQN%PI_G)/(V(5)*(EQN%GAMMA_G-1.0)) ) ! phi_g*E_g

    Q(9) = V(9)            ! phi_s

END subroutine PDEPrim2Cons

!-----------------------------------------------------------------------
! PDE flux tensor F(Q) (conservative part)
!-----------------------------------------------------------------------

subroutine PDEFlux(Q, F)
    use ader_dg, ONLY : nVar, nDim, EQN
    implicit none
    ! Argument list
    real, intent(in)  :: Q(nVar)
    real, intent(out) :: F(nVar,nDim)

    ! Local variables
    real              :: phi_s, phi_g, rho_s, rho_g, u_s, u_g, v_s, v_g, E_s, E_g, p_s, p_g

    phi_s = Q(9)
    phi_g = 1.0 - Q(9)

    rho_s = Q(1)/phi_s
    rho_g = Q(5)/phi_g

    u_s = Q(2)/Q(1); v_s = Q(3)/Q(1);
    u_g = Q(6)/Q(5); v_g = Q(7)/Q(5);

    E_s = Q(4)/phi_s
    E_g = Q(8)/phi_g

    p_s = (EQN%GAMMA_S - 1.0)*( E_s - 0.5*rho_s*(u_s*u_s + v_s*v_s) ) - EQN%GAMMA_S*EQN%PI_S
    p_g = (EQN%GAMMA_G - 1.0)*( E_g - 0.5*rho_g*(u_g*u_g + v_g*v_g) ) - EQN%GAMMA_G*EQN%PI_G

    F(1,1) = phi_s*rho_s*u_s
    F(2,1) = phi_s*(rho_s*u_s*u_s + p_s)
    F(3,1) = phi_s*rho_s*u_s*v_s
    F(4,1) = phi_s*u_s*(E_s + p_s)
    F(5,1) = phi_g*rho_g*u_g
    F(6,1) = phi_g*(rho_g*u_g*u_g + p_g)
    F(7,1) = phi_g*rho_g*u_g*v_g
    F(8,1) = phi_g*u_g*(E_g + p_g)
    F(9,1) = 0.0

    F(1,2) = phi_s*rho_s*v_s
    F(2,2) = phi_s*rho_s*u_s*v_s
    F(3,2) = phi_s*(rho_s*v_s*v_s + p_s)
    F(4,2) = phi_s*v_s*(E_s + p_s)
    F(5,2) = phi_g*rho_g*v_g
    F(6,2) = phi_g*rho_g*u_g*v_g
    F(7,2) = phi_g*(rho_g*v_g*v_g + p_g)
    F(8,2) = phi_g*v_g*(E_g + p_g)
    F(9,2) = 0.0

end subroutine PDEFlux

!-----------------------------------------------------------------------
! Eigenvalues of the PDE in the normal direction nv
!-----------------------------------------------------------------------

subroutine PDEEigenvalues(Q, nv, Lambda)
    use ader_dg, ONLY : nVar, nDim, EQN
    implicit none
    ! Argument list
    real, intent(in)  :: Q(nVar), nv(nDim)
    real, intent(out) :: Lambda(nVar)

    ! Local variables
    real              :: V(nVar), u_s, u_g, a_s, a_g

    call PDECons2Prim(Q, V)

    u_s = V(2)*nv(1) + V(3)*nv(2)                   ! normal velocity (solid)
    u_g = V(6)*nv(1) + V(7)*nv(2)                   ! normal velocity (gas)
    a_s = sqrt(EQN%GAMMA_S*(V(4) + EQN%PI_S)/V(1))  ! Sound speed (solid)
    a_g = sqrt(EQN%GAMMA_G*(V(8) + EQN%PI_G)/V(5))  ! Sound speed (gas)

    Lambda = (/ u_g-a_g, u_g+a_g, u_s-a_s, u_s+a_s, u_g, u_g, u_s, u_s, u_s/)

end subroutine PDEEigenvalues

!-----------------------------------------------------------------------
! Non-conservative matrix B in the normal direction
!-----------------------------------------------------------------------

subroutine PDEMatrixB(Bn, Q, nv)
    use ader_dg, only : nVar, nDim, EQN
    implicit none
    ! Argument list
    real, intent(in)  :: Q(nVar), nv(nDim)
    real, intent(out) :: Bn(nVar,nVar)
    ! Local variables
    real :: u_s, v_s, rho_g, u_g, v_g, E_g, phi_g, p_g, u_n

    ! Make standard assumption. Interface pressure is gas pressure and interface velocity is solid velocity

    Bn = 0.0

    phi_g = 1.0 - Q(9)

    rho_g = Q(5)/phi_g

    u_s = Q(2)/Q(1)
    v_s = Q(3)/Q(1)
    u_n = nv(1)*u_s + nv(2)*v_s

    u_g = Q(6)/Q(5)
    v_g = Q(7)/Q(5)

    E_g = Q(8)/phi_g

    p_g = (EQN%GAMMA_G - 1.0)*( E_g - 0.5*rho_g*(u_g*u_g + v_g*v_g) ) - EQN%GAMMA_G*EQN%PI_G

    Bn(2,9) = -nv(1)*p_g
    Bn(3,9) = -nv(2)*p_g
    Bn(4,9) = -p_g*u_n
    Bn(6,9) =  nv(1)*p_g
    Bn(7,9) =  nv(2)*p_g
    Bn(8,9) =  p_g*u_n
    Bn(9,9) =  u_n

 end subroutine PDEMatrixB

!-----------------------------------------------------------------------
! Nonconservative part of the PDE ( B(Q) * gradQ )
!-----------------------------------------------------------------------

subroutine PDENCP(BgradQ, Q, gradQ)
    use ader_dg, only : nVar, nDim, EQN
    implicit none
    ! Argument list
    real, intent(in)  :: Q(nVar), gradQ(nVar,nDim)
    real, intent(out) :: BgradQ(nVar)
    ! Local variables
    real :: Qx(nVar), Qy(nVar)
    real :: u_s, v_s, rho_g, u_g, v_g, E_g, phi_g, p_g

    Qx = gradQ(:,1)
    Qy = gradQ(:,2)

    ! Make standard assumption. Interface pressure is gas pressure and interface velocity is solid velocity

    phi_g = 1.0 - Q(9)

    rho_g = Q(5)/phi_g

    u_s = Q(2)/Q(1)
    v_s = Q(3)/Q(1)

    u_g = Q(6)/Q(5)
    v_g = Q(7)/Q(5)

    E_g = Q(8)/phi_g;

    p_g = (EQN%GAMMA_G - 1.0)*( E_g - 0.5*rho_g*(u_g*u_g + v_g*v_g) ) - EQN%GAMMA_G*EQN%PI_G;

    BgradQ(1) =  0.0
    BgradQ(2) = -p_g*Qx(9)
    BgradQ(3) = -p_g*Qy(9)
    BgradQ(4) = -p_g*(u_s*Qx(9) + v_s*Qy(9))
    BgradQ(5) =  0.0
    BgradQ(6) =  p_g*Qx(9)
    BgradQ(7) =  p_g*Qy(9)
    BgradQ(8) =  p_g*(u_s*Qx(9) + v_s*Qy(9))
    BgradQ(9) =  u_s*Qx(9) + v_s*Qy(9)

 end subroutine PDENCP

!-----------------------------------------------------------------------
! Make sure the input state is physically valid
!-----------------------------------------------------------------------

subroutine PDEAssurePositivity(iErr,Q, Qout)
    use ader_dg, only : nVar, EQN
    implicit none
    ! Argument list
    real, intent(in)  :: Q(nVar)
    real, intent(out) :: Qout(nVar)
    integer, intent(OUT) :: iErr
    ! Local variables
    real :: V(nVar)

    iErr = 0
    Qout = Q

    call PDECons2Prim(Q, V)

    ! 2D Baer-Nunziato equations

    if (V(1) .le. 1.0e-12 .or. V(5) .le. 1.0e-12) then ! Density
        iErr = -1
    end if

    if ( (V(4) + EQN%PI_S)  .le. 1.0e-12 .or. (V(8) + EQN%PI_G) .le. 1.0e-12) then ! Pressure
        iErr = -2
    end if

    if (V(9) .le. 0.0 .or.  V(9) .ge. 1.0) then ! Phase fraction
        iErr = -3
    end if

end subroutine PDEAssurePositivity

!---------------------------------------------------------------------------------
! No need to change anything beyond this point
!---------------------------------------------------------------------------------

!-----------------------------------------------------------------------
! Combine the source and Non conservative product term since their
! numerical treatment in many ways is very similar
!-----------------------------------------------------------------------

subroutine PDEFusedSrcNCP(Src_BgradQ, Q, gradQ)
    use ader_dg, only : nVar, nDim
    implicit none
    ! Argument list
    real, intent(in)  :: Q(nVar), gradQ(nVar,nDim)
    real, intent(out) :: Src_BgradQ(nVar)
    ! Local arguments
    real :: ncp(nVar)

    call PDENCP(ncp, Q, gradQ)

    Src_BgradQ = -ncp

end subroutine PDEFusedSrcNCP

!-----------------------------------------------------------------------
! Roe Matrix to form the path integral
!-----------------------------------------------------------------------

subroutine RoeMatrix(BRoe, Qa, Qb, nv)
    use ader_dg
    ! Argument list
    real, intent(in)  :: Qa(nVar), Qb(nVar), nv(nDim)
    real, intent(out) :: BRoe(nVar,nVar)
    ! Local variables
    integer, parameter :: nGP = 3 ! nGP = number of Gauss–Legendre points
    integer :: i
    real :: sGP(nGP), wGP(nGP)
    real :: B(nVar,nVar), Q(nVar)

    ! Definition of the Gauss–Legendre quadrature rule using 3 quadrature points
    sGP = (/ 0.5-sqrt(15.0)/10., 0.5, 0.5+sqrt(15.)/10. /)
    wGP = (/ 5./18.0, 8./18., 5./18. /)
    BRoe = 0.0

    do i = 1, nGP
        Q = Qa + sGP(i)*(Qb-Qa)
        call PDEMatrixB(B, Q, nv)
        ! Initialize Roe matrix with zero
        ! Straight-line segment path
        ! Evaluate matrix B(Q)
        BRoe = BRoe + wGP(i)*B ! Compute the path integral using numerical quadrature
    end do
end subroutine RoeMatrix

!-----------------------------------------------------------------------
! LLF (Rusanov) Riemann solver
!-----------------------------------------------------------------------

subroutine RusanovFlux(QL, FL, QR, FR, nv, Dm, Dp)
    use ader_dg
    implicit none
    ! Argument list
    real, intent(in)  :: QL(nVar), FL(nVar), QR(nVar), FR(nVar), nv(nDim)
    real, intent(out) :: Dm(nVar), Dp(nVar)

    ! Local variables
    real :: smax, LL(nVar), LR(nVar)
    real :: Bn(nVar,nVar), temp(nVar)

    call PDEEigenvalues(QL,nv,LL)
    call PDEEigenvalues(QR,nv,LR)

    smax = max( maxval(abs(LL)), maxval(abs(LR)) )

    !call PDEMatrixB(Bn, 0.50*(QL+QR), nv) ! Simpler two point average

    call RoeMatrix(Bn, QL, QR, nv)

    Dp = 0.5*(FL + FR) - 0.5*smax*(QR - QL)

    temp = matmul(Bn, (QR - QL) )

    Dm = Dp - 0.5*temp
    Dp = Dp + 0.5*temp

end subroutine RusanovFlux
