! pde.f90
! author: sunder

!-----------------------------------------------------------------------
! Routines related to the given PDE
!-----------------------------------------------------------------------

subroutine PDECons2Prim(Q, V)
    use ader_dg, only : nVar
    implicit none
    ! Argument list
    real, intent(in)     :: Q(nVar)     ! vector of conserved quantities
    real, intent(out)    :: V(nVar)     ! primitive variables

    V(1) = Q(1)

end subroutine PDECons2Prim

!-----------------------------------------------------------------------
! Convert primitive vector to conserved vector
!-----------------------------------------------------------------------

subroutine PDEPrim2Cons(V, Q)
    use ader_dg, ONLY : nVar
    implicit none
    ! Argument list
    real, intent(in)      :: V(nVar)     ! primitive variables
    real, intent(out)     :: Q(nVar)     ! vector of conserved quantities

    Q(1)   = V(1)

END subroutine PDEPrim2Cons

!-----------------------------------------------------------------------
! PDE flux tensor F(Q)
!-----------------------------------------------------------------------

subroutine PDEFlux(F, Q)
    use ader_dg, only : nVar, nDim, EQN
    implicit none
    ! Argument list
    real, intent(in)  :: Q(nVar)
    real, intent(out) :: F(nVar,nDim)

    F(1,1) = EQN%cx*Q(1)
    F(1,2) = EQN%cy*Q(1)

end subroutine PDEFlux

!-----------------------------------------------------------------------
! Viscous flux tensor F_v(Q)
!-----------------------------------------------------------------------

subroutine PDEViscFlux(F, Q, gradQ)
    use ader_dg, only : nVar, nDim, EQN
    implicit none
    ! Argument list
    real, intent(in) :: Q(nVar), gradQ(nVar,nDim)
    real, intent(out) :: F(nVar, nDim)

    F(1,1) = -EQN%a*gradQ(1,1)
    F(1,2) = -EQN%a*gradQ(1,2)

end subroutine PDEViscFlux

!-----------------------------------------------------------------------
! Eigenvalues of the PDE in the normal direction nv
!-----------------------------------------------------------------------

subroutine PDEEigenvalues(Lambda, Q, nv)
    use ader_dg, ONLY : nVar, nDim, EQN
    implicit none
    ! Argument list
    real, intent(in)  :: Q(nVar), nv(nDim)
    real, intent(out) :: Lambda(nVar)


    Lambda = (/ EQN%cx*nv(1) + EQN%cy*nv(2) /)   ! The eigenvalues of the advection-diffusion equations

end subroutine PDEEigenvalues

!-----------------------------------------------------------------------
! Eigenvalues (Viscous Part) of the PDE in the normal direction nv
!-----------------------------------------------------------------------

subroutine PDEViscEigenvalues(Lambda, Q)
    use ader_dg, ONLY : nVar, EQN
    implicit none
    ! Argument list
    real, intent(in)  :: Q(nVar)
    real, intent(out) :: Lambda(nVar)

    Lambda(1) = EQN%a

end subroutine PDEViscEigenvalues

!-----------------------------------------------------------------------
! Combine viscous and convective fluxes
!-----------------------------------------------------------------------

subroutine PDEFusedViscConvFlux(F, Q, gradQ)
    use ader_dg, only : nVar, nDim
    implicit none
    ! Argument list
    real, intent(in)  :: Q(nVar), gradQ(nVar,nDIm)
    real, intent(out) :: F(nVar,nDim)
    ! Local variables
    real :: Fc(nVar, nDim), Fv(nVar, nDim)

    call PDEFlux(Fc, Q)
    call PDEViscFlux(Fv, Q, gradQ)

    F = Fc + Fv

end subroutine PDEFusedViscConvFlux

!-----------------------------------------------------------------------
! Viscous LLF (Rusanov) Riemann solver
!-----------------------------------------------------------------------

subroutine  PDEViscLLFRiemannSolver(flux, QL, gradQL, QR, gradQR, nv)
    use ader_dg
    implicit none
    ! Argument list
    real, intent(in)  :: QL(nVar), gradQL(nVar,nDim), gradQR(nDim), QR(nVar), nv(nDim)
    real, intent(out) :: flux(nVar)
    
    ! Local variables

    real :: smax_c, smax_v, smax, factor, eta, LL(nVar), LR(nVar), LLv(nVar), LRv(nVar)
    real :: FL(nVar,nDim), FR(nVar,nDim)

    call PDEEigenvalues(LL,QL,nv)
    call PDEEigenvalues(LR,QR,nv)

    call PDEViscEigenvalues(LLv,QL)
    call PDEViscEigenvalues(LRv,QR)

    factor = real(2*N + 1)
    eta = factor/(dx(1)*sqrt(0.5*m_pi))

    smax_c = max( maxval(abs(LL)), maxval(abs(LR)) )
    smax_v = max( maxval(LLv), maxval(LRv) )
    smax   = smax_c + 2.0*eta*smax_v


    call PDEFusedViscConvFlux(FL,QL,gradQL)
    call PDEFusedViscConvFlux(FR,QR,gradQR)

    flux = 0.5*matmul(FL+FR, nv) - 0.5*smax*(QR-QL)

end subroutine  PDEViscLLFRiemannSolver


