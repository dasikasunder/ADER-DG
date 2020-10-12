! initialize.f90
! author: sunder

!-----------------------------------------------------------------------
! Initialize the ADER-DG data structures
!-----------------------------------------------------------------------

subroutine initialize
    use ader_dg
    implicit none

    ! Local variables
    integer :: i, j, k, l, iGP, VMAX(nDim)
    real :: phi(N+1), phi0(N+1), phi1(N+1), phi_xi(N+1), corner(nDim), u0(nVar), xGP(nDim)


    ! -> Setup important parameters controlling the simulation (to be read by input file)

    open(UNIT = 1, FILE = 'params.ini', STATUS = 'UNKNOWN')

    read(1,*)EQN%cx
    read(1,*)EQN%cy
    read(1,*)EQN%a
    read(1,*)IMAX
    read(1,*)JMAX
    read(1,*)xL(1)
    read(1,*)xL(2)
    read(1,*)xR(1)
    read(1,*)xR(2)
    read(1,*)tend
    read(1,*)WriteInterval

    xL(1) = 0.0
    xL(2) = 0.0

    xR(1) = 0.25*m_pi
    xR(2) = 2.0*m_pi

    close(1)

    ! -> Setup mesh and time related data

    VMAX = (/ IMAX, JMAX /)
    dx = (xR-xL)/VMAX
    timestep = 0
    dt = 0.0
    time = 0.0
    nElem = IMAX*JMAX
    nFace = JMAX*IMAX + JMAX*IMAX
    nNode = (IMAX+1)*(JMAX+1)

    BaseFile = 'Sol'

    ! -> Initialize stuff related to dofs and basis functions

    do i = 0, nDim ! i = 0 is the time dimension
        nDOF(i) = N+1
    end do

    call gauleg(0.0, 1.0, xiGPN, wGPN, N+1) ! Gauss legendre points and weights
    xin = xiGPN                                 ! Basis functions run through Gauss-Legendre points
    MM   = 0.0                                ! Element mass matrix
    Kxi  = 0.0                                ! Element stiffness matrix
    dudx = 0.0                                ! discrete derivative operator, which projects the derivatives onto the basis

    do iGP = 1, N+1
        call basis_func_1d(phi,phi_xi,xiGPN(iGP))
        do k = 1, N+1
            do l = 1, N+1
                MM(k,l)  = MM(k,l) + wGPN(iGP)*phi(k)*phi(l)      ! Mass-matrix
                Kxi(k,l) = Kxi(k,l) + wGPN(iGP)*phi_xi(k)*phi(l)  ! Stiffness-matrix
            end do
        end do
    end do

    call matrix_inverse(N+1,MM,iMM)

    dudx = matmul( iMM, transpose(Kxi) )

    call basis_func_1d(phi0, phi_xi, 0.0)   ! Compute the basis functions on the left
    call basis_func_1d(phi1, phi_xi, 1.0)   ! Compute the basis function on the right

    ! The flux matrices are all possible combinations of left and right

    do k = 1, N+1
        do l = 1, N+1
            FLm(k,l) = phi0(k)*phi1(l)   ! Left contribution to the left flux matrix    (m = left  of the interface)
            FLp(k,l) = phi0(k)*phi0(l)   ! Right contribution to the left flux matrix   (p = right of the interface)
            FRm(k,l) = phi1(k)*phi1(l)   ! Left contribution to the right flux matrix   (m = left  of the interface)
            FRp(k,l) = phi1(k)*phi0(l)   ! Right contribution to the right flux matrix  (p = right of the interface)
        end do
    end do

    ! The time flux matrices for the ADER-DG predictor method are given by the principle of upwinding in time (causality principle)

    F0 = phi0   ! upwinding in time = information comes from smaller times
    F1 = FRm    ! upwinding in time = information comes from smaller times
    K1 = F1 - Kxi

    call matrix_inverse(N+1,K1,iK1)

    FLcoeff = phi0  ! coefficients needed to extrapolate data onto the left  boundary
    FRcoeff = phi1  ! coefficients needed to extrapolate data onto the right boundary

    ! -> Main variables of the ADER-DG scheme

    allocate ( x(nDim, IMAX, JMAX) )
    allocate ( uh(nVar, nDOF(1), nDOF(2), IMAX , JMAX) )
    allocate ( wh(nVar, IMAX , JMAX) )
    allocate ( qhi(nVar, nDOF(1), nDOF(2), IMAX, JMAX) )
    allocate ( Fhi(nVar, nDim, nDOF(1), nDOF(2), IMAX, JMAX) )
    allocate ( duh(nVar, nDOF(1), nDOF(2), IMAX, JMAX) )
    allocate ( qbnd(nVar, nFac, nDOF(2), IMAX, JMAX) )
    allocate ( gradQbnd(nVar, nDim, nFac, nDOF(2), IMAX, JMAX) )
    allocate ( Fu(nVar, nDOF(2), IMAX+1 , JMAX ) )
    allocate ( Gu(nVar, nDOF(1), IMAX , JMAX+1 ) )

    ! -> Find the cell center coordinates

    do j = 1, JMAX
        do i = 1, IMAX
            x(1,i,j) = xL(1) + (real(i-1) + 0.5)*dx(1)
            x(2,i,j) = xL(2) + (real(j-1) + 0.5)*dx(2)
        end do
    end do

    ! -> Initialize the solution

    do j = 1, JMAX
        do i = 1, IMAX

            corner(:) = x(:,i,j) - 0.5*dx(:) ! coordinate of the lower left corner of the cell
            do l = 1, nDOF(2)
                do k = 1, nDOF(1)
                    xGP = corner + (/ xiGPN(k), xiGPN(l) /)*dx(:)
                    call initial_field(xGP, u0)
                    uh(:, k, l, i, j) = u0(:)
                end do
            end do

        end do
    end do

end subroutine initialize

!-----------------------------------------------------------------------
! Initial condition function
!-----------------------------------------------------------------------

subroutine initial_field(xGP, u0)
    use ader_dg
    implicit none
    ! Argument list
    real, intent(in ) :: xGP(nDim)          ! spatial position vector
    real, intent(out) :: u0(nVar)           ! initial data vector in terms of conserved variables
    ! Local variables
    real :: v0(nVar), xx, yy

    xx = xGP(1)
    yy = xGP(2)

    v0(1) = sin((yy))

    call PDEPrim2Cons(V0, u0)

end subroutine initial_field
