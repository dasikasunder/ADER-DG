! initialize.f90
! author: sunder

!-----------------------------------------------------------------------
! Initialize the ADER-DG data structures
!-----------------------------------------------------------------------

subroutine initialize
    use ader_dg
    implicit none

    ! Local variables
    integer :: ii, jj, i, j, k, l, iGP, iVar, VMAX(nDim)
    real :: phi(N+1), phi0(N+1), phi1(N+1), phi_xi(N+1), corner(nDim), u0(nVar), xGP(nDim), x0(nDim)
    real :: dxi, xi, xi1, xi2
    real, pointer  :: LSQM(:,:), iLSQM(:,:), LSQrhs(:,:)
    logical  :: dmpresult

    ! -> Setup important parameters controlling the simulation (to be read by input file)

    open(unit = 1, file = 'params.ini', status = 'UNKNOWN')

    read(1,*)EQN%GAMMA
    read(1,*)IMAX
    read(1,*)JMAX
    read(1,*)xL(1)
    read(1,*)xL(2)
    read(1,*)xR(1)
    read(1,*)xR(2)
    read(1,*)bL
    read(1,*)bR
    read(1,*)bB
    read(1,*)bT
    read(1,*)ICType
    read(1,*)tend
    read(1,*)WriteInterval

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


    if (ICType .eq. 1) then
        BaseFile = 'Sod-Shock-Tube'
    else if (ICType .eq. 2) then
        BaseFile = 'Lax-Shock-Tube'
    else if (ICType .eq. 3) then
        BaseFile = 'ShuOsher-Shock-Tube'
    else if (ICType .eq. 4) then
        BaseFile = 'ShuOsher-Vortex'
    else if (ICType .eq. 5) then
        BaseFile = 'RP2D1'
    else if (ICType .eq. 6) then
        BaseFile = 'RP2D2'
    else if (ICType .eq. 7) then
        BaseFile = 'RP2D3'
    else if (ICType .eq. 8) then
        BaseFile = 'RP2D4'
    else if (ICType .eq. 9) then
        BaseFile = 'RP2D5'
    else if (ICType .eq. 10) then
        BaseFile = 'ShockVortex'
    else if (ICType .eq. 11) then
        BaseFile = 'DMR'
    else
        BaseFile = 'Sol-'
    end if

    ! -> Initialize stuff related to dofs and basis functions

    do i = 0, nDim ! i = 0 is the time dimension
        nDOF(i) = N+1
    end do

    do i = 1, nDim
        dn(i) = 1
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
                MM(k,l) = MM(k,l) + wGPN(iGP)*phi(k)*phi(l)       ! Mass-matrix
                Kxi(k,l) = Kxi(k,l) + wGPN(iGP)*phi_xi(k)*phi(l)  ! Stiffness-matrix
            end do
        end do
    end do

    call MatrixInverse(N+1,MM,iMM)

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

    call MatrixInverse(N+1,K1,iK1)

    FLcoeff = phi0  ! coefficients needed to extrapolate data onto the left  boundary
    FRcoeff = phi1  ! coefficients needed to extrapolate data onto the right boundary

    ! -> Main variables of the ADER-DG scheme

    allocate ( x(nDim, IMAX, JMAX) )
    allocate ( uh(nVar, nDOF(1), nDOF(2), IMAX , JMAX) )
    allocate ( wh(nVar, IMAX , JMAX) )
    allocate ( olduh(nVar, nDOF(1), nDOF(2), IMAX, JMAX) )
    allocate ( qhi(nVar, nDOF(1), nDOF(2), IMAX, JMAX) )
    allocate ( Fhi(nVar, nDim, nDOF(1), nDOF(2), IMAX, JMAX) )
    allocate ( duh(nVar, nDOF(1), nDOF(2), IMAX, JMAX) )
    allocate ( qbnd(nVar, nFac, nDOF(2), IMAX, JMAX) )
    allocate ( Fbnd(nVar, 4, nDOF(2), IMAX, JMAX) )
    allocate ( Fu(nVar, nDOF(2), IMAX+1 , JMAX ) )
    allocate ( Gu(nVar, nDOF(1), IMAX , JMAX+1 ) )

    ! -> Find the cell center coordinates

    do j = 1, JMAX
        do i = 1, IMAX
            x(1,i,j) = xL(1) + (real(i-1) + 0.5)*dx(1)
            x(2,i,j) = xL(2) + (real(j-1) + 0.5)*dx(2)
        end do
    end do

    ! -> Prepare the subcell limiter

    allocate( Limiter(IMAX, JMAX), recompute(IMAX, JMAX) )

    do j = 1, JMAX
        do i = 1, IMAX
            Limiter(i,j)%oldstatus = 0
            Limiter(i,j)%status = 0
        end do
    end do

    do i = 1, nDim
        nSubLimV(i) = nSubLim
    end do

    allocate( LSQM(N+2,N+2), iLSQM(N+2,N+2), LSQrhs(N+2,nSubLim) )
    allocate( uh2lim(nSubLim,N+1), lim2uh(N+1,nSubLim)           )
    allocate( uh2lob(N+1,N+1)                                    )

    uh2lim = 0.0
    dxi = 1.0/real(nSubLim)    ! the limiter uses nSubLim subintervals in each cell

    do i = 1, nSubLim
        xi1 = real(i-1)*dxi     ! left sub-interval boundary
        xi2 = real(i  )*dxi     ! right sub-interval boundary
        do ii = 1, N+1
            xi = xi1 + dxi*xiGPN(ii)
            call basis_func_1d(phi,phi_xi,xi)
            uh2lim(i,:) = uh2lim(i,:) + wGPN(ii)*phi(:)
        end do
    end do

    ! To go back from the limited (piecewise constant sub-cell solution) to the uh solution
    ! we use a constrained least-squares reconstruction, which preserves the average exactly

    LSQM(1:N+1,1:N+1) = 2*matmul( transpose(uh2lim), uh2lim )

    LSQM(N+2,1:N+1)   =  wGPN(:)
    LSQM(1:N+1,N+2)   = -wGPN(:)

    LSQM(N+2,N+2) = 0.0
    call MatrixInverse(N+2,LSQM,iLSQM)

    LSQrhs(1:N+1,:)   = 2.0*transpose(uh2lim)
    LSQrhs(N+2,:)     = dxi
    lim2uh = matmul( iLSQM(1:N+1,:), LSQrhs )

    ! Compute the Gauss-Lobatto quadrature nodes
    call gaulob(0.0, 1.0, xiLob, wLob, N+1)

    do ii = 1, N+1
        call basis_func_1d(phi,phi_xi,xiLob(ii))
        uh2lob(ii,:) = phi(:)
    end do

    do i = 1, nSubLim
         xilimbary(i) = 0.0 + 0.5/real(nSubLim) + real(i-1)*1.0/real(nSubLim)
    end do

    deallocate( LSQM, iLSQM, LSQrhs )

    ! Compute the mapping from face to 2D neighbor index

    Face2Neigh(:,1) = (/ -1, 0 /)
    Face2Neigh(:,2) = (/ +1, 0 /)
    Face2Neigh(:,3) = (/  0,-1 /)
    Face2Neigh(:,4) = (/  0,+1 /)

    ! Find the Voronoi neighbors of a cell

    allocate(neighbor(-1:1,-1:1, nDim, IMAX, JMAX) )

    do j = 1, JMAX
        do i = 1, IMAX

            neighbor(-1,:,1,i,j) = MAX(i-1, 1)
            neighbor( 0,:,1,i,j) = i
            neighbor( 1,:,1,i,j) = MIN(i+1, IMAX)

            neighbor(:,-1,2,i,j) = MAX(j-1, 1)
            neighbor(:, 0,2,i,j) = j
            neighbor(:, 1,2,i,j) = MIN(j+1, JMAX)

            ! Add periodicity to the neighbours

            if (i .eq. 1 .and. bL .eq. 3) then
                neighbor(-1,:,1,i,j) = IMAX
            end if

            if (i .eq. IMAX .and. bR .eq. 3) then
                neighbor(1,:,1,i,j) = 1
            end if

            if (j .eq. 1 .and. bB .eq. 3) then
                neighbor(:,-1,2,i,j) = JMAX
            end if

            if (j .eq. JMAX .and. bT .eq. 3) then
                neighbor(:,1,2,i,j) = 1
            end if

        end do
    end do

    ! -> Initialize the solution

    do j = 1, JMAX
        do i = 1, IMAX

            corner(:) = x(:,i,j) - 0.5*dx(:) ! coordinate of the lower left corner of the cell
            do l = 1, nDOF(2)
                do k = 1, nDOF(1)
                    xGP = corner + (/ xiGPN(k), xiGPN(l) /)*dx(:)
                    call InitialField(xGP, u0)
                    uh(:, k, l, i, j) = u0(:)
                end do
            end do

            if (N .gt. 0) then

                do iVar = 1, nVar
                    Limiter(i,j)%lmin(iVar) = minval(uh(iVar,:,:,i,j))
                    Limiter(i,j)%lmax(iVar) = maxval(uh(iVar,:,:,i,j))
                end do

                call DMP(dmpresult,uh(:,:,:,i,j),Limiter(i,j),1.0e-1)

                if( .not. dmpresult) then

                    ! If initial condition projection does not satisfy the DMP, then activate the subcell limiter

                    if (Limiter(i,j)%oldstatus .eq. 0) then
                        allocate( Limiter(i,j)%Lh(nVar,nSubLimV(1),nSubLimV(2))    )
                        allocate( Limiter(i,j)%NewLh(nVar,nSubLimV(1),nSubLimV(2)) )
                    end if

                    Limiter(i,j)%status    = 1
                    Limiter(i,j)%oldstatus = 1

                    do jj = 1, nSubLimV(2)
                        do ii = 1, nSubLimV(1)
                            xGP = x0 + (/ xilimbary(ii), xilimbary(jj) /)*dx(:)
                            call InitialField(xGP,u0)
                            Limiter(i,j)%Lh(:,ii,jj) = u0
                            Limiter(i,j)%NewLh(:,ii,jj) = u0
                        end do
                    end do
                end if
            end if

        end do
    end do

end subroutine initialize

subroutine InitialField(xGP, u0)
    use ader_dg
    implicit none
    ! Argument list
    real, intent(in ) :: xGP(nDim)          ! spatial position vector
    real, intent(out) :: u0(nVar)           ! initial data vector in terms of conserved variables
    ! Local variables
    real :: v0(nVar), xx, yy
    real, parameter :: m_pi = 4.0*atan(1.0)
    
    xx = xGP(1)
    yy = xGP(2)

    if (xx .le. 1.0/6.0 + yy/tan(m_pi/3.0)) then

        v0(1) = 8.0
        v0(2) = 8.25*cos(m_pi/6.0)
        v0(3) = -8.25*sin(m_pi/6.0)
        v0(4) = 116.5
    else

        v0(1) = 1.4;
        v0(2) = 0.0;
        v0(3) = 0.0;
        v0(4) = 1.0;

    end if

    call PDEPrim2Cons(V0, u0)

end subroutine InitialField
