! limiter.f90
! author: sunder

!-----------------------------------------------------------------------
! Routines related to the subcell limiter
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Compute the subcell data lim from a given DG polynomial luh
!-----------------------------------------------------------------------

subroutine get_subcell_data(lim,luh)
    use ader_dg
    implicit none
    ! Argument list
    real :: lim(nVar,nSubLimV(1),nSubLimV(2)), luh(nVar,nDOF(1),nDOF(2))
    ! Local variables
    integer :: iii, jjj
    real :: limy(nVar,nDOF(1),nSubLimV(2))

     ! mapping of uh to the sub-cell limiter along y direction

     limy = 0.0

     do iii = 1, nDOF(1) ! x
         limy(:,iii,:) = matmul( luh(:,iii,:), TRANSPOSE(uh2lim) )
     end do

     ! mapping of uh to the sub-cell limiter along x direction

     lim = 0.0

     do jjj = 1, nSubLimV(2) ! y
         lim(:,:,jjj) = matmul( limy(:,:,jjj), TRANSPOSE(uh2lim) )
     end do

 end subroutine get_subcell_data

!-----------------------------------------------------------------------
! Compute the point values lob in the Gauss-Lobatto points from a given DG polynomial luh
!-----------------------------------------------------------------------

subroutine get_lobatto_data(lob,luh)
    use ader_dg
    implicit none
    ! Argument list
    real :: lob(nVar,nDOF(1),nDOF(2)), luh(nVar,nDOF(1),nDOF(2))
    ! Local variables
    integer :: iii, jjj
    real :: loby(nVar,nDOF(1),nDOF(2))

    ! mapping of uh to the Gauss-Lobatto points along y direction

    loby = 0.0

    do iii = 1, nDOF(1) ! x
        loby(:,iii,:) = matmul( luh(:,iii,:), TRANSPOSE(uh2lob) )
    end do

    ! mapping of uh to the Gauss-Lobatto points along x direction

    lob = 0.0

    do jjj = 1, nDOF(2) ! y
        lob(:,:,jjj) = matmul( loby(:,:,jjj), TRANSPOSE(uh2lob) )
    end do

end subroutine get_lobatto_data

!-----------------------------------------------------------------------
! Put subcell data back into the solution
!-----------------------------------------------------------------------

subroutine put_subcell_data(luh,lim)
    use ader_dg
    implicit none
    ! Argument list
    real :: lim(nVar,nSubLimV(1),nSubLimV(2)), luh(nVar,nDOF(1),nDOF(2))
    ! Local variables
    integer :: iii, jjj
    real :: luhy(nVar,nSubLimV(1),N+1)

    ! Mapping of the sub-cell limiter to uh along y direction
    luhy = 0.0

    do iii = 1, nSubLimV(1) ! x
        luhy(:,iii,:) = matmul( lim(:,iii,:), TRANSPOSE(lim2uh) )
    end do

    ! Mapping of the sub-cell limiter to uh along x direction

   luh = 0.0

    do jjj = 1, nDOF(2)  ! y
        luh(:,:,jjj) = matmul( luhy(:,:,jjj), TRANSPOSE(lim2uh) )
    end do

end subroutine put_subcell_data

!-----------------------------------------------------------------------
! Get the minimum and maximum in the neighbourhood of a cell
!-----------------------------------------------------------------------

subroutine get_min_max
    use ader_dg
    implicit none
    ! Local variables
    integer :: i, j, ii, jj, in, jn,  iVar
    real :: lim(nVar,nSubLimV(1),nSubLimV(2))
    real :: lob(nVar,nDOF(1),nDOF(2))
    real :: Qmin(nVar), Qmax(nVar)

    ! Step 1) Loop over all elements and get the min/max for each variable

    do j = 1, JMAX
        do i = 1, IMAX

            if(Limiter(i,j)%status .eq. 0) then

                call get_subcell_data( lim, uh(:,:,:,i,j) )
                call get_lobatto_data( lob, uh(:,:,:,i,j) )

                do iVar = 1, nVar

                    Limiter(i,j)%lmin_old(iVar) = min( minval(uh(iVar,1:nDOF(1),1:nDOF(2),i,j)), &
                                                   minval(lob(iVar,1:nDOF(1),1:nDOF(2))), &
                                                   minval(lim(iVar,1:nSubLimV(1),1:nSubLimV(2))) )

                    Limiter(i,j)%lmax_old(iVar) = max( maxval(uh(iVar,1:nDOF(1),1:nDOF(2),i,j)), &
                                                   maxval(lob(iVar,1:nDOF(1),1:nDOF(2))), &
                                                   maxval(lim(iVar,1:nSubLimV(1),1:nSubLimV(2))) )

                end do
            else
                lim = Limiter(i,j)%Lh
                do iVar = 1, nVar
                    Limiter(i,j)%lmin_old(iVar) = minval( lim(iVar,1:nSubLimV(1),1:nSubLimV(2)) )
                    Limiter(i,j)%lmax_old(iVar) = maxval( lim(iVar,1:nSubLimV(1),1:nSubLimV(2)) )
                end do
            end if

        end do
    end do

    ! Step 2) Loop over all Voronoi neighbors and get the final min/max for each variable (requires step 1 to be complete)

    do j = 1, JMAX
        do i = 1, IMAX
            Qmin = Limiter(i,j)%lmin_old(:)
            Qmax = Limiter(i,j)%lmax_old(:)

            do jj = -1, 1
                do ii = -1, 1
                    in = neighbor(ii,jj,1,i,j)
                    jn = neighbor(ii,jj,2,i,j)
                    do iVar = 1, nVar
                        Qmin(iVar) = MIN( Qmin(iVar), Limiter(in,jn)%lmin_old(iVar) )
                        Qmax(iVar) = MAX( Qmax(iVar), Limiter(in,jn)%lmax_old(iVar) )
                    end do
                end do
            end do

            Limiter(i,j)%lmin(:) = Qmin
            Limiter(i,j)%lmax(:) = Qmax

        end do
    end do

    continue

end subroutine get_min_max


!-----------------------------------------------------------------------
! Check discrete maximum principle (DMP) and physical admissibility
! in the given cell
!-----------------------------------------------------------------------

subroutine DMP(dmpresult,arguh,argLimiter,argmax)
    use ader_dg
    implicit none
    ! Argument list
    logical        :: dmpresult
    real           :: arguh(nVar,nDOF(1),nDOF(2)),argmax
    TYPE(tLimiter) :: argLimiter
    ! Local variables
    integer :: iii, jjj, iVar, iErr
    real    :: lim(nVar,nSubLimV(1),nSubLimV(2))
    real    :: lob(nVar,nDOF(1),nDOF(2))
    real    :: lmin(nVar), lmax(nVar), ldiff, Qout(nVar)


    dmpresult = .true.

    call get_subcell_data( lim, arguh )
    call get_lobatto_data( lob, arguh )

    do iVar = 1, nVar

        if (iVar .gt. 0) then

            lmin(iVar) = min( minval(arguh(iVar,1:nDOF(1),1:nDOF(2))), &
                            minval(lim(iVar,1:nSubLimV(1),1:nSubLimV(2))), &
                            minval(lob(iVar,1:nDOF(1),1:nDOF(2))) )

            lmax(iVar) = max( maxval(arguh(iVar,1:nDOF(1),1:nDOF(2))), &
                            maxval(lim(iVar,1:nSubLimV(1),1:nSubLimV(2))), &
                            maxval(lob(iVar,1:nDOF(1),1:nDOF(2))) )

            ldiff = max( 1.0e-4, 1.0e-3*(argLimiter%lmax(iVar) - argLimiter%lmin(iVar)), argmax )

            if (lmin(iVar) .lt. argLimiter%lmin(iVar)-ldiff) then
                dmpresult = .false.
                return
            end if

            if (lmax(iVar) .gt. argLimiter%lmax(iVar)+ldiff) then
                dmpresult = .false.
                return
            end if

        end if
    end do

    ! Physical detection criteria

    do jjj = 1, nSubLimV(2)
        do iii = 1, nSubLimV(1)
            call PDEAssurePositivity(iErr, lim(:,iii,jjj), Qout )
            if (iErr .lt. 0) then
                dmpresult = .false.
                return
            end if
            !do iVar = 1, nVar TODO: Add this condition for NaN later
                !if (isnan (lim(iVar,iii,jjj)) ) then
                    !dmpresult = .false.
                    !return
                !end if
            !end do
        end do
    end do

   continue

end subroutine DMP

!-----------------------------------------------------------------------
! Save uh before the limiter is applied
!-----------------------------------------------------------------------

subroutine save_olduh
    use ader_dg
    implicit none
    ! Local variables
    integer  :: i,j

    do j = 1, JMAX
        do i = 1, IMAX
            olduh(:,:,:,i,j) = uh(:,:,:,i,j)
        end do
    end do

    continue

end subroutine save_olduh

!-----------------------------------------------------------------------
! After subcell recomputation, update status of the limiter deallocate 
! uneccesary memory
!-----------------------------------------------------------------------

subroutine update_limiter
    use ader_dg
    implicit none
    ! Local variables
    integer :: i,j

    do j = 1, JMAX
        do i = 1, IMAX

            if (Limiter(i,j)%status .eq. 0) then
                if(Limiter(i,j)%oldstatus .ne. 0) then
                    deallocate( Limiter(i,j)%Lh    )
                    deallocate( Limiter(i,j)%NewLh )
                end if
            else
                Limiter(i,j)%Lh = Limiter(i,j)%NewLh
            end if

            Limiter(i,j)%oldstatus = Limiter(i,j)%status

            continue
        end do
    end do

end subroutine update_limiter

!-----------------------------------------------------------------------
! Allocate memory for the limiter
!-----------------------------------------------------------------------

subroutine allocate_limiter
    use ader_dg
    implicit none
    ! Local variables
    integer    :: i,j

    do j = 1, JMAX
        do i = 1, IMAX
            if (recompute(i,j) .eq. 1) then
                if(Limiter(i,j)%oldstatus .eq. 0) then
                    allocate( Limiter(i,j)%Lh(nVar,nSubLimV(1),nSubLimV(2))    )
                    allocate( Limiter(i,j)%NewLh(nVar,nSubLimV(1),nSubLimV(2)) )
                end if
                Limiter(i,j)%status = 1
            else
                Limiter(i,j)%status = 0
            end if
        end do
    end do

end subroutine allocate_limiter

!-----------------------------------------------------------------------
! Once a cell is marked for recomputation, also mark the face
! neighbours for recomputation
!-----------------------------------------------------------------------

subroutine spread_recompute
    use ader_dg
    implicit none
    ! Local variables
    integer  :: i, j, ii, jj, in, jn
    integer  :: iFace

    do j = 1, JMAX
        do i = 1, IMAX
            if (recompute(i,j) .eq. 1) then
                do iFace = 1, 2*nDim
                    ii = Face2Neigh(1,iFace)
                    jj = Face2Neigh(2,iFace)
                    in = neighbor(ii,jj,1,i,j)
                    jn = neighbor(ii,jj,2,i,j)
                    if (recompute(in,jn) .eq. 0) then
                        recompute(in,jn) = 2
                    end if
                end do
            end if
        end do
    end do

    continue

end subroutine spread_recompute

!-----------------------------------------------------------------------
! Rusanov flux (For the subcell limiter. This subroutine does not take 
! fluxes as input, instead computes them from the input state
!-----------------------------------------------------------------------

 subroutine RusanovFluxS(flux,QL,QR,nv)
    use ader_dg
    implicit none
    ! Argument list
    real, intent(in)  :: QL(nVar), QR(nVar), nv(nDim)
    real, intent(out) :: flux(nVar)
    ! Local variables
    real :: smax, LL(nVar), LR(nVar)
    real :: FL(nVar,nDim), FR(nVar,nDim)
    
    call PDEEigenvalues(QL,nv,LL)
    call PDEEigenvalues(QR,nv,LR)
    smax = max( maxval(abs(LL)), maxval(abs(LR)) )
    call PDEFlux(QL,FL)
    call PDEFlux(QR,FR)
    flux = 0.5*matmul(FL+FR,nv) - 0.5*smax*(QR-QL)
    
end subroutine

!-----------------------------------------------------------------------
! Min-mod limiter
!-----------------------------------------------------------------------

subroutine minmod(c,a,b,nnn)
    implicit none
    ! Argument list
    integer, intent(in)   :: nnn
    real, intent(in)      :: a(nnn), b(nnn)
    real, intent(out)     :: c(nnn)
    ! Local variables
    integer :: i

    do i = 1, nnn
        if (a(i)*b(i) .le. 0.0) then
            c(i) = 0.0
        else
            if( abs(a(i)) .lt. abs(b(i)) ) then
                c(i) = a(i)
            else
                c(i) = b(i)
            end if
        end if
    end do

end subroutine minmod

!-----------------------------------------------------------------------
! Recompute solution in the subcell
!-----------------------------------------------------------------------

subroutine subcell_recompute
    use ader_dg
    implicit none
    ! Local variables
    integer :: i, j, in, jn, ii, jj
    real    :: ldx(nDim), nv(nDim), nvi(nDim,nDim)
    real    :: subuh0(nVar,(1-nSubLimV(1)):2*nSubLimV(1),(1-nSubLimV(2)):2*nSubLimV(2))
    real    :: subuh1(nVar,nSubLimV(1),nSubLimV(2))
    real    :: slopex(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2))
    real    :: slopey(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2))
    real    :: wLx(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2))
    real    :: wLy(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2))
    real    :: wRx(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2))
    real    :: wRy(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2))
    real    :: lim(nVar,nSubLimV(1),nSubLimV(2))
    real    :: Qt(nVar)
    real    :: FLx(nVar,nDim), FRx(nVar,nDim), FLy(nVar,nDim), FRy(nVar,nDim)
    real    :: Fx(nVar,1:nSubLimV(1)+1,1:nSubLimV(2))
    real    :: Fy(nVar,1:nSubLimV(1),1:nSubLimV(2)+1)
    real    :: Va(nVar), Vb(nVar), Qbc(nVar), x_sw
    real, parameter :: m_pi = 4.0*atan(1.0)
    logical :: dmpresult

    Va(1) = 8.0
    Va(2) = 8.25*cos(m_pi/6.0)
    Va(3) = -8.25*sin(m_pi/6.0)
    Va(4) = 116.5

    Vb(1) = 1.4;
    Vb(2) = 0.0;
    Vb(3) = 0.0;
    Vb(4) = 1.0;

    call shock_position(time, x_sw)

    nvi = 0.0

    do ii = 1, nDim
        nvi(ii,ii) = 1.0
    end do

    ldx = dx/real(nSubLim)
    Fx = 0.0
    Fy = 0.0

    slopex = 0.0
    slopey = 0.0

    wLx = 0.0
    wRx = 0.0
    wLy = 0.0
    wRy = 0.0

    do j = 1, JMAX
        do i = 1, IMAX
            if (recompute(i,j) > 0) then
                ! Get the initial data, either from a previous limiter stage, or from healthy uh initial data
                ! note that here we only need the values of the face neighbors. For simplicity, we run over all the Voronoi neighbors

                do jj = -dn(2), dn(2)
                    do ii = -dn(1), dn(1)

                        in = neighbor(ii,jj,1,i,j)
                        jn = neighbor(ii,jj,2,i,j)

                        if (Limiter(in,jn)%oldstatus .eq. 0) then
                            call get_subcell_data(lim, olduh(:,:,:,in,jn))
                        else
                            lim = Limiter(in,jn)%Lh
                        end if

                        subuh0(:,1+ii*nSubLimV(1):(ii+1)*nSublimV(1),1+jj*nSubLimV(2):(jj+1)*nSublimV(2)) = lim(:,:,:)

                    end do
                end do

                !----------------------------------------------------
                ! Apply boundary conditions (reflective and inflow)
                !----------------------------------------------------

                if (i .eq. 1) then ! Left

                    call PDEPrim2Cons(Va, Qbc)

                    do jj = (1-nSubLimV(2)),2*nSubLimV(2)
                        do ii = 1-nSubLimV(1),0
                        subuh0(:, ii, jj) = Qbc
                        end do
                    end do
                end if

                if (j .eq. 1) then ! Bottom
                    if (real(i)*dx(1) .gt. 1.0/6.0) then ! Reflective boundary
                        subuh0(3,:,1-nSubLimV(2):0) = -subuh0(3,:,1-nSubLimV(2):0)
                    end if
                end if

                if (j .eq. JMAX) then ! Top

                    if (real(i)*dx(1) .le. x_sw) then
                        call PDEPrim2Cons(Va, Qbc)
                    else
                        call PDEPrim2Cons(Vb, Qbc)
                    end if

                    do jj = 1+nSubLimV(2),2*nSublimV(2)
                        do ii = 1-nSubLimV(2),2*nSubLimV(2)
                            subuh0(:,ii,jj) = Qbc
                        end do
                    end do

                end if


                ! --------------------------------------------------------------
                ! Second order TVD MUSCL reconstruction in space and time
                ! --------------------------------------------------------------

                do jj = 1, nSubLimV(2)
                    do ii = 1, nSubLimV(1)
                        ! x slopes
                        call minmod( slopex(:,ii,jj), subuh0(:,ii+1,jj)-subuh0(:,ii,jj), subuh0(:,ii,jj)-subuh0(:,ii-1,jj), nVar )
                        ! y slopes
                        call minmod( slopey(:,ii,jj), subuh0(:,ii,jj+1)-subuh0(:,ii,jj), subuh0(:,ii,jj)-subuh0(:,ii,jj-1), nVar )
                    end do
                end do

                ! Time evolution and boundary extrapolated data

                do jj = 1-dn(2), nSubLimV(2)+dn(2)
                    do ii = 1-dn(1), nSubLimV(1)+dn(1)

                        ! In the x-direction

                        wLx(:,ii,jj) = subuh0(:,ii,jj) - 0.5*slopex(:,ii,jj)
                        wRx(:,ii,jj) = subuh0(:,ii,jj) + 0.5*slopex(:,ii,jj)

                        call PDEFlux(wLx(:,ii,jj),FLx)
                        call PDEFlux(wRx(:,ii,jj),FRx)

                        Qt = - (FRx(:,1)-FLx(:,1))/ldx(1)

                        ! In the y-direction

                        wLy(:,ii,jj) = subuh0(:,ii,jj) - 0.5*slopey(:,ii,jj)
                        wRy(:,ii,jj) = subuh0(:,ii,jj) + 0.5*slopey(:,ii,jj)

                        call PDEFlux(wLy(:,ii,jj), FLy)
                        call PDEFlux(wRy(:,ii,jj), FRy)

                        Qt = Qt - (FRy(:,2)-FLy(:,2))/ldx(2)

                        wLx(:,ii,jj) = wLx(:,ii,jj) + 0.5*dt*Qt
                        wRx(:,ii,jj) = wRx(:,ii,jj) + 0.5*dt*Qt
                        wLy(:,ii,jj) = wLy(:,ii,jj) + 0.5*dt*Qt
                        wRy(:,ii,jj) = wRy(:,ii,jj) + 0.5*dt*Qt

                    end do
                end do

                nv = nvi(:,1)
                ! Solve the Riemann problems on the x-edges of the subgrid

                do jj = 1, nSubLimV(2)
                    do ii = 1, nSubLimV(1)+dn(1)
                        call RusanovFluxS(Fx(:,ii,jj),wRx(:,ii-1,jj),wLx(:,ii,jj), nv)
                    end do
                end do

                nv = nvi(:,2)

                ! Solve the Riemann problems on the y-edges of the subgrid

                do jj = 1, nSubLimV(2)+dn(2)
                    do ii = 1, nSubLimV(1)
                        call RusanovFluxS(Fy(:,ii,jj),wRy(:,ii,jj-1),wLy(:,ii,jj),nv)
                    end do
                end do

                ! Evolve the subcell data inside cell i using a simple first order finite volume scheme

                do jj = 1, nSubLimV(2)
                    do ii = 1, nSubLimV(1)
                        subuh1(:,ii,jj) = subuh0(:,ii,jj) - dt/ldx(1)*(Fx(:,ii+1,jj)-Fx(:,ii,jj)) - &
                                                            dt/ldx(2)*(Fy(:,ii,jj+1)-Fy(:,ii,jj))
                    end do
                end do

                if (recompute(i,j) .eq. 1) then ! If a cell is really troubled, then set the limiter status to 1 and set the new subcell data

                    Limiter(i,j)%status = 1

                    ! Save the evolved data in the limiter
                    do jj = 1, nSubLimV(2)
                        do ii = 1, nSubLimV(1)
                            Limiter(i,j)%NewLh(:,ii,jj) = subuh1(:,ii,jj)
                        end do
                    end do

                    ! Put the evolved cell-averaged data back into the DG polynomial
                    call put_subcell_data(uh(:,:,:,i,j),Limiter(i,j)%NewLh)

                else if (recompute(i,j) .eq. 2) then

                    ! Put the evolved cell-averaged data back into the DG polynomial
                    call put_subcell_data(uh(:,:,:,i,j),subuh1)

                    ! If we recompute neighbors of troubled cells, and the resulting reconstructed DG polynomial is troubled, then activate the limiter also in this case
                    call DMP(dmpresult,uh(:,:,:,i,j),Limiter(i,j),0.0)

                    if (dmpresult) then
                        Limiter(i,j)%status = 0
                    else
                        Limiter(i,j)%status = 1
                        if (Limiter(i,j)%oldstatus .eq. 0) then
                        allocate( Limiter(i,j)%Lh(nVar,nSubLimV(1),nSubLimV(2))    )
                        allocate( Limiter(i,j)%NewLh(nVar,nSubLimV(1),nSubLimV(2)) )
                        end if

                        ! Save the evolved data in the limiter

                        do jj = 1, nSubLimV(2)
                            do ii = 1, nSubLimV(1)
                                Limiter(i,j)%NewLh(:,ii,jj) = subuh1(:,ii,jj)
                            end do
                        end do

                    end if
                end if

            end if
        end do
    end do

    continue

end subroutine subcell_recompute
