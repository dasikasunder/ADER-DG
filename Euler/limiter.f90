! limiter.f90
! author: sunder

!-----------------------------------------------------------------------
! Routines related to the subcell limiter
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Compute the subcell data lim from a given DG polynomial luh
!-----------------------------------------------------------------------

subroutine GetSubcellData(lim,luh)
    use ader_dg
    implicit none
    ! Argument list
    double precision :: lim(nVar,nSubLimV(1),nSubLimV(2)), luh(nVar,nDOF(1),nDOF(2))
    ! Local variables
    integer :: iii, jjj
    double precision :: limy(nVar,nDOF(1),nSubLimV(2))

     ! mapping of uh to the sub-cell limiter along y direction

     limy = 0.0d0

     do iii = 1, nDOF(1) ! x
         limy(:,iii,:) = matmul( luh(:,iii,:), TRANSPOSE(uh2lim) )
     end do

     ! mapping of uh to the sub-cell limiter along x direction

     lim = 0.0d0

     do jjj = 1, nSubLimV(2) ! y
         lim(:,:,jjj) = matmul( limy(:,:,jjj), TRANSPOSE(uh2lim) )
     end do

 end subroutine GetSubcellData

!-----------------------------------------------------------------------
! Compute the point values lob in the Gauss-Lobatto points from a given DG polynomial luh
!-----------------------------------------------------------------------

subroutine GetLobattoData(lob,luh)
    use ader_dg
    implicit none
    ! Argument list
    double precision :: lob(nVar,nDOF(1),nDOF(2)), luh(nVar,nDOF(1),nDOF(2))
    ! Local variables
    integer :: iii, jjj
    double precision :: loby(nVar,nDOF(1),nDOF(2))

    ! mapping of uh to the Gauss-Lobatto points along y direction

    loby = 0.0d0

    do iii = 1, nDOF(1) ! x
        loby(:,iii,:) = matmul( luh(:,iii,:), TRANSPOSE(uh2lob) )
    end do

    ! mapping of uh to the Gauss-Lobatto points along x direction

    lob = 0.0d0

    do jjj = 1, nDOF(2) ! y
        lob(:,:,jjj) = matmul( loby(:,:,jjj), TRANSPOSE(uh2lob) )
    end do

end subroutine GetLobattoData

!-----------------------------------------------------------------------
! Put subcell data back into the solution
!-----------------------------------------------------------------------

subroutine PutSubcellData(luh,lim)
    use ader_dg
    implicit none
    ! Argument list
    double precision :: lim(nVar,nSubLimV(1),nSubLimV(2)), luh(nVar,nDOF(1),nDOF(2))
    ! Local variables
    integer :: iii, jjj
    double precision :: luhy(nVar,nSubLimV(1),N+1)

    ! Mapping of the sub-cell limiter to uh along y direction
    luhy = 0.0d0

    do iii = 1, nSubLimV(1) ! x
        luhy(:,iii,:) = matmul( lim(:,iii,:), TRANSPOSE(lim2uh) )
    end do

    ! Mapping of the sub-cell limiter to uh along x direction

   luh = 0.0d0

    do jjj = 1, nDOF(2)  ! y
        luh(:,:,jjj) = matmul( luhy(:,:,jjj), TRANSPOSE(lim2uh) )
    end do

end subroutine PutSubcellData

!-----------------------------------------------------------------------
! Get the minimum and maximum in the neighbourhood of a cell
!-----------------------------------------------------------------------

subroutine GetMinMax
    USE ader_dg
    implicit none
    ! Local variables
    integer :: i, j, ii, jj, in, jn,  iVar
    double precision :: lim(nVar,nSubLimV(1),nSubLimV(2))
    double precision :: lob(nVar,nDOF(1),nDOF(2))
    double precision :: Qmin(nVar), Qmax(nVar)

    ! Step 1) Loop over all elements and get the min/max for each variable

    DO j = 1, JMAX
        DO i = 1, IMAX

            IF(Limiter(i,j)%status.EQ.0) THEN

                CALL GetSubcellData( lim, uh(:,:,:,i,j) )
                CALL GetLobattoData( lob, uh(:,:,:,i,j) )

                DO iVar = 1, nVar

                    Limiter(i,j)%lmin_old(iVar) = MIN( MINVAL(uh(iVar,1:nDOF(1),1:nDOF(2),i,j)), &
                                                   MINVAL(lob(iVar,1:nDOF(1),1:nDOF(2))), &
                                                   MINVAL(lim(iVar,1:nSubLimV(1),1:nSubLimV(2))) )

                    Limiter(i,j)%lmax_old(iVar) = MAX( MAXVAL(uh(iVar,1:nDOF(1),1:nDOF(2),i,j)), &
                                                   MAXVAL(lob(iVar,1:nDOF(1),1:nDOF(2))), &
                                                   MAXVAL(lim(iVar,1:nSubLimV(1),1:nSubLimV(2))) )

                END DO
            ELSE
                lim = Limiter(i,j)%Lh
                DO iVar = 1, nVar
                    Limiter(i,j)%lmin_old(iVar) = MINVAL( lim(iVar,1:nSubLimV(1),1:nSubLimV(2)) )
                    Limiter(i,j)%lmax_old(iVar) = MAXVAL( lim(iVar,1:nSubLimV(1),1:nSubLimV(2)) )
                ENDDO
            END IF

        END DO
    END DO

    ! Step 2) Loop over all Voronoi neighbors and get the final min/max for each variable (requires step 1 to be complete)

    DO j = 1, JMAX
        DO i = 1, IMAX
            Qmin = Limiter(i,j)%lmin_old(:)
            Qmax = Limiter(i,j)%lmax_old(:)

            DO jj = -1, 1
                DO ii = -1, 1
                    in = neighbor(ii,jj,1,i,j)
                    jn = neighbor(ii,jj,2,i,j)
                    DO iVar = 1, nVar
                        Qmin(iVar) = MIN( Qmin(iVar), Limiter(in,jn)%lmin_old(iVar) )
                        Qmax(iVar) = MAX( Qmax(iVar), Limiter(in,jn)%lmax_old(iVar) )
                    END DO
                END DO
            END DO

            Limiter(i,j)%lmin(:) = Qmin
            Limiter(i,j)%lmax(:) = Qmax

        END DO
    END DO

    CONTINUE

END subroutine GetMinMax


!-----------------------------------------------------------------------
! Check discrete maximum principle (DMP) and physical admissibility
! in the given cell
!-----------------------------------------------------------------------

subroutine DMP(dmpresult,arguh,argLimiter,argmax)
    USE ader_dg
    implicit none
    ! Argument list
    logical        :: dmpresult
    double precision           :: arguh(nVar,nDOF(1),nDOF(2)),argmax
    TYPE(tLimiter) :: argLimiter
    ! Local variables
    integer :: iii, jjj, iVar, iErr
    double precision    :: lim(nVar,nSubLimV(1),nSubLimV(2))
    double precision    :: lob(nVar,nDOF(1),nDOF(2))
    double precision    :: lmin(nVar), lmax(nVar), ldiff, Qout(nVar)


    dmpresult = .TRUE.

    CALL GetSubcellData( lim, arguh )
    CALL GetLobattoData( lob, arguh )

    DO iVar = 1, nVar

        IF (iVar .GT. 0) THEN

            lmin(iVar) = MIN( MINVAL(arguh(iVar,1:nDOF(1),1:nDOF(2))), &
                            MINVAL(lim(iVar,1:nSubLimV(1),1:nSubLimV(2))), &
                            MINVAL(lob(iVar,1:nDOF(1),1:nDOF(2))) )

            lmax(iVar) = MAX( MAXVAL(arguh(iVar,1:nDOF(1),1:nDOF(2))), &
                            MAXVAL(lim(iVar,1:nSubLimV(1),1:nSubLimV(2))), &
                            MAXVAL(lob(iVar,1:nDOF(1),1:nDOF(2))) )

            ldiff = MAX( 1.0d-4, 1.0d-3*(argLimiter%lmax(iVar) - argLimiter%lmin(iVar)), argmax )

            IF(lmin(iVar).LT.argLimiter%lmin(iVar)-ldiff) THEN
                dmpresult = .FALSE.
                RETURN
            ENDIF

            IF(lmax(iVar).GT.argLimiter%lmax(iVar)+ldiff) THEN
                dmpresult = .FALSE.
                RETURN
            ENDIF

        END IF
    ENDDO

    ! Physical detection criteria

    DO jjj = 1, nSubLimV(2)
        DO iii = 1, nSubLimV(1)
            CALL PDEAssurePositivity(iErr, lim(:,iii,jjj), Qout )
            IF(iErr.LT.0) THEN
                dmpresult = .FALSE.
                RETURN
            ENDIF
            !DO iVar = 1, nVar TODO: Add this condition for NaN later
                !IF(ISNAN(lim(iVar,iii,jjj))) THEN
                    !dmpresult = .FALSE.
                    !RETURN
                !ENDIF
            !ENDDO
        ENDDO
    ENDDO

   CONTINUE

END subroutine DMP

!-----------------------------------------------------------------------
! Save uh before the limiter is applied
!-----------------------------------------------------------------------

subroutine Saveolduh
    USE ader_dg
    implicit none
    ! Local variables
    integer  :: i,j

    DO j = 1, JMAX
        DO i = 1, IMAX
            olduh(:,:,:,i,j) = uh(:,:,:,i,j)
        ENDDO
    END DO

    CONTINUE

END subroutine Saveolduh

subroutine UpdateLimiter
    USE ader_dg
    implicit none
    ! Local variables
    integer :: i,j

    DO j = 1, JMAX
        DO i = 1, IMAX

            IF(Limiter(i,j)%status.EQ.0) THEN
                IF(Limiter(i,j)%oldstatus.NE.0) THEN
                    DEALLOCATE( Limiter(i,j)%Lh    )
                    DEALLOCATE( Limiter(i,j)%NewLh )
                ENDIF
            ELSE
                Limiter(i,j)%Lh = Limiter(i,j)%NewLh
            ENDIF

            Limiter(i,j)%oldstatus = Limiter(i,j)%status

            CONTINUE
        ENDDO
    END DO

END subroutine UpdateLimiter

!-----------------------------------------------------------------------
! Allocate memory for the limiter
!-----------------------------------------------------------------------

subroutine AllocateLimiter
    USE ader_dg
    implicit none
    ! Local variables
    integer    :: i,j

    DO j = 1, JMAX
        DO i = 1, IMAX
            IF(recompute(i,j).EQ.1) THEN
                IF(Limiter(i,j)%oldstatus.EQ.0) THEN
                    ALLOCATE( Limiter(i,j)%Lh(nVar,nSubLimV(1),nSubLimV(2))    )
                    ALLOCATE( Limiter(i,j)%NewLh(nVar,nSubLimV(1),nSubLimV(2)) )
                ENDIF
                Limiter(i,j)%status = 1
            ELSE
                Limiter(i,j)%status = 0
            ENDIF
        END DO
    END DO

END subroutine AllocateLimiter

!-----------------------------------------------------------------------
! Once a cell is marked for recomputation, also mark the face
! neighbours for recomputation
!-----------------------------------------------------------------------

subroutine SpreadRecompute
    USE ader_dg
    implicit none
    ! Local variables
    integer  :: i, j, ii, jj, in, jn
    integer  :: iFace

    DO j = 1, JMAX
        DO i = 1, IMAX
            IF(recompute(i,j)==1) THEN
                DO iFace = 1, 2*nDim
                    ii = Face2Neigh(1,iFace)
                    jj = Face2Neigh(2,iFace)
                    in = neighbor(ii,jj,1,i,j)
                    jn = neighbor(ii,jj,2,i,j)
                    IF(recompute(in,jn)==0) THEN
                        recompute(in,jn) = 2
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
    END DO

    CONTINUE

 END subroutine SpreadRecompute

 !-----------------------------------------------------------------------
 ! Rusanov flux
 !-----------------------------------------------------------------------

 subroutine RusanovFluxS(flux,QL,QR,nv)
    USE ader_dg
    implicit none
    ! Argument list
    double precision, INTENT(IN)  :: QL(nVar), QR(nVar), nv(nDim)
    double precision, INTENT(OUT) :: flux(nVar)
    ! Local variables
    double precision :: smax, LL(nVar), LR(nVar)
    double precision :: FL(nVar,nDim), FR(nVar,nDim)
    !
    CALL PDEEigenvalues(QL,nv,LL)
    CALL PDEEigenvalues(QR,nv,LR)
    smax = MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) )
    CALL PDEFlux(QL,FL)
    CALL PDEFlux(QR,FR)
    flux = 0.50d0*matmul(FL+FR,nv) - 0.50d0*smax*(QR-QL)
    !
END subroutine

!-----------------------------------------------------------------------
! Min-mod limiter
!-----------------------------------------------------------------------

subroutine minmod(c,a,b,nnn)
    Implicit none
    ! Argument list
    integer, INTENT(IN)   :: nnn
    double precision, INTENT(IN)      :: a(nnn), b(nnn)
    double precision, INTENT(OUT)     :: c(nnn)
    ! Local variables
    integer :: i

    DO i = 1, nnn
        IF(a(i)*b(i).LE.0.0d0) THEN
            c(i) = 0.0d0
        ELSE
            IF( ABS(a(i)).LT.ABS(b(i)) ) THEN
                c(i) = a(i)
            ELSE
                c(i) = b(i)
            ENDIF
        ENDIF
    ENDDO

END subroutine minmod

!-----------------------------------------------------------------------
! Recompute solution in the subcell
!-----------------------------------------------------------------------


subroutine Subcellrecompute
    use constants, only : m_pi
    use ader_dg
    implicit none
    ! Local variables
    integer :: i, j, in, jn, ii, jj
    double precision    :: ldx(nDim), nv(nDim), nvi(nDim,nDim)
    double precision    :: subuh0(nVar,(1-nSubLimV(1)):2*nSubLimV(1),(1-nSubLimV(2)):2*nSubLimV(2))
    double precision    :: subuh1(nVar,nSubLimV(1),nSubLimV(2))
    double precision    :: slopex(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2))
    double precision    :: slopey(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2))
    double precision    :: wLx(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2))
    double precision    :: wLy(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2))
    double precision    :: wRx(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2))
    double precision    :: wRy(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2))
    double precision    :: lim(nVar,nSubLimV(1),nSubLimV(2))
    double precision    :: Qt(nVar)
    double precision    :: FLx(nVar,nDim), FRx(nVar,nDim), FLy(nVar,nDim), FRy(nVar,nDim)
    double precision    :: Fx(nVar,1:nSubLimV(1)+1,1:nSubLimV(2))
    double precision    :: Fy(nVar,1:nSubLimV(1),1:nSubLimV(2)+1)
    double precision    :: Va(nVar), Vb(nVar), Qbc(nVar), x_sw
    logical :: dmpresult

    Va(1) = 8.0d0
    Va(2) = 8.25d0*cos(m_pi/6.0d0)
    Va(3) = -8.25d0*sin(m_pi/6.0d0)
    Va(4) = 116.5d0

    Vb(1) = 1.40d0;
    Vb(2) = 0.0d0;
    Vb(3) = 0.0d0;
    Vb(4) = 1.0d0;

    call shock_position(time, x_sw)

    nvi = 0.0d0

    do ii = 1, nDim
        nvi(ii,ii) = 1.0d0
    end do

    ldx = dx/DBLE(nSubLim)
    Fx = 0.0d0
    Fy = 0.0d0

    slopex = 0.0d0
    slopey = 0.0d0

    wLx = 0.0d0
    wRx = 0.0d0
    wLy = 0.0d0
    wRy = 0.0d0

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
                        call GetSubcellData(lim, olduh(:,:,:,in,jn))
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
                if (dble(i)*dx(1) .gt. 1.0d0/6.0d0) then ! Reflective boundary
                    subuh0(3,:,1-nSubLimV(2):0) = -subuh0(3,:,1-nSubLimV(2):0)
                end if
            end if

            if (j .eq. JMAX) then ! Top

                if (dble(i)*dx(1) .le. x_sw) then
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
                call PutSubCellData(uh(:,:,:,i,j),Limiter(i,j)%NewLh)

           else if (recompute(i,j) .eq. 2) then

                ! Put the evolved cell-averaged data back into the DG polynomial
                call PutSubCellData(uh(:,:,:,i,j),subuh1)

                ! If we recompute neighbors of troubled cells, and the resulting reconstructed DG polynomial is troubled, then activate the limiter also in this case
                call DMP(dmpresult,uh(:,:,:,i,j),Limiter(i,j),0.0d0)

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

end subroutine Subcellrecompute
