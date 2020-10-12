! riemann_problem.f90
! author: sunder

!-----------------------------------------------------------------------
! Solve the riemann problem on each face of the grid. Also apply 
! appropriate boundary conditions 
!-----------------------------------------------------------------------

subroutine solve_riemann_problem
    use ader_dg
    implicit none

    ! Local variables
    real :: QL(nVar), QR(nVar), FL(nVar), FR(nVar), nv(nDim)
    integer  :: i, j, k, l

    ! Find upwind flux in x-direction

    nv(1) = 1.0; nv(2) = 0.0

    do j = 1, JMAX
        do i = 1, IMAX+1

            do l = 1, nDOF(2)

                if (i .eq. 1) then ! Left boundary

                    QL = qBnd(:,1,l,i,j);   FL = FBnd(:,1,l,i,j)
                    QR = qBnd(:,1,l,i,j);   FR = FBnd(:,1,l,i,j)

                else if (i .eq. IMAX+1) then ! Right boundary

                    QL = qBnd(:,2,l,i-1,j); FL = FBnd(:,2,l,i-1,j)
                    QR = qBnd(:,2,l,i-1,j); FR = FBnd(:,2,l,i-1,j)

                else ! Internal faces
                    QL = qBnd(:,2,l,i-1,j); FL = FBnd(:,2,l,i-1,j)
                    QR = qBnd(:,1,l,i,j);   FR = FBnd(:,1,l,i,j)

                end if

                call RusanovFlux(QL, FL, QR, FR, nv, Fm(:, l, i, j), Fp(:, l, i, j))

            end do
        end do
    end do

    ! Find upwind flux in y-direction

    nv(1) = 0.0; nv(2) = 1.0

    do j = 1, JMAX+1
        do i = 1, IMAX

            do k = 1, nDOF(1)

                if (j .eq. 1) then ! Bottom boundary

                    QL = qBnd(:,3,k,i,j);   FL = FBnd(:,3,k,i,j)
                    QR = qBnd(:,3,k,i,j);   FR = FBnd(:,3,k,i,j)

                else if (j .eq. JMAX+1) then ! Top boundary

                    QL = qBnd(:,4,k,i,j-1); FL = FBnd(:,4,k,i,j-1)
                    QR = qBnd(:,4,k,i,j-1); FR = FBnd(:,4,k,i,j-1)

                else ! Internal faces
                    QL = qBnd(:,4,k,i,j-1); FL = FBnd(:,4,k,i,j-1)
                    QR = qBnd(:,3,k,i,j);   FR = FBnd(:,3,k,i,j)
                end if

                call RusanovFlux(QL, FL, QR, FR, nv, Gm(:, k, i, j), Gp(:, k, i, j))

            end do
        end do
    end do

    continue

end subroutine solve_riemann_problem
