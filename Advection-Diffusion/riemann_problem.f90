! riemann_problem.f90
! author: sunder

!-----------------------------------------------------------------------
! Solve riemann problem on all the faces and apply boundary conditions
!-----------------------------------------------------------------------

subroutine solve_riemann_problem
    use ader_dg
    implicit none

    ! Local variables
    real    :: QL(nVar), gradQL(nVar,nDim), QR(nVar), gradQR(nVar,nDim), nv(nDim)
    integer :: i, j, k, l

    ! Find upwind flux in x-direction

    nv(1) = 1.0; nv(2) = 0.0

    do j = 1, JMAX
        do i = 1, IMAX+1

            do l = 1, nDOF(2)

                if (i .eq. 1) then ! Left boundary (Periodic)

                    QL = qBnd(:,2,l,IMAX,j); gradQL = gradqBnd(:,:,2,l,IMAX,j)
                    QR = qBnd(:,1,l,1,j);    gradQR = gradqBnd(:,:,1,l,1,j)

                else if (i .eq. IMAX+1) then ! Right boundary (Periodic)

                    QL = qBnd(:,2,l,IMAX,j); gradQL = gradqBnd(:,:,2,l,IMAX,j)
                    QR = qBnd(:,1,l,1,j);    gradQR = gradqBnd(:,:,1,l,1,j)

                else ! Internal faces
                    QL = qBnd(:,2,l,i-1,j); gradQL = gradqBnd(:,:,2,l,i-1,j)
                    QR = qBnd(:,1,l,i,j);   gradQR = gradqBnd(:,:,1,l,i,j)

                end if

                call  PDEViscLLFRiemannSolver(Fu(:, l, i, j), QL, gradQL, QR, gradQR, nv)

            end do
        end do
    end do

    ! Find upwind flux in y-direction

    nv(1) = 0.0; nv(2) = 1.0

    do j = 1, JMAX+1
        do i = 1, IMAX

            do k = 1, nDOF(1)

                if (j .eq. 1) then ! Bottom boundary (Periodic)

                    QL = qBnd(:,4,k,i,JMAX); gradQL = gradqBnd(:,:,4,k,i,JMAX)
                    QR = qBnd(:,3,k,i,1);    gradQR = gradqBnd(:,:,3,k,i,1)

                else if (j .eq. JMAX+1) then ! Top boundary (Periodic)

                    QL = qBnd(:,4,k,i,JMAX); gradQL = gradqBnd(:,:,4,k,i,JMAX)
                    QR = qBnd(:,3,k,i,1);    gradQR = gradqBnd(:,:,3,k,i,1)

                else ! Internal faces
                    QL = qBnd(:,4,k,i,j-1); gradQL = gradqBnd(:,:,4,k,i,j-1)
                    QR = qBnd(:,3,k,i,j);   gradQR = gradqBnd(:,:,3,k,i,j)

                end if

                call  PDEViscLLFRiemannSolver(Gu(:, k, i, j), QL, gradQL, QR, gradQR, nv)

            end do
        end do
    end do

    continue

end subroutine solve_riemann_problem
