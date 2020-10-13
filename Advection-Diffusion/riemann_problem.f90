! riemann_problem.f90
! author: sunder

!-----------------------------------------------------------------------
! Solve riemann problem on all the faces and apply boundary conditions
!-----------------------------------------------------------------------

subroutine solve_riemann_problem
    use ader_dg
    implicit none

    ! Local variables
    real    :: QL(nVar), gradQL(nVar,nDim), QR(nVar), gradQR(nVar,nDim), nv(nDim), y
    integer :: i, j, k, l

    ! Find upwind flux in x-direction

    nv(1) = 1.0; nv(2) = 0.0

    do j = 1, JMAX
        do i = 1, IMAX+1

            do l = 1, nDOF(2)

                if (i .eq. 1) then ! Left boundary 
    
                    if (bC(1) .eq. 1) then ! Neumann 
                        QL = qBnd(:,1,l,1,j); gradQL = gradqBnd(:,:,1,l,1,j)
                        QR = qBnd(:,1,l,1,j); gradQR = gradqBnd(:,:,1,l,1,j)
                    else if (bC(1) .eq. 2) then
                        y = (real(j-1) + xiGPN(l))*dx(2) 
                        QR = qBnd(:,2,l,IMAX,j); gradQR = gradqBnd(:,:,2,l,IMAX,j)
                        QL = -QR + 2.0*sin(m_pi*y); gradQL = gradQR
                    else if (bC(1) .eq. 3) then ! Periodic
                        QL = qBnd(:,2,l,IMAX,j); gradQL = gradqBnd(:,:,2,l,IMAX,j)
                        QR = qBnd(:,1,l,1,j);    gradQR = gradqBnd(:,:,1,l,1,j)
                    end if
                    
                else if (i .eq. IMAX+1) then ! Right boundary (Periodic)
                    
                    if (bC(2) .eq. 1) then ! Neumann 
                        QL = qBnd(:,2,l,IMAX,j); gradQL = gradqBnd(:,:,2,l,IMAX,j)
                        QR = qBnd(:,2,l,IMAX,j); gradQR = gradqBnd(:,:,2,l,IMAX,j)
                    else if (bC(2) .eq. 2) then ! Dirichlet  
                        QL = qBnd(:,2,l,IMAX,j); gradQL = gradqBnd(:,:,2,l,IMAX,j)
                        QR = -QL; gradQR = gradQL
                    else if (bC(2) .eq. 3) then ! Periodic 
                        QL = qBnd(:,2,l,IMAX,j); gradQL = gradqBnd(:,:,2,l,IMAX,j)
                        QR = qBnd(:,1,l,1,j);    gradQR = gradqBnd(:,:,1,l,1,j)
                    end if
                    
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

                if (j .eq. 1) then ! Bottom boundary
                    
                    if (bC(3) .eq. 1) then ! Neumann 
                        QL = qBnd(:,3,k,i,1); gradQL = gradqBnd(:,:,3,k,i,1)
                        QR = qBnd(:,3,k,i,1); gradQR = gradqBnd(:,:,3,k,i,1)
                    else if (bC(3) .eq. 2) then ! Dirichlet  
                        QR = qBnd(:,3,k,i,1); gradQR = gradqBnd(:,:,3,k,i,1)
                        QL = -QR; gradQL = gradQR
                    else if (bC(3) .eq. 3) then ! Periodic
                        QL = qBnd(:,4,k,i,JMAX); gradQL = gradqBnd(:,:,4,k,i,JMAX)
                        QR = qBnd(:,3,k,i,1);    gradQR = gradqBnd(:,:,3,k,i,1)
                    end if
                
                else if (j .eq. JMAX+1) then ! Top boundary 
                    
                    if (bC(4) .eq. 1) then ! Neumann 
                        QL = qBnd(:,4,k,i,JMAX); gradQL = gradqBnd(:,:,4,k,i,JMAX)
                        QR = qBnd(:,4,k,i,JMAX); gradQR = gradqBnd(:,:,4,k,i,JMAX)
                    else if (bC(4) .eq. 2) then ! Dirichlet
                        QL = qBnd(:,4,k,i,JMAX); gradQL = gradqBnd(:,:,4,k,i,JMAX)
                        QR = -QL; gradQR = gradQL
                    else if (bC(4) .eq. 3) then ! Periodic 
                        QL = qBnd(:,4,k,i,JMAX); gradQL = gradqBnd(:,:,4,k,i,JMAX)
                        QR = qBnd(:,3,k,i,1);    gradQR = gradqBnd(:,:,3,k,i,1)
                    end if
                
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
