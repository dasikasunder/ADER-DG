! analyze_soln.f90
! author: sunder

!-----------------------------------------------------------------------
! Analyzes the solution and prints L_inf error compared to the initial
! condition. Useful for smooth periodic test cases to check convergence
! rate
!-----------------------------------------------------------------------

subroutine analyze_soln
    use ader_dg
    implicit none

    real, pointer :: uh_exact(:,:,:,:,:), error_density(:, :)
    integer :: i, j, k, l
    real :: u0(nVar), xGP(nDim), corner(nDim), rho_ave, rho_ave_exact

    allocate(  uh_exact(nVar, nDOF(1), nDOF(2), IMAX, JMAX) )
    allocate(  error_density(IMAX, JMAX) )

    do j = 1, JMAX
        do i = 1, IMAX
            corner(:) = x(:,i,j) - 0.5*dx(:) ! coordinate of the lower left corner of the cell
            do l = 1, nDOF(2)
                do k = 1, nDOF(1)
                    xGP = corner + (/ xiGPN(k), xiGPN(l) /)*dx(:)
                    call InitialField(xGP, u0)
                    uh_exact(:, k, l, i, j) = u0(:)
                end do
            end do
        end do
    end do

    do j = 1, JMAX
        do i = 1, IMAX

            rho_ave= 0.0
            rho_ave_exact = 0.0

            do l = 1, nDOF(2)
                do k = 1, nDOF(1)
                    rho_ave_exact = rho_ave_exact + wGPN(l)*wGPN(k)*uh_exact(1, k, l, i, j)
                    rho_ave = rho_ave + wGPN(l)*wGPN(k)*uh(1, k, l, i, j)
                end do
            end do

            error_density(i,j) = abs(rho_ave_exact - rho_ave)

        end do
    end do

    print *, maxval(error_density)

    deallocate( uh_exact, error_density)

end subroutine analyze_soln

