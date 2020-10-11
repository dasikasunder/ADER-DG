! surface_integral.f90
! author: sunder

!-----------------------------------------------------------------------
! Find the surface integral
!-----------------------------------------------------------------------

subroutine surface_integral(lduh, lFl, lFr, lGb, lGt)

    use ader_dg
    implicit none
    ! Argument list
    real, intent(in)    :: lFl(nVar,nDOF(2)), lFr(nVar,nDOF(2)), lGb(nVar,nDOF(1)), lGt(nVar,nDOF(1))
    real, intent(inout) :: lduh(nVar, nDOF(1), nDOF(2))       ! Spatial degrees of freedom

    ! Local variables
    integer           :: i,j,iVar

    ! Now multiply the numerical fluxes on the surfaces with the test functions and compute the surface integrals

    ! x faces

    do j = 1, nDOF(2)
        do iVar = 1, nVar
            lduh(iVar,:,j) = lduh(iVar,:,j) - &
                wGPN(j)/dx(1)*( lFr(iVar,j)*FRCoeff - lFl(iVar,j)*FLCoeff )      ! left flux minus right flux
        end do
    end do

    ! y faces

    do i = 1, nDOF(1)
        do iVar = 1, nVar
            lduh(iVar,i,:) = lduh(iVar,i,:) - &
            wGPN(i)/dx(2)*( lGt(iVar,i)*FRCoeff - lGb(iVar,i)*FLCoeff )         ! bottom flux minus top flux
        end do
    end do

end subroutine surface_integral
