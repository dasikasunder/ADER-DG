! element_update.f90
! author: sunder

!-----------------------------------------------------------------------
! Update each DOF of the an element
!-----------------------------------------------------------------------

subroutine element_update(luh, lduh)
    use ader_dg
    implicit none
    ! Argument list
    real, intent(inout) :: lduh(nVar,nDOF(1),nDOF(2))      ! spatial degrees of freedom
    real, intent(out)   :: luh(nVar,nDOF(1),nDOF(2))       ! nonlinear flux tensor in each space-time DOF

    ! Local variables
    integer             :: i,j

    ! Multiply with the inverse of the mass matrix. For Gauss-Legendre nodes, we simply need to divide by the Gaussian weights

    do j = 1, nDOF(2)
        do i = 1, nDOF(1)
            lduh(:,i,j) = lduh(:,i,j)/(wGPN(i)*wGPN(j))
        end do
    end do

    ! Finally, sum the contribution to the spatial degrees of freedom

    luh = luh + dt*lduh

end subroutine element_update
