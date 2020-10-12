! volume_integral.f90
! author: sunder

!-----------------------------------------------------------------------
! Find the element volume integral
!-----------------------------------------------------------------------

subroutine volume_integral(lduh, lFhi, lShi)
    use ader_dg
    implicit none
    ! Argument list
    real, intent(in)  :: lFhi(nVar,nDim,nDOF(1),nDOF(2)) ! nonlinear flux tensor in each space-time DOF
    real, intent(in)  :: lShi(nVar,nDOF(1),nDOF(2))      ! nonlinear source/non-conservative vector in each space DOF
    real, intent(out) :: lduh(nVar,nDOF(1),nDOF(2))      ! spatial degrees of freedom
    ! Local variables
    integer           :: i,j

    ! Initialize the update DOF

    lduh = 0.0

    ! x - direction

    do j = 1, nDOF(2)
        lduh(:,:,j) = lduh(:,:,j) + matmul( lFhi(:,1,:,j), transpose(Kxi) )*wGPN(j)/dx(1)
    end do

    ! y - direction

    do i = 1, nDOF(1)
        lduh(:,i,:) = lduh(:,i,:) + matmul( lFhi(:,2,i,:), transpose(Kxi) )*wGPN(i)/dx(2)
    end do

    ! Add source/non-conservative term

    do j = 1, nDOF(2)
        do i = 1, nDOF(1)
            lduh(:,i,j) = lduh(:,i,j) + wGPN(i)*wGPN(j)*lShi(:,i,j)
        end do
    end do

    continue

end subroutine volume_integral
