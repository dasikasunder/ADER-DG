!la.f90
!author: sunder

!-----------------------------------------------------------------------
! Find the inverse of a matrix using pivoted gauss elimination
!-----------------------------------------------------------------------

subroutine MatrixInverse(N, A, iA)
        implicit none
        ! Argument list
        integer, intent(in)   :: N
        real, intent(in)  :: A(N,N)
        real, intent(out) :: iA(N,N)
        ! Local Variables
        integer  :: i,j,ml(1)
        real :: piv
        real :: temp(2*N)
        real :: C(2*N,N)

        C(1:N,:)     = transpose(A)
        C(N+1:2*N,:) = 0.0

        do i = 1, N
            C(N+i,i) = 1.0
        end do

        ! Forward elimination and row swapping (if necessary)

        do i = 1, N

            ! If pivot element is zero, then swap rows

            ml = maxloc(abs(C(i,i:N)))
            j = i - 1 + ml(1)
            temp   = C(:,j)
            C(:,j) = C(:,i)
            C(:,i) = temp

            if (abs(C(i,i)) .le. 1.0e-15) then
                print *, 'ERROR. Matrix is singular!'
                do j = 1, N
                    print *, A(j,:)
                end do
                stop
            end if

            piv    = 1.0/C(i,i)
            C(:,i) = C(:,i)*piv

            do j = i+1, N
                C(:,j) = C(:,j) - C(i,j)*C(:,i)
            end do
        end do

        ! Back substitution

        do i = N,1,-1
            do j = i-1,1,-1
                C(:,j) = C(:,j) - C(i,j)*C(:,i)
            end do
        end do

        iA = transpose( C(N+1:2*N,:) )

end subroutine MatrixInverse
