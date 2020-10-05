!la.f90
!author: sunder

!-----------------------------------------------------------------------
! Find the inverse of a matrix using pivoted gauss elimination
!-----------------------------------------------------------------------

subroutine MatrixInverse(N, A, iA)
        implicit none
        ! Argument list
        integer, intent(in)   :: N
        double precision, intent(in)  :: A(N,N)
        double precision, intent(out) :: iA(N,N)
        ! Local Variables
        integer  :: i,j,ml(1)
        double precision :: piv
        double precision :: temp(2*N)
        double precision :: C(2*N,N)

        C(1:N,:)     = transpose(A)
        C(N+1:2*N,:) = 0.0d0

        do i = 1, N
            C(N+i,i) = 1.0d0
        end do

        ! Forward elimination and row swapping (if necessary)

        do i = 1, N

            ! If pivot element is zero, then swap rows

            ml = maxloc(abs(C(i,i:N)))
            j = i - 1 + ml(1)
            temp   = C(:,j)
            C(:,j) = C(:,i)
            C(:,i) = temp

            if (C(i,i) .eq. 0.0d0) then
                print *, 'ERROR. Matrix is singular!'
                do j = 1, N
                    print *, A(j,:)
                end do
                stop
            end if

            piv    = 1.0d0/C(i,i)
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
