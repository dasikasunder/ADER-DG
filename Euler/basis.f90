! basis.f90
! author: sunder

!-----------------------------------------------------------------------
! Basis functions. Here we are using Lagrange polynomials in the
! interval [0,1]
!-----------------------------------------------------------------------

subroutine basis_func_1d(phi,phi_xi,xi)
    use ader_dg
    implicit none
    ! Argument list
    double precision, intent(in ) :: xi                    ! coordinate in [0,1] where to evaluate the basis
    double precision, intent(out) :: phi(N+1), phi_xi(N+1) ! the basis and its derivative w.r.t. xi
    ! Local variables
    integer            :: i,j,m
    double precision   :: tmp

    ! Initialize variables
    phi      = 1.0d0
    phi_xi   = 0.0d0

    ! Lagrange polynomial and its derivative
    do m = 1, N+1
       do j = 1, N+1
          if(j .eq. m) cycle
          phi(m) = phi(m)*(xi-xin(j))/(xin(m)-xin(j))
       end do
       do i = 1, N+1
          if(i .eq. m) cycle
          tmp = 1.0d0
          do j = 1, N+1
             if(j .eq. i) cycle
             if(j .eq. m) cycle
             tmp = tmp*(xi-xin(j))/(xin(m)-xin(j))
          end do
          phi_xi(m) = phi_xi(m) + tmp/(xin(m)-xin(i))
       end do
    end do

end subroutine basis_func_1d
