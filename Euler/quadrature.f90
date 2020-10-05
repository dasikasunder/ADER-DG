! quadrature.f90
! author: sunder

!-----------------------------------------------------------------------
! Various quadrature routines
!-----------------------------------------------------------------------

pure subroutine gauleg(x1, x2, x, w, n)
    use constants, ONLY : m_pi
    implicit none
    ! Argument list
    integer, intent(in)           ::  n
    double precision, intent (in) :: x1, x2
    double precision, intent(out) :: x(n), w(n)
    ! Local variables
    double precision :: EPS
    integer  :: i,j,m
    double precision :: p1, p2, p3, pp, xl, xm, z, z1

    parameter (EPS=1.0D-15)

    m  = (n+1)/2
    xm = 0.5*(x2+x1)
    xl = 0.5*(x2-x1)
    do i=1,m
        z = COS(M_PI*(i-0.25d0)/(n+0.5d0))
    1   continue
        p1 = 1.0d0
        p2 = 0.0d0
        do j = 1,n
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/dble(j)
        end do
        pp = dble(n)*(z*p1-p2)/(z*z-1.0d0)
        z1 = z
        z  = z1-p1/pp
        if (ABS(z-z1) .gt. EPS) goto 1
        x(i)    = xm-xl*z
        x(n+1-i)= xm+xl*z
        w(i)    = 2.0d0*xl/((1.-z*z)*pp*pp)
        w(n+1-i)= w(i)
    end do
    return
end subroutine gauleg


pure subroutine gaulob(x1,x2,x,w,n1)
  use constants, ONLY : m_pi
  !--------------------------------------------------------------------------
  implicit none
  !--------------------------------------------------------------------------
  ! Argument list
  integer          :: n1
  double precision :: x1,x2,x(n1),w(n1)
  ! Local variables
  double precision  :: EPS
  integer           :: i,k
  double precision  :: xold(n1)
  double precision  :: P(n1,n1)
  integer           :: n, iter, maxiter
  !--------------------------------------------------------------------------
  intent(IN)  :: x1,x2,n1
  intent(OUT) :: x,w
  !--------------------------------------------------------------------------
  PARAMETER (EPS=1.0d0-15)
  !--------------------------------------------------------------------------
  !
  n = N1-1;

  ! Use the Chebyshev-Gauss-Lobatto nodes as the first guess

  do i = 0, N
   x(i+1) = cos(M_PI*DBLE(i)/DBLE(N))
  end do
  ! The Legendre Vandermonde Matrix
  P=0.0d0
  !
  ! Compute P_(N) using the recursion relation
  ! Compute its first and second derivatives and
  ! update x using the Newton-Raphson method.
  !
  xold=2.0d0
  !
  maxiter = 100000
  do iter = 1, maxiter
    xold=x
    P(:,1)=1.0d0
    P(:,2)=x
    do k = 2, N
        P(:,k+1)=( (2*k-1)*x*P(:,k)-(k-1)*P(:,k-1) )/dble(k)
    end do
    x=xold-( x*P(:,N1)-P(:,N) )/( N1*P(:,N1) )
    IF(maxval(abs(x-xold)) .le. eps) then
      exit
    end if
  end do
  w = 2.0d0/(dble(N*N1)*P(:,N1)**2)
  w = 0.5d0*w*(x2-x1)
  x = x1 + 0.5d0*(x+1.0d0)*(x2-x1)
  !
end subroutine gaulob
