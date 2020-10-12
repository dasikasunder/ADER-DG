! quadrature.f90
! author: sunder

!-----------------------------------------------------------------------
! Various quadrature routines
!-----------------------------------------------------------------------

pure subroutine gauleg(x1, x2, x, w, n)
    implicit none
    ! Argument list
    integer, intent(in)           ::  n
    real, intent (in) :: x1, x2
    real, intent(out) :: x(n), w(n)
    ! Local variables
    real, parameter :: EPS = 1.0e-15
    real, parameter :: m_pi = 4.0*atan(1.0)
    integer  :: i,j,m
    real :: p1, p2, p3, pp, xl, xm, z, z1


    m  = (n+1)/2
    xm = 0.5*(x2+x1)
    xl = 0.5*(x2-x1)
    do i=1,m
        z = COS(M_PI*(i-0.25)/(n+0.5))
    1   continue
        p1 = 1.0
        p2 = 0.0
        do j = 1,n
            p3 = p2
            p2 = p1
            p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/real(j)
        end do
        pp = real(n)*(z*p1-p2)/(z*z-1.0)
        z1 = z
        z  = z1-p1/pp
        if (abs(z-z1) .gt. EPS) goto 1
        x(i)    = xm-xl*z
        x(n+1-i)= xm+xl*z
        w(i)    = 2.0*xl/((1.0-z*z)*pp*pp)
        w(n+1-i)= w(i)
    end do
    return
end subroutine gauleg


pure subroutine gaulob(x1,x2,x,w,n1)

  !--------------------------------------------------------------------------
  implicit none
  !--------------------------------------------------------------------------
  ! Argument list
  integer          :: n1
  real :: x1,x2,x(n1),w(n1)
  ! Local variables
  real, parameter  :: EPS = 1.0e-15
  real, parameter  :: m_pi = 4.0*atan(1.0)
  integer           :: i,k
  real  :: xold(n1)
  real  :: P(n1,n1)
  integer           :: n, iter, maxiter
  !--------------------------------------------------------------------------
  intent(IN)  :: x1,x2,n1
  intent(OUT) :: x,w
  !--------------------------------------------------------------------------
  
  !
  n = N1-1;

  ! Use the Chebyshev-Gauss-Lobatto nodes as the first guess

  do i = 0, N
   x(i+1) = cos(m_pi*real(i)/real(N))
  end do
  ! The Legendre Vandermonde Matrix
  P=0.0
  !
  ! Compute P_(N) using the recursion relation
  ! Compute its first and second derivatives and
  ! update x using the Newton-Raphson method.
  !
  xold=2.0
  !
  maxiter = 100000
  do iter = 1, maxiter
    xold=x
    P(:,1)=1.0
    P(:,2)=x
    do k = 2, N
        P(:,k+1)=( (2*k-1)*x*P(:,k)-(k-1)*P(:,k-1) )/real(k)
    end do
    x=xold-( x*P(:,N1)-P(:,N) )/( N1*P(:,N1) )
    IF(maxval(abs(x-xold)) .le. eps) then
      exit
    end if
  end do
  w = 2.0/(real(N*N1)*P(:,N1)**2)
  w = 0.5*w*(x2-x1)
  x = x1 + 0.5*(x+1.0)*(x2-x1)
  !
end subroutine gaulob
