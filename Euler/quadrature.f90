! quadrature.f90
! author: sunder

!-----------------------------------------------------------------------
! Various quadrature routines
!-----------------------------------------------------------------------

PURE SUBROUTINE gauleg(x1, x2, x, w, n)
    USE constants, ONLY : M_PI
    IMPLICIT NONE
    ! Argument list
    INTEGER, INTENT(IN)           ::  n
    DOUBLE PRECISION, INTENT (IN) :: x1, x2
    DOUBLE PRECISION, INTENT(OUT) :: x(n), w(n)
    ! Local variables
    DOUBLE PRECISION :: EPS
    INTEGER  :: i,j,m
    DOUBLE PRECISION :: p1, p2, p3, pp, xl, xm, z, z1

    PARAMETER (EPS=1.0D-15)

    m  = (n+1)/2
    xm = 0.5*(x2+x1)
    xl = 0.5*(x2-x1)
    DO i=1,m
        z = COS(M_PI*(i-0.25d0)/(n+0.5d0))
    1    CONTINUE
        p1 = 1.0d0
        p2 = 0.0d0
        DO j = 1,n
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/DBLE(j)
        END DO
        pp = DBLE(n)*(z*p1-p2)/(z*z-1.0d0)
        z1 = z
        z  = z1-p1/pp
        IF(ABS(z-z1).GT.EPS)GOTO 1
        x(i)    = xm-xl*z
        x(n+1-i)= xm+xl*z
        w(i)    = 2.0d0*xl/((1.-z*z)*pp*pp)
        w(n+1-i)= w(i)
    END DO
    RETURN
END SUBROUTINE gauleg


PURE SUBROUTINE gaulob(x1,x2,x,w,n1)
  USE constants, ONLY : M_PI
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  ! Argument list
  INTEGER          :: n1
  DOUBLE PRECISION :: x1,x2,x(n1),w(n1)
  ! Local variables
  DOUBLE PRECISION  :: EPS
  INTEGER           :: i,k
  DOUBLE PRECISION  :: xold(n1)
  DOUBLE PRECISION  :: P(n1,n1)
  INTEGER           :: n, iter, maxiter
  !--------------------------------------------------------------------------
  INTENT(IN)  :: x1,x2,n1
  INTENT(OUT) :: x,w
  !--------------------------------------------------------------------------
  PARAMETER (EPS=1.0d0-15)
  !--------------------------------------------------------------------------
  !
  n = N1-1;

  ! Use the Chebyshev-Gauss-Lobatto nodes as the first guess

  DO i = 0, N
   x(i+1) = cos(M_PI*DBLE(i)/DBLE(N))
  ENDDO
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
  DO iter = 1, maxiter
    xold=x
    P(:,1)=1.0d0
    P(:,2)=x
    DO k = 2, N
        P(:,k+1)=( (2*k-1)*x*P(:,k)-(k-1)*P(:,k-1) )/DBLE(k)
    ENDDO
    x=xold-( x*P(:,N1)-P(:,N) )/( N1*P(:,N1) )
    IF(MAXVAL(ABS(x-xold)).LE.eps) THEN
      EXIT
    ENDIF
  ENDDO
  w = 2.0d0/(DBLE(N*N1)*P(:,N1)**2)
  w = 0.5d0*w*(x2-x1)
  x = x1 + 0.5d0*(x+1.0d0)*(x2-x1)
  !
END SUBROUTINE gaulob
