subroutine shock_position(t, x_sw)
    real, intent(in)  :: t
    real, intent(out) :: x_sw
    real, parameter :: m_pi = 4.0*atan(1.0)

    x_sw = 1.0/6.0 + 1.0/tan(m_pi/3.0) + (10.0/sin(m_pi/3.0))*t

end subroutine

subroutine solve_riemann_problem
    use ader_dg
    implicit none

    ! Local variables
    real :: QL(nVar), QR(nVar), FL(nVar), FR(nVar), nv(nDim)
    real :: Va(nVar), Vb(nVar), Qbc(nVar), Fbc(nVar, nDim), x_sw
    real, parameter :: m_pi = 4.0*atan(1.0)
    integer  :: i, j, k, l

    Va(1) = 8.0
    Va(2) = 8.25*cos(m_pi/6.0)
    Va(3) = -8.25*sin(m_pi/6.0)
    Va(4) = 116.5

    Vb(1) = 1.4;
    Vb(2) = 0.0;
    Vb(3) = 0.0;
    Vb(4) = 1.0;

    call shock_position(time, x_sw)

    ! Find upwind flux in x-direction

    nv(1) = 1.0; nv(2) = 0.0

    do j = 1, JMAX
        do i = 1, IMAX+1

            do l = 1, nDOF(2)

                if (i .eq. 1) then ! Left boundary

                    call PDEPrim2Cons(Va, Qbc)
                    call PDEFlux(Qbc, Fbc)

                    QL = Qbc; FL = matmul(Fbc,nv)

                    QR = qBnd(:,1,l,i,j);   FR = FBnd(:,1,l,i,j)

                else if (i .eq. IMAX+1) then ! Right boundary

                    QL = qBnd(:,2,l,i-1,j); FL = FBnd(:,2,l,i-1,j)
                    QR = qBnd(:,2,l,i-1,j); FR = FBnd(:,2,l,i-1,j)

                else ! Internal faces
                    QL = qBnd(:,2,l,i-1,j); FL = FBnd(:,2,l,i-1,j)
                    QR = qBnd(:,1,l,i,j);   FR = FBnd(:,1,l,i,j)

                end if

                call RusanovFlux(QL, FL, QR, FR, nv, Fu(:, l, i, j))

            end do
        end do
    end do

    ! Find upwind flux in y-direction

    nv(1) = 0.0; nv(2) = 1.0

    do j = 1, JMAX+1
        do i = 1, IMAX

            do k = 1, nDOF(1)

                if (j .eq. 1) then ! Bottom boundary

                    if (real(i)*dx(1) .le. 1.0/6.0) then ! Transmissive
                        QL = qBnd(:,3,k,i,j);   FL = FBnd(:,3,k,i,j)
                    else ! Reflective
                        QL = qBnd(:,3,k,i,j);   FL = FBnd(:,3,k,i,j)
                        QL(3) = -QL(3)
                        FL(1) = -FL(1); FL(2) = -FL(2); FL(4) = -FL(4);
                    end if

                    QR = qBnd(:,3,k,i,j);   FR = FBnd(:,3,k,i,j)

                else if (j .eq. JMAX+1) then ! Top boundary

                    QL = qBnd(:,4,k,i,j-1); FL = FBnd(:,4,k,i,j-1)

                    if (real(i)*dx(1) .le. x_sw) then
                        call PDEPrim2Cons(Va, Qbc)
                    else
                        call PDEPrim2Cons(Vb, Qbc)
                    end if

                    call PDEFlux(Qbc, Fbc)

                    QR = Qbc; FR = matmul(Fbc,nv)


               else ! Internal faces
                    QL = qBnd(:,4,k,i,j-1); FL = FBnd(:,4,k,i,j-1)
                    QR = qBnd(:,3,k,i,j);   FR = FBnd(:,3,k,i,j)
                end if

                call RusanovFlux(QL, FL, QR, FR, nv, Gu(:, k, i, j))

            end do
        end do
    end do

    continue

end subroutine solve_riemann_problem
