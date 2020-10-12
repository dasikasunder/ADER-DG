! calc_time_step.f90
! author: sunder

!-----------------------------------------------------------------------
! Calculates time step for a general non-linear PDE for an explicit DG
! method
!-----------------------------------------------------------------------

subroutine calc_time_step
    use ader_dg
    implicit none

    ! Local variables

    integer  :: iDim, i, j, k, l
    real :: Lambda(nVar), Lambdav(nVar)
    real :: nv(nDim,nDim), denom, smax_c, smax_v, factor

    factor = real(2*N + 1)

    ! Normal vectors pointing into the two space dimensions

    nv = 0.0

    do i = 1, nDim
        nv(i,i) = 1.0
    end do

    ! First calculate the minimum time step on the current processor

    dt = 1.0e20

    do j = 1, JMAX
        do i = 1, IMAX

            do k = 1, nDOF(1)
                do l = 1, nDOF(2)
                    denom = 0.0
                    do iDim = 1, nDim
                        call PDEEigenvalues(Lambda,uh(:,k,l,i,j),nv(:,iDim))
                        call PDEViscEigenvalues(Lambdav,uh(:,k,l,i,j))
                        smax_c = maxval(abs(Lambda))
                        smax_v = maxval(abs(Lambdav))
                        denom = denom + (smax_c + 2.0*smax_v*factor/dx(iDim))/dx(iDim)
                    end do
                    dt = min(dt, CFL*PNPMTable(N)/denom )
                end do
            end do

        end do
    end do

    if (time+dt .gt. tend) then
        dt = tend - time
    end if

end subroutine calc_time_step
