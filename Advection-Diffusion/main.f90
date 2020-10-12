! main.f90
! author: sunder

!-----------------------------------------------------------------------
! Main program of the code. Puts everything together
!-----------------------------------------------------------------------

program ADERDG
    use ader_dg
    implicit none
    ! Local variables

    integer :: i, j

    ! Initialize important data structures and the solution

    call initialize

    ! Main loop in time

    call cpu_time(tCPU1)

    do while (time .lt. tend)

        if(mod(timestep,WriteInterval) .eq. 0) then
            !call write_data
        end if

        ! Compute the time step size according to the CFL condition

        call calc_time_step

        print*, ' t = ', time, 'dt = ', dt

        ! ADER predictor step

        do j = 1,  JMAX
            do i = 1, IMAX
                call ader_space_time_predictor(qhi(:,:,:,i,j), Fhi(:,:,:,:,i,j), qBnd(:,:,:,i,j), gradQBnd(:,:,:,:,i,j),&
                                               uh(:,:,:,i,j))
            end do
        end do

        ! Compute the element volume integral

        do j = 1, JMAX
            do i = 1, IMAX
                call  volume_integral(duh(:,:,:,i,j),Fhi(:,:,:,:,i,j))
            end do
        end do

        ! Solve Riemann problem and also apply boundary conditions

        call solve_riemann_problem

        ! Compute the surface integrals of the test function multiplied with the numerical flux

        do j = 1, JMAX
            do i = 1, IMAX
                call surface_integral(duh(:,:,:,i,j), Fu(:,:,i,j), Fu(:,:,i+1,j), Gu(:,:,i,j), Gu(:,:,i,j+1))
            end do
        end do

        ! Do the element update

        do j = 1, JMAX
            do i = 1, IMAX

                call element_update(uh(:,:,:,i,j),duh(:,:,:,i,j))

            end do
        end do


        time = time + dt
        timestep = timestep + 1

    end do

    print*, ' t = ', time, 'dt = ', dt
    call write_data

    call cpu_time(tCPU2)
    print*, ' Total CPU time = ', tCPU2-tCPU1

    call analyze_soln

end program ADERDG
