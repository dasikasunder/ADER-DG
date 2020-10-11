! main.f90
! author: sunder

!-----------------------------------------------------------------------
! Main program of the code. Puts everything together
!-----------------------------------------------------------------------

program ADERDG
    use ader_dg
    implicit none
    ! Local variables

    integer :: i, j, nRecompute
    logical       :: dmpresult

    ! Initialize important data structures and the solution

    call initialize

    ! Main loop in time

    call CPU_TIME(tCPU1)

    do while (time < tend)

        if(mod(timestep,WriteInterval) .eq. 0) then
            call write_data
            call plot_troubled_cells
        end if

        ! Compute the time step size according to the CFL condition

        call calc_time_step

        print *, ' nRecomp = ', nRecompute, ' t = ', time, 'dt = ', dt

        if (N .gt. 0) then
            call get_min_max
            recompute(:,:) = 0
            nRecompute   = 0
            call save_olduh
        end if

        ! ADER predictor step

        do j = 1,  JMAX
            do i = 1, IMAX
                call ader_space_time_predictor(qhi(:,:,:,i,j), Fhi(:,:,:,:,i,j), qBnd(:,:,:,i,j), FBnd(:,:,:,i,j), uh(:,:,:,i,j))
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

                if (N .gt. 0) then
                    call DMP(dmpresult,uh(:,:,:,i, j), Limiter(i,j), 0.0)
                    if( .not. dmpresult) then
                        recompute(i,j) = 1
                        nRecompute = nRecompute + 1
                    end if
                end if

            end do
        end do

        if ( N .gt. 0 ) then
            call spread_recompute
            call allocate_limiter
            call subcell_recompute
            call update_limiter
        end if

        time = time + dt
        timestep = timestep + 1

    end do

    print*, ' nRecomp = ', nRecompute, ' t = ', time, 'dt = ', dt
    call write_data
    call plot_troubled_cells

    call CPU_TIME(tCPU2)
    print*, ' Total CPU time = ', tCPU2-tCPU1

    call analyze_soln

end program ADERDG
