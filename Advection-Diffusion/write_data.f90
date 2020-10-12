! write_data.f90
! author: sunder

!-----------------------------------------------------------------------
! Output data in vtk format
!-----------------------------------------------------------------------

subroutine write_data
    use ader_dg
    implicit none

    integer            :: i, j, ii, jj, iVar
    real   :: x_c, y_c, Q(nVar), Vav(nVar)
    character(len=200) :: Filename
    integer, parameter :: out_unit = 20

    ! Find the cell average values of primitive variables

    do j = 1, JMAX
        do i = 1, IMAX
            Q = 0.0

            do iVar = 1, nVar
                do jj = 1, nDOF(2)
                    do ii = 1, nDOF(1)
                        Q(iVar) = Q(iVar) + wGPN(ii)*wGPN(jj)*uh(iVar, ii, jj, i, j)
                    end do
                end do
            end do

            call PDECons2Prim(Q, Vav)

            wh(:,i,j) = Vav
        end do
    end do

    write(FileName,'(a,a1,i8.8,a)') TRIM(BaseFile),'-', timestep, '.vtk'
    print *, ' Writing data to file ', TRIM(FileName)

    open (UNIT=out_unit, FILE=Filename)

    write(out_unit, '(A)') '# vtk DataFile Version 2.0'
    write(out_unit, '(A)') '2D Euler'
    write(out_unit, '(A)') 'ASCII'
    write(out_unit, '(A)') 'DATASET STRUCTURED_GRID'
    write(out_unit, '(A)') 'FIELD FieldData 1'
    write(out_unit, '(A)') 'TIME 1 1 double'
    write(out_unit, '(E15.9)')  time
    write(out_unit, '(A, I5, I5, I5)') 'DIMENSIONS ', IMAX, JMAX, 1
    write(out_unit, '(A, I10, A)') 'POINTS ', IMAX*JMAX, ' double'

    do j = 1, JMAX
        do i = 1, IMAX
            x_c = xL(1) + (real(i)-0.5)*dx(1)
            y_c = xL(2) + (real(j)-0.5)*dx(2)
            write(out_unit, '(E15.7, E15.7, E15.7)')  x_c, y_c, 0.0
        end do
    end do

    ! Plot scalar 'U'

    write(out_unit, '(A, I10)') 'POINT_DATA', JMAX*IMAX
    write(out_unit, '(A)') 'SCALARS U double 1'
    write(out_unit, '(A)') 'LOOKUP_TABLE default'

    do j = 1, JMAX
        do i = 1, IMAX
            write(out_unit, '(E15.7)')  wh(1,i,j)
        end do
    end do

    close(out_unit)

end subroutine write_data
