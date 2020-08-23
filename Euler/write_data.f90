! write_data.f90
! author: sunder

!-----------------------------------------------------------------------
! Output data in vtk format
!-----------------------------------------------------------------------

SUBROUTINE write_data
    USE ader_dg
    IMPLICIT NONE

    INTEGER            :: i, j, ii, jj, iVar
    DOUBLE PRECISION   :: x_c, y_c, Q(nVar), Vav(nVar)
    CHARACTER(LEN=200) :: Filename
    INTEGER, PARAMETER :: out_unit = 20

    ! Find the cell average values of primitive variables

    DO j = 1, JMAX
        DO i = 1, IMAX
            Q = 0.0d0

            DO iVar = 1, nVar
                DO jj = 1, nDOF(2)
                    DO ii = 1, nDOF(1)
                        Q(iVar) = Q(iVar) + wGPN(ii)*wGPN(jj)*uh(iVar, ii, jj, i, j)
                    END DO
                END DO
            END DO

            CALL PDECons2Prim(Q, Vav)

            wh(:,i,j) = Vav
        ENDDO
    ENDDO

    WRITE(FileName,'(a,a1,i8.8,a)') TRIM(BaseFile),'-', timestep, '.vtk'
    PRINT *, ' Writing data to file ', TRIM(FileName)

    OPEN (UNIT=out_unit, FILE=Filename)

    WRITE(out_unit, '(A)') '# vtk DataFile Version 2.0'
    WRITE(out_unit, '(A)') '2D Euler'
    WRITE(out_unit, '(A)') 'ASCII'
    WRITE(out_unit, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE(out_unit, '(A)') 'FIELD FieldData 1'
    WRITE(out_unit, '(A)') 'TIME 1 1 double'
    WRITE(out_unit, '(E15.9)')  time
    WRITE(out_unit, '(A, I5, I5, I5)') 'DIMENSIONS ', IMAX, JMAX, 1
    WRITE(out_unit, '(A, I10, A)') 'POINTS ', IMAX*JMAX, ' double'

    DO j = 1, JMAX
        DO i = 1, IMAX
            x_c = xL(1) + (DBLE(i)-0.5)*dx(1)
            y_c = xL(2) + (DBLE(j)-0.5)*dx(2)
            WRITE(out_unit, '(E15.7, E15.7, E15.7)')  x_c, y_c, 0.0
        ENDDO
    ENDDO

    ! Plot density

    WRITE(out_unit, '(A, I10)') 'POINT_DATA', JMAX*IMAX
    WRITE(out_unit, '(A)') 'SCALARS Density double 1'
    WRITE(out_unit, '(A)') 'LOOKUP_TABLE default'

    DO j = 1, JMAX
        DO i = 1, IMAX
            WRITE(out_unit, '(E15.7)')  wh(1,i,j)
        ENDDO
    ENDDO

    ! Plot velocity

    WRITE(out_unit, '(A)') 'VECTORS Velocity double'

    DO j = 1, JMAX
        DO i = 1, IMAX
            WRITE(out_unit, '(E15.7, E15.7, E15.7)')  wh(2,i,j), wh(3,i,j), 0.0
        ENDDO
    ENDDO

    ! Plot pressure

    WRITE(out_unit, '(A)') 'SCALARS Pressure double 1'
    WRITE(out_unit, '(A)') 'LOOKUP_TABLE default'

    DO j = 1, JMAX
        DO i = 1, IMAX
            WRITE(out_unit, '(E15.7)')  wh(4,i,j)
        ENDDO
    ENDDO

    ! Plot Troubled cells

    WRITE(out_unit, '(A)') 'SCALARS T_cells double 1'
    WRITE(out_unit, '(A)') 'LOOKUP_TABLE default'

    DO j = 1, JMAX
        DO i = 1, IMAX
            IF (recompute(i,j) .EQ. 1) THEN
                WRITE(out_unit, '(E15.7)')  1.0d0
            ELSE
                WRITE(out_unit, '(E15.7)')  0.0d0
            END IF
        END DO
    END DO

    CLOSE(out_unit)

END SUBROUTINE write_data
