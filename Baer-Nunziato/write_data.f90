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

    DO j = 1, JMAX
        DO i = 1, IMAX
            x_c = xL(1) + (real(i)-0.5)*dx(1)
            y_c = xL(2) + (real(j)-0.5)*dx(2)
            write(out_unit, '(E15.7, E15.7, E15.7)')  x_c, y_c, 0.0
        end do
    end do

    ! Plot density (Solid phase)

    write(out_unit, '(A, I10)') 'POINT_DATA', JMAX*IMAX
    write(out_unit, '(A)') 'SCALARS RHO_S double 1'
    write(out_unit, '(A)') 'LOOKUP_TABLE default'

    do j = 1, JMAX
        do i = 1, IMAX
            write(out_unit, '(E15.7)')  wh(1,i,j)
        end do
    end do

    ! Plot velocity (Solid phase)

    write(out_unit, '(A)') 'VECTORS V_S double'

    do j = 1, JMAX
        do i = 1, IMAX
            write(out_unit, '(E15.7, E15.7, E15.7)')  wh(2,i,j), wh(3,i,j), 0.0
        end do
    end do

    ! Plot pressure (Solid phase)

    write(out_unit, '(A)') 'SCALARS P_S double 1'
    write(out_unit, '(A)') 'LOOKUP_TABLE default'

    DO j = 1, JMAX
        DO i = 1, IMAX
            write(out_unit, '(E15.7)')  wh(4,i,j)
        end do
    end do

    ! Plot density (Gas phase)

    write(out_unit, '(A)') 'SCALARS RHO_G double 1'
    write(out_unit, '(A)') 'LOOKUP_TABLE default'

    do j = 1, JMAX
        do i = 1, IMAX
            write(out_unit, '(E15.7)')  wh(5,i,j)
        end do
    end do

    ! Plot velocity (Gas phase)

    write(out_unit, '(A)') 'VECTORS V_G double'

    do j = 1, JMAX
        do i = 1, IMAX
            write(out_unit, '(E15.7, E15.7, E15.7)')  wh(6,i,j), wh(7,i,j), 0.0
        end do
    end do

    ! Plot pressure (Gas phase)

    write(out_unit, '(A)') 'SCALARS P_G double 1'
    write(out_unit, '(A)') 'LOOKUP_TABLE default'

    do j = 1, JMAX
        do i = 1, IMAX
            write(out_unit, '(E15.7)')  wh(8,i,j)
        end do
    end do

   ! Phase field (Solid phase)

    write(out_unit, '(A)') 'SCALARS PHI_S double 1'
    write(out_unit, '(A)') 'LOOKUP_TABLE default'

    do j = 1, JMAX
        do i = 1, IMAX
            write(out_unit, '(E15.7)')  wh(9,i,j)
        end do
    end do

    close(out_unit)

end subroutine write_data

!-----------------------------------------------------------------------
! Plot the troubled cells in vtk (unstructured data) format
!-----------------------------------------------------------------------

subroutine plot_troubled_cells
    use ader_dg
    implicit none
    ! Local variables
    real, pointer :: xnode(:,:)
    integer, pointer :: tri(:,:), idxn(:,:), idxe(:,:)
    integer :: i, j, counter
    character(len=200) :: Filename
    integer, parameter :: out_unit = 25

    allocate( xnode(nDim, nNode)) ! Allocate the nodes
    allocate( idxn(IMAX+1,JMAX+1)  )
    allocate( tri(4, nElem)     )
    allocate( idxe(IMAX,JMAX) )

    ! Define the node coordinates and the node numbers

    counter = 0

    do j = 1, JMAX+1
        do i = 1, IMAX+1
                counter = counter + 1
                xnode(:,counter) = xL(:) + (/ i-1, j-1/)*dx(:)
                idxn(i,j) = counter
        end do
    end do

    ! Define the connectivity between the elements and the nodes

    counter = 0

    do j = 1, JMAX
        do i = 1, IMAX
            counter = counter + 1
            idxe(i,j) = counter
            tri(:,counter) = (/ idxn(i,j), idxn(i+1,j), idxn(i,j+1), idxn(i+1,j+1) /)
        end do
    end do

    write(FileName,'(a,a1,i8.8,a)') 'Tcells','-', timestep, '.vtk'

    open (unit=out_unit, file=Filename)

    write(out_unit, '(A)') '# vtk DataFile Version 2.0'
    write(out_unit, '(A)') '2D Euler'
    write(out_unit, '(A)') 'ASCII'
    write(out_unit, '(A)') 'DATASET UNSTRUCTURED_GRID'
    write(out_unit, '(A)') 'FIELD FieldData 1'
    write(out_unit, '(A)') 'TIME 1 1 double'
    write(out_unit, '(E15.9)')  time
    write(out_unit, '(A, I10, A)') 'POINTS ', nNode, ' double'

    do i = 1, nNode
        write(out_unit, '(E15.7, E15.7, E15.7)')  xnode(1,i), xnode(2,i), 0.0
    end do

    write(out_unit, '(A, I10, I10)') 'CELLS ', nElem, nElem*5

    do i = 1, nElem
        write(out_unit, '(I10, I10, I10, I10, I10)')  4, tri(1,i)-1, tri(2,i)-1, tri(4,i)-1, tri(3,i)-1
    end do

    write(out_unit, '(A, I10, I10)') 'CELL_TYPES ', nElem

    do i = 1, nElem
        write(out_unit, '(I5)')  9
    end do

    write(out_unit, '(A, I10)') 'CELL_DATA', nElem
    write(out_unit, '(A)') 'SCALARS TCells double 1'
    write(out_unit, '(A)') 'LOOKUP_TABLE default'

    do j = 1, JMAX
        do i = 1, IMAX
            if (recompute(i,j) .eq. 1) then
                write(out_unit, '(E15.7)')  1.0
            else
                write(out_unit, '(E15.7)')  0.0
            end if
       end do
    end do

    deallocate( xnode)
    deallocate( idxn  )
    deallocate( tri     )
    deallocate( idxe )

    close(out_unit)

end subroutine plot_troubled_cells

