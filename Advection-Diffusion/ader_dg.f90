! ader_dg.f90
! author: sunder

!-----------------------------------------------------------------------
! Main module containing all the necessary data structures
!-----------------------------------------------------------------------

module ader_dg

    implicit none
    public

    ! ---------------------------------------------------------------------------------------------------------
    ! Can be modified by the user
    ! ---------------------------------------------------------------------------------------------------------

    integer, parameter  :: N = 1      ! Polynomial degree of our approximation in space and time
    integer, parameter  :: nDim = 2   ! The number of space dimensions that we actually want to simulate
    real, parameter     :: CFL = 0.9  ! The Courant-Friedrichs-Lewy number < 1
    integer, parameter  :: nVar = 1   ! The number of variables of the PDE system

    ! ---------------------------------------------------------------------------------------------------------
    ! Do not change beyond this point
    ! ---------------------------------------------------------------------------------------------------------

    ! The following variables contain important information about the numerical method.

    real, parameter :: PNPMTable(0:9) = (/ 1.0, 0.33, 0.17, 0.1, 0.069, 0.045,  0.038, 0.03, 0.02, 0.015 /)   ! maximum CFL numbers for PNPM schemes
    real, parameter :: m_pi = 4.0*atan(1.0)                ! Value of Pi
    integer         :: nDOF(0:2)                           ! The number of degrees of freedom in space and time
    real            :: xiGPN(N+1), wGPN(N+1)               ! The Gauss-Legendre quadrature nodes and weights
    real            :: xin(N+1)                            ! The nodes used for our basis (can in principle be different from the Gauss-Legendre nodes, like for example also the bad Gauss-Lobatto nodes)
    real            :: MM(N+1,N+1), iMM(N+1,N+1)           ! Element mass-matrix and its inverse
    real            :: Kxi(N+1,N+1), dudx(N+1,N+1)         ! Element stiffness matrix and discrete derivative operator
    real            :: FLm(N+1,N+1), FLp(N+1,N+1)          ! Left flux matrices
    real            :: FRm(N+1,N+1), FRp(N+1,N+1)          ! Right flux matrices
    real            :: FLcoeff(N+1), FRcoeff(N+1)          ! extrapolation coefficients to the left and right boundary
    real            :: F0(N+1), F1(N+1,N+1)                ! Time flux matrices
    real            :: K1(N+1,N+1), iK1(N+1,N+1)           ! F1 - Ktau

    ! Stuff related to the problem setup, mesh and output

    integer             :: IMAX, JMAX                     ! Number of cells in each space dimension
    integer             :: nElem, nFace, nNode            ! Number of elements, faces and nodes
    integer             :: timestep                       ! Number of the current time step
    integer             :: WriteInterval                  ! Number of time steps after which to write output data
    real                :: xL(nDim), xR(nDim)             ! Computational domain
    real                :: dx(nDim), dt                   ! Vector of mesh spacings in each dimension and the time step
    real                :: time, tend                     ! Current time  and final time
    real, pointer       :: x(:,:,:)                       ! Coordinates of cell centers
    integer, parameter  :: nVtx = 2**nDim, nFac = 2*nDim  ! Number of vertices and faces per element
    integer             :: bC(4)                          ! Boundary conditions on the four faces 
    character(len=200)  :: BaseFile                       ! Basic filename to write the results

    ! Some diagnostics                                        !
    real               :: tCPU1, tCPU2                        ! CPU times
    integer(8)         :: TEU                                 ! total element updates

    ! The main variables of the ADER-DG scheme

    real, pointer  :: uh(:,:,:,:,:)               ! Coefficients of the DG polynomial
    real, pointer  :: wh(:,:,:)                   ! Cell average values of primitive variables
    real, pointer  :: duh(:,:,:,:,:)              ! Update coefficients of the DG polynomial
    real, pointer  :: qhi(:,:,:,:,:)              ! Time-averaged coefficients of the space-time predictor
    real, pointer  :: Fhi(:,:,:,:,:,:)            ! Time-averaged coefficients of the flux tensor of the space-time predictor
    real, pointer  :: qbnd(:,:,:,:,:)             ! Boundary-extrapolated data for the state vector Q in the element
    real, pointer  :: gradQbnd(:,:,:,:,:,:)       ! Boundary-extrapolated data for the state vector Q in the elem
    real, pointer  :: Fu(:,:,:,:)                 ! Upwind flux in the x-direction
    real, pointer  :: Gu(:,:,:,:)                 ! Upwind flux in the y-direction

    ! Important info and parameters concerning the governing PDE system

    type tEquations
        real :: a       ! Diffusion coefficient
        real :: cx, cy  ! Advection coefficients in x and y directions
    end type tEquations

    type(tEquations)   :: EQN

end module ader_dg

