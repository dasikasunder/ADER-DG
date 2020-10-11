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
    integer, parameter  :: nVar = 4   ! The number of variables of the PDE system

    ! ---------------------------------------------------------------------------------------------------------
    ! Do not change beyond this point
    ! ---------------------------------------------------------------------------------------------------------

    ! The following variables contain important information about the numerical method.

    real, parameter :: PNPMTable(0:9) = (/ 1.0, 0.33, 0.17, 0.1, 0.069, 0.045,  0.038, 0.03, 0.02, 0.015 /)   ! maximum CFL numbers for PNPM schemes
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

    integer            :: IMAX, JMAX                     ! Number of cells in each space dimension
    integer            :: nElem, nFace, nNode            ! Number of elements, faces and nodes
    integer            :: timestep                       ! Number of the current time step
    integer            :: WriteInterval                  ! Number of time steps after which to write output data
    real               :: xL(nDim), xR(nDim)             ! Computational domain
    real               :: dx(nDim), dt                   ! Vector of mesh spacings in each dimension and the time step
    real               :: time, tend                     ! Current time  and final time
    real, pointer      :: x(:,:,:)                       ! Coordinates of cell centers
    integer, parameter :: nVtx = 2**nDim, nFac = 2*nDim  ! Number of vertices and faces per element
    integer            :: ICType                         ! Type of initial condition to be used
    integer            :: bL, bR, bB, bT                 ! Boundary conditions on the left, right, bottom and top of the domain
    character(len=200) :: BaseFile                       ! Basic filename to write the results

    ! Some diagnostics                                   !
    real               :: tCPU1, tCPU2                   ! CPU times
    integer(8)         :: TEU                            ! total element updates

    ! Data needed for the subcell limiter

    integer, parameter :: nSubLim = 2*N+1             ! Number of subcells
    integer            :: nSubLimV(nDim)              ! Vector for number of subcells in each dimension
    real, pointer      :: uh2lim(:,:), lim2uh(:,:)    ! Mapping from DG polynomials to the subcells and back
    real, pointer      :: uh2lob(:,:)                 ! Mapping from the DG polynomials to the Gauss-Lobatto points
    real               :: xiLob(N+1), wLob(N+1)       ! Gauss lobatto points and weights
    integer, pointer   :: neighbor(:,:,:,:,:)         ! Set of Voronoi neighbors of a cell
    real, pointer      :: olduh(:,:,:,:,:)            ! for the a posteriori limiter, we need the possibility to go back
    integer, pointer   :: recompute(:,:)              ! Map containing the troubled zones that need to be recomputed
    real               :: xilimbary(nSubLim)          ! Barycenters of the subcells
    integer            :: Face2Neigh(2,4)             ! Mapping from face to neighbor index
    integer            :: dn(2)                       ! number of direct neighbors in each dimension

    type tLimiter
    integer        :: status, oldstatus
    real, pointer  :: Lh(:,:,:)
    real, pointer  :: NewLh(:,:,:)
    real           :: lmin(nVar), lmax(nVar)
    real           :: lmin_old(nVar), lmax_old(nVar)
    end type tLimiter

    type(tLimiter), pointer :: Limiter(:,:)

    ! The main variables of the ADER-DG scheme

    real, pointer  :: uh(:,:,:,:,:)               ! Coefficients of the DG polynomial
    real, pointer  :: wh(:,:,:)                   ! Cell average values of primitive variables
    real, pointer  :: duh(:,:,:,:,:)              ! Update coefficients of the DG polynomial
    real, pointer  :: qhi(:,:,:,:,:)              ! Time-averaged coefficients of the space-time predictor
    real, pointer  :: Fhi(:,:,:,:,:,:)            ! Time-averaged coefficients of the flux tensor of the space-time predictor
    real, pointer  :: qbnd(:,:,:,:,:)             ! Boundary-extrapolated data for the state vector Q in the element
    real, pointer  :: Fbnd(:,:,:,:,:)             ! Boundary-extrapolated data for the normal flux F in the element
    real, pointer  :: Fu(:,:,:,:)                 ! Upwind flux in the x-direction
    real, pointer  :: Gu(:,:,:,:)                 ! Upwind flux in the y-direction

    ! Important info and parameters concerning the governing PDE system

    type tEquations
        real :: GAMMA
    end type tEquations

    type(tEquations)   :: EQN

end module ader_dg

