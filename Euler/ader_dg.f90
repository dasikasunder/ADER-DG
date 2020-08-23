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

    integer, parameter          :: N = 2      ! Polynomial degree of our approximation in space and time
    integer, parameter          :: nDim = 2   ! The number of space dimensions that we actually want to simulate
    double precision, parameter :: CFL = 0.9  ! The Courant-Friedrichs-Lewy number < 1
    integer, parameter          :: nVar = 4   ! The number of variables of the PDE system

    ! ---------------------------------------------------------------------------------------------------------
    ! Do not change beyond this point
    ! ---------------------------------------------------------------------------------------------------------

    ! The following variables contain important information about the numerical method.

    double precision, parameter :: PNPMTable(0:9) = (/ 1.0, 0.33, 0.17, 0.1, 0.069, 0.045,  0.038, 0.03, 0.02, 0.015 /)   ! maximum CFL numbers for PNPM schemes
    integer                     :: nDOF(0:2)                           ! The number of degrees of freedom in space and time
    double precision            :: xiGPN(N+1), wGPN(N+1)               ! The Gauss-Legendre quadrature nodes and weights
    double precision            :: xin(N+1)                            ! The nodes used for our basis (can in principle be different from the Gauss-Legendre nodes, like for example also the bad Gauss-Lobatto nodes)
    double precision            :: MM(N+1,N+1), iMM(N+1,N+1)           ! Element mass-matrix and its inverse
    double precision            :: Kxi(N+1,N+1), dudx(N+1,N+1)         ! Element stiffness matrix and discrete derivative operator
    double precision            :: FLm(N+1,N+1), FLp(N+1,N+1)          ! Left flux matrices
    double precision            :: FRm(N+1,N+1), FRp(N+1,N+1)          ! Right flux matrices
    double precision            :: FLcoeff(N+1), FRcoeff(N+1)          ! extrapolation coefficients to the left and right boundary
    double precision            :: F0(N+1), F1(N+1,N+1)                ! Time flux matrices
    double precision            :: K1(N+1,N+1), iK1(N+1,N+1)           ! F1 - Ktau

    ! Stuff related to the problem setup, mesh and output

    integer                   :: IMAX, JMAX                     ! Number of cells in each space dimension
    integer                   :: nElem, nFace, nNode            ! Number of elements, faces and nodes
    integer                   :: timestep                       ! Number of the current time step
    integer                   :: WriteInterval                  ! Number of time steps after which to write output data
    double precision          :: xL(nDim), xR(nDim)             ! Computational domain
    double precision          :: dx(nDim), dt                   ! Vector of mesh spacings in each dimension and the time step
    double precision          :: time, tend                     ! Current time  and final time
    double precision, pointer :: x(:,:,:)                       ! Coordinates of cell centers
    integer, parameter        :: nVtx = 2**nDim, nFac = 2*nDim  ! Number of vertices and faces per element
    integer                   :: ICType                         ! Type of initial condition to be used
    integer                   :: bL, bR, bB, bT                     ! Boundary conditions on the left, right, bottom and top of the domain
    character(LEN=200)        :: BaseFile                       ! Basic filename to write the results

    ! Some diagnostics                                        !
    real               :: tCPU1, tCPU2                        ! CPU times
    integer(8)         :: TEU                                 ! total element updates

    ! Data needed for the subcell limiter

    integer, parameter         :: nSubLim = 2*N+1             ! Number of subcells
    integer                    :: nSubLimV(nDim)              ! Vector for number of subcells in each dimension
    double precision, pointer  :: uh2lim(:,:), lim2uh(:,:)    ! Mapping from DG polynomials to the subcells and back
    double precision, pointer  :: uh2lob(:,:)                 ! Mapping from the DG polynomials to the Gauss-Lobatto points
    double precision           :: xiLob(N+1), wLob(N+1)       ! Gauss lobatto points and weights
    integer, pointer           :: neighbor(:,:,:,:,:)         ! Set of Voronoi neighbors of a cell
    double precision, pointer  :: olduh(:,:,:,:,:)            ! for the a posteriori limiter, we need the possibility to go back
    integer, pointer           :: recompute(:,:)              ! Map containing the troubled zones that need to be recomputed
    double precision           :: xilimbary(nSubLim)          ! Barycenters of the subcells
    integer                    :: Face2Neigh(2,4)             ! Mapping from face to neighbor index
    integer                    :: dn(2)                       ! number of direct neighbors in each dimension

    type tLimiter
    integer                    :: status, oldstatus
    double precision, pointer  :: Lh(:,:,:)
    double precision, pointer  :: NewLh(:,:,:)
    double precision           :: lmin(nVar), lmax(nVar)
    double precision           :: lmin_old(nVar), lmax_old(nVar)
    end type tLimiter

    type(tLimiter), pointer :: Limiter(:,:)

    ! The main variables of the ADER-DG scheme

    double precision, pointer  :: uh(:,:,:,:,:)               ! Coefficients of the DG polynomial
    double precision, pointer  :: wh(:,:,:)                   ! Cell average values of primitive variables
    double precision, pointer  :: duh(:,:,:,:,:)              ! Update coefficients of the DG polynomial
    double precision, pointer  :: qhi(:,:,:,:,:)              ! Time-averaged coefficients of the space-time predictor
    double precision, pointer  :: Fhi(:,:,:,:,:,:)            ! Time-averaged coefficients of the flux tensor of the space-time predictor
    double precision, pointer  :: qbnd(:,:,:,:,:)             ! Boundary-extrapolated data for the state vector Q in the element
    double precision, pointer  :: Fbnd(:,:,:,:,:)             ! Boundary-extrapolated data for the normal flux F in the element
    double precision, pointer  :: Fu(:,:,:,:)                 ! Upwind flux in the x-direction
    double precision, pointer  :: Gu(:,:,:,:)                 ! Upwind flux in the y-direction

    ! Important info and parameters concerning the governing PDE system

    type tEquations
        double precision :: GAMMA
    end type tEquations

    type(tEquations)   :: EQN

end module ader_dg

