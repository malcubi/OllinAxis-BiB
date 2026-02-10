
! ****************************************
! ***   PARAMETERS ARE DECLARED HERE   ***
! ****************************************

  module param

!               VERY IMPORTANT  (PLEASE READ)
!
! Declare only one parameter per line.  Even if FORTRAN allows
! one to declare several parameters separated by commas, this
! file will be processed by a PERL script that does NOT understand
! that.
!
! Use always the format:
!
!    type :: name = value
!
! The following variations are permitted:  Character-type parameters
! that are allowed to receive multiple values at the same time
! (separated by commas) should be declared as:
!
!    character :: name = value     ! multiple
!
! Also, a range can be defined for character-type parameters as:
!
!    character :: name = value     ! range = (value1,value2,...,valuen)
!
! The range is not compulsory.  If it is not there, any value
! is allowed (for example, directory names).
!
! REMEMBER:  All parameters must be initialized.  The initial
! value should basically correspond to the code NOT doing
! anything weird.  That is, initialize to Minkowski, static
! slicing, no shift, vacuum, and all special features turned off.


! **************************
! ***   PARAMETER FILE   ***
! **************************

! The parameter file is piped to the code at run time,
! so it is not really itself a "parameter", but it is
! convenient to declare it here since its name will
! be known to all other routines (particularly) the
! checkpoint routine.

  character(1000) :: parfile = ""


! ****************
! ***   GRID   ***
! ****************

! eqsym:            Do we have equatiorial symmetry?
! dr0,dz0:          Spatial intervals.
! dr,dz:            Local spatial intervals (should not be set in parameter file).
!
! Nrtotal:          Total number of grid points in r direction for coarse grid.
! Nztotal:          Total number of grid points in z direction for coarse grid.
! Nr:               Local number of grid points in r direction (should not be set in parameter file).
! Nz:               Local number of grid points in z direction (should not be set in parameter file).
! Nrmax:            Local maximum number of grid points in r direction (should not be set in parameter file).
! Nzmax:            Local maximum number of grid points in z direction (should not be set in parameter file).
!
! Nb:               Number of refinement boxes (we allow up to 3 at the moment).
!
! Nl#:              Number of refinement levels for box #.
! Nl(#):            Array containing the above numbers to be filled at run time.
!
! Nrbox#:           Number of grid points in r direction for box #.
! Nrbox(#):         Array containing the above numbers to be filled at run time.
! Nzbox#:           Number of grid points in z direction for box #.
! Nzbox(#):         Array containing the above numbers to be filled at run time.
!
! rbox#:            r position of box #.
! rbox(#):          Array containing the above numbers to be filled at run time.
! zbox#:            z position of box #.
! zbox(#):          Array containing the above numbers to be filled at run time.
!
! ghost:            Number of ghost zones.
! interporder:      Order of interpolation.
!
! ownaxis:          Does the local grid own the axis?  (should not be set in parameter file).
! ownequator:       Does the local grid own the equator?  (should not be set in parameter file).
! ownorigin:        Does the local grid own the origin?  (should not be set in parameter file).
!
! nproctot:         Total number of processors uses in the run.  It is not really a parameter
!                   and I use it only for checkpointing (should not be set in parameter file).
!
! Notice that the number of ghost zones should not be
! given in the parameter file since it is controlled
! by the parameter "order".

  logical :: eqsym = .false.

  logical :: ownaxis = .false.       ! Should not be set in parameter file!
  logical :: ownequator = .false.    ! Should not be set in parameter file!
  logical :: ownorigin = .false.     ! Should not be set in parameter file!

  real(8) :: dr0 = 1.d0
  real(8) :: dz0 = 1.d0
  real(8) :: dr  = 0.d0              ! Should not be set in parameter file!
  real(8) :: dz  = 0.d0              ! Should not be set in parameter file!

  integer :: Nrtotal = 10
  integer :: Nztotal = 10

  integer :: Nr = 0                  ! Should not be set in parameter file!
  integer :: Nz = 0                  ! Should not be set in parameter file!
  integer :: Nrmax = 0               ! Should not be set in parameter file!
  integer :: Nzmax = 0               ! Should not be set in parameter file!

  integer :: Nb = 0

  integer :: Nl0 = 0
  integer :: Nl1 = 0
  integer :: Nl2 = 0
  integer :: Nl3 = 0
  integer :: Nl(0:3) = 0             ! Should not be set in parameter file!
  integer :: Nlmax = 1               ! Should not be set in parameter file!

  integer :: Nrbox1 = 1
  integer :: Nrbox2 = 1
  integer :: Nrbox3 = 1
  integer :: Nzbox1 = 1
  integer :: Nzbox2 = 1
  integer :: Nzbox3 = 1
  integer :: Nrbox(0:3) = 1          ! Should not be set in parameter file!
  integer :: Nzbox(0:3) = 1          ! Should not be set in parameter file!

  real(8) :: rbox1 = 0.d0
  real(8) :: rbox2 = 0.d0
  real(8) :: rbox3 = 0.d0
  real(8) :: zbox1 = 0.d0
  real(8) :: zbox2 = 0.d0
  real(8) :: zbox3 = 0.d0
  real(8) :: rbox(0:3) = 0.d0        ! Should not be set in parameter file!
  real(8) :: zbox(0:3) = 0.d0        ! Should not be set in parameter file!

  integer :: ghost = 0
  integer :: nproctot = 1            ! Should not be set in parameter file!

  character(1000) :: interporder = "bicubic"  ! range = (bilinear,bicubic)


! **************************************
! ***   REGULARIZATION FINE TUNING   ***
! **************************************

! The regularization procedure is not unique and there are several
! terms that are ambiguous and consistent as long as the sum of the
! coefficients of two contributions is equal to one.
!
! Then notation here is from the original code from J.M. Torres.
! There were in fact more parameters, but several did nothing so
! I (Miguel A.) deleted them.

! 1) Terms that only appear for evolutions with shift. These two terms
!    have to do with whether you choose to write A (KTA) or H (KTH)
!    in terms of lambda (Alambda), and appear in the shift sources
!    for these quantities in the routine sources_geometry.f90.

  real(8) :: ft1 = 0.d0     ! Ambiguous term on shift sources for lambda.
  real(8) :: ft2 = 0.d0     ! Ambiguous term on shift sources for Alambda.

! 2) Terms that appear in the routine auxiliary_geometry.f90.
!
!    Of these three parameters, ft4 and ft5 don't really affect
!    the evolutions at all since they just rewrite expressions
!    that are algrebraically entirely equivalent, their values
!    are then irrelevant.
!
!    The parameter ft3 does affect evolutions  since it switches
!    terms with KTA to terms with KTH. Taking ft3=0.5 seems optimal.

  real(8) :: ft3 = 0.5d0    ! Ambiguous regularization of KT2_lambda.
  real(8) :: ft4 = 0.d0     ! Ambiguous term in D2cov_alpha_lambda.
  real(8) :: ft5 = 0.d0     ! Ambiguous term in RIC_lambda.

! 3) Terms that appear in the routines for calculating the conformal
!    Ricci tensor:  calc_conformalRicci*.f90.
!
!    Both these parameters affect evolutions.  ft6 switches a term with
!    Delta_r/r to a term with Dr_Delta_r.  Its optimal value sems t be 1.
!    ft7 switches a term with Dz_A to a term with Dz_H.  Its optimal
!    value seems to be 0.5

  real(8) :: ft6 = 1.d0     ! Ambiguous term in RIC_lambda (calc_conformalRicci2.f90).
  real(8) :: ft7 = 0.5d0    ! Ambiguous term in RIC_lambda (calc_conformalRicci3.f90).

! 4) Terms that are only important with non-zero angular momentum.

  real(8) :: ft8 = 0.5d0    ! Ambiguous term in RIC_C1 (calc_conformalRicci4.f90).


! *************************
! ***   TIME STEPPING   ***
! *************************

! time           Local time (should not be set in parameter file).
! dt0:           Time step for base (coarsest) grid.
! dt:            Local time step (should not be set in parameter file).
! adjuststep:    Do we ajust the time step using the CFL condition?
! dtfac:         Courant parameter (dtfac=dt/dr).
! Nt:            Total number of time steps.

  logical :: adjuststep = .false.

  real(8) :: time  = 0.d0
  real(8) :: dt0   = 0.d0
  real(8) :: dt    = 0.d0
  real(8) :: dtfac = 0.5d0

  integer :: Nt = 10


! ******************
! ***   OUTPUT   ***
! ******************

! directory:       Name of directory for output.
! checkpointfile:  Name of checkpoint directory.
!
! An important point to remember is that, because of the
! way Fortran works, these strings have to be defined with a
! fixed length.  This means that before using them to identify
! a directory, they must first be trimmed to its true size.

  character(100) :: directory = "output"
  character(100) :: checkpointfile = "checkpoint"

! Output parameters.
!
! checkpoint:    Do we output checkpoint files?
! checkpointinitial:  Do we checkpoint initial data?
!
! Ninfo:         How often do we output information to screen?
! Noutput0D:     How often do we do 0D output?
! Noutput1D:     How often do we do 1D output?
! Noutput2D:     How often do we do 2D output?
! Ncheckpoint    How often do we do a checkpoint?
!
! outvars0D:     Variables that need 0D output (a list separated by commas).
! outvars1D:     Variables that need 1D output (a list separated by commas).
! outvars2D:     Variables that need 2D output (a list separated by commas).
!
! checkvars:     Variables that need to be output for a checkpoint (should not be set in parameter file).
!                It is initialized with the coordinate arrays (r,z,rr).  The evolving arrays are later
!                added to the list at run time when they are allocated.
!
! nvars0D:       Number of variables with 0D output (should not be set in parameter file).
! nvars1D:       Number of variables with 1D output (should not be set in parameter file).
! nvars2D:       Number of variables with 2D output (should not be set in parameter file).
!
! outparallel:   Type of parallel output (singlefile,fileperproc).
! commenttype:   Type of comment lines on files (xgraph,gnuplot).

  logical :: checkpointinitial = .true.
  logical :: checkpoint = .false.

  integer :: Ninfo = 10
  integer :: Noutput0D = 10
  integer :: Noutput1D = 10
  integer :: Noutput2D = 10
  integer :: Ncheckpoint = 1000
  integer :: nvars0D = 1
  integer :: nvars1D = 1
  integer :: nvars2D = 1

  character(1000) :: outvars0D = "alpha"         ! multiple
  character(1000) :: outvars1D = "alpha"         ! multiple
  character(1000) :: outvars2D = "alpha"         ! multiple
  character(1000) :: checkvars = "r,z,rr"        ! multiple

  character(1000) :: outparallel = "fileperproc" ! range = (singlefile,fileperproc)
  character(1000) :: commenttype = "gnuplot"     ! range = (xgraph,gnuplot)


! *******************
! ***   SLICING   ***
! *******************

! slicing:       Type of slicing condition.
! ilapse:        Type of initial lapse:
!                   one:        initial lapse equal to 1.
!                   shell:      Gaussian perturbation forming a shell.
!                   torus:      Gaussian perturbation forming a  torus.
!                   isotropic:  Isotropic lapse for Schwarzschild solution.
!                   quiso_kerr: Quasi-isotropic lapse for Kerr.
!                   psiminusN:  Precollapsed lapse (1/psi**N).
!
! maximalevery:  Do a maximal solve every so many steps.
!
! gauge_f:       Coefficient for Bona-Masso type slicing condition (positive).
! gauge_eta:     Coefficient of the extra term in the Generalized Harmonic slicing condition.
!
! lapsepert:     Perturbation to initial lapse.
! lapse_a0:      Amplitude of perturbation.
! lapse_r0:      Center of perturbation in r.
! lapse_z0:      Center of perturbation in z.
! lapse_rs0:     Width of perturbation in r.
! lapse_rz0:     Width of perturbation in z.

  character(30) :: slicing = "harmonic" ! range = (static,harmonic,1+log,shockavoid,harmonic2,1+log2,shockavoid2,maximal)
  character(30) :: ilapse  = "one"      ! range = (one,shell,torus,isotropic,quiso_kerr,psiminus2,psiminus4,maximal)

  integer :: maximalevery = 0

  real(8) :: gauge_f = 1.d0
  real(8) :: gauge_eta = 0.0d0

  character(30) :: lapsepert = "none"   ! range = (none,gauss)
  real(8) :: lapse_a0 = 0.0
  real(8) :: lapse_r0 = 0.0
  real(8) :: lapse_z0 = 0.0
  real(8) :: lapse_sr0 = 1.0
  real(8) :: lapse_sz0 = 1.0


! *****************
! ***   SHIFT   ***
! *****************

! shift:         Type of shift condition.
! ishift:        Type of initial shift.
! shiftafter:    After what time do we turn on the shift evolution?
!
! driverD0:      Do we add the term Delta0 to driver equations?
! drivercsi:     Coefficient of Deltar for Gammadriver.
! drivereta:     Coefficient of damping term for Gammadriver.
!
! driveradv:     Do we include advection terms in driver conditions?
!
! shiftpert:     Perturbation to initial shift.
! shift_a0:      Amplitude of perturbation.
! shift_r0:      Center of perturbation in r.
! shift_z0:      Center of perturbation in z.
! shift_rs0:     Width of perturbation in r.
! shift_rz0:     Width of perturbation in z.

  character(30) :: shift = "none"   ! range = (none,zero,static,Gammadriver0,Gammadriver1,Gammadriver2,Gammadrivershock1,Gammadrivershock2)
  character(30) :: ishift = "zero"  ! range = (zero,shell,torus,quiso_kerr)

  logical :: driverD0 = .false.
  logical :: driveradv = .false.    ! Seems to work better with advection terms off.

  real(8) :: shiftafter = 0.0
  real(8) :: drivercsi = 0.75
  real(8) :: drivereta = 0.0

  character(30) :: shiftpert = "none"   ! range = (none,gauss)

  real(8) :: shift_a0 = 0.0
  real(8) :: shift_r0 = 0.0
  real(8) :: shift_z0 = 0.0
  real(8) :: shift_sr0 = 1.0
  real(8) :: shift_sz0 = 1.0


! ************************
! ***   INITIAL DATA   ***
! ************************

! idata:         Type of initial data.

  character(30) :: idata = "minkowski"  ! range = (checkpoint,minkowski,schwarzschild,kerr,BrillLindquist,BrillWave,scalarpulse,complexpulse,testgw)


! *********************
! ***   EVOLUTION   ***
! *********************

! spacetime:     Dynamic or background (static) spacetime.
!
! integrator:    Time integration method.
! order:         Order of spatial differencing.
! icniter:       Number of iterations for icn.
!
! formulation:   Formulation of evolution equations (bssn,z4c).
! bssnflavor:    Evolution of volume element in BSSN (Eulerian/Lagrangian).
! eta:           Multiple of momentum constraint in source of Delta_i.
!
! angmom:        Turn on angular momentum terms?
! chimethod:     Use chi=1/psi**n=exp(-n*phi) in evolutions instead of phi (for BH evolutions).
! chipower:      Power of 1/psi in the definition of chi.
!
! evolveH:       Do we evolve H, or do we get it from the determinant?
! evolveKTH:     Do we evolve KTH, or do we get it from the traceless condition?
!
! shiftadvect:   Do we calculate advective derivates on the shift as one-sided?
!
! nolambda:      Do not use regularization variables.
! noDelta_r:     Do not use r BSSN Delta variable and calculate it.
! noDelta_z:     Do not use z BSSN Delta variable and calculate it.
! noDelta_p:     Do not use p BSSN Delta variable and calculate it.
!
! kappa1:        Constraint dissipation coefficient for Z4c.
! kappa2:        Constraint dissipation coefficient for Z4c.

  character(30) :: spacetime   = "dynamic"     ! range = (dynamic,background)
  character(30) :: formulation = "bssn"        ! range = (bssn,z4c)
  character(30) :: bssnflavor  = "lagrangian"  ! range = (eulerian,lagrangian)
  character(30) :: integrator  = "icn"         ! range = (icn,rk4)
  character(30) :: order       = "two"         ! range = (two,four)

  logical :: angmom = .false.

  logical :: chimethod   = .true.  ! Seems better to do this by default.
  logical :: evolveH     = .true.
  logical :: evolveKTH   = .false. ! Seems better not to evolve KTH independently.

  logical :: shiftadvect = .true.

  logical :: nolambda    = .true.  ! No regularization by default.
  logical :: noDelta_r   = .false.
  logical :: noDelta_z   = .false.
  logical :: noDelta_p   = .false.

  integer :: icniter  = 3
  integer :: chipower = 2

  real(8) :: eta = 2.0

  real(8) :: kappa1 = 0.0
  real(8) :: kappa2 = 0.0

! Kreiss-Oliger dissipation coefficients.  It is better if they are independent.

  real(8) :: geodiss     = 0.01    ! Dissipation for geometric variables.
  real(8) :: scalardiss  = 0.01    ! Dissipation for scalar fields.
  real(8) :: complexdiss = 0.01    ! Dissipation for complex scalar fields.


! **********************
! ***   BOUNDARIES   ***
! **********************

! boundary:      Type of boundary condition.

  character(30) :: boundtype = "radiative"  ! range = (none,static,flat,radiative,constraint)


! ****************************
! ***   ELLIPTIC SOLVERS   ***
! ****************************

! ELL_solver:    The elliptic solver used to get initial data:
!                     none  =  no elliptic solver active
!                     wave  =  WaveElliptic
!                     sor   =  Succesive Overrelaxation (SOR)
!
! ELL_verbose    Do we output iterations?
! ELL_Noutput    How often do we do output for iterations?
! ELL_maxiter    Maximum number of iterations for elliptic solver.
! ELL_epsilon    Tolerance for elliptic solver.
!
! WE_eta         Damping coefficient for WaveElliptic.

  character(30) :: ELL_solver = "wave"  ! range = (none,wave,sor)

  logical :: ELL_verbose = .false.
  integer :: ELL_Noutput = 1
  integer :: ELL_maxiter = 100000
  real(8) :: ELL_epsilon = 1.d-5

  real(8) :: WE_eta = 1.d-2


! ***********************
! ***   BLACK HOLES   ***
! ***********************

! N_BH:          Number of black holes (1,2).
! BHnMass:       Mass of the nth Black Hole.
! BHnSpin:       Angular momentum parameter of the nth Black Hole.
! BHnZ:          On-axis position of the nth Black Hole.

  integer :: N_BH = 1

  real(8) :: BH1mass = 1.0
  real(8) :: BH1spin = 0.0
  real(8) :: BH1Z = 0.0

  real(8) :: BH2mass = 1.0
  real(8) :: BH2spin = 0.0
  real(8) :: BH2Z = 0.0

  logical :: duplicate_bhs = .false.


! ********************
! ***   HORIZONS   ***
! ********************

! ahfind:        Do we want to find an apparent horizon?
! ahmove:        Does the center of the horizons move?
! ahshootfirst:  Always shoot on first call.
!
! ahverbose:     Do we output horizon info to screen?
! ahsaveiter:    Save data of flow algorithm?
!
! ahfind_every:  How often do we look for horizons?
! ahmethod:      Method to find horizon: shooting or spectral (relaxation) algorithm.
!
! N_ah           How many horizons do we look for (1 to 3).
! z_ah*          Initial center on axis of apparent horizons.
!
! ahzinit:       Initial z position for shooting. If set to 0 (the default)
!                the finder starts at 80% of the coarsest grid size.
!
! ahrmax:        Largest radius allowed (as a fraction of coarsest grid size).
! ahafter:       Start looking for horizons after a given time.
!
! ahNmodes:      Maximum number of Fourier modes for spectral finder.
! ahmaxiter:     Maximum number of iterations for spectral fast flow algorithm.
! ahAA:          Parameter AA for spectral fast flow algorithm.
! ahBB:          Parameter BB for spectral fast flow algorithm.
! aheccentric:   Eccentricity for initial guess in spectral fast flow algorithm.
!
! ahepsilon      Tolerance for horizon finder.

  logical :: ahfind = .false.
  logical :: ahmove = .false.
  logical :: ahshootfirst = .false.
  logical :: ahverbose = .false.
  logical :: ahsaveiter = .false.

  integer :: ahfind_every = 1
  integer :: N_ah = 1
  integer :: ahNmodes = 10
  integer :: ahmaxiter = 1500

  real(8) :: z_ah1 = 0.d0
  real(8) :: z_ah2 = 0.d0
  real(8) :: ahzinit = 0.d0
  real(8) :: ahrmax  = 0.9d0
  real(8) :: ahafter = 0.d0
  real(8) :: ahepsilon_shoot = 1.d-8
  real(8) :: ahepsilon_spectral = 1.d-5
  real(8) :: ahAA = 0.1d0
  real(8) :: ahBB = 0.5d0
  real(8) :: aheccentric = 0.d0

  character(30) :: ahmethod = "shoot"  ! range = (shoot,spectral)


! ***********************
! ***   BRILL WAVES   ***
! ***********************

! Amplitude, center and width of Brill wave data.

  real(8) :: brill_a0 = 0.0
  real(8) :: brill_sr0 = 1.0
  real(8) :: brill_sz0 = 1.0

  character(30) :: brilltype= "Holz"  ! range=(Holz)


! ***************************
! ***   TEUKOSLKY WAVES   ***
! ***************************

! Amplitude of wave.

  real(8) :: teukolsky_a0 = 0.001d0


! ******************
! ***   MATTER   ***
! ******************

! IMPORTANT:  The different matter types must have very different names,
!             and one name must NOT be contained in another.  That is,
!             don't use names such as "scalar" and "scalar-ghost" for
!             different matter types.  The reason for this is that the
!             code just checks for a given string to be contained in
!             "mattertype", and in such a case it will be contained
!             in both cases.

! mattertype (type of matter):
!
!    vacuum  = No matter.
!    scalar  = Scalar Field.
!    complex = Complex scalar field.

  character(1000) :: mattertype = "vacuum"  ! multiple, range=(vacuum,scalar,complex)


! *****************************
! ***   REAL SCALAR FIELD   ***
! *****************************

! scalarpotential:   Type of scalar field potential.
! scalar_mass:       Scalar field mass parameter.
! scalar_lambda:     Coefficient of phi**4 term in potential.
!
! scalar_a0:         Amplitude of gaussian perturbation.
! scalar_r0:         Center of gaussian perturbation.
! scalar_z0:         Center of gaussian perturbation.
! scalar_sr0:        Width of gaussian perturbation.
! scalar_sz0:        Width of gaussian perturbation.
! scalarpulse_type   Type of scalar pulse (torus,shell).
!
! scalarmethod:      Finite difference method.

  character(30) :: scalarpotential  = "none"    ! range = (none,phi2,phi4)
  character(30) :: scalarpulse_type = "shell"   ! range = (shell,torus,Baumgarte)
  character(30) :: scalarmethod     = "second"  ! range = (first,second)

  real(8) :: scalar_mass = 0.0
  real(8) :: scalar_lambda = 0.0

  real(8) :: scalar_a0 = 0.0
  real(8) :: scalar_r0 = 0.0
  real(8) :: scalar_z0 = 0.0
  real(8) :: scalar_sr0 = 1.0
  real(8) :: scalar_sz0 = 1.0


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

! complexpotential:   Type of complex field potential.
! complex_mass:       Complex field mass parameter.
! complex_lambda:     Coefficient of phi**4 term in potential.
!
! complex*_a0:        Amplitude of gaussian perturbation.
! complex*_r0:        Center of gaussian perturbation.
! complex*_z0:        Center of gaussian perturbation.
! complex*_sr0:       Width of gaussian perturbation.
! complex*_sz0:       Width of gaussian perturbation.
! complexpulse_type   Type of complex pulse (torus,shell).
!
! complexmethod:      Finite difference method.

  character(30) :: complexpotential  = "none"    ! range = (none,phi2,phi4)
  character(30) :: complexpulse_type = "shell"   ! range = (shell,torus)
  character(30) :: complexmethod     = "second"  ! range = (first,second)

  real(8) :: complex_mass = 0.0
  real(8) :: complex_lambda = 0.0

  real(8) :: complexR_a0 = 0.0
  real(8) :: complexR_r0 = 0.0
  real(8) :: complexR_z0 = 0.0
  real(8) :: complexR_sr0 = 1.0
  real(8) :: complexR_sz0 = 1.0

  real(8) :: complexI_a0 = 0.0
  real(8) :: complexI_r0 = 0.0
  real(8) :: complexI_z0 = 0.0
  real(8) :: complexI_sr0 = 1.0
  real(8) :: complexI_sz0 = 1.0


! ********************************
! ***   CURVATURE INVARIANTS   ***
! ********************************

! curvInv: Calculate cuvature invariants I,J

  logical :: curvInv = .false.


! ********************************
! ***   GRAVITATIONAL WAVES    ***
! ********************************

! wave_extract: Gravitational waves extraction
! wavextract_every : How often do we extract gravitational waves.
! rad1        : Radius of extraction
! rad2        : Radius of extraction
! rad3        : Radius of extraction

  logical :: wave_extract = .false.
  integer :: wavextract_every = 5
  real(8) :: rad1 = 0.0
  real(8) :: rad2 = 0.0
  real(8) :: rad3 = 0.0


! *******************************************
! ***   CONVERT AXISYMMETRIC DATA TO 3D   ***
! *******************************************

! convert_to_3D : Transform to 3D.
! dx,dy       : Cartesian spatial intervals.
! Nx,Ny       : Cartesian grid size.
! sgrid       : Do we stagger the cartesian grid?

  real(8) :: dxx = 0.d0
  real(8) :: dyy = 0.d0
!  real(8) :: dzz = 0.d0


  integer :: Nxx = 100
  integer :: Nyy = 100
  integer :: Nzz = 100

  logical :: convert_to_3D = .false.
  logical :: sgrid = .false.


! ***************
! ***   END   ***
! ***************

  end module param

