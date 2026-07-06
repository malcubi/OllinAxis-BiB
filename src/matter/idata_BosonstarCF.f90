
  subroutine idata_BosonstarCF

! *************************************************************
! ***   BOSON STAR INITIAL DATA IN CONFORMALLY FLAT GAUGE   ***
! *************************************************************

! Boson stars are solutions such that the spacetime is static
! and the complex scalar field has a harmonic dependence on time.
! This subroutine calculates initial data for a boson star
! in the conformally flat gauge using a shooting method.
!
! To obtain the initial data we assume that spacetime is
! static (K_ij=0), and also that the complex scalar field
! has the form:
!
! Phi(t,r) = phi(r) exp(i omega t)
!
! Notice that with this ansatz the stress-energy tensor
! is time-independent so the metric can be static.
!
! We also assume that the spatial metric is conformally flat:
!
!   2       4    2
! dl  =  psi   dl
!                flat
!
! We then need to solve the Hamiltonian constraint for the
! conformal factor "psi", the Klein-Gordon equation for
! the scalar field "phi", and the maximal slicing equation
! for the lapse "alpha".
!
! The Hamiltonian constraint takes the form:
!
! __2                    5
! \/    psi  +  2 pi psi rho  =  0
!   flat
!
! where we use the flat Laplacian, and where rho is the energy
! density of the complex scalar field given by:
!
!               /          2         2       4                         2       \
! rho  =  (1/2) | [ (d phi) + (d phi) ] / psi  +  ( omega phi / alpha )  + 2 V |
!               \     r         z                                              /
!
! and V(phi) is the self-interaction potential.
!
! On the other hand, the Klein-Gordon equation becomes:
!
! __2         ij                                    2         
! \/ phi  +  g   d phi d ln(alpha)  +  (omega/alpha) phi  -  V'  = 0 
!                 i     j  
!
! Finally, the maximal slicing equation takes the form:
!
! __2
! \/ alpha  -  4 pi (rho + trS)  =  0
!
!
! where now:
!                                        2
! rho + trS  =  2 [ ( omega phi / alpha )  -  V ]
!
! Notice that the equation for psi involves the flat Laplacian,
! while the eqautions for (phi,alpha) involve the physical
! Laplacian.  For a conformally flat metric the two Laplacian
! operators for a scalar function F are related as:
!
! __2       __2         4       ij    
! \/  F  =  \/   F / psi  +  2 g   d F  d ln(psi)
!            flat                   i    j
!
!              __2          ~ij                      4
!        =  [  \/   F  +  2 g  d F  d ln(psi) ] / psi
!               flat            i    j
!
! where the tilde indicates the flat metric in cylindrical
! coordinates, so that we simply have:
!
! ~ij
! g  d F d G  =  d F d G  +  d F d G
!     i   j       r   r       z   z
!
! The equations for phi and alpha then become:
!
!
! __2          ~ij                                              4               2
! \/   phi  +  g   [ d ln(alpha)  +  2 d ln(psi) ] d phi  +  psi [ (omega/alpha) phi  -  V' ] = 0
!  flat               i                 i           j
!
! __2              ~ij
! \/   alpha  +  2 g  d ln(psi) d alpha  -  4 pi (rho + trS)  =  0
!  flat                i         j
!
!
! This routine takes as input parameter the value of the scalar
! field at the origin "boson_phi0", and solves the coupled system
! of elliptic equations.
!
! IMPORTANT:  To avoid confusion, do notice that in fact we use the
! array "phi" for the conformal factor "psi" (this is because in the
! code "psi" is declared non-evolving).  On the other hand, the
! scalar field is called "complex_phi".

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo
  use derivatives

! Extra variables.

  implicit none

  integer box,level              ! Box number and level counters.
  integer i,j                    ! Counters.
  integer step                   ! Iteration counter.
  integer Nlmax_old              ! Original number of levels.
  real(8) NNtot

  real(8) lres,gres              ! Local and global residuals.
  real(8) waveeta                ! Damping parameter.
  real(8) cfac                   ! Courant parameter.
  real(8) one,half,smallpi       ! Numbers.

  character(3) method            ! Time integration method.

  integer s_ext(0:Nb,0:Nlmax)    ! External values of step counter.

  real(8) t_ext(0:Nb,0:Nlmax)    ! External values of time counter.
  real(8) t1_ext(0:Nb,0:Nlmax)   ! External values of t1 counter.
  real(8) t2_ext(0:Nb,0:Nlmax)   ! External values of t2 counter.


! *******************
! ***   NUMBERS   ***
! *******************

  one  = 1.d0
  half = 0.5d0

  smallpi = acos(-1.d0)


! ************************************************
! ***   SAVE EXTERNAL STEP AND TIME COUNTERS   ***
! ************************************************

! Since I just do calls to "onestep" here, I need to use the step
! and time counters.  But once we leave the routine they must be
! set back to their original external values, so I save those here.

  s_ext = s
  t_ext = t

  t1_ext = t1
  t2_ext = t2

! Now set time and time step counters to zero.

  s = 0
  t = 0.d0

  t1 = 0.d0
  t2 = 0.d0


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, 'Solving initial data for a boson star in conformally flat gauge using shooting method ...'
     print *
  end if

! Sanity check.

  if (spacetime=="minkowski") then
     print *, 'Boson star initial data is not compatible with a Minkowski background ...'
     print *, 'Aborting! (subroutine idata_BosonstarCF)'
     print *
     call die
  end if


! *************************
! ***   NORMALIZATION   ***
! *************************

! An additional factor of 1/sqrt(4*pi) can be included in
! the normalization (this is Olvier Sarbach's convention).

  if (boson_factor=="physical") then

     if (rank==0) then
        write(*,'(A)') ' Using "physical" normalization at origin: phi(r=0) = phi0'
        print *
     end if

  else

     boson_phi0 = boson_phi0/(4.d0*smallpi)

     if (rank==0) then
        write(*,'(A,ES23.16)') ' Using "harmonic" normalization at origin: phi(r=0) = phi0/sqrt(4pi) = ', boson_phi0
        print *
     end if

  end if


! ************************
! ***   SANITY CHECK   ***
! ************************

! Sanity check for potential.

  if ((complexpotential=="phi2").and.(complex_lambda/=0.d0)) then
     print *, 'You can not have a complexpotential=phi2 and complex_lambda different from zero.'
     print *, 'Aborting! (subroutine idata_BosonstarCF)'
     print *
     call die
  end if


! *******************************
! ***   INITIALIZE SOLUTION   ***
! *******************************

! Initialize boson_omega.

  boson_omega = 1.d0

! Lapse and conformal factor are initialized to 1.
! Notice that we use the array "phi" for the conformal
! factor "psi" (this is because in the code "psi" is
! declared non-evolving).
!
! The scalar field is initialized to a gaussian centered
! on the origin with the correct amplitude.
!
! The time derivatives are set initially to 0.

  do box=0,Nb
     do level=min(1,box),Nl(box)

        call currentgrid(box,level,grid(box,level))

        alpha = 1.d0
        phi   = 1.d0

!       The scalar field is initialized to a gaussian centered
!       on the origin with the correct amplitude.

        complex_phiR = boson_phi0*exp(-rr**4/(r**2*complexR_sr0**2 + z**2*complexR_sz0**2))

!       Set time derivatives to 0.

        dtalpha = 0.d0       ! Time derivative of alpha
        dtphi = 0.d0         ! Time derivatve of phi
        complex_piR = 0.d0   ! Time derivative of complex_phiR

     end do
  end do


! *****************************
! ***   SAVE INITIAL DATA   ***
! *****************************

! Save initial data to file if required.

  if (ELL_verbose) then

     do box=0,Nb
        do level=min(1,box),Nl(box)

           call currentgrid(box,level,grid(box,level))

           grabvar => alpha
           call save1Dvariable('alpha_boson',directory,box,level,outparallel,'replace')
           call save2Dvariable('alpha_boson',directory,box,level,outparallel,'replace')

           grabvar => phi
           call save1Dvariable('psi_boson',directory,box,level,outparallel,'replace')
           call save2Dvariable('psi_boson',directory,box,level,outparallel,'replace')

           grabvar => complex_phiR
           call save1Dvariable('phiR_boson',directory,box,level,outparallel,'replace')
           call save2Dvariable('phiR_boson',directory,box,level,outparallel,'replace')

        end do
     end do

  end if


! *****************************
! ***   COURANT PARAMETER   ***
! *****************************

! Internal time integration method.

  method = WE_method

! Courant parameter.

  cfac = WE_dtfac

! Time step.

  dt0 = cfac*min(dr0,dz0)

! Finer grids.

  do level=0,Nlmax
     dtl(level) = dt0/2**level
  end do


! ****************************
! ***   START ITERATIONS   ***
! ****************************

! Initialize damping parameter.

  waveeta = WE_eta

! Initialize residuals and iteration number.

  100 continue

  lres = 1.d0
  gres = 1.d0

  step = 0

! Start iterations.  Notice that we do at least 100 iterations.
! This is because since we start with ell_v=0, for the first few
! time steps the change in ell_u might be very small.

  do while ((step<100).or.((gres>ELL_epsilon).and.(step<ELL_maxiter)))


!    ******************************************
!    ***   ADVANCE ONE INTERNAL TIME STEP   ***
!    ******************************************

!    Increment step counter.

     step = step + 1

!    Advance one time step.

     call bosonstep(0,waveeta,method)


!    *************************
!    ***   FIND RESIDUAL   ***
!    *************************

!    Find local residual.

     lres = 0.d0
     gres = 0.d0

     do j=1,Nzl(0,rank)-ghost
        do i=1,Nrl(0,rank)-ghost
           lres = lres + (grid(0,0)%complex_phiR(i,j)-grid(0,0)%complex_phiR_p(i,j))**2
        end do
     end do

     NNtot = Nrl(0,rank)*Nzl(0,rank)
     lres = sqrt(lres/dble(NNtot))

!    Find global residual.

     if (size>1) then
        call MPI_Allreduce(lres,gres,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     else
        gres = lres
     end if


!    ******************
!    ***   OUTPUT   ***
!    ******************

!    Output if required.

     if ((ELL_verbose).and.(mod(step,ELL_Noutput)==0)) then

!       Data to screen.

        if (rank==0) then
           write(*,"(A,i5,A,ES15.8E2)") ' WaveElliptic:   Iteration = ',step,'         Residual = ',gres
        end if

!       Save data to file.

        do box=0,Nb
           do level=min(1,box),min(Nl(box),Nlmax)

              call currentgrid(box,level,grid(box,level))

              grabvar => alpha
              call save1Dvariable('alpha_boson',directory,box,level,outparallel,'old')
              call save2Dvariable('alpha_boson',directory,box,level,outparallel,'old')

              grabvar => phi
              call save1Dvariable('psi_boson',directory,box,level,outparallel,'old')
              call save2Dvariable('psi_boson',directory,box,level,outparallel,'old')

              grabvar => complex_phiR
              call save1Dvariable('phiR_boson',directory,box,level,outparallel,'old')
              call save2Dvariable('phiR_boson',directory,box,level,outparallel,'old')

           end do
        end do

     end if


! *************************************
! ***   END OF ITERATIONS DO LOOP   ***
! *************************************

  end do


! ****************************
! ***   DID WE CONVERGE?   ***
! ****************************

  if (rank==0) then

     if (step==ELL_maxiter) then
        write (*,'(A,i6,A)') ' BosonStarCF:   Iterations did not converge after ',ELL_maxiter,' iterations.'
        print *
     else
        if (Nlmax_old>0) then
           if (Nlmax_old/=Nlmax) then
              write (*,'(A,i5,A)') ' BosonStarCF:   Coarse grid solution converged after ',step,' iterations!'
              print *
           else
              write (*,'(A,i5,A)') ' BosonStarCF:   Finer grids solution converged after ',step,' iterations!'
              print *
           end if
        else
           write (*,'(A,i5,A)') ' BosonStarCF:   Solution converged after ',step,' iterations!'
           print *
        end if
     end if

  end if


! *****************************
! ***   REFINEMENT LEVELS   ***
! *****************************



! ************************************************************
! ***   SET BACK TO ZERO THE TIME AND TIME STEP COUNTERS   ***
! ************************************************************

! Since we advanced internally the time and time step counter arrays,
! we now need to set them back to their original external values.

  s = s_ext
  t = t_ext

  t1 = t1_ext
  t2 = t2_ext

! And fix the time step with the correct Courant factor.

  dt0 = dtfac*min(dr0,dz0)

  do level=0,Nlmax
     dtl(level) = dt0/2**level
  end do


! ********************************************
! ***   RECOVER CORRECT CONFORMAL FACTOR   ***
! ********************************************

! Remember that we where using "phi" instead of "psi",
! so copy the result back into "psi".

 do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

        psi = phi

!       Now calculate the correct "phi".

        phi  = dlog(psi)

!       Find "chi".

        if (chimethod) then
           chi  = 1.d0/psi**chipower
        end if

!       Find (psi2,psi4).

        psi2 = psi**2
        psi4 = psi**4

     end do
  end do


! **********************************************************
! ***   IMAGINARY PART OF (phi,xi) AND REAL PART OF pi   ***
! **********************************************************

 do box=0,Nb
     do level=min(1,box),Nl(box)

!       Find spatial derivatives of complex_phiR.

        diffvar => complex_phiR
        complex_xiR_r = diff1r(+1)
        complex_xiR_z = diff1z(+1)

!       Set time derivative of imaginary part to (omega/alpha)*phiR.

        complex_piI = boson_omega*complex_phiR/alpha

     end do
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine idata_BosonstarCF







  recursive subroutine bosonstep(level,waveeta,method)

! *********************************
! ***   ADVANCE ONE TIME STEP   ***
! *********************************

! Here we advance one time step in the internal evolution
! of the boson star solver.

! Include modules.

  use param
  use arrays
  use derivatives
  use procinfo

! Extra variables.

  implicit none

  logical firstcall
  logical icn,rk4

  integer box,level    ! Box number and level counters.
  integer i,j,k        ! Counters.
  integer iter         ! Counter for internal iterations.
  integer niter        ! Number of internal iterations.
  integer bmax         ! Number of boxes at this level.

  real(8) dtw          ! Internal time step.
  real(8) weight       ! Weight for rk4.
  real(8) waveeta      ! Damping parameter.

  character(*) method  ! Time integration method.

  data firstcall / .true. /

  save firstcall,icn,rk4


! *************************
! ***   LOGICAL FLAGS   ***
! *************************

  if (firstcall) then

      firstcall = .false.

      icn = (method=="icn")
      rk4 = (method=="rk4")

  end if


! ******************************************
! ***   NUMBER OF BOXES FOR THIS LEVEL   ***
! ******************************************

  if (level==0) then
     bmax = 0
  else
     bmax = Nb
  end if


! *****************************************
! ***   ITERATE OVER REFINEMENT BOXES   ***
! *****************************************

  do box=0,bmax

!    Check that the current box does indeed have
!    level l. If it doesnt't, go to next cycle.

     if (Nl(box)<level) cycle


!    *********************************
!    ***   POINT TO CURRENT GRID   ***
!    *********************************

     call currentgrid(box,level,grid(box,level))


!    ******************************
!    ***   SAVE OLD TIME STEP   ***
!    ******************************

     alpha_p = alpha
     phi_p = phi
     complex_phiR_p = complex_phiR

!    Save old values of internal boundaries.

     if (level>0) then

     end if


!    **********************************************
!    ***   FIND NUMBER OF INTERNAL ITERATIONS   ***
!    **********************************************

     if (rk4) then
        niter = 4
     else if (icn) then
        niter = icniter
     end if


!    **************************************
!    ***   ADVANCE ONE FULL TIME STEP   ***
!    **************************************

     do iter=1,niter


!       ************************
!       ***   FIND WEIGHTS   ***
!       ************************

!       Find out weights for each iteration for the
!       different time integration schemes.

!       Fourth order Runge-Kutta.

        if (rk4) then

!          In fourth order Runge-Kutta the first two iterations
!          jump half a time step and the last two a full time step.
!          Here we also set the weights with which intermediate
!          results contribute to final answer: 1/6 for first and
!          last intermediate results and 1/3 for the two middle ones.

           select case(iter)
              case(1)
                 dtw = 0.5d0*dt
                 weight = 1.d0/6.d0
              case(2)
                 dtw = 0.5d0*dt
                 weight = 1.d0/3.d0
              case(3) 
                 dtw = dt
                 weight = 1.d0/3.d0
              case(4)
                 dtw = dt
                 weight = 1.d0/6.d0
           end select

!       Iterative Crank-Nicholson (ICN).

        else if (icn) then

!          In ICN all iterations except the last one
!          jump only half a time step.

           if (iter<niter) then
              dtw = 0.5d0*dt
           else
              dtw = dt
           end if

        end if


!       ************************
!       ***   FIND SOURCES   ***
!       ************************

!       Derivatives of alpha.

        diffvar => alpha

        Dr_alpha = diff1r(+1)
        Dz_alpha = diff1z(+1)
     
        Drr_alpha = diff2r(+1)
        Dzz_alpha = diff2z(+1)
        Drz_alpha = diff2rz(+1,+1)

!       Derivatives of phi.

        diffvar => phi

        Dr_phi = diff1r(+1)
        Dz_phi = diff1z(+1)
     
        Drr_phi = diff2r(+1)
        Dzz_phi = diff2z(+1)
        Drz_phi = diff2rz(+1,+1)

!       Derivatives of complex_phiR.

        diffvar => complex_phiR

        Dr_complex_phiR = diff1r(+1)
        Dz_complex_phiR = diff1z(+1)
     
        Drr_complex_phiR = diff2r(+1)
        Dzz_complex_phiR = diff2z(+1)
        Drz_complex_phiR = diff2rz(+1,+1)

!       Scalar field potential.

!       Find frequency omega.

!       The sources for (alpha,phi,complex_phiR) are just (dtalpha,dtphi,complex_piR).

        salpha = dtalpha
        sphi = dtphi
        scomplex_phiR = complex_piR

!       The sources of (dtalpha,dtphi,complex_piR) all include the flat Laplacian.

        sdtalpha = Drr_alpha + Dzz_alpha + Dr_alpha/r
        sdtphi = Drr_phi + Dzz_phi + Dr_phi/r
        scomplex_piR = Drr_complex_phiR + Dzz_complex_phiR + Dr_complex_phiR/r

!       Extra sources for dtalpha.

!       Extra sources for dtphi.

!       Extra sources for complex_piR.

!       And add some dissipation to reduce high frequency noise.

        evolvevar => dtalpha
        sourcevar => sdtalpha
        call dissipation(+1,+1,WE_diss)

        evolvevar => dtphi
        sourcevar => sdtphi
        call dissipation(+1,+1,WE_diss)

        evolvevar => complex_piR
        sourcevar => scomplex_piR
        call dissipation(+1,+1,WE_diss)

!       Symmetries on axis.

        if (ownaxis) then
           do i=1,ghost
              sdtalpha(1-i,:) = sdtalpha(i,:)
              sdtphi(1-i,:) = sdtphi(i,:)
              scomplex_piR(1-i,:) = scomplex_piR(i,:)
           end do
        end if

!       Symmetries on equator.

        if (eqsym.and.ownequator) then
           do j=1,ghost
              sdtalpha(:,1-j) = sdtalpha(:,j)
              sdtphi(:,1-j) = sdtphi(:,j)
              scomplex_piR(:,1-j) = scomplex_piR(:,j)
           end do
        end if


!       *******************************
!       ***   BOUNDARY CONDITIONS   ***
!       *******************************

!       At the moment I use radiative boundaries that reduce
!       to a Robin boundary condition for static solutions.
!       Newmann and Diritchlet boundaries can be added later.


!       *****************************************************
!       ***   FOR RUNGE-KUTTA ADD TO ACCUMULATOR ARRAYS   ***
!       *****************************************************

        if (rk4) then
           if (iter==1) then

              alpha_a   = weight*salpha
              dtalpha_a = weight*sdtalpha

              phi_a   = weight*sphi
              dtphi_a = weight*sdtphi

              complex_phiR_a = weight*scomplex_phiR
              complex_piR_a  = weight*scomplex_piR

           else if (iter<niter) then

              alpha_a   = alpha_a + weight*salpha
              dtalpha_a = dtalpha_a + weight*sdtalpha

              phi_a   = phi_a + weight*sphi
              dtphi_a = dtphi_a + weight*sdtphi

              complex_phiR_a = complex_phiR_a + weight*scomplex_phiR
              complex_piR_a  = complex_piR_a + weight*scomplex_piR

           else

              salpha    = alpha_a + weight*salpha
              sdtalpha  = dtalpha_a + weight*sdtalpha

              sphi    = phi_a + weight*sphi
              sdtphi  = dtphi_a + weight*sdtphi

              scomplex_phiR = complex_phiR_a + weight*scomplex_phiR
              scomplex_piR  = complex_piR_a + weight*scomplex_piR

           end if
        end if


!       ****************************
!       ***   UPDATE VARIABLES   ***
!       ****************************

        alpha   = alpha_p + dtw*salpha
        dtalpha = dtalpha_p + dtw*sdtalpha

        phi   = phi_p + dtw*sphi
        dtphi = dtphi_p + dtw*sdtphi

        complex_phiR = complex_phiR_p + dtw*scomplex_phiR
        complex_piR  = complex_piR_p + dtw*scomplex_piR


!       *************************************************
!       ***   FOR FINE GRIDS INTERPOLATE BOUNDARIES   ***
!       *************************************************


!       *****************************************
!       ***   SYNCHRONIZE ACROSS PROCESSORS   ***
!       *****************************************


!       ***********************************
!       ***   END INTERNAL ITERATIONS   ***
!       ***********************************

     end do


!    ****************************************************
!    ***   ADVANCE LOCAL TIME AND TIME STEP COUNTER   ***
!    ****************************************************

!    Save old local time.

     t2(box,level) = t1(box,level)
     t1(box,level) = t (box,level)

!    Advance time step counter and local time.

     s(box,level) = s(box,level) + 1
     t(box,level) = t(box,level) + dt


!    ***********************************************
!    ***   END ITERATION OVER REFINEMENT BOXES   ***
!    ***********************************************

  end do


! **********************************
! ***   ARE THERE FINER GRIDS?   ***
! **********************************

! If there is a finer grid we need to advance it twice
! to catch up.  Notice that here I am calling the
! current subroutine "wavestep" recursively.

  if (level<Nlmax) then
     call bosonstep(level+1,waveeta,method)
     call bosonstep(level+1,waveeta,method)
  end if


! **********************
! ***   SYMMETRIES   ***
! **********************

! Apply symmetries again since we might have messed
! them up when we where on higher refinement levels
! and restricted data to the current level.


! *************************************
! ***   DO WE NEED TO SYNC BOXES?   ***
! *************************************


! ****************************************************
! ***   RESTRICT FINE GRID DATA INTO COARSE GRID   ***
! ****************************************************


! ***************
! ***   END   ***
! ***************

  end subroutine bosonstep
