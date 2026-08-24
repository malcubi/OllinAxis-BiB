
  subroutine idata_RotBoson

! ********************************************
! ***   ROTATING BOSON STAR INITIAL DATA   ***
! ********************************************

! Boson stars are solutions such that the spacetime is static and the complex
! scalar field has a harmonic dependence on time. This subroutine calculates
! initial data for a rotating boson star.
!
! We assume the spatial metric has the form:
!
!   2       4     2      2     2        2
! dl  =  psi  ( dr  +  dz  +  r  H  dphi )
!
! so we have two non-trivial metric cofficientes: (psi,H).
!
! In addition, we need the lapse function "alpha", and since the star is
! rotating we also need a non-zero angular shift component "beta_p" (index up).
!
! We assume that the space-time is static son that all metric functions are time
! independent.  Since there is a non-trivial shift vector this means that no all
! the extrinsic curvature components vanish.  We find:
!
!               4  2
! K_rp  =  ( psi  r  H / 2 alpha ) d beta_p
!                                   r
!               4  2
! K_zp  =  ( psi  r  H / 2 alpha ) d beta_p
!                                   z
!
! All other components of K (and its trace) vanish. These expressions imply:
!
!  mn           2            2                 2                2
! K   K   =  ( r  H / 2 alpha  ) [ ( d beta_p )  +  ( d beta_p )  ]
!      mn                             r                z
!
! For the complex scalar field we use the ansatz:
!
! Phi(t,r) = phi(r) exp[ - i (omega t - L phi) ]
!
! where L is an integer that corresponds to the angular momentum "quantum number".
! This ansatz guarantees that the stress-energy tensor is time-independentand
! and also independent of the angle, so that the solution will be axi-symmetric.
!
! For the case when L is not zero, the amplitude of the scalar field phi(r) has
! to behave as r**L close to the origin, so in practice we solve for a function
! F(r) such that:
!
!           L
! phi(r) = r  F(r)
!
!
! the parameter "boson_phi0" the will then correspond to the value at r=0 of F(r)
! (and not of phi(r)).
!
! The equation for alpha takes the final form:
!
!      2            2         2
!    d alpha + d alpha + d alpha / r + 2 [ d alpha d ln(psi) + d alpha d ln(psi) ]  
!     r         z         r                 r       r           z       z
!
!
!       + (1/2) [ d alpha d ln(H) + d alpha d ln(H) ]
!                 r       r         z      z
!
!           2      4
!       - (r  H psi / 2 alpha ) [ d alpha d beta_p + d alpha d beta_p ]
!                                  r       r          z       z
!                                   
!                        4                                      2
!       -  8 pi alpha psi  [ ( (omega + L beta_p ) phi / alpha )  -  V ]  = 0
!
!
! where V is the scalar field potential.
!
! The equation for beta_p takes the final form:
!
!
!
! The equation for psi takes the final form:
!
!     2       2
!    d psi + d psi 
!     r       z
!
! In the previous expressions, the physical Laplacian of a function "f" is given by:
!
!
! To solve the system we use the maximal slicing condition for the lapse, the angular
! momentum constraint for the shift, a combination of the Hamiltonian constraint and
! the equation dK_pp/dt=0 for psi and H, and finally the Klein-Gordon equation for
! the scalar field "F". (do not confuse the scalar field with the azimuthal angle).
!
! IMPORTANT:  To avoid confusion, do notice that in fact we use the array "phi" for
! the conformal factor "psi" (this is because "psi" is declared non-evolving). 
! On the other hand, for the scalar field function F(r) we use the array "complex_phiR",
! and at the end rescale with r**L.

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
  integer NNloc,NNtot            ! Total number of grid points.
  integer miniter                ! Minimum number of iterations.
  integer :: maxiter = 100000    ! Maximum number of iterations.

  real(8) lres,gres              ! Local and global residuals.
  real(8) cfac                   ! Courant parameter.
  real(8) one,half,smallpi       ! Numbers.

  character(3) method            ! Time integration method.

! Local time stepping arrays.

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
     print *, 'Solving initial data for a rotating boson star ...'
     print *
  end if

! Sanity check.

  if (spacetime=="minkowski") then
     print *, 'Rotating boson star initial data is not compatible with a Minkowski background ...'
     print *, 'Aborting! (subroutine idata_RotBoson)'
     print *
     call die
  end if


! *************************
! ***   NORMALIZATION   ***
! *************************

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


! *******************************
! ***   INITIALIZE SOLUTION   ***
! *******************************

! Initialize boson_omega.

  boson_omega = 1.d0

! The metric functions (alpha,psi,H) are initialized to 1.
! Notice that we use the array "phi" for the conformal
! factor "psi" (this is because in the code "psi" is
! declared non-evolving).
!
! The angular shift vector is initialized to 0.
!
! The time derivatives are set initially to 0.
! Notice that I use KTH as the time derivative
! of H.

  do box=0,Nb
     do level=min(1,box),Nl(box)

        call currentgrid(box,level,grid(box,level))

        alpha = one
        phi = one
        H = one

        beta_p = 0.d0

!       The scalar field is initialized to a gaussian centered
!       on the origin with the correct amplitude.  Remember
!       that this is in fact the function F(r) such that
!       phi(r) = r^L F(t).

        complex_phiR = boson_phi0*exp(-rr**2/boson_L**2)

!       Set all time derivatives to 0.

        dtalpha = 0.d0       ! Time derivative of alpha.
        dtphi = 0.d0         ! Time derivatve of phi.
        KTH = 0.d0           ! Time derivative of H.
        dtbeta_p = 0.d0      ! Time derivative of beta_p.
        complex_piR = 0.d0   ! Time derivative of complex_phiR.

     end do
  end do


! *****************************
! ***   SAVE INITIAL DATA   ***
! *****************************

! Save initial data to file if required.

  if (WE_verbose) then

     do box=0,Nb
        do level=min(1,box),Nl(box)

           call currentgrid(box,level,grid(box,level))

           grabvar => alpha
           call save1Dvariable('rotboson_alpha',directory,box,level,outparallel,'new')
           !call save2Dvariable('rotboson_alpha',directory,box,level,outparallel,'new')

           grabvar => phi
           call save1Dvariable('rotboson_psi',directory,box,level,outparallel,'new')
           !call save2Dvariable('rotboson_psi',directory,box,level,outparallel,'new')

           grabvar => H
           call save1Dvariable('rotboson_H',directory,box,level,outparallel,'new')
           !call save2Dvariable('rotboson_H',directory,box,level,outparallel,'new')

           grabvar => beta_p
           call save1Dvariable('rotboson_beta_p',directory,box,level,outparallel,'new')
           !call save2Dvariable('boson_beta_p',directory,box,level,outparallel,'new')

           grabvar => complex_phiR
           call save1Dvariable('rotboson_F',directory,box,level,outparallel,'new')
           !call save2Dvariable('rotboson_F',directory,box,level,outparallel,'new')

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

! If we have refinement boxes, we first solve on the coarse
! grid, and then we use this as initial guess for the full
! solution. This should speed up the solver considerably.

  Nlmax_old = Nlmax

  if (Nlmax>0) then
     Nlmax = 0
  end if

! Initialize residual and iteration number.

  100 continue

  lres = 1.d0
  gres = 1.d0

  step = 0

! Start iterations.

  miniter = 100

  do while ((step<miniter).or.((gres>WE_epsilon).and.(step<WE_maxiter)))

!    ******************************************
!    ***   ADVANCE ONE INTERNAL TIME STEP   ***
!    ******************************************

!    Increment step counter.

     step = step + 1

!    Advance one time step.

     call rotbosonstep(0,method)


!    *************************
!    ***   FIND RESIDUAL   ***
!    *************************

!    Find local residual.

     lres = 0.d0
     gres = 0.d0

     NNloc = 0
     NNtot = 0

     do j=1,Nzl(0,rank)-ghost
        do i=1,Nrl(0,rank)-ghost
           NNloc = NNloc + 1
           lres = lres + abs(grid(0,0)%scomplex_piR(i,j)) &
                + abs(grid(0,0)%sdtalpha(i,j)) + abs(grid(0,0)%sdtphi(i,j)) &
                + abs(grid(0,0)%sKTH(i,j)) + abs(grid(0,0)%sdtbeta_p(i,j))
        end do
     end do

!    Find global residual.

     if (size>1) then
        call MPI_Allreduce(lres,gres,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(NNloc,NNtot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierr)
     else
        gres = lres
        NNtot = NNloc
     end if

     gres = gres/dble(NNtot)


!    ******************
!    ***   OUTPUT   ***
!    ******************

!    Output if required.

     if (WE_verbose.and.(mod(step,WE_Noutput)==0)) then

!       Data to screen.

        if (rank==0) then

           write(*,"(A,i5,A,ES15.8E2,A,ES15.8E2)") ' Iteration = ',step, &
                   '     omega = ',boson_omega,'     Residual = ',gres

           !interpvar => grid(0,0)%complex_phiR
           !aux1 = interp(0,0,0.d0,0.d0,flag1)
           !print *, 'Value of phiR at origin = ',aux1

        end if

!       Save data to file.

        do box=0,Nb
           do level=min(1,box),min(Nl(box),Nlmax)

              call currentgrid(box,level,grid(box,level))

              grabvar => alpha
              call save1Dvariable('rotboson_alpha',directory,box,level,outparallel,'old')
              !call save2Dvariable('rotboson_alpha',directory,box,level,outparallel,'old')

              grabvar => phi
              call save1Dvariable('rotboson_psi',directory,box,level,outparallel,'old')
              !call save2Dvariable('rotboson_psi',directory,box,level,outparallel,'old')

              grabvar => H
              call save1Dvariable('rotboson_H',directory,box,level,outparallel,'old')
              !call save2Dvariable('rotboson_H',directory,box,level,outparallel,'old')

              grabvar => beta_p
              call save1Dvariable('rotboson_beta_p',directory,box,level,outparallel,'old')
              !call save2Dvariable('boson_beta_p',directory,box,level,outparallel,'old')

              grabvar => complex_phiR
              call save1Dvariable('rotboson_F',directory,box,level,outparallel,'old')
              !call save2Dvariable('rotboson_F',directory,box,level,outparallel,'old')

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

     if (step==WE_maxiter) then

        write (*,'(A,i6,A)') ' RotBoson: Iterations did not converge after ',WE_maxiter,' iterations.'
        print *

     else

        if (Nlmax_old>0) then

           if (Nlmax_old/=Nlmax) then
              write (*,'(A,i5,A)') ' RotBoson: Coarse grid solution converged after ',step,' iterations!'
              print *
           else
              write (*,'(A,i5,A)') ' RotBoson: Finer grids solution converged after ',step,' iterations!'
              print *
              write(*,'(A,ES23.16)') ' Final residual = ',gres
              write(*,'(A,ES23.16)') ' Omega          = ',boson_omega
              print *
           end if

        else

           write (*,'(A,i5,A)') ' RotBoson:   Solution converged after ',step,' iterations!'
           print *
           write(*,'(A,ES23.16)') ' Final residual = ',gres
           write(*,'(A,ES23.16)') ' Omega          = ',boson_omega
           print *

        end if

     end if

  end if


! *****************************
! ***   REFINEMENT LEVELS   ***
! *****************************

! Not yet implemented.


! *****************************
! ***   TOLMAN-KOMAR MASS   ***
! *****************************

! Not yet implemented.


! ***********************************
! ***   BOSON STAR PERTURBATION   ***
! ***********************************

! Perturbation not yet implemented.


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

!       Now calculate the correct "phi" and its derivatives.

        phi = dlog(psi)

        diffvar => phi

        Dr_phi = diff1r(+1)
        Dz_phi = diff1z(+1)

!       Find "chi".

        if (chimethod) then
           chi  = 1.d0/psi**chipower
        end if

!       Find (psi2,psi4).

        psi2 = psi**2
        psi4 = psi**4

     end do
  end do


! *****************************************************
! ***   RECOVER CORRECT SCALAR FIELD PHI = R**L F   ***
! ***   AND IMAGINARY PART OF TIME DERIVATIVE piI   ***
! ***            ALSO SET phiI=piR=0                ***
! *****************************************************

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Put back the factor r**L on the scalar field.

        complex_phiR = r**boson_L*complex_phiR

!       Find spatial derivatives of complex_phiR.

        diffvar => complex_phiR
        complex_xiR_r = diff1r(+1)
        complex_xiR_z = diff1z(+1)

!       From the original ansatz, we take the time derivative
!       of the imaginary part equal to:
!
!       piI = Im [ d(boson_phi)/dt - beta_p d(boson_phi)/dphi ]
!
!           = - (omega + L*beta_p)*phi/alpha

        complex_piI = - (boson_omega + boson_L*beta_p)*complex_phiR/alpha

!       Set the imaginary part of the scalar field and its spatial derivatives to zero.

        complex_phiI  = 0.d0
        complex_xiI_r = 0.d0
        complex_xiI_z = 0.d0

!       Set time derivative of the real part of phi to zero.

        complex_piR  = 0.d0

     end do
  end do


! ***************************************************************************
! ***   SET TIME DERIVATIVES BACK TO ZERO, AND FIND EXTRINSIC CURVATURE   ***
! ***************************************************************************

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Set again (dtalpha,dtphi,dtbeta_p) to 0.

        dtalpha = 0.d0
        dtphi = 0.d0
        dtbeta_p = 0.d0

!       Find the correct extrinsic curvature.  Remember that
!       for this initial data we have trK=0, and the only
!       non-zero components of the conformal extrinsic
!       curvature are:
!
!       KTC1 = H/(2*alpha) Dr_beta_p / r
!       KTC2 = H/(2*alpha) Dz_beta_p

        KTA = 0.d0
        KTB = 0.d0
        KTH = 0.d0
        KTC = 0.d0

        diffvar => beta_p
        Dr_beta_p  = diff1r(+1)
        Dz_beta_p  = diff1z(+1)

        KTC1 = 0.5d0*(H/alpha)*Dr_beta_p/r
        KTC2 = 0.5d0*(H/alpha)*Dz_beta_p

     end do
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine idata_RotBoson















  recursive subroutine rotbosonstep(level,method)

! *********************************
! ***   ADVANCE ONE TIME STEP   ***
! *********************************

! Here we advance one time step in the internal evolution
! of the rotating boson star solver.

! Include modules.

  use mpi
  use param
  use arrays
  use derivatives
  use procinfo

! Extra variables.

  implicit none

  logical firstcall
  logical icn,rk4

  integer box,level    ! Box number and level counters.
  integer i,j          ! Counters.
  integer iter         ! Counter for internal iterations.
  integer niter        ! Number of internal iterations.
  integer bmax         ! Number of boxes at this level.
  integer bbox

  real(8) dtw          ! Internal time step.
  real(8) weight       ! Weight for rk4.
  real(8) half,smallpi ! Numbers.

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


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0

  smallpi = acos(-1.d0)


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

!    Old values of lapse and time derivative.

     alpha_p   = alpha
     dtalpha_p = dtalpha

!    Old values of conformal factor and time derivative.

     phi_p   = phi
     dtphi_p = dtphi

!    Old values of H and time derivative.

     H_p   = H
     KTH_p = KTH

!    Old values of shift and time derivative.

     beta_p_p = beta_p
     dtbeta_p_p = dtbeta_p

!    Old values of complex_phiR and time derivative.

     complex_phiR_p = complex_phiR
     complex_piR_p  = complex_piR

!    Save old values of internal boundaries.

     if (level>0) then

!       THIS IS NOT YET IMPLEMENTED.

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

!       Derivatives of (alpha,dtalpha).

        diffvar => alpha

        Dr_alpha = diff1r(+1)
        Dz_alpha = diff1z(+1)
     
        Drr_alpha = diff2r(+1)
        Dzz_alpha = diff2z(+1)

        diffvar => dtalpha

        Dr_dtalpha = diff1r(+1)
        Dz_dtalpha = diff1z(+1)

!       Derivatives of (phi,dtphi).

        diffvar => phi

        Dr_phi = diff1r(+1)
        Dz_phi = diff1z(+1)
     
        Drr_phi = diff2r(+1)
        Dzz_phi = diff2z(+1)

        diffvar => dtphi

        Dr_dtphi = diff1r(+1)
        Dz_dtphi = diff1z(+1)

!       Derivatives of (H,KTH).

        diffvar => H

        Dr_H = diff1r(+1)
        Dz_H = diff1z(+1)
     
        Drr_H = diff2r(+1)
        Dzz_H = diff2z(+1)

        diffvar => KTH

        Dr_KTH = diff1r(+1)
        Dz_KTH = diff1z(+1)

!       Derivatives of (beta_p,dtbeta_p).

        diffvar => beta_p

        Dr_beta_p = diff1r(+1)
        Dz_beta_p = diff1z(+1)
     
        Drr_beta_p = diff2r(+1)
        Dzz_beta_p = diff2z(+1)

        diffvar => dtbeta_p

        Dr_dtbeta_p = diff1r(+1)
        Dz_dtbeta_p = diff1z(+1)

!       Derivatives of (complex_phiR,complex_piR).

        diffvar => complex_phiR

        Dr_complex_phiR = diff1r(+1)
        Dz_complex_phiR = diff1z(+1)
     
        Drr_complex_phiR = diff2r(+1)
        Dzz_complex_phiR = diff2z(+1)

        diffvar => complex_piR

        Dr_complex_piR = diff1r(+1)
        Dz_complex_piR = diff1z(+1)

!       Scalar field potential (at the moment I assume a pure mass term).
!       Remember that we need to include the factor of r**L.

        complex_V = half*complex_mass**2*(complex_phiR*r**boson_L)**2

!       Find frequency omega.  We solve for omega from the
!       Klein-Gordon equation at the first grid point.

        if (level==Nlmax) then

           if (ownaxis.and.ownequator) then

!             NOT YET IMPLEMENTED.

           end if

        end if

        call MPI_BCAST(boson_omega,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!       The sources for (alpha,phi,H,BETA_Pcomplex_phiR) are just (dtalpha,dtphi,KTH,DTBETA_P,complex_piR).

        salpha = dtalpha

        sphi = dtphi

        sH = KTH

        sbeta_p = dtbeta_p

        scomplex_phiR = complex_piR

        psi4 = phi**4

!       Sources for dtalpha.

         sdtalpha = Drr_alpha + Dzz_alpha + Dr_alpha/r &
                  + 2.d0*(Dr_alpha*Dr_phi + Dz_alpha*Dz_phi)/phi &
                  + half*(Dr_alpha*Dr_H + Dz_alpha*Dz_H)/H &
                  - half*(Dr_beta_p**2 + Dz_beta_p**2)*H*psi4*r**2/alpha &
                  - 8.d0*smallpi*alpha*psi4*(((boson_omega + boson_L*beta_p) &
                  *(r**boson_L*complex_phiR)/alpha)**2 - complex_V)

!       Sources for dtphi.

        sdtphi = 0.d0

!       Sources for KTH.

        sKTH = 0.d0

!       Sources for dtbeta_p.

        sdtbeta_p = Drr_beta_p + Dzz_beta_p + 3.d0*Dr_beta_p/r &
                  + 6.d0*(Dr_beta_p*Dr_phi + Dz_beta_p*Dz_phi)/phi &
                  + 1.5d0*(Dr_beta_p*Dr_H + Dz_beta_p*Dz_H)/H &
                  - (Dr_beta_p*Dr_alpha + Dz_beta_p*Dz_alpha)/alpha &
                  - 16.d0*smallpi*boson_L*(r**boson_L*complex_phiR)**2 &
                  *(boson_omega + boson_L*beta_p)/(r**2*H)

!       Sources for complex_piR.

        scomplex_piR = 0.d0

!       Damping term.  This is needed in order to avoid large
!       oscillations and make the iterations stable.  But it
!       seems we only need it for the Klein-Gordon equation.

        !scomplex_piR = scomplex_piR - complex_piR

!       And add some dissipation to reduce high frequency noise.

        evolvevar => dtalpha
        sourcevar => sdtalpha
        call dissipation(+1,+1,WE_diss)

        evolvevar => dtphi
        sourcevar => sdtphi
        call dissipation(+1,+1,WE_diss)

        evolvevar => KTH
        sourcevar => sKTH
        call dissipation(+1,+1,WE_diss)

        !evolvevar => dtbeta_p
        !sourcevar => sdtbeta_p
        !call dissipation(+1,+1,WE_diss)

        evolvevar => complex_piR
        sourcevar => scomplex_piR
        call dissipation(+1,+1,WE_diss)

!       But set the source at point (1,1) such that the
!       value at r=0 does not change.


!       Symmetries on axis.

        if (ownaxis) then
           do i=1,ghost
              sdtalpha(1-i,:) = sdtalpha(i,:)
              sdtphi(1-i,:) = sdtphi(i,:)
              sKTH(1-i,:) = sKTH(i,:)
              sdtbeta_p(1-i,:) = sdtbeta_p(i,:)
              scomplex_piR(1-i,:) = scomplex_piR(i,:)
           end do
        end if

!       Symmetries on equator.

        if (eqsym.and.ownequator) then
           do j=1,ghost
              sdtalpha(:,1-j) = sdtalpha(:,j)
              sdtphi(:,1-j)   = sdtphi(:,j)
              sKTH(:,1-j) = sKTH(:,j)
              sdtbeta_p(:,1-j) = sdtbeta_p(:,j)
              scomplex_piR(:,1-j) = scomplex_piR(:,j)
           end do
        end if


!       *******************************
!       ***   BOUNDARY CONDITIONS   ***
!       *******************************

!       Simple radiative boundaries for all three equations.

        if (level==0) then

!          Radiative conditions at r boundary.

           if (mod(rank+1,nprocr)==0) then
              i = Nr
              sdtalpha(i,:)     = - (r(i,:)*Dr_dtalpha(i,:)     + z(i,:)*Dz_dtalpha(i,:)     + dtalpha(i,:))/rr(i,:)
              sdtphi(i,:)       = - (r(i,:)*Dr_dtphi(i,:)       + z(i,:)*Dz_dtphi(i,:)       + dtphi(i,:))/rr(i,:)
              sKTH(i,:)         = - (r(i,:)*Dr_KTH(i,:)         + z(i,:)*Dz_KTH(i,:)         + KTH(i,:))/rr(i,:)
              sdtbeta_p(i,:)    = - (r(i,:)*Dr_dtbeta_p(i,:)    + z(i,:)*Dz_dtbeta_p(i,:)    + dtbeta_p(i,:))/rr(i,:)
              scomplex_piR(i,:) = - (r(i,:)*Dr_complex_piR(i,:) + z(i,:)*Dz_complex_piR(i,:) + complex_piR(i,:))/rr(i,:)
           end if

!          Radiative conditions at z boundaries.

           if (rank>=size-nprocr) then
              j = Nz
              sdtalpha(:,j)     = - (r(:,j)*Dr_dtalpha(:,j)     + z(:,j)*Dz_dtalpha(:,j)     + dtalpha(:,j))/rr(:,j)
              sdtphi(:,j)       = - (r(:,j)*Dr_dtphi(:,j)       + z(:,j)*Dz_dtphi(:,j)       + dtphi(:,j  ))/rr(:,j)
              sKTH(:,j)         = - (r(:,j)*Dr_KTH(:,j)         + z(:,j)*Dz_KTH(:,j)         + KTH(:,j  ))/rr(:,j)
              sdtbeta_p(:,j)    = - (r(:,j)*Dr_dtbeta_p(:,j)    + z(:,j)*Dz_dtbeta_p(:,j)    + dtbeta_p(:,j  ))/rr(:,j)
              scomplex_piR(:,j) = - (r(:,j)*Dr_complex_piR(:,j) + z(:,j)*Dz_complex_piR(:,j) + complex_piR(:,j))/rr(:,j)
           end if

        end if


!       *****************************************************
!       ***   FOR RUNGE-KUTTA ADD TO ACCUMULATOR ARRAYS   ***
!       *****************************************************

        if (rk4) then

        end if


!       ****************************
!       ***   UPDATE VARIABLES   ***
!       ****************************

        alpha   = alpha_p   + dtw*salpha
        dtalpha = dtalpha_p + dtw*sdtalpha

        phi   = phi_p   + dtw*sphi
        dtphi = dtphi_p + dtw*sdtphi

        H   = H_p   + dtw*sH
        KTH = KTH_p + dtw*sKTH

        beta_p   = beta_p_p   + dtw*sbeta_p
        dtbeta_p = dtbeta_p_p + dtw*sdtbeta_p

        complex_phiR = complex_phiR_p + dtw*scomplex_phiR
        complex_piR  = complex_piR_p  + dtw*scomplex_piR


!       *************************************************
!       ***   FOR FINE GRIDS INTERPOLATE BOUNDARIES   ***
!       *************************************************

!       For fine grids we need to interpolate from the new
!       time level of the coarse grid to get boundary data.
!
!       Remember that the coarse grid has already advanced
!       to the next time level.

        if (level>0) then

        end if


!       **********************
!       ***   SYMMETRIES   ***
!       **********************

!       Do we own axis and/or equator?

        ownaxis = (axis(box,rank)/=-1)
        ownequator = (eqz(box,rank)/=-1)

!       Symmetries on axis.

        if (ownaxis) then
           do i=1,ghost

              alpha(1-i,:)   = alpha(i,:)
              dtalpha(1-i,:) = dtalpha(i,:)

              phi(1-i,:)   = phi(i,:)
              dtphi(1-i,:) = dtphi(i,:)

              H(1-i,:)   = H(i,:)
              KTH(1-i,:) = KTH(i,:)

              beta_p(1-i,:)   = beta_p(i,:)
              dtbeta_p(1-i,:) = dtbeta_p(i,:)

              complex_phiR(1-i,:) = complex_phiR(i,:)
              complex_piR(1-i,:)  = complex_piR(i,:)

           end do
        end if

!       Symmetries on equator.

        if (eqsym.and.ownequator) then
           do j=1,ghost

              alpha(:,1-j)   = alpha(:,j)
              dtalpha(:,1-j) = dtalpha(:,j)

              phi(:,1-j)   = phi(:,j)
              dtphi(:,1-j) = dtphi(:,j)

              H(:,1-j)   = H(:,j)
              KTH(:,1-j) = KTH(:,j)

              beta_p(:,1-j)   = beta_p(:,j)
              dtbeta_p(:,1-j) = dtbeta_p(:,j)

              complex_phiR(:,1-j) = complex_phiR(:,j)
              complex_piR(:,1-j)  = complex_piR(:,j)

           end do
        end if


!       ***********************
!       ***   SYNCHRONIZE   ***
!       ***********************

!       If we have more than one processor we must now
!       synchronize ghost zones.

        if (size>1) then

           call sync(alpha)
           call sync(dtalpha)

           call sync(phi)
           call sync(dtphi)

           call sync(H)
           call sync(KTH)

           call sync(beta_p)
           call sync(dtbeta_p)

           call sync(complex_phiR)
           call sync(complex_piR)

        end if


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
     t(box,level) = t(box,level) + 1.d0/2**level


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
     call bosonstep(level+1,method)
     call bosonstep(level+1,method)
  end if


! *************************************
! ***   DO WE NEED TO SYNC BOXES?   ***
! *************************************

! Check if refinement boxes at this level intersect,
! and if they do make sure they agree. We basically
! just copy data from the interior of one box to
! the other.  This is similar to synchronization
! across inter-processor boundaries, but in this
! case it is across different refinement boxes on
! the same time level.


! ****************************************************
! ***   RESTRICT FINE GRID DATA INTO COARSE GRID   ***
! ****************************************************

! Restrict the data from the fine to the coarse grid
! across all boxes when both levels coincide in time.
! This restriction does not change data in the current
! grid level, but rather in the coarser level.


! ***************
! ***   END   ***
! ***************

  end subroutine rotbosonstep
