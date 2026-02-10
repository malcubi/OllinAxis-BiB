
  subroutine initial

! ************************
! ***   INITIAL DATA   ***
! ************************

! This subroutine finds initial data.
!
! For simple cases the data can be calculated directly here,
! and for more complicated cases other initial data routines
! are called from here.
!
! NOTE:  Do remember to fill in intial data at all refinement
! boxes and grid levels!

! Include modules.

  use procinfo
  use param
  use arrays
  use derivatives
  use derivadvect

! Extra variables.

  implicit none

  logical contains

  integer box,level        ! Box number and level counters.

  real(8) zero,half,third
  real(8) one,two


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  half  = 0.5d0
  third = 1.d0/3.d0

  one = 1.d0
  two = 2.d0


! ***********************************
! ***   INITIALIZE TO MINKOWSKI   ***
! ***********************************

! By default, initial data is always set to Minkowski first.
! This is in order to avoid the code crashing because the
! metric componets are all zero.

! Loop over boxes and grids levels.

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Lapse.

        alpha = one

!       Shift.

        if (shift/="none") then

           beta_r = zero
           beta_z = zero

           dtbeta_r = zero
           dtbeta_z = zero

           fdriver = one

        end if

!       Conformal factor.

        psi = one
        phi = zero
        chi = one

!       Metric Functions.

        A = one
        B = one
        H = one
        C = zero

        if (angmom) then
           C1 = zero
           C2 = zero
        end if

!       Extrinsic curvature.

        trK = zero

        KTA = zero
        KTB = zero
        KTC = zero
        KTH = zero

        if (angmom) then
           KTC1 = zero
           KTC2 = zero
        end if

!       Deltas.

        Delta_r = zero
        Delta_z = zero

        if (angmom) then
           Delta_p = zero
        end if

!       Regularization variables.

        lambda  = zero
        Alambda = zero

!       Z4c theta variable (for BSSN it does not evolve).

        z4theta = zero

     end do
  end do


! ****************************************
! ***   RESTART FROM CHECKPOINT FILE   ***
! ****************************************

! If we are restarting from a checkpoint file, call
! the restart routine and then return.

  if (idata=="checkpoint") then
     call checkpointrestart
     return
  end if


! *************************
! ***   SCHWARZSCHILD   ***
! *************************

! Initial data for a Schwarzschild black hole
! in isotropic (conformally flat) coordinates.
! In this case the conformal factor is:
!
! psi = 1 + M/(2r)
!
! Notice that for this metric the horizon is
! located (initially) at r=M/2.

! The isotropic lapse (for which the solution
! is static) is given by:
!
! alpha  =  (1 - M/(2r)) / (1 + M/(2r))  =  (2 - psi) / psi
!
! But notice that this lapse will cause numerical
! problems at the origin.
!
! Also, for this metric the horizon is located
! initially at r=M/2.

  if (idata=="schwarzschild") then

!    Loop over boxes and grids levels.

     do box=0,Nb
        do level=min(1,box),Nl(box)

!          Point to current grid.

           call currentgrid(box,level,grid(box,level))

!          Conformal factor.

           psi = one + half*BH1mass/abs(rr)
           phi = log(psi)
           chi = one/psi**dble(chipower)

           psi2 = psi**2
           psi4 = psi**4

!          Lapse.

           if (ilapse=="isotropic") then
              alpha = (one - half*BH1mass/abs(rr))/(one + half*BH1mass/abs(rr))
           end if

        end do
     end do

  end if


! ****************
! ***   KERR   ***
! ****************

! Initial data for a Kerr black hole in quasi-isotropic coordinates.
! See: S.R. Brandt and E. Seidel, Phys.Rev.D 54, 1403.

  if (idata=="kerr") then

!    Sanity check.

     if (.not.angmom) then
        if (rank==0) then
           print *
           print *, 'Kerr data has angular momentum.'
           print *, 'Try setting:  angmom = .true.'
           print *, 'Aborting (subroutine initial.f90)'
           print *
        end if
        call die
     end if

     if (shift=="none") then
        if (rank==0) then
           print *
           print *, 'Kerr data has requires a shift.'
           print *, 'Aborting (subroutine initial.f90)'
           print *
        end if
        call die
     end if

!    Loop over boxes and grids levels.

     do box=0,Nb
        do level=min(1,box),Nl(box)

!          Point to current grid.

           call currentgrid(box,level,grid(box,level))

!          Construct Kerr data.

           call idata_kerr

        end do
     end do

  end if


! ***************************
! ***   BRILL-LINDQUIST   ***
! ***************************

! Initial data for Brill-Lindquist data for two black
! holes in isotropic (conformally flat) coordinates.
! In this case the conformal factor is:
!
! psi = 1 + SUM(M_i/2 abs(rr-rr_i))
!
! where M_i are the bare mass parameters.

  if (idata=="BrillLindquist") then

!    Loop over boxes and grids levels.

     do box=0,Nb
        do level=min(1,box),Nl(box)

!          Point to current grid.

           call currentgrid(box,level,grid(box,level))

!          Conformal factor.

           if (N_BH==1) then

              psi = one + half*BH1mass/sqrt(r**2 + (z-BH1Z)**2)

!             Add reflected black holes if required.

              if (eqsym.and.duplicate_bhs) then
                 psi = psi + half*BH1mass/sqrt(r**2 + (z+BH1Z)**2)
              end if

           else if (N_BH==2) then

              psi = one + half*(BH1mass/sqrt(r**2 + (z-BH1Z)**2) + BH2mass/sqrt(r**2 + (z-BH2Z)**2))

!             Add reflected black holes if required.

              if (eqsym.and.duplicate_bhs) then
                 psi = psi + half*(BH1mass/sqrt(r**2 + (z+BH1Z)**2) + BH2mass/sqrt(r**2 + (z+BH2Z)**2))
              end if

           end if

!          Phi,chi,psi2,spi4.

           phi = log(psi)
           chi = one/psi**dble(chipower)

           psi2 = psi**2
           psi4 = psi**4

        end do
     end do

  end if


! ***********************
! ***   BRILL WAVES   ***
! ***********************

! Initial data for nonlinear Brill waves.

  if (idata=="BrillWave") then
     call idata_brillwave
  end if


! ************************
! ***   TEUKOLSKY GW   ***
! ************************

! Teukolsky gravitational wave.

  if (idata=="testgw") then
     !call testgw

     do box=0,Nb
        do level=min(1,box),Nl(box)

!          Point to current grid.

           call currentgrid(box,level,grid(box,level))

!          Conformal metric coefficients.

           psi = 1.d0

           A = 1.d0 - 24.d0*exp(-rr**2)*teukolsky_a0*(-1.d0 + z**2*r**2)
           B = 1.d0 - 24.d0*exp(-rr**2)*teukolsky_a0*(2.d0 - 4.d0*r**2 + r**4)
           H = (1.d0 + 24.d0*exp(-rr**2)*teukolsky_a0*(1.d0 + (-4.d0 + z**2)*r**2 + r**4))
           C = 24.d0*exp(-rr**2)*z*teukolsky_a0*(-2.d0 + r**2)

        end do
     end do
     !return
  end if


! ************************
! ***   SCALAR PULSE   ***
! ************************

! This initial data corresponds to a small initial pulse
! in the real scalar field.

  if (idata=="scalarpulse") then
     if (contains(mattertype,"scalar")) then
        call idata_scalarpulse
     else
        if (rank==0) then
           print *, 'Scalar pulse initial data needs scalar field type matter.'
           print *, 'Aborting (subroutine initial.f90)'
           print *
        end if
        call die
     end if
  end if


! *************************
! ***   COMPLEX PULSE   ***
! *************************

! This initial data corresponds to a small initial pulse
! in the complex scalar field.

  if (idata=="complexpulse") then
     if (contains(mattertype,"complex")) then
        call idata_complexpulse
     else
        if (rank==0) then
           print *, 'Complex pulse initial data needs complex scalar field type matter.'
           print *, 'Aborting (subroutine initial.f90)'
           print *
        end if
        call die
     end if
  end if


! *************************************************
! ***   NO MORE INITIAL DATA AFTER THIS POINT   ***
! *************************************************

! IMPORTANT:  NO MORE INITIAL DATA ROUTINES MUST
! BE ADDED AFTER THIS POINT!
!
! WHAT COMES AFTER APPLIES TO ALL DIFFERENT TYPES
! OF INITIAL DATA.


! *************************************
! ***   FACTOR METRIC DETERMINANT   ***
! *************************************

! For BSSN one should have the determinant
! of the conformal metric equal to its flat
! value (r**2).  But some initial data don't
! do this (Brill waves for example), so here
! we fix it.
!
! NOTE: In fact, the version of BSSN we use
! here should work in any case as long as the
! determinant is time-independent (for lagrangian
! evolutions), but in practice having a non-trivial
! determinant that remains fixed in time can lead
! to late-time instabilities.

! Loop over boxes and grids levels.

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Find metric determinant (divided by r**2).

        if (angmom) then
           hdet = A*B*H - r**2*(H*C**2 + A*C2**2 + r**2*C1*(B*C1 - two*C*C2))
        else
           hdet = (A*B - (r*C)**2)*H
        end if

!       Find new conformal factor.

        phi = phi + log(hdet)/12.d0

!       Find new metric.

        auxarray = hdet**third

        A = A/auxarray
        B = B/auxarray
        C = C/auxarray
        H = H/auxarray

        if (angmom) then
           C1 = C1/auxarray
           C2 = C2/auxarray
        end if

!       psi, psi2, psi4.

        psi  = exp(phi)
        psi2 = psi**2
        psi4 = psi**4

     end do
  end do


! *************************************************
! ***   SUBTRACT TRACE OF EXTRINSIC CURVATURE   ***
! *************************************************

! For all initial data we must have KT traceless.
! Here we just make sure this is true.

! Loop over boxes and grids levels.

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Calculate the trace of KT which should be zero.

        if (.not.angmom) then
           auxarray = g_A*KTA + g_B*KTB + g_H*KTH + two*r**2*g_C*KTC
        else
           auxarray = g_A*KTA + g_B*KTB + g_H*KTH &
                    + two*r**2*(g_C*KTC + r**2*g_C1*KTC1 + g_C2*KTC2)
        end if

!       Subtract trace.

        KTA = KTA - third*A*auxarray
        KTB = KTB - third*B*auxarray
        KTH = KTH - third*H*auxarray
        KTC = KTC - third*C*auxarray

        if (angmom) then
           KTC1 = KTC1 - third*C1*auxarray
           KTC2 = KTC2 - third*C2*auxarray
        end if

     end do
  end do


! ******************************
! ***   PRECOLLAPSED LAPSE   ***
! ******************************

! Initialize lapse to an inverse power of psi if needed.

  if (ilapse=="psiminus2") then

     do box=0,Nb
        do level=min(1,box),Nl(box)
           call currentgrid(box,level,grid(box,level))
           alpha = one/psi**2
        end do
     end do

  else if (ilapse=="psiminus4") then

     do box=0,Nb
        do level=min(1,box),Nl(box)
           call currentgrid(box,level,grid(box,level))
           alpha = one/psi**4
        end do
     end do

  end if


! ******************************
! ***   LAPSE PERTURBATION   ***
! ******************************

! Add perturbation to initial lapse.

  if (lapsepert=="gauss") then

!    Loop over boxes and levels.

     do box=0,Nb
        do level=min(1,box),Nl(box)

!          Point to current grid.

           call currentgrid(box,level,grid(box,level))

!          Shell type perturbation.  This is a spheroidal
!          shell around the origin.  The idea here is to
!          have a gaussian of the form:
!
!                                 2               2       2
!          exp[ - (R  -  (rho0 cos(theta) + z0 sin(theta)) )
!
!                           2   2                2   2
!               / (sigma_rho cos(theta) + sigma_z sin(theta))
!
!
!          and remember that:  sin(theta)=z/R, cos(theta)=rho/R.
!
!          We add a reflected gaussian on R to guarantee symmetry
!          on the origin.

           if (ilapse=='shell') then

              if ((lapse_r0==0.d0).and.(lapse_z0==0.d0)) then
                 alpha = alpha + lapse_a0*exp(-(rr**4/(r**2*lapse_sr0**2 + z**2*lapse_sz0**2)))
              else
                 alpha = alpha + lapse_a0 &
                       *(exp(-(rr-(r**2*lapse_r0+z**2*lapse_z0)/rr**2)**2 &
                       *rr**2/(r**2*lapse_sr0**2 + z**2*lapse_sz0**2)) &
                       + exp(-(rr+(r**2*lapse_r0+z**2*lapse_z0)/rr**2)**2 &
                       *rr**2/(r**2*lapse_sr0**2 + z**2*lapse_sz0**2)))
              end if

!          Torus type perturbation.  This is basically
!          just the product of two gaussians in r and z
!          centered at r0 and z0 respectively.
!
!          We add reflected gaussians on r (and z for
!          equatorial symmetry) to guarantee the symmetry
!          on the axis (and equator).

           else if (ilapse=='torus') then

              if ((lapse_r0==0.d0).and.(lapse_z0==0.d0)) then
                 alpha = alpha + lapse_a0*exp(-(r/lapse_sr0)**2)*exp(-(z/lapse_sz0)**2)
              else
                 if (eqsym) then
                    alpha = alpha + lapse_a0 &
                          *(exp(-((r-lapse_r0)/lapse_sr0)**2) + exp(-((r+lapse_r0)/lapse_sr0)**2)) &
                          *(exp(-((z-lapse_z0)/lapse_sz0)**2) + exp(-((z+lapse_z0)/lapse_sz0)**2))
                 else
                    alpha = alpha + lapse_a0*exp(-((z-lapse_z0)/lapse_sz0)**2) &
                          *(exp(-((r-lapse_r0)/lapse_sr0)**2) + exp(-((r+lapse_r0)/lapse_sr0)**2))
                 end if
              end if

           end if

        end do
     end do

  end if


! ******************************
! ***   SHIFT PERTURBATION   ***
! ******************************

! Add perturbation to initial shift.

  if (shiftpert=="gauss") then

!    Sanity check.

     if ((shift=="none").or.(shift=="zero")) then
        print *
        print *, 'shift=(none,zero) is incompatible with adding a perturbation.'
        print *
        call die
     end if

!    Loop over boxes and levels.

     do box=0,Nb
        do level=min(1,box),Nl(box)

!          Point to current grid.

           call currentgrid(box,level,grid(box,level))

!          Shell type perturbation.

           if (ishift=='shell') then

              if ((shift_r0==0.d0).and.(shift_z0==0.d0)) then
                 beta_r = beta_r + (r/rr)*shift_a0*exp(-(rr**4/(r**2*shift_sr0**2 + z**2*shift_sz0**2)))
                 beta_z = beta_z + (z/rr)*shift_a0*exp(-(rr**4/(r**2*shift_sr0**2 + z**2*shift_sz0**2)))
              else
                 beta_r = beta_r + (r/rr)*shift_a0 &
                       *(exp(-(rr-(r**2*shift_r0+z**2*shift_z0)/rr**2)**2 &
                       *rr**2/(r**2*shift_sr0**2 + z**2*shift_sz0**2)) &
                       - exp(-(rr+(r**2*shift_r0+z**2*shift_z0)/rr**2)**2 &
                       *rr**2/(r**2*shift_sr0**2 + z**2*shift_sz0**2)))
                 beta_z = beta_z + (z/rr)*shift_a0 &
                       *(exp(-(rr-(r**2*shift_r0+z**2*shift_z0)/rr**2)**2 &
                       *rr**2/(r**2*shift_sr0**2 + z**2*shift_sz0**2)) &
                       - exp(-(rr+(r**2*shift_r0+z**2*shift_z0)/rr**2)**2 &
                       *rr**2/(r**2*shift_sr0**2 + z**2*shift_sz0**2)))
              end if

!          Torus type perturbation.

           else if (ishift=='torus') then

           end if

        end do
     end do

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine initial

