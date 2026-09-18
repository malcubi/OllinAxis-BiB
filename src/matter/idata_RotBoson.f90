
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
!   2          2     2      2       2
! dl  =  A ( dr  + dz  ) + r  H dphi
!
! so we have two non-trivial metric cofficientes: (A,H).  We use this form
! of the metric to be consistent with the paper by Ontañon and Alcubierre:
! Class. Quatum. Grav. 38 (2021) 154003.
!
! In addition, we also need the lapse function "alpha", and since the star is
! rotating we also need a non-zero angular shift component "beta_p" (index up).
!
! We assume that the space-time is static so that all metric functions are
! time independent.  Since there is a non-trivial shift vector this means
! that not all the extrinsic curvature components vanish.  We find:
!
! KTC1 = H/(2*alpha) Dr_beta_p / r
! KTC2 = H/(2*alpha) Dz_beta_p
!
! All other components of K (and its trace) vanish.
!
! For the complex scalar field we use the ansatz (the sign convention here
! is opposite to the one in the routine idat_BosonstarCF.f90):
!
! Phi(t,r,z,phi) = phi(r,z) exp[ - i (omega t - L phi) ]
!
! where L is an integer that corresponds to the angular momentum "quantum number".
! This ansatz guarantees that the stress-energy tensor is time-independent and
! and also independent of the angle, so that the solution will be axi-symmetric.
!
! DO NOT CONFUSE the scalar field Phi(r,z) with the angular coordinate "phi".
!
! For the case when L is not zero, the amplitude of the scalar field Phi(r,z)
! has to behave as r**L close to the origin, so in practice we solve for a
! function F(r,z) such that:
!
!             L
! Phi(r,z) = r  F(r,z)
!
!
! the parameter "boson_phi0" the will then correspond to the value at r=z=0 of
! F(r,z) (and not of Phi(r,z)).
!
! The equations to be solved are found in the following way (I will not write
! the full equations here):
!
! 1) The equation for the scalar field Phi (in fact F) is the Klein-Gordon equation.
!
! 2) The equation for alpha it obtained from the maximal slicing condition.
!
! 3) The equation for beta_p is ontained from the angular component of the
!    momentum constraint.
!
! 4) The equation for the metric function H comes from the angular component
!    of the ADM evolution equations, by asking for dK_[phi,phi]/dt=0.
!
! 5) The equation for the metric function A is obtained from a combination of
!    the Hamiltonian constraint and the equation for H above.
!
! Finally, we meed to introduce a regularization variable "lambda" defined as:
!
! lambda := (A-H)/r**2
!
! This is needed to make sure the system is regular at the axis of symmetry r=0.
!
! 6) The equation for lambda is then obtained from a combination of the ADM
!    equations for dK_[r,r]/dt=0 and dK_[phi,phi]/dt=0 that is chosen in order
!    to have purely regular terms.
!
! We solve the system by turning the elliptic equations into hyperbolic wave
! equations and "evolving" to a steady state.
!
! This method is very slow (specially at high resolution), but seems robust. 
! It can be improved by using a good initial guess (I'll try this later),
! and maybe something like multi-grid since it converges much faster.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo
  use derivatives

! Extra variables.

  implicit none

  logical flag1,flag2            ! Interpolation flags.

  integer box,level              ! Box number and level counters.
  integer i,j,m,n                ! Counters.
  integer step                   ! Iteration counter.
  integer Nlmax_old              ! Original number of levels.
  integer NNloc,NNtot            ! Total number of grid points.
  integer miniter                ! Minimum number of iterations.
  integer :: maxiter = 100000    ! Maximum number of iterations.

  real(8) lres,gres              ! Local and global residuals.
  real(8) r0,z0,interp           ! For interpolation.
  real(8) cfac                   ! Courant parameter.
  real(8) one,half,smallpi       ! Numbers.
  real(8) integral               ! Integrating function.
  real(8) aux1,aux2

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

! This initial data is such that hdet is not one,
! so you should not factor it.

  if (factorhdet) then
     print *, 'Rotating boson star initial data is not compatible with factorhdet=.true.'
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

! The metric functions (alpha,A,H) are initialized to 1.
!
! The angular shift vector is initialized to 0.
! The time derivatives are set initially to 0.

  do box=0,Nb
     do level=min(1,box),Nl(box)

        call currentgrid(box,level,grid(box,level))

        alpha = one
        A     = one
        H     = one

        lambda = 0.d0
        beta_p = 0.d0

!       The scalar field is initialized to a gaussian centered
!       on the origin with the correct amplitude.  Remember that
!       phi(r) = r^L F(t).

        complex_phiR = boson_phi0*exp(-rr**2/(complex_l+1.d0)**2)

!       Set all time derivatives to 0.

        dtalpha  = 0.d0      ! Time derivative of alpha.
        KTA      = 0.d0      ! Time derivatve of A.
        KTH      = 0.d0      ! Time derivatve of H.
        Alambda  = 0.d0      ! Time derivative of lambda.
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

           grabvar => A
           call save1Dvariable('rotboson_A',directory,box,level,outparallel,'new')
           !call save2Dvariable('rotboson_A',directory,box,level,outparallel,'new')

           grabvar => H
           call save1Dvariable('rotboson_H',directory,box,level,outparallel,'new')
           !call save2Dvariable('rotboson_H',directory,box,level,outparallel,'new')

           grabvar => lambda
           call save1Dvariable('rotboson_lambda',directory,box,level,outparallel,'new')
           !call save2Dvariable('rotboson_lambda',directory,box,level,outparallel,'new')

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
           lres = lres + abs(grid(0,0)%complex_phiR(i,j)-grid(0,0)%complex_phiR_p(i,j)) &
                + abs(grid(0,0)%alpha(i,j) - grid(0,0)%alpha_p(i,j)) &
                + abs(grid(0,0)%beta_p(i,j) - grid(0,0)%beta_p_p(i,j)) &
                + abs(grid(0,0)%A(i,j) - grid(0,0)%A_p(i,j)) &
                + abs(grid(0,0)%H(i,j) - grid(0,0)%H_p(i,j)) &
                + abs(grid(0,0)%lambda(i,j) - grid(0,0)%lambda_p(i,j))
        End do
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

              grabvar => A
              call save1Dvariable('rotboson_A',directory,box,level,outparallel,'old')
              !call save2Dvariable('rotboson_A',directory,box,level,outparallel,'old')

              grabvar => H
              call save1Dvariable('rotboson_H',directory,box,level,outparallel,'old')
              !call save2Dvariable('rotboson_H',directory,box,level,outparallel,'old')

              grabvar => lambda
              call save1Dvariable('rotboson_lambda',directory,box,level,outparallel,'old')
              !call save2Dvariable('rotboson_lambda',directory,box,level,outparallel,'old')

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

! If we converged and we have refinement levels, inject the
! coarse solution into fine grids and solve again.

  !if ((step/=WE_maxiter).and.(Nlmax/=Nlmax_old)) then
  if (Nlmax/=Nlmax_old) then

!    Set again Nlmax to its original value.

     Nlmax = Nlmax_old

!    Set time derivatives in the base grid back to 0.

     grid(0,0)%dtalpha = 0.d0
     grid(0,0)%KTA = 0.d0
     grid(0,0)%KTH = 0.d0
     grid(0,0)%dtbeta_p = 0.d0
     grid(0,0)%Alambda = 0.d0
     grid(0,0)%complex_piR = 0.d0

!    Iterate over boxes and levels higher than 0.

     do box=0,Nb
        do level=1,Nl(box)

!          Point to current grid.

           call currentgrid(box,level,grid(box,level))

!          Set time derivatives to 0.

           dtalpha  = 0.d0
           KTA = 0.d0
           KTH = 0.d0
           Alambda = 0.d0
           dtbeta_p = 0.d0
           complex_piR = 0.d0

!          Loop over ALL points in the current grid (and I do mean
!          all of them, not just the ones in the current processor).

           do m=0,Nrbox(box)+ghost-1
              do n=0,Nzbox(box)+ghost-1

!                Figure out (r,z) position for interpolation.

                 r0 = rminl(box,level) + dble(m)*drl(level)
                 z0 = zminl(box,level) + dble(n)*dzl(level)

!                Figure out to which grid point this (r0,z0) values
!                would correspond in the local processor at the
!                fine grid level.

                 i = nint((r0-r(1-ghost,0))/drl(level)) + 1 - ghost
                 j = nint((z0-z(0,1-ghost))/dzl(level)) + 1 - ghost

!                Notice now that the (i,j) values above might be outside
!                the range of the current processor.  This means that
!                we should not try to access this location as it belongs
!                to another processor.

                 if ((i>=1-ghost).and.(i<=Nrl(box,rank)).and.(j>=1-ghost).and.(j<=Nzl(box,rank))) then
                    flag1 = .true.
                 else
                    flag1 = .false.
                 end if

!                Interpolate variables from coarse grid level. But we only update
!                the values if the location belongs to us.

                 interpvar => grid(0,0)%alpha
                 aux1 = interp(0,0,r0,z0,flag2)
                 call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

                 if (flag1) then
                    alpha(i,j) = aux2
                 end if

                 interpvar => grid(0,0)%A
                 aux1 = interp(0,0,r0,z0,flag2)
                 call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

                 if (flag1) then
                    A(i,j) = aux2
                 end if

                 interpvar => grid(0,0)%H
                 aux1 = interp(0,0,r0,z0,flag2)
                 call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

                 if (flag1) then
                    H(i,j) = aux2
                 end if

                 interpvar => grid(0,0)%lambda
                 aux1 = interp(0,0,r0,z0,flag2)
                 call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

                 if (flag1) then
                    lambda(i,j) = aux2
                 end if

                 interpvar => grid(0,0)%beta_p
                 aux1 = interp(0,0,r0,z0,flag2)
                 call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

                 if (flag1) then
                    beta_p(i,j) = aux2
                 end if

                 interpvar => grid(0,0)%complex_phiR
                 aux1 = interp(0,0,r0,z0,flag2)
                 call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

                 if (flag1) then
                    complex_phiR(i,j) = aux2
                 end if

              end do
           end do

        end do
     end do

!    Set time and time step counters back to zero,
!    and restart iterations.

     s = 0
     t = 0.d0

     t1 = 0.d0
     t2 = 0.d0

     goto 100

  end if


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

        complex_phiR = r**complex_l*complex_phiR

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

        complex_piI = - (boson_omega + complex_l*beta_p)*complex_phiR/alpha

!       Set the imaginary part of the scalar field and its spatial derivatives to zero.

        complex_phiI  = 0.d0
        complex_xiI_r = 0.d0
        complex_xiI_z = 0.d0

!       Set time derivative of the real part of phi to zero.

        complex_piR  = 0.d0

     end do
  end do


! ******************************************************
! ***   SET B=A, PUT TIME DERIVATIVES BACK TO ZERO   ***
! ***        AND FIND EXTRINSIC CURVATURE            ***
! ******************************************************

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Remember that for this initial data we have B = A.

        B = A

!       Set again (dtalpha,KTA,KTH,dtbeta_p,Alambda) to 0.

        dtalpha = 0.d0
        KTA = 0.d0
        KTH = 0.d0
        dtbeta_p = 0.d0
        Alambda = 0.d0

        sdtalpha = 0.d0
        sKTA = 0.d0
        sKTH = 0.d0
        sdtbeta_p = 0.d0
        sAlambda = 0.d0

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
        Dr_beta_p = diff1r(+1)
        Dz_beta_p = diff1z(+1)

        KTC1 = 0.5d0*(H/alpha)*Dr_beta_p/r
        KTC2 = 0.5d0*(H/alpha)*Dz_beta_p

     end do
  end do


! *****************************
! ***   TOLMAN-KOMAR MASS   ***
! *****************************

! The Tolman-Komar mass only makes sense for static solutions.
! It is based on the existence on a Killing field, but can be
! expressed as a volume integral that depends on the lapse
! function and the stress-energy of matter.
!
! The general expression is:
!
!             /                               m
! mass_TK  =  | [ alpha (rho + trS) - 2 J beta  ] dV
!             /                          m
!
! with alpha the lapse, rho the energy density, trS the trace
! of the stress tensor, beta^m the contravariant shift vector,
! J_m the covariant momentum density, and dV the physical
! volume element:
!
! dV = 2 pi r psi**6 sqrt(hdet) dr dz
!
! The factor 2*pi comes from the integral over the angle.
!
! Notice that for a complex scalar field we have:
!
! rho + trS  =  2 ( Pi**2 - V )
!
! so that the integral becomes:
!
!                  /                               p        6
! mass_TK  =  4 pi | [ alpha (PiI**2 - V) - J  beta  ] r psi sqrt(hdet) dr dz
!                  /                         p

! At the moment I only do the integral in the coarse grid.

  call currentgrid(0,0,grid(0,0))

! Scalar field potential.

  call potential

! Conformal metric determinant (divided by r**2). In this
! case we just have hdet = A**2*H.  Notice also that psi=1.

  hdet = A**2*H

! The covariant angular momentum densityis given by: - L PiI phiR.
! Notice: I use the array J_p which normally is the contravariant
! component, but it will be calculated correctly later.

  J_p = - complex_l*complex_piI*complex_phiR

! Integrate.

  auxarray = (4.d0*smallpi*r*sqrt(hdet))*(alpha*(complex_piI**2 - complex_V) - J_p*beta_p)

  mass_TK = integral(0,0,auxarray)

  if (rank==0) then
     write (*,'(A,ES13.6)') ' Tolman-Komar mass (mass_TK) = ',mass_TK
     print *
  end if


! *****************************************
! ***   TOLMAN-KOMAR ANGULAR MOMENTUM   ***
! *****************************************


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
  integer i,j,k        ! Counters.
  integer iter         ! Counter for internal iterations.
  integer niter        ! Number of internal iterations.
  integer bmax         ! Number of boxes at this level.
  integer bbox

  real(8) dtw          ! Internal time step.
  real(8) weight       ! Weight for rk4.
  real(8) smallpi      ! Numbers.
  real(8) aux

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

!    Old values of (A,H) and time derivatives.

     A_p = A
     KTA_p = KTA

     H_p   = H
     KTH_p = KTH

!    Old values of lambda and Alambda.

     lambda_p  = lambda
     Alambda_p = Alambda

!    Old values of shift and time derivative.

     beta_p_p = beta_p
     dtbeta_p_p = dtbeta_p

!    Old values of complex_phiR and time derivative.

     complex_phiR_p = complex_phiR
     complex_piR_p  = complex_piR

!    Save old values of internal boundaries.

     if (level>0) then

        do i=0,ghost-1

!          alpha.

           alpha_bound_rL(i,:,2) = alpha_bound_rL(i,:,1)
           alpha_bound_rL(i,:,1) = alpha_bound_rL(i,:,0)
           alpha_bound_rL(i,:,0) = alpha(1-ghost+i,:)
           alpha_bound_rR(i,:,2) = alpha_bound_rR(i,:,1)
           alpha_bound_rR(i,:,1) = alpha_bound_rR(i,:,0)
           alpha_bound_rR(i,:,0) = alpha(Nr-i,:)

           alpha_bound_zL(:,i,2) = alpha_bound_zL(:,i,1)
           alpha_bound_zL(:,i,1) = alpha_bound_zL(:,i,0)
           alpha_bound_zL(:,i,0) = alpha(:,1-ghost+i)
           alpha_bound_zR(:,i,2) = alpha_bound_zR(:,i,1)
           alpha_bound_zR(:,i,1) = alpha_bound_zR(:,i,0)
           alpha_bound_zR(:,i,0) = alpha(:,Nz-i)

!          dtalpha.

           dtalpha_bound_rL(i,:,2) = dtalpha_bound_rL(i,:,1)
           dtalpha_bound_rL(i,:,1) = dtalpha_bound_rL(i,:,0)
           dtalpha_bound_rL(i,:,0) = dtalpha(1-ghost+i,:)
           dtalpha_bound_rR(i,:,2) = dtalpha_bound_rR(i,:,1)
           dtalpha_bound_rR(i,:,1) = dtalpha_bound_rR(i,:,0)
           dtalpha_bound_rR(i,:,0) = dtalpha(Nr-i,:)

           dtalpha_bound_zL(:,i,2) = dtalpha_bound_zL(:,i,1)
           dtalpha_bound_zL(:,i,1) = dtalpha_bound_zL(:,i,0)
           dtalpha_bound_zL(:,i,0) = dtalpha(:,1-ghost+i)
           dtalpha_bound_zR(:,i,2) = dtalpha_bound_zR(:,i,1)
           dtalpha_bound_zR(:,i,1) = dtalpha_bound_zR(:,i,0)
           dtalpha_bound_zR(:,i,0) = dtalpha(:,Nz-i)

!          A.

           A_bound_rL(i,:,2) = A_bound_rL(i,:,1)
           A_bound_rL(i,:,1) = A_bound_rL(i,:,0)
           A_bound_rL(i,:,0) = A(1-ghost+i,:)
           A_bound_rR(i,:,2) = A_bound_rR(i,:,1)
           A_bound_rR(i,:,1) = A_bound_rR(i,:,0)
           A_bound_rR(i,:,0) = A(Nr-i,:)

           A_bound_zL(:,i,2) = A_bound_zL(:,i,1)
           A_bound_zL(:,i,1) = A_bound_zL(:,i,0)
           A_bound_zL(:,i,0) = A(:,1-ghost+i)
           A_bound_zR(:,i,2) = A_bound_zR(:,i,1)
           A_bound_zR(:,i,1) = A_bound_zR(:,i,0)
           A_bound_zR(:,i,0) = A(:,Nz-i)

!          KTA.

           KTA_bound_rL(i,:,2) = KTA_bound_rL(i,:,1)
           KTA_bound_rL(i,:,1) = KTA_bound_rL(i,:,0)
           KTA_bound_rL(i,:,0) = KTA(1-ghost+i,:)
           KTA_bound_rR(i,:,2) = KTA_bound_rR(i,:,1)
           KTA_bound_rR(i,:,1) = KTA_bound_rR(i,:,0)
           KTA_bound_rR(i,:,0) = KTA(Nr-i,:)

           KTA_bound_zL(:,i,2) = KTA_bound_zL(:,i,1)
           KTA_bound_zL(:,i,1) = KTA_bound_zL(:,i,0)
           KTA_bound_zL(:,i,0) = KTA(:,1-ghost+i)
           KTA_bound_zR(:,i,2) = KTA_bound_zR(:,i,1)
           KTA_bound_zR(:,i,1) = KTA_bound_zR(:,i,0)
           KTA_bound_zR(:,i,0) = KTA(:,Nz-i)

!          H.

           H_bound_rL(i,:,2) = H_bound_rL(i,:,1)
           H_bound_rL(i,:,1) = H_bound_rL(i,:,0)
           H_bound_rL(i,:,0) = H(1-ghost+i,:)
           H_bound_rR(i,:,2) = H_bound_rR(i,:,1)
           H_bound_rR(i,:,1) = H_bound_rR(i,:,0)
           H_bound_rR(i,:,0) = H(Nr-i,:)

           H_bound_zL(:,i,2) = H_bound_zL(:,i,1)
           H_bound_zL(:,i,1) = H_bound_zL(:,i,0)
           H_bound_zL(:,i,0) = H(:,1-ghost+i)
           H_bound_zR(:,i,2) = H_bound_zR(:,i,1)
           H_bound_zR(:,i,1) = H_bound_zR(:,i,0)
           H_bound_zR(:,i,0) = H(:,Nz-i)

!          KTH.

           KTH_bound_rL(i,:,2) = KTH_bound_rL(i,:,1)
           KTH_bound_rL(i,:,1) = KTH_bound_rL(i,:,0)
           KTH_bound_rL(i,:,0) = KTH(1-ghost+i,:)
           KTH_bound_rR(i,:,2) = KTH_bound_rR(i,:,1)
           KTH_bound_rR(i,:,1) = KTH_bound_rR(i,:,0)
           KTH_bound_rR(i,:,0) = KTH(Nr-i,:)

           KTH_bound_zL(:,i,2) = KTH_bound_zL(:,i,1)
           KTH_bound_zL(:,i,1) = KTH_bound_zL(:,i,0)
           KTH_bound_zL(:,i,0) = KTH(:,1-ghost+i)
           KTH_bound_zR(:,i,2) = KTH_bound_zR(:,i,1)
           KTH_bound_zR(:,i,1) = KTH_bound_zR(:,i,0)
           KTH_bound_zR(:,i,0) = KTH(:,Nz-i)

!          beta_p.

           beta_p_bound_rL(i,:,2) = beta_p_bound_rL(i,:,1)
           beta_p_bound_rL(i,:,1) = beta_p_bound_rL(i,:,0)
           beta_p_bound_rL(i,:,0) = beta_p(1-ghost+i,:)
           beta_p_bound_rR(i,:,2) = beta_p_bound_rR(i,:,1)
           beta_p_bound_rR(i,:,1) = beta_p_bound_rR(i,:,0)
           beta_p_bound_rR(i,:,0) = beta_p(Nr-i,:)

           beta_p_bound_zL(:,i,2) = beta_p_bound_zL(:,i,1)
           beta_p_bound_zL(:,i,1) = beta_p_bound_zL(:,i,0)
           beta_p_bound_zL(:,i,0) = beta_p(:,1-ghost+i)
           beta_p_bound_zR(:,i,2) = beta_p_bound_zR(:,i,1)
           beta_p_bound_zR(:,i,1) = beta_p_bound_zR(:,i,0)
           beta_p_bound_zR(:,i,0) = beta_p(:,Nz-i)

!          dtbeta_p.

           dtbeta_p_bound_rL(i,:,2) = dtbeta_p_bound_rL(i,:,1)
           dtbeta_p_bound_rL(i,:,1) = dtbeta_p_bound_rL(i,:,0)
           dtbeta_p_bound_rL(i,:,0) = dtbeta_p(1-ghost+i,:)
           dtbeta_p_bound_rR(i,:,2) = dtbeta_p_bound_rR(i,:,1)
           dtbeta_p_bound_rR(i,:,1) = dtbeta_p_bound_rR(i,:,0)
           dtbeta_p_bound_rR(i,:,0) = dtbeta_p(Nr-i,:)

           dtbeta_p_bound_zL(:,i,2) = dtbeta_p_bound_zL(:,i,1)
           dtbeta_p_bound_zL(:,i,1) = dtbeta_p_bound_zL(:,i,0)
           dtbeta_p_bound_zL(:,i,0) = dtbeta_p(:,1-ghost+i)
           dtbeta_p_bound_zR(:,i,2) = dtbeta_p_bound_zR(:,i,1)
           dtbeta_p_bound_zR(:,i,1) = dtbeta_p_bound_zR(:,i,0)
           dtbeta_p_bound_zR(:,i,0) = dtbeta_p(:,Nz-i)

!          lambda.

           lambda_bound_rL(i,:,2) = lambda_bound_rL(i,:,1)
           lambda_bound_rL(i,:,1) = lambda_bound_rL(i,:,0)
           lambda_bound_rL(i,:,0) = lambda(1-ghost+i,:)
           lambda_bound_rR(i,:,2) = lambda_bound_rR(i,:,1)
           lambda_bound_rR(i,:,1) = lambda_bound_rR(i,:,0)
           lambda_bound_rR(i,:,0) = lambda(Nr-i,:)

           lambda_bound_zL(:,i,2) = lambda_bound_zL(:,i,1)
           lambda_bound_zL(:,i,1) = lambda_bound_zL(:,i,0)
           lambda_bound_zL(:,i,0) = lambda(:,1-ghost+i)
           lambda_bound_zR(:,i,2) = lambda_bound_zR(:,i,1)
           lambda_bound_zR(:,i,1) = lambda_bound_zR(:,i,0)
           lambda_bound_zR(:,i,0) = lambda(:,Nz-i)

!          Alambda.

           Alambda_bound_rL(i,:,2) = Alambda_bound_rL(i,:,1)
           Alambda_bound_rL(i,:,1) = Alambda_bound_rL(i,:,0)
           Alambda_bound_rL(i,:,0) = Alambda(1-ghost+i,:)
           Alambda_bound_rR(i,:,2) = Alambda_bound_rR(i,:,1)
           Alambda_bound_rR(i,:,1) = Alambda_bound_rR(i,:,0)
           Alambda_bound_rR(i,:,0) = Alambda(Nr-i,:)

           Alambda_bound_zL(:,i,2) = Alambda_bound_zL(:,i,1)
           Alambda_bound_zL(:,i,1) = Alambda_bound_zL(:,i,0)
           Alambda_bound_zL(:,i,0) = Alambda(:,1-ghost+i)
           Alambda_bound_zR(:,i,2) = Alambda_bound_zR(:,i,1)
           Alambda_bound_zR(:,i,1) = Alambda_bound_zR(:,i,0)
           Alambda_bound_zR(:,i,0) = Alambda(:,Nz-i)

!          complex_phiR.

           complex_phiR_bound_rL(i,:,2) = complex_phiR_bound_rL(i,:,1)
           complex_phiR_bound_rL(i,:,1) = complex_phiR_bound_rL(i,:,0)
           complex_phiR_bound_rL(i,:,0) = complex_phiR(1-ghost+i,:)
           complex_phiR_bound_rR(i,:,2) = complex_phiR_bound_rR(i,:,1)
           complex_phiR_bound_rR(i,:,1) = complex_phiR_bound_rR(i,:,0)
           complex_phiR_bound_rR(i,:,0) = complex_phiR(Nr-i,:)

           complex_phiR_bound_zL(:,i,2) = complex_phiR_bound_zL(:,i,1)
           complex_phiR_bound_zL(:,i,1) = complex_phiR_bound_zL(:,i,0)
           complex_phiR_bound_zL(:,i,0) = complex_phiR(:,1-ghost+i)
           complex_phiR_bound_zR(:,i,2) = complex_phiR_bound_zR(:,i,1)
           complex_phiR_bound_zR(:,i,1) = complex_phiR_bound_zR(:,i,0)
           complex_phiR_bound_zR(:,i,0) = complex_phiR(:,Nz-i)

!          complex_piR.

           complex_piR_bound_rL(i,:,2) = complex_piR_bound_rL(i,:,1)
           complex_piR_bound_rL(i,:,1) = complex_piR_bound_rL(i,:,0)
           complex_piR_bound_rL(i,:,0) = complex_piR(1-ghost+i,:)
           complex_piR_bound_rR(i,:,2) = complex_piR_bound_rR(i,:,1)
           complex_piR_bound_rR(i,:,1) = complex_piR_bound_rR(i,:,0)
           complex_piR_bound_rR(i,:,0) = complex_piR(Nr-i,:)

           complex_piR_bound_zL(:,i,2) = complex_piR_bound_zL(:,i,1)
           complex_piR_bound_zL(:,i,1) = complex_piR_bound_zL(:,i,0)
           complex_piR_bound_zL(:,i,0) = complex_piR(:,1-ghost+i)
           complex_piR_bound_zR(:,i,2) = complex_piR_bound_zR(:,i,1)
           complex_piR_bound_zR(:,i,1) = complex_piR_bound_zR(:,i,0)
           complex_piR_bound_zR(:,i,0) = complex_piR(:,Nz-i)

        end do

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

!       Derivatives of (A,KTA).

        diffvar => A

        Dr_A = diff1r(+1)
        Dz_A = diff1z(+1)

        Drr_A = diff2r(+1)
        Dzz_A = diff2z(+1)

        diffvar => KTA

        Dr_KTA = diff1r(+1)
        Dz_KTA = diff1z(+1)

!       Derivatives of (H,KTH).

        diffvar => H

        Dr_H = diff1r(+1)
        Dz_H = diff1z(+1)

        Drr_H = diff2r(+1)
        Dzz_H = diff2z(+1)

        diffvar => KTH

        Dr_KTH = diff1r(+1)
        Dz_KTH = diff1z(+1)

!       Derivatives of (lambda,Alambda).

        diffvar => lambda

        Dr_lambda = diff1r(+1)
        Dz_lambda = diff1z(+1)

        Drr_lambda = diff2r(+1)
        Dzz_lambda = diff2z(+1)

        diffvar => Alambda

        Dr_Alambda = diff1r(+1)
        Dz_Alambda = diff1z(+1)

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

        complex_V   = 0.5d0*complex_mass**2*(complex_phiR*r**complex_l)**2
        complex_VPR = complex_mass**2*(complex_phiR*r**complex_l)

!       Find frequency omega.  We solve for omega from the
!       Klein-Gordon equation at the first grid point.

        if (level==Nlmax) then

           if (ownaxis.and.ownequator) then

              aux = - alpha(1,1)**2/A(1,1)/complex_phiR(1,1) &
                  *(Drr_complex_phiR(1,1) + Dzz_complex_phiR(1,1) + (2.d0*complex_l + 1.d0)*Dr_complex_phiR(1,1)/r(1,1) &
                  + (Dr_complex_phiR(1,1)*Dr_alpha(1,1) + Dz_complex_phiR(1,1)*Dz_alpha(1,1))/alpha(1,1) &
                  + 0.5d0*(Dr_complex_phiR(1,1)*Dr_H(1,1) + Dz_complex_phiR(1,1)*Dz_H(1,1))/H(1,1) &
                  + complex_l*complex_phiR(1,1)/r(1,1)*(Dr_alpha(1,1)/alpha(1,1) + 0.5d0*Dr_H(1,1)/H(1,1)) &
                  - A(1,1)*complex_VPR(1,1)/r(1,1)**complex_l - complex_l**2*complex_phiR(1,1)*lambda(1,1)/H(1,1))

              boson_omega = sqrt(abs(aux)) - complex_l*beta_p(1,1)
 
           end if

        end if

        call MPI_BCAST(boson_omega,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!       The sources for (alpha,A,H,lambda,beta_p,complex_phiR)
!       are just (dtalpha,KTA,KTH,Alambda,dtbeta_p,complex_piR).

        salpha = dtalpha

        sA = KTA

        sH = KTH

        slambda = Alambda

        sbeta_p = dtbeta_p

        scomplex_phiR = complex_piR

!       Sources for dtalpha.

        sdtalpha = Drr_alpha + Dzz_alpha + Dr_alpha/r &
                 + 0.5d0/H*(Dr_alpha*Dr_H + Dz_alpha*Dz_H) &
                 - 0.5d0*r**2*H/alpha*(Dr_beta_p**2 + Dz_beta_p**2) &
                 + 8.d0*smallpi*A*(alpha*complex_V &
                 - (r**complex_l*complex_phiR*(complex_l*beta_p + boson_omega))**2/alpha)

!       Sources for dtbeta_p.

        sdtbeta_p = Drr_beta_p + Dzz_beta_p + 3.d0*Dr_beta_p/r &
                  + 1.5d0/H*(Dr_beta_p*Dr_H + Dz_beta_p*Dz_H) &
                  - 1.d0/alpha*(Dr_beta_p*Dr_alpha + Dz_beta_p*Dz_alpha) &
                  - 16.d0*smallpi*complex_l*A/H*(complex_l*beta_p + boson_omega) &
                  *(r**complex_l*complex_phiR/r)**2

!       Sources for KTA.

        sKTA = Drr_A + Dzz_A - (Dr_A**2 + Dz_A**2)/A &
             - A/(alpha*H)*(Dr_alpha*Dr_H + Dz_alpha*Dz_H) &
             - 0.5d0*r**2*A*H/alpha**2*(Dr_beta_p**2 + Dz_beta_p**2) &
             - 2.d0*A*Dr_alpha/(r*alpha) &
             + 8.d0*smallpi*A*r**(2*complex_l)*(Dr_complex_phiR**2 + Dz_complex_phiR**2 &
             + 2.d0*complex_l*complex_phiR*Dr_complex_phiR/r  &
             + A*(complex_phiR*(complex_l*beta_p + boson_omega)/alpha)**2 &
             - lambda*(complex_l*complex_phiR)**2/H)

!       Sources for KTH.

        sKTH = Drr_H + Dzz_H + 2.d0*Dr_H/r - 0.5d0*(Dr_H**2 + Dz_H**2)/H &
             + (Dr_alpha*Dr_H + Dz_alpha*Dz_H)/alpha + 2.d0*H*Dr_alpha/(r*alpha) &
             + (r*H/alpha)**2*(Dr_beta_p**2 + Dz_beta_p**2) &
             + 16.d0*smallpi*A*(H*complex_V + (complex_l*r**complex_l*complex_phiR/r)**2)

!       Sources for complex_piR.

        scomplex_piR = Drr_complex_phiR + Dzz_complex_phiR + (2.d0*complex_l + 1.d0)*Dr_complex_phiR/r &
                     + (Dr_complex_phiR*Dr_alpha + Dz_complex_phiR*Dz_alpha)/alpha &
                     + 0.5d0*(Dr_complex_phiR*Dr_H + Dz_complex_phiR*Dz_H)/H &
                     + complex_l*complex_phiR/r*(Dr_alpha/alpha + 0.5d0*Dr_H/H) &
                     + A*(complex_phiR*(complex_l*beta_p + boson_omega)**2/alpha**2 - complex_VPR/r**complex_l) &
                     - complex_l**2*complex_phiR*lambda/H

!       Sources for Alambda.

        sAlambda = Drr_lambda + Dzz_lambda + 3.d0*Dr_lambda/r &
                 - 0.5d0*(Dr_lambda*Dr_H - Dz_lambda*Dz_H)/H &
                 - 2.d0 *(Dr_lambda*Dr_H + Dz_lambda*Dz_H)/A &
                 - r**2*(Dr_lambda**2 + Dz_lambda**2)/A &
                 - (Dr_lambda*Dr_alpha - Dz_lambda*Dz_alpha)/alpha &
                 - 4.d0*lambda*(lambda + r*Dr_lambda)/A &
                 - 2.d0*lambda*Dr_alpha/(r*alpha) &
                 + lambda*Dr_H/r*(1.d0/H - 4.d0/A) &
                 - (Dr_H/r)**2*(1.d0/A + 0.5d0/H) &
                 + lambda*(Dz_H**2/(A*H) - 0.5d0*Dr_H**2/H**2) &
                 - 2.d0*Dr_alpha*Dr_H/(r**2*alpha) &
                 + 2.d0*lambda*Drr_alpha/alpha + lambda*Drr_H/H &
                 - (H/alpha)**2*(Dr_beta_p**2 + Dz_beta_p**2) - A*H/alpha**2*Dr_beta_p**2 &
                 + 16.d0*smallpi*A*(lambda*complex_V &
                 + r**(2*(complex_l-1))*(Dr_complex_phiR/r) &
                 *(2.d0*complex_l*complex_phiR + r*Dr_complex_phiR))

        auxarray = Dr_alpha/r
        diffvar => auxarray
        DD_alphar = diff1r(+1)

        sAlambda = sAlambda + 2.d0*H*DD_alphar/(r*alpha)

        auxarray = Dr_H/r
        diffvar => auxarray
        DD_alphar = diff1r(+1)

        sAlambda = sAlambda + DD_alphar/r

!       Damping term.  We need it for complex_piR always or the iterations
!       go unstable.  And it needs to be larger for larger complex_l.
!
!       For the other quantities I only add it for fine grids since by
!       then the coarse grid has essentially converged.

        scomplex_piR = scomplex_piR - (2.d0*complex_l + 1.d0)*complex_piR

        if (Nlmax>0) then
           sdtalpha = sdtalpha - dtalpha
           sKTA = sKTA - KTA
           sKTH = sKTH - KTH
           sAlambda = sAlambda - Alambda
           sdtbeta_p = sdtbeta_p - dtbeta_p
        end if

!       And add some dissipation to reduce high frequency noise.

        evolvevar => dtalpha
        sourcevar => sdtalpha
        call dissipation(+1,+1,WE_diss)

        evolvevar => KTA
        sourcevar => sKTA
        call dissipation(+1,+1,WE_diss)

        evolvevar => KTH
        sourcevar => sKTH
        call dissipation(+1,+1,WE_diss)

        evolvevar => Alambda
        sourcevar => sAlambda
        call dissipation(+1,+1,WE_diss)

        evolvevar => dtbeta_p
        sourcevar => sdtbeta_p
        call dissipation(+1,+1,WE_diss)

        evolvevar => complex_piR
        sourcevar => scomplex_piR
        call dissipation(+1,+1,WE_diss)

!       But set the source for complex_phiR at point (1,1)
!       such that the value at r=0 does not change.

        if ((rank==0).and.(level==Nlmax)) then
           !scomplex_piR(1,1) = 0.d0
           scomplex_piR(1,1) = scomplex_piR(2,2)/9.d0
        end if

!       I also fix the source for KTA at the grid points closer
!       to the origin to guarantee that:  A = H + r**2 lambda.

        if ((rank==0).and.(level==Nlmax)) then
           sKTA(1,0) = sKTH(1,0) + r(1,0)**2*sAlambda(1,0)
           sKTA(0,1) = sKTH(0,1) + r(0,1)**2*sAlambda(0,1)
           sKTA(1,1) = sKTH(1,1) + r(1,1)**2*sAlambda(1,1)
        end if

!       Symmetries on axis.

        if (ownaxis) then
           do i=1,ghost
              sdtalpha(1-i,:)  = sdtalpha(i,:)
              sKTA(1-i,:)      = sKTA(i,:)
              sKTH(1-i,:)      = sKTH(i,:)
              sAlambda(1-i,:)  = sAlambda(i,:)
              sdtbeta_p(1-i,:) = sdtbeta_p(i,:)
              scomplex_piR(1-i,:) = scomplex_piR(i,:)
           end do
        end if

!       Symmetries on equator.

        if (eqsym.and.ownequator) then
           do j=1,ghost
              sdtalpha(:,1-j)  = sdtalpha(:,j)
              sKTA(:,1-j)      = sKTA(:,j)
              sKTH(:,1-j)      = sKTH(:,j)
              sAlambda(:,1-j)  = sAlambda(:,j)
              sdtbeta_p(:,1-j) = sdtbeta_p(:,j)
              scomplex_piR(:,1-j) = scomplex_piR(:,j)
           end do
        end if


!       *******************************
!       ***   BOUNDARY CONDITIONS   ***
!       *******************************

!       Radiative boundaries for all equations. However, we use the fact that
!       asymtotically (alpha,A,H) decay as 1/r, while beta_p decays as 1/r**3
!       (see paper by Ontañon and Alcubierre).

        if (level==0) then

!          Radiative conditions at r boundary.

           if (mod(rank+1,nprocr)==0) then
              i = Nr
              sdtalpha(i,:)     = - (r(i,:)*Dr_dtalpha(i,:)     + z(i,:)*Dz_dtalpha(i,:)     + dtalpha(i,:))/rr(i,:)
              sKTA(i,:)         = - (r(i,:)*Dr_KTA(i,:)         + z(i,:)*Dz_KTA(i,:)         + KTA(i,:))/rr(i,:)
              sKTH(i,:)         = - (r(i,:)*Dr_KTH(i,:)         + z(i,:)*Dz_KTH(i,:)         + KTH(i,:))/rr(i,:)
              sAlambda(i,:)     = - (r(i,:)*Dr_Alambda(i,:)     + z(i,:)*Dz_Alambda(i,:)     + Alambda(i,:))/rr(i,:)
              sdtbeta_p(i,:)    = - (r(i,:)*Dr_dtbeta_p(i,:)    + z(i,:)*Dz_dtbeta_p(i,:)    + 3.d0*dtbeta_p(i,:))/rr(i,:)
              scomplex_piR(i,:) = - (r(i,:)*Dr_complex_piR(i,:) + z(i,:)*Dz_complex_piR(i,:) + complex_piR(i,:))/rr(i,:)
           end if

!          Radiative conditions at z boundaries.

           if (rank>=size-nprocr) then
              j = Nz
              sdtalpha(:,j)     = - (r(:,j)*Dr_dtalpha(:,j)     + z(:,j)*Dz_dtalpha(:,j)     + dtalpha(:,j))/rr(:,j)
              sKTA(:,j)         = - (r(:,j)*Dr_KTA(:,j)         + z(:,j)*Dz_KTA(:,j)         + KTA(:,j))/rr(:,j)
              sKTH(:,j)         = - (r(:,j)*Dr_KTH(:,j)         + z(:,j)*Dz_KTH(:,j)         + KTH(:,j))/rr(:,j)
              sAlambda(:,j)     = - (r(:,j)*Dr_Alambda(:,j)     + z(:,j)*Dz_Alambda(:,j)     + Alambda(:,j))/rr(:,j)
              sdtbeta_p(:,j)    = - (r(:,j)*Dr_dtbeta_p(:,j)    + z(:,j)*Dz_dtbeta_p(:,j)    + 3.d0*dtbeta_p(:,j))/rr(:,j)
              scomplex_piR(:,j) = - (r(:,j)*Dr_complex_piR(:,j) + z(:,j)*Dz_complex_piR(:,j) + complex_piR(:,j))/rr(:,j)
           end if

           if ((.not.eqsym).and.(rank<nprocr)) then
              j = 1-ghost
              sdtalpha(:,j)     = - (r(:,j)*Dr_dtalpha(:,j)     + z(:,j)*Dz_dtalpha(:,j)     + dtalpha(:,j))/rr(:,j)
              sKTA(:,j)         = - (r(:,j)*Dr_KTA(:,j)         + z(:,j)*Dz_KTA(:,j)         + KTA(:,j))/rr(:,j)
              sKTH(:,j)         = - (r(:,j)*Dr_KTH(:,j)         + z(:,j)*Dz_KTH(:,j)         + KTH(:,j))/rr(:,j)
              sAlambda(:,j)     = - (r(:,j)*Dr_Alambda(:,j)     + z(:,j)*Dz_Alambda(:,j)     + Alambda(:,j))/rr(:,j)
              sdtbeta_p(:,j)    = - (r(:,j)*Dr_dtbeta_p(:,j)    + z(:,j)*Dz_dtbeta_p(:,j)    + 3.d0*dtbeta_p(:,j))/rr(:,j)
              scomplex_piR(:,j) = - (r(:,j)*Dr_complex_piR(:,j) + z(:,j)*Dz_complex_piR(:,j) + complex_piR(:,j))/rr(:,j)
           end if

        end if


!       *****************************************************
!       ***   FOR RUNGE-KUTTA ADD TO ACCUMULATOR ARRAYS   ***
!       *****************************************************

!       Not yet implemented.

        if (rk4) then

        end if


!       ****************************
!       ***   UPDATE VARIABLES   ***
!       ****************************

        alpha   = alpha_p   + dtw*salpha
        dtalpha = dtalpha_p + dtw*sdtalpha

        A   = A_p   + dtw*sA
        KTA = KTA_p + dtw*sKTA

        H   = H_p   + dtw*sH
        KTH = KTH_p + dtw*sKTH

        lambda  = lambda_p  + dtw*slambda
        Alambda = Alambda_p + dtw*sAlambda

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

!          Boundaries for alpha.

           finevar   => alpha
           finevar_p => alpha_p
           finevar_bound_rR => alpha_bound_rR
           finevar_bound_rL => alpha_bound_rL
           finevar_bound_zR => alpha_bound_zR
           finevar_bound_zL => alpha_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%alpha
           else
              coarsevar => grid(box,level-1)%alpha
           end if

           call finebound(box,level,dtw,.false.)

!          Boundaries for dtalpha.

           finevar   => dtalpha
           finevar_p => dtalpha_p
           finevar_bound_rR => dtalpha_bound_rR
           finevar_bound_rL => dtalpha_bound_rL
           finevar_bound_zR => dtalpha_bound_zR
           finevar_bound_zL => dtalpha_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%dtalpha
           else
              coarsevar => grid(box,level-1)%dtalpha
           end if

           call finebound(box,level,dtw,.false.)

!          Boundaries for A.

           finevar   => A
           finevar_p => A_p
           finevar_bound_rR => A_bound_rR
           finevar_bound_rL => A_bound_rL
           finevar_bound_zR => A_bound_zR
           finevar_bound_zL => A_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%A
           else
              coarsevar => grid(box,level-1)%A
           end if

           call finebound(box,level,dtw,.false.)

!          Boundaries for KTA.

           finevar   => KTA
           finevar_p => KTA_p
           finevar_bound_rR => KTA_bound_rR
           finevar_bound_rL => KTA_bound_rL
           finevar_bound_zR => KTA_bound_zR
           finevar_bound_zL => KTA_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%KTA
           else
              coarsevar => grid(box,level-1)%KTA
           end if

           call finebound(box,level,dtw,.false.)

!          Boundaries for H.

           finevar   => H
           finevar_p => H_p
           finevar_bound_rR => H_bound_rR
           finevar_bound_rL => H_bound_rL
           finevar_bound_zR => H_bound_zR
           finevar_bound_zL => H_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%H
           else
              coarsevar => grid(box,level-1)%H
           end if

           call finebound(box,level,dtw,.false.)

!          Boundaries for KTH.

           finevar   => KTH
           finevar_p => KTH_p
           finevar_bound_rR => KTH_bound_rR
           finevar_bound_rL => KTH_bound_rL
           finevar_bound_zR => KTH_bound_zR
           finevar_bound_zL => KTH_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%KTH
           else
              coarsevar => grid(box,level-1)%KTH
           end if

           call finebound(box,level,dtw,.false.)

!          Boundaries for lambda.

           finevar   => lambda
           finevar_p => lambda_p
           finevar_bound_rR => lambda_bound_rR
           finevar_bound_rL => lambda_bound_rL
           finevar_bound_zR => lambda_bound_zR
           finevar_bound_zL => lambda_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%lambda
           else
              coarsevar => grid(box,level-1)%lambda
           end if

           call finebound(box,level,dtw,.false.)

!          Boundaries for Alambda.

           finevar   => Alambda
           finevar_p => Alambda_p
           finevar_bound_rR => Alambda_bound_rR
           finevar_bound_rL => Alambda_bound_rL
           finevar_bound_zR => Alambda_bound_zR
           finevar_bound_zL => Alambda_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%Alambda
           else
              coarsevar => grid(box,level-1)%Alambda
           end if

           call finebound(box,level,dtw,.false.)

!          Boundaries for beta_p.

           finevar   => beta_p
           finevar_p => beta_p_p
           finevar_bound_rR => beta_p_bound_rR
           finevar_bound_rL => beta_p_bound_rL
           finevar_bound_zR => beta_p_bound_zR
           finevar_bound_zL => beta_p_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%beta_p
           else
              coarsevar => grid(box,level-1)%beta_p
           end if

           call finebound(box,level,dtw,.false.)

!          Boundaries for dtbeta_p.

           finevar   => dtbeta_p
           finevar_p => dtbeta_p_p
           finevar_bound_rR => dtbeta_p_bound_rR
           finevar_bound_rL => dtbeta_p_bound_rL
           finevar_bound_zR => dtbeta_p_bound_zR
           finevar_bound_zL => dtbeta_p_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%dtbeta_p
           else
              coarsevar => grid(box,level-1)%dtbeta_p
           end if

           call finebound(box,level,dtw,.false.)

!          Boundaries for complex_phiR.

           finevar   => complex_phiR
           finevar_p => complex_phiR_p
           finevar_bound_rR => complex_phiR_bound_rR
           finevar_bound_rL => complex_phiR_bound_rL
           finevar_bound_zR => complex_phiR_bound_zR
           finevar_bound_zL => complex_phiR_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%complex_phiR
           else
              coarsevar => grid(box,level-1)%complex_phiR
           end if

           call finebound(box,level,dtw,.false.)

!          Boundaries for complex_piR.

           finevar   => complex_piR
           finevar_p => complex_piR_p
           finevar_bound_rR => complex_piR_bound_rR
           finevar_bound_rL => complex_piR_bound_rL
           finevar_bound_zR => complex_piR_bound_zR
           finevar_bound_zL => complex_piR_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%complex_piR
           else
              coarsevar => grid(box,level-1)%complex_piR
           end if

           call finebound(box,level,dtw,.false.)

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

              A(1-i,:)   = A(i,:)
              KTA(1-i,:) = KTA(i,:)

              H(1-i,:)   = H(i,:)
              KTH(1-i,:) = KTH(i,:)

              lambda(1-i,:)  = lambda(i,:)
              Alambda(1-i,:) = Alambda(i,:)

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

              A(:,1-j)   = A(:,j)
              KTA(:,1-j) = KTA(:,j)

              H(:,1-j)   = H(:,j)
              KTH(:,1-j) = KTH(:,j)

              lambda(:,1-j)  = lambda(:,j)
              Alambda(:,1-j) = Alambda_p(:,j)

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

           call sync(A)
           call sync(KTA)

           call sync(H)
           call sync(KTH)

           call sync(lambda)
           call sync(Alambda)

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
! current subroutine recursively.

  if (level<Nlmax) then
     call rotbosonstep(level+1,method)
     call rotbosonstep(level+1,method)
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

  if (level>0) then

!    Loop over all pairs of boxes in this level.

     do box=0,Nb

        if (Nl(box)<level) cycle

        do k=0,Nb

           if ((Nl(k)<level).or.(k==box)) cycle

           print *, 'Sync boxes not yet implemented'

        end do

     end do

  end if


! ****************************************************
! ***   RESTRICT FINE GRID DATA INTO COARSE GRID   ***
! ****************************************************

! Restrict the data from the fine to the coarse grid
! across all boxes when both levels coincide in time.
! This restriction does not change data in the current
! grid level, but rather in the coarser level.

  if (level>0) then
     do box=0,bmax

!       If this level does not exist for this box cycle.

        if (Nl(box)<level) cycle

!       We only restrict when the fine grid catches
!       up with the coarse grid.

        if (mod(s(box,level),2)==0) then

!          Figure out on which box is level-1.

           if (level==1) then
              bbox = 0
           else
              bbox = box
           end if

!          Restriction for alpha.

           finevar => grid(box,level)%alpha
           coarsevar => grid(bbox,level-1)%alpha
           call restrict(box,level,.false.)

!          Restriction for dtalpha.

           finevar => grid(box,level)%dtalpha
           coarsevar => grid(bbox,level-1)%dtalpha
           call restrict(box,level,.false.)

!          Restriction for A.

           finevar => grid(box,level)%A
           coarsevar => grid(bbox,level-1)%A
           call restrict(box,level,.false.)

!          Restriction for KTA.

           finevar => grid(box,level)%KTA
           coarsevar => grid(bbox,level-1)%KTA
           call restrict(box,level,.false.)

!          Restriction for H.

           finevar => grid(box,level)%H
           coarsevar => grid(bbox,level-1)%H
           call restrict(box,level,.false.)

!          Restriction for KTH.

           finevar => grid(box,level)%KTH
           coarsevar => grid(bbox,level-1)%KTH
           call restrict(box,level,.false.)

!          Restriction for lambda.

           finevar => grid(box,level)%lambda
           coarsevar => grid(bbox,level-1)%lambda
           call restrict(box,level,.false.)

!          Restriction for Alambda.

           finevar => grid(box,level)%Alambda
           coarsevar => grid(bbox,level-1)%Alambda
           call restrict(box,level,.false.)

!          Restriction for beta_p.

           finevar => grid(box,level)%beta_p
           coarsevar => grid(bbox,level-1)%beta_p
           call restrict(box,level,.false.)

!          Restriction for dtbeta_p.

           finevar => grid(box,level)%dtbeta_p
           coarsevar => grid(bbox,level-1)%dtbeta_p
           call restrict(box,level,.false.)

!          Restriction for complex_phiR.

           finevar => grid(box,level)%complex_phiR
           coarsevar => grid(bbox,level-1)%complex_phiR
           call restrict(box,level,.false.)

!          Restriction for complex_piR.

           finevar => grid(box,level)%complex_piR
           coarsevar => grid(bbox,level-1)%complex_piR
           call restrict(box,level,.false.)

!          Point to grid on level-1.

           call currentgrid(bbox,level-1,grid(bbox,level-1))

!          Symmetries.

           if (ownaxis) then
              do i=1,ghost

                 alpha(1-i,:)   = alpha(i,:)
                 dtalpha(1-i,:) = dtalpha(i,:)

                 A(1-i,:)   = A(i,:)
                 KTA(1-i,:) = KTA(i,:)

                 H(1-i,:)   = H(i,:)
                 KTH(1-i,:) = KTH(i,:)

                 lambda(1-i,:)  = lambda(i,:)
                 Alambda(1-i,:) = Alambda(i,:)

                 beta_p(1-i,:)   = beta_p(i,:)
                 dtbeta_p(1-i,:) = dtbeta_p(i,:)

                 complex_phiR(1-i,:) = complex_phiR(i,:)
                 complex_piR(1-i,:)  = complex_piR(i,:)

              end do
           end if

           if (eqsym.and.ownequator) then
              do j=1,ghost

                 alpha(:,1-j)   = alpha(:,j)
                 dtalpha(:,1-j) = dtalpha(:,j)

                 A(:,1-j)   = A(:,j)
                 KTA(:,1-j) = KTA(:,j)

                 H(:,1-j)   = H(:,j)
                 KTH(:,1-j) = KTH(:,j)

                 lambda(:,1-j)  = lambda(:,j)
                 Alambda(:,1-j) = Alambda(:,j)

                 beta_p(:,1-j)   = beta_p(:,j)
                 dtbeta_p(:,1-j) = dtbeta_p(:,j)

                 complex_phiR(:,1-j) = complex_phiR(:,j)
                 complex_piR(:,1-j)  = complex_piR(:,j)

              end do
           end if

!          Sync.

           if (size>1) then

              call sync(alpha)
              call sync(dtalpha)

              call sync(A)
              call sync(KTA)

              call sync(H)
              call sync(KTH)

              call sync(lambda)
              call sync(Alambda)

              call sync(beta_p)
              call sync(dtbeta_p)

              call sync(complex_phiR)
              call sync(complex_piR)

           end if

        end if

     end do
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine rotbosonstep
