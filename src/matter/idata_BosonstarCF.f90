
  subroutine idata_BosonstarCF

! *************************************************************
! ***   BOSON STAR INITIAL DATA IN CONFORMALLY FLAT GAUGE   ***
! *************************************************************

! Boson stars are solutions such that the spacetime is static
! and the complex scalar field has a harmonic dependence on time.
! This subroutine calculates initial data for a boson star
! in the conformally flat gauge.
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
! __2              ~ij                              4
! \/   alpha  +  2 g  d ln(psi) d alpha  -  4 pi psi (rho + trS)  =  0
!  flat                i         j
!
!
! This routine takes as input parameter the value of the scalar
! field at the origin "boson_phi0", and solves the coupled system
! of elliptic equations.
!
! We solve the system by turning the elliptic equations into
! hyperbolic wave equations and "evolving" to a steady state.
!
! This method is slow (specially at high resolution), but seems quite robust. 
! It can be improved by using a good initial guess (I'll try this later),
! and maybe something like multi-grid since it converges much faster

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
     print *, 'Solving initial data for a boson star in conformally flat gauge ...'
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

        complex_phiR = boson_phi0*exp(-rr**2)

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

  if (WE_verbose) then

     do box=0,Nb
        do level=min(1,box),Nl(box)

           call currentgrid(box,level,grid(box,level))

           grabvar => alpha
           call save1Dvariable('boson_alpha',directory,box,level,outparallel,'replace')
           !call save2Dvariable('boson_alpha',directory,box,level,outparallel,'replace')

           grabvar => phi
           call save1Dvariable('boson_psi',directory,box,level,outparallel,'replace')
           !call save2Dvariable('boson_psi',directory,box,level,outparallel,'replace')

           grabvar => complex_phiR
           call save1Dvariable('boson_phiR',directory,box,level,outparallel,'replace')
           !call save2Dvariable('boson_phiR',directory,box,level,outparallel,'replace')

           grabvar => dtalpha
           call save1Dvariable('boson_dtalpha',directory,box,level,outparallel,'replace')
           !call save2Dvariable('boson_dtalpha',directory,box,level,outparallel,'replace')

           grabvar => dtphi
           call save1Dvariable('boson_dtpsi',directory,box,level,outparallel,'replace')
           !call save2Dvariable('boson_dtpsi',directory,box,level,outparallel,'replace')

           grabvar => complex_piR
           call save1Dvariable('boson_piR',directory,box,level,outparallel,'replace')
           !call save2Dvariable('boson_piR',directory,box,level,outparallel,'replace')

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

     call bosonstep(0,method)


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
                + abs(grid(0,0)%sdtalpha(i,j)) + abs(grid(0,0)%sdtphi(i,j)) 
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
              call save1Dvariable('boson_alpha',directory,box,level,outparallel,'old')
              !call save2Dvariable('boson_alpha',directory,box,level,outparallel,'old')

              grabvar => phi
              call save1Dvariable('boson_psi',directory,box,level,outparallel,'old')
              !call save2Dvariable('boson_psi',directory,box,level,outparallel,'old')

              grabvar => complex_phiR
              call save1Dvariable('boson_phiR',directory,box,level,outparallel,'old')
              !call save2Dvariable('boson_phiR',directory,box,level,outparallel,'old')

              grabvar => dtalpha
              call save1Dvariable('boson_dtalpha',directory,box,level,outparallel,'old')
              !call save2Dvariable('boson_dtalpha',directory,box,level,outparallel,'old')

              grabvar => dtphi
              call save1Dvariable('boson_dtpsi',directory,box,level,outparallel,'old')
              !call save2Dvariable('boson_dtpsi',directory,box,level,outparallel,'old')

              grabvar => complex_piR
              call save1Dvariable('boson_piR',directory,box,level,outparallel,'old')
              !call save2Dvariable('boson_piR',directory,box,level,outparallel,'old')

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

        write (*,'(A,i6,A)') ' BosonStarCF: Iterations did not converge after ',WE_maxiter,' iterations.'
        print *

     else

        if (Nlmax_old>0) then

           if (Nlmax_old/=Nlmax) then
              write (*,'(A,i5,A)') ' BosonStarCF: Coarse grid solution converged after ',step,' iterations!'
              print *
           else
              write (*,'(A,i5,A)') ' BosonStarCF: Finer grids solution converged after ',step,' iterations!'
              print *
              write(*,'(A,ES23.16)') ' Final residual = ',gres
              write(*,'(A,ES23.16)') ' Omega          = ',boson_omega
              print *
           end if

        else

           write (*,'(A,i5,A)') ' BosonStarCF:   Solution converged after ',step,' iterations!'
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

  if ((step/=WE_maxiter).and.(Nlmax/=Nlmax_old)) then

!    Set again Nlmax to its original value.

     Nlmax = Nlmax_old

!    Set time derivatives in the base grid back to 0.

     grid(0,0)%dtalpha = 0.d0
     grid(0,0)%dtphi = 0.d0
     grid(0,0)%complex_piR = 0.d0

!    Iterate over boxes and levels higher than 0.

     do box=0,Nb
        do level=1,Nl(box)

!          Point to current grid.

           call currentgrid(box,level,grid(box,level))

           dtalpha = 0.d0
           dtphi = 0.d0
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

                 interpvar => grid(0,0)%phi
                 aux1 = interp(0,0,r0,z0,flag2)
                 call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

                 if (flag1) then
                    phi(i,j) = aux2
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


! ********************************
! ***   IMAGINARY PART OF pi   ***
! ********************************

! From the original ansatz, we take the time derivative of
! the imaginary part equal to (omega/alpha)*phi.

  do box=0,Nb
     do level=min(1,box),Nl(box)
        call currentgrid(box,level,grid(box,level))
        complex_piI = boson_omega*complex_phiR/alpha
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
!            /
! massTK  =  | alpha (rho + trS) dV
!            /
!
! with alpha the lapse, rho the energy density, trS
! the trace of the stress tensor, and dV the physical
! volume element:
!
! dV = 2 pi r sqrt(hdet) psi**6 dr dz
!
! The factor 2*pi comes from the integral over the angle.
!
! Notice that for a boson star we have:
!
! rho + trS  =  2 [ PiI**2 - V ]
!
! Remember that since we have conformally flat data hdet=1,
! and also that we have been using "phi" instead of "psi"
! (this is fixed at the end of this routine).
!
! At the moment I only do the integral in the coarse grid.

  call currentgrid(0,0,grid(0,0))

  call potential

  auxarray = (4.d0*smallpi*r*phi**6)*alpha*(complex_piI**2 - complex_V)
  mass_TK  = integral(0,0,auxarray)

  if (rank==0) then
     write (*,'(A,ES13.6)') ' Tolman-Komar mass (mass_TK) = ',mass_TK
     print *
  end if



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

!       Now calculate the correct "phi".

        phi = dlog(psi)

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

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Find spatial derivatives of complex_phiR.

        diffvar => complex_phiR
        complex_xiR_r = diff1r(+1)
        complex_xiR_z = diff1z(+1)

!       Set the imaginary part of the scalar field and its spatial derivatives to zero.

        complex_phiI  = 0.d0
        complex_xiI_r = 0.d0
        complex_xiI_z = 0.d0

!       Set time derivative of the real part of phi to zero.

        complex_piR  = 0.d0

     end do
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine idata_BosonstarCF













  recursive subroutine bosonstep(level,method)

! *********************************
! ***   ADVANCE ONE TIME STEP   ***
! *********************************

! Here we advance one time step in the internal evolution
! of the boson star solver.

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

!    Old values of conformal factor and time derivative.

     phi_p   = phi
     dtphi_p = dtphi

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

!          phi.

           phi_bound_rL(i,:,2) = phi_bound_rL(i,:,1)
           phi_bound_rL(i,:,1) = phi_bound_rL(i,:,0)
           phi_bound_rL(i,:,0) = phi(1-ghost+i,:)
           phi_bound_rR(i,:,2) = phi_bound_rR(i,:,1)
           phi_bound_rR(i,:,1) = phi_bound_rR(i,:,0)
           phi_bound_rR(i,:,0) = phi(Nr-i,:)

           phi_bound_zL(:,i,2) = phi_bound_zL(:,i,1)
           phi_bound_zL(:,i,1) = phi_bound_zL(:,i,0)
           phi_bound_zL(:,i,0) = phi(:,1-ghost+i)
           phi_bound_zR(:,i,2) = phi_bound_zR(:,i,1)
           phi_bound_zR(:,i,1) = phi_bound_zR(:,i,0)
           phi_bound_zR(:,i,0) = phi(:,Nz-i)

!          dtphi.

           dtphi_bound_rL(i,:,2) = dtphi_bound_rL(i,:,1)
           dtphi_bound_rL(i,:,1) = dtphi_bound_rL(i,:,0)
           dtphi_bound_rL(i,:,0) = dtphi(1-ghost+i,:)
           dtphi_bound_rR(i,:,2) = dtphi_bound_rR(i,:,1)
           dtphi_bound_rR(i,:,1) = dtphi_bound_rR(i,:,0)
           dtphi_bound_rR(i,:,0) = dtphi(Nr-i,:)

           dtphi_bound_zL(:,i,2) = dtphi_bound_zL(:,i,1)
           dtphi_bound_zL(:,i,1) = dtphi_bound_zL(:,i,0)
           dtphi_bound_zL(:,i,0) = dtphi(:,1-ghost+i)
           dtphi_bound_zR(:,i,2) = dtphi_bound_zR(:,i,1)
           dtphi_bound_zR(:,i,1) = dtphi_bound_zR(:,i,0)
           dtphi_bound_zR(:,i,0) = dtphi(:,Nz-i)

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

!       Derivatives of (phi,dtphi).

        diffvar => phi

        Dr_phi = diff1r(+1)
        Dz_phi = diff1z(+1)
     
        Drr_phi = diff2r(+1)
        Dzz_phi = diff2z(+1)

        diffvar => dtphi

        Dr_dtphi = diff1r(+1)
        Dz_dtphi = diff1z(+1)

!       Derivatives of (complex_phiR,complex_piR).

        diffvar => complex_phiR

        Dr_complex_phiR = diff1r(+1)
        Dz_complex_phiR = diff1z(+1)
     
        Drr_complex_phiR = diff2r(+1)
        Dzz_complex_phiR = diff2z(+1)

        diffvar => complex_piR

        Dr_complex_piR = diff1r(+1)
        Dz_complex_piR = diff1z(+1)

!       Scalar field potential.

        call potential

!       Find frequency omega.  We solve for omega from the
!       Klein-Gordon equation at the first grid point.

        if (level==Nlmax) then

           if (ownaxis.and.ownequator) then

              aux = alpha(1,1)**2/complex_phiR(1,1)*(complex_VPR(1,1) - 1.d0/phi(1,1)**4 &
                  *(Drr_complex_phiR(1,1) + Dzz_complex_phiR(1,1) + Dr_complex_phiR(1,1)/r(1,1) &
                  + Dr_complex_phiR(1,1)*(Dr_alpha(1,1)/alpha(1,1) + 2.d0*Dr_phi(1,1)/phi(1,1)) &
                  + Dz_complex_phiR(1,1)*(Dz_alpha(1,1)/alpha(1,1) + 2.d0*Dz_phi(1,1)/phi(1,1))))

              boson_omega = sqrt(abs(aux))

           end if

        end if

        call MPI_BCAST(boson_omega,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!       The sources for (alpha,phi,complex_phiR) are just (dtalpha,dtphi,complex_piR).

        salpha = dtalpha

        sphi = dtphi

        scomplex_phiR = complex_piR

!       The sources of (dtalpha,dtphi,complex_piR) all include the flat Laplacian.

        sdtalpha = Drr_alpha + Dzz_alpha + Dr_alpha/r

        sdtphi = Drr_phi + Dzz_phi + Dr_phi/r

        scomplex_piR = Drr_complex_phiR + Dzz_complex_phiR + Dr_complex_phiR/r

!       Extra sources for dtalpha.

        sdtalpha = sdtalpha + 2.d0*(Dr_alpha*Dr_phi/phi + Dz_alpha*Dz_phi/phi) &
                 - 8.d0*smallpi*alpha*phi**4*((boson_omega*complex_phiR/alpha)**2 - complex_V)

!       Extra sources for dtphi.

        sdtphi = sdtphi + smallpi*phi**5 &
               *((Dr_complex_phiR**2 + Dz_complex_phiR**2)/phi**4 &
               + (boson_omega*complex_phiR/alpha)**2 + 2.d0*complex_V)

!       Extra sources for complex_piR.

        scomplex_piR = scomplex_piR + phi**4*(complex_phiR*(boson_omega/alpha)**2 - complex_VPR) &
                     + Dr_complex_phiR*(Dr_alpha/alpha + 2.d0*Dr_phi/phi) &
                     + Dz_complex_phiR*(Dz_alpha/alpha + 2.d0*Dz_phi/phi)

!       Damping term.  This is needed in order to avoid large
!       oscillations and make the iterations stable.  But it
!       seems we only need it for the Klein-Gordon equation.

        scomplex_piR = scomplex_piR - WE_eta*complex_piR

!       And add some dissipation to reduce high frequency noise.

        evolvevar => complex_piR
        sourcevar => scomplex_piR
        call dissipation(+1,+1,WE_diss)

        evolvevar => dtalpha
        sourcevar => sdtalpha
        call dissipation(+1,+1,WE_diss)

        evolvevar => dtphi
        sourcevar => sdtphi
        call dissipation(+1,+1,WE_diss)

!       But set the source at point (1,1) such that the
!       cubic interpolated value to r=0 does not change.
!
!       In order to do this we notice that a cubic
!       interpolation to r=0 (assuming that phiR is
!       symmetric) implies that:
!
!       f0 = (9*f(1,1) - f(2,2))/8
!
!       from which we find:
!
!       f(1,1) = (8*f0 + f(2,2))/9
!
!       and its time derivative:
!
!       df(1,1)/dt = (df(2,2)/dt)/9

        if ((rank==0).and.(level==Nlmax)) then
           !scomplex_piR(1,1) = 0.d0
           scomplex_piR(1,1) = scomplex_piR(2,2)/9.d0
        end if

!       Symmetries on axis.

        if (ownaxis) then
           do i=1,ghost
              sdtalpha(1-i,:) = sdtalpha(i,:)
              sdtphi(1-i,:)   = sdtphi(i,:)
              scomplex_piR(1-i,:) = scomplex_piR(i,:)
           end do
        end if

!       Symmetries on equator.

        if (eqsym.and.ownequator) then
           do j=1,ghost
              sdtalpha(:,1-j) = sdtalpha(:,j)
              sdtphi(:,1-j)   = sdtphi(:,j)
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
              scomplex_piR(i,:) = - (r(i,:)*Dr_complex_piR(i,:) + z(i,:)*Dz_complex_piR(i,:) + complex_piR(i,:))/rr(i,:)
           end if

!          Radiative conditions at z boundaries.

           if (rank>=size-nprocr) then
              j = Nz
              sdtalpha(:,j)     = - (r(:,j)*Dr_dtalpha(:,j)     + z(:,j)*Dz_dtalpha(:,j)     + dtalpha(:,j))/rr(:,j)
              sdtphi(:,j)       = - (r(:,j)*Dr_dtphi(:,j)       + z(:,j)*Dz_dtphi(:,j)       + dtphi(:,j  ))/rr(:,j)
              scomplex_piR(:,j) = - (r(:,j)*Dr_complex_piR(:,j) + z(:,j)*Dz_complex_piR(:,j) + complex_piR(:,j))/rr(:,j)
           end if

           if ((.not.eqsym).and.(rank<nprocr)) then
              j = 1-ghost
              sdtalpha(:,j)     = - (r(:,j)*Dr_dtalpha(:,j) + z(:,j)*Dz_dtalpha(:,j) + dtalpha(:,j))/rr(:,j)
              sdtphi(:,j)       = - (r(:,j)*Dr_dtphi(:,j) + z(:,j)*Dz_dtphi(:,j) + dtphi(:,j))/rr(:,j)
              scomplex_piR(:,j) = - (r(:,j)*Dr_complex_piR(:,j) + z(:,j)*Dz_complex_piR(:,j) + complex_piR(:,j))/rr(:,j)
           end if

        end if


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

        alpha   = alpha_p   + dtw*salpha
        dtalpha = dtalpha_p + dtw*sdtalpha

        phi   = phi_p   + dtw*sphi
        dtphi = dtphi_p + dtw*sdtphi

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

!          Boundaries for phi.

           finevar   => phi
           finevar_p => phi_p
           finevar_bound_rR => phi_bound_rR
           finevar_bound_rL => phi_bound_rL
           finevar_bound_zR => phi_bound_zR
           finevar_bound_zL => phi_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%phi
           else
              coarsevar => grid(box,level-1)%phi
           end if

           call finebound(box,level,dtw,.false.)

!          Boundaries for dtphi.

           finevar   => dtphi
           finevar_p => dtphi_p
           finevar_bound_rR => dtphi_bound_rR
           finevar_bound_rL => dtphi_bound_rL
           finevar_bound_zR => dtphi_bound_zR
           finevar_bound_zL => dtphi_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%dtphi
           else
              coarsevar => grid(box,level-1)%dtphi
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

              phi(1-i,:)   = phi(i,:)
              dtphi(1-i,:) = dtphi(i,:)

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

  if (level>0) then

!    Loop over all pairs of boxes in this level.

     do box=0,Nb

        if (Nl(box)<level) cycle

        do k=0,Nb

           if ((Nl(k)<level).or.(k==box)) cycle

           box1var => grid(box,level)%alpha
           box2var => grid(k  ,level)%alpha
           call syncboxes(box,k,level,.false.)

           box1var => grid(box,level)%dtalpha
           box2var => grid(k  ,level)%dtalpha
           call syncboxes(box,k,level,.false.)

           box1var => grid(box,level)%phi
           box2var => grid(k  ,level)%phi
           call syncboxes(box,k,level,.false.)

           box1var => grid(box,level)%dtphi
           box2var => grid(k  ,level)%dtphi
           call syncboxes(box,k,level,.false.)

           box1var => grid(box,level)%complex_phiR
           box2var => grid(k  ,level)%complex_phiR
           call syncboxes(box,k,level,.false.)

           box1var => grid(box,level)%complex_piR
           box2var => grid(k  ,level)%complex_piR
           call syncboxes(box,k,level,.false.)

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

!          Restriction for phi.

           finevar => grid(box,level)%phi
           coarsevar => grid(bbox,level-1)%phi
           call restrict(box,level,.false.)

!          Restriction for dtphi.

           finevar => grid(box,level)%dtphi
           coarsevar => grid(bbox,level-1)%dtphi
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

                 phi(1-i,:)   = phi(i,:)
                 dtphi(1-i,:) = dtphi(i,:)

                 complex_phiR(1-i,:) = complex_phiR(i,:)
                 complex_piR(1-i,:)  = complex_piR(i,:)

              end do
           end if

           if (eqsym.and.ownequator) then
              do j=1,ghost

                 alpha(:,1-j)   = alpha(:,j)
                 dtalpha(:,1-j) = dtalpha(:,j)

                 phi(:,1-j)   = phi(:,j)
                 dtphi(:,1-j) = dtphi(:,j)

                 complex_phiR(:,1-j) = complex_phiR(:,j)
                 complex_piR(:,1-j)  = complex_piR(:,j)

              end do
           end if

!          Sync.

           if (size>1) then

              call sync(alpha)
              call sync(dtalpha)

              call sync(phi)
              call sync(dtphi)

              call sync(complex_phiR)
              call sync(complex_piR)

           end if

        end if

     end do
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine bosonstep
