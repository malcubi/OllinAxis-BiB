!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/elliptic/wave_elliptic.f90,v 1.35 2021/02/26 18:35:59 malcubi Exp $

  subroutine wave_elliptic(type,init)

! ****************************************************
! ***   ELLIPTIC SOLVER VIA HYPERBOLIC EVOLUTION   ***
! ****************************************************

! This is a simple elliptic solver for the Poisson equation:
!
!  __2                   5
!  \/ u  +  S1 u  +  S5 u   =  S0
!
!
! Here S0 is a source term, S1 is the coefficient of the 
! linear term, and S5 the coefficient of the u^5 term
! (this term often appears in the Hamiltonian constraint).
!
! The equation is solved by transforming it into:
!
!   2      2    __2                  5
!  d u / dt  =  \/ u  +  S1 u  + S5 u  -  S0  -  k du/dt
!
!
! and looking for stationary solutions. Here k is a damping
! coefficient.
!
! The main routine just calls recursively a simplified version
! of src/base/onestep.f90.

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
  integer i,j,n,m                ! Counters.
  integer step                   ! Iteration counter.
  integer Nlmax_old

  real(8) lres,gres              ! Local and global residuals.
  real(8) r0,z0,interp           ! For interpolation.
  real(8) vrmax,vzmax            ! Local maximum speeds.
  real(8) waveeta                ! Damping parameter.
  real(8) cfac                   ! Courant parameter.
  real(8) zero,one               ! Numbers.
  real(8) aux1,aux2

  character(*) type              ! Type of Laplacian (flat,conformal,physical).
  character(*) init              ! Initial guess.
  character(3) method            ! Time integration method.

  integer s_ext(0:Nb,0:Nlmax)    ! External values of step counter.
  real(8) t_ext(0:Nb,0:Nlmax)    ! External values of time counter.
  real(8) t1_ext(0:Nb,0:Nlmax)   ! External values of t1 counter.
  real(8) t2_ext(0:Nb,0:Nlmax)   ! External values of t2 counter.

  real(8) vr_max(0:Nb,0:Nlmax),vz_max(0:Nb,0:Nlmax) ! Arrays with maximum speeds for different boxes and levels.


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  one  = 1.d0


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


! *******************************
! ***   INITIALIZE SOLUTION   ***
! *******************************

! We always initialize ell_v to 0.
!
! By default we initialize ell_u to 1, but we can
! overrdide this depending on the value of 'init'.

  do box=0,Nb
     do level=min(1,box),Nl(box)

        call currentgrid(box,level,grid(box,level))

        ell_v = zero

        if (init=='zero') then
           ell_u = zero
        else if (init=='one') then
           ell_u = one
        else
           if (rank==0) then
              print *
              print *, 'Unknow initial data for wave_elliptic'
              print *, 'Setting to 1.'
              print *
           end if
           ell_u = one
        end if

!       Derivatives of u.

        diffvar => ell_u

        Dr_ell_u = diff1r(+1)
        Dz_ell_u = diff1z(+1)
     
        Drr_ell_u = diff2r(+1)
        Dzz_ell_u = diff2z(+1)
        Drz_ell_u = diff2rz(+1,+1)

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

           grabvar => ell_u
           call save1Dvariable('ell_u',directory,box,level,outparallel,'replace')
           call save2Dvariable('ell_u',directory,box,level,outparallel,'replace')

           grabvar => ell_v
           call save1Dvariable('ell_v',directory,box,level,outparallel,'replace')
           call save2Dvariable('ell_v',directory,box,level,outparallel,'replace')

        end do
     end do

  end if


! *****************************
! ***   COURANT PARAMETER   ***
! *****************************

! Internal time integration method.  Since we evolve to
! stationary solution, it doesn't really matter if the
! time evolution is second (icn) or fourth (rk4) order.
! By default I set rk4.

  method = "rk4"

! Courant parameter:  For rk4 we set the Courant
! parameter to 0.7, and for icn to 0.5 (it goes
! unstable otherwise).

  if (method=="icn") then
     cfac = 0.5d0
  else if (method=="rk4") then
     cfac = 0.7d0
  end if

! Choose time step.

  if (type=='flat') then

     dt0 = cfac*min(dr0,dz0)

  else

! For a curved Laplacian we need to adapt the Courant
! parameter to the physical metric.
!
! First we find the maximum speeds in each box and level.

     do box=0,Nb
        do level=min(1,box),Nl(box)

!          Point to current grid.

           call currentgrid(box,level,grid(box,level))

!          Wave speeds.

           vl_rp = sqrt(abs(g_A/psi4))
           vl_zp = sqrt(abs(g_B/psi4))

!          Find maximum speeds in r and z.

           vrmax = 0.d0
           vzmax = 0.d0

           do j=1-ghost,Nz
              do i=1-ghost,Nr

                 if (abs(vl_rp(i,j))>vrmax) vrmax=abs(vl_rp(i,j))
                 if (abs(vl_zp(i,j))>vzmax) vzmax=abs(vl_zp(i,j))

                 if (abs(va_rp(i,j))>vrmax) vrmax=abs(va_rp(i,j))
                 if (abs(va_zp(i,j))>vzmax) vzmax=abs(va_zp(i,j))

              end do
           end do

!          For parallel runs find the global maximum.

           if (size>1) then

              call MPI_ALLREDUCE(vrmax,aux1,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
              vrmax = aux1

              call MPI_ALLREDUCE(vzmax,aux1,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
              vzmax = aux1

           end if

!          Save maximum speeds for current box and level.

           vr_max(box,level) = vrmax
           vz_max(box,level) = vzmax

        end do
     end do

!    Now find maximum speeds across boxes and levels.

     vrmax = 0.d0
     vzmax = 0.d0

     do box=0,Nb
        do level=min(1,box),Nl(box)
           if (abs(vr_max(box,level))>vrmax) vrmax=abs(vr_max(box,level))
           if (abs(vz_max(box,level))>vzmax) vzmax=abs(vz_max(box,level))
        end do
     end do

!    Modify time step.

     dt0 = cfac*min(dr0,dz0,dr0/vrmax,dz0/vzmax)

  end if

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

! Initialize damping parameter.

  waveeta = WE_eta

! Initialize residuals and iteration number.

  100 continue

  lres = 1.d0
  gres = 1.d0

  step = 0

! Start iterations.  Notice that we do at least 100 iterations.
! This is because since we start with ell_v=0, the first few
! time steps the change in ell_u might be very small.

  do while ((step<100).or.((gres>ELL_epsilon).and.(step<ELL_maxiter)))


!    ******************************************
!    ***   ADVANCE ONE INTERNAL TIME STEP   ***
!    ******************************************

!    Increment step counter.

     step = step + 1

!    Advance one time step.

     call wavestep(0,type,waveeta,method)


!    *************************
!    ***   FIND RESIDUAL   ***
!    *************************

!    Find local residual.

     lres = 0.d0
     gres = 0.d0

     do j=1,Nzl(0,rank)-ghost
        do i=1,Nrl(0,rank)-ghost
           lres = lres + abs(grid(0,0)%ell_u(i,j)-grid(0,0)%ell_u_p(i,j))
        end do
     end do

!    Find global residual.

     call MPI_Allreduce(lres,gres,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)


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

              grabvar => ell_u
              call save1Dvariable('ell_u',directory,box,level,outparallel,'old')
              call save2Dvariable('ell_u',directory,box,level,outparallel,'old')

              grabvar => ell_v
              call save1Dvariable('ell_v',directory,box,level,outparallel,'old')
              call save2Dvariable('ell_v',directory,box,level,outparallel,'old')

           end do
        end do

     end if

  end do


! ****************************
! ***   DID WE CONVERGE?   ***
! ****************************

  if ((rank==0).and.(ELL_verbose)) then

     if (step==ELL_maxiter) then
        write (*,'(A,i6,A)') ' WaveElliptic:   Iterations did not converge after ',ELL_maxiter,' iterations.'
     else
        if (Nlmax_old>0) then
           if (Nlmax_old/=Nlmax) then
              write (*,'(A,i5,A)') ' WaveElliptic:   Coarse grid solution converged after ',step,' iterations!'
           else
              write (*,'(A,i5,A)') ' WaveElliptic:   Finer grids solution converged after ',step,' iterations!'
           end if
        else
           write (*,'(A,i5,A)') ' WaveElliptic:   Solution converged after ',step,' iterations!'
        end if
     end if

  end if


! *****************************
! ***   REFINEMENT LEVELS   ***
! *****************************

! If we converged and we have refinement levels, inject the
! coarse solution into fine grids and solve again.

  if ((step/=ELL_maxiter).and.(Nlmax/=Nlmax_old)) then

!    Set damping parameter to zero. This parameter is required
!    to avoid large oscilations when we start far from the
!    true solution.  But once we get here we have a good solution
!    on the coarse level, so we don't need it anymore and it
!    only slows the code.

     waveeta = 0.d0

!    Set again Nlmax to its original value.

     Nlmax = Nlmax_old

!    Set ell_v in the base grid back to 0.

     grid(0,0)%ell_v = 0.d0

!    Iterate over boxes and levels higher than 0.

     do box=0,Nb
        do level=1,Nl(box)

!          Point to current grid.

           call currentgrid(box,level,grid(box,level))

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

!                Interpolate variable from coarse grid level.

                 interpvar => grid(0,0)%ell_u
                 aux1 = interp(0,0,r0,z0,flag2)
                 call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!                But we only update the value if the location belongs to us.

                 if (flag1) then
                    ell_u(i,j) = aux2
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


! ***************
! ***   END   ***
! ***************

  end subroutine wave_elliptic









  recursive subroutine wavestep(level,type,waveeta,method)

! ***************************************************
! ***   EVOLUTION OF WAVE EQUATION WITH SOURCES   ***
! ***************************************************

! This routine does one time step of a wave equation
! with sources in first order form:
!
! du/dt  =  v
!           __2                 5
! dv/dt  =  \/ u  +  S1 u  +  S5 u  -  S0  -  k v
!
! This routine is is just a simplified version of src/base/onestep.f90.

! Include modules.

  use param
  use arrays
  use derivatives
  use procinfo

! Extra variables.

  implicit none

  integer box,level    ! Box number and level counters.
  integer i,j,k        ! Counters.
  integer iter         ! Counter for internal iterations.
  integer niter        ! Number of internal iterations.
  integer bmax         ! Number of boxes at this level.
  integer bbox

  real(8) dtw          ! Internal time step.
  real(8) weight       ! Weight for rk4.
  real(8) waveeta      ! Damping parameter.

  character(*) type    ! Type of Laplacian (flat,conformal,physical)
  character(*) method  ! Time integration method.


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

!    Save old value of (u,v).

     ell_u_p = ell_u
     ell_v_p = ell_v

!    Save old values of boundaries.

     if (level>0) then

        do i=0,ghost-1
           ell_u_bound_rL(i,:,2) = ell_u_bound_rL(i,:,1)
           ell_u_bound_rL(i,:,1) = ell_u_bound_rL(i,:,0)
           ell_u_bound_rL(i,:,0) = ell_u(1-ghost+i,:)
           ell_u_bound_rR(i,:,2) = ell_u_bound_rR(i,:,1)
           ell_u_bound_rR(i,:,1) = ell_u_bound_rR(i,:,0)
           ell_u_bound_rR(i,:,0) = ell_u(Nr-i,:)
           ell_u_bound_zL(:,i,2) = ell_u_bound_zL(:,i,1)
           ell_u_bound_zL(:,i,1) = ell_u_bound_zL(:,i,0)
           ell_u_bound_zL(:,i,0) = ell_u(:,1-ghost+i)
           ell_u_bound_zR(:,i,2) = ell_u_bound_zR(:,i,1)
           ell_u_bound_zR(:,i,1) = ell_u_bound_zR(:,i,0)
           ell_u_bound_zR(:,i,0) = ell_u(:,Nz-i)
        end do

        do i=0,ghost-1
           ell_v_bound_rL(i,:,2) = ell_v_bound_rL(i,:,1)
           ell_v_bound_rL(i,:,1) = ell_v_bound_rL(i,:,0)
           ell_v_bound_rL(i,:,0) = ell_v(1-ghost+i,:)
           ell_v_bound_rR(i,:,2) = ell_v_bound_rR(i,:,1)
           ell_v_bound_rR(i,:,1) = ell_v_bound_rR(i,:,0)
           ell_v_bound_rR(i,:,0) = ell_v(Nr-i,:)
           ell_v_bound_zL(:,i,2) = ell_v_bound_zL(:,i,1)
           ell_v_bound_zL(:,i,1) = ell_v_bound_zL(:,i,0)
           ell_v_bound_zL(:,i,0) = ell_v(:,1-ghost+i)
           ell_v_bound_zR(:,i,2) = ell_v_bound_zR(:,i,1)
           ell_v_bound_zR(:,i,1) = ell_v_bound_zR(:,i,0)
           ell_v_bound_zR(:,i,0) = ell_v(:,Nz-i)
        end do

     end if


!    **********************************************
!    ***   FIND NUMBER OF INTERNAL ITERATIONS   ***
!    **********************************************

     if (method=="icn") then
        niter = 3
     else if (method=="rk4") then
        niter = 4
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

!       Iterative Crank-Nicholson (ICN).

        if (method=="icn") then

!          In ICN all iterations except the last one
!          jump only half a time step.

           if (iter<niter) then
              dtw = 0.5d0*dt
           else
              dtw = dt
           end if

!       Fourth order Runge-Kutta.

        else if (method=="rk4") then

!          In fourth order Runge-Kutta the first two iterations
!          jump half a time step and the last two a full time step.
!          Here we also set the weights with which intermediate
!          results contribute to final answer: 1/6 for first and
!          last intermediate results and 1/3 for the two middle ones.

           if (iter==1) then
              dtw = 0.5d0*dt
              weight = 1.d0/6.d0
           else if (iter==2) then
              dtw = 0.5d0*dt
              weight = 1.d0/3.d0
           else if (iter==3) then
              dtw = dt
              weight = 1.d0/3.d0
           else
              dtw = dt
              weight = 1.d0/6.d0
           end if

        end if


!       ************************
!       ***   FIND SOURCES   ***
!       ************************

!       Derivatives of u.

        diffvar => ell_u

        Dr_ell_u = diff1r(+1)
        Dz_ell_u = diff1z(+1)
     
        Drr_ell_u = diff2r(+1)
        Dzz_ell_u = diff2z(+1)
        Drz_ell_u = diff2rz(+1,+1)

!       Derivatives of v.

        diffvar => ell_v

        Dr_ell_v = diff1r(+1)
        Dz_ell_v = diff1z(+1)

!       The source for u is just v.

        sell_u = ell_v

!       Source for v:  Laplacian term.

        if (type=='flat') then

!          Flat Laplacian.

           sell_v = Drr_ell_u + Dzz_ell_u + Dr_ell_u/r

        else if ((type=='physical').or.(type=='conformal')) then

!          Physical Laplacian: second derivative terms.

           sell_v = g_A*Drr_ell_u + g_B*Dzz_ell_u + 2.d0*r*g_C*Drz_ell_u

!          Physical laplacian: derivatives of conformal metric.

           sell_v = sell_v + Dr_ell_u*(Dr_g_A + r*Dz_g_C) &
                           + Dz_ell_u*(Dz_g_B + r*Dr_g_C + g_C)

!          Physical laplacian: derivatives of determinant of conformal metric.

           sell_v = sell_v + g_A*Dr_ell_u/r + g_C*Dz_ell_u + 0.5d0*ihdet &
                  *(Dr_hdet*(g_A*Dr_ell_u + r*g_C*Dz_ell_u) &
                  + Dz_hdet*(g_B*Dz_ell_u + r*g_C*Dr_ell_u))

!          Physical Laplacian: derivatives of conformal factor.

           if (type/='conformal') then
              sell_v = sell_v + 2.d0*(Dr_phi*(g_A*Dr_ell_u + r*g_C*Dz_ell_u) &
                                    + Dz_phi*(g_B*Dz_ell_u + r*g_C*Dr_ell_u))
           end if

!          Rescale by conformal factor since all the above terms
!          were calculated with the conformal metric.

           if (type/='conformal') then
              sell_v = sell_v/psi4
           end if

        else

!          Unknown Laplacian type.

           if (rank==0) then
              print *
              print *, 'Unknow Laplacian type'
              print *, 'Aborting! (Subroutine wave_elliptic.f90)'
              print *
           end if

           call die

        end if

!       Source for v: source terms

        sell_v = sell_v + ell_S1*ell_u + ell_S5*ell_u**5 - ell_S0

!       Source for v: damping term.

        sell_v = sell_v - waveeta*ell_v/dtl(level)

!       And add some dissipation.  This is mostly because when
!       we have refinement levels the boundaries of the fine
!       grids introduce high frequency reflections that we
!       need to dissipate away.

        if (Nlmax>0) then
           evolvevar => ell_v
           sourcevar => sell_v
           call dissipation(+1,+1,0.01d0)
        end if

!       Symmetries on axis.

        if (ownaxis) then
           do i=1,ghost
              sell_u(1-i,:) = sell_u(i,:)
              sell_v(1-i,:) = sell_v(i,:)
           end do
        end if

!       Symmetries on equator.

        if (eqsym.and.ownequator) then
           do j=1,ghost
              sell_u(:,1-j) = sell_u(:,j)
              sell_v(:,1-j) = sell_v(:,j)
           end do
        end if


!       *******************************
!       ***   BOUNDARY CONDITIONS   ***
!       *******************************

!       At the moment I use radiative boundaries that reduce
!       to a Robin boundary condition for static solutions.
!       Newmann and Diritchlet boundaries can be added later.

        evolvevar => ell_v
        sourcevar => sell_v
        Dr_var => Dr_ell_v
        Dz_var => Dz_ell_v
        call radbound(1.d0,0.d0)


!       *****************************************************
!       ***   FOR RUNGE-KUTTA ADD TO ACCUMULATOR ARRAYS   ***
!       *****************************************************

        if (method=="rk4") then
           if (iter==1) then
              ell_u_a = weight*sell_u
              ell_v_a = weight*sell_v
           else if (iter<niter) then
              ell_u_a = ell_u_a + weight*sell_u
              ell_v_a = ell_v_a + weight*sell_v
           else
              sell_u  = ell_u_a + weight*sell_u
              sell_v  = ell_v_a + weight*sell_v
           end if
        end if


!       ****************************
!       ***   UPDATE VARIABLES   ***
!       ****************************

        ell_u = ell_u_p + dtw*sell_u
        ell_v = ell_v_p + dtw*sell_v


!       *************************************************
!       ***   FOR FINE GRIDS INTERPOLATE BOUNDARIES   ***
!       *************************************************

!       For fine grids we need to interpolate from the new
!       time level of the coarse grid to get boundary data.
!
!       Remember that the coarse grid has already advanced
!       to the next time level.

        if (level>0) then

!          Boundaries for ell_u.

           finevar   => ell_u
           finevar_p => ell_u_p
           finevar_bound_rR => ell_u_bound_rR
           finevar_bound_rL => ell_u_bound_rL
           finevar_bound_zR => ell_u_bound_zR
           finevar_bound_zL => ell_u_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%ell_u
           else
              coarsevar => grid(box,level-1)%ell_u
           end if

           call finebound(box,level,dtw,.false.)

!          Boundaries for ell_v.

           finevar   => ell_v
           finevar_p => ell_v_p
           finevar_bound_rR => ell_v_bound_rR
           finevar_bound_rL => ell_v_bound_rL
           finevar_bound_zR => ell_v_bound_zR
           finevar_bound_zL => ell_v_bound_zL

           if (level==1) then
              coarsevar => grid(0,level-1)%ell_v
           else
              coarsevar => grid(box,level-1)%ell_v
           end if

           call finebound(box,level,dtw,.false.)

        end if


!       *****************************************
!       ***   SYNCHRONIZE ACROSS PROCESSORS   ***
!       *****************************************

!       If we have more than one processor we must now
!       synchronize ghost zones.

        if (size>1) then
           call sync(ell_u)
           call sync(ell_v)
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
     call wavestep(level+1,type,waveeta,method)
     call wavestep(level+1,type,waveeta,method)
  end if


! **********************
! ***   SYMMETRIES   ***
! **********************

! Apply symmetries again since we might have messed
! them up when we where on higher refinement levels
! and restricted data to the current level.

  if (Nlmax>0) then
     do box=0,bmax

!       If this level does not exist for this box cycle.

        if (Nl(box)<level) cycle

!       Do we own axis and/or equator?

        ownaxis = (axis(box,rank)/=-1)
        ownequator = (eqz(box,rank)/=-1)

!       Symmetries on axis.

        if (ownaxis) then
           do i=1,ghost
              grid(box,level)%ell_u(1-i,:) = grid(box,level)%ell_u(i,:)
              grid(box,level)%ell_v(1-i,:) = grid(box,level)%ell_v(i,:)
           end do
        end if

!       Symmetries on equator.

        if (eqsym.and.ownequator) then
           do j=1,ghost
              grid(box,level)%ell_u(:,1-j) = grid(box,level)%ell_u(:,j)
              grid(box,level)%ell_v(:,1-j) = grid(box,level)%ell_v(:,j)
           end do
        end if

     end do
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

           box1var => grid(box,level)%ell_u
           box2var => grid(k  ,level)%ell_u
           call syncboxes(box,k,level,.false.)

           box1var => grid(box,level)%ell_v
           box2var => grid(k  ,level)%ell_v
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

!          Restriction for ell_u.

           finevar => grid(box,level)%ell_u
           coarsevar => grid(bbox,level-1)%ell_u
           call restrict(box,level,.false.)

!          Restriction for ell_v.

           finevar => grid(box,level)%ell_v
           coarsevar => grid(bbox,level-1)%ell_v

           call restrict(box,level,.false.)

!          Point to grid on level-1.

           call currentgrid(bbox,level-1,grid(bbox,level-1))


           if (ownaxis) then
              do i=1,ghost
                 sell_u(1-i,:) = sell_u(i,:)
                 sell_v(1-i,:) = sell_v(i,:)
              end do
           end if

           if (eqsym.and.ownequator) then
              do j=1,ghost
                 sell_u(:,1-j) = sell_u(:,j)
                 sell_v(:,1-j) = sell_v(:,j)
              end do
           end if

           if (size>1) then
              call sync(ell_u)
              call sync(ell_v)
           end if

        end if

     end do
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine wavestep

