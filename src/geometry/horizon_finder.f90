
  subroutine horizon_finder

! ***********************************
! ***   APPARENT HORIZON FINDER   ***
! ***********************************

! Routine to find apparent horizons in axisymmetry.
! In order to find the apparent horizon we need to
! solve an ODE (see subroutines below).
!
! At the moment we allow to search for up to three
! apparent horizons (0,1,2), all centered on the axis
! (horizon 0 is always centered on the origin).
!
! The routine is not really parallel.  All processors
! do the exact same calculation, but only proc 0 does
! the output.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical interpflag         ! Interpolation flag.
  logical bracket            ! Did we bracket the horizon center?

  integer k                  ! Counter.
  integer box,level,bb,ll    ! Box and level counters.

  real(8) r0,z0,z0_old       ! Interpolation point.
  real(8) interp             ! Interpolation function.
  real(8) Dzphi,Dzphi_old    ! Interpolated value of Dz_phi.
  real(8) aux                ! Auxiliary.
  real(8) zah(0:2)           ! Array with horizon centers.


! ************************
! ***   SANITY CHECK   ***
! ************************

! Sanity check.

  if (N_ah>3) then
     if (rank==0) then
        print *
        print *, 'The value of N_ah can not be larger than 3.'
        print *, 'Aborting! (subroutine horizon_finder)'
        print *
     end if
     call die
  end if

! If we have equatorial symmetry we can only look for
! two horizons.

  if (eqsym.and.(N_ah==3)) then
     if (rank==0) then
        print *
        print *, 'For equatorial symmetry we can only look for 2 horizons.'
        print *, 'Aborting! (subroutine horizon_finder)'
        print *
     end if
     call die
  end if

! By construction, horizon 0 is centered on the origin,
! while horizons 1 and 2 must be above and below the
! equator respectively.

  if (z_ah1<0.0d0) then
     if (rank==0) then
        print *
        print *, 'Horizon 1 must be above the equator.'
        print *, 'Aborting! (subroutine horizon_finder)'
        print *
     end if
     call die
  end if

  if (z_ah2>0.d0) then
     if (rank==0) then
        print *
        print *, 'Horizon 2 must be below the equator.'
        print *, 'Aborting! (subroutine horizon_finder)'
        print *
     end if
     call die
  end if


! ********************************
! ***   FIND HORIZON CENTERS   ***
! ********************************

! Horizon 0 is always centered on the origin.

  zah(0) = 0.d0

! Initialize positions of horizons 1 and 2.

  if (N_ah>1) then
     zah(1) = z_ah1
     zah(2) = z_ah2
  end if

! If the centers do move, locate the closest
! maximum of the conformal factor to z_ah.

  if ((ahmove).and.(N_ah>1)) then

!    Initialize r0.

     r0 = 0.d0

!    Loop over horizons.

     do k=1,N_ah-1

!       Initialize z0.

        z0 = zah(k)

!       We need to interpolate at the highest level
!       available at the current location.

        box = 0
        level = 0

        if (Nlmax>0) then
           do bb=0,Nb
              do ll=1,Nl(bb)
                 if ((r0>rminl(bb,ll)+drl(ll)).and.(r0<rmaxl(bb,ll)-drl(ll)).and. &
                     (z0>zminl(bb,ll)+dzl(ll)).and.(z0<zmaxl(bb,ll)-dzl(ll))) then
                    box = bb
                    level = ll
                 end if
              end do
           end do
        end if

!       Find local maximum of phi.

        bracket = .false.

        Dzphi     = 0.d0
        Dzphi_old = 0.d0

        do while (.not.bracket)

!          Find value of Dz_phi at current location.

           interpvar => grid(box,level)%Dz_phi
           aux = interp(box,level,r0,z0,interpflag)
           call MPI_ALLREDUCE(aux,Dzphi,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!          Now move along axis to find local maximum of phi.

           if ((.not.bracket).and.(Dzphi*Dzphi_old>=0.d0)) then
              z0_old = z0
              Dzphi_old = Dzphi
              if (Dzphi>0.d0) then
                 z0 = z0 + dz0
              else
                 z0 = z0 - dz0
              end if
           else if (.not.bracket) then
              bracket = .true.
              zah(k) = z0_old - Dzphi_old*(z0-z0_old)/(Dzphi-Dzphi_old)
           end if

        end do

     end do

  end if


! *********************************
! ***   NOW LOOK FOR HORIZONS   ***
! *********************************

! Call routine that integrates horizon equation.
! There are two methods at the moment: a shooting
! method and specral fast flow method.

  do k=0,N_ah-1

     if (ahverbose.and.(rank==0)) then
        if (N_ah==1) then
           print *
           print *, 'Looking for apparent horizon'
        else
           print *
           print *, 'Looking for apparent horizon centered on:',zah(k)
        end if
     end if

     if (ahmethod=="shoot") then
        call ahshoot(k,zah(k))
     else if (ahmethod=="spectral") then
        call ahspectral(k,zah(k))
     end if

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine horizon_finder








  subroutine ahshoot(k,zcenter)

! ****************************************************************
! ***   FIND APPARENT HORIZON EQUATION USING SHOOTING METHOD   ***
! ****************************************************************

! This subroutine solves the apparent horizon ODE
! in axisymmetry using a shooting method and a
! bisection refinement to obtain the correct initial
! condition.
!
! Originally written by Erik Jimenez for the old code.
! It has been adapted to the new code.
!
! This version of the integrate_ah routine uses an
! adaptive stepsize fourth order Runge-Kutta. The
! adaptive parte is quite simple, but it seems to
! do the job.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  common /ahfinder/ ack,ahfound

  logical ahfound(0:2)        ! Did we find a horizon?
  logical firstcall(0:2)      ! Is this the first call?
  logical bracket             ! Did we bracket the horizon?
  logical firstiter           ! First iteration?
  logical errorflag           ! Variable for error control.
  logical errorflag_old       ! Error flag on last iteration.

  integer i,j                 ! Counters.
  integer k                   ! Horizon counter (0,1,2).
  integer Ntheta,Nthetamax    ! Number of points in theta, and maximum allowed.
  integer box,level,bb,ll     ! Box and level counters.
  integer iaux

  real(8) zcenter             ! Center of horizon.
  real(8) dzz,dth,dth2        ! z and theta intervals.
  real(8) thetafinal          ! Final value of theta.
  real(8) dzlocal             ! Local grid size.
  real(8) zinit,zinit_old     ! Initial z position for shooting method.
  real(8) za,zb,drada,dradb   ! Values of z and drad for bracketing.
  real(8) dradlast            ! Derivative dr/dtheta at last point.
  real(8) dradlast_old        ! Old value of dradlast.
  real(8) k_1,k_2,k_3,k_4     ! Runge-Kutta steps for rad.
  real(8) l_1,l_2,l_3,l_4     ! Runge-Kutta steps for dhdth.
  real(8) m_1,m_2,m_3,m_4     ! Runge-Kutta steps for horizon area.
  real(8) n_1,n_2,n_3,n_4     ! Runge-Kutta steps for angular momentum.
  real(8) sdh                 ! Source for horizon function.
  real(8) sarea,sJmom         ! Sources for area and angular momentum.
  real(8) ah_area             ! Area of horizon.
  real(8) ah_J                ! Angular momentum of horizon.
  real(8) ah_echarge          ! Charge of horizon.
  real(8) MAH,MI,RAH          ! Apparent horizon mass, irreducible mass, and areal radius.
  real(8) eps                 ! Small number.
  real(8) zero,one,two,half,smallpi
  real(8) rad_aux,drad_aux,rad_error,drad_error

  real(8), dimension(0:2,0:50) :: ack  ! Fourier coefficients for each horizon.

  real(8), allocatable, dimension(:) :: rad,theta   ! Coordinates of horizon as arrays.
  real(8), allocatable, dimension(:) :: drad        ! Array for derivative dr/dtheta.
  real(8), allocatable, dimension(:) :: area        ! Array for area integration.
  real(8), allocatable, dimension(:) :: Jmom        ! Array for angular momentum integration.

  character(4)  filen
  character(20) filestatus

  data firstcall /3*.true./
  save firstcall


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  one  = 1.d0
  two  = 2.0d0
  half = 0.5d0

  smallpi = acos(-1.d0)


! **********************
! ***   INITIALIZE   ***
! **********************

! Give Nthetamax.

  Nthetamax = 50000

! Allocate arrays (rad,theta).

  allocate(rad(0:Nthetamax),theta(0:Nthetamax),drad(0:Nthetamax))

! Allocate arrays for area and angular momentum.

  allocate(area(0:Nthetamax),Jmom(0:Nthetamax))


! **************************************
! ***   INITIALIZE SHOOTING METHOD   ***
! **************************************

! Sanity check.

  if (ahrmax>0.9d0) then
      if (rank==0) then
        print *, 'Value of ahrmax is too large, setting it to 0.9'
        print *
     end if
     ahrmax = 0.9d0
  end if

! Initial z position for shooting.

  zinit_old = zero

  if (ahzinit==0.d0) then
     if (zcenter>=zero) then
        zinit = 0.8d0*zmaxl(0,0)
     else
        zinit = 0.8d0*zminl(0,0)
     end if
  else if ((ahzinit>0.d0).and.(ahzinit<=0.9d0*zmaxl(0,0))) then
     if (zcenter>=zero) then
        zinit = + ahzinit
     else
        zinit = - ahzinit
     end if
  else
     if (rank==0) then
        print *, 'Value of ahzinit out of bounds, setting it to 80% of coarsest grid size'
        print *
     end if
     if (zcenter>=zero) then
        zinit = 0.8d0*zmaxl(0,0)
     else
        zinit = 0.8d0*zminl(0,0)
     end if
  end if

! Set initial value of dzlocal.

  dzlocal = dzl(0)


! ***************************
! ***   MOVE ALONG AXIS   ***
! ***************************

! Open file for iterations data (for testing).

  if (ahsaveiter) then

     if (N_ah==1) then
        write(filen,'(A)') '.tl'
     else
        write(filen,'(i1,A)') k,'.tl'
     end if

     open(1,file=trim(directory)//'/ah_iter'//trim(filen),form='formatted',status='replace')

  end if

! The boundary conditions are dh/dtheta(0)=0 and dh/dtheta(pi)=0.
! Only in case of equatorial symmetry we can stop the integration
! at theta=pi/2 and fix the boundary condition as dh/dtheta(pi/2)=0.
!
! The shooting method uses a 4th Runge-Kutta integration.
! If we obtain a sign change in the derivative at the last
! point we continue with a bisection refinement.

  50 continue

! Flags.

  bracket   = .false.
  firstiter = .true.
  errorflag = .false.

! Initialize (za,zb) to clearly different things.
! (they are used in the bisection method below).

  za = +1.d10
  zb = -1.d10

! Move along axis. We iterate until we bracket a horizon,
! or until we get too close to the horizon center.
! The number of iterations is saved in the counter 'j'.

  j = 0

  do while ((abs(zinit-zcenter)>5.d0*dzlocal).and.(abs(zinit_old-zinit)>ahepsilon_shoot))

     errorflag_old = errorflag

!    Find local highest resolution on axis,
!    and modify (dzlocal,dzz) accordingly.

     box = 0
     level = 0

     if (Nlmax>0) then
        do bb=0,Nb
           do ll=1,Nl(bb)
              if ((0.d0>rminl(bb,ll)+drl(ll)).and.(0.d0<rmaxl(bb,ll)-drl(ll)).and. &
                  (zinit>zminl(bb,ll)+dzl(ll)).and.(zinit<zmaxl(bb,ll)-dzl(ll))) then
                 box = bb
                 level = ll
              end if
           end do
        end do
     end if

     dzlocal = dzl(level)
     dzz = 0.5d0*dzlocal

!    Initialize rad and drad.

     rad  = 0.d0
     drad = 0.d0

     rad(0) = abs(zinit - zcenter)

!    Store the old value of dradlast. On first
!    iteration set it to something large.

     if (.not.firstiter) then
        dradlast_old = dradlast
     else
        dradlast_old = 1.d10
     end if

!    Initialize area and angular momentum.

     area(0) = zero
     Jmom(0) = zero


!    **************************************
!    ***   INTEGRATE HORIZON EQUATION   ***
!    **************************************

!    Initialize i.

     i = 0

!    Initialize theta and dth.  Notice that we
!    start a little bit away from the axis since
!    the horizon source term is singular there.

     dth = smallpi/100.d0

     if (zcenter>=zero) then
        theta(0) = 1.d-3
     else
        dth = - dth
        theta(0) = smallpi - 1.d-3
     end if

!    Find final value of theta. Again, for the cases
!    where we should integrate all the way to the
!    axis we stop a little bit away.

     if (eqsym.and.(zcenter==zero)) then
        thetafinal = half*smallpi
     else if (zcenter>=zero) then
        thetafinal = smallpi - 1.d-3
     else
        thetafinal = 1.d-3
     end if

!    Loop while we reach 'thetafinal'.

     do while (abs(theta(i)-thetafinal)>1.d-3)

!       Make sure we don't overshoot.

        if (zcenter>=zero) then
           if (theta(i)+dth>thetafinal) dth = thetafinal - theta(i)
        else
           if (theta(i)+dth<thetafinal) dth = thetafinal - theta(i)
        end if

        dth2 = half*dth

!       First take one large Rung-Kutta step.

        call ah_source(rad(i),drad(i),theta(i),zcenter,sdh,sarea,sJmom,errorflag)
        if (errorflag) goto 100
        k_1 = dth*drad(i)
        l_1 = dth*sdh
        m_1 = abs(dth2)*sarea
        n_1 = abs(dth2)*sJmom

        call ah_source(rad(i)+0.5d0*k_1,drad(i)+0.5d0*l_1,theta(i)+0.5d0*dth,zcenter,sdh,sarea,sJmom,errorflag)
        if (errorflag) goto 100
        k_2 = dth*(drad(i)+0.5d0*l_1)
        l_2 = dth*sdh

        call ah_source(rad(i)+0.5d0*k_2,drad(i)+0.5d0*l_2,theta(i)+0.5d0*dth,zcenter,sdh,sarea,sJmom,errorflag)
        if (errorflag) goto 100
        k_3 = dth*(drad(i)+0.5d0*l_2)
        l_3 = dth*sdh

        call ah_source(rad(i)+k_3,drad(i)+l_3,theta(i)+dth,zcenter,sdh,sarea,sJmom,errorflag)
        if (errorflag) goto 100
        k_4 = dth*(drad(i)+l_3)
        l_4 = dth*sdh

        rad_aux  = rad(i)  + 1.d0/6.d0*(k_1 + 2.d0*k_2 + 2.d0*k_3 + k_4)
        drad_aux = drad(i) + 1.d0/6.d0*(l_1 + 2.d0*l_2 + 2.d0*l_3 + l_4)

!       Take first small Runge-Kutta step. But notice that we don't need
!       to repeat the call to calculate (l1,k1,m1,n1) since it 

        !call ah_source(rad(i),drad(i),theta(i),zcenter,sdh,sarea,sJmom,errorflag)
        !if (errorflag) goto 100
        !k_1 = dth2*drad(i)
        !l_1 = dth2*sdh
        !m_1 = abs(dth2)*sarea
        !n_1 = abs(dth2)*sJmom

        k_1 = half*k_1
        l_1 = half*l_1

        call ah_source(rad(i)+0.5d0*k_1,drad(i)+0.5d0*l_1,theta(i)+0.5d0*dth2,zcenter,sdh,sarea,sJmom,errorflag)
        if (errorflag) goto 100
        k_2 = dth2*(drad(i)+0.5d0*l_1)
        l_2 = dth2*sdh
        m_2 = abs(dth2)*sarea
        n_2 = abs(dth2)*sJmom

        call ah_source(rad(i)+0.5d0*k_2,drad(i)+0.5d0*l_2,theta(i)+0.5d0*dth2,zcenter,sdh,sarea,sJmom,errorflag)
        if (errorflag) goto 100
        k_3 = dth2*(drad(i)+0.5d0*l_2)
        l_3 = dth2*sdh
        m_3 = abs(dth2)*sarea
        n_3 = abs(dth2)*sJmom

        call ah_source(rad(i)+k_3,drad(i)+l_3,theta(i)+dth2,zcenter,sdh,sarea,sJmom,errorflag)
        if (errorflag) goto 100
        k_4 = dth2*(drad(i)+l_3)
        l_4 = dth2*sdh
        m_4 = abs(dth2)*sarea
        n_4 = abs(dth2)*sJmom

        rad(i+1)  = rad(i)  + 1.d0/6.d0*(k_1 + 2.d0*k_2 + 2.d0*k_3 + k_4)
        drad(i+1) = drad(i) + 1.d0/6.d0*(l_1 + 2.d0*l_2 + 2.d0*l_3 + l_4)
        area(i+1) = area(i) + 1.d0/6.d0*(m_1 + 2.d0*m_2 + 2.d0*m_3 + m_4)
        Jmom(i+1) = Jmom(i) + 1.d0/6.d0*(n_1 + 2.d0*n_2 + 2.d0*n_3 + n_4)

!       Take second small Runge-Kutta step.

        call ah_source(rad(i+1),drad(i+1),theta(i)+dth,zcenter,sdh,sarea,sJmom,errorflag)
        if (errorflag) goto 100
        k_1 = dth2*drad(i+1)
        l_1 = dth2*sdh
        m_1 = abs(dth2)*sarea
        n_1 = abs(dth2)*sJmom

        call ah_source(rad(i+1)+0.5d0*k_1,drad(i+1)+0.5d0*l_1,theta(i)+1.5d0*dth2,zcenter,sdh,sarea,sJmom,errorflag)
        if (errorflag) goto 100
        k_2 = dth2*(drad(i+1)+0.5d0*l_1)
        l_2 = dth2*sdh
        m_2 = abs(dth2)*sarea
        n_2 = abs(dth2)*sJmom

        call ah_source(rad(i+1)+0.5d0*k_2,drad(i+1)+0.5d0*l_2,theta(i)+1.5d0*dth2,zcenter,sdh,sarea,sJmom,errorflag)
        if (errorflag) goto 100
        k_3 = dth2*(drad(i+1)+0.5d0*l_2)
        l_3 = dth2*sdh
        m_3 = abs(dth2)*sarea
        n_3 = abs(dth2)*sJmom

        call ah_source(rad(i+1)+k_3,drad(i+1)+l_3,theta(i)+dth,zcenter,sdh,sarea,sJmom,errorflag)
        if (errorflag) goto 100
        k_4 = dth2*(drad(i+1)+l_3)
        l_4 = dth2*sdh
        m_4 = abs(dth2)*sarea
        n_4 = abs(dth2)*sJmom

        rad(i+1)  = rad(i+1)  + 1.d0/6.d0*(k_1 + 2.d0*k_2 + 2.d0*k_3 + k_4)
        drad(i+1) = drad(i+1) + 1.d0/6.d0*(l_1 + 2.d0*l_2 + 2.d0*l_3 + l_4)
        area(i+1) = area(i+1) + 1.d0/6.d0*(m_1 + 2.d0*m_2 + 2.d0*m_3 + m_4)
        Jmom(i+1) = Jmom(i+1) + 1.d0/6.d0*(n_1 + 2.d0*n_2 + 2.d0*n_3 + n_4)

!       Find difference between one large step and two small ones.

        rad_error  = abs(rad(i+1)-rad_aux)
        drad_error = abs(drad(i+1)-drad_aux)

!       Now check if the errors are too large, in which case
!       we reduce dth by half and ignore this step.

        eps = 1.d-4*(abs(drad(i+1))+abs(drad_aux))

        if (drad_error>max(eps,1.d-8)) then

           dth = half*dth

!       If the error was small enough we accept this step.

        else

!          Advance theta and i (in that order).

           theta(i+1) = theta(i) + dth
           i = i+1

!          Local extrapolation.  We now have two fourth order
!          approximations, one with a large time step and one
!          with two small time steps.  We can use them to
!          obtain a fifth order extrapolation.

           rad(i+1)  = rad(i+1)  + (rad(i+1)-rad_aux)/15.d0
           drad(i+1) = drad(i+1) + (drad(i+1)-drad_aux)/15.d0

!          Check that 'i' has not exceeded Nthetamax.
!          If it did give up.

           if (i>=Nthetamax) then
              !print *, 'Nthetamax exceeded'
              errorflag = .true.
              goto 100
           end if

!          If the error is too small, we can safely increase dth.

           if (drad_error<0.1d0*eps) then
              dth = two*dth
           end if

        end if

     end do

!    Once we leave the loop, set Ntheta to the last value of i,
!    set dradlast=drad(Ntheta), and increase the iteration counter j.

     100 continue

     Ntheta = i 
     dradlast = drad(Ntheta)

     j = j+1

!    Save iterations data if needed (for testing).

     if (ahsaveiter.and.(mod(j,10)==0)) then

        !print *, j,theta(Ntheta),zinit,zinit_old,dradlast,dradlast_old

        write(1,"(A,ES20.10)") '"Time   = ',dble(j)
        write(1,"(A,ES20.10)") '"Center = ',zcenter

        do i=0,Ntheta
           write(1,*) theta(i),rad(i)
        end do

        write(1,*)

     end if


!    ***********************************************
!    ****   CHECK IF WE HAVE BRACKETED HORIZON   ***
!    ***********************************************

!    When we get here we have integrated the horizon equation
!    starting at some point along the axis. If this is not the
!    first time we do it, we can compare the value of the derivative
!    at the last point with the previous integration to see if
!    we have bracketed a horizon.

!    If we still have not bracketed a horizon move along the axis.

     if (firstiter.or. &
        ((.not.bracket).and.(dradlast*dradlast_old>=zero))) then

        if (firstiter) firstiter = .false.

        zinit_old = zinit

        if (zcenter>=zero) then
           zinit = zinit - dzz
        else
           zinit = zinit + dzz
        end if

!    If we found a sign change in dradlast then set the
!    'bracket' flag to true, and start the bisections.

     else if (.not.bracket) then

        bracket = .true.

        za = zinit_old
        zb = zinit

        drada = dradlast_old
        dradb = dradlast

        zinit_old = zinit
        zinit = half*(za+zb)

        if ((rank==0).and.ahverbose) then
           print *, 'horizon bracketed between z values: ', za,zb
           print *, 'with derivatives: ',drada,dradb
        end if

!    If we have already bracketed the horizon,
!    then continue with the bisections.

     else

        if (drada*dradlast<=zero) then
           zb = zinit
           dradb = dradlast
        else
           za = zinit
           drada = dradlast
        end if

        zinit_old = zinit
        zinit = half*(za+zb)

        !print *,za,zb,drada,dradb

     end if

  end do


! ************************************
! ***   HAVE WE FOUND A HORIZON?   ***
! ************************************

! When we get here we are out of the loop that moves
! along the axis. This can be because one of two
! conditions was satisfied:

! 1) The bisections converged.

  if (abs(zinit-zinit_old)<=ahepsilon_shoot) then

!    If there was no error in either this step or the last,
!    and the derivative is small, we have found a horizon.

     if ((.not.errorflag).and.(abs(dradlast)<1.d-1)) then

        ahfound(k) = .true.

        !print *, j,zinit,theta(Ntheta),dradlast

        if ((rank==0).and.ahverbose) then
           print *, 'Apparent horizon found with initial position in axis: ',zinit
           print *
        end if

!    If there was an error, the bisections converged to something
!    that is not really a horizon. Check now if we still have some
!    way to go along the axis. In that case move further along and
!    restart the search.

     else if (abs(zinit-zcenter)>dzlocal) then

        if ((rank==0).and.ahverbose) then
           print *, 'Bracketed a possible solution but it was not a horizon.'
           print *, 'errorflag =', errorflag
           print *, 'za,zb,drada,dradb =',za,zb,drada,dradb
           print *, 'Restarting the search ...'
        end if

        if (zcenter>=zero) then
           zinit = zinit - 2.d0*dzz
        else
           zinit = zinit + 2.d0*dzz
        end if

        goto 50

     end if

! 2) We finished moving along the axis and found no horizon.

  else

     ahfound(k) = .false.

     if ((rank==0).and.ahverbose) then
        print *, 'No horizon found.'
        print *
     end if

  end if

! Close file for iterations data (for testing).

  if (ahsaveiter) then
     close(1)
  end if


! ********************************
! ***   FOURIER COEFFICIENTS   ***
! ********************************

! If we found a horizon, calculate the Fourier (cosine)
! coefficients for its radius.  This is useful if we
! want to use this information later as an initial
! guess for the spectral finder (otherwise we can
! just ignore this).
!
! Notice that, since the shooting method uses adaptive
! steps, the integral is done over unevenly spaced
! points.  I use the trapezoidal rule, which is not very
! accurate.  For this reason I only use a few Fourier
! modes, as higher modes are very inaccurate. This shouldn't
! be a problem as I only need a "reasonable" initial guess
! for the spectral finder.

  if (ahfound(k)) then

     if (eqsym.and.(zcenter==0.d0)) then
        iaux = 2
     else
        iaux = 1
     end if

!    Loop over Fourier modes.

     ack(k,:) = zero

     do j=0,20,iaux

!       Integral over cosines.

        do i=0,Ntheta-1
           ack(k,j) = ack(k,j) + half*(theta(i+1) - theta(i)) &
                 *(rad(i+1)*cos(dble(j)*theta(i+1)) + rad(i)*cos(dble(j)*theta(i)))
        end do

!       Rescale.

        ack(k,j) = two*dble(iaux)*ack(k,j)/smallpi

        !if (rank==0) print *, "ac(",j,") =",ac(k,j)

     end do

  end if


! ************************
! ***   OUTPUT FILES   ***
! ************************

! Determine file status and extension.

  if (firstcall(k)) then
     firstcall(k) = .false.
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

  if (N_ah==1) then
     write(filen,'(A)') '.tl'
  else
     write(filen,'(i1,A)') k,'.tl'
  end if


! *********************************
! ***   SAVE HORIZON POSITION   ***
! *********************************

! Save horizon position.

  if (rank==0) then

     if (filestatus=='replace') then
        open(1,file=trim(directory)//'/ah_radius'//trim(filen),form='formatted', &
        status=filestatus)
     else
        open(1,file=trim(directory)//'/ah_radius'//trim(filen),form='formatted', &
        status=filestatus,position='append')
     end if

!    Write current time, and horizon center.

     write(1,"(A,ES20.8)") '"Time   = ',t(0,0)
     write(1,"(A,ES20.8)") '"Center = ',zcenter

!    Save horizon position.

     if (ahfound(k)) then
        do i=0,Ntheta
           write(1,"(2ES20.8E3)") theta(i),rad(i)
        end do
     else
        do i=0,Ntheta
           write(1,"(2ES20.8E3)") theta(i),0.d0
        end do
     end if

!    Leave blank space before next time.

     write(1,*)

!    Close file.

     close(1)

  end if


! **************************************
! ***   SAVE HORIZON AREA AND MASS   ***
! **************************************

  if (rank==0) then

!    Save horizon area and horizon areal radius defined as:
!
!    R = sqrt(Area/pi) / 2

     if (.not.ahfound(k)) then
        ah_area = zero
     else
        ah_area = area(Ntheta)
        if (eqsym.and.(zcenter==0.d0)) ah_area = two*ah_area
     end if

     RAH = half*sqrt(ah_area/smallpi)

     if (ahverbose.and.ahfound(k)) then
        print *, 'Horizon area:',ah_area
        print *, 'Areal radius:',RAH
     end if

     if (filestatus=='replace') then
        open(1,file=trim(directory)//'/ah_area'//trim(filen),form='formatted', &
        status=filestatus)
     else
        open(1,file=trim(directory)//'/ah_area'//trim(filen),form='formatted', &
        status=filestatus,position='append')
     end if

     write(1,"(2ES20.8E3)") t(0,0),ah_area

     close(1)

!    Save horizon angular momentum.

     if (angmom) then

        if (.not.ahfound(k)) then
           ah_J = zero
        else
           ah_J = Jmom(Ntheta)
           if (eqsym.and.(zcenter==0.d0)) ah_J = two*ah_J
        end if

        if (ahverbose.and.ahfound(k)) then
           print *, 'Horizon angular momentum:',ah_J
        end if

        if (filestatus == 'replace') then
           open(1,file=trim(directory)//'/ah_J'//trim(filen),form='formatted', &
           status=filestatus)
        else
           open(1,file=trim(directory)//'/ah_J'//trim(filen),form='formatted', &
           status=filestatus,position='append')
        end if

        write(1,"(2ES20.8E3)") t(0,0),ah_J

        close(1)

     else

        ah_J = zero

     end if

!    Find and save horizon mass.
!
!    The irreducible mass is then defined as:
!
!    MI  =  R/2  =  sqrt(Area/pi) / 4
!
!    with R the horizon areal radius.
!
!    The total mass, including the angular momentum
!    contribution is (still missing the charge
!    contribution):
!
!                 4       2
!    M  =  sqrt( R  +  4 J ) / 2 R
!
!    with J the horizon angular momentum.

     if (ahfound(k)) then
        MI = half*RAH
        if (angmom) then
           MAH = half*sqrt(RAH**4 + 4.d0*ah_J**2)/RAH
        else
           MAH = MI
        end if
     else
        MI  = zero
        MAH = zero
     end if

     if (ahverbose.and.ahfound(k)) then
        print *, 'Irreducible horizon mass:',MI
        print *, 'Total horizon mass:      ',MAH
        print *
     end if

     if (filestatus == 'replace') then
        open(1,file=trim(directory)//'/ah_mass'//trim(filen),form='formatted', &
        status=filestatus)
     else
        open(1,file=trim(directory)//'/ah_mass'//trim(filen),form='formatted', &
        status=filestatus,position='append')
     end if

     write(1,"(2ES20.8E3)") t(0,0),MAH

     close(1)

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine ahshoot









  subroutine ahspectral(k,zcenter)

! ****************************************************************
! ***   FIND APPARENT HORIZON EQUATION USING SPECTRAL METHOD   ***
! ****************************************************************

! For this method we use a cosine expansion in order to satisfy
! the correct boundary conditions (zero derivative at both theta=0
! and theta=pi):
!
!                              k_max
! h(theta)  =  A_0 / 2   +  sum       [ A_k  cos (k theta) ]
!                              k=1
!
! Remember that we parametrize the horizon as the surface r=h(theta)
! (but notice that the code uses "rad" for the function h).
!
! Notice that for the case of equatorial symmetry we must have zero
! derivative at theta=pi/2 also, which implies that we should only
! use even values of n in the above expansion.
!
! We then use an iterative procedure to find the expansion coefficients.
! For this we use the "fast flow" algorithm of Gundlach	(Phys.Rev.D57-863,
! arXiv:gr-qc/9707050), but modified for axisymmetry.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  common /ahfinder/ ack,ahfound

  logical ahfound(0:2)        ! Did we find a horizon?
  logical firstcall(0:2)      ! Is this the first call?
  logical errorflag           ! Variable for error control.

  integer i,j,q,iter          ! Counters.
  integer step                ! 1 for eqsym=.false., 2 for eqsym=.true.
  integer k                   ! Horizon counter (0,1,2).
  integer Ntheta              ! Number of points in theta (power of 2).
  integer Nmode               ! Number of Fourier modes.
  integer maxiter             ! Maximum number of iterations.

  real(8) zcenter             ! Center of horizon.
  real(8) dtheta              ! Interval in theta.
  real(8) res                 ! Residual.
  real(8) AA,BB               ! Parameters for flow algorithm.
  real(8) sdh                 ! Source for horizon function.
  real(8) sarea,sJmom         ! Source for area and angular momentum integrals.
  real(8) ah_area             ! Area of horizon.
  real(8) ah_J                ! Angular momentum of horizon.
  real(8) ah_echarge          ! Charge of horizon.
  real(8) MAH,MI,RAH          ! Apparent horizon mass, irreducible mass, and areal radius.
  real(8) zero,one,two,half,smallpi
  real(8) aux1,aux2

  real(8), dimension(0:50) :: ac,cc         ! Arrays for Fourier coefficients.
  real(8), dimension(0:50) :: ac_old,cc_old ! Old values of arrays.
  real(8), dimension(0:50) :: tc            ! Test step.
  real(8), dimension(0:2,0:50) :: ack       ! Fourier coefficients for each horizon.

  real(8), allocatable, dimension(:) :: theta          ! Array for angular coordinate.
  real(8), allocatable, dimension(:) :: rad,drad       ! Horizon radius and derivative.
  real(8), allocatable, dimension(:) :: rad_p,drad_p   ! Previous values of rad and drad.
  real(8), allocatable, dimension(:) :: Ifunc          ! Integrand function (see below).
  real(8), allocatable, dimension(:) :: area           ! Array for area integration.
  real(8), allocatable, dimension(:) :: Jmom           ! Array for angular momentum integration.

  character(4)  filen
  character(20) filestatus

  data firstcall /3*.true./

  save firstcall


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  one  = 1.d0
  two  = 2.0d0
  half = 0.5d0

  smallpi = acos(-1.d0)


! ***************************
! ***   ALLOCATE ARRAYS   ***
! ***************************

! Maximum number of iterations.

  maxiter = ahmaxiter

! Give Nmode and Ntheta. Notice that Nmode
! is copied from the global parameter "ahNmodes".

  if (ahNmodes<=20) then
     Nmode = ahNmodes
  else
     if (rank==0) then
         print *, 'Maximum number of Fourier modes allowed is 20.'
     end if
     Nmode = 20
  end if

  Nmode = ahNmodes

  Ntheta = 10*Nmode

! Allocate arrays (theta,rad,drad,Ifunc).

  allocate(theta(0:Ntheta))
  allocate(rad(0:Ntheta),drad(0:Ntheta),rad_p(0:Ntheta),drad_p(0:Ntheta))
  allocate(Ifunc(0:Ntheta))

! Allocate arrays for area and angular momentum.

  allocate(area(0:Ntheta),Jmom(0:Ntheta))


! **********************
! ***   INITIALIZE   ***
! **********************

! Find "step" and "dtheta" (step=2 for
! equatorial symmetry and step=1 otherwise).

  if (eqsym.and.(zcenter==0.d0)) then
     step = 2
  else
     step = 1
  end if

  dtheta = smallpi/dble(step*Ntheta)

! Initialize parameters AA and BB.

  AA = ahAA
  BB = ahBB

! Initialize theta.

  do i=0,Ntheta
     theta(i) = dble(i)*dtheta
  end do

! Sanity check.

  if (ahrmax>0.9d0) then
      if (rank==0) then
        print *, 'Value of ahrmax is too large, setting it to 0.9'
        print *
     end if
     ahrmax = 0.9d0
  end if

! Initialize Fourier coefficients to zero.

  ac = zero

  if (firstcall(k)) then
     ack(k,:) = zero
  end if


! *******************************************
! ***   CALL SHOOTING FOR INITIAL GUESS   ***
! *******************************************

  if (firstcall(k).and.ahshootfirst) then
     call ahshoot(k,zcenter)
     ack(k,0) = 1.05d0*ack(k,0)   ! We increase the size a bit for the flow.
     !do j=0,Nmode,step
     !   if (rank==0) print *, j,ack(k,j)
     !end do
  end if


! ********************************************************
! ***   CHOOSE INITIAL SURFACE AS LAST HORIZON FOUND   ***
! ********************************************************

! If we found a previous horizon we use it as the initial guess.

  if (ahfound(k)) then

     if ((rank==0).and.ahverbose) then
        print *, 'Using previous horizon as initial guess'
     end if

     ac(:) = ack(k,:)

     do i=0,Ntheta
        rad(i) = half*ac(0)
        drad(i) = zero
        do j=step,Nmode,step
           rad(i)  = rad(i) + ac(j)*cos(dble(j)*theta(i))
           drad(i) = drad(i) - dble(j)*ac(j)*sin(dble(j)*theta(i))
        end do
     end do

     goto 50

  end if


! **************************************************
! ***   CHOOSE INITIAL SURFACE AS LARGE CIRCLE   ***
! **************************************************

! Initialize horizon function rad to a large
! circle, and derivative drad to zero.

  if (zcenter==0.d0) then
     if (ahzinit==0.d0) then
        ac(0) = 1.8d0*ahrmax*zmaxl(0,0)
     else
        ac(0) = 2.d0*ahzinit
     end if
  else
     ac(0) = 1.95d0*abs(zcenter)
  end if

  rad = half*ac(0)
  drad = zero


! ************************************************
! ***   CHOOSE INITIAL SURFACE AS AN ELLIPSE   ***
! ************************************************

! If aheccentric is not zero we initialize the
! surface to an ellipse.  This is only done
! if the horizon center is not on the origin.

  if ((zcenter/=0.d0).and.(aheccentric/=0.d0)) then

!    Initialize rad to an ellipse.

     do i=0,Ntheta
        rad(i) = half*ac(0)/sqrt(one - (aheccentric*sin(theta(i)))**2)
     end do

!    Find Fourier coefficients.

     do j=0,Nmode

!       First point.

        ac(j) = half*dtheta*rad(0)

!       Interior points.

        do i=1,Ntheta-1
           ac(j) = ac(j) + dtheta*rad(i)*cos(dble(j)*theta(i))
        end do

!       Last point.

        ac(j) = ac(j) + half*dtheta*rad(Ntheta)*cos(dble(j)*theta(Ntheta))

!       Normalization.

        ac(j) = two*dble(step)*ac(j)/smallpi
        !if (rank==0) print *, j,ac(j)

     end do

!    Reconstruct rad and its derivative
!    from cosine expansion.

     do i=0,Ntheta
        rad(i) = half*ac(0)
        drad(i) = zero
        do j=step,Nmode,step
           rad(i)  = rad(i) + ac(j)*cos(dble(j)*theta(i))
           drad(i) = drad(i) - dble(j)*ac(j)*sin(dble(j)*theta(i))
        end do
     end do

  end if


! ****************************
! ***   START ITERATIONS   ***
! ****************************

  50 continue

! Open file for iterations data and save initial data (for testing).

  if ((rank==0).and.(ahsaveiter)) then

     if (N_ah==1) then
        write(filen,'(A)') '.tl'
     else
        write(filen,'(i1,A)') k,'.tl'
     end if

     open(1,file=trim(directory)//'/ah_iter'//trim(filen),form='formatted',status='replace')

     write(1,"(A,i5)") '"Time   = ',0
     write(1,"(A,ES20.10)") '"Center = ',zcenter

     do i=0,Ntheta
        write(1,*) theta(i),rad(i)
     end do

     write(1,*)

  end if

! Initialize iter and residual.

  res = 1.d0
  iter = 0

  do while ((res>ahepsilon_spectral).and.(iter<maxiter))

     iter = iter + 1

!    Save old Fourier coefficients.

     ac_old = ac
     cc_old = cc

!    Save previous values of rad and drad.

     rad_p = rad
     drad_p = drad


!    **************************
!    ***   INTERNAL STEPS   ***
!    **************************

!    For the surface update we use a two step
!    Runge-Kutta method.

     do q=1,2


!       ******************************************
!       ***   FIND FUNCTION TO BE INTEGRATED   ***
!       ******************************************

!       Integrand function.  This is a multiple of the
!       expansion H and is given by:
!
!       I  =  - ( Q  +  h'' )  ~  H
!
!       where Q is the (negative) of the source function "sdh"
!       calculated in the routine ah_source.
!
!       Notice that the routine ah_source not only calculates
!       "sdh" (-Q), but it also returns the area element and the
!       local contribution to the angular momentum integral.
!       The errorflag returns "false" is everything went OK,
!       and "true" if there was a problem (we got out of the
!       computational domain, or have a negative radius).

        do i=0,Ntheta
           call ah_source(rad(i),drad(i),theta(i),zcenter,sdh,sarea,sJmom,errorflag)
           if (errorflag) goto 100
           Ifunc(i) = sdh
           area(i) = sarea
           Jmom(i) = sJmom
        end do

!       The function "sdh" is singular at theta=(0,pi), so at those
!       points I need to correct.

        Ifunc(0) = Ifunc(1)
        area(0) = zero
        Jmom(0) = zero

        if (step==1) then
           Ifunc(Ntheta) = Ifunc(Ntheta-1)
           area(Ntheta) = zero
           Jmom(Ntheta) = zero
        end if


!       ***********************************************
!       ***   FIND FOURIER COEFFICIENTS FOR Ifunc   ***
!       ***********************************************

!       The Fourier coefficients are defined as:
!
!                      / pi/s
!       c   =  (2s/pi) |      F(theta) cos(j*theta) dtheta
!        j             / 0
!
!       with "s" a constant whose value depends on the symmetry:
!
!       s=1  if there is no equatorial symmetry.
!       s=2  if there is equatorial symmetry.
!
!       The above integrals are done using the trapezoidal rule,
!       which is second order accurate, but works extremely
!       well (it is in fact exact) for integrals of pure
!       cosines (and sines).

        do j=0,Nmode,step

!          First point.

           cc(j) = half*dtheta*Ifunc(0)

!          Interior points.

           do i=1,Ntheta-1
              cc(j) = cc(j) + dtheta*Ifunc(i)*cos(dble(j)*theta(i))
           end do

!          Last point.

           cc(j) = cc(j) + half*dtheta*Ifunc(Ntheta)*cos(dble(j)*theta(Ntheta))

!          Normalization.

           cc(j) = two*dble(step)*cc(j)/smallpi

!          We still need to add the contribution of minus h'' to Ifunc.
!          we do this by taking the second derivative of the expansion
!          in cosines.

           cc(j) = cc(j) + dble(j)**2*ac(j)

!          On first RK step save value of coefficients.

           if (q==1) cc_old(j) = cc(j)

        end do

!       Reconstruct Ifunc (for testing).

        !do i=0,Ntheta
        !   Ifunc(i) = half*cc(0)
        !   do j=step,Nmode,step
        !      Ifunc(i) = Ifunc(i) + (cc(j)-dble(j)**2*ac(j))*cos(dble(j)*theta(i))
        !   end do
        !end do


!       **************************************************************
!       ***   FIND NEW FOURIER COEFFICIENTS FOR HORIZON FUNCTION   ***
!       **************************************************************

!       Here I use essentially the method of Gundlach, but instead of
!       spherical harmonics (in 3D) we use a cosine expansion in theta,
!       so there is a slight modification.  The iteration takes the form:
!
!        n+1     n                   2
!       a    =  a   -  AA / (1 + BB j )  c
!        j       j                        j
!
!       where c_j are the Fourier coefficients of the function Ifunc defined
!       above.  The values of AA and BB given here were found empirically
!       and seem to work quite well.
!
!       The factor q/2 is there for the two step Runge-Kutta update.

        do j=0,Nmode,step
           ac(j) = ac_old(j) - half*AA*dble(q)/(one + BB*dble(j)**2)*cc(j)
        end do


!       ****************************************
!       ***   RECONSTRUCT HORIZON FUNCTION   ***
!       ****************************************

!       Reconstruct rad and its first derivative
!       using the expansion in cosines.

        do i=0,Ntheta
           rad(i)  = half*ac(0)
           drad(i) = zero
           do j=step,Nmode,step
              rad(i)  = rad(i) + ac(j)*cos(dble(j)*theta(i))
              drad(i) = drad(i) - dble(j)*ac(j)*sin(dble(j)*theta(i))
           end do
        end do

     end do

!    Save iterations data if needed (for testing).

     if ((rank==0).and.(ahsaveiter).and.(mod(iter,10)==0)) then

        write(1,"(A,i5)") '"Time   = ',iter
        write(1,"(A,ES20.10)") '"Center = ',zcenter

        do i=0,Ntheta
           write(1,*) theta(i),rad(i)
        end do

        write(1,*)

     end if


!    *************************
!    ***   FIND RESIDUAL   ***
!    *************************

!    Remember that this function Ifunc is proportional to
!    the expansion, so it should go to zero (and hence also
!    its Fourier coefficients) at a horizon.

     res = zero

     do j=0,Nmode,step
        res = max(res,abs(cc(j)))
     end do

!    Output to screen.

     !if ((rank==0).and.(mod(iter,10)==0)) print *, iter,res


!    ***********************************************
!    ***   DO WE NEED TO CHANGE THE STEP SIZE?   ***
!    ***********************************************

!    Find the difference between the two step Runge-Kutta
!    result and a single large Euler step.

     aux1 = zero
     aux2 = zero

     do j=0,Nmode,step
        tc(j) = ac_old(j) - AA/(one + BB*dble(j)**2)*cc_old(j)
        aux1 = aux1 + abs(ac(j)-tc(j))
        aux2 = aux2 + abs(ac_old(j))
     end do
 
     aux1 = aux1
     aux2 = aux2

!    If the difference between two small steps and one
!    big step is too large we ignore this iteration
!    and make AA smaller.

     if (aux1>0.05d0*aux2) then

         AA = half*AA
         !if (rank==0) print *, 'decrease AA: ',AA

         do j=0,Nmode,step
            ac(j) = ac_old(j)
         end do

         rad = rad_p
         drad = drad_p

!    If the difference is very small, then we can increase
!    the value of AA for the next iteration. But we don't
!    allow it to grow too much.

     else if (aux1<0.001d0*aux2) then

        AA = min(1.5d0*AA,ahAA)
        !if (rank==0) print *, 'increase AA: ',AA

     end if

  end do


! ************************************
! ***   HAVE WE FOUND A HORIZON?   ***
! ************************************

  100 continue

! If we got here before the maximum number of iterations,
! and we had no error, then we found a horizon.

  if ((iter<maxiter).and.(.not.errorflag)) then

     ahfound(k) = .true.
     ack(k,:) = ac(:)

     if ((rank==0).and.ahverbose) then
        write(*,'(A,i5,A)') ' Apparent horizon found after ',iter,' interations'
        print *, 'with Fourier coefficients (a,c):'
        do j=0,Nmode,step
           print *, j,ac(j),cc(j)
        end do
        print *
     end if

! If we reached the maximum number of iterations
! or we had an error then we did not find a horizon.

  else

     ahfound(k) = .false.

     if ((rank==0).and.ahverbose) then
        if (iter>=maxiter) then
           print *, 'No horizon found after  ',iter,' interations'
           print *, 'Residual = ',res
        else
           print *, 'No horizon found.'
        end if
        !do j=0,Nmode,step
        !   print *, j,ac(j)
        !end do
        print *
     end if

  end if

! Close file for iterations data (for testing).

  if ((rank==0).and.ahsaveiter) then
     close(1)
  end if


! ************************
! ***   OUTPUT FILES   ***
! ************************

! Determine file status and extension.

  if (firstcall(k)) then
     firstcall(k) = .false.
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

  if (N_ah==1) then
     write(filen,'(A)') '.tl'
  else
     write(filen,'(i1,A)') k,'.tl'
  end if


! *********************************
! ***   SAVE HORIZON POSITION   ***
! *********************************

! Save horizon position.

  if (rank==0) then

     if (filestatus=='replace') then
        open(1,file=trim(directory)//'/ah_radius'//trim(filen),form='formatted', &
        status=filestatus)
     else
        open(1,file=trim(directory)//'/ah_radius'//trim(filen),form='formatted', &
        status=filestatus,position='append')
     end if

!    Write current time, and horizon center.

     write(1,"(A,ES20.8)") '"Time   = ',t(0,0)
     write(1,"(A,ES20.8)") '"Center = ',zcenter

!    Save horizon position.

     if (ahfound(k)) then
        do i=0,Ntheta
           write(1,"(2ES20.8E3)") theta(i),rad(i)
        end do
     else
        do i=0,Ntheta
           write(1,"(2ES20.8E3)") theta(i),0.d0
        end do
     end if

!    Leave blank space before next time.

     write(1,*)

!    Close file.

     close(1)

  end if


! **************************************
! ***   SAVE HORIZON AREA AND MASS   ***
! **************************************

  if (rank==0) then

!    Save horizon area and horizon areal radius defined as:
!
!    R = sqrt(Area/pi) / 2

!    No horizon, area=0.

     if (.not.ahfound(k)) then

        ah_area = zero

!    Found a horizon, so we integrate the area
!    using the trapezoidal rule.

     else

        ah_area = half*dtheta*area(0)

        do i=1,Ntheta-1
           ah_area = ah_area + dtheta*area(i)
        end do

        ah_area = ah_area + half*dtheta*area(Ntheta)

        if (eqsym.and.(zcenter==0.d0)) ah_area = two*ah_area

     end if

     RAH = half*sqrt(ah_area/smallpi)

     if (ahverbose.and.ahfound(k)) then
        print *, 'Horizon area:',ah_area
        print *, 'Areal radius:',RAH
     end if

     if (filestatus=='replace') then
        open(1,file=trim(directory)//'/ah_area'//trim(filen),form='formatted', &
        status=filestatus)
     else
        open(1,file=trim(directory)//'/ah_area'//trim(filen),form='formatted', &
        status=filestatus,position='append')
     end if

     write(1,"(2ES20.8E3)") t(0,0),ah_area

     close(1)

!    Save horizon angular momentum.

     if (angmom) then

!       No horizon, so J=0.

        if (.not.ahfound(k)) then

           ah_J = zero

!       Found a horizon, so we integrate angular momentum
!       using the trapezoidal rule.

        else

           ah_J = half*dtheta*Jmom(0)

           do i=1,Ntheta-1
              ah_J = ah_J + dtheta*Jmom(i)
           end do

           ah_J = ah_J + half*dtheta*Jmom(Ntheta)

           if (eqsym.and.(zcenter==0.d0)) ah_J = two*ah_J

        end if

        if (ahverbose.and.ahfound(k)) then
           print *, 'Horizon angular momentum:',ah_J
        end if

        if (filestatus == 'replace') then
           open(1,file=trim(directory)//'/ah_J'//trim(filen),form='formatted', &
           status=filestatus)
        else
           open(1,file=trim(directory)//'/ah_J'//trim(filen),form='formatted', &
           status=filestatus,position='append')
        end if

        write(1,"(2ES20.8E3)") t(0,0),ah_J

        close(1)

     else

        ah_J = zero

     end if

!    Find and save horizon mass.
!
!    The irreducible mass is then defined as:
!
!    MI  =  R/2  =  sqrt(Area/pi) / 4
!
!    with R the horizon areal radius.
!
!    The total mass, including the angular momentum
!    contribution is (still missing the charge
!    contribution):
!
!                 4       2
!    M  =  sqrt( R  +  4 J ) / 2 R
!
!    with J the horizon angular momentum.

     if (ahfound(k)) then
        MI = half*RAH
        if (angmom) then
           MAH = half*sqrt(RAH**4 + 4.d0*ah_J**2)/RAH
        else
           MAH = MI
        end if
     else
        MI  = zero
        MAH = zero
     end if

     if (ahverbose.and.ahfound(k)) then
        print *, 'Irreducible horizon mass:',MI
        print *, 'Total horizon mass:      ',MAH
        print *
     end if

     if (filestatus == 'replace') then
        open(1,file=trim(directory)//'/ah_mass'//trim(filen),form='formatted', &
        status=filestatus)
     else
        open(1,file=trim(directory)//'/ah_mass'//trim(filen),form='formatted', &
        status=filestatus,position='append')
     end if

     write(1,"(2ES20.8E3)") t(0,0),MAH

     close(1)

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine ahspectral










  subroutine ah_source(rad,dh,th,zcenter,sdh,sarea,sJmom,errorflag)

! ***********************************************
! ***   SOURCES FOR APPARENT HORIZON FINDER   ***
! ***********************************************

! This subroutines calculates the sources for the apparent horizon ODE.
! For more details check eq. (6.7.21) from Alcubierre's book.
!
! We assume the horizon is a ray-body that can be expressed a r = h(theta).
! We then find a second order ODE for the function h(theta):
!
!              2  ij     i   j           k
! h''  =  - ( u  g   -  d F d F ) ( Gamma   d F  +  u K  )
!                                        ij  k         ij
!
!            rr  theta theta      r theta 2
!      /  [ g   g            -  (g       )  ]
!
!
! Here g is the spatial metric, F(r,theta) = r - h(theta), the Gammas
! are the Christoffel symbols, and u is the norm of dF/di:
!
!  2         i       rr       r theta         theta theta    2
! u  =  d F d F  =  g   -  2 g        h'  +  g           (h')
!        i
!
! But do notice that the whole expression has to be calculated in
! spherical coordinates (not in the cylindrical the code uses).
! The code does calculate the conformal spherical metric, but
! from the origin of the coordinate system, so they won't
! work if "zcenter" is not zero and we will need to recalculate
! them here.
!
! Notice also that the source is in fact singular at the axis
! (theta=0,pi), due to a term that goes as 1/rho. This term
! comes from:
!
!  pp      t
! g   Gamma   ~ cos(theta)/sin(theta)    (even in flat space)
!          pp
!
! (with p=phi,t=theta). This could be regularized, but I have
! not done so yet. So avoid calling it for those angles!
!
! An original version of this was coded by Erik Jimenez in the old
! code, but this version is very different.
!
! Note for parellel runs:  This routine is not really parallel,
! all processors in fact do the same calculation.  But when
! running in parallel, the interpolations need to be done on
! all processors at the same time, and the result is then
! reduced among all of them (processors that don't own
! the interpolating location return 0, and the results from
! all processors it added together).

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical errorflag          ! Variable for error control.
  logical interpflag         ! Interpolation flag.

  integer box,level,bb,ll    ! Box and level counters.
  integer i,j,k,m            ! Counters.

  real(8) rad,dh,th          ! Radius, its derivative, and theta.
  real(8) zcenter            ! z position of horizon center.
  real(8) sdh,sarea          ! Sources for horizon function and area.
  real(8) sJmom              ! Source for angular momentum.
  real(8) rc,zc              ! Coordinates from horizon center.
  real(8) r0,z0              ! Grid coordinates.
  real(8) interp             ! Interpolation function.
  real(8) u,u2               ! Norm of gradient of F.
  real(8) denom              ! Denominator in expression for h''.
  real(8) n_r,n_t            ! Components of unit normal vector.
  real(8) zero,half,one,two  ! Numbers.
  real(8) smallpi            ! Pi.
  real(8) epsilon            ! Small number.
  real(8) aux,diffr,diffz    ! Auxiliary.

! Interpolated conformal factor and derivatives.

  real(8) interp_phi,interp_psi,interp_psi4
  real(8) interp_Drphi,interp_Dzphi

! Interpolated metric components and derivatives.

  real(8) interp_A,interp_DrA,interp_DzA
  real(8) interp_B,interp_DrB,interp_DzB
  real(8) interp_C,interp_DrC,interp_DzC
  real(8) interp_H,interp_DrH,interp_DzH

  real(8) interp_C1,interp_DrC1,interp_DzC1
  real(8) interp_C2,interp_DrC2,interp_DzC2

! Interpolated extrinsic curvature.

  real(8) interp_trK

  real(8) interp_KA,interp_KB,interp_KC
  real(8) interp_KH,interp_KC1,interp_KC2

! Array for gradient of F in spherical coordinates.

  real(8) DF(1:3)

! Array for Kroenecker delta.

  real(8) delta(1:3,1:3)

! Array for derivatives of phi in spherica coordinates.

  real(8) Dphi(1:3)

! Arrays for spherical conformal metric and inverse.

  real(8) gd(1:3,1:3)
  real(8) gu(1:3,1:3)

! Array for spherical metric derivatives
! (first index indicates derivative).

  real(8) Dg(1:3,1:3,1:3)

! Array for spherical traceless extrinsic curvature.

  real(8) Kd(1:3,1:3)

! Arrays for spherical conformal Christoffel symbols
! (first index is up, the last two down).

  real(8) Gamma(1:3,1:3,1:3)


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  half = 0.5d0
  one  = 1.d0
  two  = 2.d0

  smallpi = acos(-1.d0)


! *****************************************
! ***   COORDINATES FOR INTERPOLATION   ***
! *****************************************

! Initialize errorflag.

  errorflag = .false.

! Check if the radius is negative.

  if (rad<0.d0) then
     !if (rank==0) print *, 'Negative radius!'
     sdh = 0.d0
     errorflag = .true.
     return
  end if

! Cylindrical coordinates:
! (rc,zc) from horizon center.
! (r0,z0) from origin.

  rc = rad*sin(th)
  zc = rad*cos(th)

  r0 = rc
  z0 = zc + zcenter

  if (eqsym.and.(z0<0.d0)) then
     z0 = - z0
  end if

! Check if the horizon radius is larger than the maximum allowed.

  if ((r0>ahrmax*rmaxl(0,0)).or.(abs(z0)>ahrmax*zmaxl(0,0))) then
     !if (rank==0) print *, 'Horizon radius is too large!',ahrmax,r0,z0
     sdh = 0.d0
     errorflag = .true.
     return
  end if

! Check if we are too close to the grid edges.

  if ((r0<rminl(0,0)+drl(0)).or.(r0>rmaxl(0,0)-drl(0)).or. &
      (z0<zminl(0,0)+dzl(0)).or.(z0>zmaxl(0,0)-drl(0))) then
     !if (rank==0) print *, 'Too close to grid edge!',r0,z0
     sdh = 0.d0
     errorflag = .true.
     return
  end if

! Check if horizons 1 and 2 have crossed the equator.
! For the moment, if we are using the "shooting" method,
! we will force them not to cross the equator, as otherwise
! the finder settles on the common horizon once it forms.

  if ((z0*zcenter<0.d0).and.(ahmethod=="shoot")) then
     !if (rank==0) print *, 'Horizon 1 or 2 crossed the equator: ',zcenter,zc
     sdh = 0.d0
     errorflag = .true.
     return
  end if


! *****************************************************
! ***   FIND GRID BOX AND LEVEL FOR INTERPOLATION   ***
! *****************************************************

! We need to interpolate at the highest level
! available at the current location.

  box = 0
  level = 0

  if (Nlmax>0) then
     do bb=0,Nb
        do ll=1,Nl(bb)
           if ((r0>rminl(bb,ll)+drl(ll)).and.(r0<rmaxl(bb,ll)-drl(ll)).and. &
               (z0>zminl(bb,ll)+dzl(ll)).and.(z0<zmaxl(bb,ll)-dzl(ll))) then
              box = bb
              level = ll
           end if
        end do
     end do
  end if


! ****************************
! ***   CONFORMAL FACTOR   ***
! ****************************

! Interpolate conformal factor and derivatives.

  interpvar => grid(box,level)%phi
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_phi,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  interpvar => grid(box,level)%psi
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_psi,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  interp_psi4 = interp_psi**4

  interpvar => grid(box,level)%Dr_phi
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_Drphi,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  interpvar => grid(box,level)%Dz_phi
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_Dzphi,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

! Transform derivatives to spherical coordinates.
! Remember that:
!
! d_R  =  (r d_r  +  z d_z ) / R
! d_T  =  (z d_r  -  r d_z )
!
! where (r,z) are the cylindrical coordinates
! and (R,T) the spherical coordinates.

  Dphi(1) = (rc*interp_Drphi + zc*interp_Dzphi)/rad
  Dphi(2) = (zc*interp_Drphi - rc*interp_Dzphi)
  Dphi(3) = zero


! **************************************
! ***   CONFORMAL SPHERICAL METRIC   ***
! **************************************

! Interpolate metric in cylindrical coordinates.

  interpvar => grid(box,level)%A
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_A,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  interpvar => grid(box,level)%Dr_A
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_DrA,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  interpvar => grid(box,level)%Dz_A
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_DzA,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  interpvar => grid(box,level)%B
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_B,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  interpvar => grid(box,level)%Dr_B
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_DrB,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  interpvar => grid(box,level)%Dz_B
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_DzB,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  interpvar => grid(box,level)%C
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_C,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  interpvar => grid(box,level)%Dr_C
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_DrC,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  interpvar => grid(box,level)%Dz_C
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_DzC,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  interpvar => grid(box,level)%H
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_H,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  interpvar => grid(box,level)%Dr_H
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_DrH,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  interpvar => grid(box,level)%Dz_H
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_DzH,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  if (angmom) then

     interpvar => grid(box,level)%C1
     aux = interp(box,level,r0,z0,interpflag)
     call MPI_ALLREDUCE(aux,interp_C1,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     interpvar => grid(box,level)%Dr_C1
     aux = interp(box,level,r0,z0,interpflag)
     call MPI_ALLREDUCE(aux,interp_DrC1,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     interpvar => grid(box,level)%Dz_C1
     aux = interp(box,level,r0,z0,interpflag)
     call MPI_ALLREDUCE(aux,interp_DzC1,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

     interpvar => grid(box,level)%C2
     aux = interp(box,level,r0,z0,interpflag)
     call MPI_ALLREDUCE(aux,interp_C2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     interpvar => grid(box,level)%Dr_C2
     aux = interp(box,level,r0,z0,interpflag)
     call MPI_ALLREDUCE(aux,interp_DrC2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     interpvar => grid(box,level)%Dz_C2
     aux = interp(box,level,r0,z0,interpflag)
     call MPI_ALLREDUCE(aux,interp_DzC2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  end if

! Transform to spherical coordinates.

  gd(1,1) = (rc**2*interp_A + zc**2*interp_B + two*rc**2*zc*interp_C)/rad**2
  gd(2,2) = (zc**2*interp_A + rc**2*interp_B - two*rc**2*zc*interp_C)
  gd(1,2) = rc*(zc*(interp_A - interp_B) + (zc**2 - rc**2)*interp_C)/rad
  gd(3,3) = rc**2*interp_H

  if (angmom) then
     gd(1,3) = rc**2*(rc**2*interp_C1 + zc*interp_C2)/rad
     gd(2,3) = rc**3*(zc*interp_C1 - interp_C2)
  else
     gd(1,3) = zero
     gd(2,3) = zero
  end if

  gd(2,1) = gd(1,2)
  gd(3,1) = gd(1,3)
  gd(3,2) = gd(2,3)

! Inverse spherical metric. When theta=0 we need to take care.
! Several expressions simplify since the metric becomes diagonal,
! but gd(3,3) goes to 0 and gu(3,3) becomes infinite. To
! avoid NaN's I set gu(3,3) to 0, but then we need to
! take care with expressions below that might use it.

  if ((th/=zero).and.(th/=smallpi)) then
     if (angmom) then
        aux = gd(1,1)*gd(2,2)*gd(3,3) + two*gd(1,2)*gd(1,3)*gd(2,3) &
            - gd(1,1)*gd(2,3)**2 - gd(2,2)*gd(1,3)**2 - gd(3,3)*gd(1,2)**2
        gu(1,1) = (gd(2,2)*gd(3,3) - gd(2,3)**2)/aux
        gu(2,2) = (gd(1,1)*gd(3,3) - gd(1,3)**2)/aux
        gu(3,3) = (gd(1,1)*gd(2,2) - gd(1,2)**2)/aux
        gu(1,2) = (gd(1,3)*gd(2,3) - gd(1,2)*gd(3,3))/aux
        gu(1,3) = (gd(1,2)*gd(2,3) - gd(1,3)*gd(2,2))/aux
        gu(2,3) = (gd(1,2)*gd(1,3) - gd(2,3)*gd(1,1))/aux
     else
        aux = gd(1,1)*gd(2,2) - gd(1,2)**2
        gu(1,1) = gd(2,2)/aux
        gu(2,2) = gd(1,1)/aux
        gu(3,3) = one/gd(3,3)
        gu(1,2) = - gd(1,2)/aux
        gu(1,3) = zero
        gu(2,3) = zero
     end if
  else
     gu(1,1) = one/gd(1,1)
     gu(2,2) = one/gd(2,2)
     gu(3,3) = zero ! This is really infinite, but I set it to zero to avoid NaNs.
     gu(1,2) = zero
     gu(1,3) = zero
     gu(2,3) = zero
  end if

  gu(2,1) = gu(1,2)
  gu(3,1) = gu(1,3)
  gu(3,2) = gu(2,3)


! ******************************************************
! ***   METRIC DERIVATIVES IN SPHERICAL COORDINATES   **
! ******************************************************

! Metric derivatives. Remember that:
!
! d_R  =  (r d_r  +  z d_z ) / R
! d_T  =  (z d_r  -  r d_z )
!
! where (r,z) are the cylindrical coordinates
! and (R,T) the spherical coordinates.

  diffr = (two*rc*interp_A + rc**2*interp_DrA + zc**2*interp_DrB &
        + two*rc*zc*(two*interp_C + rc*interp_DrC) - two*rc*gd(1,1))/rad**2
  diffz = (rc**2*interp_DzA + two*zc*interp_B + zc**2*interp_DzB &
        + two*rc**2*(interp_C + zc*interp_DzC) - two*zc*gd(1,1))/rad**2
  Dg(1,1,1) = (rc*diffr + zc*diffz)/rad
  Dg(2,1,1) = (zc*diffr - rc*diffz)

  diffr = (two*rc*interp_B + rc**2*interp_DrB + zc**2*interp_DrA &
        - two*rc*zc*(two*interp_C + rc*interp_DrC))
  diffz = (rc**2*interp_DzB + two*zc*interp_A + zc**2*interp_DzA &
        - two*rc**2*(interp_C + zc*interp_DzC))
  Dg(1,2,2) = (rc*diffr + zc*diffz)/rad
  Dg(2,2,2) = (zc*diffr - rc*diffz)

  diffr = ((one - (rc/rad)**2)*(zc*(interp_A-interp_B) + (zc**2-rc**2)*interp_C) &
        + rc*(zc*(interp_DrA-interp_DrB) + (zc**2-rc**2)*interp_DrC - two*rc*interp_C))/rad
  diffz = (rc/rad)*(-(zc/rad**2)*(zc*(interp_A-interp_B) + (zc**2-rc**2)*interp_C) &
        + (interp_A-interp_B) + zc*(interp_DzA-interp_DzB) + (zc**2-rc**2)*interp_DzC + two*zc*interp_C)
  Dg(1,1,2) = (rc*diffr + zc*diffz)/rad
  Dg(2,1,2) = (zc*diffr - rc*diffz)

  diffr = rc*(two*interp_H + rc*interp_DrH)
  diffz = rc**2*interp_DzH
  Dg(1,3,3) = (rc*diffr + zc*diffz)/rad
  Dg(2,3,3) = (zc*diffr - rc*diffz)

  if (angmom) then

     diffr = rc/rad*((two - (rc/rad)**2)*(rc**2*interp_C1 + zc*interp_C2) &
           + rc*(two*rc*interp_C1 + rc**2*interp_DrC1 + zc*interp_DrC2))
     diffz = rc**2/rad*(rc**2*interp_DzC1 + interp_C2 + zc*interp_DzC2 &
           - zc/rad**2*(rc**2*interp_C1 + zc*interp_C2))
     Dg(1,1,3) = (rc*diffr + zc*diffz)/rad
     Dg(2,1,3) = (zc*diffr - rc*diffz)

     diffr = rc**2*(3.d0*(zc*interp_C1 - interp_C2) + rc*(zc*interp_DrC1 - interp_DrC2))
     diffz = rc**3*(interp_C1 + zc*interp_DzC1 - interp_DzC2)
     Dg(1,2,3) = (rc*diffr + zc*diffz)/rad
     Dg(2,2,3) = (zc*diffr - rc*diffz)

  else

     Dg(1,1,3) = zero
     Dg(2,1,3) = zero

     Dg(1,2,3) = zero
     Dg(2,2,3) = zero

  end if

  Dg(1,2,1) = Dg(1,1,2)
  Dg(2,2,1) = Dg(2,1,2)

  Dg(1,3,1) = Dg(1,1,3)
  Dg(2,3,1) = Dg(2,1,3)

  Dg(1,3,2) = Dg(1,2,3)
  Dg(2,3,2) = Dg(2,2,3)

! Remember that nothing depends on phi.

  Dg(3,:,:) = zero


! *******************************
! ***   EXTRINSIC CURVATURE   ***
! *******************************

! Interpolate trK.

  interpvar => grid(box,level)%trK
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_trK,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

! Interpolate physical extrinsic curvature.

  interpvar => grid(box,level)%KA
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_KA,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  interpvar => grid(box,level)%KB
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_KB,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  interpvar => grid(box,level)%KC
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_KC,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  interpvar => grid(box,level)%KH
  aux = interp(box,level,r0,z0,interpflag)
  call MPI_ALLREDUCE(aux,interp_KH,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  if (angmom) then

     interpvar => grid(box,level)%KC1
     aux = interp(box,level,r0,z0,interpflag)
     call MPI_ALLREDUCE(aux,interp_KC1,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

     interpvar => grid(box,level)%KC2
     aux = interp(box,level,r0,z0,interpflag)
     call MPI_ALLREDUCE(aux,interp_KC2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  end if

! Transform to spherical coordinates.

  Kd(1,1) = (rc**2*interp_KA + zc**2*interp_KB + two*rc**2*zc*interp_KC)/rad**2
  Kd(2,2) = (zc**2*interp_KA + rc**2*interp_KB - two*rc**2*zc*interp_KC)
  Kd(1,2) = rc*(zc*(interp_KA - interp_KB) + (zc**2 - rc**2)*interp_KC)/rad
  Kd(3,3) = rc**2*interp_KH

  if (angmom) then
     Kd(1,3) = rc**2*(rc**2*interp_KC1 + zc*interp_KC2)/rad
     Kd(2,3) = rc**3*(zc*interp_KC1 - interp_KC2)
  else
     Kd(1,3) = zero
     Kd(2,3) = zero
  end if

  Kd(2,1) = Kd(1,2)
  Kd(3,1) = Kd(1,3)
  kd(3,2) = Kd(2,3)


! ************************************************************
! ***   FIND GRADIENT OF F, ITS NORM AND THE DENOMINATOR   ***
! ************************************************************

! Gradient of F (index up).

  DF(1) = (gu(1,1) - gu(2,1)*dh)/interp_psi4
  DF(2) = (gu(1,2) - gu(2,2)*dh)/interp_psi4
  DF(3) = (gu(1,3) - gu(2,3)*dh)/interp_psi4

! Norm of gradient of F.

  u2 = (gu(1,1) - two*gu(1,2)*dh + gu(2,2)*dh**2)/interp_psi4
  u  = dsqrt(abs(u2))

! Denominator: g^rr*g^tt - (g^rt)^2.

  denom = (gu(1,1)*gu(2,2) - gu(1,2)**2)/interp_psi4**2


! *******************************
! ***   CHRISTOFFEL SYMBOLS   ***
! *******************************

! Kroenecker delta.

  do i=1,3
     do j=1,3
        delta(i,j) = int(min(i,j)/max(i,j))
     end do
  end do

! Initialize Christoffels to zero.

  Gamma = zero

! Conformally flat Christoffels (for testing).
! I only consider the non-zero contribution.
! I do not consider the case when the first (up)
! index is equal to 3 since this is never used.
! This should be commented out for production runs.

  !Gamma(1,2,2) = - rad
  !Gamma(1,3,3) = - rad*sin(th)**2
  !Gamma(2,1,2) = one/rad
  !Gamma(2,3,3) = - sin(th)*cos(th)
  !Gamma(2,2,1) = Gamma(2,1,2)

! Terms from conformal Christoffel symbols.
! I do not consider the case when the first (up)
! index is equal to 3 since it is never used.

  do k=1,2
     do i=1,3
        do j=1,3
           do m=1,3
              Gamma(k,i,j) = Gamma(k,i,j) &
                           + half*gu(k,m)*(Dg(i,j,m) + Dg(j,i,m) - Dg(m,i,j))
           end do
        end do
     end do
  end do

! Terms from derivatives of conformal factor.
! I do not consider the case when the first (up)
! index is equal to 3 since it is never used.

  do k=1,2
     do i=1,3
        do j=1,3
           Gamma(k,i,j) = Gamma(k,i,j) + two*(delta(k,i)*Dphi(j) + delta(k,j)*Dphi(i) &
                        - gd(i,j)*(gu(k,1)*Dphi(1) + gu(k,2)*Dphi(2)))
        end do
     end do
  end do


! ***************************************
! ***   SOURCE FOR HORIZON EQUATION   ***
! ***************************************

! Now construct the source term by term.

! Initialize source to 0.

  sdh = zero

! 1) Term:  u^2 g^ij Gamma^k_ij (d_k F).

  do i=1,3
     do j=1,3
        sdh = sdh + u2*gu(i,j)*(Gamma(1,i,j) - Gamma(2,i,j)*dh)/interp_psi4
     end do
  end do

! 2) Term:  (d^i F) (d^j F) Gamma^k_ij (d_k F).

  do i=1,3
     do j=1,3
        sdh = sdh - DF(i)*DF(j)*(Gamma(1,i,j) - Gamma(2,i,j)*dh)
     end do
  end do

! 3) Term:  u^3 trK.

  sdh = sdh + u**3*interp_trK

! 4) Term:  u (d^i F) (d^j F) K_ij.

  do i=1,3
     do j=1,3
        sdh = sdh - u*DF(i)*DF(j)*Kd(i,j)
     end do
  end do

! 5) Finally change sign and divide by denom.

  sdh = - sdh/denom

! If theta=0 (or theta=pi) we need to correct.  At the moment
! I just set sdh=0. This can be fixed by regularizing the source
! at theta=0 (or theta=pi) but I don't think this is really needed.
!
! Notice that the source is only 0 at theta={0,pi} if dh=0.
! Otherwise it is actually singular.  This means that for the
! shhoting method we can set it to 0 when we start integrating
! from the axis, but when we return to the axis on the other
! side it actually blows up unless we do have a horizon.

  epsilon = 1.d-5

  if ((abs(th)<epsilon).or.(abs(th-smallpi)<epsilon)) then
     sdh = zero
  end if

! Check for NaN's.

  if ((sdh/=sdh).or.(abs(sdh)>1.d50)) then
     !print *, 'Source diverged',sdh
     sdh = 1.d10
     errorflag = .true.
     return
  end if


! *****************************
! ***   FIND AREA ELEMENT   ***
! *****************************

! In order to calculate the area we start from
! the line element:
!
! dl^2  =  grr dr^2 + gtt dtheta^2 + gpp dphi^2
!
!       + 2 ( grt dr dtheta + grp dr dphi + gtp dtheta dphi )
!
! We then take r = h(theta), so that dr = h'*dtheta.
! Substituting we find:
!
! dl^2  =  [ grr (h')^2 + 2 grt h' + gtt ] dtheta^2
!
!       + 2 [ grp h' + gtp ] dtheta dphi + gpp dphi^2
!
! This is a 2D metric with area element equal to
! the root of the determinant:
!
! g  =  gpp [ grr (h')^2 + 2 grt h' + gtt ]
!
!    -  [ grp h' + gtp ]^2 
!
! The integral over phi is trivially equal to 2*pi.
! But do remember that we are using the conformal
! spherical metric, so we need to multiple with psi^4.

  aux = two*smallpi*interp_psi4

  sarea = 0.d0

  if (angmom) then
     sarea = aux*sqrt(abs(gd(3,3)*(dh**2*gd(1,1) + two*dh*gd(1,2) + gd(2,2)) &
           - (dh*gd(1,3) + gd(2,3))**2))
  else
     sarea = aux*sqrt(abs(gd(3,3)*(dh**2*gd(1,1) + two*dh*gd(1,2) + gd(2,2))))
  end if


! **********************************************
! ***   FIND ANGULAR MOMENTUM CONTRIBUTION   ***
! **********************************************

! Find unit normal vector to surface r=h(theta).
! The vector is found by asking for it for:
!
! n.n = 1
! n.v = 0
!
! where v is the tangent vector to the surface r=r(theta):
!
! v = (dh,1)
!
! Notice that for Kerr we just have n^r=1/sqrt(grr) and n^t=0.

  aux = - (gd(1,1)*dh + gd(1,2))/(gd(1,2)*dh + gd(2,2))

  n_r = one/sqrt(interp_psi4*(gd(1,1) + gd(2,2)*aux**2 + two*gd(1,2)*aux))
  n_t = n_r*aux

  !print *, interp_psi4*(gd(1,1)*n_r**2 + gd(2,2)*n_t**2 + two*gd(1,2)*n_r*n_t)
  !print *, interp_psi4*(gd(1,1)*n_r*dh + gd(2,2)*n_t + gd(1,2)*(n_r + n_t*dh))
  !print *

! The angular momentum contribution is given by:
!
!           r         t
! dJ  =  ( n  K   +  n  K  ) dA / 8pi
!              rp        tp

  if (angmom) then
     sJmom = (n_r*Kd(1,3) + n_t*Kd(2,3))*sarea/(8.d0*smallpi)
  else
     sJmom = zero
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine ah_source

