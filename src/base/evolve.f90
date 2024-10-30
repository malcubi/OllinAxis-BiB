!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/evolve.f90,v 1.59 2022/10/02 04:36:54 erik Exp $

  subroutine evolve

! **********************************
! ***   MAIN ITERATION ROUTINE   ***
! **********************************

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical contains

  integer i,j,k,step      ! Grid point and time step counters.
  integer box,level,proc  ! Box, level and processor counters.

  real(8) zero,half       ! Numbers.
  real(8) aux1,aux2       ! Auxiliary.


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  half  = 0.5d0


! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

! Loop over boxes and grid levels.

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Find r coordinate.  We always stagger the axis and
!       take a number of ghost zones on the negative side.

        if (rbox(box)==0.d0) then
           do i=1-ghost,Nr
              r(i,:) = (dble(Nminl_r(box,rank) + i) - 0.5d0)*dr
           end do
        else
           do i=1-ghost,Nr
              r(i,:) = rbox(box) + (dble(Nminl_r(box,rank) + i) - 0.5d0*dble(Nrbox(box) - ghost + 1))*dr
           end do
        end if

!       Find z coordinate.

        if (eqsym.and.(zbox(box)==0.d0)) then
           do j=1-ghost,Nz
              z(:,j) = (dble(Nminl_z(box,rank) + j) - 0.5d0)*dz
           end do
        else
           do j=1-ghost,Nz
              z(:,j) = zbox(box) + (dble(Nminl_z(box,rank) + j) - 0.5d0*dble(Nzbox(box) - ghost + 1))*dz
           end do
        end if

!       Distance to origin.

        rr = sqrt(r**2 + z**2)

!       Minimum and maximum values of (r,z) across processors.

        aux1 = r(1-ghost,0)
        call MPI_Allreduce(aux1,aux2,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
        rminl(box,level) = aux2

        aux1 = r(Nr,0)
        call MPI_Allreduce(aux1,aux2,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
        rmaxl(box,level) = aux2

        aux1 = z(0,1-ghost)
        call MPI_Allreduce(aux1,aux2,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
        zminl(box,level) = aux2

        aux1 = z(0,Nz)
        call MPI_Allreduce(aux1,aux2,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
        zmaxl(box,level) = aux2

    end do
  end do


! *****************************************************
! ***   SET TIME AND NUMBER OF TIME STEPS TO ZERO   ***
! *****************************************************

! Step counter and current time.

  s = 0
  t = 0.d0

! Old times.

  t1 = 0.d0
  t2 = 0.d0


! *****************************
! ***   FIND INITIAL DATA   ***
! *****************************

! Call initial data routine.  We do not loop over
! boxes and grid levels here since the individual
! initial data routines need to be aware of the
! whole grid structure.

  call initial

! Make sure we have correct symmetries at axis and equator.

  do box=0,Nb
     do level=min(1,box),Nl(box)

        call currentgrid(box,level,grid(box,level))

        if (ownaxis) then
           call symmetries_r
        end if

        if (eqsym.and.ownequator) then
           call symmetries_z
        end if

     end do
  end do

! If we have more than one processor synchronize ghost zones.

  if (size>1) then
     do box=0,Nb
        do level=min(1,box),Nl(box)
           call currentgrid(box,level,grid(box,level))
           call syncall
        end do
     end do
  end if


! ********************************
! ***   AUXILIARY QUANTITIES   ***
! ********************************

! Calculate auxiliary quantities for matter and geometry.

  do box=0,Nb
     do level=min(1,box),Nl(box)
        call currentgrid(box,level,grid(box,level))
        call auxiliary
     end do
  end do


! *****************************************
! ***   INITIALIZE GAMMA DRIVER SHIFT   ***
! *****************************************

  if ((shift=="Gammadriver2").or.(shift=="Gammadrivershock2")) then

     if (.not.driverD0) then

        do box=0,Nb
           do level=min(1,box),Nl(box)

              call currentgrid(box,level,grid(box,level))

              dtbeta_r = drivercsi*Delta_r
              dtbeta_z = drivercsi*Delta_z

              if (angmom) then
                 dtbeta_p = drivercsi*Delta_p
              end if

           end do
        end do

     end if

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
! the same time level. But notice that here we sync
! also all derivatives.

  do level=1,Nlmax
     do box=0,Nb
        if (Nl(box)<level) cycle
        do k=0,Nb
           if ((Nl(k)<level).or.(k==box)) cycle
           call syncboxes(box,k,level,.true.)
        end do
     end do
  end do


! ********************************************
! ***   MAXIMAL SLICING FOR INITIAL DATA   ***
! ********************************************

! If we want an initial maximal solve, do it here.

  if ((slicing=='maximal').or.(ilapse=='maximal').or.(maximalevery/=0)) then
     call maximalslicing
  end if


! *****************************
! ***   ANALYSIS ROUTINES   ***
! *****************************

! Call the different analysis routines.

  do box=0,Nb
     do level=min(1,box),Nl(box)
        call currentgrid(box,level,grid(box,level))
        call analysis(box,level)
     end do
  end do


! ********************
! ***   HORIZONS   ***
! ********************

  if (ahfind.and.(t(0,0)>=ahafter)) then
     call horizon_finder
  end if


! ****************************************
! ***   FIND SOURCES AT INITIAL TIME   ***
! ****************************************

! Here we calculate the sources for all evolution
! equations once.  This is not really necessary,
! it is only done so that we can output them at
! the initial time if needed.

  do box=0,Nb
     do level=min(1,box),Nl(box)

        if (spacetime=="dynamic") then

!          Point to current grid.

           call currentgrid(box,level,grid(box,level))

!          Sources for geometry.

           call sources_geometry

!          Sources for lapse. The sources for the lapse should
!          be calculated after the sources for the geometry,
!          since some lapse conditions might require strK.

           call sources_lapse

!          Sources for shift. The sources for the shift should
!          be calculated after the sources for the geometry
!          and the sources for the lapse, since the shift
!          conditions use sDeltar and some use salpha.

           if (shift/="none") then
              call sources_shift
           end if

        end if

!       Sources for matter.

        if (mattertype/="vacuum") then
           call sources_matter
        end if

     end do
  end do


! *********************************
! ***   SAVE THE INITIAL DATA   ***
! *********************************

! Save 0D, 1D and 2D data.

  call save0Ddata
  call save1Ddata
  call save2Ddata

  if (convert_to_3D) then

     call save3D ('psi',directory)

  end if
! Checkpoint initial data, except when we are
! restarting from a checkpoint file (since we
! already have it).

  if (checkpointinitial.and.(idata/="checkpoint")) then
     call checkpointsave
  end if


! ****************************
! ***   OUTPUT TO SCREEN   ***
! ****************************

  if (rank==0) then
     print *,'------------------------------'
     print *,'|  Time step  |     Time     |'
     print *,'------------------------------'
     write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',0,'   | ',t(0,0),'  | '
  end if


! *************************************
! ***   START MAIN EVOLUTION LOOP   ***
! *************************************

  do step=1,Nt


!    ***************************************
!    ***   DO WE ADJUST THE TIME STEP?   ***
!    ***************************************

!    If needed, adjust the time step according to
!    the values of the propagation speeds on the
!    base grid.

     if (adjuststep) then
        call cfl(step)
     end if


!    ****************************************
!    ***   ADVANCE ONE GLOBAL TIME STEP   ***
!    ****************************************

!    We only call the time stepping routine for
!    the coarse grid.  It is then called recursively
!    from inside for the other boxes and levels.

     call onestep(0)


!    ***************************
!    ***   MAXIMAL SLICING   ***
!    ***************************

!    If we want an initial maximal solve, do it here.

     if ((slicing=='maximal').or. &
        ((maximalevery/=0).and.(mod(step,maximalevery)==0))) then
        call maximalslicing
     end if


!    ****************************
!    ***   ANALYSIS ROUTINES  ***
!    ****************************

!    Call the different analysis routines. They are
!    called only when we need output.
!
!    Notice that if any of the Noutput parameters
!    is zero we don't do any output here (this is
!    a special case for testing purposes only).

     if (Noutput0D*Noutput1D*Noutput2D>0) then
        if ((mod(step,Noutput0D)==0).or. &
            (mod(step,Noutput1D)==0).or. &
            (mod(step,Noutput2D)==0)) then
           do box=0,Nb
              do level=min(1,box),Nl(box)
                 call currentgrid(box,level,grid(box,level))
                 call analysis(box,level)
              end do
           end do
        end if
     end if


!    ********************
!    ***   HORIZONS   ***
!    ********************

!    We only look for apparent horizons every "ahfind_every"
!    time steps.  Notice that even if ahfind_every=1, we
!    only look for horizons at the coarse time steps.

     if (ahfind.and.(t(0,0)>=ahafter)) then
        if (mod(step,ahfind_every).eq.0) then
           call horizon_finder
        end if
     end if


!    ***************************
!    ***   WAVE EXTRACTION   ***
!    ***************************

!    Call the wave extraction subroutine only when we need output.

     if (wave_extract.and.(mod(step,wavextract_every).eq.0)) then
        call wavextract
     end if


!    *****************************
!    ***   SAVE DATA TO FILE   ***
!    *****************************

!    Save 0D data (only every Noutput0D time steps).

     if (Noutput0D>0) then
        if (mod(step,Noutput0D)==0) then
           call save0Ddata
        end if
     end if

!    Save 1D data (only every Noutput1D time steps).

     if (Noutput1D>0) then
        if (mod(step,Noutput1D)==0) then
           call save1Ddata
        end if
     end if

!    Save 2D data (only every Noutput2D time steps).

     if (Noutput2D>0) then
        if (mod(step,Noutput2D)==0) then
           call save2Ddata
        end if
     end if

!    Save checkpoint data.

     if (checkpoint.and.(mod(step,Ncheckpoint)==0)) then
        call checkpointsave
     end if


!    ***********************************
!    ***   END MAIN EVOLUTION LOOP   ***
!    ***********************************

!    Time step information to screen.

     if (rank==0) then
        if (mod(step,Ninfo).eq.0) then
           write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',step,'   | ',t(0,0),'  | '
        end if
     end if

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine evolve

