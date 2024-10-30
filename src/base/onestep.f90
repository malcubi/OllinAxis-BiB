!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/onestep.f90,v 1.64 2021/03/27 20:55:10 malcubi Exp $

  recursive subroutine onestep(level)

! *********************************
! ***   ADVANCE ONE TIME STEP   ***
! *********************************

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical contains

  integer i            ! Counter.
  integer box,level    ! Box number and level counters.
  integer k            ! Auxiliary counter.
  integer iter         ! Counter for internal iterations.
  integer niter        ! Number of internal iterations.
  integer bmax         ! Number of boxes at this level.
  integer bbox         ! Auxiliary for restriction.

  real(8) dtw          ! Internal time step.
  real(8) weight       ! Weight for rk4.


! ******************************************
! ***   NUMBER OF BOXES FOR THIS LEVEL   ***
! ******************************************

! Iterate over refinement boxes at the current
! grid level l. But notice that level 0 only
! exists on box 0.

  if (level==0) then
     bmax = 0
  else
     bmax = Nb
  end if


! *****************************************
! ***   ITERATE OVER REFINEMENT BOXES   ***
! *****************************************

! Iterate over boxes.

  do box=0,bmax

!    Check that the current box does indeed have
!    this level. If it doesn't, go to next cycle.

     if (Nl(box)<level) cycle


!    *********************************
!    ***   POINT TO CURRENT GRID   ***
!    *********************************

     call currentgrid(box,level,grid(box,level))


!    ******************************
!    ***   SAVE OLD TIME STEP   ***
!    ******************************

     call saveold


!    *******************************
!    ***   AUXILIARY VARIABLES   ***
!    *******************************

!    Before calculating sources, make sure we
!    have the correct auxiliary quantities.
!    This step mighti n fact not by needed, and
!    it makes the code slower, but I need to check.

     call auxiliary


!    **********************************************
!    ***   FIND NUMBER OF INTERNAL ITERATIONS   ***
!    **********************************************

     if (integrator=="icn") then
        niter = icniter
     else if (integrator=="rk4") then
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

        if (integrator=="icn") then

!          In ICN all iterations except the last one
!          jump only half a time step.

           if (iter<niter) then
              dtw = 0.5d0*dt
           else
              dtw = dt
           end if

!       Fourth order Runge-Kutta.

        else if (integrator=="rk4") then

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


!       *********************************
!       ***   SOURCES FOR SPACETIME   ***
!       *********************************

        if (spacetime=="dynamic") then

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


!       ******************************
!       ***   SOURCES FOR MATTER   ***
!       ******************************

        if (mattertype/="vacuum") then
           call sources_matter
        end if


!       *******************************
!       ***   BOUNDARY CONDITIONS   ***
!       *******************************

!       Outer boundaries are only needed for level 0.
!       The boundary conditions are applied on all processors,
!       the synchronization step below should fix this.

        if (level==0) then

!          Static and flat boundaries.  The routine "simpleboundary"
!          applies either static or flat boundaries to the sources
!          of all those arrays that are not declared with the
!          attribute NOBOUND.

           if ((boundtype=="static").or.(boundtype=="flat")) then

              call simpleboundary

!          Radiative boundaries for geometry. The routine "radiative_geometry"
!          applies outgoing wave boundary conditions to the geometric variables.
!          These boundary conditions take into account the characteristic structure,
!          but not the constraints. Notice that matter variables should define
!          their own radiative boundaries in their source routines.

           else if (boundtype=="radiative") then

              call radiative_geometry

!          Constraint preserving boundary conditions. STILL NOT IMPLEMENTED.

           else if (boundtype=="constraint") then

              call radiative_geometry
              call constraintbound

           end if

        end if


!       *****************************************************
!       ***   FOR RUNGE-KUTTA ADD TO ACCUMULATOR ARRAYS   ***
!       *****************************************************

!       The accumulator arrays add the contributions from the
!       different Runge-Kutta iterations with their corresponding
!       weights.  The routine "accumulate" (generated by perl)
!       works in the following way:
!
!       1) On the first step (iter=1) it just copies the source array times
!          its weight to the accumulator array.
!       2) On all subsequent steps except the last (1<iter<niter), the routine
!          adds the source times its weight to the accumulator array.
!       3) On the last call (iter=niter), it adds the final source times its
!          weight, but stores the result back into the source arrays so
!          that the last update will work correctly.

        if (integrator=="rk4") then
           call accumulate(iter,niter,weight)
        end if


!       ****************************
!       ***   UPDATE VARIABLES   ***
!       ****************************

        call update(dtw)


!       *************************************************
!       ***   FOR FINE GRIDS INTERPOLATE BOUNDARIES   ***
!       *************************************************

!       For fine grids we need to interpolate from the new
!       time level of the coarse grid to get boundary data.
!
!       Remember that the coarse grid has already advanced
!       to the next time level.

        if (level>0) then
           call finebound(box,level,dtw,.true.)
        end if


!       **********************
!       ***   SYMMETRIES   ***
!       **********************

!       Symmetries on axis and equator.

        if (eqsym.and.ownequator) then
           call symmetries_z
        end if

        if (ownaxis) then
           call symmetries_r
        end if


!       *****************************************
!       ***   SYNCHRONIZE ACROSS PROCESSORS   ***
!       *****************************************

!       If we have more than one processor we must now
!       synchronize ghost zones.

        if (size>1) then
           call syncall
        end if


!       *******************************
!       ***   AUXILIARY VARIABLES   ***
!       *******************************

!       Auxiliary variables for geometry and matter.

        call auxiliary


!       ***********************************
!       ***   END INTERNAL ITERATIONS   ***
!       ***********************************

     end do


!    ****************************************************
!    ***   ADVANCE LOCAL TIME AND TIME STEP COUNTER   ***
!    ****************************************************

!    Save old local times.

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
! current subroutine "onestep" recursively.

  if (level<Nlmax) then
     call onestep(level+1)
     call onestep(level+1)
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
! the same time level. But notice that here we also
! sync all derivatives.

  if (level>0) then

!    Loop over boxes at this level.

     do box=0,Nb
        if (Nl(box)<level) cycle
        do k=0,Nb
           if ((Nl(k)<level).or.(k==box)) cycle
           call syncboxes(box,k,level,.true.)
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
!
! Remember to fix the symmetries for level-1
! and synchronize again for multi-processor runs!

  if (level>0) then

     do box=0,bmax

        if (Nl(box)<level) cycle

        if (mod(s(box,level),2)==0) then

!          Restrict data on coarse level.

           call restrict(box,level,.true.)

!          Figure out on which box is level-1.

           if (level==1) then
              bbox = 0
           else
              bbox = box
           end if

!          Point to grid on level-1.

           call currentgrid(bbox,level-1,grid(bbox,level-1))

           if (ownaxis) then
              call symmetries_r
           end if

           if (eqsym.and.ownequator) then
              call symmetries_z
           end if

           if (size>1) then
              call syncall
           end if

!          Auxiliary quantities for matter and geometry on level-1.

           call auxiliary

        end if

     end do

  end if


! **************************
! ***   SPECIAL OUTPUT   ***
! **************************

! Here we do output of ALL time steps for the current box and level,
! and not only at the coarse level. Notice that the output times
! will then not coincide, as we will have more output for fine
! levels than for coarse levels.
!
! This is really intended only for testing, and is done when the
! corresponding parameter Noutput is equal to 0.

  if (Noutput0D*Noutput1D*Noutput2D==0) then

     do box=0,bmax

        if (Nl(box)<level) cycle

        call currentgrid(box,level,grid(box,level))

        if (Noutput0D==0) then
           do i=1,nvars0D
              call grabarray(trim(outvars0Darray(i)))
              call save0Dvariable(trim(outvars0Darray(i)),directory,box,level,s(box,level),t(box,level),'old')
           end do
        end if

        if (Noutput1D==0) then
           do i=1,nvars1D
              call grabarray(trim(outvars1Darray(i)))
              call save1Dvariable(trim(outvars1Darray(i)),directory,box,level,outparallel,'old')
           end do
        end if

        if (Noutput2D==0) then
           do i=1,nvars2D
              call grabarray(trim(outvars2Darray(i)))
              call save2Dvariable(trim(outvars2Darray(i)),directory,box,level,outparallel,'old')
           end do
        end if

     end do

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine onestep
