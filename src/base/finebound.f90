!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/finebound.f90,v 1.28 2021/02/18 15:44:40 malcubi Exp $

  subroutine finebound(box,level,dtw,all)

! ****************************************
! ****   FINE BOUNDARY INTERPOLATION   ***
! ****************************************

! For fine grids we need to interpolate from the new
! time level of the coarse grid to get boundary data.
!
! The logical flag "all", when true, means the boundary
! conditions are applied to all evolving arrays.  When
! false the condition is applied only to a single array.
!
! When we apply the condition to a single array, the
! corresponding coarse array, fine array, and fine array
! in the previous time level, and boundaries on previos
! time levels must be assigned to the pointers "coarsevar",
! "finevar", "finevar_p" and "finebound_bound_*" BEFORE
! calling the routine.
!
! At the moment, the interpolation in time is a simple
! Lagrange quadratic interpolation, except for the
! first two timesteps for which we don't have enough
! points to the past, in which case it is only linear.
!
! I have tried cubic Lagrange interpolation (it is commented
! out in a couple of places below), but it seems to cause
! an instability!  If we want to go to cubic interpolation
! we will probably need a monotone cubic spline.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical contains
  logical all              ! Do we apply the boundary to all evolving arrays?
  logical flag1,flag2      ! Interpolator flags.

  integer i,j              ! Grid point counters.
  integer imax,jmax        ! Size of boundary region.
  integer box,level        ! Box number and level counters.
  integer m,n              ! Auxiliary counters.
  integer bbox             ! Auxiliary box identificator.
  integer border           ! Order of interpolation at fine boundary.

  real(8) dtw              ! Internal time step.
  real(8) tp,t0,tm1,tm2,tl ! Local time at different time steps (for interpolation).
  real(8) r0,z0,interp     ! For interpolation.
  real(8) aux1,aux2

  character(2) bd          ! Identifies the boundary in which we are: (rL,rR,zL,zR)


! ************************************************************
! ***   FIGURE OUR FROM WHICH COARSE GRID TO INTERPOLATE   ***
! ************************************************************

! Normally we interpolate the boundaries from the coarse
! level on the same box we are in.  But if we are at
! level 1 we always intepolate from box 0 (the base grid).

  if (level==1) then
     bbox = 0
  else
     bbox = box
  end if


! ************************************
! ***   TIMES FOR INTEROPOLATION   ***
! ************************************

! Figure out times for interpolation:
!
!    tp  = Local time on coarse grid (which has already advanced).
!    t0  = Local time on fine grid.
!    tm1 = Time 1 step  to the past for fine grid.
!    tm2 = Time 2 steps to the past for fine grid.
!    tl  = t0 + dtw (time to which we want to interpolate).
!
! Notice that I am not assuming here that the time steps are uniform.

  tp  = t(bbox,level-1)
  t0  = t(box ,level  )

  tm1 = t1(box,level)
  tm2 = t2(box,level)

  tl = t0 + dtw


! ***************************************************
! ***   INTERPOLATE BOUNDARIES FROM COARSE GRID   ***
! ***************************************************

! How many boundary points get data from interpolation?

  imax = ghost-1
  jmax = ghost-1

! Order of interpolation at boundary.  If we are at the
! first or second time steps the order is 1, otherwise
! the order is 2.

  if (tm1==0.d0) then
     border = 1
  else
     border = 2
  end if

! Apply boundary condition. Since here I don't know how many
! variables are evolving, I generate the necessary code
! automatically at compile time and include it here.
!
! Notice that for multi-processor runs all processors need to
! call the interpolating routine since: 1) We don't know who
! owns the point and, 2) if we don't do this the code
! will hang.  Processors that don't own the point will just
! return 0, and the MPI_ALLREDUCE call will later make sure
! all processors end up with the same value (aux2).
!
! But then we need to be careful. The interpolation is done
! at grid level l-1, so even if a given processor owns the
! point al level l, it might not own it at level l-1 (and
! viceversa). We then need to make sure we don't end up
! injecting the interpolated value to the wrong processor.
! That is the purpose of flag1 (flag2 is just a return
! diagnostic from the interpolation routine that is not
! used here).

! BOUNDARIES ON R: We iterate over ghost zones on r,
! and all values of z for the given refinement box.

! Upper boundary.

  do m=0,imax
     do n=0,Nzbox(box)+ghost-1

!       Figure out (r,z) position for interpolation.

        r0 = rmaxl(box,level) - dble(m)*drl(level)
        z0 = zminl(box,level) + dble(n)*dzl(level)

!       Figure out to which grid point this (r0,z0) values
!       would correspond in the local processor at the
!       fine grid level.

        i = nint((r0-r(1-ghost,0))/drl(level)) + 1 - ghost
        j = nint((z0-z(0,1-ghost))/dzl(level)) + 1 - ghost

!       Notice now that the (i,j) values above might be outside
!       the range of the current processor.  This means that
!       we should not try to access this location as it belongs
!       to another processor.

        if ((i>=1-ghost).and.(i<=Nrl(box,rank)).and.(j>=1-ghost).and.(j<=Nzl(box,rank))) then
           flag1 = .true.
        else
           flag1 = .false.
        end if

!       Now we interpolate all evolving variables at grid level l-1
!       to inject them to grid level l.

        bd = 'rR'

        if (all) then
           include '../auto/bound_interp.inc'
        else
           interpvar => coarsevar
           aux1 = interp(bbox,level-1,r0,z0,flag2)
           call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
           if (flag1) then
              if (border==1) then
                 finevar(i,j) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*finevar_p(i,j)
              else
                 finevar(i,j) = (tl-t0)*(tl-tm1)/((tp -t0)*(tp-tm1))*aux2 &
                              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*finevar_bound_rR(m,j,0) &
                              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*finevar_bound_rR(m,j,1)
              end if
           end if
        end if

     end do
  end do

! Lower boundary when applicable. Same comments
! apply as above, so I don't repeat them here.

  if ((box>0).and.(rbox(box)/=0.d0)) then

     do m=0,imax
        do n=0,Nzbox(box)+ghost-1

           r0 = rminl(box,level) + dble(m)*drl(level)
           z0 = zminl(box,level) + dble(n)*dzl(level)

           i = nint((r0-r(1-ghost,0))/drl(level)) + 1 - ghost
           j = nint((z0-z(0,1-ghost))/dzl(level)) + 1 - ghost

           if ((i>=1-ghost).and.(i<=Nrl(box,rank)).and.(j>=1-ghost).and.(j<=Nzl(box,rank))) then
              flag1 = .true.
           else
              flag1 = .false.
           end if

           bd = 'rL'

           if (all) then
              include '../auto/bound_interp.inc'
           else
              interpvar => coarsevar
              aux1 = interp(bbox,level-1,r0,z0,flag2)
              call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
              if (flag1) then
                 if (border==1) then
                    finevar(i,j) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*finevar_p(i,j)
                 else
                    finevar(i,j) = (tl-t0)*(tl-tm1)/((tp -t0)*(tp-tm1))*aux2 &
                                 + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*finevar_bound_rL(m,j,0) &
                                 + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*finevar_bound_rL(m,j,1)
                 end if
              end if
           end if

        end do
     end do

  end if

! BOUNDARIES ON Z. We iterate over ghost zones on r,
! and all values of z for the given refinement box.

! Upper boundary.

  do n=0,jmax
     do m=0,Nrbox(box)+ghost-1

        r0 = rminl(box,level) + dble(m)*drl(level)
        z0 = zmaxl(box,level) - dble(n)*dzl(level)

        i = nint((r0-r(1-ghost,0))/drl(level)) + 1 - ghost
        j = nint((z0-z(0,1-ghost))/dzl(level)) + 1 - ghost

        if ((i>=1-ghost).and.(i<=Nrl(box,rank)).and.(j>=1-ghost).and.(j<=Nzl(box,rank))) then
           flag1 = .true.
        else
           flag1 = .false.
        end if

        bd = 'zR'

        if (all) then
           include '../auto/bound_interp.inc'
        else
           interpvar => coarsevar
           aux1 = interp(bbox,level-1,r0,z0,flag2)
           call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
           if (flag1) then
              if (border==1) then
                 finevar(i,j) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*finevar_p(i,j)
              else
                 finevar(i,j) = (tl-t0)*(tl-tm1)/((tp -t0)*(tp-tm1))*aux2 &
                              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*finevar_bound_zR(i,n,0) &
                              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*finevar_bound_zR(i,n,1)
              end if
           end if
        end if

     end do
  end do

! Lower boundary when applicable.

  if ((.not.eqsym).or.((box>0).and.(zbox(box)/=0.d0))) then

     do n=0,jmax
        do m=0,Nrbox(box)+ghost-1

           r0 = rminl(box,level) + dble(m)*drl(level)
           z0 = zminl(box,level) + dble(n)*dzl(level)

           i = nint((r0-r(1-ghost,0))/drl(level)) + 1 - ghost
           j = nint((z0-z(0,1-ghost))/dzl(level)) + 1 - ghost

           if ((i>=1-ghost).and.(i<=Nrl(box,rank)).and.(j>=1-ghost).and.(j<=Nzl(box,rank))) then
              flag1 = .true.
           else
              flag1 = .false.
           end if

           bd = 'zL'

           if (all) then
              include '../auto/bound_interp.inc'
           else
              interpvar => coarsevar
              aux1 = interp(bbox,level-1,r0,z0,flag2)
              call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
              if (flag1) then
                 if (border==1) then
                    finevar(i,j) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*finevar_p(i,j)
                 else
                    finevar(i,j) = (tl-t0)*(tl-tm1)/((tp -t0)*(tp-tm1))*aux2 &
                                 + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*finevar_bound_zL(i,n,0) &
                                 + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*finevar_bound_zL(i,n,1)
                 end if
              end if
           end if

        end do
     end do

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine finebound
