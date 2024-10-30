!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/restrict.f90,v 1.22 2021/02/19 20:27:39 malcubi Exp $

  subroutine restrict(box,level,all)

! ********************
! ***   RESTRICT   ***
! ********************

! Restrict data from fine grid into coarse grid.
!
! The logical flag "all", when true, means we restrict
! all evolving arrays.  When false we restrict a single array.
!
! The logical flag "all", when true, means we restrict all
! evolving arrays.  When false we restrict a single array.
!
! When we apply restriction to a single array, the corresponding
! coarse array and fine array must be assigned to the pointers
! "coarsevar" and "finevar" BEFORE calling the routine.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Declare variables.
  
  implicit none

  logical contains
  logical all                ! Do we apply the boundary to all evolving arrays?
  logical flag1,flag2        ! Interpolating flags.

  integer i,j,m,n            ! Grid point counters.
  integer box,level          ! Box number and refinement level counters.
  integer bbox               ! Auxiliary box identificator.

  real(8) r0,z0,interp       ! For interpolation.
  real(8) deltar,deltaz      ! Safe distance to boundaries.
  real(8) aux1,aux2


! ************************
! ***   SANITY CHECK   ***
! ************************

! This subroutine should never be called for level=0.

  if (level==0) then
     print *
     print *, 'The subroutine "restrict" should never be called with level=0!'
     print *, 'Aborting! (subroutine restrict)'
     print *
     call die
  end if

  ! write(*,"(A,I2,A,I2)") ' Restricting level ',l,' to level ',level-1


! *******************************************************
! ***   FIGURE OUT ON WHICH COARSE GRID TO RESTRICT   ***
! *******************************************************

! Normally we restrict data to the coarse grid on the
! same box we are in.  But if we are at level 1 we
! always restrict to box 0 (the base grid).

  if (level==1) then
     bbox = 0
  else
     bbox = box
  end if


! ****************************************
! ***   RESTRICT DATA TO COARSE GRID   ***
! ****************************************

! Safe distance to boundaries on fine grid.
! We need to be suficiently far, as otherwise
! an instability seems to appear.
!
! We add a small number to prevent changes
! from round-off error.

  deltar = dble(ghost)*drl(level) + 1.d-10
  deltaz = dble(ghost)*dzl(level) + 1.d-10

! All processors loop over ALL points in the coarse time level,
! even if they don't own them.  This is to make sure that all
! call the interpolation routine for the same point at the
! same time.

  do n=0,Nzbox(bbox)+ghost-1
     do m=0,Nrbox(bbox)+ghost-1

!       Figure out (r,z) position for interpolation.

        r0 = rminl(bbox,level-1) + dble(m)*drl(level-1)
        z0 = zminl(bbox,level-1) + dble(n)*dzl(level-1)

!       We will only interpolate if point (r0,z0) is well
!       inside the fine grid.

        if ((r0>rminl(box,level)+deltar).and.(r0<rmaxl(box,level)-deltar).and. &
            (z0>zminl(box,level)+deltaz).and.(z0<zmaxl(box,level)-deltaz)) then

!          Figure out to which grid point this (r0,z0) values
!          would correspond in the local processor at the
!          coarse grid level.

           i = nint((r0-grid(bbox,level-1)%r(1-ghost,0))/drl(level-1)) + 1 - ghost
           j = nint((z0-grid(bbox,level-1)%z(0,1-ghost))/dzl(level-1)) + 1 - ghost

!          Notice now that the (i,j) values above might be outside
!          the range of the current processor in the coarse grid.
!          This means that we should not try to access this location
!          as it belongs to another processor.

           if ((i>=1-ghost).and.(i<=Nrl(bbox,rank)).and.(j>=1-ghost).and.(j<=Nzl(bbox,rank))) then
              flag1 = .true.
           else
              flag1 = .false.
           end if

!          Now we interpolate all evolving variables at grid level l to
!          inject them into grid level-1. Since here I don't know how
!          many variables are evolving, I generate the necessary code
!          automatically at compile time and include it here.

           if (all) then
              include '../auto/restrict_interp.inc'
           else
              interpvar => finevar
              aux1 = interp(box,level,r0,z0,flag2)
              call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
              if (flag1) then
                 coarsevar(i,j) = aux2
              end if
           end if

        end if

     end do
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine restrict
