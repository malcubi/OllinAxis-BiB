!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/syncboxes.f90,v 1.9 2019/10/29 19:39:52 malcubi Exp $

  subroutine syncboxes(b1,b2,level,all)

! ***************************************
! ***   CHECK FOR BOX INTERSECTIONS   ***
! ***************************************

! Check if boxes on the same level intersect, and if they
! do synchronize the data.
!
! When we get here we have selected a pair of boxes (b1,b2)
! that exist at the same grid level to see if they intersect.
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

  integer b1,b2              ! Pair of boxes.
  integer level              ! Grid level.
  integer i,j,m,n            ! Grid point counters.
  integer imax,jmax          ! Size of boundary region.

  real(8) r0,z0,interp       ! For interpolation.
  real(8) deltar,deltaz      ! Safe distance to boundaries.
  real(8) aux1,aux2


! ***********************************
! ***   CHECK FOR INTERSECTIONS   ***
! ***********************************

! Simple check of cases where no intersection is possible.

  if ((rmaxl(b1,level)<rminl(b2,level)).or.(rmaxl(b2,level)<rminl(b1,level)).or. &
      (zmaxl(b1,level)<zminl(b2,level)).or.(zmaxl(b2,level)<zminl(b1,level))) return

! Size of region to synchronize. It can be as large as
! 3*ghost because of the internal iterations of rk4.
 
  imax = 3*ghost
  jmax = 3*ghost

! Safe distance to boundaries of box b2.

  deltar = dble(3*ghost)*drl(level)
  deltaz = dble(3*ghost)*dzl(level)

! Upper r boundary.

  do m=0,imax
     do n=0,Nzbox(b1)+ghost-1

!       Figure out (r,z) position for interpolation.

        r0 = rmaxl(b1,level) - dble(m)*drl(level)
        z0 = zminl(b1,level) + dble(n)*dzl(level)

!       We will only interpolate if point (r0,z0)
!       is well inside box b2.

        if ((r0>=rminl(b2,level)+deltar).and.(r0<=rmaxl(b2,level)-deltar).and. &
            (z0>=zminl(b2,level)+deltaz).and.(z0<=zmaxl(b2,level)-deltaz)) then

!          Figure out to which grid point this (r0,z0) values
!          would correspond in the local processor for box b1.

           i = nint((r0-grid(b1,level)%r(1-ghost,0))/drl(level)) + 1 - ghost
           j = nint((z0-grid(b1,level)%z(0,1-ghost))/dzl(level)) + 1 - ghost

!          Notice now that the (i,j) values above might be outside
!          the range of the current processor in box b1.
!          This means that we should not try to access this location
!          as it belongs to another processor.

           if ((i>=1-ghost).and.(i<=Nrl(b1,rank)).and.(j>=1-ghost).and.(j<=Nzl(b1,rank))) then
              flag1 = .true.
           else
              flag1 = .false.
           end if

!          Now we interpolate all evolving variables from box b2 to
!          inject them into box b1. Since here I don't know how many
!          variables are evolving, I generate the necessary code
!          automatically at compile time and include it here.

           if (all) then
              include '../auto/intersect_interp.inc'
           else
              interpvar => box2var
              aux1 = interp(b2,level,r0,z0,flag2)
              call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
              if (flag1) then
                 box1var(i,j) = aux2
              end if
           end if

        end if

     end do
  end do

! Lower r boundary when applicable.  Same comments
! apply as above, so I don't repeat them here.

  if ((b1>0).and.(rbox(b1)/=0.d0)) then

     do m=0,imax
        do n=0,Nzbox(b1)+ghost-1

           r0 = rminl(b1,level) + dble(m)*drl(level)
           z0 = zminl(b1,level) + dble(n)*dzl(level)

           if ((r0>=rminl(b2,level)+deltar).and.(r0<=rmaxl(b2,level)-deltar).and. &
               (z0>=zminl(b2,level)+deltaz).and.(z0<=zmaxl(b2,level)-deltaz)) then

              i = nint((r0-grid(b1,level)%r(1-ghost,0))/drl(level)) + 1 - ghost
              j = nint((z0-grid(b1,level)%z(0,1-ghost))/dzl(level)) + 1 - ghost

              if ((i>=1-ghost).and.(i<=Nrl(b1,rank)).and.(j>=1-ghost).and.(j<=Nzl(b1,rank))) then
                 flag1 = .true.
              else
                 flag1 = .false.
              end if

              if (all) then
                 include '../auto/intersect_interp.inc'
              else
                 interpvar => box2var
                 aux1 = interp(b2,level,r0,z0,flag2)
                 call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
                 if (flag1) then
                    box1var(i,j) = aux2
                 end if
              end if

           end if

        end do
     end do

  end if

! Upper z boundary when applicable.

  do n=0,jmax
     do m=0,Nrbox(b1)+ghost-1

        r0 = rminl(b1,level) + dble(m)*drl(level)
        z0 = zmaxl(b1,level) - dble(n)*dzl(level)

        if ((r0>=rminl(b2,level)+deltar).and.(r0<=rmaxl(b2,level)-deltar).and. &
            (z0>=zminl(b2,level)+deltaz).and.(z0<=zmaxl(b2,level)-deltaz)) then

           i = nint((r0-grid(b1,level)%r(1-ghost,0))/drl(level)) + 1 - ghost
           j = nint((z0-grid(b1,level)%z(0,1-ghost))/dzl(level)) + 1 - ghost

           if ((i>=1-ghost).and.(i<=Nrl(b1,rank)).and.(j>=1-ghost).and.(j<=Nzl(b1,rank))) then
              flag1 = .true.
           else
              flag1 = .false.
           end if

           if (all) then
              include '../auto/intersect_interp.inc'
           else
              interpvar => box2var
              aux1 = interp(b2,level,r0,z0,flag2)
              call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
              if (flag1) then
                 box1var(i,j) = aux2
              end if
           end if

        end if

     end do
  end do

! Lower z boundary.

  if ((.not.eqsym).or.((b1>0).and.(zbox(b1)/=0.d0))) then

     do n=0,jmax
        do m=0,Nrbox(b1)+ghost-1

           r0 = rminl(b1,level) + dble(m)*drl(level)
           z0 = zminl(b1,level) + dble(n)*dzl(level)

           if ((r0>=rminl(b2,level)+deltar).and.(r0<=rmaxl(b2,level)-deltar).and. &
               (z0>=zminl(b2,level)+deltaz).and.(z0<=zmaxl(b2,level)-deltaz)) then

              i = nint((r0-grid(b1,level)%r(1-ghost,0))/drl(level)) + 1 - ghost
              j = nint((z0-grid(b1,level)%z(0,1-ghost))/dzl(level)) + 1 - ghost

              if ((i>=1-ghost).and.(i<=Nrl(b1,rank)).and.(j>=1-ghost).and.(j<=Nzl(b1,rank))) then
                 flag1 = .true.
              else
                 flag1 = .false.
              end if

              if (all) then
                 include '../auto/intersect_interp.inc'
              else
                 interpvar => box2var
                 aux1 = interp(b2,level,r0,z0,flag2)
                 call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
                 if (flag1) then
                    box1var(i,j) = aux2
                 end if
              end if

           end if

        end do
     end do

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine syncboxes
