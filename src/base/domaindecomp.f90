!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/domaindecomp.f90,v 1.7 2019/08/28 18:03:39 malcubi Exp $

  subroutine domaindecomp

! ********************************
! ***   DOMAIN DECOMPOSITION   ***
! ********************************

! This routine does the domain decomposition for
! parallel runs. It first needs to figure out how
! many processors to use in each direction.
!
! Then it figures out, for each refinement box, how
! many grid points to allocate to each processor
! in each direction.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer i,j                ! Grid point counters.
  integer box,proc           ! Box and processor counters.
  integer Naux               ! Auxiliary.

  real(8) aux                ! Auxiliary.


! *******************************************************
! ***   FIND NUMBER OF PROCESSORS IN EACH DIRECTION   ***
! *******************************************************

! Output total number of processors.

  if (rank==0) then
     if (size==1) then
        print *, 'Running on 1 processor'
        print *
     else
        write(*,'(A,I3,A)') ' Running on',size,' processors'
        print *
     end if
  end if

! Find out the square of the number of processors.

  aux = sqrt(real(size))

! Special case of one processor.

  if (size==1) then

     nprocr = 1
     nprocz = 1

! General case. Here we find the most efficient way to
! factor the total number of processors in 2 dimensions.
! Notice that if the number of processors is, for example,
! a prime number, this will NOT be very efficient.

  else

     do i=1,int(aux)
        if (real(size)/real(i)==size/i) then
           nprocr=i
           nprocz=size/i
        end if
     end do

  end if

! Message to screen.

  if ((rank==0).and.(size>1)) then
     write(*,'(A,I3,A)') ' with ',nprocr,' processors in r direction'
     write(*,'(A,I3,A)') ' and  ',nprocz,' processors in z direction'
     print *
  end if


! ********************************
! ***   DOMAIN DECOMPOSITION   ***
! ********************************

! Allocate arrays from number of grid points on each box and processor.

  allocate(Nrmaxl(0:Nb),Nzmaxl(0:Nb))

  allocate(Nminl_r(0:Nb,0:size-1),Nmaxl_r(0:Nb,0:size-1),Nrl(0:Nb,0:size-1))
  allocate(Nminl_z(0:Nb,0:size-1),Nmaxl_z(0:Nb,0:size-1),Nzl(0:Nb,0:size-1))

! Arrays for minimum and maximum values of (r,z) for each box and processor.

  allocate(rminl(0:Nb,0:Nlmax),rmaxl(0:Nb,0:Nlmax))
  allocate(zminl(0:Nb,0:Nlmax),zmaxl(0:Nb,0:Nlmax))

! Single processor run.

  if (size==1) then

!    Local size is equal to global size for all boxes.

     Nrl(:,0) = Nrbox(:)
     Nzl(:,0) = Nzbox(:)

!    First point is zero and last point is the local size
!    for all boxes.

     Nminl_r = 0
     Nminl_z = 0

     Nmaxl_r = Nrl
     Nmaxl_z = Nzl

! Parallel run.

! The idea here is that each processor owns grid points with:
!
! i,j = 1,...,N-ghost      (interior points)
!
! with N={Nr,Nz}.
!
! Boundary points, or points belonging to a different processor
! that need to be synchronized, are points such that:
!
!    i,j = 1-ghost,...,0      (lower boundary)
!    i,j = N+1-ghost,...,N    (upper boundary)
!
! For example, for ghost=1, the boundary points are i,j=0
! and i,j=N.  For ghost=2, the boundary points are i,j=-1,0
! and i,j=N-1,N.

! Now figure out the coordinates of the corners for each processor.
! Notice that mod(rank,nprocr) tells me on which "r" bin I am,
! while rank/nprocr tells me on which "z" bin I am.

   else

!    First loop over processors.

     do proc=0,size-1

        i = mod(proc,nprocr) ! r bin
        j = proc/nprocr      ! z bin

!       Now loop over refinement boxes.

        do box=0,Nb

           Naux = Nrbox(box)/nprocr

           Nmaxl_r(box,proc) = (i+1)*Naux
           Nminl_r(box,proc) = Nmaxl_r(box,proc) - Naux

           Naux = Nzbox(box)/nprocz

           Nmaxl_z(box,proc) = (j+1)*Naux
           Nminl_z(box,proc) = Nmaxl_z(box,proc) - Naux

!          Take ghost zones into account, and make sure
!          that the last processor in a given direction
!          always has Nmax=Nbox(box).

           if (i==nprocr-1) then
              Nmaxl_r(box,proc) = Nrbox(box)
           else
              Nmaxl_r(box,proc) = Nmaxl_r(box,proc) + ghost
           end if

           if (j==nprocz-1) then
              Nmaxl_z(box,proc) = Nzbox(box)
           else
              Nmaxl_z(box,proc) = Nmaxl_z(box,proc) + ghost
           end if

!          Now find local number of grid points in each direction.

           Nrl(box,proc) = Nmaxl_r(box,proc) - Nminl_r(box,proc)
           Nzl(box,proc) = Nmaxl_z(box,proc) - Nminl_z(box,proc)

        end do

     end do

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine domaindecomp

