!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/sync.f90,v 1.4 2019/10/04 15:49:36 malcubi Exp $

  subroutine sync(syncvar)

! ********************************
! ***   SYNCHRONIZE VARIABLE   ***
! ********************************

! Load modules.

  use mpi
  use param,  only: ghost,Nr,Nz,Nrmax,Nzmax
  use procinfo

! Extra variables.

  implicit none

  integer i,j,p,Naux
  integer i0
  integer status(MPI_STATUS_SIZE)

  real(8), dimension (1-ghost:Nrmax,1-ghost:Nzmax) :: syncvar
  real(8), dimension (1-ghost:Nrmax) :: varr
  real(8), dimension (1-ghost:Nzmax) :: varz


! ***************************
! ***   SYNCHRONIZATION   ***
! ***************************

! We only need to synchronize if we have
! more than one processor.

  if (size>1) then

!    ***********************
!    ***   R DIRECTION   ***
!    ***********************

     if (nprocr>1) then

        Naux = Nzmax + ghost

!       Send information to processor on the left
!       (only if we are not at left boundary).

        if (mod(rank,nprocr)/=0) then
           p = rank-1
           do i=1,ghost
              varz = syncvar(i,:)
              call MPI_SEND(varz,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
           end do
        end if

!       Send information to processor on the right
!       (only if we are not at right boundary).

        if (mod(rank+1,nprocr)/=0) then
           p = rank+1
           do i=1,ghost
              varz = syncvar(Nr-2*ghost+i,:)
              call MPI_SEND(varz,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
           end do
        end if

!       Receive information from processor on the right
!       (only if we are not at right boundary).

        if (mod(rank+1,nprocr)/=0) then
           p = rank+1
           do i=1,ghost
              call MPI_RECV(varz,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              syncvar(Nr-ghost+i,:) = varz
           end do
        end if

!       Receive information from processor on the left
!       (only if we are not at left boundary).

        if (mod(rank,nprocr)/=0) then
           p = rank-1
           do i=1,ghost
              call MPI_RECV(varz,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              syncvar(i-ghost,:) = varz
           end do
        end if

     end if


!    ***********************
!    ***   Z DIRECTION   ***
!    ***********************

     if (nprocz>1) then

        Naux = Nrmax + ghost

!       Send information to processor below us
!       (only if we are not at the lower boundary).

        if (rank>=nprocr) then
           p = rank - nprocr
           do j=1,ghost
              varr = syncvar(:,j)
              call MPI_SEND(varr,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
           end do
        end if

!       Send information to processor above us
!       (only if we are not at the upper boundary).

        if (rank<size-nprocr) then
           p = rank + nprocr
           do j=1,ghost
              varr = syncvar(:,Nz-2*ghost+j)
              call MPI_SEND(varr,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
           end do
        end if

!       Receive information from processor above us
!       (only if we are not at the upper boundary).

        if (rank<size-nprocr) then
           p = rank + nprocr
           do j=1,ghost
              call MPI_RECV(varr,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              syncvar(:,Nz-ghost+j) = varr
           end do
        end if

!       Receive information from processor below us
!       (only if we are not at the lower boundary).

        if (rank>=nprocr) then
           p = rank - nprocr
           do j=1,ghost
              call MPI_RECV(varr,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              syncvar(:,j-ghost) = varr
           end do
        end if

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sync







  subroutine syncr(syncvar)

! *******************************************
! ***   SYNCHRONIZE VARIABLE ACROSS RHO   ***
! *******************************************

! Load modules.

  use mpi
  use param,  only: ghost,Nr,Nrmax,Nzmax
  use procinfo

! Extra variables.

  implicit none

  integer i,p,Naux
  integer status(MPI_STATUS_SIZE)

  real(8), dimension (1-ghost:Nrmax,1-ghost:Nzmax) :: syncvar
  real(8), dimension (1-ghost:Nzmax) :: varz


! ***********************
! ***   R DIRECTION   ***
! ***********************

! We only need to synchronize if we have
! more than one processor across rho.

  if (nprocr>1) then

     Naux = Nzmax + ghost

!    Send information to processor on the left
!    (only if we are not on axis).

     if (mod(rank,nprocr)/=0) then
        p = rank-1
        do i=1,ghost
           varz = syncvar(i,:)
           call MPI_SEND(varz,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
        end do
     end if

!    Send information to processor on the right
!    (only if we are not at outer boundary).

     if (mod(rank+1,nprocr)/=0) then
        p = rank+1
        do i=1,ghost
           varz = syncvar(Nr-2*ghost+i,:)
           call MPI_SEND(varz,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
        end do
     end if

!    Receive information from processor on the right
!    (only if we are not at outer boundary).

     if (mod(rank+1,nprocr)/=0) then
        p = rank+1
        do i=1,ghost
           call MPI_RECV(varz,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           syncvar(Nr-ghost+i,:) = varz
        end do
     end if

!    Receive information from processor on the left
!    (only if we are not on axis).

     if (mod(rank,nprocr)/=0) then
        p = rank-1
        do i=1,ghost
           call MPI_RECV(varz,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           syncvar(i-ghost,:) = varz
        end do
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine syncr







  subroutine syncz(syncvar)

! *****************************************
! ***   SYNCHRONIZE VARIABLE ACROSS Z   ***
! *****************************************

! Load modules.

  use mpi
  use param,  only: ghost,Nz,Nrmax,Nzmax
  use procinfo

! Extra variables.

  implicit none

  integer j,p,Naux
  integer status(MPI_STATUS_SIZE)

  real(8), dimension (1-ghost:Nrmax,1-ghost:Nzmax) :: syncvar
  real(8), dimension (1-ghost:Nrmax) :: varr


! ***********************
! ***   Z DIRECTION   ***
! ***********************

! We only need to synchronize if we have
! more than one processor across z.

  if (nprocz>1) then

     Naux = Nrmax + ghost

!    Send information to processor below us
!    (only if we are not at the lower boundary).

     if (rank>=nprocr) then
        p = rank - nprocr
        do j=1,ghost
           varr = syncvar(:,j)
           call MPI_SEND(varr,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
        end do
     end if

!    Send information to processor above us
!    (only if we are not at the upper boundary).

     if (rank<size-nprocr) then
        p = rank + nprocr
        do j=1,ghost
           varr = syncvar(:,Nz-2*ghost+j)
           call MPI_SEND(varr,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
        end do
     end if

!    Receive information from processor above us
!    (only if we are not at the upper boundary).

     if (rank<size-nprocr) then
        p = rank + nprocr
        do j=1,ghost
           call MPI_RECV(varr,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           syncvar(:,Nz-ghost+j) = varr
        end do
     end if

!    Receive information from processor below us
!    (only if we are not at the lower boundary).

     if (rank>=nprocr) then
        p = rank - nprocr
        do j=1,ghost
           call MPI_RECV(varr,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           syncvar(:,j-ghost) = varr
        end do
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine syncz


