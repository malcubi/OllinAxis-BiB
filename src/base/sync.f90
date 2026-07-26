
  subroutine sync(syncvar)

! ********************************
! ***   SYNCHRONIZE VARIABLE   ***
! ********************************

! This routine has been improved with ChatGPT.

! Load modules.

  use mpi
  use param
  use procinfo

! Extra variables.

  implicit none

  integer :: i,j,p,Naux
  integer :: nreq
  integer :: request(4*ghost)
  integer :: status(MPI_STATUS_SIZE,4*ghost)
  integer, parameter :: TAG = 1

  real(8), intent(inout) :: syncvar(1-ghost:Nrmax,1-ghost:Nzmax)
  real(8) :: varz(1-ghost:Nzmax,4*ghost)


! ***************************
! ***   SYNCHRONIZATION   ***
! ***************************

! We only need to synchronize if we have
! more than one processor.


! ***********************
! ***   R DIRECTION   ***
! ***********************

  if (nprocr>1) then

     Naux = Nzmax + ghost
     nreq = 0


!    --------------------------------------------------
!    Pack send buffers BEFORE posting nonblocking sends
!    --------------------------------------------------

!    Left boundary.

     if (mod(rank,nprocr) /= 0) then
        do i=1,ghost
           varz(:,i) = syncvar(i,:)
        end do
     end if

!       Right boundary.

     if (mod(rank+1,nprocr) /= 0) then
        do i=1,ghost
           varz(:,ghost+i) = syncvar(Nr-2*ghost+i,:)
        end do
     end if


!    -------------
!    Post receives
!    -------------

!    Receive from right.

     if (mod(rank+1,nprocr) /= 0) then
        p = rank + 1
        do i=1,ghost
           nreq = nreq + 1
           call MPI_IRECV(varz(:,2*ghost+i), Naux, MPI_REAL8, &
                          p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
        end do
     end if

!    Receive from left.

     if (mod(rank,nprocr) /= 0) then
        p = rank - 1
        do i=1,ghost
           nreq = nreq + 1
           call MPI_IRECV(varz(:,3*ghost+i), Naux, MPI_REAL8, &
                       p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
           end do
     end if


!    ----------
!    Post sends
!    ----------

!    Send to left.

     if (mod(rank,nprocr) /= 0) then
        p = rank - 1
        do i=1,ghost
           nreq = nreq + 1
           call MPI_ISEND(varz(:,i), Naux, MPI_REAL8, &
                          p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
        end do
     end if

!    Send to right.

     if (mod(rank+1,nprocr) /= 0) then
        p = rank + 1
        do i=1,ghost
           nreq = nreq + 1
           call MPI_ISEND(varz(:,ghost+i), Naux, MPI_REAL8, &
                          p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
        end do
     end if


!    --------------------
!    Finish communication
!    --------------------

     if (nreq > 0) then
        call MPI_WAITALL(nreq, request, status, ierr)
     end if


!    --------------------
!    Unpack received data
!    --------------------

!    From right.

     if (mod(rank+1,nprocr) /= 0) then
        do i=1,ghost
           syncvar(Nr-ghost+i,:) = varz(:,2*ghost+i)
        end do
     end if

!    From left.

    if (mod(rank,nprocr) /= 0) then
        do i=1,ghost
           syncvar(i-ghost,:) = varz(:,3*ghost+i)
        end do
     end if

  end if


! ***********************
! ***   Z DIRECTION   ***
! ***********************

  if (nprocz>1) then

     Naux = Nrmax + ghost
     nreq = 0


!    -------------------
!    Post receives first
!    -------------------

!    Receive from processor above.

     if (rank < size-nprocr) then
        p = rank + nprocr
        do j=1,ghost
           nreq = nreq + 1
           call MPI_IRECV(syncvar(:,Nz-ghost+j), Naux, MPI_REAL8, &
                          p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
        end do
     end if

!    Receive from processor below.

     if (rank >= nprocr) then
        p = rank - nprocr
        do j=1,ghost
           nreq = nreq + 1
           call MPI_IRECV(syncvar(:,j-ghost), Naux, MPI_REAL8, &
                          p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
        end do
     end if


!    -----------
!    Post sends
!    -----------

!    Send to processor below.

     if (rank >= nprocr) then
        p = rank - nprocr
        do j=1,ghost
           nreq = nreq + 1
           call MPI_ISEND(syncvar(:,j), Naux, MPI_REAL8, &
                          p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
        end do
     end if

!    Send to processor above.

     if (rank < size-nprocr) then
        p = rank + nprocr
        do j=1,ghost
           nreq = nreq + 1
           call MPI_ISEND(syncvar(:,Nz-2*ghost+j), Naux, MPI_REAL8, &
                          p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
        end do
     end if


!    ---------------------------
!    Wait for all communications
!    ---------------------------

     if (nreq > 0) then
        call MPI_WAITALL(nreq, request, status, ierr)
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

! This routine has been improved mith ChatGPT.

! Load modules.

  use mpi
  use param
  use procinfo

! Extra variables.

  implicit none

  integer :: i,p,Naux
  integer :: nreq
  integer :: request(4*ghost)
  integer :: status(MPI_STATUS_SIZE,4*ghost)
  integer, parameter :: TAG = 1

  real(8), intent(inout) :: syncvar(1-ghost:Nrmax,1-ghost:Nzmax)
  real(8) :: varz(1-ghost:Nzmax,4*ghost)


! ***********************
! ***   R DIRECTION   ***
! ***********************

  if (nprocr<=1) return

  Naux = Nzmax + ghost
  nreq = 0


! --------------------------------------------------
! Pack send buffers BEFORE posting nonblocking sends
! --------------------------------------------------

! Left boundary.

  if (mod(rank,nprocr) /= 0) then
     do i=1,ghost
        varz(:,i) = syncvar(i,:)
     end do
  end if

! Right boundary.

  if (mod(rank+1,nprocr) /= 0) then
     do i=1,ghost
        varz(:,ghost+i) = syncvar(Nr-2*ghost+i,:)
     end do
  end if


! -------------
! Post receives
! -------------

! Receive from right.

  if (mod(rank+1,nprocr) /= 0) then
     p = rank + 1
     do i=1,ghost
        nreq = nreq + 1
        call MPI_IRECV(varz(:,2*ghost+i), Naux, MPI_REAL8, &
                       p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
     end do
  end if

! Receive from left.

  if (mod(rank,nprocr) /= 0) then
     p = rank - 1
     do i=1,ghost
        nreq = nreq + 1
        call MPI_IRECV(varz(:,3*ghost+i), Naux, MPI_REAL8, &
                       p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
     end do
  end if


! ----------
! Post sends
! ----------

! Send to left.

  if (mod(rank,nprocr) /= 0) then
     p = rank - 1
     do i=1,ghost
        nreq = nreq + 1
        call MPI_ISEND(varz(:,i), Naux, MPI_REAL8, &
                       p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
     end do
  end if

! Send to right.

  if (mod(rank+1,nprocr) /= 0) then
     p = rank + 1
     do i=1,ghost
        nreq = nreq + 1
        call MPI_ISEND(varz(:,ghost+i), Naux, MPI_REAL8, &
                       p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
     end do
  end if


! --------------------
! Finish communication
! --------------------

  if (nreq > 0) then
     call MPI_WAITALL(nreq, request, status, ierr)
  end if


! --------------------
! Unpack received data
! --------------------

! From right.

  if (mod(rank+1,nprocr) /= 0) then
     do i=1,ghost
        syncvar(Nr-ghost+i,:) = varz(:,2*ghost+i)
     end do
  end if

! From left.

  if (mod(rank,nprocr) /= 0) then
     do i=1,ghost
        syncvar(i-ghost,:) = varz(:,3*ghost+i)
     end do
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine syncr








  subroutine syncz(syncvar)

! *****************************************
! ***   SYNCHRONIZE VARIABLE ACROSS Z   ***
! *****************************************

! This routine has been improved with ChatGPT.

! Load modules.

  use mpi
  use param
  use procinfo

! Extra variables.

  implicit none

  integer, parameter :: TAG = 1

  integer :: Naux
  integer :: nreq
  integer :: j, p

  integer :: request(4*ghost)
  integer :: status(MPI_STATUS_SIZE,4*ghost)

  real(8), intent(inout) :: syncvar(1-ghost:Nrmax,1-ghost:Nzmax)


! ***********************
! ***   Z DIRECTION   ***
! ***********************

  if (nprocz <= 1) return

  Naux = Nrmax + ghost
  nreq = 0


! -------------------
! Post receives first
! -------------------

! Receive from processor above.

  if (rank < size-nprocr) then
     p = rank + nprocr
     do j=1,ghost
        nreq = nreq + 1
        call MPI_IRECV(syncvar(:,Nz-ghost+j), Naux, MPI_REAL8, &
                       p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
     end do
  end if

! Receive from processor below.

  if (rank >= nprocr) then
     p = rank - nprocr
     do j=1,ghost
        nreq = nreq + 1
        call MPI_IRECV(syncvar(:,j-ghost), Naux, MPI_REAL8, &
                       p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
     end do
  end if


! ----------
! Post sends
! ----------

! Send to processor below.

  if (rank >= nprocr) then
     p = rank - nprocr
     do j=1,ghost
        nreq = nreq + 1
        call MPI_ISEND(syncvar(:,j), Naux, MPI_REAL8, &
                       p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
     end do
  end if

! Send to processor above.

  if (rank < size-nprocr) then
     p = rank + nprocr
     do j=1,ghost
        nreq = nreq + 1
        call MPI_ISEND(syncvar(:,Nz-2*ghost+j), Naux, MPI_REAL8, &
                       p, TAG, MPI_COMM_WORLD, request(nreq), ierr)
     end do
  end if


! ---------------------------
! Wait for all communications
! ---------------------------

  if (nreq > 0) then
     call MPI_WAITALL(nreq, request, status, ierr)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine syncz








  subroutine syncrold(syncvar)

! *******************************************
! ***   SYNCHRONIZE VARIABLE ACROSS RHO   ***
! *******************************************

! Original version of the routine.
! I fimproved it with ChatGPT.

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

  end subroutine syncrold








  subroutine synczold(syncvar)

! *****************************************
! ***   SYNCHRONIZE VARIABLE ACROSS Z   ***
! *****************************************

! Original version of the routine.
! It hangs for very large Nrmax,Nzmax.
! I fixed it and improved with ChatGPT.

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

  end subroutine synczold






