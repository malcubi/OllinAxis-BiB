
! *********************************************
! ***   STUBS FOR SOME BASIC MPI ROUTINES   ***
! *********************************************

! This file should only be compiled when there is no MPI installed
! in the system.  It will just make sure that the code compiles
! in serial mode without complaining, basically by bypassing all
! MPI calls.

! ***   MPI_INIT   ***

  subroutine MPI_INIT(ierr)
  implicit none
  integer ierr
  end subroutine MPI_INIT

! ***   MPI_FINALIZE   ***

  subroutine MPI_FINALIZE(ierr)
  implicit none
  integer ierr
  end subroutine MPI_FINALIZE

! ***   MPI_COMM_SIZE   ***

  subroutine MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  implicit none
  integer MPI_COMM_WORLD
  integer size,ierr
  size = 1
  end subroutine MPI_COMM_SIZE

! ***   MPI_COMM_RANK   ***

  subroutine MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
  implicit none
  integer MPI_COMM_WORLD
  integer rank,ierr
  rank = 0
  end subroutine MPI_COMM_RANK

! ***   MPI_Allreduce   ***

  subroutine MPI_Allreduce(local,global,tag,MPI_TYPE,MPI_OP,MPI_COMM_WORLD,ierr)
  implicit none
  integer MPI_COMM_WORLD
  integer MPI_TYPE,MPI_OP
  integer tag,ierr
  real(8) local,global
  global = local
  end subroutine MPI_Allreduce

! ***   MPI_SEND   ***

  subroutine MPI_SEND(buffer,N,MPI_TYPE,proc,tag,MPI_COMM_WORLD,ierr)
  implicit none
  integer MPI_COMM_WORLD
  integer MPI_TYPE
  integer N,proc,tag,ierr
  real(8) buffer(1:N)
  end subroutine MPI_SEND

! ***   MPI_RECV   ***

  subroutine MPI_RECV(buffer,N,MPI_TYPE,proc,tag,MPI_COMM_WORLD,ierr)
  implicit none
  integer MPI_COMM_WORLD
  integer MPI_TYPE
  integer N,proc,tag,ierr
  real(8) buffer(1:N)
  end subroutine MPI_RECV

! ***   MPI_ISEND   ***

  subroutine MPI_ISEND(buffer,N,MPI_TYPE,proc,tag,MPI_COMM_WORLD,request,ierr)
  implicit none
  integer MPI_COMM_WORLD
  integer MPI_TYPE
  integer N,proc,tag,ierr
  integer request(:)
  real(8) buffer(1:N)
  end subroutine MPI_ISEND

! ***   MPI_ISEND   ***

  subroutine MPI_IRECV(buffer,N,MPI_TYPE,proc,tag,MPI_COMM_WORLD,request,ierr)
  implicit none
  integer MPI_COMM_WORLD
  integer MPI_TYPE
  integer N,proc,tag,ierr
  integer request(:)
  real(8) buffer(1:N)
  end subroutine MPI_IRECV

! ***   MPI_BCAST   ***

  subroutine MPI_BCAST(buffer,N,MPI_TYPE,proc,MPI_COMM_WORLD,ierr)
  implicit none
  integer MPI_COMM_WORLD
  integer MPI_TYPE
  integer N,proc,ierr
  real(8) buffer(1:N)
  end subroutine MPI_BCAST

! ***   MPI_BARRIER   ***

  subroutine MPI_BARRIER(MPI_COMM_WORLD,ierr)
  implicit none
  integer MPI_COMM_WORLD
  integer ierr
  end subroutine MPI_BARRIER

! ***   MPI_WAITALL   ***

  subroutine MPI_WAITALL(N,request,MPI_STATUS,ierr)
  implicit none
  integer MPI_STATUS
  integer N,ierr
  integer request(1:N)
  end subroutine MPI_WAITALL
