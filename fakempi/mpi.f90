!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/fakempi/mpi.f90,v 1.3 2019/08/28 02:23:20 malcubi Exp $

! *********************************
! ***   MPI FAKE DECLARATIONS   ***
! *********************************

! This file should only be compiled when there is no MPI installed
! in the system.  It will just make sure that the code compiles
! in serial mode without complaining, basically by bypassing all
! MPI calls.

  module mpi

  integer MPI_COMM_WORLD
  integer MPI_STATUS_SIZE

  integer MPI_LOGICAL

  integer MPI_INT
  integer MPI_INTEGER

  integer MPI_REAL8
  integer MPI_DOUBLE_PRECISION

  integer MPI_SUM
  integer MPI_MAX
  integer MPI_MIN
  integer MPI_LOR


! ***************
! ***   END   ***
! ***************

  end module mpi



