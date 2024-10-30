!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/distribute.f90,v 1.2 2019/12/03 20:19:18 malcubi Exp $

  subroutine distribute(varl,varg)

! *******************************************************
! ***   DISTRIBUTE ONE LARGE ARRAY AMONG PROCESSORS   ***
! *******************************************************

! This is a small routine that only makes sense for multi-processor
! runs.  It distributes a large array "varg" of dimension( Nrtotal,Nztotal)
! owned by processor 0 among local arrays "varl" on all processors.

! Include modules.

  use mpi
  use param
  use procinfo

! Extra variables.

  implicit none

  real(8), dimension (1-ghost:Nrmax,1-ghost:Nzmax) :: varl       ! Local array.
  real(8), dimension (1-ghost:Nrtotal,1-ghost:Nztotal) :: varg   ! Global array.

  if (rank==0) then
     print *, 'Distributing 2D data from processor 0 to other processors not yet implemented.'
     print *, 'Aborting (subroutine distribute.f90)'
     print *
  end if
  call die


! ***********************
! ***   PROCESSOR 0   ***
! ***********************


! ****************************
! ***   OTHER PROCESSORS   ***
! ****************************


! ***********************
! ***   SYNCHRONIZE   ***
! ***********************


! ***************
! ***   END   ***
! ***************

  end subroutine distribute
