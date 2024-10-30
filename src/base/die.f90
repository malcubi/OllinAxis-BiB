!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/die.f90,v 1.1 2018/08/31 18:02:14 malcubi Exp $

  subroutine die

  use mpi
  use procinfo

! Subroutine for a clean code abort.


! ************************
! ***   FINALIZE MPI   ***
! ************************

  call MPI_FINALIZE(ierr)


! ***************
! ***  STOP   ***
! ***************

  stop


! ***************
! ***   END   ***
! ***************

  end subroutine die

