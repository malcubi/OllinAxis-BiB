
  subroutine analysis_matter

! ********************************************************
! ***   CALCULATION OF ANALYSIS VARIABLES FOR OUTPUT   ***
! ********************************************************

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical contains


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

  if (contains(mattertype,"complex")) then

!    Norm of complex scalar field: sqrt(phiR**2 + phiI**2).

     if (allocated(complex_phi_norm)) then
        complex_phi_norm = dsqrt(complex_phiR**2 + complex_phiI**2)
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine analysis_matter

