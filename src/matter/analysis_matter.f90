
  subroutine analysis_matter(box,level)

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

  integer box,level


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

  if (contains(mattertype,"complex")) then

!    Norm of complex scalar field: sqrt(phiR**2 + phiI**2).

     if (associated(grid(box,level)%complex_phi_norm)) then
        complex_phi_norm = dsqrt(complex_phiR**2 + complex_phiI**2)
     end if

!    The integrated boson charge is:
!
!                   /             
!    complex_NB  =  |  complex_Bdens dV
!                   /
!
!    with dV the physical volumen element.

     if (associated(grid(box,level)%complex_NB)) then
        !call bosonintegral
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine analysis_matter

