
! *****************************************
! ***   GRAVITATIONAL WAVE EXTRACTION   ***
! *****************************************

! Here we define the main gravitational wave extraction method:
!
! wavextract    Calls the appropriate wave extraction subroutines
!               (Moncrief-Zerilli, Weyl scalar Psi4, or both)
!               based on user configuration.


  subroutine wavextract

! Routine for gravitational wave extraction. It dispatches execution
! to the Moncrief-Zerilli formalism, the Newman-Penrose (Weyl scalar Psi4)
! formalism, or both, according to the wavextract_method parameter.

  use mpi
  use param
  use arrays
  use procinfo
  use harmonix

! Extra variables.

  implicit none
  
  logical contains


! ************************************************
! ***   CALL SELECTED WAVE EXTRACTION METHOD   ***
! ************************************************

! Execute both Moncrief-Zerilli and Weyl scalar Psi4 extractions.

  if (contains(wavextract_method,"both")) then

     call extract_moncrief
     call extract_weyl

! Execute polar Moncrief-Zerilli wave extraction only.

  else if (contains(wavextract_method,"moncrief")) THEN

     call extract_moncrief

! Execute Weyl scalar Psi4 (Newman-Penrose) wave extraction only.

  else if (contains(wavextract_method,"weyl")) THEN

     call extract_weyl

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine wavextract
