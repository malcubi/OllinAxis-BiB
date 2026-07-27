! *******************************************************
! ***   GRAVITATIONAL WAVE EXTRACTION                 ***
! *******************************************************

! Here I define the main gravitational wave extraction:
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


  implicit none
  
  logical contains

! Check selected wave extraction method
  
  IF (contains(wavextract_method, "both")) THEN

!    Execute both Moncrief-Zerilli and Weyl scalar Psi4 extractions.

     CALL extract_moncrief
     CALL extract_weyl

  ELSE IF (contains(wavextract_method, "moncrief")) THEN

!    Execute polar Moncrief-Zerilli wave extraction only.

     CALL extract_moncrief

  ELSE IF (contains(wavextract_method, "weyl")) THEN

!    Execute Weyl scalar Psi4 (Newman-Penrose) wave extraction only.

     CALL extract_weyl

  ELSE IF (contains(wavextract_method, "none")) THEN

!    No wave extraction requested; do nothing.

  END IF

! ***************
! ***   END   ***
! ***************

end subroutine wavextract
