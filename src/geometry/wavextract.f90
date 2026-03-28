
  subroutine wavextract

! *****************************************
! ***   GRAVITATIONAL WAVE EXTRACTION   ***
! *****************************************

! Routine for gravitational wave extraction using
! the Newmann-Penrose formalism and the multipolar
! decomposition of the Weyl scalar Psi4.
! Or the Moncrief-Zerilli formalism and the multipolar
! decomposition of the Q factor.

  use mpi
  use param
  use arrays
  use procinfo

  implicit none
  
  logical contains

   ! Extracción con Moncrief
   if (contains(wavextract_method, "moncrief") .or. contains(wavextract_method, "both")) then
      call extract_moncrief
   end if

   ! Extracción con Weyl
   if (contains(wavextract_method, "weyl") .or. contains(wavextract_method, "both")) then
      call extract_weyl
   end if


! ***************
! ***   END   ***
! ***************

end subroutine wavextract