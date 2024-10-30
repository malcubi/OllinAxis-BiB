!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/matter/sources_matter.f90,v 1.2 2019/11/13 18:59:37 malcubi Exp $

  subroutine sources_matter

! *****************************
! ***  SOURCES FOR MATTER   ***
! *****************************

! This routine calculates the sources for the matter.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  logical contains


! *****************************
! ***   REAL SCALAR FIELD   ***
! *****************************

  if (contains(mattertype,"scalar")) then
     call sources_scalar
  end if


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

  if (contains(mattertype,"complex")) then
     call sources_complex
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_matter
