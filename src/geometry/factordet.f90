!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/factordet.f90,v 1.1 2020/04/08 19:17:41 malcubi Exp $

  subroutine factordet

! **************************************************
! ***   FACTOR DETERMINANT OF CONFORMAL METRIC   ***
! **************************************************

! For Lagrangian evolutions, the determinant of the
! conformal metric must be time independent, but it
! can change due to numerical error.  Here we rescale
! the conformal metric to make sure the determinant
! remains constant in time.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  real(8) third,two


! *******************
! ***   NUMBERS   ***
! *******************

  third = 1.d0/3.d0
  two   = 2.d0


! ****************************
! ***   FIND DETERMINANT   ***
! ****************************

! Find metric determinant (divided by r**2).

  if (angmom) then
     auxarray = A*B*H - r**2*(H*C**2 + A*C2**2 + r**2*C1*(B*C1 - two*C*C2))
  else
     auxarray = (A*B - (r*C)**2)*H
  end if


! **************************
! ***   RESCALE METRIC   ***
! **************************

  auxarray = auxarray**third

  A = A/auxarray
  B = B/auxarray
  C = C/auxarray
  H = H/auxarray

  if (angmom) then
     C1 = C1/auxarray
     C2 = C2/auxarray
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine factordet

