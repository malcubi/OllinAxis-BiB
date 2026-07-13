
  subroutine bosonintegral(complex_NB)

! **************************************************
! ***   CALCULATION OF INTEGRATED BOSON CHARGE   ***
! **************************************************

! This is the integrated boson number (boson charge):
!
!                /
! complex_NB  =  | complex_Bdens dV
!                /
!
! with dV the physical volume element.
!
! Notice that this assumes that the spacetime is REGULAR at the
! origin, so it will not work for eternal black holes such as
! Schwarzschild or Kerr.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  real(8) integral
  real(8) smallpi
  real(8) complex_NB


! *******************
! ***   NUMBERS   ***
! *******************

  smallpi = acos(-1.d0)


! ***********************************
! ***   INTEGRATE BOSON DENSITY   ***
! ***********************************

  auxarray = 2.d0*smallpi*complex_Bdens*dsqrt(hdet)*psi**6*r

  complex_NB = integral(auxarray)


! ***************
! ***   END   ***
! ***************

  end subroutine bosonintegral
