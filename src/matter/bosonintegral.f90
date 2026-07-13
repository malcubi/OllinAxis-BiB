
  subroutine bosonintegral(box,level)

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

  integer box,level

  real(8) integral
  real(8) smallpi

  real(8) complex_NB(0:Nb,0:Nlmax)


! *******************
! ***   NUMBERS   ***
! *******************

  smallpi = acos(-1.d0)


! ***********************************
! ***   INTEGRATE BOSON DENSITY   ***
! ***********************************

  auxarray = 2.d0*smallpi*complex_Bdens*dsqrt(hdet)*psi**6*r

! At the moment we only integrate it at t=0.

  if (t(0,0)==0.d0) then

     complex_NB(box,level) = integral(auxarray)

!    We only output the integrated boson number for the coarse grid.

     if ((box==0).and.(level==0)) then
        if (rank==0) then
           write (*,'(A,ES12.5)') ' Total boson number NB = ',complex_NB(0,0)
           print *
        end if
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine bosonintegral
