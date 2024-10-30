!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/constraintbound.f90,v 1.4 2020/01/31 22:38:47 malcubi Exp $

  subroutine constraintbound

! ********************************************
! ***   IMPOSE CONSTRAINTS AT BOUNDARIES   ***
! ********************************************

! This routine imposes the constraints at the boundaries.
!
! Notice that the boundary conditions here are applied at the level
! of the source terms, and not directly to the dynamical quantities
! themselves.
!
! Also , this routine only applies boundary conditions to a subset
! of the geometric variables.  When we get here the radiative
! boundary conditions must have already been applied to the rest.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  real(8) zero,one,two,third
  real(8) aux1,aux2


! **********************************
! ***   DO WE WANT TO BE HERE?   ***
! **********************************

  if (rank==0) then
     print *
     print *, 'boundtype=constraint not yet implemented.'
     print *, 'Aborting! (subroutine constraintbound)'
     print *
  end if
  call die

  if (spacetime=="background") return


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  one   = 1.d0
  two   = 2.d0
  third = 1.d0/3.d0


! ************************
! ***   SANITY CHECK   ***
! ************************

! Only allow standard BSSN.

  if (eta/=2.d0) then
     if (rank==0) then
        print *
        print *, 'boundtype=constraint only works for eta=2.'
        print *, 'Aborting! (subroutine constraintbound)'
        print *
     end if
     call die
  end if


! **********************************************
! ***   BOUNDARY FOR DIAGONAL KT COMPONENTS  ***
! **********************************************

  !evolvevar => KTA
  !sourcevar => sKTA
  !Dr_var => Dr_KTA
  !Dz_var => Dz_KTA
  !call radbound_r(1.d0,0.d0)

  aux1 = 0.d0
  aux2 = 0.d0

  if (mod(rank+1,nprocr)==0) then
     sKTA(Nr,:) = sKTA(Nr,:) &
             - aux1*alpha(Nr,:)*(one/g_A(Nr,:) - third*A(Nr,:))*(RSCAL(Nr,:) - KT2(Nr,:) + two*third*trK(Nr,:)**2) &
             - aux2*(DIV_KT_r(Nr,:) + 6.d0*(KT_up_rr(Nr,:)*Dr_phi(Nr,:) + KT_up_rz(Nr,:)*Dz_phi(Nr,:)) &
             - two*third*(g_A(Nr,:)*Dr_trK(Nr,:) + r(Nr,:)*g_C(Nr,:)*Dz_trK(Nr,:)))
  end if


! **************************************************
! ***   BOUNDARY FOR OFF-DIAGONAL KT COMPONENTS  ***
! **************************************************


! ********************************
! ***   BOUNDARY FOR ALAMBDA   ***
! ********************************


! ***************
! ***   END   ***
! ***************

  end subroutine constraintbound

