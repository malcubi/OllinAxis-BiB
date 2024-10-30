!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/weyl.f90,v 1.4 2020/10/27 23:29:01 malcubi Exp $

  subroutine weyl

! **********************************************************
! ***   ELECTRIC AND MAGNETIC COMPONENTS OF WEYL TENSOR  ***
! **********************************************************

! This routine calculates the elevtric and magnetic parts
! of the Weyl tensor, the Weyls scalars psi0 to psi4, and
! the curvature invariantes I and J.
!
! Originally written by Erik Jimenez.

! Include modules.

  use arrays
  use param
  use derivatives

! Extra variables.

  implicit none

  real(8) third,one,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  third = 1.d0/3.d0
  one = 1.d0

  smallpi = acos(-one)


! *************************
! ***   ELECTRIC PART   ***
! *************************

! Electric part of Weyl:
!                           m
! E   =  R  + trK K   - K  K  -  4 pi [ S  + 1/3 gamma  (4 rho - trS) ]
!  ij     ij       ij    im  j           ij           ij

  WEYL_E_A = RIC_A + psi4*(-KT2_A + third*trK*KTA + 2.d0/9.d0*A*trK**2)
  WEYL_E_B = RIC_B + psi4*(-KT2_B + third*trK*KTB + 2.d0/9.d0*B*trK**2)
  WEYL_E_H = RIC_H + psi4*(-KT2_H + third*trK*KTH + 2.d0/9.d0*H*trK**2)
  WEYL_E_C = RIC_C + psi4*(-KT2_C + third*trK*KTC + 2.d0/9.d0*C*trK**2)

  if (mattertype/='vacuum') then
     WEYL_E_A = WEYL_E_A - 4.d0*smallpi*(S_A + third*psi4*A*(4.d0*rho - trS))
     WEYL_E_B = WEYL_E_B - 4.d0*smallpi*(S_B + third*psi4*B*(4.d0*rho - trS))
     WEYL_E_H = WEYL_E_H - 4.d0*smallpi*(S_H + third*psi4*H*(4.d0*rho - trS))
     WEYL_E_C = WEYL_E_C - 4.d0*smallpi*(S_C + third*psi4*C*(4.d0*rho - trS))
  end if

! Non-zero angular momentum.

  if (angmom) then

     WEYL_E_C1 = RIC_C1 + psi4*(-KT2_C1 + third*trK*KTC1 + 2.d0/9.d0*C1*trK**2)
     WEYL_E_C2 = RIC_C2 + psi4*(-KT2_C2 + third*trK*KTC2 + 2.d0/9.d0*C2*trK**2)

     if (mattertype/='vacuum') then
        WEYL_E_C1 = WEYL_E_C1 - 4.d0*smallpi*(S_C1 + third*psi4*C1*(4.d0*rho - trS))
        WEYL_E_C2 = WEYL_E_C2 - 4.d0*smallpi*(S_C2 + third*psi4*C2*(4.d0*rho - trS))
     end if

  end if


! *************************
! ***   MAGNETIC PART   ***
! *************************

! Magnetic part of Weyl:
!
!         mn
! B  =  e   [ D  K   - 4 pi gamma   J  ]
!  ij    i     m  nj             jm  n
!
! where D_i is the 3D covariant derivative,
! and e_ijk is the Levi-Civita tensor.
!

! Case with non-zero angular momentum (STILL MISSING).

  if (angmom) then

!    Matter terms (STILL MISSING).

     if (mattertype/='vacuum') then

     end if

! Case with no angular momentum.  Because of the
! Levi-Civita symbol, when there is NO angular
! momentum the only non-trivial components turn
! out to be the (r,phi) and (z,phi) components.

  else

     WEYL_B_A = 0.d0
     WEYL_B_B = 0.d0
     WEYL_B_H = 0.d0
     WEYL_B_C = 0.d0

     WEYL_B_C1 = -(0.25*psi2*(A*H**2*KTB*Dz_A - A**2*H*KTB*Dz_H + &
          r**2*C*(2.*C**2*H*KTH + 2.*H*KTC*(-H*Dz_A + A*Dz_H) + C*(A*KTH*Dz_H + &
          2.*H**2*(-KTC + Dz_KTA + 2.*KTA*Dz_phi) - H*(KTA*Dz_H + 2.*A*(Dz_KTH &
          + 2.*KTH*Dz_phi)))) + r*H*(A*H*(KTC*Dr_B - 2.*KTB*Dr_C) + &
          C*(H*(KTB*Dr_A - KTA*Dr_B) + A*KTB*Dr_H)) + r**3*C*(-C**2*KTH*Dr_H + &
          2.*H**2*(-C*Dr_KTC + KTC*(Dr_C - 2.*C*Dr_phi)) + C*H*(-KTC*Dr_H + &
          2.*C*(Dr_KTH + 2.*KTH*Dr_phi))) + B*(-A**2*KTH*Dz_H + &
          A*H*(2.*A*(Dz_KTH + 2.*KTH*Dz_phi) - r*KTC*Dr_H) + H**2*(KTA*(Dz_A - &
          4.*A*Dz_phi) - r*KTC*Dr_A + 2.*A*(-Dz_KTA + r*(Dr_KTC + &
          2.*KTC*Dr_phi))) + C*(2.*H**2*KTA + r*A*KTH*Dr_H + H*(r*KTA*Dr_H - &
          2.*A*(r*Dr_KTH + KTH*(1. + 2.*r*Dr_phi)))))))/(r**2*sqrt(hdet)**3)

     WEYL_B_C2 = -(0.25*psi2*(r*(r**3*C**3*(-2.*H*Dz_KTH + KTH*(Dz_H - 4.*H*Dz_phi)) - &
          r*C*H*(A*KTB*Dz_H + H*(-KTB*Dz_A + KTA*Dz_B + 2.*r*KTC*(r*Dz_C - &
          Dr_B))) + A*H**2*(r*KTC*Dz_B - KTB*Dr_B) + r*C**2*H*(r*(r*KTC*Dz_H + &
          KTB*Dr_H) - 2.*H*(r*(-r*Dz_KTC - 2.*r*KTC*Dz_phi + Dr_KTB) + KTB*(-1. &
          + 2.*r*Dr_phi)))) + B**2*(2.*H**2*KTA + r*A*KTH*Dr_H + H*(r*KTA*Dr_H &
          - 2.*A*(r*Dr_KTH + KTH*(1. + 2.*r*Dr_phi)))) + r*B*(-r*C*(4.*H**2*KTC &
          + A*KTH*Dz_H + H*(KTA*Dz_H - 2.*A*(Dz_KTH + 2.*KTH*Dz_phi) + &
          2.*r*KTC*Dr_H)) + r*C**2*(-r*KTH*Dr_H + 2.*H*(r*Dr_KTH + KTH*(1. + &
          2.*r*Dr_phi))) + H*(r*A*KTC*Dz_H - H*(r*KTC*(Dz_A + 4.*A*Dz_phi) + &
          KTA*(-2.*r*Dz_C + Dr_B) + 2.*A*(r*Dz_KTC - Dr_KTB - &
          2.*KTB*Dr_phi))))))/(r**2*sqrt(hdet)**3)

!    Matter terms (STILL MISSING).

     if (mattertype/='vacuum') then

     end if

  end if


! ********************************
! ***   WEYL SCALARS (PSI's)   ***
! ********************************

! At the moment we only have expressions for the case with no
! angular momentum.  But I need to check them (Miguel).

  if (angmom) then

  else

!    As this subroutine is called also for gravitational wave extraction
!    we need to calculate Psi4, but the other scalars are not always necessary.

     WEYL_psi4 = (0.5*((2.*r**2*(-z*A + z*B + (-z**2 + r**2)*C)*(z*WEYL_E_A - &
          z*WEYL_E_B + (z**2 - r**2)*WEYL_E_C))/(rr**4*(A*B - &
          r**2*C**2)) + ((z**2*B + r**2*(A + 2.*z*C))*(z**2*WEYL_E_A + &
          r**2*(WEYL_E_B - 2.*z*WEYL_E_C)))/(rr**4*(A*B - &
          r**2*C**2)) + (r**2*(-z*A + z*B + (-z**2 + &
          r**2)*C)**2*(r**2*WEYL_E_A + z*(z*WEYL_E_B + &
          2.*r**2*WEYL_E_C)))/(rr**4*(A*B - r**2*C**2)*(z**2*B + r**2*(A &
          + 2.*z*C))) - (WEYL_E_H)/H + (2.*r**2*(r**2*WEYL_B_C1 + &
          z*WEYL_B_C2)*(z*A - z*B + (z**2 - r**2)*C))/(rr**2*Sqrt(r**2*A &
          + z*(z*B + 2.*r**2*C))*Sqrt(hdet)) + (2.*r**2*(-z*WEYL_B_C1 + &
          WEYL_B_C2)*Sqrt(r**2*A + z*(z*B + &
          2.*r**2*C)))/(rr**2*Sqrt(hdet))))/psi4

     if (curvInv) then

        WEYL_psi0 = (0.5*((2.*r**2*(-z*A + z*B + (-z**2 + r**2)*C)*(z*WEYL_E_A - &
             z*WEYL_E_B + (z - r)*(z + r)*WEYL_E_C))/(rr**4*(A*B - &
             r**2*C**2)) + ((z**2*B + r**2*(A + 2.*z*C))*(z**2*WEYL_E_A + &
             r**2*(WEYL_E_B - 2.*z*WEYL_E_C)))/(rr**4*(A*B - &
             r**2*C**2)) - (r**2*(-z*A + z*B + (-z + r)*(z + &
             r)*C)**2*(z**2*WEYL_E_B + r**2*(WEYL_E_A + &
             2.*z*WEYL_E_C)))/(rr**4*(-A*B + r**2*C**2)*(z**2*B + r**2*(A + &
             2.*z*C))) - (WEYL_E_H)/H + (2.*r**2*(r**2*WEYL_B_C1 + &
             z*WEYL_B_C2)*(-z*A + z*B + (-z**2 + &
             r**2)*C))/(rr**2*Sqrt(z**2*B + r**2*(A + 2.*z*C))*Sqrt(hdet)) + &
             (2.*r**2*(z*WEYL_B_C1 - WEYL_B_C2)*Sqrt(z**2*B + &
             r**2*(A + 2.*z*C)))/(rr**2*Sqrt(hdet))))/psi4

        WEYL_psi1 = (0.5*r*((-z*WEYL_E_A + z*WEYL_E_B + (-z**2 + &
             r**2)*WEYL_E_C)/(rr**2*Sqrt(A*B - r**2*C**2)) - ((-z*A + z*B + &
             (-z**2 + r**2)*C)*(r**2*WEYL_E_A + z*(z*WEYL_E_B + &
             2.*r**2*WEYL_E_C)))/(rr**2*Sqrt(A*B - r**2*C**2)*(z**2*B + &
             r**2*(A + 2.*z*C))) - ((r**2*WEYL_B_C1 + &
             z*WEYL_B_C2))/Sqrt((r**2*A + z*(z*B + 2.*r**2*C))*H)))/psi4

        WEYL_psi2 = (0.5*(r**2*WEYL_E_A + z*(z*WEYL_E_B + &
             2.*r**2*WEYL_E_C)))/(psi4*(z**2*B + r**2*(A + 2.*z*C)))

        WEYL_psi3 = (0.5*r*((z*WEYL_E_A - z*WEYL_E_B + (z**2 - &
             r**2)*WEYL_E_C)/(rr**2*Sqrt(A*B - r**2*C**2)) + ((-z*A + z*B + &
             (-z**2 + r**2)*C)*(r**2*WEYL_E_A + z*(z*WEYL_E_B + &
             2.*r**2*WEYL_E_C)))/(rr**2*Sqrt(A*B - r**2*C**2)*(z**2*B + &
             r**2*(A + 2.*z*C))) - ((r**2*WEYL_B_C1 + &
             z*WEYL_B_C2))/Sqrt((r**2*A + z*(z*B + 2.*r**2*C))*H)))/psi4

     end if

  end if


! ********************************
! ***   CURVATURE INVARIANTS   ***
! ********************************

  if (curvInv) then

!    I = 3 psi2^2 - 4 psi1 psi3 + psi0 psi4.

     curvI = 3.d0*WEYL_psi2**2 - 4.d0*WEYL_psi1*WEYL_psi3 + WEYL_psi0*WEYL_psi4

!    J = psi0 psi2 psi4 + 2 psi1 psi2 psi3 - psi0 psi3^2 - psi1^2 psi4 - psi2^3.

     curvJ = WEYL_psi0*WEYL_psi2*WEYL_psi4 + 2.d0*WEYL_psi1*WEYL_psi2*WEYL_psi3 &
           - WEYL_psi0*WEYL_psi3**2 - WEYL_psi1**2*WEYL_psi4 - WEYL_psi2**3

!    Invariant S := 27 J^2/I^3.

     curvS = 27.d0*curvJ**2/curvI**3

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine weyl

