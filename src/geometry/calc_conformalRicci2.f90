!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/calc_conformalRicci2.f90,v 1.10 2021/03/23 17:44:57 malcubi Exp $

  subroutine calc_conformalRicci2

! ****************************************************************
! ***   CALCULATES SECOND TERM OF THE CONFORMAL RICCI TENSOR   ***
! ****************************************************************

! Originally written by Jose Manuel Torres.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  real(8) half,third,one,two


! *******************
! ***   NUMBERS   ***
! *******************

  one = 1.d0
  two = 2.d0

  half  = 0.5d0
  third = 1.d0/3.d0


! **************************************************
! ***   SECOND TERM FOR RICCI CONFORMAL TENSOR   ***
! **************************************************

! The expressions are rather long, we calculate the components using:
!
! _           _mn. . _    _           m          mn                  mn
! R  = - half*g  D D g  + g   D  Delta  + 2 Delta    Delta    + Delta   Delta
!  ab             m n ab   m(a b)                  (a     b)mn         a     mnb
!
! Here we calculate the second term.

  if (angmom) then

     RIC_A = RIC_A + C*Dr_Delta_z*r + A*Dr_Delta_r + C1*Dr_Delta_p*r**3 &
           + half*(Delta_z*Dz_A + Delta_r*Dr_A)

     RIC_B = RIC_B + C*Dz_Delta_r*r + B*Dz_Delta_z + Dz_Delta_p*r**2*C2 &
           + half*(Delta_r*Dr_B + Delta_z*Dz_B)

     RIC_H = RIC_H + half*(Delta_r*Dr_H + Delta_z*Dz_H) + Delta_r*H/r

     RIC_C = RIC_C + half*C*(Dz_Delta_z + Dr_Delta_r) &
           + half*(Delta_r*Dr_C + Delta_z*Dz_C) &
           + half*(C*Delta_r + A*Dz_Delta_r + B*Dr_Delta_z)/r &
           + half*(C1*Dz_Delta_p*r**2 + C2*Dr_Delta_p*r)

     RIC_C1 = RIC_C1 + half*(Delta_r*Dr_C1 + Delta_z*Dz_C1 + C1*Dr_Delta_r) &
           + half*(3.0*Delta_r*C1 + H*Dr_Delta_p + Dr_Delta_z*C2)/r

     RIC_C2 = RIC_C2 + half*(Delta_r*Dr_C2 + C2*Dz_Delta_z + C1*Dz_Delta_r*r) &
           + half*(H*Dz_Delta_p + Delta_z*Dz_C2) + Delta_r*C2/r

     RIC_lambda = RIC_lambda + (Dr_Delta_z*C + Delta_r*lambda)/r + Dr_Delta_p*C1*r &
           + half*(Delta_r*Dr_lambda + Delta_z*Dz_lambda)  &
           +       ft6*(A*DD_Delta_rr + lambda*Delta_r)/r &
           + (1.0-ft6)*(H*DD_Delta_rr/r + lambda*Dr_Delta_r)

  else

     RIC_A = RIC_A + C*Dr_Delta_z*r + A*Dr_Delta_r &
           + half*(Delta_r*Dr_A + Delta_z*Dz_A)

     RIC_B = RIC_B + C*Dz_Delta_r*r + B*Dz_Delta_z &
           + half*(Delta_r*Dr_B + Delta_z*Dz_B)

     RIC_H = RIC_H + half*(Delta_r*Dr_H + Delta_z*Dz_H) + Delta_r*H/r

     RIC_C = RIC_C + half*C*(Dr_Delta_r + Dz_Delta_z) &
           + half*(Delta_r*Dr_C + Delta_z*Dz_C) &
           + half*(C*Delta_r + A*Dz_Delta_r + B*Dr_Delta_z)/r

     RIC_lambda = RIC_lambda + (Dr_Delta_z*C + Delta_r*lambda)/r &
           + half*(Delta_r*Dr_lambda + Delta_z*Dz_lambda)  &
           +       ft6*(A*DD_Delta_rr + lambda*Delta_r)/r &
           + (1.0-ft6)*(H*DD_Delta_rr/r + lambda*Dr_Delta_r)

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine calc_conformalRicci2
