!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/calc_conformalRicci1.f90,v 1.8 2021/03/16 20:46:08 malcubi Exp $

  subroutine calc_conformalRicci1

! ***************************************************************
! ***   CALCULATES FIRST TERM OF THE CONFORMAL RICCI TENSOR   ***
! ***************************************************************

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


! *************************************************
! ***   FIRST TERM FOR RICCI CONFORMAL TENSOR   ***
! *************************************************

! The expressions are rather long, we calculate the components using:
!
! _           _mn. . _    _           m          mn                  mn
! R  = - half*g  D D g  + g   D  Delta  + 2 Delta    Delta    + Delta   Delta
!  ab             m n ab   m(a b)                  (a     b)mn         a     mnb
!
! Here we calculate the first term.

  if (angmom) then

     RIC_A = RIC_A + 2.d0*g_C1*Dr_C1*r**3 + (2.0*g_C2*Dz_C1 + 2.0*g_C1*C1)*r**2 &
           - g_C*Drz_A*r - half*g_A*Drr_A - half*g_B*Dzz_A &
           - half*g_H*Dr_A/r + lambda*g_H

     RIC_B = RIC_B - g_C*Drz_B*r-half*g_A*Drr_B - half*g_B*Dzz_B &
           - half*g_H*Dr_B/r

     RIC_H = RIC_H - 2.d0*g_C1*Dr_C1*r**3 + (-2.0*g_C2*Dz_C1 - 2.0*g_C1*C1)*r**2 &
           - g_C*Drz_H*r - half*g_A*Drr_H - half*g_B*Dzz_H &
           - half*g_H*Dr_H/r - lambda*g_H

     RIC_C = RIC_C + (g_C1*Dr_C2-g_C*Drz_C)*r - half*g_A*Drr_C + g_C2*Dz_C2 &
           - half*g_B*Dzz_C - g_C*Dz_C + (-g_A*Dr_C - half*g_H*Dr_C)/r

     RIC_C1= RIC_C1 - g_C*Drz_C1*r - half*g_B*Dzz_C1 - 2.0*g_C*Dz_C1 &
           - half*g_A*Drr_C1 + (-2.0*g_A*Dr_C1 + g_C1*Dr_H - half*g_H*Dr_C1 &
           - g_C1*Dr_A)/r + g_C1*lambda - C1*g_lambda - g_C2*Dz_lambda

     RIC_C2= RIC_C2 + (-g_C1*Dr_C - g_C*Drz_C2)*r - half*g_A*Drr_C2 &
           - g_C*Dz_C2 - half*g_B*Dzz_C2 - g_C2*Dz_C &
           + (-g_A*Dr_C2 - half*g_H*Dr_C2)/r

     RIC_lambda = RIC_lambda + 4.0*g_C1*Dr_C1*r - half*g_B*Dzz_lambda &
           + 4.0*g_C2*Dz_C1 - half*g_A*Drr_lambda + 4.0*g_C1*C1 &
           - g_lambda*lambda - (2.0*g_A + half*g_H)*Dr_lambda/r &
           - g_C*(r*Drz_lambda + 2.d0*Dz_lambda)

  else

     RIC_A = RIC_A - g_C*Drz_A*r - half*g_B*Dzz_A - half*g_A*Drr_A &
           - half*g_H*Dr_A/r + g_H*lambda

     RIC_B = RIC_B - g_C*Drz_B*r - half*g_A*Drr_B - half*g_B*Dzz_B &
           - half*g_H*Dr_B/r

     RIC_H = RIC_H - g_C*Drz_H*r - half*g_B*Dzz_H - half*g_A*Drr_H &
           - half*g_H*Dr_H/r - g_H*lambda

     RIC_C = RIC_C - g_C*Drz_C*r - half*g_A*Drr_C - g_C*Dz_C &
           - half*g_B*Dzz_C + (-g_A*Dr_C - half*g_H*Dr_C)/r

     RIC_lambda = RIC_lambda - g_C*(r*Drz_lambda + 2.d0*Dz_lambda) &
          - half*g_B*Dzz_lambda - half*g_A*Drr_lambda &
          - (2.d0*g_A + half*g_H)*Dr_lambda/r ! -g_lambda*lambda cancels with term in ricci3

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine calc_conformalRicci1

