!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/calc_DIVKTA.f90,v 1.5 2021/03/22 19:20:10 malcubi Exp $

  subroutine calc_DIVKTA

! *****************************
! ***   DIVERGENCE OF KTA   ***
! *****************************

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


! **********************************
! ***   FIND DIVERGENCE OF KTA   ***
! **********************************

! The divergence of the traceless extrinsic curvature
! appears as source of the Delta vector.
!
! Since we use mainly the covariant version of the extrinsic
! curvature we rewrite the divergence as
!
!      ij       im     jn .         i       m        ij  mn .                i   mn
!  D KT  = gamma  gamma   D KT  - KT   Delta  - gamma  KT   D gamma   + Delta  KT
!   j                      j  mn     m                       m     jn        mn
!
! We also define a vector DIVA_KTA = DIVKT - Delta*KT, which
! corresponds to the three first terms of the previous equation
! and is the one that truly appears on the evolution equations
! when no constraints are added.
!
! NOTE ON OPTIMAL PERFORMANCE (J.M. Torres):
!
! I know that I am basically calculating it twice and it should be
! better to calculate just the first three terms, assign them to
! DIVA_KT, then complete DIV_KT with the last one. I'm leaving it for
! testing pourposes and since the regularization seems a little more
! natural. And because even then, these are not the most time
! consuming calculations (e.g. Ricci).
!
! Side note, the irregular term on DIVA_KT_r:
!
! g_A*g_H*(KTA - g_H*KTH*A + g_H*KTH*H - KTH)/r
!
! is regularized for simplicity as:
!
! r*g_A*g_H*(Alambda - lambda*g_H*KTH)
!
! But it can be rewritten as an equivalent two parameter-free
! expressions involving KTA*(1-A*g_A) and KTH*(1-H*g_H).

  if (angmom) then

     DIV_KT_r =-g_C1**2*Dz_C1*g_C*KTC1*r**9-g_C1**2*Dr_C1*g_A*KTC1*r&
          &**8+(-4.0*g_C1**2*C1*g_C2*KTC+2.0*g_C1**2*C1*g_C*KTC2-2.0&
          &*g_C1**3*C1*KTA+3.0*g_C1**2*C2*g_C*KTC1-4.0*g_C1**2*C1*g_H&
          &*KTC1+g_C1**3*C1*KTH-g_C1**2*C*g_C2*KTC1+2.0*g_C1**2*C1&
          &*g_A*KTC1-g_C**2*Dz_C*g_C1*KTC1-g_C1*Dz_C2*g_C*g_C2*KTC1&
          &-g_C1**3*KTC1*A-g_C1*Dz_C1*g_C**2*KTC-g_C1*Dz_C1*g_C2*g_C&
          &*KTC2+3.0*g_C1**3*H*KTC1)*r**7+(-g_C*Dr_C*g_A*g_C1*KTC1&
          &-g_C1*Dr_C1*g_A*g_C*KTC-g_C1*Dr_C2*g_C2*g_A*KTC1-g_A*g_C2&
          &*g_C1*KTC2*Dr_C1)*r**6+(-2.0*g_C1*C1*g_C2**2*KTB-g_C1**2&
          &*g_C2*KTC2*A-2.0*g_C*C*g_C1**2*KTA-g_C1**2*C*g_B*KTC+g_C1&
          &**2*C2*g_B*KTC2-g_C1*Dz_C2*g_C*g_B*KTC+2.0*g_C1*C1*g_A*g_C&
          &*KTC-g_C**3*Dz_C*KTC-3.0*g_C*C*g_C1*g_H*KTC1+3.0*g_C1*C2&
          &*g_C**2*KTC-half*g_C1**2*Dz_H*g_C*KTA-g_C**2*Dz_C2*g_C2&
          &*KTC-g_C**2*Dz_C2*g_H*KTC1-g_C**2*Dz_C*g_C2*KTC2+3.0*g_C1&
          &**2*H*g_C*KTC+3.0*g_C1**2*H*g_C2*KTC2-g_C**2*Dz_B*g_C2&
          &*KTC1-3.0*g_C*C*g_C1*g_C2*KTC-4.0*g_C1*C1*g_C2*g_H*KTC2&
          &-g_C1*Dz_H*g_H*g_C*KTC1-g_C1*Dz_H*g_C2*g_C*KTC-g_C1*Dz_C1&
          &*g_C*g_A*KTA-g_C*Dz_C1*g_C2*g_A*KTC-g_C*Dz_C1*g_H*g_A*KTC1&
          &-g_A*Dz_C*g_C*g_C2*KTC1-g_C1*Dz_C2*g_C**2*KTA-g_C*Dz_C*g_B&
          &*g_C1*KTC2-g_C*Dz_C*g_C2*g_C1*KTH-g_C*Dz_C2*g_C2**2*KTC2&
          &-g_C1**2*g_C*KTC*A+g_C1*C1*g_C**2*KTB+g_C1**2*C2*g_C2*KTH&
          &-g_C1*C*g_C2**2*KTC2-g_C*C*g_A*g_C1*KTC1-g_C*Dz_A*g_A*g_C1&
          &*KTC1+3.0*g_C1*C2*g_C2*g_C*KTC2+g_C1*g_C*Dz_KTC1-g_C**2&
          &*Dz_C1*g_C2*KTB-g_C*Dz_C1*g_H*g_C1*KTH-g_C**2*Dz_C1*g_H&
          &*KTC2-half*g_C*Dz_A*g_C1**2*KTH-g_C**2*Dz_A*g_C1*KTC2)*r&
          &**5+(2.0*g_C*g_C1*Dr_KTC2+2.0*g_A*g_C1*Dr_KTC1+g_C1**2&
          &*Dr_KTH-g_C2*g_A**2*KTC1*Dr_C-g_A*g_C2*g_C1*KTH*Dr_C-g_A&
          &*g_C2*g_C*KTB*Dr_C1-g_C*Dr_B*g_C2*g_A*KTC1-g_C2*g_A**2*KTC&
          &*Dr_C1-g_H*g_A**2*KTC1*Dr_C1-g_A**2*Dr_A*g_C1*KTC1-g_C**2&
          &*Dr_C*g_A*KTC-half*g_A*Dr_A*g_C1**2*KTH-g_C1*Dr_C2*g_C*g_A&
          &*KTA-g_C1*Dr_C2*g_B*g_A*KTC-g_A*g_H*g_C*KTC2*Dr_C1-g_A*g_H&
          &*g_C1*KTH*Dr_C1-g_A*g_C*g_C2*KTC2*Dr_C-g_C1*Dr_C1*g_A**2&
          &*KTA-half*g_C1**2*Dr_H*g_A*KTA-g_C1*Dr_H*g_H*g_A*KTC1-g_C1&
          &*Dr_H*g_C2*g_A*KTC-g_C*Dr_C2*g_H*g_A*KTC1-g_C*Dr_C2*g_C2&
          &*g_A*KTC-g_A*Dr_A*g_C*g_C1*KTC2-g_A*g_B*g_C1*KTC2*Dr_C-g_A&
          &*Dr_C2*g_C2**2*KTC2-Delta_r*g_C1*KTC1)*r**4+(2.0*g_C1**2&
          &*KTA-Delta_p*g_C1*KTH+4.0*g_C1*g_H*KTC1+2.0*g_A*g_C1*KTC1&
          &+g_C*g_C2*Dz_KTC2+3.0*g_C1*g_C2*KTC+g_C1*g_C2*Dz_KTH+g_C&
          &*g_C1*KTC2+g_C1*g_B*Dz_KTC2+g_A*g_C2*Dz_KTC1-g_C*Dz_H*g_C2&
          &*g_H*KTC2-g_C1**2*KTH+g_C*C2*g_C2*g_A*KTC+g_C*C2*g_H*g_A&
          &*KTC1+3.0*g_C1*H*g_H*g_C*KTC2+3.0*g_C1*H*g_H*g_A*KTC1+3.0&
          &*g_C1*H*g_C2*g_A*KTC+3.0*g_C1*H*g_C2*g_C*KTB+3.0*g_C1**2*H&
          &*g_H*KTH+3.0*g_C1**2*H*g_A*KTA-2.0*g_C1*C1*g_H**2*KTH-g_A&
          &*C2*g_C2**2*KTC2-g_C*C*g_C2**2*KTB-g_C1**2*g_H*KTH*A+2.0&
          &*g_C*C2*g_H*g_C1*KTH+2.0*g_C**2*C2*g_H*KTC2-half*g_C*Dz_H&
          &*g_C2**2*KTB-g_A*g_C*g_C2*KTC2*C-g_A*g_C2*g_C1*KTH*C-2.0&
          &*g_C*C*g_C2*g_H*KTC2-g_C1*g_C*g_C2*KTB*A-g_C1*g_C*g_H*KTC2&
          &*A-3.0*g_A*g_C1*g_C2*KTC*A-3.0*g_A*g_C1*g_H*KTC1*A+2.0*g_C&
          &**2*C2*g_C2*KTB-g_C2*g_A**2*KTC1*C-g_C**2*Dz_C*g_B*KTB-2.0&
          &*g_A*g_C1**2*KTA*A-g_A*Dz_C*g_C**2*KTA-half*g_C*Dz_B*g_C2&
          &**2*KTH-g_C**2*Dz_B*g_B*KTC-g_C**2*C*g_A*KTC+g_C1*C1*g_A&
          &**2*KTA-g_C**2*Dz_A*g_A*KTC+g_C**2*Dz_KTC-g_C1*C*g_C2*g_B&
          &*KTB-g_C1*C*g_H*g_C2*KTH-g_C1*C*g_H*g_B*KTC2-g_A*Dz_C*g_C&
          &*g_B*KTC-g_C*Dz_B*g_B*g_C2*KTC2-half*g_C**3*Dz_A*KTB-g_C&
          &*Dz_C2*g_C2*g_B*KTB-g_C*Dz_C2*g_H*g_B*KTC2-g_A*g_B*g_C1&
          &*KTC2*C-g_C*Dz_C2*g_H*g_C2*KTH+g_C1*C2*g_B*g_C*KTB+2.0*g_C&
          &*C2*g_C1*g_A*KTA-Delta_p*g_A*KTC1-Delta_z*g_C1*KTC2-half&
          &*g_C**3*Dz_B*KTA-Delta_p*g_C*KTC2)*r**3+(-g_A*Dr_B*g_B&
          &*g_C2*KTC2-half*g_A*Dr_H*g_C2**2*KTB-half*g_A*Dr_A*g_C**2&
          &*KTB-g_A*Dr_C2*g_H*g_B*KTC2-g_A*Dr_C2*g_C2*g_B*KTB-half&
          &*g_A*Dr_B*g_C2**2*KTH-g_A*g_B*g_C*KTB*Dr_C-Delta_r*g_C*KTC&
          &-g_B*g_A**2*KTC*Dr_C-g_A**2*Dr_A*g_C*KTC-g_A*Dr_C2*g_H&
          &*g_C2*KTH+g_C**2*Dr_KTB-g_C*Dr_B*g_B*g_A*KTC-g_A*Dr_H*g_C2&
          &*g_H*KTC2+2.0*g_A*g_C*Dr_KTC-g_C*Dr_C*g_A**2*KTA-half*g_C&
          &**2*Dr_B*g_A*KTA)*r**2+(-Delta_z*g_C*KTB+2.0*g_A*g_C*KTC&
          &-half*g_C*Dz_B*g_B**2*KTB-g_A*C2*g_H*g_C2*KTH-g_A*C2*g_C2&
          &*g_B*KTB+g_C*g_B*Dz_KTB-half*g_A**2*Dz_A*g_C*KTA-g_A*g_C2&
          &*KTC2-g_A*C2*g_H*g_B*KTC2-g_A*g_B*g_C*KTB*C+g_A*g_B*Dz_KTC&
          &+2.0*g_C2*g_H*KTC2+g_C*g_H*KTC+g_C2**2*KTB-Delta_z*g_A*KTC&
          &-g_A*g_C2**2*KTB*A-g_C*C*g_H**2*KTH-2.0*g_A*g_C2*g_H*KTC2&
          &*A-g_C*C*g_A**2*KTA-g_B*g_A**2*KTC*C+g_A*g_C*Dz_KTA-half&
          &*g_C*Dz_H*g_H**2*KTH)*r-half*g_A*Dr_B*g_B**2*KTB-half*g_A&
          &*Dr_H*g_H**2*KTH-half*g_A**3*Dr_A*KTA+g_A**2*Dr_KTA&
          &-Delta_r*g_A*KTA+r*g_H*(g_A*Alambda+g_H*KTH*(r**2*C1*g_C1&
          &+C*g_C))

     DIV_KT_z = -g_C*r**9*g_C1**2*KTC1*Dr_C1-r**8*g_B*Dz_C1*g_C1**2&
          &*KTC1+(-g_C2*Dr_C2*g_C*g_C1*KTC1-g_C2*Dr_C1*g_C*g_C1*KTC2&
          &-g_C**2*g_C1*KTC1*Dr_C-g_C1*g_C**2*KTC*Dr_C1)*r**7+(-4.0&
          &*g_C2*C1*g_C1*g_H*KTC1+3.0*g_C2*H*g_C1**2*KTC1+2.0*g_B*C2&
          &*g_C1**2*KTC1-g_B*Dz_C1*g_C1*g_C*KTC-g_B*Dz_C*g_C*g_C1&
          &*KTC1-g_C2*Dz_C1*g_B*g_C1*KTC2-g_C**2*g_C1*KTC1*C-g_C2&
          &*g_C1**2*KTC1*A-g_B*Dz_C2*g_C1*g_C2*KTC1+g_C2*C1*g_C1**2&
          &*KTH+2.0*g_C2*C1*g_C*g_C1*KTC2-C*g_C1*g_C2**2*KTC1-2.0&
          &*g_C2*C1*g_C1**2*KTA-4.0*C1*g_C1*g_C2**2*KTC+2.0*g_C2*C1&
          &*g_A*g_C1*KTC1)*r**6+(-g_C**2*g_C2*KTC2*Dr_C-g_C2*Dr_H&
          &*g_C1*g_C*KTC-g_C*g_C2*g_C1*KTH*Dr_C-g_C**2*Dr_A*g_C1*KTC2&
          &-g_C**2*Dr_C2*g_C1*KTA-g_C2*Dr_C1*g_C**2*KTB-g_C**2*Dr_C2&
          &*g_H*KTC1-g_C2*Dr_C2*g_C**2*KTC-g_H*g_C**2*KTC2*Dr_C1-g_C&
          &*Dr_A*g_A*g_C1*KTC1-g_C*Dr_H*g_C1*g_H*KTC1-g_C*g_C2*g_A&
          &*KTC1*Dr_C-g_B*Dr_C*g_C*g_C1*KTC2-g_C*g_C1*g_A*KTA*Dr_C1&
          &-half*g_C*Dr_A*g_C1**2*KTH-g_C*g_H*g_A*KTC1*Dr_C1+g_C*g_C1&
          &*Dr_KTC1-g_C**3*KTC*Dr_C-g_C2*Dr_C1*g_A*g_C*KTC-g_B*Dr_C2&
          &*g_C1*g_C*KTC-g_C*g_H*g_C1*KTH*Dr_C1-half*g_C*Dr_H*g_C1**2&
          &*KTA-g_C**2*Dr_B*g_C2*KTC1-g_C2**2*Dr_C2*g_C*KTC2)*r**5+(&
          &-g_C**3*KTC*C-g_C2*C*g_C1*g_C*KTA-g_C*g_C2*g_A*KTC1*C-g_C2&
          &*Dz_C1*g_B*g_A*KTC-g_C2*Dz_C1*g_B*g_C*KTB+2.0*g_C2*C1*g_A&
          &*g_C*KTC-g_B*Dz_A*g_C*g_C1*KTC2-g_B*C*g_C*g_C1*KTC2+3.0&
          &*g_C2*H*g_C1*g_C*KTC-g_C*Dz_C*g_B*g_C2*KTC2-g_B*Dz_C*g_C2&
          &*g_A*KTC1-g_B*Dz_C*g_C2*g_C1*KTH-g_B*Dz_C2*g_C1*g_C*KTA&
          &-g_B*Dz_C2*g_H*g_C*KTC1-g_C*g_C2*g_C1*KTH*C+g_B*C2*g_C1&
          &*g_C*KTC-3.0*g_C*g_C1*g_C2*KTC*A-2.0*g_C*g_C1*g_H*KTC1*A&
          &-g_C**2*C2*g_C1*KTA-g_C**2*C2*g_H*KTC1-g_B*C*g_C1**2*KTA&
          &-g_C**2*g_C2*KTC2*C+g_C2**2*C2*g_C1*KTH+g_C2**2*C2*g_A&
          &*KTC1-g_C*g_C1**2*KTA*A-g_C2**2*Dz_C2*g_B*KTC2-4.0*C1*g_C2&
          &**2*g_H*KTC2+g_C2*C1*g_C**2*KTB-g_B*Dz_H*g_C1*g_H*KTC1-g_B&
          &*Dz_A*g_A*g_C1*KTC1-g_B*Dz_C1*g_H*g_A*KTC1-g_B*Dz_C1*g_C1&
          &*g_A*KTA-g_B*Dz_C1*g_H*g_C1*KTH-g_C2*C*g_H*g_C*KTC1-g_B&
          &*Dz_C1*g_H*g_C*KTC2+3.0*g_C2*C2*g_B*g_C1*KTC2-g_B*Dz_B*g_C&
          &*g_C2*KTC1-g_C2*Dz_C2*g_C*g_B*KTC-g_C2*Dz_H*g_C1*g_B*KTC&
          &-3.0*g_B*C*g_C1*g_C2*KTC-2.0*g_B*C*g_C1*g_H*KTC1-g_C1*g_C2&
          &**2*KTC2*A-Dz_C2*g_C1*g_B**2*KTC-C*g_C2**2*g_C*KTC-g_C**2&
          &*Dz_C*g_B*KTC+3.0*g_C2**2*H*g_C1*KTC2-Dz_C*g_B**2*g_C1&
          &*KTC2+2.0*g_C*g_C2*Dz_KTC1-2.0*C1*g_C2**3*KTB-Delta_p*g_C&
          &*KTC1-C*g_C2**3*KTC2-half*g_B*Dz_H*g_C1**2*KTA-half*g_B&
          &*Dz_A*g_C1**2*KTH)*r**4+(-g_B*Dr_B*g_C**2*KTC-g_C**2*g_A&
          &*KTA*Dr_C+g_B*g_C1*Dr_KTC2-g_C*Dr_C2*g_H*g_C2*KTH-g_C2&
          &*Dr_C2*g_B*g_C*KTB-g_B*Dr_C2*g_H*g_C*KTC2+g_C2*g_C*Dr_KTC2&
          &-half*g_C**3*Dr_A*KTB-Delta_r*g_C2*KTC1+g_C**2*Dr_KTC-g_B&
          &*Dr_B*g_C2*g_C*KTC2-half*g_C2**2*Dr_H*g_C*KTB-half*g_C&
          &*Dr_B*g_C2**2*KTH-g_C2*Dr_H*g_H*g_C*KTC2+g_C2*g_C1*Dr_KTH&
          &-g_B*Dr_C*g_C**2*KTB+g_C2*g_A*Dr_KTC1-half*g_C**3*Dr_B*KTA&
          &-g_C**2*Dr_A*g_A*KTC-g_B*Dr_C*g_A*g_C*KTC)*r**3+(-g_C2&
          &*g_C1*KTH+2.0*g_C2*g_A*KTC1+2.0*g_C*g_B*Dz_KTC+g_B*C2*g_H&
          &*g_C*KTC2+3.0*g_B*C2*g_C2*g_A*KTC-g_B*C*g_A*g_C*KTC-g_B&
          &*Dz_C*g_C*g_A*KTA+3.0*g_C2*H*g_C1*g_A*KTA+3.0*g_C2*H*g_H&
          &*g_A*KTC1+3.0*g_C2*H*g_H*g_C*KTC2+3.0*g_C2*H*g_H*g_C1*KTH&
          &+g_C2*C2*g_C*g_A*KTA-3.0*g_C*g_C2*g_H*KTC2*A-g_B*Dz_C2*g_H&
          &*g_C2*KTH-g_C*C2*g_H*g_C2*KTH-3.0*g_B*C*g_C2*g_H*KTC2&
          &-Delta_p*g_B*KTC2+g_C**2*KTC+g_C2**2*KTC+g_C**2*Dz_KTA+2.0&
          &*g_B*C2*g_H*g_C1*KTH+2.0*g_B*C2*g_C2*g_C*KTB+2.0*g_B*C2&
          &*g_H*g_A*KTC1+2.0*g_B*C2*g_C1*g_A*KTA+g_C2*g_C1*KTA-g_A&
          &*g_C2**2*KTC*A-Dz_C2*g_H*g_B**2*KTC2-half*g_C2**2*Dz_H*g_B&
          &*KTB-Dz_C*g_B**2*g_A*KTC-g_C*Dz_C*g_B**2*KTB-g_C2*Dz_H*g_H&
          &*g_B*KTC2+2.0*g_B*g_C2*Dz_KTC2+2.0*g_C2*g_H*KTC1-g_C2*g_A&
          &*g_C1*KTA*A-g_C2*g_A*g_H*KTC1*A-half*g_B*Dz_A*g_C**2*KTB&
          &+g_C2**2*Dz_KTH-Delta_p*g_C2*KTH-Delta_z*g_C2*KTC2-g_C2&
          &*g_C1*g_H*KTH*A-g_B*Dz_A*g_A*g_C*KTC-g_B*C*g_C**2*KTB-C&
          &*g_H*g_C2**2*KTH-g_C**2*g_A*KTA*C-2.0*g_B*C*g_C2**2*KTB&
          &-half*g_B*Dz_B*g_C**2*KTA-g_B**2*Dz_B*g_C*KTC-g_B**2*Dz_B&
          &*g_C2*KTC2-half*g_B*Dz_B*g_C2**2*KTH+3.0*g_C2**2*H*g_A*KTC&
          &+3.0*g_C2**2*H*g_C*KTB-2.0*g_C2*C1*g_H**2*KTH-g_C2*Dz_C2&
          &*g_B**2*KTB-2.0*g_C*g_C2**2*KTB*A+g_C2*C1*g_A**2*KTA&
          &-Delta_z*g_C*KTC)*r**2+(g_B*g_C*Dr_KTB-half*g_C*Dr_H*g_H&
          &**2*KTH-Delta_r*g_C*KTA-Delta_r*g_B*KTC+g_C*g_A*Dr_KTA+g_B&
          &*g_A*Dr_KTC-half*g_B**2*Dr_B*g_C*KTB-half*g_C*Dr_A*g_A**2&
          &*KTA)*r-g_B*C*g_H**2*KTH+g_B**2*Dz_KTB+g_B*g_H*KTC-g_C*g_H&
          &*KTH+g_C*g_H*KTA-Delta_z*g_B*KTB-half*Dz_B*g_B**3*KTB+g_B&
          &*g_A*KTC-half*g_B*Dz_H*g_H**2*KTH-g_C*g_H**2*KTH*A-half&
          &*g_B*Dz_A*g_A**2*KTA    

     DIV_KT_p = -r**8*Dz_C1*g_C1**2*g_C2*KTC1+g_H*g_B*g_A*KTC*C2-g_H*C&
          &*g_C2*g_B*KTB-g_H*g_C*g_C2*KTB*A-g_H*g_A*g_C2*KTC*A-g_H&
          &*g_A*g_C1*KTA*A+g_H*g_C*g_A*KTA*C2+3.0*g_H*g_A*g_C1*KTA*H&
          &+3.0*g_H*g_A*g_C2*KTC*H+3.0*g_H*g_C*g_C2*KTB*H+g_H*g_A**2&
          &*KTA*C1-half*g_C2*Dz_B*g_B**2*KTB+3.0*g_C*g_H**2*KTC2*H-C&
          &*g_H**2*g_B*KTC2-g_A*g_H**2*KTC1*A-g_C*g_H**2*KTC2*A-half&
          &*Dz_H*g_H**2*g_C2*KTH+3.0*g_A*g_H**2*KTC1*H+g_H*g_B*g_C&
          &*KTB*C2-2.0*g_C2*C*g_H**2*KTH+3.0*g_C1*g_H**2*KTH*H-2.0&
          &*g_C1*g_H**2*KTH*A-r**9*g_C1**3*KTC1*Dr_C1-2.0*C1*g_H**3&
          &*KTH+(-g_H*g_C1**2*KTH*Dr_C1-g_C1**2*Dr_C2*g_C*KTA-g_C1&
          &*g_C2*g_C*KTB*Dr_C1-g_C1*g_C*g_C2*KTC2*Dr_C-g_C1**2*g_A&
          &*KTA*Dr_C1-g_C2*g_C1**2*KTH*Dr_C-g_C1**2*Dr_C2*g_B*KTC&
          &-g_C2**2*g_C1*KTC2*Dr_C2-g_H*Dr_H*g_C1**2*KTC1-g_H*g_C&
          &*g_C1*KTC1*Dr_C2-g_C1*g_C**2*KTC*Dr_C-g_C2*Dr_B*g_C*g_C1&
          &*KTC1-g_C1*g_C2*g_A*KTC*Dr_C1-half*g_C1**3*Dr_H*KTA+g_C1&
          &**2*Dr_KTC1-g_B*g_C1**2*KTC2*Dr_C-g_C1*g_H*g_A*KTC1*Dr_C1&
          &-g_C1*g_H*g_C*KTC2*Dr_C1-g_C2*g_C1*g_C*KTC*Dr_C2-g_C1**2&
          &*Dr_H*g_C2*KTC-half*Dr_A*g_C1**3*KTH-Dr_A*g_C*g_C1**2*KTC2&
          &-g_C1*g_C2*g_A*KTC1*Dr_C-Dr_A*g_A*g_C1**2*KTC1)*r**5+(-2.0&
          &*g_C1**2*KTC1-4.0*g_H*C1*g_C1*g_C2*KTC-g_C2*Dz_C2*g_H*g_C&
          &*KTC1-3.0*g_C2*C*g_C1*g_H*KTC1-g_C2*Dz_C1*g_C1*g_A*KTA+2.0&
          &*g_H*g_C*g_C1*KTC2*C1+2.0*g_H*g_A*g_C1*KTC1*C1-g_H*Dz_H&
          &*g_C1*g_C2*KTC1-g_C2*C*g_C1**2*KTA-g_C2*Dz_C*g_C**2*KTC&
          &-half*g_C2*Dz_H*g_C1**2*KTA-g_C2**2*Dz_H*g_C1*KTC-g_C2**2&
          &*Dz_C1*g_A*KTC-g_C2**2*Dz_C1*g_C*KTB-g_C2*g_C1**2*KTH*C-C2&
          &*g_C1**2*g_B*KTC-half*Dz_A*g_C2*g_C1**2*KTH-g_B*g_C1**2&
          &*KTC2*C+g_C1*C2*g_C2**2*KTC2-g_C1*Dz_C*g_C2**2*KTH-g_C1&
          &*g_C**2*KTC*C-2.0*g_C1**2*g_C2*KTC*A+3.0*g_C1**2*g_H*KTC1&
          &*H-3.0*g_C1**2*g_H*KTC1*A+g_C1*g_C2*Dz_KTC1-Delta_p*g_C1&
          &*KTC1-g_C1**3*KTA*A-Dz_C2*g_C2**3*KTC2-g_C1*g_C*g_C2*KTC2&
          &*C-g_C1*Dz_A*g_C2*g_A*KTC1-g_C1*Dz_A*g_C*g_C2*KTC2-g_C1&
          &*g_C2*g_A*KTC1*C-g_C1*Dz_C*g_B*g_C2*KTC2-g_C1*Dz_C1*g_H&
          &*g_C2*KTH-g_C2*Dz_C2*g_C1*g_C*KTA-g_C2*Dz_C2*g_C1*g_B*KTC&
          &+g_C1*C2*g_C2*g_C*KTC-g_H*Dz_C1*g_C2*g_A*KTC1-g_H*Dz_C1&
          &*g_C2*g_C*KTC2-C2*g_C1**2*g_C*KTA-Dz_C2*g_C2**2*g_C*KTC&
          &-2.0*C*g_C1*g_C2**2*KTC-Dz_B*g_C*g_C2**2*KTC1-Dz_C*g_C2**2&
          &*g_A*KTC1-Dz_C*g_C2**2*g_C*KTC2-2.0*g_H*C1*g_C1**2*KTA+g_H&
          &*g_C1**2*KTH*C1-4.0*C1*g_C1*g_H**2*KTC1)*r**4+(-g_C1*Dr_C2&
          &*g_C2*g_B*KTB+g_C1*g_C*Dr_KTC-g_C1*g_B*g_C*KTB*Dr_C-g_C2&
          &*g_H*g_C1*KTH*Dr_C2-g_C2*Dr_B*g_B*g_C1*KTC2-half*Dr_B*g_C2&
          &**2*g_C1*KTH-half*g_C1*Dr_H*g_C2**2*KTB-half*g_C1*Dr_B*g_C&
          &**2*KTA-half*g_C1*Dr_A*g_C**2*KTB-g_H*g_B*g_C1*KTC2*Dr_C2&
          &+g_C2*g_C1*Dr_KTC2-g_C1*g_B*g_A*KTC*Dr_C-g_H*Dr_H*g_C2&
          &*g_C1*KTC2-g_C1*g_C*g_A*KTA*Dr_C-g_C1*Dr_A*g_A*g_C*KTC&
          &-g_C1*Dr_B*g_C*g_B*KTC)*r**3+(-g_C1*g_B*g_C*KTB*C-g_C1*g_B&
          &*g_A*KTC*C-g_C1*g_C*g_A*KTA*C+2.0*g_C2*C2*g_C1*g_A*KTA&
          &-g_C2*Dz_B*g_C*g_B*KTC-g_C2*Dz_C*g_B*g_C*KTB-g_C2*Dz_C*g_B&
          &*g_A*KTC-3.0*g_C1*g_C2*g_H*KTC2*A+3.0*g_C1*g_C2*g_H*KTC2*H&
          &-g_C1*C2*g_C2*g_B*KTB+2.0*g_C1*C2*g_H*g_C2*KTH-g_H*C*g_C1&
          &*g_B*KTC+3.0*g_H*g_C2*g_C*KTC2*C2+3.0*g_H*g_C2*g_A*KTC1*C2&
          &-g_H*C*g_C1*g_C*KTA-g_C2*Dz_C*g_C*g_A*KTA-g_C2*Dz_C2*g_H&
          &*g_B*KTC2+3.0*g_H*g_C*g_C1*KTC*H+2.0*g_H*g_A*g_C*KTC*C1&
          &-g_C2*Dz_A*g_A*g_C*KTC-g_H*C*g_C2*g_C*KTC-g_H*g_C*g_C1*KTC&
          &*A-g_C1*g_C*KTC-3.0*g_C2*g_C1*KTC2-half*g_C2**3*Dz_H*KTB&
          &+g_C1*g_C*Dz_KTA+g_C1*g_B*Dz_KTC+g_C2*g_C*Dz_KTC+g_H*g_C&
          &*Dz_KTC1-Delta_z*g_C1*KTC-Delta_p*g_C2*KTC2-half*Dz_B*g_C2&
          &**3*KTH-C*g_C2**3*KTB-half*g_C2*Dz_B*g_C**2*KTA+2.0*g_C2&
          &**2*C2*g_A*KTC+2.0*g_C2**2*C2*g_C*KTB-half*g_C2*Dz_A*g_C&
          &**2*KTB-g_C1*g_C2**2*KTB*A-Dz_C2*g_H*g_C2**2*KTH-Dz_C2&
          &*g_C2**2*g_B*KTB-Dz_B*g_B*g_C2**2*KTC2-3.0*C*g_C2**2*g_H&
          &*KTC2-2.0*g_H*C1*g_C2**2*KTB+g_H*g_C**2*KTB*C1-g_H*Dz_H&
          &*g_C2**2*KTC2+g_H*g_C**2*KTC*C2-C*g_H**2*g_C*KTC1-4.0*C1&
          &*g_C2*g_H**2*KTC2+g_C2**2*Dz_KTC2)*r**2+(g_C1*g_A*Dr_KTA&
          &+g_C2*g_A*Dr_KTC-Delta_r*g_H*KTC1+g_H*g_C1*Dr_KTH+g_H*g_A&
          &*Dr_KTC1-Delta_r*g_C1*KTA+g_H*g_C*Dr_KTC2-Delta_r*g_C2*KTC&
          &-half*Dr_H*g_H**2*g_C1*KTH-half*g_C1*Dr_B*g_B**2*KTB+g_C2&
          &*g_C*Dr_KTB-half*g_C1*Dr_A*g_A**2*KTA)*r-2.0*g_C1*g_A*KTA&
          &-g_C2*g_A*KTC-2.0*g_C2*g_C*KTB-g_H*g_C*KTC2-4.0*g_H*g_C1&
          &*KTH+2.0*g_C1*g_H*KTA+g_C2*g_B*Dz_KTB+2.0*g_C2*g_H*KTC+g_H&
          &*g_B*Dz_KTC2+g_H*g_C2*Dz_KTH-half*g_C2*Dz_A*g_A**2*KTA+(&
          &-g_C1**2*g_C*KTC*Dr_C1-g_C2*g_C1**2*KTC1*Dr_C2-g_C2*g_C1&
          &**2*KTC2*Dr_C1-g_C*g_C1**2*KTC1*Dr_C)*r**7+(-g_C1*Dz_C*g_C&
          &*g_C2*KTC1+C2*g_C1**2*g_C2*KTC1-g_C1*Dz_C1*g_C2**2*KTC2&
          &-Dz_C2*g_C1*g_C2**2*KTC1-g_C1*Dz_C1*g_C2*g_C*KTC-g_C*g_C1&
          &**2*KTC1*C)*r**6+2.0*g_H**2*KTC1-Delta_z*g_C2*KTB-Delta_z&
          &*g_H*KTC2-Delta_p*g_H*KTH

     DIVA_KT_r =-g_C1**2*Dz_C1*g_C*KTC1*r**9+(-g_C1**3*Dr_C1*KTH-2.0&
          &*g_C1**2*Dr_C1*g_C*KTC2-g_C1**3*Dr_H*KTC1-3.0*g_C1**2&
          &*Dr_C1*g_A*KTC1-2.0*g_C1**2*Dr_C2*g_C*KTC1)*r**8+(-2.0&
          &*g_C1*Dz_C1*g_C2*g_A*KTC1-g_C1**2*Dz_H*g_C2*KTC1+g_C1**3*H&
          &*KTC1-2.0*g_C1**3*C1*KTA-g_C1**2*Dz_C1*g_C2*KTH-4.0*g_C1&
          &**2*C1*g_C*KTC2-3.0*g_C1*Dz_C2*g_C*g_C2*KTC1-2.0*g_C1**3&
          &*C1*KTH-4.0*g_C1**2*C1*g_A*KTC1-g_C1**2*C2*g_C*KTC1-g_C**2&
          &*Dz_C*g_C1*KTC1-g_C1**2*C*g_C2*KTC1-g_C1**2*Dz_C1*g_B*KTC2&
          &-4.0*g_C1**2*C1*g_H*KTC1-g_C1*Dz_C1*g_C2*g_C*KTC2-g_C1**3&
          &*KTC1*A-g_C1*Dz_C1*g_C**2*KTC-4.0*g_C1**2*C1*g_C2*KTC)*r&
          &**7+(-g_C1**2*Dr_C2*g_C2*KTH-g_C*Dr_C*g_C1**2*KTH-2.0*g_C&
          &**2*Dr_C*g_C1*KTC2-g_C1*Dr_C2*g_C2*g_A*KTC1-3.0*g_C*Dr_C&
          &*g_A*g_C1*KTC1-3.0*g_C1*Dr_C1*g_A*g_C*KTC-g_C1*Dr_C1*g_C&
          &**2*KTB-g_C1**2*Dr_H*g_C*KTC-g_C**2*Dr_B*g_C1*KTC1-2.0&
          &*g_C1*Dr_C2*g_C**2*KTC-g_C1**2*Dr_H*g_C2*KTC2-g_A*g_C2&
          &*g_C1*KTC2*Dr_C1-g_C1**2*Dr_C2*g_B*KTC2-2.0*g_C1*Dr_C2&
          &*g_C2*g_C*KTC2)*r**6+(-2.0*g_C1*C1*g_C2**2*KTB-g_C1**2&
          &*g_C2*KTC2*A-2.0*g_C*C*g_C1**2*KTA-g_C1**2*C*g_B*KTC-g_C1&
          &**2*C2*g_B*KTC2-3.0*g_C1*Dz_C2*g_C*g_B*KTC-g_C1*Dz_C2*g_C2&
          &**2*KTH-4.0*g_C1*C1*g_A*g_C*KTC-g_C**3*Dz_C*KTC-3.0*g_C*C&
          &*g_C1*g_H*KTC1-g_C1*C2*g_C**2*KTC-g_C1**2*Dz_H*g_C*KTA&
          &-g_C1*Dz_H*g_C2**2*KTC2-g_C1**2*Dz_H*g_B*KTC-g_C**2*Dz_C2&
          &*g_C2*KTC-g_C**2*Dz_C2*g_H*KTC1-g_C**2*Dz_C*g_C2*KTC2-2.0&
          &*g_C**2*C*g_C1*KTC2+g_C1**2*H*g_C*KTC+g_C1**2*H*g_C2*KTC2&
          &-2.0*g_C**2*Dz_B*g_C2*KTC1-3.0*g_C*C*g_C1*g_C2*KTC-4.0&
          &*g_C1*C1*g_C2*g_H*KTC2-g_C1*Dz_H*g_H*g_C*KTC1-g_C1*Dz_C1&
          &*g_B*g_C*KTB-g_C1*Dz_H*g_C2*g_C*KTC-2.0*g_C1*Dz_C1*g_C*g_A&
          &*KTA-g_C*Dz_C1*g_C2*g_A*KTC-g_C*Dz_C1*g_H*g_A*KTC1-3.0*g_A&
          &*Dz_C*g_C*g_C2*KTC1-2.0*g_C1*Dz_C2*g_C**2*KTA-2.0*g_C1&
          &*Dz_C1*g_B*g_A*KTC-g_C*Dz_C*g_B*g_C1*KTC2-g_C*Dz_C*g_C2&
          &*g_C1*KTH-g_C*C*g_C1**2*KTH-g_C*Dz_C2*g_C2**2*KTC2-g_C1**2&
          &*g_C*KTC*A-g_A*Dz_C1*g_C2**2*KTC2-2.0*g_C1*C1*g_C**2*KTB&
          &-g_C1**2*C2*g_C2*KTH-g_C1*C*g_C2**2*KTC2-2.0*g_C1*Dz_C2&
          &*g_B*g_C2*KTC2-3.0*g_C*C*g_A*g_C1*KTC1-g_C*Dz_A*g_A*g_C1&
          &*KTC1-g_C1*C2*g_C2*g_C*KTC2+g_C1*g_C*Dz_KTC1)*r**5+(2.0&
          &*g_C*g_C1*Dr_KTC2+2.0*g_A*g_C1*Dr_KTC1+g_C1**2*Dr_KTH-g_C&
          &**3*Dr_C*KTB-g_C**2*Dr_B*g_C2*KTC2-g_C2*g_A**2*KTC1*Dr_C&
          &-g_A*g_C2*g_C1*KTH*Dr_C-g_A*g_C2*g_C*KTB*Dr_C1-g_C1*Dr_H&
          &*g_C2*g_C*KTB-g_C*Dr_B*g_C2*g_A*KTC1-g_C2*g_A**2*KTC*Dr_C1&
          &-g_H*g_A**2*KTC1*Dr_C1-2.0*g_A**2*Dr_A*g_C1*KTC1-3.0*g_C&
          &**2*Dr_C*g_A*KTC-g_A*Dr_A*g_C1**2*KTH-g_C1**2*Dr_H*g_H*KTH&
          &-g_C**2*Dr_C2*g_C2*KTB-g_C**3*Dr_B*KTC-2.0*g_C1*Dr_C2*g_C&
          &*g_A*KTA-g_C1*Dr_C2*g_B*g_A*KTC-g_C1*Dr_C2*g_B*g_C*KTB-g_A&
          &*g_H*g_C*KTC2*Dr_C1-g_A*g_H*g_C1*KTH*Dr_C1-g_A*g_C*g_C2&
          &*KTC2*Dr_C-2.0*g_C1*Dr_C1*g_A**2*KTA-g_C1**2*Dr_H*g_A*KTA&
          &-g_C**2*Dr_C2*g_H*KTC2-g_C1*Dr_H*g_H*g_A*KTC1-g_C*Dr_B*g_B&
          &*g_C1*KTC2-g_C*Dr_B*g_C2*g_C1*KTH-g_C1*Dr_H*g_C2*g_A*KTC&
          &-g_C*Dr_C2*g_H*g_A*KTC1-g_C*Dr_C2*g_C2*g_A*KTC-g_C*Dr_C2&
          &*g_H*g_C1*KTH-g_C1*Dr_H*g_H*g_C*KTC2-2.0*g_A*Dr_A*g_C*g_C1&
          &*KTC2-g_A*g_B*g_C1*KTC2*Dr_C-Delta_r*g_C1*KTC1)*r**4+(g_C1&
          &**2*KTA-Delta_p*g_C1*KTH+2.0*g_C1*g_H*KTC1+2.0*g_A*g_C1&
          &*KTC1+g_C*g_C2*Dz_KTC2+g_C1*g_C2*KTC+g_C1*g_C2*Dz_KTH+g_C&
          &*g_C1*KTC2+g_C1*g_B*Dz_KTC2+g_A*g_C2*Dz_KTC1-g_C**3*C*KTB&
          &-g_C1**2*KTH+g_C*C2*g_C2*g_A*KTC+g_C*C2*g_H*g_A*KTC1-g_A&
          &*Dz_C1*g_C2*g_B*KTB+g_C1*H*g_H*g_C*KTC2+3.0*g_C1*H*g_H*g_A&
          &*KTC1+3.0*g_C1*H*g_C2*g_A*KTC+g_C1*H*g_C2*g_C*KTB-g_A*Dz_A&
          &*g_C2*g_C1*KTH+g_C1**2*H*g_H*KTH+2.0*g_C1**2*H*g_A*KTA-2.0&
          &*g_C1*C1*g_H**2*KTH+g_A*C2*g_C2**2*KTC2-g_C*C*g_C2**2*KTB&
          &-g_C1**2*g_H*KTH*A-g_A*g_C*g_C2*KTC2*C-g_A*g_C2*g_C1*KTH*C&
          &-2.0*g_A*Dz_C*g_B*g_C2*KTC2-2.0*g_C*C*g_C2*g_H*KTC2-g_A&
          &*Dz_C1*g_H*g_B*KTC2-g_A*Dz_C1*g_H*g_C2*KTH-g_C1*g_C*g_C2&
          &*KTB*A-g_C1*g_C*g_H*KTC2*A-3.0*g_A*g_C1*g_C2*KTC*A-3.0*g_A&
          &*g_C1*g_H*KTC1*A-g_C2*g_A**2*KTC1*C-g_C**2*Dz_C*g_B*KTB&
          &-g_A**2*Dz_A*g_C2*KTC1-2.0*g_A*g_C1**2*KTA*A-g_A*Dz_C*g_C2&
          &**2*KTH-2.0*g_A*Dz_C*g_C**2*KTA-g_C*Dz_B*g_C2**2*KTH-g_C1&
          &*Dz_C2*g_B**2*KTB-2.0*g_C**2*Dz_B*g_B*KTC-3.0*g_C**2*C*g_A&
          &*KTC-2.0*g_C1*C1*g_A**2*KTA-g_C**2*Dz_A*g_A*KTC+g_C**2&
          &*Dz_KTC-g_C1*C*g_C2*g_B*KTB-g_C1*C*g_H*g_C2*KTH-g_A*Dz_A&
          &*g_B*g_C1*KTC2-g_C1*C*g_H*g_B*KTC2-g_C1*Dz_H*g_H*g_B*KTC2&
          &-g_C1*Dz_H*g_H*g_C2*KTH-g_A*Dz_A*g_C2*g_C*KTC2-3.0*g_A&
          &*Dz_C*g_C*g_B*KTC-2.0*g_C*Dz_B*g_B*g_C2*KTC2-g_C1*Dz_H&
          &*g_C2*g_B*KTB-g_C*Dz_C2*g_C2*g_B*KTB-g_C*Dz_C2*g_H*g_B&
          &*KTC2-g_A*g_B*g_C1*KTC2*C-g_C*Dz_C2*g_H*g_C2*KTH-g_C1*C2&
          &*g_B*g_C*KTB-Delta_p*g_A*KTC1-Delta_z*g_C1*KTC2-g_C**3&
          &*Dz_B*KTA-Delta_p*g_C*KTC2)*r**3+(-g_C**2*Dr_B*g_B*KTB-2.0&
          &*g_C*Dr_C*g_A**2*KTA-g_A*Dr_A*g_C**2*KTB-2.0*g_A**2*Dr_A&
          &*g_C*KTC-g_B*g_A**2*KTC*Dr_C+2.0*g_A*g_C*Dr_KTC-g_C**2&
          &*Dr_B*g_A*KTA-g_C*Dr_B*g_B*g_A*KTC+g_C**2*Dr_KTB-Delta_r&
          &*g_C*KTC-g_A*g_B*g_C*KTB*Dr_C)*r**2+(-Delta_z*g_C*KTB-g_A&
          &*g_C2*KTC2-Delta_z*g_A*KTC+g_A*C2*g_C2*g_B*KTB+g_A*H*g_C2&
          &**2*KTB-g_A*g_B*g_C*KTB*C+g_A*g_B*Dz_KTC-g_A**2*Dz_A*g_B&
          &*KTC+2.0*g_A*g_C*KTC-g_C*C*g_H**2*KTH-g_A*Dz_A*g_B*g_C*KTB&
          &+g_A*g_C*Dz_KTA-g_A**2*Dz_A*g_C*KTA-2.0*g_C*C*g_A**2*KTA&
          &-2.0*g_A*g_C2*g_H*KTC2*A-g_C*Dz_B*g_B**2*KTB-g_A*g_C2**2&
          &*KTB*A+g_C*g_H*KTC-g_B*g_A**2*KTC*C+g_C*g_B*Dz_KTB-g_A&
          &*Dz_C*g_B**2*KTB+2.0*g_A*H*g_C2*g_H*KTC2+g_A*C2*g_H*g_C2&
          &*KTH+g_A*C2*g_H*g_B*KTC2)*r+g_A**2*Dr_KTA-Delta_r*g_A*KTA&
          &-g_A**3*Dr_A*KTA+r*g_A*g_H*(Alambda-lambda*g_H*KTH)

     DIVA_KT_z =-g_C*r**9*g_C1**2*KTC1*Dr_C1-2.0*g_C*r**8*Dz_C1*g_C1&
          &*g_C2*KTC1+(-g_C1*g_C**2*KTC*Dr_C1-2.0*g_C2*Dr_C1*g_A*g_C1&
          &*KTC1-g_C**2*g_C1*KTC1*Dr_C-g_C2*Dr_C1*g_C1**2*KTH-g_C2&
          &*Dr_H*g_C1**2*KTC1-g_C2*Dr_C2*g_C*g_C1*KTC1-3.0*g_C2*Dr_C1&
          &*g_C*g_C1*KTC2-g_B*Dr_C2*g_C1**2*KTC1)*r**7+(-4.0*g_C2*C1&
          &*g_A*g_C1*KTC1-g_C2*g_C1**2*KTC1*A-2.0*g_C2*Dz_C1*g_C**2&
          &*KTC-g_C**2*Dz_C1*g_C1*KTA-2.0*g_C2*C1*g_C1**2*KTH-4.0&
          &*g_C2*C1*g_C1*g_H*KTC1-g_C2**2*Dz_H*g_C1*KTC1-2.0*g_C**2&
          &*Dz_C*g_C2*KTC1-2.0*g_C2*C1*g_C1**2*KTA-4.0*C1*g_C1*g_C2&
          &**2*KTC-g_C2**2*Dz_C1*g_C1*KTH+g_C2*H*g_C1**2*KTC1-g_B&
          &*Dz_C*g_C*g_C1*KTC1-g_C**2*Dz_C1*g_H*KTC1-4.0*g_C2*C1*g_C&
          &*g_C1*KTC2-g_B*Dz_C2*g_C1*g_C2*KTC1-2.0*g_C2**2*Dz_C1*g_C&
          &*KTC2-g_C**2*Dz_A*g_C1*KTC1-g_C2*Dz_C1*g_B*g_C1*KTC2-g_C2&
          &**2*Dz_C1*g_A*KTC1-2.0*g_C2**2*Dz_C2*g_C*KTC1-C*g_C1*g_C2&
          &**2*KTC1-g_B*Dz_C1*g_C1*g_C*KTC-g_C**2*g_C1*KTC1*C)*r**6+(&
          &-2.0*g_C2*Dr_C2*g_B*g_C1*KTC2-g_C*g_C2*g_A*KTC1*Dr_C-g_C&
          &*g_C2*g_C1*KTH*Dr_C-g_C*g_C1*g_A*KTA*Dr_C1-3.0*g_C2*Dr_C1&
          &*g_A*g_C*KTC-g_C*Dr_A*g_C1**2*KTH-2.0*g_B*Dr_C*g_A*g_C1&
          &*KTC1-g_B*Dr_C2*g_C1*g_C*KTC-3.0*g_B*Dr_C*g_C*g_C1*KTC2&
          &-g_B*Dr_B*g_C*g_C1*KTC1-g_C*g_H*g_A*KTC1*Dr_C1-g_C*g_H&
          &*g_C1*KTH*Dr_C1-g_C2*Dr_H*g_C1*g_C*KTC-2.0*g_C*Dr_A*g_A&
          &*g_C1*KTC1-2.0*g_C2*Dr_C1*g_C**2*KTB+g_C*g_C1*Dr_KTC1-g_C&
          &**3*KTC*Dr_C-g_C2**2*Dr_C2*g_A*KTC1-g_C2*Dr_C2*g_C**2*KTC&
          &-g_C2**2*Dr_H*g_C1*KTC2-g_B*Dr_C*g_C1**2*KTH-g_C**2*g_C2&
          &*KTC2*Dr_C-g_H*g_C**2*KTC2*Dr_C1-g_C2**2*Dr_C2*g_C*KTC2&
          &-g_C2**2*Dr_C2*g_C1*KTH-2.0*g_C**2*Dr_A*g_C1*KTC2)*r**5+(&
          &-g_C**3*KTC*C-g_C**3*Dz_C*KTA-g_C2**3*Dz_H*KTC2-g_C2*Dz_C1&
          &*g_C*g_A*KTA-g_C2*Dz_H*g_H*g_C*KTC1-g_C2*C*g_C1*g_C*KTA&
          &-g_C*g_C2*g_A*KTC1*C-g_C2*Dz_C1*g_B*g_A*KTC-2.0*g_C2*Dz_C1&
          &*g_B*g_C*KTB-4.0*g_C2*C1*g_A*g_C*KTC-g_B*Dz_A*g_C*g_C1&
          &*KTC2-2.0*g_B*C*g_A*g_C1*KTC1-3.0*g_B*C*g_C*g_C1*KTC2+3.0&
          &*g_C2*H*g_C1*g_C*KTC-3.0*g_C*Dz_C*g_B*g_C2*KTC2-g_B*Dz_C&
          &*g_C2*g_A*KTC1-g_B*Dz_C*g_C2*g_C1*KTH-g_B*Dz_C2*g_C1*g_C&
          &*KTA-g_B*Dz_C2*g_H*g_C*KTC1-g_C*g_C2*g_C1*KTH*C+g_B*C2&
          &*g_C1*g_C*KTC-g_C*Dz_C1*g_H*g_C2*KTH-g_C*Dz_A*g_C2*g_A&
          &*KTC1-3.0*g_C*g_C1*g_C2*KTC*A-2.0*g_C*g_C1*g_H*KTC1*A+2.0&
          &*g_C*H*g_C1*g_H*KTC1-g_B*C*g_C1**2*KTH+g_C**2*C2*g_C1*KTA&
          &+g_C**2*C2*g_H*KTC1-g_B*C*g_C1**2*KTA-g_C**2*g_C2*KTC2*C&
          &-g_C2**2*C2*g_C1*KTH-g_C**2*Dz_A*g_C2*KTC2-g_C2**2*C2*g_A&
          &*KTC1-g_C*g_C1**2*KTA*A-g_C2*Dz_C2*g_C**2*KTA-3.0*g_C2**2&
          &*Dz_C2*g_B*KTC2-4.0*C1*g_C2**2*g_H*KTC2+g_C*H*g_C1**2*KTA&
          &-2.0*g_C2*C1*g_C**2*KTB-g_C2*C*g_H*g_C*KTC1-g_B*Dz_C1*g_H&
          &*g_C*KTC2-g_C*Dz_A*g_C2*g_C1*KTH-g_C2*C2*g_B*g_C1*KTC2-2.0&
          &*g_B*Dz_B*g_C*g_C2*KTC1-3.0*g_C2*Dz_C2*g_C*g_B*KTC-g_C2&
          &*Dz_H*g_C1*g_C*KTA-g_C2*Dz_H*g_C1*g_B*KTC-3.0*g_B*C*g_C1&
          &*g_C2*KTC-2.0*g_B*C*g_C1*g_H*KTC1-g_C2**2*Dz_H*g_C*KTC&
          &-g_C1*g_C2**2*KTC2*A-Dz_C2*g_C1*g_B**2*KTC-C*g_C2**2*g_C&
          &*KTC-3.0*g_C**2*Dz_C*g_B*KTC+g_C2**2*H*g_C1*KTC2-Dz_C*g_B&
          &**2*g_C1*KTC2-g_C*Dz_C*g_C2**2*KTH+2.0*g_C*g_C2*Dz_KTC1&
          &-g_C**3*Dz_A*KTC-2.0*C1*g_C2**3*KTB-g_C2**3*Dz_C2*KTH&
          &-Delta_p*g_C*KTC1-C*g_C2**3*KTC2)*r**4+(g_C2*g_C1*Dr_KTH&
          &-g_C2**2*Dr_H*g_C*KTB-g_C2**2*Dr_H*g_A*KTC-g_B**2*Dr_B&
          &*g_C1*KTC2-2.0*g_C**2*Dr_A*g_A*KTC-g_B*Dr_C2*g_H*g_C1*KTH&
          &-g_C2*Dr_H*g_H*g_C1*KTH-g_C2*Dr_C2*g_C*g_A*KTA-g_B*Dr_B&
          &*g_C2*g_C1*KTH-g_B*Dr_B*g_C2*g_A*KTC1-g_B*Dr_C2*g_C1*g_A&
          &*KTA-2.0*g_C2*Dr_C2*g_B*g_A*KTC-g_B*Dr_B*g_C2*g_C*KTC2&
          &-g_C2*Dr_H*g_C1*g_A*KTA-g_C2*Dr_H*g_H*g_A*KTC1-g_B*Dr_B&
          &*g_C**2*KTC-2.0*g_B*Dr_C*g_C**2*KTB-g_C**2*g_A*KTA*Dr_C&
          &-g_C2*Dr_C1*g_A**2*KTA-g_B*Dr_C2*g_H*g_C*KTC2+g_B*g_C1&
          &*Dr_KTC2+g_C2*g_A*Dr_KTC1-3.0*g_B*Dr_C*g_A*g_C*KTC-g_C2&
          &*Dr_H*g_H*g_C*KTC2-g_B*Dr_C2*g_H*g_A*KTC1-2.0*g_C2*Dr_C2&
          &*g_B*g_C*KTB+g_C**2*Dr_KTC-g_C**3*Dr_A*KTB+g_C2*g_C&
          &*Dr_KTC2-Delta_r*g_C2*KTC1)*r**3+(-g_C2*g_C1*KTH+2.0*g_C2&
          &*g_A*KTC1+2.0*g_C*g_B*Dz_KTC+g_B*C2*g_H*g_C*KTC2-g_B*C2&
          &*g_C2*g_A*KTC-3.0*g_B*C*g_A*g_C*KTC-g_B*Dz_C*g_C*g_A*KTA&
          &+g_C2*H*g_C1*g_A*KTA+g_C2*H*g_H*g_A*KTC1+3.0*g_C2*H*g_H&
          &*g_C*KTC2+g_C2*H*g_H*g_C1*KTH-g_C2*C2*g_C*g_A*KTA-3.0*g_C&
          &*g_C2*g_H*KTC2*A-g_B*Dz_C2*g_H*g_C2*KTH+g_C*C2*g_H*g_C2&
          &*KTH-3.0*g_B*C*g_C2*g_H*KTC2-Delta_p*g_B*KTC2+g_C**2*KTC&
          &+g_C2**2*KTC+g_C**2*Dz_KTA+g_C2*g_C1*KTA-g_A*g_C2**2*KTC*A&
          &-Dz_C2*g_H*g_B**2*KTC2-g_C2**2*Dz_H*g_B*KTB-Dz_C*g_B**2&
          &*g_A*KTC-2.0*g_C*Dz_C*g_B**2*KTB-g_C2*Dz_H*g_H*g_B*KTC2&
          &+2.0*g_B*g_C2*Dz_KTC2+2.0*g_C2*g_H*KTC1-g_C2*g_A*g_C1*KTA&
          &*A-g_C2*g_A*g_H*KTC1*A-g_B*Dz_A*g_C**2*KTB+g_C2**2*Dz_KTH&
          &-Delta_p*g_C2*KTH-Delta_z*g_C2*KTC2-g_C2*g_C1*g_H*KTH*A&
          &-g_B*Dz_A*g_A*g_C*KTC-2.0*g_B*C*g_C**2*KTB-C*g_H*g_C2**2&
          &*KTH-g_C**2*g_A*KTA*C-2.0*g_B*C*g_C2**2*KTB-g_C**2*Dz_A&
          &*g_A*KTA-g_B*Dz_B*g_C**2*KTA-2.0*g_B**2*Dz_B*g_C*KTC-2.0&
          &*g_B**2*Dz_B*g_C2*KTC2-g_B*Dz_B*g_C2**2*KTH+g_C2**2*H*g_A&
          &*KTC+2.0*g_C2**2*H*g_C*KTB-2.0*g_C2*C1*g_H**2*KTH-2.0*g_C2&
          &*Dz_C2*g_B**2*KTB-2.0*g_C*g_C2**2*KTB*A-g_C2**2*Dz_H*g_H&
          &*KTH-2.0*g_C2*C1*g_A**2*KTA-Delta_z*g_C*KTC)*r**2+(g_B*g_C&
          &*Dr_KTB-g_B*Dr_C*g_A**2*KTA-g_B**2*Dr_B*g_C*KTB-Delta_r&
          &*g_C*KTA+g_C*g_A*Dr_KTA+g_B*g_A*Dr_KTC-Delta_r*g_B*KTC-g_C&
          &*Dr_A*g_A**2*KTA-g_B**2*Dr_B*g_A*KTC-g_B*Dr_B*g_C*g_A*KTA)&
          &*r-Dz_B*g_B**3*KTB+g_B*g_H*KTC+g_C*g_H*KTA+g_C*g_H**2*KTH&
          &*H+g_B**2*Dz_KTB-Delta_z*g_B*KTB-g_B*g_A**2*KTA*C+g_B*g_A&
          &*KTC-g_C*g_H**2*KTH*A-g_B*C*g_H**2*KTH-g_C*g_H*KTH

     DIVA_KT_p =(-g_H*Dz_C1*g_C2*g_A*KTC1-g_H*Dz_C1*g_C2*g_C*KTC2-g_H&
          &*Dz_H*g_C1*g_C2*KTC1-4.0*g_H*C1*g_C1*g_C2*KTC-3.0*g_C1&
          &*g_C2*g_A*KTC1*C-3.0*g_C1*g_C*g_C2*KTC2*C-g_C1*Dz_A*g_C2&
          &*g_A*KTC1-g_C1*Dz_A*g_C*g_C2*KTC2-2.0*g_C1*Dz_C*g_C*g_B&
          &*KTC-3.0*g_C1*Dz_C*g_B*g_C2*KTC2-g_C1*Dz_C1*g_C2*g_B*KTB&
          &-g_C2*Dz_C2*g_C1*g_C*KTA-g_C2*Dz_C2*g_C1*g_B*KTC-3.0*g_C2&
          &*Dz_C2*g_H*g_C*KTC1-3.0*g_C2*C*g_C1*g_H*KTC1-2.0*g_C1&
          &*Dz_C1*g_H*g_B*KTC2-2.0*g_C1*Dz_C1*g_H*g_C2*KTH+g_C1*C2&
          &*g_C2*g_C*KTC-4.0*g_H*g_A*g_C1*KTC1*C1-4.0*g_H*g_C*g_C1&
          &*KTC2*C1+2.0*g_C1**2*g_C2*KTC*H-2.0*g_C1**2*g_C2*KTC*A&
          &-Dz_C*g_C2**2*g_A*KTC1-Dz_C*g_C2**2*g_C*KTC2-2.0*Dz_B*g_C&
          &*g_C2**2*KTC1-Dz_C2*g_C2**2*g_C*KTC-2.0*C*g_C1*g_C2**2*KTC&
          &-g_C2*C*g_C1**2*KTA-4.0*C1*g_C1*g_H**2*KTC1-2.0*g_H*g_C1&
          &**2*KTH*C1+3.0*g_C1**2*g_H*KTC1*H-3.0*g_C1**2*g_H*KTC1*A&
          &-g_C1*g_C**2*KTC*C-g_C1*Dz_A*g_C**2*KTC-2.0*g_C1*Dz_C*g_C2&
          &**2*KTH-g_H*Dz_C1*g_C**2*KTC+g_C1*C2*g_C2**2*KTC2-2.0*g_H&
          &*C1*g_C1**2*KTA-Dz_A*g_B*g_C1**2*KTC2-Dz_A*g_C2*g_C1**2&
          &*KTH+C2*g_C1**2*g_C*KTA+C2*g_C1**2*g_B*KTC+g_C1*g_C2&
          &*Dz_KTC1-Delta_p*g_C1*KTC1+g_C1**3*KTA*H-g_C1**3*KTA*A&
          &-Dz_C2*g_C2**3*KTC2-g_C2*Dz_C*g_C**2*KTC-g_C1*Dz_C*g_C**2&
          &*KTA-g_B*g_C1**2*KTC2*C-2.0*g_C2*g_C1**2*KTH*C)*r**4-g_H*C&
          &*g_C2*g_B*KTB+g_C2*g_A*KTC+2.0*g_H*g_A*KTC1+g_H*g_C*KTC2&
          &-2.0*g_H*g_C1*KTH+(-g_C1**2*g_C*KTC*Dr_C1-g_C*g_C1**2*KTC1&
          &*Dr_C-g_C2*g_C1**2*KTC2*Dr_C1-g_C2*g_C1**2*KTC1*Dr_C2)*r&
          &**7+2.0*g_H**2*KTC1+g_C2*g_B*Dz_KTB+2.0*g_C2*g_H*KTC+g_H&
          &*g_B*Dz_KTC2+g_H*g_C2*Dz_KTH+2.0*g_C1*g_H*KTA+(-g_H*g_B&
          &*g_C1*KTC2*Dr_C2-g_H*Dr_H*g_C2*g_C1*KTC2-2.0*g_C2*g_H*g_A&
          &*KTC1*Dr_C2-2.0*g_H*g_A*g_C*KTC*Dr_C1-2.0*g_C1*Dr_A*g_A&
          &*g_C*KTC-g_C1*g_C*g_A*KTA*Dr_C-Dr_B*g_C2**2*g_C*KTC2-Dr_B&
          &*g_C2**2*g_C1*KTH-g_C2**2*g_A*KTC*Dr_C2-g_C2**2*g_C*KTB&
          &*Dr_C2-g_C2*g_C**2*KTB*Dr_C-g_C2*Dr_B*g_C**2*KTC-g_H*g_C&
          &**2*KTB*Dr_C1-g_H*g_C**2*KTC*Dr_C2-g_C1*Dr_A*g_C**2*KTB&
          &-Dr_B*g_C2**2*g_A*KTC1+g_C2*g_C1*Dr_KTC2+g_C1*g_C*Dr_KTC&
          &-g_C1*g_B*g_A*KTC*Dr_C-2.0*g_C2*g_A*g_C*KTC*Dr_C-g_C2*Dr_B&
          &*g_B*g_C1*KTC2-g_C2*g_C1*g_A*KTA*Dr_C2-2.0*g_C2*g_H*g_C&
          &*KTC2*Dr_C2-2.0*g_C2*g_H*g_C1*KTH*Dr_C2-g_C1*g_B*g_C*KTB&
          &*Dr_C-g_H*Dr_H*g_C1*g_C*KTC)*r**3+(-g_H*g_C2*g_A*KTC1*C2&
          &-g_H*g_C2*g_C*KTC2*C2-g_H*Dz_C1*g_C*g_A*KTA-g_H*Dz_C1*g_B&
          &*g_A*KTC-g_H*Dz_C1*g_B*g_C*KTB-g_H*Dz_H*g_C1*g_C*KTA-g_H&
          &*Dz_H*g_C1*g_B*KTC-g_H*Dz_H*g_C2*g_C*KTC+g_H*g_C*g_C1*KTC&
          &*H-g_H*g_C*g_C1*KTC*A-g_H*C*g_C1*g_C*KTA-g_H*C*g_C1*g_B&
          &*KTC-g_H*C*g_C2*g_C*KTC-g_C1*g_B*g_C*KTB*C-g_C1*Dz_A*g_B&
          &*g_C*KTB-g_C1*Dz_A*g_B*g_A*KTC-g_C1*Dz_A*g_C*g_A*KTA-g_C1&
          &*g_C*g_A*KTA*C-g_C2*Dz_C*g_C*g_A*KTA-2.0*g_C2*g_A*g_C*KTC&
          &*C-g_C2*Dz_C*g_B*g_A*KTC-g_C2*Dz_C*g_B*g_C*KTB-2.0*g_C2&
          &*Dz_B*g_C*g_B*KTC-3.0*g_C2*Dz_C2*g_H*g_B*KTC2-g_C1*g_B*g_A&
          &*KTC*C+g_C1*C2*g_C2*g_B*KTB+3.0*g_C1*g_C2*g_H*KTC2*H-3.0&
          &*g_C1*g_C2*g_H*KTC2*A-4.0*g_H*g_A*g_C*KTC*C1-2.0*g_H*Dz_C2&
          &*g_C*g_B*KTC-g_C2*g_C1*KTC2+g_C1*g_C2**2*KTB*H-2.0*Dz_B&
          &*g_B*g_C2**2*KTC2-Dz_C2*g_C2**2*g_B*KTB-2.0*Dz_C2*g_H*g_C2&
          &**2*KTH-3.0*C*g_C2**2*g_H*KTC2-g_C2*g_C**2*KTB*C-g_C2*Dz_B&
          &*g_C**2*KTA-Dz_H*g_H**2*g_C*KTC1-C*g_H**2*g_C*KTC1+g_C1&
          &*g_C*KTC+g_C2*g_C*Dz_KTC+g_H*g_C*Dz_KTC1+g_C1*g_C*Dz_KTA&
          &+g_C1*g_B*Dz_KTC-Delta_z*g_C1*KTC-Delta_p*g_C2*KTC2-Dz_B&
          &*g_C2**3*KTH-C*g_C2**3*KTB+g_C2**2*Dz_KTC2-4.0*C1*g_C2*g_H&
          &**2*KTC2-g_H*g_C**2*KTC*C2-g_C1*g_C2**2*KTB*A-g_C1*Dz_C&
          &*g_B**2*KTB-g_H*Dz_C2*g_C**2*KTA-g_H*Dz_H*g_C2**2*KTC2-2.0&
          &*g_H*C1*g_C2**2*KTB-2.0*g_H*g_C**2*KTB*C1)*r**2+(g_H*g_C&
          &*Dr_KTC2-g_H*Dr_H*g_C2*g_C*KTB+g_C1*g_A*Dr_KTA-g_H*g_B*g_C&
          &*KTB*Dr_C2+g_C2*g_C*Dr_KTB-g_H*g_B*g_A*KTC*Dr_C2-Delta_r&
          &*g_H*KTC1-Delta_r*g_C1*KTA+g_C2*g_A*Dr_KTC-g_H*g_A**2*KTA&
          &*Dr_C1-g_H*Dr_H*g_C1*g_A*KTA-Dr_H*g_H**2*g_A*KTC1-Delta_r&
          &*g_C2*KTC-g_C2*Dr_B*g_B*g_C*KTB+g_H*g_A*Dr_KTC1-g_H*Dr_H&
          &*g_C2*g_A*KTC-g_C1*Dr_A*g_A**2*KTA-g_C2*g_A**2*KTA*Dr_C&
          &-Dr_H*g_H**2*g_C*KTC2-g_C2*Dr_B*g_C*g_A*KTA-g_H*g_C*g_A&
          &*KTA*Dr_C2-Dr_H*g_H**2*g_C1*KTH-g_C2*Dr_B*g_B*g_A*KTC+g_H&
          &*g_C1*Dr_KTH)*r+(-3.0*g_C1*g_C*g_C2*KTC2*Dr_C+g_C1**2&
          &*Dr_KTC1-2.0*Dr_A*g_A*g_C1**2*KTC1-g_C2**2*g_C1*KTC2*Dr_C2&
          &-g_H*Dr_H*g_C1**2*KTC1-2.0*g_C2*g_C1**2*KTH*Dr_C-3.0*g_C1&
          &*g_C2*g_A*KTC1*Dr_C-Dr_A*g_C1**3*KTH-g_C2*Dr_B*g_C*g_C1&
          &*KTC1-2.0*g_H*g_C1**2*KTH*Dr_C1-g_C1*g_C2*g_C*KTB*Dr_C1&
          &-3.0*g_C1*g_H*g_A*KTC1*Dr_C1-g_C2*g_C1*g_C*KTC*Dr_C2-2.0&
          &*Dr_A*g_C*g_C1**2*KTC2-g_C1*g_C2*g_A*KTC*Dr_C1-g_C1**2*g_A&
          &*KTA*Dr_C1-g_C1*g_C**2*KTC*Dr_C-g_B*g_C1**2*KTC2*Dr_C-g_H&
          &*g_C*g_C1*KTC1*Dr_C2-3.0*g_C1*g_H*g_C*KTC2*Dr_C1)*r**5+(&
          &-Dz_C2*g_C1*g_C2**2*KTC1-Dz_A*g_C*g_C1**2*KTC1-g_C1*Dz_C1&
          &*g_C2*g_C*KTC-Dz_C1*g_C1**2*g_B*KTC-g_C1*Dz_C1*g_C2**2&
          &*KTC2-3.0*g_C1*Dz_C*g_C*g_C2*KTC1-2.0*g_C1*Dz_C1*g_H*g_C&
          &*KTC1-g_C*g_C1**2*KTC1*C+C2*g_C1**2*g_C2*KTC1-Dz_C1*g_C1&
          &**2*g_C*KTA)*r**6-r**9*g_C1**3*KTC1*Dr_C1+2.0*g_C1*g_H**2&
          &*KTH*H-2.0*g_C1*g_H**2*KTH*A-g_C2*g_A**2*KTA*C-g_C2*Dz_B&
          &*g_B**2*KTB-2.0*g_C2*C*g_H**2*KTH-Dz_H*g_H**2*g_B*KTC2&
          &-Dz_H*g_H**2*g_C2*KTH+g_A*g_H**2*KTC1*H-g_A*g_H**2*KTC1*A&
          &+g_C*g_H**2*KTC2*H-g_C*g_H**2*KTC2*A-C*g_H**2*g_B*KTC2-2.0&
          &*g_H*g_A**2*KTA*C1-g_H*Dz_C2*g_B**2*KTB+g_H*g_A*g_C2*KTC*H&
          &-r**8*Dz_C1*g_C1**2*g_C2*KTC1+g_H*g_A*g_C1*KTA*H-g_H*g_A&
          &*g_C1*KTA*A-g_H*g_C*g_A*KTA*C2-g_H*g_B*g_A*KTC*C2-g_H*g_B&
          &*g_C*KTB*C2-2.0*C1*g_H**3*KTH-Delta_z*g_C2*KTB-Delta_z*g_H&
          &*KTC2-Delta_p*g_H*KTH-g_H*Dz_H*g_C2*g_B*KTB-g_H*g_A*g_C2&
          &*KTC*A+g_H*g_C*g_C2*KTB*H-g_H*g_C*g_C2*KTB*A

  else

     DIV_KT_r =-g_C**3*Dz_C*KTC*r**5-g_A*g_C**2*KTC*Dr_C*r**4+(g_C**2&
          &*Dz_KTC-half*g_C**3*Dz_B*KTA-g_A*g_C**2*KTC*C-half*g_C**3&
          &*Dz_A*KTB-g_A*Dz_C*g_C**2*KTA-g_A*Dz_C*g_C*g_B*KTC-g_C**2&
          &*Dz_C*g_B*KTB-g_C**2*Dz_B*g_B*KTC-g_A*Dz_A*g_C**2*KTC)*r&
          &**3+(2.0*g_A*g_C*Dr_KTC-Delta_r*g_C*KTC-g_C*Dr_B*g_B*g_A&
          &*KTC-g_B*g_A**2*KTC*Dr_C-Dr_A*g_A**2*g_C*KTC-g_C*g_A**2&
          &*KTA*Dr_C+g_C**2*Dr_KTB-g_A*g_B*g_C*KTB*Dr_C-half*g_A*Dr_A&
          &*g_C**2*KTB-half*g_C**2*Dr_B*g_A*KTA)*r**2+(2.0*g_A*g_C&
          &*KTC+g_A*g_C*Dz_KTA-g_C*g_A**2*KTA*C+g_C*g_H*KTC-half*g_H&
          &**2*KTH*g_C*Dz_H-Delta_z*g_C*KTB-g_B*g_A**2*KTC*C-half*g_A&
          &**2*Dz_A*g_C*KTA+g_C*g_B*Dz_KTB+g_A*g_B*Dz_KTC-Delta_z*g_A&
          &*KTC-g_A*g_B*g_C*KTB*C-g_C*g_H**2*KTH*C-half*g_C*Dz_B*g_B&
          &**2*KTB)*r-half*g_A*Dr_B*g_B**2*KTB-half*g_H**2*KTH*g_A&
          &*Dr_H-Delta_r*g_A*KTA+g_A**2*Dr_KTA-half*g_A**3*Dr_A*KTA+r&
          &*g_H*(g_A*Alambda+g_H*KTH*C*g_C)

     DIV_KT_z = -r**5*g_C**3*KTC*Dr_C+(-g_C**3*KTC*C-Dz_C*g_C**2*g_B&
          &*KTC)*r**4+(g_C**2*Dr_KTC-g_B*g_C**2*KTB*Dr_C-g_C*g_B*g_A&
          &*KTC*Dr_C-half*g_C**3*Dr_B*KTA-g_C**2*g_A*KTA*Dr_C-half&
          &*Dr_A*g_C**3*KTB-Dr_A*g_A*g_C**2*KTC-g_B*Dr_B*g_C**2*KTC)&
          &*r**3+(-g_B*Dz_C*g_C*g_A*KTA+g_C**2*KTC-Dz_B*g_C*g_B**2&
          &*KTC-Delta_z*g_C*KTC-g_C*Dz_C*g_B**2*KTB+g_C**2*Dz_KTA+2.0&
          &*g_C*g_B*Dz_KTC-g_C**2*g_A*KTA*C-Dz_C*g_B**2*g_A*KTC-half&
          &*Dz_A*g_B*g_C**2*KTB-g_C*Dz_A*g_B*g_A*KTC-g_B*g_C**2*KTB*C&
          &-half*g_B*Dz_B*g_C**2*KTA-g_C*g_B*g_A*KTC*C)*r**2+(g_B*g_C&
          &*Dr_KTB+g_B*g_A*Dr_KTC+g_C*g_A*Dr_KTA-half*g_H**2*KTH*g_C&
          &*Dr_H-Delta_r*g_B*KTC-Delta_r*g_C*KTA-half*g_B**2*Dr_B*g_C&
          &*KTB-half*g_C*Dr_A*g_A**2*KTA)*r+g_C*g_H*KTA+g_B*g_A*KTC&
          &-half*Dz_B*g_B**3*KTB+g_B**2*Dz_KTB-Delta_z*g_B*KTB-g_C&
          &*g_H**2*KTH*A-half*g_B*Dz_A*g_A**2*KTA-half*g_H**2*KTH*g_B&
          &*Dz_H-g_C*g_H*KTH+g_B*g_H*KTC-g_B*g_H**2*KTH*C

     DIVA_KT_r =-r**5*g_C**3*Dz_C*KTC+(-g_C**3*Dr_C*KTB-g_C**3*Dr_B&
          &*KTC-3.0*g_C**2*Dr_C*g_A*KTC)*r**4+(-g_C**3*C*KTB-2.0*g_A&
          &*Dz_C*g_C**2*KTA-3.0*g_C**2*C*g_A*KTC+g_C**2*Dz_KTC-g_C**2&
          &*Dz_A*g_A*KTC-3.0*g_A*Dz_C*g_C*g_B*KTC-2.0*g_C**2*Dz_B*g_B&
          &*KTC-g_C**3*Dz_B*KTA-g_C**2*Dz_C*g_B*KTB)*r**3+(-g_A*Dr_B&
          &*g_C**2*KTA-g_A*g_B*g_C*KTB*Dr_C+2.0*g_A*g_C*Dr_KTC-2.0&
          &*g_C*Dr_C*g_A**2*KTA-Delta_r*g_C*KTC+g_C**2*Dr_KTB-g_C**2&
          &*Dr_B*g_B*KTB-2.0*g_A**2*Dr_A*g_C*KTC-g_A*Dr_A*g_C**2*KTB&
          &-g_B*g_A**2*KTC*Dr_C-g_A*Dr_B*g_C*g_B*KTC)*r**2+(-g_A*Dz_C&
          &*g_B**2*KTB-g_B*g_A**2*KTC*C-g_C*Dz_B*g_B**2*KTB+2.0*g_A&
          &*g_C*KTC+g_C*g_B*Dz_KTB+g_C*g_H*KTC+g_A*g_C*Dz_KTA+g_A*g_B&
          &*Dz_KTC-g_A*Dz_A*g_B*g_C*KTB-g_C*g_H**2*KTH*C-Delta_z*g_C&
          &*KTB-g_A**2*Dz_A*g_B*KTC-2.0*g_C*C*g_A**2*KTA-Delta_z*g_A&
          &*KTC-g_A**2*Dz_A*g_C*KTA-g_A*g_B*g_C*KTB*C)*r-Dr_A*g_A**3&
          &*KTA-Delta_r*g_A*KTA+g_A**2*Dr_KTA+r*g_A*g_H*(Alambda&
          &-lambda*g_H*KTH)

     DIVA_KT_z =-r**5*g_C**3*KTC*Dr_C+(-3.0*Dz_C*g_C**2*g_B*KTC-Dz_A*g_C&
          &**3*KTC-Dz_C*g_C**3*KTA-g_C**3*KTC*C)*r**4+(-3.0*g_C*g_B*g_A&
          &*KTC*Dr_C-2.0*Dr_A*g_A*g_C**2*KTC+g_C**2*Dr_KTC-g_B*Dr_B*g_C**2&
          &*KTC-Dr_A*g_C**3*KTB-2.0*g_B*g_C**2*KTB*Dr_C-g_C**2*g_A*KTA&
          &*Dr_C)*r**3+(-g_C*Dz_A*g_B*g_A*KTC-g_C**2*g_A*KTA*C-2.0*Dz_B&
          &*g_C*g_B**2*KTC-Dz_A*g_C**2*g_A*KTA-Dz_A*g_B*g_C**2*KTB-3.0*g_C&
          &*g_B*g_A*KTC*C-2.0*g_C*Dz_C*g_B**2*KTB-Delta_z*g_C*KTC-Dz_C*g_B&
          &**2*g_A*KTC+g_C**2*KTC+g_C**2*Dz_KTA+2.0*g_C*g_B*Dz_KTC-g_B&
          &*Dz_C*g_C*g_A*KTA-g_B*Dz_B*g_C**2*KTA-2.0*g_B*g_C**2*KTB*C)*r&
          &**2+(-Delta_r*g_B*KTC+g_B*g_A*Dr_KTC+g_C*g_A*Dr_KTA+g_B*g_C&
          &*Dr_KTB-g_B**2*Dr_B*g_C*KTB-Delta_r*g_C*KTA-g_B**2*Dr_B*g_A*KTC&
          &-g_C*Dr_A*g_A**2*KTA-g_B*Dr_B*g_C*g_A*KTA-g_B*g_A**2*KTA&
          &*Dr_C)*r-Delta_z*g_B*KTB-g_B*g_H**2*KTH*C+g_B**2*Dz_KTB-g_C*g_H**2&
          &*KTH*A-g_B*g_A**2*KTA*C+g_B*g_H*KTC+g_C*g_H*KTA+g_B*g_A&
          &*KTC-Dz_B*g_B**3*KTB

  end if


! ***************
! ***   END   ***
! ***************
  
  end subroutine calc_DIVKTA

