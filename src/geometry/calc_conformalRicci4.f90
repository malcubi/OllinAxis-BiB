!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/calc_conformalRicci4.f90,v 1.5 2021/03/10 18:19:50 malcubi Exp $

  subroutine calc_conformalRicci4

! ****************************************************************
! ***   CALCULATES FOURTH TERM OF THE CONFORMAL RICCI TENSOR   ***
! ****************************************************************

! Originally written by Jose Manuel Torres.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  real(8) half,third,one,two,fourth


! *******************
! ***   NUMBERS   ***
! *******************

  one = 1.d0
  two = 2.d0

  half  = 0.5d0
  third = 1.d0/3.d0
  fourth = half*half


! **************************************************
! ***   FOURTH TERM FOR RICCI CONFORMAL TENSOR   ***
! **************************************************

! The expressions are rather long, we calculate the components using:
!
! _           _mn. . _    _           m          mn                  mn
! R  = - half*g  D D g  + g   D  Delta  + 2 Delta    Delta    + Delta   Delta
!  ab             m n ab   m(a b)                  (a     b)mn         a     mnb
!
! Here we calculate the fourth (last) term.

  if (angmom) then

     RIC_A = RIC_A +(-C1*g_C1**2*Dr_C1-g_C2*g_C1*Dr_C1*Dz_C1)*r**7+(&
          &-7.0*half*C1*g_C2*g_C1*Dz_C1-half*g_C2**2*Dz_C1**2-3.0*C1&
          &**2*g_C1**2+g_C2*g_C1*Dr_C1*Dr_C2)*r**6+(-C1*g_C*g_C1*Dr_C&
          &+5.0*half*C1*g_C2*g_C1*Dr_C2+g_C2*g_C*Dr_C*Dz_C1+Dz_A*g_C1&
          &*g_C*Dr_C1-g_B*g_C1*Dr_C*Dz_C1+g_C2*g_C1*Dr_C1*C2+g_H*g_C&
          &*Dr_C1*Dz_C1)*r**5+(7.0*half*C1*g_C1*g_C*Dz_A-C*g_B*g_C1&
          &*Dz_C1+Dr_A*g_C*g_C1*Dr_C2+half*Dr_A*g_C1**2*Dr_H+half&
          &*g_C2**2*Dr_C2**2+g_C2*g_C*Dr_C*Dr_C2+g_C2*g_C1*Dr_C*Dr_H&
          &+half*g_H*g_B*Dz_C1**2-C1*g_C*g_C1*C+Dz_A*g_B*g_C1*Dz_C1&
          &+7.0*half*C1*g_H*g_C*Dz_C1-Dz_A*g_C2*g_C*Dz_C1+C*g_C2*g_C&
          &*Dz_C1+g_B*g_C1*Dr_C*Dr_C2+2.0*g_C2*g_A*Dr_C1*Dr_C-half&
          &*g_C2**2*C2*Dz_C1+g_H*g_C*Dr_C1*Dr_C2+g_C2*g_C*Dr_C1*Dr_B&
          &+g_H*g_C1*Dr_C1*Dr_H+Dr_A*g_C1*g_A*Dr_C1+2.0*C1*g_C2*g_C1&
          &*C2+g_H*g_A*Dr_C1**2)*r**4+(2.0*g_C2*g_C*Dr_C*C2+C*g_C2&
          &*g_C*Dr_C2+half*Dr_A*g_C1**2*H+Dz_A*g_C**2*Dr_C+2.0*g_H&
          &*g_C*Dr_C1*C2-g_C1*Dr_C1+2.0*g_C2*g_A*Dr_C1*C+6.0*C1*g_H&
          &*g_A*Dr_C1+5.0*half*C1*g_C1*g_A*Dr_A+5.0*half*C1*g_C2*g_C&
          &*Dr_B+g_H*g_C1*Dr_C1*H-Dz_A*g_C2*g_A*Dr_C1+g_C2*g_C1*Dr_C&
          &*H+5.0*half*C1*g_H*g_C1*Dr_H+3.0*half*g_C2**2*C2*Dr_C2+C&
          &*g_B*g_C1*Dr_C2+3.0*half*Dr_A*g_C*g_C1*C2+g_B*g_C1*Dr_C*C2&
          &+6.0*C1*g_C2*g_A*Dr_C+C*g_C2*g_C1*Dr_H+5.0*half*C1*g_H*g_C&
          &*Dr_C2)*r**3+(half*g_H*g_B*Dr_C2**2+2.0*C*g_C2*g_C*C2-3.0&
          &*g_C1*C1+g_B*g_A*Dr_C**2+half*Dr_A*g_C**2*Dr_B+g_C2*g_C1*C&
          &*H-7.0*half*C1*g_C2*g_A*Dz_A+half*Dz_A*g_B*g_C1*C2-half&
          &*Dz_A*g_C2*g_C*C2+Dr_A*g_C*g_A*Dr_C+g_C2**2*C2**2+2.0*C1&
          &*g_H*g_C1*H-half*Dz_A**2*g_C**2+half*g_H*g_B*Dz_C1*C2+g_H&
          &*g_C2*Dr_H*Dr_C2+g_B*g_C*Dr_C*Dr_B+Dz_A*g_C**2*C+6.0*C1&
          &*g_C2*g_A*C+5.0*C1*g_H*g_C*C2+9.0*C1**2*g_H*g_A+g_B*g_C2&
          &*Dr_C2*Dr_B+g_B*g_C1*C*C2-half*g_C2*Dz_C1+half*g_C2**2&
          &*Dr_H*Dr_B)*r**2+(-half*g_C2*Dr_C2+g_B*g_C*C*Dr_B+3.0*half&
          &*g_H*g_C2*Dr_H*C2+3.0*half*g_B*g_C2*C2*Dr_B+Dr_A*g_C*g_A*C&
          &-Dz_A*g_B*g_A*Dr_C+2.0*g_B*g_A*Dr_C*C+half*g_C2**2*H*Dr_B&
          &+3.0*half*g_H*g_B*C2*Dr_C2+g_H*g_C2*H*Dr_C2)*r+half*Dz_A&
          &**2*g_B*g_A+g_H*g_B*C2**2+1.0/4.0*Dr_A**2*g_A**2+g_B*g_A*C&
          &**2-Dz_A*g_B*g_A*C-g_C2*C2+g_H*g_C2*H*C2+1.0/4.0*g_H**2&
          &*Dr_H**2+1.0/4.0*g_B**2*Dr_B**2+(half*g_H**2*Dr_H*H-half&
          &*g_H*Dr_H)/r
     
     RIC_B = RIC_B +half*r**8*g_C1**2*Dz_C1**2+(-half*g_C1**2*Dr_C2&
          &**2+g_C1*g_C*Dz_C*Dz_C1+g_C1*g_C2*Dz_C2*Dz_C1)*r**6+(g_C1&
          &*g_C*Dz_C*Dr_C2-g_C1*g_C2*Dz_C2*Dr_C2-2.0*C2*g_C1**2&
          &*Dr_C2)*r**5+(g_H*g_C*Dz_C2*Dz_C1+g_A*g_C1*Dz_C1*Dz_A+2.0&
          &*g_C1*g_C*Dz_C*C2+half*g_C1**2*Dz_H*Dz_A+g_C2*g_C*Dz_B&
          &*Dz_C1+2.0*g_C1*g_B*Dz_C2*Dz_C-2.0*g_C1**2*C2**2+g_H*g_C1&
          &*Dz_H*Dz_C1+half*g_H*g_A*Dz_C1**2-2.0*g_C1*g_C2*Dz_C2*C2&
          &+Dz_H*g_C1*g_C2*Dz_C-g_C*g_C1*Dr_C2*Dr_B+g_A*g_C2*Dz_C&
          &*Dz_C1+g_C1*g_C*Dz_C2*Dz_A)*r**4+(-g_A*g_C2*Dz_C*Dr_C2&
          &-Dr_B*g_C1*g_B*Dz_C2+Dr_B*g_C**2*Dz_C+g_H*g_C*Dz_C2*Dr_C2&
          &+Dr_B*g_C2*g_C*Dz_C2-2.0*C2*g_C1*g_C*Dr_B)*r**3+(Dz_H*g_H&
          &*g_C2*Dz_C2+half*g_H*g_A*Dr_C2**2+Dz_B*g_C*g_B*Dz_C+Dz_B&
          &*g_C2*g_B*Dz_C2+g_C2*g_A*Dr_C2*Dr_B+half*Dz_B*g_C2**2*Dz_H&
          &-2.0*g_A*g_C2*Dz_C*C2+g_A*g_C*Dz_C*Dz_A+Dz_C2**2*g_H*g_B&
          &+half*Dz_A*g_C**2*Dz_B+2.0*g_H*g_C*Dz_C2*C2-half*g_C**2&
          &*Dr_B**2+g_A*g_B*Dz_C**2)*r**2+(2.0*C2*g_H*g_A*Dr_C2-Dr_B&
          &*g_A*g_B*Dz_C+2.0*C2*g_C2*g_A*Dr_B)*r+half*g_B*g_A*Dr_B**2&
          &+1.0/4.0*g_B**2*Dz_B**2+1.0/4.0*g_A**2*Dz_A**2+2.0*g_H*g_A&
          &*C2**2+1.0/4.0*g_H**2*Dz_H**2

     RIC_H = RIC_H -half*g_C**2*Dz_C1**2*r**6+(g_C**2*Dz_C1*Dr_C2&
          &+g_C*g_C1*Dr_H*Dz_C1)*r**5+(-g_C*g_C1*Dr_H*Dr_C2-half*C&
          &*g_B*g_C1*Dz_C1-half*A*g_C1*g_C*Dz_C1-half*g_C**2*Dr_C2**2&
          &+3.0*half*g_C**2*Dz_C1*C2-g_C2*g_C*Dz_H*Dz_C1+g_B*g_C1&
          &*Dz_H*Dz_C1+half*C*g_C2*g_C*Dz_C1+half*g_B*g_A*Dz_C1**2&
          &-half*g_C1**2*Dr_H**2+3.0*half*g_C*g_C1*H*Dz_C1)*r**4+(&
          &-g_B*g_A*Dr_C2*Dz_C1+half*C*g_B*g_C1*Dr_C2-3.0*half*g_C&
          &*g_C1*Dr_H*C2-g_B*g_C1*Dz_H*Dr_C2-3.0*half*g_C**2*C2*Dr_C2&
          &+half*A*g_C1*g_C*Dr_C2-half*C*g_C2*g_C*Dr_C2-g_C2*g_A*Dr_H&
          &*Dz_C1-Dz_H*g_C2*g_C1*Dr_H+g_C2*g_C*Dz_H*Dr_C2-3.0*half&
          &*g_C*g_C1*H*Dr_C2-3.0*half*g_C1**2*Dr_H*H+half*A*g_C1**2&
          &*Dr_H+half*C*g_C2*g_C1*Dr_H)*r**3+(-2.0*g_C1*C1+A*g_C1**2&
          &*H+g_C2*g_C1*C*H-g_C**2*C2**2+A*g_C1*g_C*C2-half*Dz_H**2&
          &*g_C2**2-C*g_C2*g_C*C2-3.0*half*g_B*g_C1*Dz_H*C2-2.0*g_C&
          &*g_C1*H*C2+half*g_C2*Dz_C1+half*A*g_C2*g_C1*Dz_H-g_C1**2*H&
          &**2+half*C*g_C2**2*Dz_H-3.0*half*g_C2*g_A*H*Dz_C1-3.0*half&
          &*Dz_H*g_C2*g_C1*H+g_C2*g_A*Dr_H*Dr_C2+g_B*g_C1*C*C2-3.0&
          &*half*g_B*g_A*C2*Dz_C1+half*A*g_C2*g_A*Dz_C1+half*g_B*g_A&
          &*Dr_C2**2+3.0*half*g_C2*g_C*Dz_H*C2)*r**2+(-half*C*g_H*g_C&
          &*Dr_H+Dz_H*g_H*g_C*Dr_H-half*g_C2*Dr_C2+3.0*half*g_B*g_A&
          &*C2*Dr_C2+3.0*half*g_C2*g_A*Dr_H*C2+3.0*half*g_C2*g_A*H&
          &*Dr_C2-half*A*g_C2*g_A*Dr_C2)*r-half*C*g_H*g_B*Dz_H-C*g_H&
          &*g_C*H-half*g_C*Dz_H+g_B*g_A*C2**2+half*g_H*g_A*Dr_H**2&
          &-g_C2*C2-half*A*g_H*g_C*Dz_H-A*g_C2*g_A*C2+3.0*half*Dz_H&
          &*g_H*g_C*H+half*Dz_H**2*g_H*g_B+2.0*g_C2*g_A*H*C2+(-half&
          &*g_H*Dr_H-half*g_A*Dr_H-half*A*g_H*g_A*Dr_H+3.0*half*g_H&
          &*g_A*Dr_H*H)/r+g_H*lambda*(1.0-g_A*H)
     
     RIC_C = RIC_C +half*g_C1**2*Dr_C1*Dz_C1*r**7+(3.0*half*C1*g_C1&
          &**2*Dz_C1-half*g_C1**2*Dr_C1*Dr_C2)*r**6+(g_C1*g_C*Dr_C1&
          &*Dz_C-g_C1**2*Dr_C1*C2+half*g_C*g_C1*Dr_C*Dz_C1-3.0*half&
          &*C1*g_C1**2*Dr_C2)*r**5+(-half*g_C*g_C1*Dr_C2*Dr_C+half&
          &*g_C1*g_B*Dz_C*Dz_C1-half*g_C2*g_C*Dz_C1*Dz_C-half*Dr_B&
          &*g_C1*g_C*Dr_C1+3.0*C1*g_C1*g_C*Dz_C-half*g_C2**2*Dz_C2&
          &*Dz_C1-3.0*g_C1**2*C1*C2+half*g_C*g_C1*C*Dz_C1)*r**4+(half&
          &*g_C2**2*Dz_C2*Dr_C2+1.0/4.0*g_C1**2*Dr_H*Dz_A+1.0/4.0&
          &*Dz_H*g_C1**2*Dr_A+g_C**2*Dr_C*Dz_C+Dz_C2*g_C2*g_C*Dr_C&
          &+Dz_C2*g_H*g_C*Dr_C1-half*g_B*g_C1*Dz_C1*Dr_B+half*g_C1&
          &*g_A*Dr_A*Dz_C1+half*Dz_H*g_H*g_C1*Dr_C1-g_C*g_C1*C2*Dr_C&
          &+half*Dz_H*g_C2*g_C1*Dr_C+half*g_C1*g_C2*Dz_C*Dr_H-half&
          &*g_C*g_C1*Dr_C2*C+g_C*g_C1*Dr_C2*Dz_A+half*g_C1*g_A*Dr_C1&
          &*Dz_A-3.0*half*C1*g_C1*g_C*Dr_B+Dr_B*g_C2*g_C*Dz_C1+half&
          &*g_C2*g_A*Dr_C*Dz_C1+g_H*g_C*Dz_C1*Dr_C2+half*Dr_A*g_C1&
          &*g_C*Dz_C2+half*g_H*g_A*Dr_C1*Dz_C1+half*g_H*g_C1*Dr_H&
          &*Dz_C1+half*g_C1*g_B*Dz_C*Dr_C2+half*g_C*g_C2*Dz_C*Dr_C2&
          &+half*Dz_B*g_C2*g_C*Dr_C1)*r**3+(Dz_C2*g_C2*g_C*C+g_C1*g_B&
          &*Dz_C*C2+half*Dz_H*g_C2*g_C1*C+3.0*half*C1*g_C1*g_A*Dz_A&
          &+3.0*half*C1*g_H*g_A*Dz_C1+3.0*half*C1*g_C2*g_C*Dz_B+half&
          &*g_C2*g_A*Dr_C*Dr_C2+3.0*C1*g_H*g_C*Dz_C2-half*Dz_A*g_C2&
          &*g_C*Dz_C2+C1*g_H*g_C1*Dz_H+C2*g_C*g_C2*Dz_C+2.0*C2*g_C1&
          &*g_C*Dz_A+2.0*C2*g_H*g_C*Dz_C1+half*g_C2*g_A*C*Dz_C1+half&
          &*g_H*g_A*Dr_C1*Dr_C2+H*g_C1*g_C2*Dz_C+H*g_C1*g_H*Dz_C1&
          &+half*g_H*g_B*Dz_C2*Dz_C1+half*Dz_A*g_C1*g_B*Dz_C2+half&
          &*Dr_B*g_C2*g_A*Dr_C1-half*g_C1*Dz_C1-half*Dz_A*g_C**2*Dz_C&
          &+half*H*g_C1**2*Dz_A+C2*g_C2**2*Dz_C2-half*g_C**2*Dr_B&
          &*Dr_C+g_C**2*C*Dz_C-g_C*g_C1*C2*C)*r**2+(3.0*half*C1*g_H&
          &*g_A*Dr_C2+half*Dz_H*g_C2*g_H*Dr_C2+g_H*g_A*Dr_C1*C2+3.0&
          &/4.0*g_C**2*Dr_B*Dz_A+1.0/4.0*Dz_B*g_C2**2*Dr_H+half*Dz_B&
          &*g_B*g_C*Dr_C-half*g_C1*Dr_C2+1.0/4.0*Dz_H*g_C2**2*Dr_B&
          &+half*g_C*g_A*Dr_C*Dz_A+half*Dr_A*g_A*g_C*Dz_C+1.0/4.0&
          &*Dz_B*g_C**2*Dr_A+g_C2*g_A*Dr_C*C2+half*g_C2*g_B*Dz_C2&
          &*Dr_B+half*g_H*g_C2*Dz_C2*Dr_H-half*g_C**2*Dr_B*C+half&
          &*g_C2*g_A*C*Dr_C2+3.0*half*C1*g_C2*g_A*Dr_B+half*g_C*g_B&
          &*Dz_C*Dr_B-half*g_C2*g_A*Dz_A*Dr_C2+half*g_H*g_B*Dz_C2&
          &*Dr_C2+half*Dz_B*g_B*g_C2*Dr_C2)*r+half*C2*g_H*g_C2*Dz_H&
          &+g_H*g_B*Dz_C2*C2+half*Dz_B*g_C2**2*H-C2*g_C2*g_A*Dz_A+3.0&
          &*g_H*g_A*C1*C2+Dz_B*g_B*g_C2*C2+g_C2*g_A*C*C2-g_C1*C2+half&
          &*g_C*g_A*C*Dz_A+half*Dr_B*g_B*g_A*Dr_C+H*g_H*g_C2*Dz_C2&
          &+half*Dz_B*g_B*g_C*C+half*Dz_A*g_A*g_B*Dz_C-g_C2*Dz_C2&
          &+(half*Dr_B*g_B*g_A*C+1.0/4.0*g_A**2*Dr_A*Dz_A-half*g_B&
          &*g_A*Dz_A*Dr_B+1.0/4.0*Dz_B*g_B**2*Dr_B+1.0/4.0*Dz_H*g_H&
          &**2*Dr_H)/r

     RIC_C1 = RIC_C1+A*g_C*g_C1*C-half*A*g_C1*g_C*Dz_A-half*A*g_H*g_C&
          &*Dz_C1+A*g_C2*g_C1*C2+half*C*g_C2*g_C*Dz_A-half*C*g_C1*g_B&
          &*Dz_A-half*C*g_H*g_B*Dz_C1+half*H*g_H*g_C*Dz_C1-g_C*g_C1*H&
          &*C-half*Dz_A*g_C2*g_C*Dz_H+half*Dz_A*g_B*g_A*Dz_C1+half&
          &*Dz_A*g_B*g_C1*Dz_H+half*H*g_C1*g_C*Dz_A+half*g_H*g_B*Dz_H&
          &*Dz_C1+(-half*g_C2*g_C*Dz_C1**2+half*g_B*g_C1*Dz_C1**2&
          &-half*g_C1**2*Dr_H*Dr_C1+3.0*half*C1*g_C1*g_C*Dz_C1-half&
          &*g_C1*g_C*Dr_C2*Dr_C1)*r**4+(C*g_C2*g_C1*Dr_C1-half*g_C2&
          &*g_A*Dz_C1*Dr_C1-half*g_B*g_C1*Dz_C1*Dr_C2-g_C1*g_C*C2&
          &*Dr_C1+half*g_C2*g_C*Dr_C2*Dz_C1+A*g_C1**2*Dr_C1-5.0*half&
          &*g_C*g_C1*Dr_C2*C1-2.0*g_C1**2*Dr_H*C1-half*g_C2*g_C1*Dz_H&
          &*Dr_C1+half*C1*g_C1**2*Dr_A-g_C1*g_A*C1*Dr_C1+half*g_C**2&
          &*Dz_C1*Dr_C+C1*g_H*g_C1*Dr_C1+C1*g_C2*g_C1*Dr_C-half*g_C2&
          &*g_C1*Dr_H*Dz_C1-g_C1**2*H*Dr_C1)*r**3+(-3.0*half*C1*g_C2&
          &*g_C1*Dz_H+half*g_C2*g_C*Dz_C1*C2+half*g_C**2*Dz_C1*C-3.0&
          &*g_C1*g_A*C1**2-3.0*half*C1*g_C2*g_A*Dz_C1+3.0*A*g_C1**2&
          &*C1-5.0*g_C*g_C1*C2*C1+4.0*C1*g_C2*g_C1*C+half*g_C2*g_A&
          &*Dr_C2*Dr_C1+half*C*g_C2**2*Dz_C1-half*Dz_A*g_C**2*Dz_C1&
          &-half*g_B*g_C1*C2*Dz_C1-half*g_C*g_C1*Dr_H*Dr_C+3.0*C1**2&
          &*g_H*g_C1-4.0*H*g_C1**2*C1-half*H*g_C2*g_C1*Dz_C1+half*A&
          &*g_C2*g_C1*Dz_C1-half*g_C2**2*Dz_H*Dz_C1-half*g_C**2*Dr_C2&
          &*Dr_C)*r**2+(A*g_C*g_C1*Dr_C+half*A*g_C2*g_C1*Dr_C2+half*C&
          &*g_C1*g_C*Dr_A+C*g_B*g_C1*Dr_C-g_C*g_C1*H*Dr_C+half*Dz_A&
          &*g_C*g_C1*Dr_H+half*g_H*g_C*Dr_H*Dz_C1+half*Dz_A*g_C**2&
          &*Dr_C2+half*C*g_C2**2*Dr_C2+half*C1*g_C2**2*Dr_B-half*g_C&
          &**2*Dr_B*C1-half*g_C**2*Dr_C2*C-g_C**2*C2*Dr_C+half*g_C2&
          &*g_C*Dz_H*Dr_C-half*g_B*g_A*Dz_C1*Dr_C-half*g_C*g_C1*Dr_H&
          &*C-half*H*g_C2*g_C1*Dr_C2+C1*g_C2*g_H*Dr_C2-half*g_B*g_C1&
          &*C2*Dr_C2-half*g_C2*g_C*C2*Dr_C2-half*g_C2*g_C1*Dr_H*C2&
          &-g_C*g_A*Dr_C*C1+3.0*half*C1*g_C2*g_A*Dr_C2+half*g_H*g_C&
          &*Dz_H*Dr_C1-half*g_B*g_C1*Dz_H*Dr_C)*r+(3.0*half*C1*g_H&
          &*g_A*Dr_H+half*C*g_C2*g_B*Dr_B-half*H*g_H*g_C*Dr_C2-half*H&
          &*g_H*g_C1*Dr_H+half*C*g_H*g_C2*Dr_H-half*g_C*g_A*Dr_A*C2&
          &-half*g_A**2*Dr_A*C1+half*A*g_H*g_C1*Dr_H-half*g_C1*Dr_H&
          &+half*C1*g_H**2*Dr_H+half*g_B*g_A*Dr_C2*C+half*g_C2*g_A&
          &*Dr_H*C-half*H*g_C1*g_A*Dr_A+half*A*g_C2*g_C*Dr_B-half*H&
          &*g_C2*g_C*Dr_B+half*A*g_H*g_C*Dr_C2-half*Dz_A*g_C2*g_A&
          &*Dr_H+half*A*g_C1*g_A*Dr_A-half*g_B*g_C*Dr_B*C2+half*C*g_H&
          &*g_B*Dr_C2-half*Dz_A*g_B*g_A*Dr_C2)/r+half*g_C**2*Dz_A*C2&
          &+C*g_C2**2*C2+C**2*g_B*g_C1-g_H*C1*(r**2*C1*g_C1+C2*g_C2+C&
          &*g_C+g_A*lambda)+half*g_A*g_C2*lambda*Dz_A-g_C**2*C2*C-g_B&
          &*g_C1*C2**2-g_C2*g_C*C2**2+half*g_C2*g_C*Dz_H*C-half*g_B&
          &*g_A*Dz_C1*C-half*g_B*g_C1*Dz_H*C-2.0*H*g_C2*g_C1*C2+2.0&
          &*C1*g_C2*g_H*C2-g_C*g_A*C*C1+3.0*half*C1*g_H*g_C*Dz_H+half&
          &*g_H*g_A*Dr_H*Dr_C1+half*g_B*g_A*Dr_C2*Dr_C+half*g_C2*g_A&
          &*Dr_H*Dr_C-half*(1.0-ft8)*(-g_C2*Dz_lambda+(g_C2*(g_A&
          &*lambda+r**2*C1*g_C1+C*g_C)-C1*g_A*g_C)*Dz_A)-half&
          &*ft8*(Dz_lambda*g_A*g_B*C2+(g_C2*(g_A*lambda+r**2*C1&
          &*g_C1+C*g_C)-C1*g_A*g_C)*Dz_H)+half*g_C1*g_C*Dz_C1*Dr_C1*r&
          &**5

     RIC_C2 = RIC_C2+half*C1*g_H**2*Dz_H-half*g_A**2*Dz_A*C1-half*g_A&
          &*g_C*Dz_A*C2+half*Dz_H*g_H*g_B*Dz_C2-half*g_B*g_C*Dz_B*C2&
          &-g_A*g_B*Dz_C*C2+half*g_C*g_C1*Dz_C1*Dr_C2*r**5+(half*C&
          &*g_C2*g_C1*Dz_C1-g_A*g_C1*Dz_C1*C1-half*g_C1**2*Dr_C2*Dr_H&
          &+half*g_C*g_C1*Dz_C1*C2-half*H*g_C1**2*Dz_C1-half*g_C**2&
          &*Dz_C*Dz_C1-half*g_C*g_C1*Dr_C2**2-g_C1*g_C*Dz_C2*C1+C1&
          &*g_C1*g_H*Dz_C1+half*C1*g_C1**2*Dz_A-half*g_C1**2*Dz_H*C1&
          &+C1*g_C1*g_C2*Dz_C-half*g_C2*g_C*Dz_C2*Dz_C1+half*g_C1*g_B&
          &*Dz_C2*Dz_C1+half*A*g_C1**2*Dz_C1)*r**4+(half*Dz_H*g_C1&
          &*g_B*Dz_C-3.0*half*g_C2*g_C1*Dz_H*C2+half*C1*g_C2**2*Dz_B&
          &-3.0*half*g_C2*g_A*Dz_C1*C2-half*Dz_C2*g_C2**2*Dz_H-H*g_C1&
          &*g_C2*Dz_C2+A*g_C1*g_C2*Dz_C2-H*g_C1**2*C2+A*g_C1**2*C2&
          &-g_A*g_C*Dz_C*C1+C*g_C2*g_C1*C2-half*Dr_B*g_C**2*Dr_C2+C&
          &*g_C*g_C2*Dz_C-half*Dr_B*g_C*g_C1*Dr_H-g_C1*g_B*Dz_C2*C2&
          &+half*g_C2*g_A*Dr_C2**2-half*g_C**2*Dz_B*C1-g_C*g_C1*C2**2&
          &+half*C*g_H*g_C*Dz_C1+half*C*g_C1*g_C*Dz_A+half*g_A*g_B&
          &*Dz_C*Dz_C1+C*g_C2**2*Dz_C2-half*Dz_H*g_C*g_C2*Dz_C+C1*g_H&
          &*g_C2*Dz_C2)*r**2+(-half*A*g_C2*g_A*Dr_B+half*H*g_H*g_A&
          &*Dr_C2+half*Dr_B*g_C2*g_A*H-half*A*g_H*g_A*Dr_C2+half*Dr_B&
          &*g_B*g_A*C2+g_H*g_A*C2*Dr_H)/r+(3.0*half*g_C2*g_A*C2*Dr_C2&
          &-half*C*g_C2*g_C*Dr_B+half*A*g_C1*g_C*Dr_B-half*g_A*g_C2&
          &*Dz_C*Dr_H+half*Dz_H*g_H*g_C*Dr_C2+half*C*g_C1*g_B*Dr_B&
          &-half*Dr_B*g_C*g_C1*H+half*g_H*g_C*Dz_C2*Dr_H-half*Dr_B&
          &*g_C**2*C2+half*Dr_B*g_C2*g_C*Dz_H-half*Dr_B*g_B*g_C1*Dz_H&
          &-half*g_A*g_B*Dz_C*Dr_C2-half*Dr_B*g_B*g_A*Dz_C1-half*C&
          &*g_H*g_C*Dr_C2)*r+(half*g_C2*g_C*Dz_C2*Dr_C2+half*C*g_C2&
          &*g_C1*Dr_C2+half*g_C**2*Dz_C*Dr_C2-half*H*g_C1**2*Dr_C2&
          &+half*A*g_C1**2*Dr_C2-3.0*half*g_C*g_C1*C2*Dr_C2-half*g_C2&
          &*g_A*Dz_C1*Dr_C2-half*g_C1*g_C2*Dz_C2*Dr_H-g_C1**2*C2*Dr_H&
          &+half*g_C1*g_C*Dz_C*Dr_H+half*Dr_B*g_C**2*Dz_C1-half*g_C1&
          &*g_B*Dz_C2*Dr_C2-half*g_C2*g_C1*Dz_H*Dr_C2)*r**3+half*g_H&
          &*g_A*Dr_C2*Dr_H+Dz_H*g_H*g_C*C2-half*H*g_H*g_C1*Dz_H-half&
          &*H*g_C1*g_A*Dz_A-half*H*g_C2*g_C*Dz_B+half*A*g_C1*g_A*Dz_A&
          &+half*A*g_C2*g_C*Dz_B+half*C*g_H*g_C2*Dz_H+A*g_A*g_C2*Dz_C&
          &+half*C*g_C2*g_B*Dz_B+half*A*g_H*g_C1*Dz_H+half*Dr_B*g_B&
          &*g_A*Dr_C2+half*A*g_H*g_A*Dz_C1-C*g_H*g_C*C2+half*Dr_B&
          &*g_C2*g_A*Dr_H-half*H*g_H*g_A*Dz_C1-H*g_A*g_C2*Dz_C+g_C2&
          &*g_A*C2**2-lambda*g_H*g_A*C2

  else

    RIC_A = RIC_A + Dz_A*g_C**2*Dr_C*r**3+(-half*Dz_A**2*g_C**2+Dz_A&
          &*g_C**2*C+g_B*g_A*Dr_C**2+half*Dr_A*g_C**2*Dr_B+Dr_A*g_C&
          &*g_A*Dr_C+g_B*g_C*Dr_C*Dr_B)*r**2+(-Dz_A*g_B*g_A*Dr_C+2.0&
          &*g_B*g_A*Dr_C*C+g_B*g_C*C*Dr_B+Dr_A*g_C*g_A*C)*r+1.0&
          &*fourth*g_H**2*Dr_H**2-g_B*g_A*C*Dz_A+fourth*g_B**2&
          &*Dr_B**2+g_B*g_A*C**2+fourth*Dr_A**2*g_A**2+half*g_B&
          &*g_A*Dz_A**2

     RIC_B = RIC_B + Dr_B*g_C**2*r**3*Dz_C+(g_A*g_B*Dz_C**2+half*Dz_A&
          &*g_C**2*Dz_B-half*g_C**2*Dr_B**2+Dz_B*g_C*g_B*Dz_C+g_A*g_C&
          &*Dz_C*Dz_A)*r**2-Dr_B*g_A*g_B*Dz_C*r+fourth*g_A**2*Dz_A**2&
          &+half*g_B*g_A*Dr_B**2+fourth*g_H**2*Dz_H**2+fourth*g_B**2&
          &*Dz_B**2

     RIC_H = RIC_H + (-half*C*g_C*g_H*Dr_H+g_H*g_C*Dz_H*Dr_H)*r+half&
          &*g_H*Dz_H**2*g_B+3.0*half*g_C*Dz_H+half*g_H*g_A*Dr_H&
          &**2-half*C*g_B*g_H*Dz_H-half*g_H*A*g_C*Dz_H-C*g_C&
          &-half*g_C*Dz_H+half*(2.0*g_A-g_H-g_H*A*g_A)*Dr_H/r&
          &+g_H*lambda*(1.0-g_A*H)

     RIC_C = RIC_C + g_C**2*Dr_C*Dz_C*r**3+(-half*Dz_A*g_C**2*Dz_C+g_C&
          &**2*C*Dz_C-half*g_C**2*Dr_B*Dr_C)*r**2+(-half*g_C**2*Dr_B&
          &*C+fourth*Dz_B*g_C**2*Dr_A+half*Dr_A*g_A*g_C*Dz_C+half*g_C&
          &*g_B*Dz_C*Dr_B+3.0*fourth*g_C**2*Dr_B*Dz_A+half*Dz_B*g_B&
          &*g_C*Dr_C+half*g_C*g_A*Dr_C*Dz_A)*r+half*Dz_A*g_A*g_B*Dz_C&
          &+half*Dz_B*g_B*g_C*C+half*g_C*g_A*C*Dz_A+half*Dr_B*g_B*g_A&
          &*Dr_C+(fourth*g_H**2*Dz_H*Dr_H+half*g_B*g_A*Dr_B*C+fourth&
          &*g_A**2*Dr_A*Dz_A+fourth*Dz_B*g_B**2*Dr_B-half*g_B*g_A&
          &*Dz_A*Dr_B)/r

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine calc_conformalRicci4

