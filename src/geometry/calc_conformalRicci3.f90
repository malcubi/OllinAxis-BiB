!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/calc_conformalRicci3.f90,v 1.4 2021/03/11 16:05:19 malcubi Exp $

  subroutine calc_conformalRicci3

! ***************************************************************
! ***   CALCULATES THIRD TERM OF THE CONFORMAL RICCI TENSOR   ***
! ***************************************************************

! Originally written by Jose Manuel Torres.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  real(8) half,third,fourth,one,two


! *******************
! ***   NUMBERS   ***
! *******************

  one = 1.d0
  two = 2.d0

  half  = 0.5d0
  third = 1.d0/3.d0
  fourth = half*half


! *************************************************
! ***   THIRD TERM FOR RICCI CONFORMAL TENSOR   ***
! *************************************************

! The expressions are rather long, we calculate the components using:
!
! _           _mn. . _    _           m          mn                  mn
! R  = - half*g  D D g  + g   D  Delta  + 2 Delta    Delta    + Delta   Delta
!  ab             m n ab   m(a b)                  (a     b)mn         a     mnb
!
! Here we calculate the third term.

  if (angmom) then

     RIC_A = RIC_A +(g_C2*g_C1*Dr_C1*Dz_C1-2.0*C1*g_C1**2*Dr_C1)*r**7&
          &+(-6.0*C1**2*g_C1**2-g_C2*g_C1*Dr_C1*Dr_C2+3.0*g_C2*g_C1&
          &*C1*Dz_C1)*r**6+(g_C2*g_C*Dr_C*Dz_C1+g_H*g_C*Dr_C1*Dz_C1&
          &+Dz_A*g_C1*g_C*Dr_C1-2.0*C1*g_C*g_C1*Dr_C-5.0*g_C2*g_C1*C1&
          &*Dr_C2+g_C1*g_C*Dr_A*Dz_C1+2.0*g_C2*g_C*Dr_C1*Dz_C+g_B&
          &*g_C1*Dr_C*Dz_C1+g_C2**2*Dr_C2*Dz_C1-2.0*g_C2*g_C1*Dr_C1&
          &*C2)*r**5+(Dr_A*g_C1*g_A*Dr_C1-g_C2**2*Dr_C2**2-g_C2*g_C&
          &*Dr_C*Dr_C2-g_C2*g_C*Dr_C1*Dr_B+2.0*g_C2**2*C2*Dz_C1-2.0&
          &*C1*g_C*g_C1*C+3.0*Dz_A*g_C1*g_C*C1+g_B*g_C1*C*Dz_C1&
          &-10.0*g_C2*g_C1*C1*C2+g_C2*g_C*C*Dz_C1-g_H*g_C*Dr_C1&
          &*Dr_C2+6.0*g_C2*g_C*C1*Dz_C-g_C2*g_C1*Dr_C*Dr_H+3.0*g_H&
          &*g_C*C1*Dz_C1-g_B*g_C1*Dr_C*Dr_C2-g_H*g_C1*Dr_C1*Dr_H)*r&
          &**4+(-g_C2*g_C1*C*Dr_H-5.0*g_H*g_C1*C1*Dr_H-2.0*g_C2*g_C1&
          &*Dr_C*H+Dr_A*g_C1*g_A*C1-5.0*g_C2*g_C*C1*Dr_B+Dz_A*g_C2&
          &*g_A*Dr_C1+2.0*g_B*g_C*Dr_C*Dz_C+g_C2*g_B*Dr_B*Dz_C1+2.0*A&
          &*g_C2*g_C1*Dr_C+g_H*g_B*Dz_C1*Dr_C2+Dz_A*g_B*g_C1*Dr_C2&
          &+Dz_A*g_C2*g_C*Dr_C2+Dz_A*g_C2*g_C1*Dr_H+2.0*A*g_H*g_C1&
          &*Dr_C1-5.0*g_H*g_C*C1*Dr_C2-2.0*g_H*g_C1*Dr_C1*H+2.0*g_B&
          &*g_C2*Dr_C2*Dz_C-2.0*C1*g_H*g_A*Dr_C1-2.0*C1*g_C2*g_A*Dr_C&
          &-g_C2*g_C*C*Dr_C2-g_B*g_C1*C*Dr_C2+g_H*g_C2*Dr_H*Dz_C1-2.0&
          &*g_B*g_C1*Dr_C*C2-2.0*g_C2*g_C*Dr_C*C2-2.0*g_H*g_C*Dr_C1&
          &*C2+A*g_C1**2*Dr_A+Dz_A*g_C**2*Dr_C-4.0*g_C2**2*C2*Dr_C2&
          &+g_C**2*Dr_A*Dz_C+g_C2**2*Dr_H*Dz_C)*r**3+(2.0*g_C1*C1&
          &+Dr_A*g_C*g_A*Dr_C-2.0*g_B*g_C2*Dr_C2*Dr_B-2.0*g_C2*g_C*C&
          &*C2+2.0*g_H*g_B*Dz_C1*C2+3.0*Dz_A*g_C2*g_A*C1+2.0*g_B*g_C&
          &*C*Dz_C-g_C2*Dz_C1+2.0*Dz_A*g_C2*g_C1*H+2.0*Dz_A*g_B*g_C1&
          &*C2+2.0*A*g_C2*g_C1*C-2.0*g_C2*g_C1*C*H-10.0*g_H*g_C1*C1&
          &*H-g_B*g_C*Dr_C*Dr_B+6.0*A*g_H*g_C1*C1+4.0*g_B*g_C2*C2&
          &*Dz_C+2.0*Dz_A*g_C2*g_C*C2-2.0*C1*g_C2*g_A*C+2.0*g_H*g_C2&
          &*H*Dz_C1-2.0*g_B*g_C1*C*C2-10.0*g_H*g_C*C1*C2-2.0*g_H&
          &*g_C2*Dr_H*Dr_C2-4.0*g_C2**2*C2**2+2.0*g_C2**2*H*Dz_C-6.0&
          &*C1**2*g_H*g_A+Dz_A*g_C**2*C-g_C2**2*Dr_H*Dr_B-g_H*g_B&
          &*Dr_C2**2)*r**2+(2.0*A*g_C2*g_H*Dr_C2+A*g_C2**2*Dr_B+Dz_A&
          &*g_B*g_C*Dr_B-4.0*g_B*g_C2*C2*Dr_B+Dr_A*g_C*g_A*C-4.0*g_H&
          &*g_B*C2*Dr_C2-4.0*g_H*g_C2*H*Dr_C2+Dz_A*g_B*g_A*Dr_C-2.0&
          &*g_C2**2*H*Dr_B-4.0*g_H*g_C2*Dr_H*C2+Dz_A*g_C*g_A*Dr_A-g_B&
          &*g_C*C*Dr_B+g_B**2*Dr_B*Dz_C+g_C2*Dr_C2)*r+Dz_A*g_B*g_A*C&
          &-4.0*g_H*g_B*C2**2-8.0*g_H*g_C2*H*C2-half*g_B**2*Dr_B**2&
          &+2.0*g_C2*C2-half*g_H**2*Dr_H**2+4.0*A*g_C2*g_H*C2+half&
          &*Dr_A**2*g_A**2+(g_H*Dr_H+A*g_H**2*Dr_H-2.0*g_H**2*Dr_H*H)&
          &/r-2.0*g_H*lambda*(1.0-g_H*H)

     RIC_B = RIC_B -r**8*g_C1**2*Dz_C1**2+r**7*g_C1**2*Dz_C1*Dr_C2+(&
          &-g_C1*g_C2*Dz_C2*Dz_C1-g_C1*g_C*Dz_C*Dz_C1)*r**6+(Dr_B*g_C&
          &*g_C1*Dz_C1+g_C1**2*Dz_H*Dr_C+2.0*g_C1*g_C*Dz_C2*Dr_C+g_C1&
          &*g_C2*Dz_C2*Dr_C2+2.0*g_A*g_C1*Dz_C1*Dr_C+g_C1*g_C*Dz_C&
          &*Dr_C2)*r**5+(2.0*C*g_C1*g_C2*Dz_C-g_C1**2*Dz_H*Dz_A+g_C1&
          &**2*Dz_H*C-Dz_H*g_C1*g_C2*Dz_C-g_C1*g_C*Dz_C2*Dz_A-g_H*g_A&
          &*Dz_C1**2+2.0*C*g_C1*g_H*Dz_C1-g_H*g_C*Dz_C2*Dz_C1-2.0*g_H&
          &*g_C1*Dz_H*Dz_C1-g_A*g_C2*Dz_C*Dz_C1+2.0*g_A*g_C1*Dz_C1*C&
          &-2.0*g_A*g_C1*Dz_C1*Dz_A+2.0*g_C1*g_C*Dz_C2*C+C*g_C1**2&
          &*Dz_A)*r**4+(g_A*g_C2*Dz_C*Dr_C2+g_C2*g_C*Dz_B*Dr_C2+Dr_B&
          &*g_C1*g_B*Dz_C2+Dr_B*g_C**2*Dz_C+g_H*g_C1*Dz_H*Dr_C2+2.0&
          &*g_A*g_C*Dz_C*Dr_C+g_H*g_A*Dz_C1*Dr_C2+Dr_B*g_C2*g_A*Dz_C1&
          &+Dr_B*g_C2*g_C*Dz_C2+Dr_B*g_C2*g_C1*Dz_H+g_C1*g_A*Dz_A&
          &*Dr_C2+g_H*g_C*Dz_C2*Dr_C2+g_C**2*Dz_B*Dr_C)*r**3+(2.0*g_A&
          &*g_C*Dz_C*C+Dz_B*g_C*g_B*Dz_C+Dz_B*g_C2*g_B*Dz_C2+g_C**2&
          &*Dz_B*C+C*g_C2**2*Dz_B+2.0*C*g_H*g_C2*Dz_C2-g_A*g_C*Dz_C&
          &*Dz_A-Dz_H*g_H*g_C2*Dz_C2)*r**2+(Dr_B*g_A*g_B*Dz_C+Dr_B&
          &*g_B*g_C*Dz_B+g_A**2*Dz_A*Dr_C+Dr_B*g_A*g_C*Dz_A)*r+g_A**2&
          &*Dz_A*C-half*g_H**2*Dz_H**2-half*g_A**2*Dz_A**2+half*g_B&
          &**2*Dz_B**2+C*g_H**2*Dz_H

     RIC_H = RIC_H +2.0*g_C1*Dr_C1*r**3+(4.0*g_C1*C1+g_C2*Dz_C1)*r**2&
          &+g_C2*Dr_C2*r+2.0*g_C2*C2-g_C*Dz_H+(g_H*Dr_H-g_A*Dr_H)/r

     RIC_C = RIC_C +half*Dz_A*g_B*g_C*Dz_B+A*g_H*g_C2*Dz_C2-half*g_C&
          &*g_A*C*Dz_A+half*Dr_B*g_B*g_A*Dr_C+C*g_C2*g_H*C2-H*g_H&
          &*g_C2*Dz_C2-2.0*C2*g_H*g_C2*Dz_H+half*Dr_B*g_C*g_A*Dr_A&
          &+half*Dz_B*g_B*g_C*C+half*Dz_A*g_A*g_B*Dz_C-g_H*g_B*Dz_C2&
          &*C2-half*g_C1**2*Dr_C1*Dz_C1*r**7+half*g_A**2*Dr_A*Dr_C&
          &+half*g_B**2*Dz_B*Dz_C+half*Dz_A**2*g_A*g_C+g_C*g_A*C**2&
          &+half*Dr_B**2*g_B*g_C+half*A*g_C2**2*Dz_B+(half*Dr_B*g_B&
          &*g_A*C+half*g_A**2*Dr_A*C-half*Dz_H*g_H**2*Dr_H+half*C*g_H&
          &**2*Dr_H)/r+(-5.0*half*C1*g_C1**2*Dz_C1+half*g_C2*g_C1&
          &*Dz_C1**2+half*g_C1**2*Dr_C1*Dr_C2)*r**6+(3.0*half*C1*g_C1&
          &**2*Dr_C2-g_C2*g_C1*Dr_C2*Dz_C1-half*g_C*g_C1*Dr_C*Dz_C1)&
          &*r**5+(-C1*g_C1*g_C2*Dz_C2-half*g_C*g_C1*C*Dz_C1+g_C1*g_A&
          &*Dr_C1*Dr_C+half*g_H*g_C*Dz_C1**2+3.0*half*g_C2*g_C*Dz_C1&
          &*Dz_C+half*g_C1*g_B*Dz_C*Dz_C1+half*Dr_B*g_C1*g_C*Dr_C1&
          &+half*g_C2*g_C1*Dr_C2**2+half*g_C2**2*Dz_C2*Dz_C1+3.0*half&
          &*g_C*g_C1*Dr_C2*Dr_C-2.0*C2*g_C2*g_C1*Dz_C1-C1*g_C1*g_C&
          &*Dz_C+half*g_C1**2*Dr_H*Dr_C+g_C1*g_C*Dz_A*Dz_C1)*r**4+(&
          &-half*g_C2**2*Dz_C2*Dr_C2+g_C1**2*H*Dr_C+half*g_C1**2*Dr_H&
          &*C+half*C*g_C1**2*Dr_A-half*g_C1**2*Dr_H*Dz_A-half*Dz_H&
          &*g_H*g_C1*Dr_C1+2.0*g_C*g_C1*C2*Dr_C+C*g_C2*g_C1*Dr_C+C&
          &*g_H*g_C1*Dr_C1-half*Dz_H*g_C2*g_C1*Dr_C-half*g_C1*g_C2&
          &*Dz_C*Dr_H+3.0*half*g_C*g_C1*Dr_C2*C-g_C*g_C1*Dr_C2*Dz_A&
          &-half*g_C1*g_A*Dr_C1*Dz_A+3.0*half*C1*g_C1*g_C*Dr_B-Dr_B&
          &*g_C2*g_C*Dz_C1-half*g_C2*g_A*Dr_C*Dz_C1-g_H*g_C*Dz_C1&
          &*Dr_C2+half*Dr_A*g_C1*g_C*Dz_C2+C2*g_C2*g_C1*Dr_C2-half&
          &*g_H*g_A*Dr_C1*Dz_C1-g_H*g_C1*Dr_H*Dz_C1+3.0*g_C1*g_A*C1&
          &*Dr_C+g_C1*g_A*Dr_C1*C-half*g_C1*g_B*Dz_C*Dr_C2-half*g_C&
          &*g_C2*Dz_C*Dr_C2+half*Dz_B*g_C2*g_C*Dr_C1)*r**3+(Dr_B*g_C2&
          &*g_C1*H+3.0*half*C1*g_H*g_A*Dr_C2+half*Dr_A*g_A*g_C*Dz_C&
          &-half*g_C1*Dr_C2-half*g_H*g_C2*Dz_C2*Dr_H+2.0*g_C*g_A*Dr_C&
          &*C-Dz_H*g_C2*g_H*Dr_C2+C2*g_H*g_C*Dr_C2+C*g_C2*g_H*Dr_C2&
          &-half*g_C*g_B*Dz_C*Dr_B-half*Dz_H*g_C2**2*Dr_B+g_H*g_C1*H&
          &*Dr_C2-half*g_C2*g_B*Dz_C2*Dr_B+C2*g_C1*g_B*Dr_B+3.0*half&
          &*C1*g_C2*g_A*Dr_B+half*C*g_C2**2*Dr_B-half*g_C**2*Dr_B&
          &*Dz_A-half*g_H*g_B*Dz_C2*Dr_C2-half*g_C*g_A*Dr_C*Dz_A+C2&
          &*g_C2*g_C*Dr_B+half*g_C2*g_A*C*Dr_C2+half*Dz_B*g_C**2*Dr_A&
          &+g_C**2*Dr_B*C+half*Dz_B*g_B*g_C*Dr_C)*r+half*g_H*Dz_H&
          &*(lambda*g_H+g_C2*C2)+(2.0*C*g_H*g_C1*C1+half*g_C1*g_A&
          &*Dr_A*Dr_C2+half*g_H*g_C1*Dr_H*Dr_C2+half*g_C2*g_B*Dz_B&
          &*Dz_C1+half*g_H*g_C2*Dz_H*Dz_C1+half*Dz_A*g_C2*g_C1*Dz_H&
          &-g_C1*g_B*Dz_C*C2-half*Dz_H*g_C2*g_C1*C-C1*g_A*g_C2*Dz_C&
          &-5.0*half*C1*g_C1*g_A*Dz_A-5.0*half*C1*g_H*g_A*Dz_C1+half&
          &*C1*g_C2*g_C*Dz_B+half*g_C2*g_A*Dr_C*Dr_C2-C1*g_H*g_C&
          &*Dz_C2+half*Dz_A*g_C2*g_C*Dz_C2-2.0*C1*g_H*g_C1*Dz_H-C2&
          &*g_C*g_C2*Dz_C-2.0*C2*g_C1*g_C*Dz_A-2.0*C2*g_H*g_C*Dz_C1&
          &-half*g_C2*g_A*C*Dz_C1+half*g_H*g_A*Dr_C1*Dr_C2-H*g_C1&
          &*g_C2*Dz_C-2.0*H*g_C1*g_H*Dz_C1+3.0*g_C1*g_A*C1*C+Dr_B&
          &*g_C2*g_C*Dr_C2+half*g_H*g_B*Dz_C2*Dz_C1+half*Dz_A*g_C1&
          &*g_B*Dz_C2+A*g_C1*g_C2*Dz_C+A*g_C1*g_H*Dz_C1+g_C2*g_B&
          &*Dz_C2*Dz_C+half*Dr_B*g_B*g_C1*Dr_C2+half*Dr_B*g_C2*g_C1&
          &*Dr_H+half*Dr_B*g_C2*g_A*Dr_C1+half*g_C1*Dz_C1+Dz_A*g_C**2&
          &*Dz_C+g_C2*g_C1*C**2+g_C*g_B*Dz_C**2+g_C*g_A*Dr_C**2-H&
          &*g_C1**2*Dz_A+half*A*g_C1**2*Dz_A+half*g_C2**2*Dz_H*Dz_C&
          &+half*g_H*g_C*Dr_C2**2-C2*g_C2**2*Dz_C2+g_C**2*Dr_B*Dr_C&
          &+g_C1**2*H*C+2.0*g_C*g_C1*C2*C+half*g_C2*g_A*Dz_A*Dz_C1)*r&
          &**2

     RIC_C1= RIC_C1+C2*(g_C*(r**2*C1*g_C1+C*g_C)-C1*g_A*g_C2)-g_C1&
          &*lambda+(3.0*half*C1*g_C1*g_C*Dz_C1+3.0*half*g_C1*g_C&
          &*Dr_C2*Dr_C1+g_C1**2*Dr_H*Dr_C1+g_C1*g_A*Dr_C1**2)*r**4&
          &+(half*g_C2*g_C*Dr_C2*Dz_C1+3.0*g_C1*g_C*C2*Dr_C1+6.0*g_C1&
          &*g_A*C1*Dr_C1+C1*g_C2*g_C1*Dr_C+half*g_C2*g_C1*Dr_H*Dz_C1&
          &+C1*g_H*g_C1*Dr_C1+3.0*g_C1**2*Dr_H*C1+g_C1**2*H*Dr_C1&
          &+half*g_C2*g_A*Dz_C1*Dr_C1+Dz_C2*g_C2*g_C*Dr_C1+half*g_C2&
          &*g_C1*Dz_H*Dr_C1+9.0*half*g_C*g_C1*Dr_C2*C1+half*g_C**2&
          &*Dz_C1*Dr_C+half*g_B*g_C1*Dz_C1*Dr_C2+half*C1*g_C1**2&
          &*Dr_A)*r**3+(H*g_C2*g_C1*Dz_C1+g_C2*g_C*Dz_C1*C2+3.0*half&
          &*C1*g_C2*g_A*Dz_C1+half*g_C2*g_C*Dr_C2**2+half*g_C**2*Dr_B&
          &*Dr_C1+g_C2*g_C1*Dr_H*Dr_C2+9.0*g_C1*g_A*C1**2+3.0*half*C1&
          &*g_C2*g_C1*Dz_H+half*g_C*g_C1*Dr_H*Dr_C+3.0*Dz_C2*g_C2*g_C&
          &*C1+half*g_C2*g_A*Dr_C2*Dr_C1+half*g_C**2*Dz_C1*C+9.0*g_C&
          &*g_C1*C2*C1+g_C*g_A*Dr_C*Dr_C1+half*g_C**2*Dr_C2*Dr_C+g_B&
          &*g_C1*C2*Dz_C1+3.0*C1**2*g_H*g_C1+half*g_B*g_C1*Dr_C2**2&
          &+3.0*H*g_C1**2*C1+C1*g_C2*g_C1*C)*r**2+(half*Dz_C2*g_C**2&
          &*Dr_A+half*Dz_C2*g_C2**2*Dr_H+half*Dz_H*g_C2**2*Dr_C2+half&
          &*C1*g_C2**2*Dr_B+3.0*half*g_C**2*Dr_B*C1+half*g_C**2*Dr_C2&
          &*C+g_C**2*C2*Dr_C+half*g_B*g_C*Dr_B*Dz_C1+Dz_C2*g_B*g_C2&
          &*Dr_C2+half*g_C2*g_C*Dz_H*Dr_C+half*g_B*g_A*Dz_C1*Dr_C&
          &+Dz_C2*g_B*g_C*Dr_C+half*g_C*g_C1*Dr_H*C+H*g_C2*g_C1*Dr_C2&
          &+half*Dz_H*g_C1*g_C*Dr_A+g_C*g_A*C*Dr_C1+C1*g_C2*g_H*Dr_C2&
          &+2.0*g_B*g_C1*C2*Dr_C2+g_C2*g_A*Dr_C1*C2+2.0*g_C2*g_C*C2&
          &*Dr_C2+2.0*g_C2*g_C1*Dr_H*C2+3.0*g_C*g_A*Dr_C*C1+3.0*half&
          &*C1*g_C2*g_A*Dr_C2+half*g_H*g_C*Dz_H*Dr_C1+half*g_B*g_C1&
          &*Dz_H*Dr_C+half*g_C*g_A*Dr_A*Dz_C1)*r+(half*Dz_H*g_H*g_C2&
          &*Dr_H+3.0*half*C1*g_H*g_A*Dr_H+half*g_C1*Dr_A+half*g_B*g_A&
          &*Dr_C2*C+half*g_C2*g_A*Dr_H*C+g_B*g_C*Dr_B*C2+half*Dz_H&
          &*g_H*g_B*Dr_C2+half*Dz_C2*g_B**2*Dr_B+C2*g_H*g_C*Dr_H+g_B&
          &*g_A*Dr_C*C2+H*g_H*g_C1*Dr_H+half*g_C*Dr_C2+half*Dz_H*g_C2&
          &*g_B*Dr_B+g_C*g_A*Dr_A*C2+half*C1*g_H**2*Dr_H+3.0*half*g_A&
          &**2*Dr_A*C1)/r+half*g_C2*Dz_lambda+C1*g_lambda+C2*g_C2**2&
          &*Dz_H+Dz_C2*g_C2**2*H+half*g_A**2*Dr_A*Dr_C1+g_C**2*C2*C&
          &+2.0*g_B*g_C1*C2**2+2.0*g_C2*g_C*C2**2+half*g_H*g_C1*Dr_H&
          &**2-C1*g_H*(r**2*C1*g_C1+C2*g_C2)+half*g_C2*g_C*Dz_H*C&
          &+Dz_C2*g_B*g_C*C+half*g_B*g_A*Dz_C1*C+half*g_B*g_C1*Dz_H*C&
          &+2.0*H*g_C2*g_C1*C2+2.0*C1*g_C2*g_H*C2+3.0*g_C2*g_A*C1*C2&
          &+3.0*g_C*g_A*C*C1+half*g_C2*g_C*Dr_B*Dr_H+half*C1*g_H*g_C&
          &*Dz_H+half*g_H*g_A*Dr_H*Dr_C1+half*g_H*g_C*Dr_C2*Dr_H+2.0&
          &*Dz_C2*g_B*g_C2*C2+half*g_C1*g_A*Dr_A*Dr_H+half*g_C*g_A&
          &*Dr_A*Dr_C2+half*g_B*g_C*Dr_B*Dr_C2+half*g_B*g_A*Dr_C2&
          &*Dr_C+half*g_C2*g_A*Dr_H*Dr_C-half*g_C*Dz_C1+half*g_C1*g_C&
          &*Dz_C1*Dr_C1*r**5
     
     RIC_C2= RIC_C2+half*g_C*g_C1*Dz_C1**2*r**6+(half*g_C*g_C1*Dz_C1&
          &*Dr_C2+half*g_C1**2*Dz_C1*Dr_H+g_C1*g_C*Dz_C2*Dr_C1+half&
          &*g_C1**2*Dz_H*Dr_C1+g_A*g_C1*Dz_C1*Dr_C1)*r**5+(g_C2*g_C1&
          &*Dz_H*Dz_C1+C1*g_C1*g_H*Dz_C1+3.0*g_C1*g_C*Dz_C2*C1+half&
          &*g_C**2*Dz_C*Dz_C1+half*g_C2*g_A*Dz_C1**2+C1*g_C1*g_C2&
          &*Dz_C+3.0*half*g_C2*g_C*Dz_C2*Dz_C1+3.0*g_A*g_C1*Dz_C1*C1&
          &+3.0*half*g_C1**2*Dz_H*C1+g_C*g_C1*Dz_C1*C2+half*g_C1*g_B&
          &*Dz_C2*Dz_C1+half*C1*g_C1**2*Dz_A)*r**4+(g_A*g_C*Dz_C&
          &*Dr_C1+half*g_C**2*Dz_B*Dr_C1+half*g_C2*g_C*Dz_C2*Dr_C2&
          &+half*g_C1*g_C*Dz_C*Dr_H+half*g_C1*g_B*Dz_C2*Dr_C2+half&
          &*g_C2*g_A*Dz_C1*Dr_C2+half*g_C1*g_C2*Dz_C2*Dr_H+half*g_C&
          &**2*Dz_C*Dr_C2+half*g_C2*g_C1*Dz_H*Dr_C2)*r**3+(g_C2*g_A&
          &*Dz_C1*C2+3.0*half*g_C**2*Dz_B*C1+g_C2*g_C*Dz_C2*C2+Dz_C2&
          &*g_C*g_B*Dz_C+half*Dz_C2*g_C**2*Dz_A+g_C1*g_B*Dz_C2*C2&
          &+half*Dz_H*g_C*g_C2*Dz_C+g_C2*g_C1*Dz_H*C2+half*g_B*g_C&
          &*Dz_B*Dz_C1+half*g_A*g_B*Dz_C*Dz_C1+C1*g_H*g_C2*Dz_C2+g_C&
          &**2*Dz_C*C2+half*Dz_H*g_C1*g_B*Dz_C+half*g_A*g_C*Dz_A&
          &*Dz_C1+Dz_C2*g_C2**2*Dz_H+3.0*g_A*g_C*Dz_C*C1+half*Dz_H&
          &*g_H*g_C*Dz_C1+Dz_C2**2*g_C2*g_B+half*Dz_H*g_C1*g_C*Dz_A&
          &+half*C1*g_C2**2*Dz_B)*r**2+(half*g_H*g_A*Dz_C1*Dr_H+half&
          &*g_C2*g_C*Dz_B*Dr_H+half*g_A*g_C2*Dz_C*Dr_H+half*g_A*g_B&
          &*Dz_C*Dr_C2+half*g_H*g_C1*Dz_H*Dr_H+g_C1*Dr_C+half*g_B*g_C&
          &*Dz_B*Dr_C2+half*g_H*g_C*Dz_C2*Dr_H+half*g_C1*g_A*Dz_A&
          &*Dr_H+half*g_A**2*Dz_A*Dr_C1+half*g_A*g_C*Dz_A*Dr_C2)*r&
          &-half*g_H*Dz_C1+half*C1*g_H**2*Dz_H+3.0*half*g_A**2*Dz_A&
          &*C1+half*Dz_H*g_C2*g_B*Dz_B+g_A*g_C*Dz_A*C2+half*Dz_H*g_H&
          &*g_B*Dz_C2+g_B*g_C*Dz_B*C2+half*Dz_C2*g_B**2*Dz_B+half&
          &*Dz_H**2*g_H*g_C2+g_A*g_B*Dz_C*C2+half*g_C1*Dz_H-half*g_C1&
          &*Dz_A+half*g_A*Dz_C1+(half*g_H*Dr_C2-half*g_A*Dr_C2+half&
          &*g_C2*Dr_B)/r

     RIC_lambda = RIC_lambda + half*Dr_H/r*(3.0*g_A*(C1*g_C1*r**2+C2&
          &*g_C2)+2.0*lambda*g_H**2+g_H*(C2*g_C2-C*g_C))-g_H*lambda&
          &*(3.0*C1*g_C1*r**2+3.0*C2*g_C2-g_lambda*H)+(Dz_A*g_B*g_C&
          &*Dr_B+2.0*g_B*g_A*Dr_C*C-3.0*half*g_B*g_A*C2*Dr_C2+A*g_C2&
          &**2*Dr_B+g_B**2*Dr_B*Dz_C-5.0*half*g_H*g_C2*Dr_H*C2-3.0&
          &*g_H*g_C2*H*Dr_C2+2.0*Dr_A*g_C*g_A*C+Dz_A*g_C*g_A*Dr_A+2.0&
          &*A*g_C2*g_H*Dr_C2-3.0*half*g_C2*g_A*Dr_H*C2-g_C*g_H*Dz_H&
          &*Dr_H+half*C*g_C*g_H*Dr_H+half*A*g_A*g_C2*Dr_C2-3.0*half&
          &*g_C2*g_A*H*Dr_C2-5.0*half*g_B*g_C2*C2*Dr_B-5.0*half*g_H&
          &*g_B*C2*Dr_C2-3.0*half*g_C2**2*H*Dr_B)/r-3.0*C1*g_C1**2&
          &*Dr_C1*r**5-3.0*g_C2**2*C2**2+half*Dz_H**2*g_C2**2-3.0&
          &*g_C2*Dz_C1-half*g_C2**2*Dr_H*Dr_B+half*Dr_A*g_C**2*Dr_B&
          &-half*g_H*g_B*Dr_C2**2-half*g_B*g_A*Dr_C2**2+3.0*C1**2*g_H&
          &*g_A-A*g_C1**2*H-half*C*g_C2**2*Dz_H+g_B*g_A*Dr_C**2+2.0&
          &*Dz_A*g_C**2*C-half*Dz_A**2*g_C**2+half*g_B*(ft7&
          &*(g_lambda*Dz_A**2+g_H*Dz_lambda*(Dz_A+Dz_H))+(1.0-ft7)&
          &*(g_A*Dz_lambda*(Dz_A+Dz_H)+g_lambda*Dz_H**2))+(half*g_C1&
          &**2*Dr_H**2-3.0*C1*g_C*g_C1*C+6.0*g_C2*g_C*C1*Dz_C+half&
          &*Dr_A*g_C1**2*Dr_H-8.0*g_C2*g_C1*C1*C2+g_C2*g_C*Dz_H*Dz_C1&
          &+g_C*g_C1*Dr_H*Dr_C2-3.0*half*g_C*g_C1*H*Dz_C1+3.0*half&
          &*g_C2*g_C*C*Dz_C1+13.0*half*g_H*g_C*C1*Dz_C1+13.0*half&
          &*Dz_A*g_C1*g_C*C1+half*A*g_C*g_C1*Dz_C1-g_B*g_C1*Dz_H&
          &*Dz_C1+2.0*Dr_A*g_C1*g_A*Dr_C1+3.0*half*g_C2**2*C2*Dz_C1&
          &-half*g_B*g_A*Dz_C1**2+g_H*g_A*Dr_C1**2-3.0*half*g_C**2*C2&
          &*Dz_C1-half*g_C2**2*Dr_C2**2+half*g_C**2*Dr_C2**2+half*g_B&
          &*g_H*Dz_C1**2+half*g_B*g_C1*C*Dz_C1+Dz_A*g_B*g_C1*Dz_C1&
          &+Dr_A*g_C*g_C1*Dr_C2+2.0*g_C2*g_A*Dr_C1*Dr_C-Dz_A*g_C*g_C2&
          &*Dz_C1)*r**2-2.0*g_B*g_C1*C*C2-3.0*half*g_C2*g_C*Dz_H*C2&
          &+2.0*g_C*g_C1*H*C2+3.0*half*g_B*g_C1*Dz_H*C2-g_C2*g_A*Dr_H&
          &*Dr_C2+4.0*C1*g_C2*g_A*C+6.0*A*g_H*g_C1*C1-half*A*g_A*g_C2&
          &*Dz_C1-A*g_C*g_C1*C2+2.0*Dr_A*g_C*g_A*Dr_C-2.0*g_C2*g_C1*C&
          &*H-8.0*g_H*g_C1*C1*H+3.0*half*g_B*g_A*Dz_C1*C2-g_H*g_C2&
          &*Dr_H*Dr_C2+2.0*g_H*g_C2*H*Dz_C1-5.0*g_H*g_C*C1*C2+3.0&
          &*half*g_C*(C1*g_C1*r**2+C2*g_C2)*Dz_H+(g_B*g_A*Dz_C1*Dr_C2&
          &+g_C2*g_A*Dr_H*Dz_C1+3.0*half*g_C*g_C1*Dr_H*C2+2.0*g_C2&
          &*g_A*Dr_C1*C+4.0*C1*g_C2*g_A*Dr_C+g_C1*g_C2*Dz_H*Dr_H-half&
          &*A*g_C*g_C1*Dr_C2+3.0*half*g_C*g_C1*H*Dr_C2+2.0*A*g_H*g_C1&
          &*Dr_C1+half*g_C2*g_C*C*Dr_C2-5.0*half*g_H*g_C*C1*Dr_C2+7.0&
          &*half*Dr_A*g_C1*g_A*C1+g_B*g_C1*Dz_H*Dr_C2-g_H*g_C1*Dr_C1&
          &*H-5.0*half*g_H*g_C1*C1*Dr_H-g_C2*g_C1*Dr_C*H-g_B*g_C1&
          &*Dr_C*C2+g_C2**2*Dr_H*Dz_C+A*g_C1**2*Dr_A+half*Dr_A*g_C1&
          &**2*H+3.0*half*g_C**2*C2*Dr_C2+3.0*half*g_C1**2*Dr_H*H-5.0&
          &*half*g_C2**2*C2*Dr_C2-half*A*g_C1**2*Dr_H+2.0*Dz_A*g_C**2&
          &*Dr_C+g_C**2*Dr_A*Dz_C-half*g_B*g_C1*C*Dr_C2+Dz_A*g_B*g_C1&
          &*Dr_C2-half*g_C2*g_C1*C*Dr_H+2.0*g_B*g_C*Dr_C*Dz_C+g_C2&
          &*g_B*Dr_B*Dz_C1+4.0*C1*g_H*g_A*Dr_C1+2.0*A*g_C2*g_C1*Dr_C&
          &+2.0*g_B*g_C2*Dr_C2*Dz_C+g_H*g_B*Dz_C1*Dr_C2+g_H*g_C2*Dr_H&
          &*Dz_C1+Dz_A*g_C2*g_C*Dr_C2+Dz_A*g_C2*g_C1*Dr_H-3.0*g_C1&
          &*Dr_C1-5.0*half*g_C2*g_C*C1*Dr_B+3.0*half*Dr_A*g_C*g_C1*C2&
          &-g_C2*g_C*Dz_H*Dr_C2)*r+(-3.0*C1*g_C*g_C1*Dr_C+2.0*g_C2&
          &*g_C*Dr_C*Dz_C1+2.0*g_H*g_C*Dr_C1*Dz_C1-g_C2*g_C1*Dr_C1*C2&
          &+2.0*Dz_A*g_C1*g_C*Dr_C1-g_C*g_C1*Dr_H*Dz_C1+g_C1*g_C*Dr_A&
          &*Dz_C1+2.0*g_C2*g_C*Dr_C1*Dz_C+g_C2**2*Dr_C2*Dz_C1-5.0&
          &*half*g_C2*g_C1*C1*Dr_C2-g_C**2*Dr_C2*Dz_C1)*r**3+(3.0/4.0&
          &*Dr_A**2*g_A**2-half*g_A*g_H*Dr_H**2-1.0/4.0*g_H**2*Dr_H&
          &**2-1.0/4.0*g_B**2*Dr_B**2)/r**2+C*(C*g_C**2+H*g_C1*g_C2)&
          &+g_C1**2*H**2+3.0*half*Dz_A*g_C2*g_C*C2+2.0*Dz_A*g_C2*g_C1&
          &*H+5.0*half*g_H*g_B*Dz_C1*C2+2.0*g_B*g_C*C*Dz_C+5.0*half&
          &*Dz_A*g_B*g_C1*C2-half*Dz_A*g_C2*g_A*C1+2.0*A*g_C2*g_C1*C&
          &+4.0*g_B*g_C2*C2*Dz_C-g_B*g_C2*Dr_C2*Dr_B+g_C2*g_C*C*C2&
          &+2.0*g_C2**2*H*Dz_C+g_C**2*C2**2-3.0*g_C1*C1+(-9.0*C1**2&
          &*g_C1**2-half*g_C2**2*Dz_C1**2-half*g_C2*g_C1*C1*Dz_C1&
          &+half*g_C**2*Dz_C1**2)*r**4+(g_A+3.0*g_H)*C2*C1*g_C+(4.0&
          &*g_H+g_A)*C2*g_C2*lambda-half*g_H*C1*g_C2*Dz_H-half*A*g_C1&
          &*g_C2*Dz_H+3.0*half*g_C1*g_C2*Dz_H*H+3.0*half*g_C2*g_A*H&
          &*Dz_C1

  else

     RIC_A = RIC_A + (2.0*g_B*g_C*Dr_C*Dz_C+g_C**2*Dr_A*Dz_C+Dz_A*g_C&
          &**2*Dr_C)*r**3+(2.0*g_B*g_C*C*Dz_C-g_B*g_C*Dr_C*Dr_B+Dr_A&
          &*g_C*g_A*Dr_C+Dz_A*g_C**2*C)*r**2+(Dz_A*g_B*g_C*Dr_B+Dz_A&
          &*g_B*g_A*Dr_C+Dr_A*g_C*g_A*C-g_B*g_C*C*Dr_B+Dz_A*g_C*g_A&
          &*Dr_A+g_B**2*Dr_B*Dz_C)*r-half*g_H**2*Dr_H**2-half*g_B**2&
          &*Dr_B**2+g_B*g_A*C*Dz_A+half*Dr_A**2*g_A**2+(g_H*Dr_H+g_H&
          &**2*A*Dr_H-2.0*g_H**2*Dr_H*H)/r

     RIC_B = RIC_B + (Dr_B*g_C**2*Dz_C+g_C**2*Dz_B*Dr_C+2.0*g_A*g_C&
          &*Dz_C*Dr_C)*r**3+(2.0*g_A*g_C*Dz_C*C-g_A*g_C*Dz_C*Dz_A+g_C&
          &**2*Dz_B*C+Dz_B*g_C*g_B*Dz_C)*r**2+(Dr_B*g_A*g_C*Dz_A+g_A&
          &**2*Dz_A*Dr_C+Dr_B*g_A*g_B*Dz_C+Dr_B*g_B*g_C*Dz_B)*r+half&
          &*g_B**2*Dz_B**2+g_A**2*Dz_A*C-half*g_A**2*Dz_A**2+g_H**2&
          &*Dz_H*C-half*g_H**2*Dz_H**2

     !RIC_H = RIC_H - g_C*Dz_H + (g_H-g_A)*Dr_H/r
     RIC_H = RIC_H - g_C*Dz_H - r*g_lambda*Dr_H

     RIC_C = RIC_C + (g_C**2*Dr_B*Dr_C+g_C*g_B*Dz_C**2+Dz_A*g_C**2&
          &*Dz_C+g_C*g_A*Dr_C**2)*r**2+(-half*g_C*g_A*Dr_C*Dz_A+2.0&
          &*g_C*g_A*Dr_C*C-half*g_C**2*Dr_B*Dz_A-half*g_C*g_B*Dz_C&
          &*Dr_B+half*Dz_B*g_B*g_C*Dr_C+g_C**2*Dr_B*C+half*Dz_B*g_C&
          &**2*Dr_A+half*Dr_A*g_A*g_C*Dz_C)*r+half*Dr_B*g_C*g_A*Dr_A&
          &+half*Dr_B*g_B*g_A*Dr_C+half*g_B**2*Dz_B*Dz_C+half*Dz_A&
          &*g_A*g_B*Dz_C+half*g_A**2*Dr_A*Dr_C+half*Dr_B**2*g_B*g_C&
          &+g_C*g_A*C**2-half*g_C*g_A*C*Dz_A+half*Dz_A*g_B*g_C*Dz_B&
          &+half*Dz_A**2*g_A*g_C+half*Dz_B*g_B*g_C*C+(half*g_H**2*C&
          &*Dr_H+half*g_B*g_A*Dr_B*C-half*g_H**2*Dz_H*Dr_H+half*g_A&
          &**2*Dr_A*C)/r+half*Dz_H*g_H**2*lambda

     RIC_lambda = RIC_lambda + (two*g_B*g_C*Dr_C*Dz_C+g_C**2*Dr_A*Dz_C+two*Dz_A&
          &*g_C**2*Dr_C)*r+two*Dr_A*g_C*g_A*Dr_C+g_B*g_A*Dr_C**2+half&
          &*Dr_A*g_C**2*Dr_B+two*Dz_A*g_C**2*C+two*g_B*g_C*C*Dz_C&
          &-half*Dz_A**2*g_C**2+(g_B**2*Dr_B*Dz_C+Dz_A*g_C*g_A*Dr_A&
          &+half*C*g_C*g_H*Dr_H+Dz_A*g_B*g_C*Dr_B+two*g_B*g_A*Dr_C*C&
          &+two*Dr_A*g_C*g_A*C-g_H*g_C*Dz_H*Dr_H)/r+(-fourth*g_H**2&
          &*Dr_H**2-fourth*g_B**2*Dr_B**2+3.0*fourth*Dr_A**2*g_A**2&
          &-half*g_H*g_A*Dr_H**2)/r**2+half*g_B*(ft7*(g_lambda*Dz_A&
          &**2+g_H*Dz_lambda*(Dz_A+Dz_H))+(1.0-ft7)*(g_A*Dz_lambda&
          &*(Dz_A+Dz_H)+g_lambda*Dz_H**2))+C**2*g_C**2-g_H*(half*C&
          &*g_C-g_H*lambda)*Dr_H/r

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine calc_conformalRicci3

