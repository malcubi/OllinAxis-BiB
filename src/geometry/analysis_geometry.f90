!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/analysis_geometry.f90,v 1.17 2020/12/05 00:51:19 malcubi Exp $

  subroutine analysis_geometry(box,level)

! ********************************************************
! ***   CALCULATION OF ANALYSIS VARIABLES FOR OUTPUT   ***
! ********************************************************

! Include modules.

  use arrays
  use param
  use derivatives

! Extra variables.

  implicit none

  integer box,level

  real(8) third,half,one,two,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  third = 1.d0/3.d0
  half = 0.5d0

  one = 1.d0
  two = 2.d0
  smallpi = acos(-one)


! ************************************
! ***   PSEUDO SCWARZSCHILD MASS   ***
! ************************************

! The pseudo-Schwarzschild mass is obtained by noticing
! that for the Schwarzschild metric we have:
!
! g   =  1 / (1 - 2M/R)
!  RR
!
! with R the Schwarzschild or areal radius. This implies:
!
!                                             2         4
! M  =  R/2 ( 1 - 1/g  )  =  R/2 ( 1 - (dR/dr) / grr psi )
!                    RR
!
! with r the coordinate radius.  Notice also that in general:
!
!
!  d   =  ( rho d    +  z d  ) / r
!   r            rho       z
!
! But remember that the code uses r for "rho" and "rr" for r.
!
! We now take the areal radius as:
!
! R  =  psi**2 sqrt(gtt)
!
! where gtt is the angular conformal metric: gtt := g_(theta,theta).
! (The code calculates both grr and gtt in auxiliary_geometry.f90)
!
! Note that taking this value of R we are assuming spherical symmetry,
! which will only be true sufficiently far away.

  if (associated(grid(box,level)%mass_sch)) then

!    Find areal radius and save it in "auxarray".

     auxarray = psi2*sqrt(gtt)
     !auxarray = psi2*sqrt(abs(rr/r))*(gtt*gpp)**0.25d0 ! This is another possible way of defining R.

!    Find derivatives of areal radius.

     diffvar => auxarray
     Dr_auxarray = diff1r(+1)
     Dz_auxarray = diff1z(+1)

!    Calculate mass.

     mass_sch = half*auxarray*(one - ((r*Dr_auxarray + z*Dz_auxarray)/rr)**2/(grr*psi4))

  end if


! ********************
! ***   ADM MASS   ***
! ********************

  if (associated(grid(box,level)%mass_ADM)) then

  end if


! *********************************
! ***   CHARACTERISTIC SPEEDS   ***
! *********************************

! Light speed.

  vl_rp = + abs(alpha)*sqrt(abs(g_A/psi4))
  vl_rm = - abs(alpha)*sqrt(abs(g_A/psi4))

  vl_zp = + abs(alpha)*sqrt(abs(g_B/psi4))
  vl_zm = - abs(alpha)*sqrt(abs(g_B/psi4))

! Slicing speeds.

  if (slicing=='maximal') then

     va_rp = vl_rp
     va_rm = vl_rm

     va_zp = vl_zp
     va_zm = vl_zm

  else if ((index(slicing,"harmonic")/=0).or. &
           (index(slicing,"1+log")/=0).or. &
           (index(slicing,"shockavoid")/=0)) then

     va_rp = + sqrt(abs(falpha*g_A/psi4))
     va_rm = - sqrt(abs(falpha*g_A/psi4))

     va_zp = + sqrt(abs(falpha*g_B/psi4))
     va_zm = - sqrt(abs(falpha*g_B/psi4))

  end if

! Shift speeds.

  if ((shift(1:11)=="Gammadriver").and.(drivercsi/=0.d0)) then

     vs_rp = + sqrt(abs(4.d0/3.d0*drivercsi*g_A/psi4))
     vs_rm = - sqrt(abs(4.d0/3.d0*drivercsi*g_A/psi4))

     vs_zp = + sqrt(abs(4.d0/3.d0*drivercsi*g_B/psi4))
     vs_zm = - sqrt(abs(4.d0/3.d0*drivercsi*g_B/psi4))

  end if

! Add shift contribution.

  if (shift/="none") then

     vl_rp = - beta_r + vl_rp
     vl_rm = - beta_r + vl_rm
     vl_zp = - beta_z + vl_zp
     vl_zm = - beta_z + vl_zm

     va_rp = - beta_r + va_rp
     va_rm = - beta_r + va_rm
     va_zp = - beta_z + va_zp
     va_zm = - beta_z + va_zm

     vs_rp = - beta_r + vs_rp
     vs_rm = - beta_r + vs_rm
     vs_zp = - beta_z + vs_zp
     vs_zm = - beta_z + vs_zm

  end if


! ***********************
! ***   WEYL TENSOR   ***
! ***********************

  if (curvInv) then
     call weyl
  end if


! ***************************
! ***   WAVE EXTRACTION   ***
! ***************************


! ***************
! ***   END   ***
! ***************

  end subroutine analysis_geometry

