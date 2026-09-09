
  subroutine sources_complexproca

! *******************************************
! ***   SOURCES FOR COMPLEX PROCA FIELD   ***
! *******************************************

! This routine calculates the sources for
! a complex Proca field. 

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  real(8) vl,var0
  real(8) zero,half,one,two


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  half = 0.5d0
  one  = 1.d0
  two  = 2.d0


! *******************
! ***   SOURCES   ***
! *******************

! Source for scalar potential.

  sproc_phiR = alpha * proc_phiR * trK                                           &
               - proc_AR_r*Dr_alpha - proc_AR_z*Dz_alpha                         &
               - alpha*proc_ARu_r * (8.d0*Dr_phi + (1.d0/hdet)*Dr_hdet)          &
               - alpha*proc_ARu_z * (8.d0*Dz_phi + (1.d0/hdet)*Dz_hdet)          &
               - (alpha/psi4) * ( r * ( proc_AR_r*Dz_g_C + proc_AR_z*Dr_g_C)     &
                     + proc_AR_r*Dr_g_A + proc_AR_z*g_C + proc_AR_z*Dz_g_B)      &
               - (alpha/psi4) * (g_A*Dr_proc_AR_r + r*g_C*Dz_proc_AR_r +         &
                                 g_B*Dz_proc_AR_z + r*g_C*Dr_proc_AR_z )
               
  sproc_phiI = alpha * proc_phiI * trK                                           &
               - proc_AI_r*Dr_alpha - proc_AI_z*Dz_alpha                         &
               - alpha*proc_AIu_r * (8.d0*Dr_phi + (1.d0/hdet)*Dr_hdet)          &
               - alpha*proc_AIu_z * (8.d0*Dz_phi + (1.d0/hdet)*Dz_hdet)          &
               - (alpha/psi4) * ( r * ( proc_AI_r*Dz_g_C + proc_AI_z*Dr_g_C)     &
                     + proc_AI_r*Dr_g_A + proc_AI_z*g_C + proc_AI_z*Dz_g_B)      &
               - (alpha/psi4) * (g_A*Dr_proc_AI_r + r*g_C*Dz_proc_AI_r +         &
                                 g_B*Dz_proc_AI_z + r*g_C*Dr_proc_AI_z )

!  Notice that no new terms appear in the case of angular momentum.

  if (shift/="none") then 
     sproc_phiR = sproc_phiR + zero   ! Pongo cero provisionalmente
  end if

! Source for vector potential.

  sproc_AR_r = - alpha * Dr_proc_phiR - proc_phiR * Dr_alpha      &
               - alpha * psi4 * (A*proc_ER_r + r*C*proc_ER_z)

  sproc_AI_r = - alpha * Dr_proc_phiI - proc_phiI * Dr_alpha      &
               - alpha * psi4 * (A*proc_EI_r + r*C*proc_EI_z)

  sproc_AR_z = - alpha * Dz_proc_phiR - proc_phiR * Dz_alpha      &
               - alpha * psi4 * (r*C*proc_ER_r + B*proc_ER_z)

  sproc_AI_z = - alpha * Dz_proc_phiI - proc_phiI * Dz_alpha      &
               - alpha * psi4 * (r*C*proc_EI_r + B*proc_EI_z)

  if (shift/="none") then 
     sproc_AR_r = sproc_AR_r + zero   ! Pongo cero provisionalmente
     sproc_AR_z = sproc_AR_z + zero   ! Pongo cero provisionalmente
  end if

! Source for electric field.

  sproc_ER_r = alpha*trK*proc_ER_r + alpha * proc_mass**2 * (1/ psi4) * &
                                   ( g_A*proc_AR_r + r*g_C*proc_AR_z )

  sproc_EI_r = alpha*trK*proc_EI_r + alpha * proc_mass**2 * (1/ psi4) * &
                                   ( g_A*proc_AI_r + r*g_C*proc_AI_z )

  sproc_ER_z = alpha*trK*proc_ER_z + alpha * proc_mass**2 * (1/ psi4) * &
                                   ( r*g_C*proc_AR_r + g_B*proc_AR_z )

  sproc_EI_z = alpha*trK*proc_EI_z + alpha * proc_mass**2 * (1/ psi4) * &
                                   ( r*g_C*proc_AI_r + g_B*proc_AI_z )

  if (shift/="none") then 
     sproc_ER_r = sproc_ER_r + zero   ! Pongo cero provisionalmente
     sproc_ER_z = sproc_ER_z + zero   ! Pongo cero provisionalmente
  end if

! Source for magnetic field.

  sproc_BR_r = alpha*trK*proc_BR_r
  sproc_BI_r = alpha*trK*proc_BI_r
  sproc_BR_z = alpha*trK*proc_BR_z
  sproc_BI_z = alpha*trK*proc_BI_z

  if (shift/="none") then 
     sproc_BR_r = sproc_BR_r + zero   ! Pongo cero provisionalmente
     sproc_BR_z = sproc_BR_z + zero   ! Pongo cero provisionalmente
  end if

! Dissipation.

  if (complexprocadiss/=zero) then

!    Real part.

     evolvevar => proc_phiR
     sourcevar => sproc_phiR
     call dissipation(+1,+1,complexprocadiss)

     evolvevar => proc_AR_r
     sourcevar => sproc_AR_r
     call dissipation(-1,+1,complexprocadiss)

     evolvevar => proc_AR_z
     sourcevar => sproc_AR_z
     call dissipation(+1,-1,complexprocadiss)

     evolvevar => proc_ER_r
     sourcevar => sproc_ER_r
     call dissipation(-1,+1,complexprocadiss)

     evolvevar => proc_ER_z
     sourcevar => sproc_ER_z
     call dissipation(+1,-1,complexprocadiss)

     evolvevar => proc_BR_r
     sourcevar => sproc_BR_r
     call dissipation(-1,+1,complexprocadiss)

     evolvevar => proc_BR_z
     sourcevar => sproc_BR_z
     call dissipation(+1,-1,complexprocadiss)

!    Imaginary part.

     evolvevar => proc_phiI 
     sourcevar => sproc_phiI 
     call dissipation(+1,+1,complexprocadiss)

     evolvevar => proc_AI_r
     sourcevar => sproc_AI_r
     call dissipation(-1,+1,complexprocadiss)

     evolvevar => proc_AI_z
     sourcevar => sproc_AI_z
     call dissipation(+1,-1,complexprocadiss)

     evolvevar => proc_EI_r
     sourcevar => sproc_EI_r
     call dissipation(-1,+1,complexprocadiss)

     evolvevar => proc_EI_z
     sourcevar => sproc_EI_z
     call dissipation(+1,-1,complexprocadiss)

     evolvevar => proc_BI_r
     sourcevar => sproc_BI_r
     call dissipation(-1,+1,complexprocadiss)

     evolvevar => proc_BI_z
     sourcevar => sproc_BI_z
     call dissipation(+1,-1,complexprocadiss)

  end if


! ********************************
! ***   RADIATIVE BOUNDARIES   ***
! ********************************

  vl = one
  var0 = zero

! Real part.

  evolvevar => proc_phiR
  sourcevar => sproc_phiR
  Dr_var => Dr_proc_phiR
  Dz_var => Dz_proc_phiR
  call radbound(vl,var0)

  evolvevar => proc_AR_r
  sourcevar => sproc_AR_r
  Dr_var => Dr_proc_AR_r
  Dz_var => Dz_proc_AR_r
  call radbound(vl,var0)

  evolvevar => proc_AR_z
  sourcevar => sproc_AR_z
  Dr_var => Dr_proc_AR_z
  Dz_var => Dz_proc_AR_z
  call radbound(vl,var0)

  evolvevar => proc_ER_r
  sourcevar => sproc_ER_r
  Dr_var => Dr_proc_ER_r
  Dz_var => Dz_proc_ER_r
  call radbound(vl,var0)

  evolvevar => proc_ER_z
  sourcevar => sproc_ER_z
  Dr_var => Dr_proc_ER_z
  Dz_var => Dz_proc_ER_z
  call radbound(vl,var0)

  evolvevar => proc_BR_r
  sourcevar => sproc_BR_r
  Dr_var => Dr_proc_BR_r
  Dz_var => Dz_proc_BR_r
  call radbound(vl,var0)

  evolvevar => proc_BR_z
  sourcevar => sproc_BR_z
  Dr_var => Dr_proc_BR_z
  Dz_var => Dz_proc_BR_z
  call radbound(vl,var0)

! Imaginary part.

  evolvevar => proc_phiI
  sourcevar => sproc_phiI
  Dr_var => Dr_proc_phiI
  Dz_var => Dz_proc_phiI
  call radbound(vl,var0)

  evolvevar => proc_AI_r
  sourcevar => sproc_AI_r
  Dr_var => Dr_proc_AI_r
  Dz_var => Dz_proc_AI_r
  call radbound(vl,var0)

  evolvevar => proc_AI_z
  sourcevar => sproc_AI_z
  Dr_var => Dr_proc_AI_z
  Dz_var => Dz_proc_AI_z
  call radbound(vl,var0)

  evolvevar => proc_EI_r
  sourcevar => sproc_EI_r
  Dr_var => Dr_proc_EI_r
  Dz_var => Dz_proc_EI_r
  call radbound(vl,var0)

  evolvevar => proc_EI_z
  sourcevar => sproc_EI_z
  Dr_var => Dr_proc_EI_z
  Dz_var => Dz_proc_EI_z
  call radbound(vl,var0)

  evolvevar => proc_BI_r
  sourcevar => sproc_BI_r
  Dr_var => Dr_proc_BI_r
  Dz_var => Dz_proc_BI_r
  call radbound(vl,var0)

  evolvevar => proc_BI_z
  sourcevar => sproc_BI_z
  Dr_var => Dr_proc_BI_z
  Dz_var => Dz_proc_BI_z
  call radbound(vl,var0)


! ***************
! ***   END   ***
! ***************

  end subroutine sources_complexproca
