
  subroutine sources_complexproca

! ********************************************
! ***   SOURCES FOR COMPLEX SCALAR FIELD   ***
! ********************************************

! This routine calculates the sources for a complex Proca field.
! Here varphi is the scalar potential.
!
! The Proca equation in 3+1 formalism has the form:
!
!                                             i
! d varphi  -  L  varphi = - nabla (alpha * Ap ) + alpha * varphi * K
!  t            beta              i
!
!
! d x  -  L    x = - nabla (alpha * varphi) - alpha * Ep
!  t i     beta i         i                             i
!
!     i           i                          i                i              2    i
! d Ep  -  L    Ep = + [nabla X (alpha * Bp)] + alpha * K * Ep + alpha * mass * Ap
!  t       beta

!     i           i                          i                i
! d Bp  -  L    Bp = - [nabla X (alpha * Ep)] + alpha * K * Bp
!  t       beta

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
               - alpha*proc_ARu_r * (2.d0*Dr_phi + (0.5d0/hdet)*Dr_hdet)         &
               - alpha*proc_ARu_z * (2.d0*Dz_phi + (0.5d0/hdet)*Dz_hdet)         &
               - (alpha/psi4) * ( r * ( proc_AR_r*Dz_g_C + proc_AR_z*Dr_g_C)     &
                     + proc_AR_r*Dr_g_A + proc_AR_z*g_C + proc_AR_z*Dz_g_B)      &
               - (alpha/psi4) * (g_A*Dr_proc_AR_r + r*g_C*Dz_proc_AR_r +         &
                                 g_B*Dz_proc_AR_z + r*g_C*Dr_proc_AR_z )

  sproc_phiI = alpha * proc_phiI * trK                                           &
               - proc_AI_r*Dr_alpha - proc_AI_z*Dz_alpha                         &
               - alpha*proc_AIu_r * (2.d0*Dr_phi + (0.5d0/hdet)*Dr_hdet)         &
               - alpha*proc_AIu_z * (2.d0*Dz_phi + (0.5d0/hdet)*Dz_hdet)         &
               - (alpha/psi4) * ( r * ( proc_AI_r*Dz_g_C + proc_AI_z*Dr_g_C)     &
                     + proc_AI_r*Dr_g_A + proc_AI_z*g_C + proc_AI_z*Dz_g_B)      &
               - (alpha/psi4) * (g_A*Dr_proc_AI_r + r*g_C*Dz_proc_AI_r +         &
                                 g_B*Dz_proc_AI_z + r*g_C*Dr_proc_AI_z )

!  Notice that adding angular momentum does not introduce any additional terms.

  if (shift/="none") then 
     sproc_phiR = sproc_phiR + beta_r*DAr_proc_phiR + beta_z*DAz_proc_phiR
     sproc_phiI = sproc_phiI + beta_r*DAr_proc_phiI + beta_z*DAz_proc_phiI
  end if

! Source for vector potential.

  sproc_AR_r = - alpha * Dr_proc_phiR - proc_phiR * Dr_alpha - alpha * proc_ERd_r
  sproc_AI_r = - alpha * Dr_proc_phiI - proc_phiI * Dr_alpha - alpha * proc_EId_r
  sproc_AR_z = - alpha * Dz_proc_phiR - proc_phiR * Dz_alpha - alpha * proc_ERd_z
  sproc_AI_z = - alpha * Dz_proc_phiI - proc_phiI * Dz_alpha - alpha * proc_EId_z

  if (angmom) then
     sproc_AR_p = - alpha * proc_ERd_p
     sproc_AI_p = - alpha * proc_EId_p
  end if

  if (shift/="none") then

     sproc_AR_r = sproc_AR_r + beta_r*DAr_proc_AR_r + beta_z*DAz_proc_AR_r    &
                             + proc_AR_r*Dr_beta_r  + proc_AR_z*Dr_beta_z
     sproc_AI_r = sproc_AI_r + beta_r*DAr_proc_AI_r + beta_z*DAz_proc_AI_r    &
                             + proc_AI_r*Dr_beta_r  + proc_AI_z*Dr_beta_z
     sproc_AR_z = sproc_AR_z + beta_r*DAr_proc_AR_z + beta_z*DAz_proc_AR_z    &
                             + proc_AR_r*Dz_beta_r  + proc_AR_z*Dz_beta_z
     sproc_AI_z = sproc_AI_z + beta_r*DAr_proc_AI_z + beta_z*DAz_proc_AI_z    &
                             + proc_AI_r*Dz_beta_r  + proc_AI_z*Dz_beta_z

     if (angmom) then
        sproc_AR_p = sproc_AR_p + beta_r*DAr_proc_AR_p + beta_z*DAz_proc_AR_p
        sproc_AI_p = sproc_AI_p + beta_r*DAr_proc_AI_p + beta_z*DAz_proc_AI_p
     end if

  end if

! Source for electric field.

  sproc_ER_r = alpha * trK * proc_ER_r + alpha * proc_mass**2 * proc_ARu_r
  sproc_EI_r = alpha * trK * proc_EI_r + alpha * proc_mass**2 * proc_AIu_r
  sproc_ER_z = alpha * trK * proc_ER_z + alpha * proc_mass**2 * proc_ARu_z
  sproc_EI_z = alpha * trK * proc_EI_z + alpha * proc_mass**2 * proc_AIu_z

  if (angmom) then

     sproc_ER_r = sproc_ER_r - proc_BRd_p * Dz_alpha + alpha *           &
                  ( r*C*CovDp_proc_BR_r     + B*CovDp_proc_BR_z          &
                  + r**2*C2*CovDp_proc_BR_p - r**3*C1*CovDz_proc_BR_r    & ! Warning. Remember that CovDp_proc_BR_p contents singular christoffel
                  - r**2*C2*CovDz_proc_BR_z - r**2*H*CovDz_proc_BR_p )

     sproc_EI_r = sproc_EI_r - proc_BId_p * Dz_alpha + alpha *           &
                  ( r*C*CovDp_proc_BI_r     + B*CovDp_proc_BI_z          &
                  + r**2*C2*CovDp_proc_BI_p - r**3*C1*CovDz_proc_BI_r    & ! Warning. Remember that CovDp_proc_BI_p contents singular christoffel
                  - r**2*C2*CovDz_proc_BI_z - r**2*H*CovDz_proc_BI_p )

     sproc_ER_z = sproc_ER_z + proc_BRd_p * Dr_alpha + alpha *           &
                  ( r**3*C1*CovDr_proc_BR_r + r**2*C2*CovDr_proc_BR_z    &
                  + r**2*H*CovDr_proc_BR_p  - A*CovDp_proc_BR_r          & ! Warning. Remember that CovDr_proc_BR_p contents singular christoffel
                  - r*C*CovDp_proc_BR_z     - r**3*C1*CovDp_proc_BR_p )    ! Warning. Remember that CovDp_proc_BR_p contents singular christoffel

     sproc_EI_z = sproc_EI_z + proc_BId_p * Dr_alpha + alpha *           &
                  ( r**3*C1*CovDr_proc_BI_r + r**2*C2*CovDr_proc_BI_z    &
                  + r**2*H*CovDr_proc_BI_p  - A*CovDp_proc_BI_r          & ! Warning. Remember that CovDr_proc_BI_p contents singular christoffel
                  - r*C*CovDp_proc_BI_z     - r**3*C1*CovDp_proc_BI_p )    ! Warning. Remember that CovDp_proc_BI_p contents singular christoffel

     sproc_ER_p = alpha*trK*proc_ER_p + alpha*proc_mass**2*proc_ARu_p    &
                  + proc_BRd_r*Dz_alpha - proc_BRd_z*Dr_alpha + alpha*   &
                  ( A*CovDz_proc_BR_r       + r*C*CovDz_proc_BR_z        &
                  + r**3*C1*CovDz_proc_BR_p - r*C*CovDr_proc_BR_r        &
                  - B*CovDr_proc_BR_z       - r**2*C2*CovDr_proc_BR_p )    ! Warning. Remember that CovDr_proc_BR_p contents singular christoffel

     sproc_EI_p = alpha*trK*proc_EI_p + alpha*proc_mass**2*proc_AIu_p    &
                  + proc_BId_r*Dz_alpha - proc_BId_z*Dr_alpha + alpha*   &
                  ( A*CovDz_proc_BI_r       + r*C*CovDz_proc_BI_z        &
                  + r**3*C1*CovDz_proc_BI_p - r*C*CovDr_proc_BI_r        &
                  - B*CovDr_proc_BI_z       - r**2*C2*CovDr_proc_BI_p )    ! Warning. Remember that CovDr_proc_BI_p contents singular christoffel

  end if

  if (shift/="none") then

     sproc_ER_r = sproc_ER_r + beta_r*DAr_proc_ER_r + beta_z*DAz_proc_ER_r    &
                             - proc_ER_r*Dr_beta_r  - proc_ER_z*Dz_beta_r

     sproc_EI_r = sproc_EI_r + beta_r*DAr_proc_EI_r + beta_z*DAz_proc_EI_r    &
                             - proc_EI_r*Dr_beta_r  - proc_EI_z*Dz_beta_r

     sproc_ER_z = sproc_ER_z + beta_r*DAr_proc_ER_z + beta_z*DAz_proc_ER_z    &
                             - proc_ER_r*Dr_beta_z  - proc_ER_z*Dz_beta_z

     sproc_EI_z = sproc_EI_z + beta_r*DAr_proc_EI_z + beta_z*DAz_proc_EI_z    &
                             - proc_EI_r*Dr_beta_z  - proc_EI_z*Dz_beta_z

     if (angmom) then
        sproc_ER_p = sproc_ER_p + beta_r*DAr_proc_ER_p + beta_z*DAz_proc_ER_p    &
                                - proc_ER_r*Dr_beta_p  - proc_ER_z*Dz_beta_p
        sproc_EI_p = sproc_EI_p + beta_r*DAr_proc_EI_p + beta_z*DAz_proc_EI_p    &
                                - proc_EI_r*Dr_beta_p  - proc_EI_z*Dz_beta_p
     end if

  end if

! Source for magnetic field.

  sproc_BR_r = alpha * trK * proc_BR_r
  sproc_BI_r = alpha * trK * proc_BI_r
  sproc_BR_z = alpha * trK * proc_BR_z
  sproc_BI_z = alpha * trK * proc_BI_z

  if (angmom) then

     sproc_BR_r = sproc_BR_r + proc_ERd_p * Dz_alpha - alpha *           &
                  ( r*C*CovDp_proc_ER_r     + B*CovDp_proc_ER_z          &
                  + r**2*C2*CovDp_proc_ER_p - r**3*C1*CovDz_proc_ER_r    & ! Warning. Remember that CovDp_proc_ER_p contents  singular christoffel
                  - r**2*C2*CovDz_proc_ER_z - r**2*H*CovDz_proc_ER_p )

     sproc_BI_r = sproc_BI_r + proc_EId_p * Dz_alpha - alpha *           &
                  ( r*C*CovDp_proc_EI_r     + B*CovDp_proc_EI_z          &
                  + r**2*C2*CovDp_proc_EI_p - r**3*C1*CovDz_proc_EI_r    & ! Warning. Remember that CovDp_proc_EI_p contents  singular christoffel
                  - r**2*C2*CovDz_proc_EI_z - r**2*H*CovDz_proc_EI_p )

     sproc_BR_z = sproc_BR_z - proc_ERd_p * Dr_alpha - alpha *           &
                  ( r**3*C1*CovDr_proc_ER_r + r**2*C2*CovDr_proc_ER_z    &
                  + r**2*H*CovDr_proc_ER_p  - A*CovDp_proc_ER_r          & ! Warning. Remember that CovDr_proc_ER_p contents  singular christoffel
                  - r*C*CovDp_proc_ER_z     - r**3*C1*CovDp_proc_ER_p )    ! Warning. Remember that CovDp_proc_ER_p contents  singular christoffel

     sproc_BI_z = sproc_BI_z - proc_EId_p * Dr_alpha - alpha *           &
                  ( r**3*C1*CovDr_proc_EI_r + r**2*C2*CovDr_proc_EI_z    &
                  + r**2*H*CovDr_proc_EI_p  - A*CovDp_proc_EI_r          & ! Warning. Remember that CovDr_proc_EI_p contents  singular christoffel
                  - r*C*CovDp_proc_EI_z     - r**3*C1*CovDp_proc_EI_p )    ! Warning. Remember that CovDp_proc_EI_p contents  singular christoffel

     sproc_BR_p = alpha*trK*proc_BR_p                                    &
                  - proc_ERd_r*Dz_alpha + proc_ERd_z*Dr_alpha - alpha*   &
                  ( A*CovDz_proc_ER_r       + r*C*CovDz_proc_ER_z        &
                  + r**3*C1*CovDz_proc_ER_p - r*C*CovDr_proc_ER_r        &
                  - B*CovDr_proc_ER_z       - r**2*C2*CovDr_proc_ER_p )    ! Warning. Remember that CovDr_proc_ER_p contents  singular christoffel

     sproc_BI_p = alpha*trK*proc_BI_p                                    &
                  - proc_EId_r*Dz_alpha + proc_EId_z*Dr_alpha - alpha*   &
                  ( A*CovDz_proc_EI_r       + r*C*CovDz_proc_EI_z        &
                  + r**3*C1*CovDz_proc_EI_p - r*C*CovDr_proc_EI_r        &
                  - B*CovDr_proc_EI_z       - r**2*C2*CovDr_proc_EI_p )    ! Warning. Remember that CovDr_proc_EI_p contents  singular christoffel

  end if

  if (shift/="none") then

     sproc_BR_r = sproc_BR_r + beta_r*DAr_proc_BR_r + beta_z*DAz_proc_BR_r    &
                             - proc_BR_r*Dr_beta_r  - proc_BR_z*Dz_beta_r

     sproc_BI_r = sproc_BI_r + beta_r*DAr_proc_BI_r + beta_z*DAz_proc_BI_r    &
                             - proc_BI_r*Dr_beta_r  - proc_BI_z*Dz_beta_r

     sproc_BR_z = sproc_BR_z + beta_r*DAr_proc_BR_z + beta_z*DAz_proc_BR_z    &
                             - proc_BR_r*Dr_beta_z  - proc_BR_z*Dz_beta_z

     sproc_BI_z = sproc_BI_z + beta_r*DAr_proc_BI_z + beta_z*DAz_proc_BI_z    &
                             - proc_BI_r*Dr_beta_z  - proc_BI_z*Dz_beta_z

     if (angmom) then
        sproc_BR_p = sproc_BR_p + beta_r*DAr_proc_BR_p + beta_z*DAz_proc_BR_p    &
                                - proc_BR_r*Dr_beta_p  - proc_BR_z*Dz_beta_p
        sproc_BI_p = sproc_BI_p + beta_r*DAr_proc_BI_p + beta_z*DAz_proc_BI_p    &
                                - proc_BI_r*Dr_beta_p  - proc_BI_z*Dz_beta_p
     end if

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

! Radiative boundaries are only applied to Pi.

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
