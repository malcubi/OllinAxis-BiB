
  subroutine sources_maxwell

! ******************************************
! ***   SOURCES FOR REAK MAXWELL FIELD   ***
! ******************************************

! This routine calculates the sources for a Maxwell field.
! Here phi is the scalar potential.
!
! The Maxwell equation in 3+1 formalism has the form
!
!                                            i
! d phi  -  L  phi = - nabla (alpha * Am ) + alpha * phi * K
!  t         beta              i
!
!
! d Am  -  L    Am = - nabla (alpha * phi ) - alpha * Ep
!  t  i     beta  i         i                           i
!
!     i           i                          i                i                   i
! d Em  -  L    Em = + [nabla X (alpha * Bm)] + alpha * K * Em  - 4 pi * alpha * J
!  t       beta
!
!     i           i                           i                i
! d Bm  -  L    Bm = - [nabla X (alpha * Em )] + alpha * K * Bm
!  t       beta


! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  real(8) vl,var0
  real(8) zero,half,one,two,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  half = 0.5d0
  one  = 1.d0
  two  = 2.d0
  smallpi = acos(-one)


! *******************
! ***   SOURCES   ***
! *******************

! Source for scalar potential.

  smaxw_phi = alpha * maxw_phi * trK                                           &
               - maxw_A_r*Dr_alpha - maxw_A_z*Dz_alpha                         &
               - alpha*maxw_Au_r * (2.d0*Dr_phi + (0.5d0/hdet)*Dr_hdet)        &
               - alpha*maxw_Au_z * (2.d0*Dz_phi + (0.5d0/hdet)*Dz_hdet)        &
               - (alpha/psi4) * ( r * ( maxw_A_r*Dz_g_C + maxw_A_z*Dr_g_C)     &
                     + maxw_A_r*Dr_g_A + maxw_A_z*g_C + maxw_A_z*Dz_g_B)       &
               - (alpha/psi4) * (g_A*Dr_maxw_A_r + r*g_C*Dz_maxw_A_r +         &
                                 g_B*Dz_maxw_A_z + r*g_C*Dr_maxw_A_z )

  if (shift/="none") then 
     smaxw_phi = smaxw_phi + beta_r*DAr_maxw_phi + beta_z*DAz_maxw_phi
  end if

! Source for vector potential.

  smaxw_A_r = - alpha * Dr_maxw_phi - maxw_phi * Dr_alpha - alpha * maxw_Ed_r
  smaxw_A_z = - alpha * Dz_maxw_phi - maxw_phi * Dz_alpha - alpha * maxw_Ed_z

  if (angmom) then
     smaxw_A_p = - alpha * maxw_Ed_p
  end if

  if (shift/="none") then 

     smaxw_A_r = smaxw_A_r + beta_r*DAr_maxw_A_r + beta_z*DAz_maxw_A_r    &
                           + maxw_A_r*Dr_beta_r  + maxw_A_z*Dr_beta_z
     smaxw_A_z = smaxw_A_z + beta_r*DAr_maxw_A_z + beta_z*DAz_maxw_A_z    &
                           + maxw_A_r*Dz_beta_r  + maxw_A_z*Dz_beta_z

     if (angmom) then
        smaxw_A_p = smaxw_A_p + beta_r*DAr_maxw_A_p + beta_z*DAz_maxw_A_p
     end if

  end if

! Source for electric field.

  smaxw_E_r = alpha*trK*maxw_E_r - 4.d0*smallpi*alpha*J_r
  smaxw_E_z = alpha*trK*maxw_E_z - 4.d0*smallpi*alpha*J_z

  if (angmom) then

     smaxw_E_r = smaxw_E_r - maxw_Bd_p * Dz_alpha + alpha *           &
                 ( r*C*CovDp_maxw_B_r     + B*CovDp_maxw_B_z          &
                 + r**2*C2*CovDp_maxw_B_p - r**3*C1*CovDz_maxw_B_r    & ! Warning. Remember that CovDp_maxw_B_p contents singular christoffel
                 - r**2*C2*CovDz_maxw_B_z - r**2*H*CovDz_maxw_B_p )

     smaxw_E_z = smaxw_E_z + maxw_Bd_p * Dr_alpha + alpha *           &
                 ( r**3*C1*CovDr_maxw_B_r + r**2*C2*CovDr_maxw_B_z    &
                 + r**2*H*CovDr_maxw_B_p  - A*CovDp_maxw_B_r          & ! Warning. Remember that CovDr_maxw_B_p contents singular christoffel
                 - r*C*CovDp_maxw_B_z     - r**3*C1*CovDp_maxw_B_p )    ! Warning. Remember that CovDp_maxw_B_p contents singular christoffel

     smaxw_E_p = alpha*trK*maxw_E_p                                   &
                 + maxw_Bd_r*Dz_alpha - maxw_Bd_z*Dr_alpha + alpha*   &
                 ( A*CovDz_maxw_B_r       + r*C*CovDz_maxw_B_z        &
                 + r**3*C1*CovDz_maxw_B_p - r*C*CovDr_maxw_B_r        &
                 - B*CovDr_maxw_B_z       - r**2*C2*CovDr_maxw_B_p )    ! Warning. Remember that CovDr_maxw_B_p contents singular christoffel

  end if

  if (shift/="none") then 

     smaxw_E_r = smaxw_E_r + beta_r*DAr_maxw_E_r + beta_z*DAz_maxw_E_r    &
                           - maxw_E_r*Dr_beta_r  - maxw_E_z*Dz_beta_r

     smaxw_E_z = smaxw_E_z + beta_r*DAr_maxw_E_z + beta_z*DAz_maxw_E_z    &
                           - maxw_E_r*Dr_beta_z  - maxw_E_z*Dz_beta_z

     if (angmom) then
        smaxw_E_p = smaxw_E_p + beta_r*DAr_maxw_E_p + beta_z*DAz_maxw_E_p &
                              - maxw_E_r*Dr_beta_p  - maxw_E_z*Dz_beta_p
     end if

  end if

! Source for magnetic field.

  smaxw_B_r = alpha*trK*maxw_B_r
  smaxw_B_z = alpha*trK*maxw_B_z

  if (angmom) then

     smaxw_B_r = smaxw_B_r + maxw_Ed_p * Dz_alpha - alpha *           &
                 ( r*C*CovDp_maxw_E_r     + B*CovDp_maxw_E_z          &
                 + r**2*C2*CovDp_maxw_E_p - r**3*C1*CovDz_maxw_E_r    & ! Warning. Remember that CovDp_maxw_E_p contents  singular christoffel
                 - r**2*C2*CovDz_maxw_E_z - r**2*H*CovDz_maxw_E_p )

     smaxw_B_z = smaxw_B_z - maxw_Ed_p * Dr_alpha - alpha *           &
                 ( r**3*C1*CovDr_maxw_E_r + r**2*C2*CovDr_maxw_E_z    &
                 + r**2*H*CovDr_maxw_E_p  - A*CovDp_maxw_E_r          & ! Warning. Remember that CovDr_maxw_E_p contents  singular christoffel
                 - r*C*CovDp_maxw_E_z     - r**3*C1*CovDp_maxw_E_p )    ! Warning. Remember that CovDp_maxw_E_p contents  singular christoffel

     smaxw_B_p = alpha*trK*maxw_B_p                                   &
                 - maxw_Ed_r*Dz_alpha + maxw_Ed_z*Dr_alpha - alpha*   &
                 ( A*CovDz_maxw_E_r       + r*C*CovDz_maxw_E_z        &
                 + r**3*C1*CovDz_maxw_E_p - r*C*CovDr_maxw_E_r        &
                 - B*CovDr_maxw_E_z       - r**2*C2*CovDr_maxw_E_p )    ! Warning. Remember that CovDr_maxw_E_p contents  singular christoffel

  end if

  if (shift/="none") then

     smaxw_B_r = smaxw_B_r + beta_r*DAr_maxw_B_r + beta_z*DAz_maxw_B_r    &
                           - maxw_B_r*Dr_beta_r  - maxw_B_z*Dz_beta_r

     smaxw_B_z = smaxw_B_z + beta_r*DAr_maxw_B_z + beta_z*DAz_maxw_B_z    &
                           - maxw_B_r*Dr_beta_z  - maxw_B_z*Dz_beta_z

     if (angmom) then
        smaxw_B_p = smaxw_B_p + beta_r*DAr_maxw_B_p + beta_z*DAz_maxw_B_p &
                              - maxw_B_r*Dr_beta_p  - maxw_B_z*Dz_beta_p
     end if

  end if

! Dissipation.

  if (maxwelldiss/=zero) then

!    Real part.

     evolvevar => maxw_phi
     sourcevar => smaxw_phi
     call dissipation(+1,+1,maxwelldiss)

     evolvevar => maxw_A_r
     sourcevar => smaxw_A_r
     call dissipation(-1,+1,maxwelldiss)

     evolvevar => maxw_A_z
     sourcevar => smaxw_A_z
     call dissipation(+1,-1,maxwelldiss)

     evolvevar => maxw_E_r
     sourcevar => smaxw_E_r
     call dissipation(-1,+1,maxwelldiss)

     evolvevar => maxw_E_z
     sourcevar => smaxw_E_z
     call dissipation(+1,-1,maxwelldiss)

     evolvevar => maxw_B_r
     sourcevar => smaxw_B_r
     call dissipation(-1,+1,maxwelldiss)

     evolvevar => maxw_B_z
     sourcevar => smaxw_B_z
     call dissipation(+1,-1,maxwelldiss)

  end if


! ********************************
! ***   RADIATIVE BOUNDARIES   ***
! ********************************

  vl = one
  var0 = zero

  evolvevar => maxw_phi
  sourcevar => smaxw_phi
  Dr_var => Dr_maxw_phi
  Dz_var => Dz_maxw_phi
  call radbound(vl,var0)

  evolvevar => maxw_A_r
  sourcevar => smaxw_A_r
  Dr_var => Dr_maxw_A_r
  Dz_var => Dz_maxw_A_r
  call radbound(vl,var0)

  evolvevar => maxw_A_z
  sourcevar => smaxw_A_z
  Dr_var => Dr_maxw_A_z
  Dz_var => Dz_maxw_A_z
  call radbound(vl,var0)

  evolvevar => maxw_E_r
  sourcevar => smaxw_E_r
  Dr_var => Dr_maxw_E_r
  Dz_var => Dz_maxw_E_r
  call radbound(vl,var0)

  evolvevar => maxw_E_z
  sourcevar => smaxw_E_z
  Dr_var => Dr_maxw_E_z
  Dz_var => Dz_maxw_E_z
  call radbound(vl,var0)

  evolvevar => maxw_B_r
  sourcevar => smaxw_B_r
  Dr_var => Dr_maxw_B_r
  Dz_var => Dz_maxw_B_r
  call radbound(vl,var0)

  evolvevar => maxw_B_z
  sourcevar => smaxw_B_z
  Dr_var => Dr_maxw_B_z
  Dz_var => Dz_maxw_B_z
  call radbound(vl,var0)


! ***************
! ***   END   ***
! ***************

  end subroutine sources_maxwell
