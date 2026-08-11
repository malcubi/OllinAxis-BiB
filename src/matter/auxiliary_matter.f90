
  subroutine auxiliary_matter

! ******************************************
! ***   AUXILIARY VARIABLES FOR MATTER   ***
! ******************************************

! Auxiliary quantities for matter: derivatives, etc.

! Include modules.

  use param
  use arrays
  use derivatives
  use derivadvect

! Extra variables.

  implicit none

  logical contains


! *****************************
! ***   REAL SCALAR FIELD   ***
! *****************************

  if (contains(mattertype,"scalar")) then

!    Derivatives of phi.

     diffvar => scalar_phi
     Dr_scalar_phi = diff1r(+1)
     Dz_scalar_phi = diff1z(+1)
     Drr_scalar_phi = diff2r(+1)
     Dzz_scalar_phi = diff2z(+1)
     Drz_scalar_phi = diff2rz(+1,+1)

!    Derivatives of pi.

     diffvar => scalar_pi
     Dr_scalar_pi = diff1r(+1)
     Dz_scalar_pi = diff1z(+1)

!    Derivatives of xi.

     diffvar => scalar_xi_r
     Dr_scalar_xi_r = diff1r(-1)
     Dz_scalar_xi_r = diff1z(+1)

     diffvar => scalar_xi_z
     Dr_scalar_xi_z = diff1r(+1)
     Dz_scalar_xi_z = diff1z(-1)

!    Advective derivatives.

     if (shift/="none") then

!       Derivatives of pi.

        diffvar => scalar_pi
        DAr_scalar_pi = diffadvr(+1)
        DAz_scalar_pi = diffadvz(+1)

!       Derivatives of xi.

        diffvar => scalar_xi_r
        DAr_scalar_xi_r = diffadvr(-1)
        DAz_scalar_xi_r = diffadvz(+1)

        diffvar => scalar_xi_z
        DAr_scalar_xi_z = diffadvr(+1)
        DAz_scalar_xi_z = diffadvz(-1)

     end if

  end if


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

  if (contains(mattertype,"complex")) then

!    Derivatives of phi.

     diffvar => complex_phiR
     Dr_complex_phiR = diff1r(+1)
     Dz_complex_phiR = diff1z(+1)
     Drr_complex_phiR = diff2r(+1)
     Dzz_complex_phiR = diff2z(+1)
     Drz_complex_phiR = diff2rz(+1,+1)

     diffvar => complex_phiI
     Dr_complex_phiI = diff1r(+1)
     Dz_complex_phiI = diff1z(+1)
     Drr_complex_phiI = diff2r(+1)
     Dzz_complex_phiI = diff2z(+1)
     Drz_complex_phiI = diff2rz(+1,+1)

!    Derivatives of pi.

     diffvar => complex_piR
     Dr_complex_piR = diff1r(+1)
     Dz_complex_piR = diff1z(+1)

     diffvar => complex_piI
     Dr_complex_piI = diff1r(+1)
     Dz_complex_piI = diff1z(+1)

!    Derivatives of xi.

     diffvar => complex_xiR_r
     Dr_complex_xiR_r = diff1r(-1)
     Dz_complex_xiR_r = diff1z(+1)

     diffvar => complex_xiR_z
     Dr_complex_xiR_z = diff1r(-1)
     Dz_complex_xiR_z = diff1z(+1)

     diffvar => complex_xiI_r
     Dr_complex_xiI_r = diff1r(-1)
     Dz_complex_xiI_r = diff1z(+1)

     diffvar => complex_xiI_z
     Dr_complex_xiI_z = diff1r(-1)
     Dz_complex_xiI_z = diff1z(+1)

!    Advective derivatives.

     if (shift/="none") then

!       Derivatives of pi.

        diffvar => complex_piR
        DAr_complex_piR = diffadvr(+1)
        DAz_complex_piR = diffadvz(+1)

        diffvar => complex_piI
        DAr_complex_piI = diffadvr(+1)
        DAz_complex_piI = diffadvz(+1)

!       Derivatives of xi.

        diffvar => complex_xiR_r
        DAr_complex_xiR_r = diffadvr(-1)
        DAz_complex_xiR_r = diffadvz(+1)

        diffvar => complex_xiR_z
        DAr_complex_xiR_z = diffadvr(+1)
        DAz_complex_xiR_z = diffadvz(-1)

        diffvar => complex_xiI_r
        DAr_complex_xiI_r = diffadvr(-1)
        DAz_complex_xiI_r = diffadvz(+1)

        diffvar => complex_xiI_z
        DAr_complex_xiI_z = diffadvr(+1)
        DAz_complex_xiI_z = diffadvz(-1)

     end if

  end if


! *************************
! ***   MAXWELL FIELD   ***
! *************************

  if (contains(mattertype,"electric")) then

!    Derivatives of phi.

     diffvar => maxw_phi
     Dr_maxw_phi = diff1r(+1)
     Dz_maxw_phi = diff1z(+1)

!    Drr_maxw_phi = diff2r(+1)
!    Dzz_maxw_phi = diff2z(+1)
!    Drz_maxw_phi = diff2rz(+1,+1)

!    Derivatives of Ai.

     diffvar => maxw_A_r
     Dr_maxw_A_r = diff1r(-1)
     Dz_maxw_A_r = diff1z(+1)

     diffvar => maxw_A_z
     Dr_maxw_A_z = diff1r(+1)
     Dz_maxw_A_z = diff1z(-1)

!    Derivatives of Ei.

     diffvar => maxw_E_r
     Dr_maxw_E_r = diff1r(-1)
     Dz_maxw_E_r = diff1z(+1)

     diffvar => maxw_E_z
     Dr_maxw_E_z = diff1r(+1)
     Dz_maxw_E_z = diff1z(-1)

!    Derivatives of Bi.

     diffvar => maxw_B_r
     Dr_maxw_B_r = diff1r(-1)
     Dz_maxw_B_r = diff1z(+1)

     diffvar => maxw_B_z
     Dr_maxw_B_z = diff1r(+1)
     Dz_maxw_B_z = diff1z(-1)

!    Advective derivatives.

     if (shift/="none") then

!       Derivatives of phi.

        diffvar => maxw_phi
        DAr_maxw_phi = diffadvr(+1)
        DAz_maxw_phi = diffadvz(+1)

!       Derivatives of A.

        diffvar => maxw_A_r
        DAr_maxw_A_r = diffadvr(-1)
        DAz_maxw_A_r = diffadvz(+1)

        diffvar => maxw_A_z
        DAr_maxw_A_z = diffadvr(+1)
        DAz_maxw_A_z = diffadvz(-1)

!       Derivatives of E.

        diffvar => maxw_E_r
        DAr_maxw_E_r = diffadvr(-1)
        DAz_maxw_E_r = diffadvz(+1)

        diffvar => maxw_E_z
        DAr_maxw_E_z = diffadvr(+1)
        DAz_maxw_E_z = diffadvz(-1)

!       Derivatives of B.

        diffvar => maxw_B_r
        DAr_maxw_B_r = diffadvr(-1)
        DAz_maxw_B_r = diffadvz(+1)

        diffvar => maxw_B_z
        DAr_maxw_B_z = diffadvr(+1)
        DAz_maxw_B_z = diffadvz(-1)

     end if

!    Covariant derivatives (without angular momentum at the moment).

     CovDr_maxw_A_r = Dr_maxw_A_r - chris_rrr*maxw_A_r - chris_zrr*maxw_A_z
     CovDz_maxw_A_r = Dz_maxw_A_r - chris_rzr*maxw_A_r - chris_zzr*maxw_A_z
     CovDr_maxw_A_z = Dr_maxw_A_z - chris_rrz*maxw_A_r - chris_zrz*maxw_A_z
     CovDz_maxw_A_z = Dz_maxw_A_z - chris_rzz*maxw_A_r - chris_zzz*maxw_A_z

     CovDr_maxw_E_r = Dr_maxw_E_r + chris_rrr*maxw_E_r + chris_rrz*maxw_E_z
     CovDr_maxw_E_z = Dr_maxw_E_z + chris_zrr*maxw_E_r + chris_zrz*maxw_E_z
     CovDz_maxw_E_r = Dz_maxw_E_r + chris_rzr*maxw_E_r + chris_rzz*maxw_E_z
     CovDz_maxw_E_z = Dz_maxw_E_z + chris_zzr*maxw_E_r + chris_zzz*maxw_E_z

     CovDr_maxw_B_r = Dr_maxw_B_r + chris_rrr*maxw_B_r + chris_rrz*maxw_B_z
     CovDr_maxw_B_z = Dr_maxw_B_z + chris_zrr*maxw_B_r + chris_zrz*maxw_B_z
     CovDz_maxw_B_r = Dz_maxw_B_r + chris_rzr*maxw_B_r + chris_rzz*maxw_B_z
     CovDz_maxw_B_z = Dz_maxw_B_z + chris_zzr*maxw_B_r + chris_zzz*maxw_B_z

  end if


! *******************************
! ***   COMPLEX PROCA FIELD   ***
! *******************************

  if (contains(mattertype,"complexproca")) then

!    Derivatives of phi.

     diffvar => proc_phiR
     Dr_proc_phiR = diff1r(+1)
     Dz_proc_phiR = diff1z(+1)

     diffvar => proc_phiI
     Dr_proc_phiI = diff1r(+1)
     Dz_proc_phiI = diff1z(+1)

!    Derivatives of A.

     diffvar => proc_AR_r
     Dr_proc_AR_r = diff1r(-1)
     Dz_proc_AR_r = diff1z(+1)

     diffvar => proc_AR_z
     Dr_proc_AR_z = diff1r(-1)
     Dz_proc_AR_z = diff1z(+1)

     diffvar => proc_AI_r
     Dr_proc_AI_r = diff1r(-1)
     Dz_proc_AI_r = diff1z(+1)

     diffvar => proc_AI_z
     Dr_proc_AI_z = diff1r(-1)
     Dz_proc_AI_z = diff1z(+1)

!    Derivatives of E.

     diffvar => proc_ER_r
     Dr_proc_ER_r = diff1r(-1)
     Dz_proc_ER_r = diff1z(+1)

     diffvar => proc_ER_z
     Dr_proc_ER_z = diff1r(-1)
     Dz_proc_ER_z = diff1z(+1)

     diffvar => proc_EI_r
     Dr_proc_EI_r = diff1r(-1)
     Dz_proc_EI_r = diff1z(+1)

     diffvar => proc_EI_z
     Dr_proc_EI_z = diff1r(-1)
     Dz_proc_EI_z = diff1z(+1)

!    Derivatives of B.

     diffvar => proc_BR_r
     Dr_proc_BR_r = diff1r(-1)
     Dz_proc_BR_r = diff1z(+1)

     diffvar => proc_BR_z
     Dr_proc_BR_z = diff1r(-1)
     Dz_proc_BR_z = diff1z(+1)

     diffvar => proc_BI_r
     Dr_proc_BI_r = diff1r(-1)
     Dz_proc_BI_r = diff1z(+1)

     diffvar => proc_BI_z
     Dr_proc_BI_z = diff1r(-1)
     Dz_proc_BI_z = diff1z(+1)

!    Advective derivatives.

     if (shift/="none") then

     ! Still missing.

     end if

!    Covariant derivatives (without angular momentum at the moment).

     CovDr_proc_AR_r = Dr_proc_AR_r - chris_rrr*proc_AR_r - chris_zrr*proc_AR_z
     CovDz_proc_AR_r = Dz_proc_AR_r - chris_rzr*proc_AR_r - chris_zzr*proc_AR_z
     CovDr_proc_AR_z = Dr_proc_AR_z - chris_rrz*proc_AR_r - chris_zrz*proc_AR_z
     CovDz_proc_AR_z = Dz_proc_AR_z - chris_rzz*proc_AR_r - chris_zzz*proc_AR_z

     CovDr_proc_ER_r = Dr_proc_ER_r + chris_rrr*proc_ER_r + chris_rrz*proc_ER_z
     CovDr_proc_ER_z = Dr_proc_ER_z + chris_zrr*proc_ER_r + chris_zrz*proc_ER_z
     CovDz_proc_ER_r = Dz_proc_ER_r + chris_rzr*proc_ER_r + chris_rzz*proc_ER_z
     CovDz_proc_ER_z = Dz_proc_ER_z + chris_zzr*proc_ER_r + chris_zzz*proc_ER_z

     CovDr_proc_BR_r = Dr_proc_BR_r + chris_rrr*proc_BR_r + chris_rrz*proc_BR_z
     CovDr_proc_BR_z = Dr_proc_BR_z + chris_zrr*proc_BR_r + chris_zrz*proc_BR_z
     CovDz_proc_BR_r = Dz_proc_BR_r + chris_rzr*proc_BR_r + chris_rzz*proc_BR_z
     CovDz_proc_BR_z = Dz_proc_BR_z + chris_zzr*proc_BR_r + chris_zzz*proc_BR_z

     CovDr_proc_AI_r = Dr_proc_AI_r - chris_rrr*proc_AI_r - chris_zrr*proc_AI_z
     CovDz_proc_AI_r = Dz_proc_AI_r - chris_rzr*proc_AI_r - chris_zzr*proc_AI_z
     CovDr_proc_AI_z = Dr_proc_AI_z - chris_rrz*proc_AI_r - chris_zrz*proc_AI_z
     CovDz_proc_AI_z = Dz_proc_AI_z - chris_rzz*proc_AI_r - chris_zzz*proc_AI_z

     CovDr_proc_EI_r = Dr_proc_EI_r + chris_rrr*proc_EI_r + chris_rrz*proc_EI_z
     CovDr_proc_EI_z = Dr_proc_EI_z + chris_zrr*proc_EI_r + chris_zrz*proc_EI_z
     CovDz_proc_EI_r = Dz_proc_EI_r + chris_rzr*proc_EI_r + chris_rzz*proc_EI_z
     CovDz_proc_EI_z = Dz_proc_EI_z + chris_zzr*proc_EI_r + chris_zzz*proc_EI_z

     CovDr_proc_BI_r = Dr_proc_BI_r + chris_rrr*proc_BI_r + chris_rrz*proc_BI_z
     CovDr_proc_BI_z = Dr_proc_BI_z + chris_zrr*proc_BI_r + chris_zrz*proc_BI_z
     CovDz_proc_BI_r = Dz_proc_BI_r + chris_rzr*proc_BI_r + chris_rzz*proc_BI_z
     CovDz_proc_BI_z = Dz_proc_BI_z + chris_zzr*proc_BI_r + chris_zzz*proc_BI_z

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine auxiliary_matter

