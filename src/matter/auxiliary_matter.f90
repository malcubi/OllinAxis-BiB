
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

!    We have defined vector quantities with indices as follows:
!       maxw_A: index down
!       maxw_E: index up
!       maxw_B: index up

!    Define the Maxwell vector potential (index up).

     maxw_Au_r = (1.d0/psi4) * (g_A*maxw_A_r + r*g_C*maxw_A_z)
     maxw_Au_z = (1.d0/psi4) * (r*g_C*maxw_A_r + g_B*maxw_A_z)

!    Define the Maxwell electric field (index down).

     maxw_Ed_r = psi4 * (A*maxw_E_r + r*C*maxw_E_z)
     maxw_Ed_z = psi4 * (r*C*maxw_E_r + B*maxw_E_z)

!    Define the Maxwell magnetic field (index down).

     maxw_Bd_r = psi4 * (A*maxw_B_r + r*C*maxw_B_z)
     maxw_Bd_z = psi4 * (r*C*maxw_B_r + B*maxw_B_z)

     if (angmom) then

!       Extra terms to previous ones:

!       Maxwell vector potential (index up).

        maxw_Au_r = maxw_Au_r + (1.d0/psi4) * (r*g_C1*maxw_A_p)
        maxw_Au_z = maxw_Au_z + (1.d0/psi4) * (g_C2*maxw_A_p)

!       Maxwell electric field (index down).

        maxw_Ed_r = maxw_Ed_r + psi4 * (r**3*C1*maxw_E_p)
        maxw_Ed_z = maxw_Ed_z + psi4 * (r**2*C2*maxw_E_p)

!       Maxwell magnetic field (index down).

        maxw_Bd_r = maxw_Bd_r + psi4 * (r**3*C1*maxw_B_p)
        maxw_Bd_z = maxw_Bd_z + psi4 * (r**2*C2*maxw_B_p)

!       Pure angular momentum:

        maxw_Au_p = (1.d0/psi4) * (r*g_C1*maxw_A_r + g_C2*maxw_A_z + (g_H/r**2)*maxw_A_p)
        maxw_Ed_p = psi4 * (r**3*C1*maxw_E_r + r**2*C2*maxw_E_z + r**2*H*maxw_E_p)
        maxw_Bd_p = psi4 * (r**3*C1*maxw_B_r + r**2*C2*maxw_B_z + r**2*H*maxw_B_p)

     end if

!    Covariant derivatives. Although the index p is present in chris_rpp,
!    chris_zpp, chris_prp, chris_pzp, remember that they are nonzero with
!    even without angular momentum.

     CovDr_maxw_A_r = Dr_maxw_A_r - chris_rrr*maxw_A_r - chris_zrr*maxw_A_z
     CovDz_maxw_A_r = Dz_maxw_A_r - chris_rrz*maxw_A_r - chris_zrz*maxw_A_z
     CovDr_maxw_A_z = Dr_maxw_A_z - chris_rrz*maxw_A_r - chris_zrz*maxw_A_z
     CovDz_maxw_A_z = Dz_maxw_A_z - chris_rzz*maxw_A_r - chris_zzz*maxw_A_z
     CovDp_maxw_A_p =             - chris_rpp*maxw_A_r - chris_zpp*maxw_A_z

     CovDr_maxw_E_r = Dr_maxw_E_r + chris_rrr*maxw_E_r + chris_rrz*maxw_E_z
     CovDz_maxw_E_r = Dz_maxw_E_r + chris_rrz*maxw_E_r + chris_rzz*maxw_E_z
     CovDr_maxw_E_z = Dr_maxw_E_z + chris_zrr*maxw_E_r + chris_zrz*maxw_E_z
     CovDz_maxw_E_z = Dz_maxw_E_z + chris_zrz*maxw_E_r + chris_zzz*maxw_E_z
     CovDp_maxw_E_p =             + chris_prp*maxw_E_r + chris_pzp*maxw_E_z   ! WARNING! Here we have chris_prp singular

     CovDr_maxw_B_r = Dr_maxw_B_r + chris_rrr*maxw_B_r + chris_rrz*maxw_B_z
     CovDz_maxw_B_r = Dz_maxw_B_r + chris_rrz*maxw_B_r + chris_rzz*maxw_B_z
     CovDr_maxw_B_z = Dr_maxw_B_z + chris_zrr*maxw_B_r + chris_zrz*maxw_B_z
     CovDz_maxw_B_z = Dz_maxw_B_z + chris_zrz*maxw_B_r + chris_zzz*maxw_B_z
     CovDp_maxw_B_p =             + chris_prp*maxw_B_r + chris_pzp*maxw_B_z   ! WARNING! Here we have chris_prp singular

!    We add terms related to angular momentum.

     if (angmom) then

!       Extra terms to previous ones.

        CovDr_maxw_A_r = CovDr_maxw_A_r - chris_prr*maxw_A_p
        CovDz_maxw_A_r = CovDz_maxw_A_r - chris_prz*maxw_A_p   ! WARNING! Here we have chris_prz singular
        CovDr_maxw_A_z = CovDr_maxw_A_z - chris_prz*maxw_A_p   ! WARNING! Here we have chris_prz singular
        CovDz_maxw_A_z = CovDz_maxw_A_z - chris_pzz*maxw_A_p
        CovDp_maxw_A_p = CovDp_maxw_A_p - chris_ppp*maxw_A_p

        CovDr_maxw_E_r = CovDr_maxw_E_r + chris_rrp*maxw_E_p
        CovDz_maxw_E_r = CovDz_maxw_E_r + chris_rzp*maxw_E_p
        CovDr_maxw_E_z = CovDr_maxw_E_z + chris_zrp*maxw_E_p
        CovDz_maxw_E_z = CovDz_maxw_E_z + chris_zzp*maxw_E_p
        CovDp_maxw_E_p = CovDp_maxw_E_p + chris_ppp*maxw_E_p

        CovDr_maxw_B_r = CovDr_maxw_B_r + chris_rrp*maxw_B_p
        CovDz_maxw_B_r = CovDz_maxw_B_r + chris_rzp*maxw_B_p
        CovDr_maxw_B_z = CovDr_maxw_B_z + chris_zrp*maxw_B_p
        CovDz_maxw_B_z = CovDz_maxw_B_z + chris_zzp*maxw_B_p
        CovDp_maxw_B_p = CovDp_maxw_B_p + chris_ppp*maxw_B_p

!       Pure angular momentum.

        CovDp_maxw_A_r = - chris_rrp*maxw_A_r - chris_zrp*maxw_A_z - chris_prp*maxw_A_p   ! WARNING! Here we have chris_prp singular
        CovDp_maxw_A_z = - chris_rzp*maxw_A_r - chris_zzp*maxw_A_z - chris_pzp*maxw_A_p
        CovDr_maxw_A_p = CovDp_maxw_A_r
        CovDz_maxw_A_p = CovDp_maxw_A_z

        CovDp_maxw_E_r = + chris_rrp*maxw_E_r + chris_rzp*maxw_E_z + chris_rpp*maxw_E_p
        CovDp_maxw_E_z = + chris_zrp*maxw_E_r + chris_zzp*maxw_E_z + chris_zpp*maxw_E_p
        CovDr_maxw_E_p = + chris_prr*maxw_E_r + chris_prz*maxw_E_z + chris_prp*maxw_E_p   ! WARNING! Here we have chris_prp and chris_prz singular
        CovDz_maxw_E_p = + chris_prz*maxw_E_r + chris_pzz*maxw_E_z + chris_pzp*maxw_E_p

        CovDp_maxw_B_r = + chris_rrp*maxw_B_r + chris_rzp*maxw_B_z + chris_rpp*maxw_B_p
        CovDp_maxw_B_z = + chris_zrp*maxw_B_r + chris_zzp*maxw_B_z + chris_zpp*maxw_B_p
        CovDr_maxw_B_p = + chris_prr*maxw_B_r + chris_prz*maxw_B_z + chris_prp*maxw_B_p   ! WARNING! Here we have chris_prp and chris_prz singular
        CovDz_maxw_B_p = + chris_prz*maxw_B_r + chris_pzz*maxw_B_z + chris_pzp*maxw_B_p

     end if

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
     Dr_proc_AR_z = diff1r(+1)
     Dz_proc_AR_z = diff1z(-1)

     diffvar => proc_AI_r
     Dr_proc_AI_r = diff1r(-1)
     Dz_proc_AI_r = diff1z(+1)

     diffvar => proc_AI_z
     Dr_proc_AI_z = diff1r(+1)
     Dz_proc_AI_z = diff1z(-1)

!    Derivatives of E.

     diffvar => proc_ER_r
     Dr_proc_ER_r = diff1r(-1)
     Dz_proc_ER_r = diff1z(+1)

     diffvar => proc_ER_z
     Dr_proc_ER_z = diff1r(+1)
     Dz_proc_ER_z = diff1z(-1)

     diffvar => proc_EI_r
     Dr_proc_EI_r = diff1r(-1)
     Dz_proc_EI_r = diff1z(+1)

     diffvar => proc_EI_z
     Dr_proc_EI_z = diff1r(+1)
     Dz_proc_EI_z = diff1z(-1)

!    Derivatives of B.

     diffvar => proc_BR_r
     Dr_proc_BR_r = diff1r(-1)
     Dz_proc_BR_r = diff1z(+1)

     diffvar => proc_BR_z
     Dr_proc_BR_z = diff1r(+1)
     Dz_proc_BR_z = diff1z(-1)

     diffvar => proc_BI_r
     Dr_proc_BI_r = diff1r(-1)
     Dz_proc_BI_r = diff1z(+1)

     diffvar => proc_BI_z
     Dr_proc_BI_z = diff1r(+1)
     Dz_proc_BI_z = diff1z(-1)

!    Advective derivatives.

     if (shift/="none") then

     ! Still missing.

     end if

!    We have defined vector quantities with index as follow
!       proc_A: index down
!       proc_E: index up
!       proc_B: index up

!    Define the Proca vector potential (index up).

     proc_ARu_r = (1/psi4) * (g_A*proc_AR_r + r*g_C*proc_AR_z)
     proc_AIu_r = (1/psi4) * (g_A*proc_AI_r + r*g_C*proc_AI_z)
     proc_ARu_z = (1/psi4) * (r*g_C*proc_AR_r + g_B*proc_AR_z)
     proc_AIu_z = (1/psi4) * (r*g_C*proc_AI_r + g_B*proc_AI_z)

!    Define the Proca electric field (index down).

     proc_ERd_r = psi4 * (A*proc_ER_r + r*C*proc_ER_z)
     proc_EId_r = psi4 * (A*proc_EI_r + r*C*proc_EI_z)
     proc_ERd_z = psi4 * (r*C*proc_ER_r + B*proc_ER_z)
     proc_EId_z = psi4 * (r*C*proc_EI_r + B*proc_EI_z)

!    Define the Proca magnetic field (index down).

     proc_BRd_r = psi4 * (A*proc_BR_r + r*C*proc_BR_z)
     proc_BId_r = psi4 * (A*proc_BI_r + r*C*proc_BI_z)
     proc_BRd_z = psi4 * (r*C*proc_BR_r + B*proc_BR_z)
     proc_BId_z = psi4 * (r*C*proc_BI_r + B*proc_BI_z)

     if (angmom) then

!       Extra terms to previous one:

!       Proca vector potential (index up).

        proc_ARu_r = proc_ARu_r + (1/psi4) * (r*g_C1*proc_AR_p)
        proc_AIu_r = proc_AIu_r + (1/psi4) * (r*g_C1*proc_AI_p)
        proc_ARu_z = proc_ARu_z + (1/psi4) * (g_C2*proc_AR_p)
        proc_AIu_z = proc_AIu_z + (1/psi4) * (g_C2*proc_AI_p)

!       Proca electric field (index down).

        proc_ERd_r = proc_ERd_r + psi4 * (r**3*C1*proc_ER_p)
        proc_EId_r = proc_EId_r + psi4 * (r**3*C1*proc_EI_p)
        proc_ERd_z = proc_ERd_z + psi4 * (r**2*C2*proc_ER_p)
        proc_EId_z = proc_EId_z + psi4 * (r**2*C2*proc_EI_p)

!       Proca magnetic field (index down).

        proc_BRd_r = proc_BRd_r + psi4 * (r**3*C1*proc_BR_p)
        proc_BId_r = proc_BId_r + psi4 * (r**3*C1*proc_BI_p)
        proc_BRd_z = proc_BRd_z + psi4 * (r**2*C2*proc_BR_p)
        proc_BId_z = proc_BId_z + psi4 * (r**2*C2*proc_BI_p)

!       Pure angular momentum:

        proc_ARu_p = (1/psi4) * (r*g_C1*proc_AR_r + g_C2*proc_AR_z + (g_H/r**2)*proc_AR_p)
        proc_AIu_p = (1/psi4) * (r*g_C1*proc_AI_r + g_C2*proc_AI_z + (g_H/r**2)*proc_AI_p)
        proc_ERd_p = psi4 * (r**3*C1*proc_ER_r + r**2*C2*proc_ER_z + r**2*H*proc_ER_p)
        proc_EId_p = psi4 * (r**3*C1*proc_EI_r + r**2*C2*proc_EI_z + r**2*H*proc_EI_p)
        proc_BRd_p = psi4 * (r**3*C1*proc_BR_r + r**2*C2*proc_BR_z + r**2*H*proc_BR_p)
        proc_BId_p = psi4 * (r**3*C1*proc_BI_r + r**2*C2*proc_BI_z + r**2*H*proc_BI_p)

     end if

!    Covariant derivatives. Although the index p is present in chris_rpp,
!    chris_zpp, chris_prp, chris_pzp, remember that they are nonzero with
!    even without angular momentum.

!    Real part.

     CovDr_proc_AR_r = Dr_proc_AR_r - chris_rrr*proc_AR_r - chris_zrr*proc_AR_z
     CovDz_proc_AR_r = Dz_proc_AR_r - chris_rrz*proc_AR_r - chris_zrz*proc_AR_z
     CovDr_proc_AR_z = Dr_proc_AR_z - chris_rrz*proc_AR_r - chris_zrz*proc_AR_z
     CovDz_proc_AR_z = Dz_proc_AR_z - chris_rzz*proc_AR_r - chris_zzz*proc_AR_z
     CovDp_proc_AR_p =              - chris_rpp*proc_AR_r - chris_zpp*proc_AR_z

     CovDr_proc_ER_r = Dr_proc_ER_r + chris_rrr*proc_ER_r + chris_rrz*proc_ER_z
     CovDz_proc_ER_r = Dz_proc_ER_r + chris_rrz*proc_ER_r + chris_rzz*proc_ER_z
     CovDr_proc_ER_z = Dr_proc_ER_z + chris_zrr*proc_ER_r + chris_zrz*proc_ER_z
     CovDz_proc_ER_z = Dz_proc_ER_z + chris_zrz*proc_ER_r + chris_zzz*proc_ER_z
     CovDp_proc_ER_p =              + chris_prp*proc_ER_r + chris_pzp*proc_ER_z  ! WARNING! Here we have chris_prp singular

     CovDr_proc_BR_r = Dr_proc_BR_r + chris_rrr*proc_BR_r + chris_rrz*proc_BR_z
     CovDz_proc_BR_r = Dz_proc_BR_r + chris_rrz*proc_BR_r + chris_rzz*proc_BR_z
     CovDr_proc_BR_z = Dr_proc_BR_z + chris_zrr*proc_BR_r + chris_zrz*proc_BR_z
     CovDz_proc_BR_z = Dz_proc_BR_z + chris_zrz*proc_BR_r + chris_zzz*proc_BR_z
     CovDp_proc_BR_p =              + chris_prp*proc_BR_r + chris_pzp*proc_BR_z  ! WARNING! Here we have chris_prp singular

!    Imaginary part.

     CovDr_proc_AI_r = Dr_proc_AI_r - chris_rrr*proc_AI_r - chris_zrr*proc_AI_z
     CovDz_proc_AI_r = Dz_proc_AI_r - chris_rrz*proc_AI_r - chris_zrz*proc_AI_z
     CovDr_proc_AI_z = Dr_proc_AI_z - chris_rrz*proc_AI_r - chris_zrz*proc_AI_z
     CovDz_proc_AI_z = Dz_proc_AI_z - chris_rzz*proc_AI_r - chris_zzz*proc_AI_z
     CovDp_proc_AI_p =              - chris_rpp*proc_AI_r - chris_zpp*proc_AI_z

     CovDr_proc_EI_r = Dr_proc_EI_r + chris_rrr*proc_EI_r + chris_rrz*proc_EI_z
     CovDz_proc_EI_r = Dz_proc_EI_r + chris_rrz*proc_EI_r + chris_rzz*proc_EI_z
     CovDr_proc_EI_z = Dr_proc_EI_z + chris_zrr*proc_EI_r + chris_zrz*proc_EI_z
     CovDz_proc_EI_z = Dz_proc_EI_z + chris_zrz*proc_EI_r + chris_zzz*proc_EI_z
     CovDp_proc_EI_p =              + chris_prp*proc_EI_r + chris_pzp*proc_EI_z  ! WARNING! Here we have chris_prp singular

     CovDr_proc_BI_r = Dr_proc_BI_r + chris_rrr*proc_BI_r + chris_rrz*proc_BI_z
     CovDz_proc_BI_r = Dz_proc_BI_r + chris_rrz*proc_BI_r + chris_rzz*proc_BI_z
     CovDr_proc_BI_z = Dr_proc_BI_z + chris_zrr*proc_BI_r + chris_zrz*proc_BI_z
     CovDz_proc_BI_z = Dz_proc_BI_z + chris_zrz*proc_BI_r + chris_zzz*proc_BI_z
     CovDp_proc_BI_p =              + chris_prp*proc_BI_r + chris_pzp*proc_BI_z  ! WARNING! Here we have chris_prp singular

!    We add terms related to angular momentum.

     if (angmom) then

!       Extra terms to previous ones.

!       Real part.

        CovDr_proc_AR_r = CovDr_proc_AR_r - chris_prr*proc_AR_p
        CovDz_proc_AR_r = CovDz_proc_AR_r - chris_prz*proc_AR_p   ! WARNING! Here we have chris_prz singular
        CovDr_proc_AR_z = CovDr_proc_AR_z - chris_prz*proc_AR_p   ! WARNING! Here we have chris_prz singular
        CovDz_proc_AR_z = CovDz_proc_AR_z - chris_pzz*proc_AR_p
        CovDp_proc_AR_p = CovDp_proc_AR_p - chris_ppp*proc_AR_p

        CovDr_proc_ER_r = CovDr_proc_ER_r + chris_rrp*proc_ER_p
        CovDz_proc_ER_r = CovDz_proc_ER_r + chris_rzp*proc_ER_p
        CovDr_proc_ER_z = CovDr_proc_ER_z + chris_zrp*proc_ER_p
        CovDz_proc_ER_z = CovDz_proc_ER_z + chris_zzp*proc_ER_p
        CovDp_proc_ER_p = CovDp_proc_ER_p + chris_ppp*proc_ER_p

        CovDr_proc_BR_r = CovDr_proc_BR_r + chris_rrp*proc_BR_p
        CovDz_proc_BR_r = CovDz_proc_BR_r + chris_rzp*proc_BR_p
        CovDr_proc_BR_z = CovDr_proc_BR_z + chris_zrp*proc_BR_p
        CovDz_proc_BR_z = CovDz_proc_BR_z + chris_zzp*proc_BR_p
        CovDp_proc_BR_p = CovDp_proc_BR_p + chris_ppp*proc_BR_p

!       Imaginary part.

        CovDr_proc_AI_r = CovDr_proc_AI_r - chris_prr*proc_AI_p
        CovDz_proc_AI_r = CovDz_proc_AI_r - chris_prz*proc_AI_p   ! WARNING! Here we have chris_prz singular
        CovDr_proc_AI_z = CovDr_proc_AI_z - chris_prz*proc_AI_p   ! WARNING! Here we have chris_prz singular
        CovDz_proc_AI_z = CovDz_proc_AI_z - chris_pzz*proc_AI_p
        CovDp_proc_AI_p = CovDp_proc_AI_p - chris_ppp*proc_AI_p

        CovDr_proc_EI_r = CovDr_proc_EI_r + chris_rrp*proc_EI_p
        CovDz_proc_EI_r = CovDz_proc_EI_r + chris_rzp*proc_EI_p
        CovDr_proc_EI_z = CovDr_proc_EI_z + chris_zrp*proc_EI_p
        CovDz_proc_EI_z = CovDz_proc_EI_z + chris_zzp*proc_EI_p
        CovDp_proc_EI_p = CovDp_proc_EI_p + chris_ppp*proc_EI_p

        CovDr_proc_BI_r = CovDr_proc_BI_r + chris_rrp*proc_BI_p
        CovDz_proc_BI_r = CovDz_proc_BI_r + chris_rzp*proc_BI_p
        CovDr_proc_BI_z = CovDr_proc_BI_z + chris_zrp*proc_BI_p
        CovDz_proc_BI_z = CovDz_proc_BI_z + chris_zzp*proc_BI_p
        CovDp_proc_BI_p = CovDp_proc_BI_p + chris_ppp*proc_BI_p

!       Pure angular momentum.

!       Real part.

        CovDp_proc_AR_r = - chris_rrp*proc_AR_r - chris_zrp*proc_AR_z - chris_prp*proc_AR_p  ! WARNING! Here we have chris_prp singular
        CovDp_proc_AR_z = - chris_rzp*proc_AR_r - chris_zzp*proc_AR_z - chris_pzp*proc_AR_p
        CovDr_proc_AR_p = CovDp_proc_AR_r
        CovDz_proc_AR_p = CovDp_proc_AR_z

        CovDp_proc_ER_r = + chris_rrp*proc_ER_r + chris_rzp*proc_ER_z + chris_rpp*proc_ER_p
        CovDp_proc_ER_z = + chris_zrp*proc_ER_r + chris_zzp*proc_ER_z + chris_zpp*proc_ER_p
        CovDr_proc_ER_p = + chris_prr*proc_ER_r + chris_prz*proc_ER_z + chris_prp*proc_ER_p  ! WARNING! Here we have chris_prp, chris_prp singular
        CovDz_proc_ER_p = + chris_prz*proc_ER_r + chris_pzz*proc_ER_z + chris_pzp*proc_ER_p

        CovDp_proc_BR_r = + chris_rrp*proc_BR_r + chris_rzp*proc_BR_z + chris_rpp*proc_BR_p
        CovDp_proc_BR_z = + chris_zrp*proc_BR_r + chris_zzp*proc_BR_z + chris_zpp*proc_BR_p
        CovDr_proc_BR_p = + chris_prr*proc_BR_r + chris_prz*proc_BR_z + chris_prp*proc_BR_p  ! WARNING! Here we have chris_prp, chris_prp singular
        CovDz_proc_BR_p = + chris_prz*proc_BR_r + chris_pzz*proc_BR_z + chris_pzp*proc_BR_p

!       Imaginary part.

        CovDp_proc_AI_r = - chris_rrp*proc_AI_r - chris_zrp*proc_AI_z - chris_prp*proc_AI_p  ! WARNING! Here we have chris_prp singular
        CovDp_proc_AI_z = - chris_rzp*proc_AI_r - chris_zzp*proc_AI_z - chris_pzp*proc_AI_p
        CovDr_proc_AI_p = CovDp_proc_AI_r
        CovDz_proc_AI_p = CovDp_proc_AI_z

        CovDp_proc_EI_r = + chris_rrp*proc_EI_r + chris_rzp*proc_EI_z + chris_rpp*proc_EI_p
        CovDp_proc_EI_z = + chris_zrp*proc_EI_r + chris_zzp*proc_EI_z + chris_zpp*proc_EI_p
        CovDr_proc_EI_p = + chris_prr*proc_EI_r + chris_prz*proc_EI_z + chris_prp*proc_EI_p  ! WARNING! Here we have chris_prp, chris_prp singular
        CovDz_proc_EI_p = + chris_prz*proc_EI_r + chris_pzz*proc_EI_z + chris_pzp*proc_EI_p

        CovDp_proc_BI_r = + chris_rrp*proc_BI_r + chris_rzp*proc_BI_z + chris_rpp*proc_BI_p
        CovDp_proc_BI_z = + chris_zrp*proc_BI_r + chris_zzp*proc_BI_z + chris_zpp*proc_BI_p
        CovDr_proc_BI_p = + chris_prr*proc_BI_r + chris_prz*proc_BI_z + chris_prp*proc_BI_p  ! WARNING! Here we have chris_prp, chris_prp singular
        CovDz_proc_BI_p = + chris_prz*proc_BI_r + chris_pzz*proc_BI_z + chris_pzp*proc_BI_p

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine auxiliary_matter

