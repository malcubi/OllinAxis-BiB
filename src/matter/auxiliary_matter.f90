!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/matter/auxiliary_matter.f90,v 1.9 2021/08/20 17:39:10 malcubi Exp $

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


! ***************
! ***   END   ***
! ***************

  end subroutine auxiliary_matter

