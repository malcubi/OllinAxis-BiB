!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/sources_geometry.f90,v 1.33 2021/03/23 17:44:43 malcubi Exp $

  subroutine sources_geometry

! *******************************************
! ***   SOURCES FOR EVOLUTION EQUATIONS   ***
! *******************************************

! This routine calculates the sources for the evolution
! equations of the different geometrical variables.
!
! Originally written by Jose Manuel Torres.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  real(8) sigma     ! sigma:  BSSN flavor, 0 for eulerian, 1 for lagrangian.

  real(8) zero,half,third,sixth
  real(8) one,two,three,four,six,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  half  = 0.5d0
  third = 1.d0/3.d0
  sixth = 1.d0/6.d0

  one   = 1.d0
  two   = 2.d0
  three = 3.d0
  four  = 4.d0
  six   = 6.d0

  smallpi = acos(-one)


! ***********************
! ***   BSSN FLAVOR   ***
! ***********************

  !print *, 'Entering sources_geometry'

! Depending on the choice for the evolution equation
! of the conformal volume elements, we have two
! different BSSN schemes, Eulerian and Lagrangian,
! that change the form of the shift terms in the
! sources of the geometric variables.

  if (bssnflavor=='eulerian') then
     sigma = 0.d0
  else if (bssnflavor=='lagrangian') then
     sigma = 1.d0
  end if


! ****************************************
! ***   SOURCES FOR CONFORMAL FACTOR   ***
! ****************************************

! Sources for conformal factor.  Notice that if
! chimethod is true we evolve chi instead of phi.
!
! Notice also that sphi must always be defined, even
! if we are using chimethod=true.  This is because
! it might later be used to construct other sources.

  sphi = - sixth*alpha*trK

  if (formulation=="z4c") then
     sphi = sphi - two*sixth*alpha*z4theta
  end if

  if (shift/="none") then
     sphi = sphi + beta_r*DAr_phi + beta_z*DAz_phi &
          + sixth*sigma*DIV_beta
  end if

  if (chimethod) then

     schi = dble(chipower)*sixth*alpha*chi*trK

     if (formulation=="z4c") then
        schi = schi + two*dble(chipower)*sixth*alpha*chi*z4theta
     end if

     if (shift/="none") then
        schi = schi + beta_r*DAr_chi + beta_z*DAz_chi &
             - sigma*dble(chipower)*sixth*chi*DIV_beta
     end if

  end if


! ****************************************
! ***   SOURCES FOR CONFORMAL METRIC   ***
! ****************************************

  sA = - two*alpha*KTA
  sB = - two*alpha*KTB
  sH = - two*alpha*KTH
  sC = - two*alpha*KTC

! Shift terms.

  if (shift/="none") then

     sA = sA + beta_r*DAr_A + beta_z*DAz_A &
        + two*(A*Dr_beta_r + r*C*Dr_beta_z) - two*third*sigma*A*DIV_beta

     sB = sB + beta_r*DAr_B + beta_z*DAz_B &
        + two*(B*Dz_beta_z + r*C*Dz_beta_r) - two*third*sigma*B*DIV_beta

     sH = sH + beta_r*DAr_H + beta_z*DAz_H &
        + two*H*beta_r/r - two*third*sigma*H*DIV_beta

     sC = sC + beta_r*DAr_C + beta_z*DAz_C &
        + C*beta_r/r + C*(Dr_beta_r + Dz_beta_z) + (B*Dr_beta_z + A*Dz_beta_r)/r &
        - two*third*sigma*C*DIV_beta

     if (angmom) then
        sA = sA + two*r**3*C1*Dr_beta_p
        sB = sB + two*r**2*C2*Dz_beta_p
        sC = sC + r*(C2*Dr_beta_p +r*C1*Dz_beta_p)
     end if

  end if

! Sources for rotational metric functions.

  if (angmom) then

     sC1 = - two*alpha*KTC1
     sC2 = - two*alpha*KTC2

!    Shift terms.

     if (shift/="none") then
        sC1 = sC1 + beta_r*DAr_C1 + beta_z*DAz_C1 &
            + three*C1*beta_r/r + C1*Dr_beta_r + (C2*Dr_beta_z + H*Dr_beta_p)/r &
            - two*third*sigma*C1*DIV_beta
        sC2 = sC2 + beta_r*DAr_C2 + beta_z*DAz_C2 &
            + two*C2*beta_r/r + r*C1*Dz_beta_r + C2*Dz_beta_z + H*Dz_beta_p &
            - two*third*sigma*C2*DIV_beta
     end if

  end if


! ******************************
! ***   SOURCES FOR LAMBDA   ***
! ******************************

! Sources for lambda = (A - H)/r**2

  if (.not.nolambda) then

     slambda = - two*alpha*Alambda

!    Shift terms.

     if (shift/="none") then

        slambda = slambda + beta_r*DAr_lambda + beta_z*DAz_lambda &
                + two*(lambda*beta_r + C*Dr_beta_z)/r &
                + two*(ft1*(lambda*Dr_beta_r + H*DD_beta_rr/r) &
                + (one-ft1)*(lambda*beta_r + A*DD_beta_rr)/r ) &
                - two*third*sigma*lambda*DIV_beta

        if (angmom) then
           slambda = slambda + two*r*C1*Dr_beta_p 
        end if

     end if

!    Dissipation.

     !if (geodiss/=0.d0) then
     !   evolvevar => lambda
     !   sourcevar => slambda
     !   call dissipation(+1,+1,geodiss)
     !end if

  end if


! ****************************************************
! ***   SOURCES FOR TRACE OF EXTRINSIC CURVATURE   ***
! ****************************************************

  if (slicing=="maximal") then

!    For maximal slicing the source for trK is set to zero.

     strK = zero

  else

!    For any other slicing we calculate the full source term:
!    which is given by:
!
!                 2                 ij
!    strK  =  -  D alpha  +  alpha K  K   +  4 pi alpha ( rho + S )
!                                      ij
!
!           2
!    where D  is the physical Laplacian.

!    Terms coming from Laplacian of lapse.

     strK = - Lapla_alpha

!    Quadratic terms.

     strK = strK + alpha*K2

!    Damping term for Z4c.

     if (formulation=="z4c") then
        strK = strK + alpha*kappa1*(one - kappa2)*z4theta
     end if

!    Shift terms.

     if (shift/="none") then
        strK = strK + beta_r*DAr_trK + beta_z*DAz_trK
     end if

!    Matter terms.

     if (mattertype/="vacuum") then
        strK = strK + four*smallpi*alpha*(rho + trS)
     end if

!    Dissipation.

     if (geodiss/=0.d0) then
        evolvevar => trK
        sourcevar => strK
        call dissipation(+1,+1,geodiss)
     end if

  end if


! *****************************************************
! ***   SOURCES FOR TRACELESS EXTRINSIC CURVATURE   ***
! *****************************************************

! Quadratic terms.

  sKTA = alpha*((trK + two*z4theta)*KTA - two*KT2_A)
  sKTB = alpha*((trK + two*z4theta)*KTB - two*KT2_B)
  sKTH = alpha*((trK + two*z4theta)*KTH - two*KT2_H)
  sKTC = alpha*((trK + two*z4theta)*KTC - two*KT2_C)

! Ricci tensor and derivatives of lapse.

  sKTA = sKTA + (-D2cov_alpha_A + alpha*RIC_A)/psi4 &
       - third*A*(-Lapla_alpha + alpha*RSCAL)
  sKTB = sKTB + (-D2cov_alpha_B + alpha*RIC_B)/psi4 &
       - third*B*(-Lapla_alpha + alpha*RSCAL)
  sKTH = sKTH + (-D2cov_alpha_H + alpha*RIC_H)/psi4 &
       - third*H*(-Lapla_alpha + alpha*RSCAL)
  sKTC = sKTC + (-D2cov_alpha_C + alpha*RIC_C)/psi4 &
       - third*C*(-Lapla_alpha + alpha*RSCAL)

! Shift terms identical to metric terms taking gamma->KT.

  if (shift/="none") then

     sKTA = sKTA + beta_r*DAr_KTA + beta_z*DAz_KTA + two*(KTA*Dr_beta_r + r*KTC*Dr_beta_z) &
          - two*third*sigma*KTA*DIV_beta
     sKTB = sKTB + beta_r*DAr_KTB + beta_z*DAz_KTB + two*(KTB*Dz_beta_z + r*KTC*Dz_beta_r) &
          - two*third*sigma*KTB*DIV_beta
     sKTH = sKTH + beta_r*DAr_KTH + beta_z*DAz_KTH + two*KTH*beta_r/r &
          - two*third*sigma*KTH*DIV_beta
     sKTC = sKTC + beta_r*DAr_KTC + beta_z*DAz_KTC + KTC*beta_r/r &
          + KTC*(Dr_beta_r + Dz_beta_z) + (KTB*Dr_beta_z + KTA*Dz_beta_r)/r &
          - two*third*sigma*KTC*DIV_beta

     if (angmom) then
        sKTA = sKTA + two*r**3*KTC1*Dr_beta_p
        sKTB = sKTB + two*r**2*KTC2*Dz_beta_p
        sKTC = sKTC + r*(KTC2*Dr_beta_p + r*KTC1*Dz_beta_p)
     end if

  end if

! Sources for angular momentum coefficients.

  if (angmom) then

!    Quadratic terms.

     sKTC1 = alpha*((trK + two*z4theta)*KTC1 - two*KT2_C1)
     sKTC2 = alpha*((trK + two*z4theta)*KTC2 - two*KT2_C2)

!    Derivatives.

     sKTC1 = sKTC1 + (-D2cov_alpha_C1 + alpha*RIC_C1)/psi4 &
           - third*C1*(-Lapla_alpha + alpha*RSCAL)
     sKTC2 = sKTC2 + (-D2cov_alpha_C2 + alpha*RIC_C2)/psi4 &
           - third*C2*(-Lapla_alpha + alpha*RSCAL)

     if (shift/="none") then
        sKTC1 = sKTC1 + beta_r*DAr_KTC1 + beta_z*DAz_KTC1 + three*KTC1*beta_r/r &
              + KTC1*Dr_beta_r + (KTC2*Dr_beta_z + KTH*Dr_beta_p)/r &
              - two*third*sigma*KTC1*DIV_beta
        sKTC2 = sKTC2 + beta_r*DAr_KTC2 + beta_z*DAz_KTC2 + two*KTC2*beta_r/r &
              + r*KTC1*Dz_beta_r + KTC2*Dz_beta_z + KTH*Dz_beta_p &
              - two*third*sigma*KTC2*DIV_beta
     end if

  end if

! Matter terms.

  if (mattertype/="vacuum") then

     sKTA = sKTA - 8.d0*smallpi*alpha*(S_A/psi4 - third*A*trS)
     sKTB = sKTB - 8.d0*smallpi*alpha*(S_B/psi4 - third*B*trS)
     sKTH = sKTH - 8.d0*smallpi*alpha*(S_H/psi4 - third*H*trS)
     sKTC = sKTC - 8.d0*smallpi*alpha*(S_C/psi4 - third*C*trS)

     if (angmom) then
        sKTC1= sKTC1 - 8.d0*smallpi*alpha*(S_C1/psi4 - third*C1*trS)
        sKTC2= sKTC2 - 8.d0*smallpi*alpha*(S_C2/psi4 - third*C2*trS)
     end if

  end if

! Dissipation.

  if (geodiss/=0.d0) then

     evolvevar => KTA
     sourcevar => sKTA
     call dissipation(+1,+1,geodiss)

     evolvevar => KTB
     sourcevar => sKTB
     call dissipation(+1,+1,geodiss)

     evolvevar => KTH
     sourcevar => sKTH
     call dissipation(+1,+1,geodiss)

     evolvevar => KTC
     sourcevar => sKTC
     call dissipation(+1,-1,geodiss)

     if (angmom) then

        evolvevar => KTC1
        sourcevar => sKTC1
        call dissipation(+1,+1,geodiss)

        evolvevar => KTC2
        sourcevar => sKTC2
        call dissipation(+1,-1,geodiss)

     end if

  end if


! *******************************
! ***   SOURCES FOR ALAMBDA   ***
! *******************************

! Sources for Alambda = (KTA - KTH)/r**2

  if (.not.nolambda) then

!    Quadratic terms.

     sAlambda = alpha*(trK*Alambda - two*KT2_lambda)

!    Derivatives.

     sAlambda = sAlambda + (-D2cov_alpha_lambda + alpha*RIC_lambda)/psi4 &
              - third*lambda*(-Lapla_alpha + alpha*RSCAL)

!    Shift terms.

     if (shift/="none") then

        sAlambda = sAlambda + beta_r*DAr_Alambda + beta_z*DAz_Alambda &
                 + two*(Alambda*beta_r + KTC*Dr_beta_z)/r &
                 + two*(ft2*(Alambda*Dr_beta_r + KTH*DD_beta_rr/r) &
                 + (one-ft2)*(Alambda*beta_r + KTA*DD_beta_rr)/r) &
                 - two*third*sigma*Alambda*DIV_beta

        if (angmom) then
           sAlambda = sAlambda + two*r*KTC1*Dr_beta_p 
        end if

     end if

!    Matter terms.

     if (mattertype/="vacuum") then
        sAlambda = sAlambda - 8.d0*smallpi*alpha*(S_lambda/psi4 - third*lambda*trS)
     end if

!    Dissipation.

     if (geodiss/=0.d0) then
        evolvevar => Alambda
        sourcevar => sAlambda
        call dissipation(+1,+1,geodiss)
     end if

  end if


! **********************************
! ***   SOURCES FOR BSSN DELTA   ***
! **********************************

! Sources for Delta_r.

  if (.not.noDelta_r) then

!    Standard terms.

     sDelta_r = - two*(KT_up_rr*Dr_alpha + KT_up_rz*Dz_alpha) &
              - two*alpha*DIVA_KT_r

!    Multiple of momentum constraint.

     sDelta_r = sDelta_r + eta*alpha*(DIV_KT_r &
              + six*(KT_up_rr*Dr_phi + KT_up_rz*Dz_phi) &
              - two*third*(g_A*Dr_trK + r*g_C*Dz_trK))

!    Terms for Z4c formulation.

     if (formulation=="z4c") then
        sDelta_r = sDelta_r - third*eta*alpha*(g_A*Dr_z4theta + r*g_C*Dz_z4theta) &
                 - two*alpha*kappa1*(Delta_r + (Dr_g_A + r*Dz_g_C + r*g_lambda &
                 + half*(g_A*Dr_hdet + r*g_C*Dz_hdet)*ihdet))
     end if

!    Shift terms.

     if (shift/="none") then

!       Lie derivative along shift.

        sDelta_r = sDelta_r &
                 + beta_r*DAr_Delta_r + beta_z*DAz_Delta_r &
                 - Delta_r*Dr_beta_r - Delta_z*Dz_beta_r

!       Terms from the flat laplacian.

        sDelta_r = sDelta_r + two*r*g_C*Drz_beta_r &
                 + g_A*Drr_beta_r + g_B*Dzz_beta_r + g_H*DD_beta_rr

!       Angular momentum terms.

        if (angmom) then
           sDelta_r = sDelta_r - two*r*(r*g_C1*Dr_beta_p + g_C2*Dz_beta_p)
        end if

!       Sigma terms.

        sDelta_r = sDelta_r + third*sigma*(g_A*Dr_DIV_beta &
                 + r*g_C*Dz_DIV_beta + two*Delta_r*DIV_beta)

     end if

!    Matter terms.

     if (mattertype/="vacuum") then
        sDelta_r = sDelta_r - 8.d0*smallpi*eta*alpha*J_r*psi4
     end if

!    Dissipation.

     if (geodiss/=0.d0) then
        evolvevar => Delta_r
        sourcevar => sDelta_r
        call dissipation(-1,+1,geodiss)
     end if

  end if

! Sources for Delta_z.

  if (.not.noDelta_z) then

!    Standard terms.

     sDelta_z = - two*(KT_up_rz*Dr_alpha + KT_up_zz*Dz_alpha) &
              - two*alpha*DIVA_KT_z

!    Multiple of momentum constraint.

     sDelta_z = sDelta_z + eta*alpha*(DIV_KT_z &
              + six*(KT_up_rz*Dr_phi + KT_up_zz*Dz_phi) &
              - two*third*(r*g_C*Dr_trK + g_B*Dz_trK))

!    Terms for Z4c formulation.

     if (formulation=="z4c") then
        sDelta_z = sDelta_z - third*eta*alpha*(r*g_C*Dr_z4theta + g_B*Dz_z4theta) &
                 - two*alpha*kappa1*(DelDef_z + (Dr_g_C*r + Dz_g_B + two*g_C &
                 + half*(r*g_C*Dr_hdet + g_B*Dz_hdet)*ihdet))
     end if

!    Shift terms.

     if (shift/="none") then

!       Lie derivative along shift.

        sDelta_z = sDelta_z &
                 + beta_r*DAr_Delta_z + beta_z*DAz_Delta_z &
                 - Delta_r*Dr_beta_z - Delta_z*Dz_beta_z

!       Terms from the flat laplacian.

        sDelta_z = sDelta_z + two*r*g_C*Drz_beta_z &
                 + g_A*Drr_beta_z + g_B*Dzz_beta_z + g_H*Dr_beta_z/r

!       Sigma terms.

        sDelta_z = sDelta_z + third*sigma*(r*g_C*Dr_DIV_beta &
                 + g_B*Dz_DIV_beta + two*Delta_z*DIV_beta)

     end if

!    Matter terms.

     if (mattertype/="vacuum") then
        sDelta_z = sDelta_z - 8.d0*smallpi*eta*alpha*J_z*psi4
     end if

!    Dissipation.

     if (geodiss/=0.d0) then
        evolvevar => Delta_z
        sourcevar => sDelta_z
        call dissipation(+1,-1,geodiss)
     end if

  end if

! Sources for Delta_p.

  if ((.not.noDelta_p).and.angmom) then

!    Standard terms.

     sDelta_p = - two*(KT_up_rp*Dr_alpha + KT_up_zp*Dz_alpha) &
              - two*alpha*DIVA_KT_p

!    Multiple of momentum constraint.

     sDelta_p = sDelta_p + eta*alpha*(DIV_KT_p &
              + six*(KT_up_rp*Dr_phi + KT_up_zp*Dz_phi) &
              - two*third*(r*g_C1*Dr_trK + g_C2*Dz_trK))

!    Terms for Z4c formulation.

     if (formulation=="z4c") then
        sDelta_p = sDelta_p - third*eta*alpha*(r*g_C1*Dr_z4theta + g_C2*Dz_z4theta) &
                 - two*alpha*kappa1*(DelDef_p + (Dr_g_C1*r + Dz_g_C2 + four*g_C1 &
                 + half*(r*g_C1*Dr_hdet + g_C2*Dz_hdet)*ihdet))
     end if

!    Shift terms.

     if (shift/="none") then

!       Lie derivative along shift.

        sDelta_p = sDelta_p &
                 + beta_r*DAr_Delta_p + beta_z*DAz_Delta_p &
                 - Delta_r*Dr_beta_p - Delta_z*Dz_beta_p

!       Terms from the flat laplacian.

        sDelta_p = sDelta_p + two*r*g_C*Drz_beta_p + g_A*Drr_beta_p+g_B*Dzz_beta_p &
                 + two*(g_C1*Dr_beta_r + g_C*Dz_beta_p) &
                 + (Dr_beta_p*(two*g_A + g_H) + two*(g_C2*Dz_beta_r - g_C1*beta_r))/r

!       Sigma terms.

        sDelta_p = sDelta_p + third*sigma*(r*g_C1*Dr_DIV_beta &
                 + g_C2*Dz_DIV_beta + two*Delta_p*DIV_beta)

     end if

!    Matter terms.

     if (mattertype/="vacuum") then
        sDelta_p = sDelta_p - 8.d0*smallpi*eta*alpha*J_p*psi4
     end if

!    Dissipation.

     if (geodiss/=0.d0) then
        evolvevar => Delta_p
        sourcevar => sDelta_p
        call dissipation(+1,+1,geodiss)
     end if

  end if


! ********************************
! ***   SOURCES FOR Z4 THETA   ***
! ********************************

! Sources for z4theta. The source is essentially the
! the Hamiltonian constraint, plus a shift term and
! a damping term.  The only special thing is that
! trK is replaced with: trK + 2*z4theta.

  if (formulation=="z4c") then

!    Hamiltonian constraint term.

     sz4theta = half*alpha*(RSCAL - KT2 + two*third*(trK + two*z4theta)**2)

     if (mattertype/="vacuum") then
        sz4theta = sz4theta - 8.d0*smallpi*alpha*rho
     end if

!    Shift terms.

     if (shift/="none") then
        sz4theta = sz4theta + beta_r*DAr_z4theta + beta_z*DAz_z4theta
     end if

!    Damping terms.

     sz4theta = sz4theta - alpha*kappa1*(two + kappa2)*z4theta

!    Dissipation.

     if (geodiss/=0.d0) then
        evolvevar => z4theta
        sourcevar => sz4theta
        call dissipation(+1,+1,geodiss)
     end if

! For BSSN the source is set to 0.

  else

     sz4theta = zero

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_geometry

