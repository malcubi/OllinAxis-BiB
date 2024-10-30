!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/auxiliary_geometry.f90,v 1.57 2021/03/10 19:18:40 malcubi Exp $

  subroutine auxiliary_geometry

! *******************************************
! ***   AUXLIARY VARIABLES FOR GEOMETRY   ***
! *******************************************

! Auxiliary quantities for geometry: derivatives, inverse metric,
! Ricci scalar, Laplacian of the lapse, etc.
!
! This routine is essentially copied (with some adjustments) from the one
! in the old version of the code OllinAxis.  It was written originally
! by Jose Manuel Torres.
!
! The changes have to do mainly with the way derivative routines are
! calculated (we now use pointers), and comments here and there.

! Include modules.

  use param
  use arrays
  use derivatives
  use derivadvect

! Extra variables.

  implicit none

  integer i,j

  real(8) zero,half,third
  real(8) one,two,three,four


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  one   = 1.d0
  two   = 2.d0
  three = 3.d0
  four  = 4.d0

  half  = 0.5d0
  third = 1.d0/3.d0


! *****************
! ***   LAPSE   ***
! *****************

! Derivatives of lapse.

  diffvar => alpha

  Dr_alpha = diff1r(+1)
  Dz_alpha = diff1z(+1)

  Drr_alpha = diff2r(+1)
  Dzz_alpha = diff2z(+1)
  Drz_alpha = diff2rz(+1,+1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_alpha = diffadvr(+1)
        DAz_alpha = diffadvz(+1)
     else
        DAr_alpha = Dr_alpha
        DAz_alpha = Dz_alpha
     end if
  end if

! Derivatives of dtalpha.

  diffvar => dtalpha

  Dr_dtalpha = diff1r(+1)
  Dz_dtalpha = diff1z(+1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_dtalpha = diffadvr(+1)
        DAz_dtalpha = diffadvz(+1)
     else
        DAr_dtalpha = Dr_dtalpha
        DAz_dtalpha = Dz_dtalpha
     end if
  end if

! Regularization.
  
  auxarray = Dr_alpha/r

  diffvar => auxarray
  DD_alphar = diff1r(+1)


! **************************************
! ***   BONA-MASSO GAUGE FUNCTIONS   ***
! **************************************

! Here we calculate the Bona-Masso gauge function
! which is used for the slicing condition, and
! also to calculate gauge speeds. Notice that I
! actually define falpha as:
!
! falpha = alpha**2 f(alpha)
!
! This is in order to avoid divisions with
! alpha for small alpha.


! GEODESIC OR MAXIMAL.

  if ((slicing=="static").or.(slicing=="maximal")) then

!    In this case we just set falpha it to zero.

     falpha = 0.d0


! HARMONIC FAMILY.

  else if (index(slicing,"harmonic")/=0) then

!    In this case the gauge function is:
!
!    f  =  gauge_f   =>  falpha = alpha**2*gauge_f
!
!    True harmonic slicing corresponds to the case when
!    the constant gauge_f is equal to 1.

     falpha = alpha**2*gauge_f


! 1+LOG FAMILY.

  else if (index(slicing,"1+log")/=0) then

!    In this case the gauge function is:
!
!    f  =  gauge_f/alpha  =>  falpha = alpha*gauge_f
!
!    Standard 1+log slicing corresponds to the case when
!    the constant gauge_f is equal to 2.

     falpha = alpha*gauge_f


! SHOCK AVOIDING FAMILY.

  else if (index(slicing,"shockavoid")/=0) then

!    In this case the gauge function is:
!
!    f = 1 + k/alpha^2  =>  falpha = alpha**2 + k
!
!    with k=constant. The whole family avoids shocks.
!    (See: Class.Quant.Grav. 20 (2003) 607-624; gr-qc/0210050)
!
!    Below we take k=gauge_f-1.  This is just to guarantee
!    that for alpha=1 we always have falpha=gauge_f (as in
!    the harmonic and 1+log cases above).

     falpha = alpha**2 + (gauge_f-1.d0)

  end if


! ****************************
! ***   CONFORMAL FACTOR   ***
! ****************************

! Calculate phi, chi, their derivatives,
! and psi = exp(phi).

  if (.not.chimethod.or.(time==zero)) then

     chi = exp(-dble(chipower)*phi)

     diffvar => phi

     Dr_phi = diff1r(+1)
     Dz_phi = diff1z(+1)
     
     Drr_phi = diff2r(+1)
     Dzz_phi = diff2z(+1)
     Drz_phi = diff2rz(+1,+1)

     if (shift/="none") then
        if (shiftadvect) then
           DAr_phi = diffadvr(+1)
           DAz_phi = diffadvz(+1)
        else
           DAr_phi = Dr_phi
           DAz_phi = Dz_phi
        end if
     end if

  else

     phi = - log(chi)/dble(chipower)

     diffvar => chi

     Dr_chi = diff1r(+1)
     Dz_chi = diff1z(+1)

     Drr_chi = diff2r(+1)
     Dzz_chi = diff2z(+1)
     Drz_chi = diff2rz(+1,+1)

     Dr_phi = - Dr_chi/chi/dble(chipower)
     Dz_phi = - Dz_chi/chi/dble(chipower)

     Drr_phi = - (Drr_chi/chi - (Dr_chi/chi)**2)/dble(chipower)
     Dzz_phi = - (Dzz_chi/chi - (Dz_chi/chi)**2)/dble(chipower)
     Drz_phi = - (Drz_chi/chi - Dr_chi*Dz_chi/chi**2)/dble(chipower)

     if (shift/="none") then

        if (shiftadvect) then
           DAr_chi = diffadvr(+1)
           DAz_chi = diffadvz(+1)
        else
           DAr_chi = Dr_chi
           DAz_chi = Dz_chi
        end if

        DAr_phi = - DAr_chi/chi/dble(chipower)
        DAz_phi = - DAz_chi/chi/dble(chipower)

     end if

  end if

! psi, psi2, psi4.

  psi  = exp(phi)
  psi2 = psi**2
  psi4 = psi**4

! DD_phir = d ( d phi / r )
!            r   r

  auxarray = Dr_phi/r

  diffvar => auxarray
  DD_phir = diff1r(+1)


! *******************************************
! ***   DETERMINANT OF CONFORMAL METRIC   ***
! *******************************************

! Notice that hdet is actually the ratio between the
! conformal determinant and the flat one, so we are
! dividing a factor of r**2. Also, for Lagrangian
! evolutions (the default), hdet must remain equal
! to its initial value.
!
! At t=0 the routine initial.f90 already made sure
! that the determinant of the conformal metric is
! equal to the flat one, so we must set hdet=1.

  if (time==zero) then
     hdet  = one
     ihdet = one
     Dr_hdet = zero
     Dz_hdet = zero
  end if

! Eulerian evolution.

  if (bssnflavor=="eulerian") then

!    Determinant.

     if (angmom) then
        hdet = A*B*H - r**2*(H*C**2 + A*C2**2 + r**2*C1*(B*C1 - two*C*C2))
     else
        hdet = (A*B - (r*C)**2)*H
     end if

!    Inverse.

     ihdet = one/hdet

!    Derivatives.

     diffvar => hdet
     Dr_hdet = diff1r(+1)
     Dz_hdet = diff1z(+1)

  end if

! For Lagrangian evolutions make sure that the determinant
! of the conformal metric remains time indepedent.

  if (bssnflavor=="lagrangian") then
     if (.not.evolveH) then
        if (.not.angmom) then
           H = hdet/(A*B - (r*C)**2)
        else
           H = (hdet + r**2*(A*C2**2 + r**2*C1*(B*C1 - two*C*C2)))/(A*B - (r*C)**2)
        end if
     else
        call factordet
     end if
  end if


! **************************
! ***   SPATIAL METRIC   ***
! **************************

! Derivatives of conformal metric component A.

  diffvar => A
  Dr_A  = diff1r(+1)
  Dz_A  = diff1z(+1)
  Drr_A = diff2r(+1)
  Dzz_A = diff2z(+1)
  Drz_A = diff2rz(+1,+1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_A = diffadvr(+1)
        DAz_A = diffadvz(+1)
     else
        DAr_A = Dr_A
        DAz_A = Dz_A
     end if
  end if

! Derivatives of conformal metric component B.

  diffvar => B
  Dr_B  = diff1r(+1)
  Dz_B  = diff1z(+1)
  Drr_B = diff2r(+1)
  Dzz_B = diff2z(+1)
  Drz_B = diff2rz(+1,+1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_B = diffadvr(+1)
        DAz_B = diffadvz(+1)
     else
        DAr_B = Dr_B
        DAz_B = Dz_B
     end if
  end if

! Derivatives of conformal metric component H.

  diffvar => H
  Dr_H  = diff1r(+1)
  Dz_H  = diff1z(+1)
  Drr_H = diff2r(+1)
  Dzz_H = diff2z(+1)
  Drz_H = diff2rz(+1,+1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_H = diffadvr(+1)
        DAz_H = diffadvz(+1)
     else
        DAr_H = Dr_H
        DAz_H = Dz_H
     end if
  end if

! Derivatives of conformal metric component C.

  diffvar => C
  Dr_C  = diff1r(+1)
  Dz_C  = diff1z(-1)
  Drr_C = diff2r(+1)
  Dzz_C = diff2z(-1)
  Drz_C = diff2rz(+1,-1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_C = diffadvr(+1)
        DAz_C = diffadvz(-1)
     else
        DAr_C = Dr_C
        DAz_C = Dz_C
     end if
  end if

! Derivatives of conformal metric components C1 and C2.

  if (angmom) then

     diffvar => C1
     Dr_C1  = diff1r(+1)
     Dz_C1  = diff1z(+1)
     Drr_C1 = diff2r(+1)
     Dzz_C1 = diff2z(+1)
     Drz_C1 = diff2rz(+1,+1)

     if (shift/="none") then
        if (shiftadvect) then
           DAr_C1 = diffadvr(+1)
           DAz_C1 = diffadvz(+1)
        else
           DAr_C1 = Dr_C1
           DAz_C1 = Dz_C1
        end if
     end if

     diffvar => C2
     Dr_C2  = diff1r(+1)
     Dz_C2  = diff1z(-1)
     Drr_C2 = diff2r(+1)
     Dzz_C2 = diff2z(-1)
     Drz_C2 = diff2rz(+1,-1)

     if (shift/="none") then
        if (shiftadvect) then
           DAr_C2 = diffadvr(+1)
           DAz_C2 = diffadvz(-1)
        else
           DAr_C2 = Dr_C2
           DAz_C2 = Dz_C2
        end if
     end if

  end if

! Physical spatial metric.

  A_PHYS = psi4*A
  B_PHYS = psi4*B
  H_PHYS = psi4*H
  C_PHYS = psi4*C

  if (angmom) then
     C1_PHYS = psi4*C1
     C2_PHYS = psi4*C2
  end if


! ************************************
! ***   INVERSE CONFORMAL METRIC   ***
! ************************************

! The components of the conformal inverse metric
! are given in terms of the functions g_f as:
!
! gamma^rr  =  g_A
!
! gamma^zz  =  g_B
!
! gamma^pp  =  g_H/r**2
!
! gamma^rz  =  r g_C
!
! gamma^rp  =  r g_C1
!
! gamma^zp  =  g_C2

  if (angmom) then

     g_A  = (B*H - (r*C2)**2)*ihdet
     g_B  = (A*H - (r**2*C1)**2)*ihdet
     g_H  = (A*B - (r*C)**2)*ihdet
     g_C  = (r**2*C1*C2 - C*H)*ihdet
     g_C1 = (C*C2 - B*C1)*ihdet
     g_C2 = (r**2*C*C1 - A*C2)*ihdet

  else

     g_A = B*H*ihdet
     g_B = A*H*ihdet
     g_H = one/H
     g_C = - C*H*ihdet

  end if

! Derivatives.

  diffvar => g_A
  Dr_g_A = diff1r(+1)
  Dz_g_A = diff1z(+1)

  diffvar => g_B
  Dr_g_B = diff1r(+1)
  Dz_g_B = diff1z(+1)

  diffvar => g_H
  Dr_g_H = diff1r(+1)
  Dz_g_H = diff1z(+1)

  diffvar => g_C
  Dr_g_C = diff1r(+1)
  Dz_g_C = diff1z(-1)

  if (angmom) then

     diffvar => g_C1
     Dr_g_C1 = diff1r(+1)
     Dz_g_C1 = diff1z(+1)

     diffvar => g_C2
     Dr_g_C2 = diff1r(+1)
     Dz_g_C2 = diff1z(-1)

  end if


! *****************
! ***   SHIFT   ***
! *****************

  if (shift/="none")then

!    Derivatives of beta_r.

     diffvar => beta_r
     Dr_beta_r  = diff1r(-1)
     Dz_beta_r  = diff1z(+1)
     Drr_beta_r = diff2r(-1)
     Dzz_beta_r = diff2z(+1)
     Drz_beta_r = diff2rz(-1,+1)

     if (shiftadvect) then
        DAr_beta_r = diffadvr(-1)
        DAz_beta_r = diffadvz(+1)
     else
        DAr_beta_r = Dr_beta_r
        DAz_beta_r = Dz_beta_r
     end if

!    Derivatives of beta_z.

     diffvar => beta_z
     Dr_beta_z  = diff1r(+1)
     Dz_beta_z  = diff1z(-1)
     Drr_beta_z = diff2r(+1)
     Dzz_beta_z = diff2z(-1)
     Drz_beta_z = diff2rz(+1,-1)

     if (shiftadvect) then
        DAr_beta_z = diffadvr(+1)
        DAz_beta_z = diffadvz(-1)
     else
        DAr_beta_z = Dr_beta_z
        DAz_beta_z = Dz_beta_z
     end if

!    Derivatives of beta_p.

     if (angmom) then

        diffvar => beta_p
        Dr_beta_p  = diff1r(+1)
        Dz_beta_p  = diff1z(+1)
        Drr_beta_p = diff2r(+1)
        Dzz_beta_p = diff2z(+1)
        Drz_beta_p = diff2rz(+1,+1)

        if (shiftadvect) then
           DAr_beta_p = diffadvr(+1)
           DAz_beta_p = diffadvz(+1)
        else
           DAr_beta_p = Dr_beta_p
           DAz_beta_p = Dz_beta_p
        end if

     end if

!    Derivatives of fdriver.

     diffvar => fdriver

     Dr_fdriver = diff1r(+1)
     Dz_fdriver = diff1z(+1)

     if (shiftadvect) then
        DAr_fdriver = diffadvr(+1)
        DAz_fdriver = diffadvz(+1)
     else
        DAr_fdriver = Dr_fdriver
        DAz_fdriver = Dz_fdriver
     end if

!                         r
!    DD_beta_rr = d ( beta / r )
!                  r

     auxarray = beta_r/r

     diffvar => auxarray
     DD_beta_rr = diff1r(+1)

!    Derivatives of dtbeta_r.

     diffvar => dtbeta_r

     Dr_dtbeta_r  = diff1r(-1)
     Dz_dtbeta_r  = diff1z(+1)

     if (shiftadvect) then
        DAr_dtbeta_r = diffadvr(-1)
        DAz_dtbeta_r = diffadvz(+1)
     else
        DAr_dtbeta_r = Dr_dtbeta_r
        DAz_dtbeta_r = Dz_dtbeta_r
     end if

!    Derivatives of dtbeta_z.

     diffvar => dtbeta_z

     Dr_dtbeta_z  = diff1r(+1)
     Dz_dtbeta_z  = diff1z(-1)

     if (shiftadvect) then
        DAr_dtbeta_z = diffadvr(+1)
        DAz_dtbeta_z = diffadvz(-1)
     else
        DAr_dtbeta_z = Dr_dtbeta_z
        DAz_dtbeta_z = Dz_dtbeta_z
     end if

!    Derivatives of dtbeta_p.

     if (angmom) then

        diffvar => dtbeta_p

        Dr_dtbeta_p  = diff1r(+1)
        Dz_dtbeta_p  = diff1z(+1)

        if (shiftadvect) then
           DAr_dtbeta_p = diffadvr(+1)
           DAz_dtbeta_p = diffadvz(+1)
        else
           DAr_dtbeta_p = Dr_dtbeta_p
           DAz_dtbeta_p = Dz_dtbeta_p
        end if

     end if

!    Conformal divergence of the shift vector.
!    Remember that for lagrangian evolutions we must
!    have hdet=1.

     DIV_beta = Dr_beta_r + Dz_beta_z + beta_r/r

     Dr_DIV_beta = Drr_beta_r + Drz_beta_z + DD_beta_rr
     Dz_DIV_beta = Drz_beta_r + Dzz_beta_z + Dz_beta_r/r

     if (bssnflavor=="eulerian") then

        auxarray = half*(beta_r*Dr_hdet + beta_z*Dz_hdet)*ihdet
        DIV_beta = DIV_beta + auxarray

        diffvar => auxarray
        Dr_DIV_beta = Dr_DIV_beta + diff1r(+1)
        Dz_DIV_beta = Dz_DIV_beta + diff1z(+1)

     end if

  end if


! *********************************
! ***   SUBTRACT TRACE OF KTA   ***
! *********************************

! Subtract trace of KTA to make sure that it is traceless.
! Must be done before the regularization variables.

  if (.not.evolveKTH) then
     if (.not.angmom) then
        KTH = - (g_A*KTA + g_B*KTB + two*r**2*g_C*KTC)/g_H
     else
        KTH = - (g_A*KTA + g_B*KTB + two*r**2*(g_C*KTC &
            + r**2*g_C1*KTC1 + g_C2*KTC2))/g_H
     end if
  else
     call detrace_KT
  end if


! *******************
! ***   LAMBDAS   ***
! *******************

! If we don't want to introduce the regularization variables
! lambda and Alambda, then here we calculate them in terms
! of metric and extrinsic curvature.

  lamDef  = (A - H)/r**2
  AlamDef = (KTA - KTH)/r**2

  if (nolambda.or.(time==zero)) then
     lambda  = lamDef
     Alambda = AlamDef
  end if

! Derivatives of lambda.

  diffvar => lambda
  Dr_lambda  = diff1r(+1)
  Dz_lambda  = diff1z(+1)
  Drr_lambda = diff2r(+1)
  Dzz_lambda = diff2z(+1)
  Drz_lambda = diff2rz(+1,+1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_lambda = diffadvr(+1)
        DAz_lambda = diffadvz(+1)
     else
        DAr_lambda = Dr_lambda
        DAz_lambda = Dz_lambda
     end if
  end if

! Derivatives of Alambda.

  diffvar => Alambda
  Dr_Alambda = diff1r(+1)
  Dz_Alambda = diff1z(+1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_Alambda = diffadvr(+1)
        DAz_Alambda = diffadvz(+1)
     else
        DAr_Alambda = Dr_Alambda
        DAz_Alambda = Dz_Alambda
     end if
  end if

! g_lambda := (g_A - g_H)/r**2 = (C**2 - B*lambda - C2**2)*ihdet.

  if (angmom) then
     g_lambda = (C**2 - B*lambda - C2**2)*ihdet
  else
     g_lambda = (C**2 - B*lambda)*ihdet
  end if


! *******************************
! ***   EXTRINSIC CURVATURE   ***
! *******************************

! Derivatives of trK.

  diffvar => trK
  Dr_trK = diff1r(+1)
  Dz_trK = diff1z(+1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_trK = diffadvr(+1)
        DAz_trK = diffadvz(+1)
     else
        DAr_trK = Dr_trK
        DAz_trK = Dz_trK
     end if
  end if

! Derivatives of KTA.

  diffvar => KTA
  Dr_KTA = diff1r(+1)
  Dz_KTA = diff1z(+1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_KTA = diffadvr(+1)
        DAz_KTA = diffadvz(+1)
     else
        DAr_KTA = Dr_KTA
        DAz_KTA = Dz_KTA
     end if
  end if

! Derivatives of KTB.

  diffvar => KTB
  Dr_KTB = diff1r(+1)
  Dz_KTB = diff1z(+1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_KTB = diffadvr(+1)
        DAz_KTB = diffadvz(+1)
     else
        DAr_KTB = Dr_KTB
        DAz_KTB = Dz_KTB
     end if
  end if

! Derivatives of KTH.

  diffvar => KTH
  Dr_KTH = diff1r(+1)
  Dz_KTH = diff1z(+1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_KTH = diffadvr(+1)
        DAz_KTH = diffadvz(+1)
     else
        DAr_KTH = Dr_KTH
        DAz_KTH = Dz_KTH
     end if
  end if

! Derivatives of KTC.

  diffvar => KTC
  Dr_KTC = diff1r(+1)
  Dz_KTC = diff1z(-1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_KTC = diffadvr(+1)
        DAz_KTC = diffadvz(-1)
     else
        DAr_KTC = Dr_KTC
        DAz_KTC = Dz_KTC
     end if
  end if

! Derivatives of KTC1 and KTC2.

  if (angmom) then

     diffvar => KTC1
     Dr_KTC1  = diff1r(+1)
     Dz_KTC1  = diff1z(+1)

     if (shift/="none") then
        if (shiftadvect) then
           DAr_KTC1 = diffadvr(+1)
           DAz_KTC1 = diffadvz(+1)
        else
           DAr_KTC1 = Dr_KTC1
           DAz_KTC1 = Dz_KTC1
        end if
     end if

     diffvar => KTC2
     Dr_KTC2 = diff1r(+1)
     Dz_KTC2 = diff1z(-1)

     if (shift/="none") then
        if (shiftadvect) then
           DAr_KTC2 = diffadvr(+1)
           DAz_KTC2 = diffadvz(-1)
        else
           DAr_KTC2 = Dr_KTC2
           DAz_KTC2 = Dz_KTC2
        end if
     end if

  end if

! Contravariant extrinsic curvature.

  if (angmom) then

     KT_up_rr = (2.0*g_A*g_C1*KTC1 + 2.0*g_C*g_C1*KTC2 + g_C1**2*KTH)*r**4 &
              + (2.0*g_A*g_C*KTC + g_C**2*KTB)*r**2 + g_A**2*KTA
     KT_up_zz = 2.0*g_C*r**4*g_C2*KTC1 + (g_C**2*KTA + 2.0*g_C*g_B*KTC &
              + 2.0*g_B*g_C2*KTC2 + g_C2**2*KTH)*r**2 + g_B**2*KTB

     !***************************************************************
     !The next one is for completeness, since it has leading negative
     ! power of r, and shouldn't appear in equations
     KT_up_pp = (g_C1**2*KTA+2.0*g_C1*g_C2*KTC+2.0*g_C1*g_H*KTC1)*r**2 &
              + 2.0*g_C2*g_H*KTC2+g_C2**2*KTB+g_H**2*KTH/r**2
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     KT_up_rz = g_C*r**5*g_C1*KTC1 + (g_C**2*KTC + g_C2*g_A*KTC1 &
              + g_B*g_C1*KTC2 + g_C2*g_C*KTC2 + g_C2*g_C1*KTH)*r**3 &
              + (g_C*g_A*KTA + g_B*g_A*KTC + g_B*g_C*KTB)*r
     KT_up_rp = r**5*g_C1**2*KTC1 + (g_C1*g_C*KTC + g_C2*g_C1*KTC2)*r**3 &
              + (g_C1*g_A*KTA + g_H*g_A*KTC1 + g_C2*g_A*KTC + g_C2*g_C*KTB &
              + g_H*g_C*KTC2+g_H*g_C1*KTH)*r
     KT_up_zp = r**4*g_C1*g_C2*KTC1 + (g_C1*g_B*KTC + g_H*g_C*KTC1 &
              + g_C2**2*KTC2 + g_C2*g_C*KTC + g_C1*g_C*KTA)*r**2 &
              + g_C2*g_B*KTB + g_H*g_B*KTC2 + g_H*g_C2*KTH

  else

     KT_up_rr = (2.0*g_A*g_C*KTC + g_C**2*KTB)*r**2 + g_A**2*KTA
     KT_up_zz = (g_C**2*KTA + 2.0*g_C*g_B*KTC)*r**2 + g_B**2*KTB
     KT_up_pp = g_H**2/r**2*KTH
     KT_up_rz = g_C**2*r**3*KTC + (g_C*g_A*KTA + g_B*g_A*KTC + g_B*g_C*KTB)*r

  end if

! Squared traceless extrinsic curvature:
!
! KT2 = KT_im gamma^mn KT_nj

  if (angmom) then

     KT2_A = (2.0*KTA*g_C1*KTC1 + 2.0*KTC*g_C2*KTC1 + KTC1**2*g_H)*r**4 &
           + (2.0*KTA*g_C*KTC + KTC**2*g_B)*r**2 + g_A*KTA**2

     KT2_B = 2.0*KTC*r**4*g_C1*KTC2 + (KTC**2*g_A + 2.0*KTC*g_C*KTB &
           +2.0*KTB*g_C2*KTC2 + KTC2**2*g_H)*r**2 + g_B*KTB**2
     
     KT2_H = (2.0*KTC1*g_C1*KTH + KTC1**2*g_A + 2.0*KTC1*g_C*KTC2)*r**4 &
           + (2.0*KTC2*g_C2*KTH + KTC2**2*g_B)*r**2 + KTH**2*g_H

     KT2_C = KTC1*g_C1*KTC*r**4 + (KTA*g_C1*KTC2 + KTC**2*g_C &
           + KTC*g_C2*KTC2+KTC1*g_C2*KTB + KTC1*g_H*KTC2)*r**2 &
           + KTA*g_A*KTC + KTA*g_C*KTB+KTC*g_B*KTB

     KT2_C1 = KTC1**2*g_C1*r**4 + (KTC*g_C*KTC1 + KTC1*g_C2*KTC2)*r**2 &
            + KTA*g_A*KTC1 + KTA*g_C*KTC2 + KTA*g_C1*KTH + KTC*g_B*KTC2 &
            + KTC*g_C2*KTH + KTC1*g_H*KTH

     KT2_C2 = KTC2*g_C1*KTC1*r**4 + (KTC*g_C1*KTH + KTB*g_C*KTC1 &
            + KTC*g_C*KTC2 + KTC2**2*g_C2 + KTC*g_A*KTC1)*r**2 &
            + KTB*g_B*KTC2 + KTB*g_C2*KTH + KTC2*g_H*KTH

     KT2_lambda = (2.0*KTC*g_C2*KTC1 + KTC1**2*g_H + 2.0*KTA*g_C1*KTC1 &
                - 2.0*KTC1*g_C1*KTH - KTC1**2*g_A - 2.0*KTC1*g_C*KTC2)*r**2 &
                + 2.0*KTA*g_C*KTC+KTC**2*g_B - 2.0*KTC2*g_C2*KTH - KTC2**2*g_B &
                +       ft3*(g_lambda*KTA**2 + g_H*Alambda*(KTA + KTH)) &
                + (1.0-ft3)*(g_lambda*KTH**2 + g_A*Alambda*(KTA + KTH))

!    KT2 = KT^ij KT_ij

     KT2 = g_A*KT2_A + g_B*KT2_B + g_H*KT2_H &
         + two*r**2*(g_C*KT2_C + r**2*g_C1*KT2_C1 + g_C2*KT2_C2)

  else

     KT2_A = (2.0*KTA*g_C*KTC + KTC**2*g_B)*r**2 + g_A*KTA**2

     KT2_B = (KTC**2*g_A + 2.0*KTC*g_C*KTB)*r**2 + g_B*KTB**2

     KT2_H = KTH**2*g_H

     KT2_C = KTC**2*g_C*r**2 + KTA*g_A*KTC + KTA*g_C*KTB + KTC*g_B*KTB

     KT2_lambda = 2.0*KTA*g_C*KTC + KTC**2*g_B &
                +       ft3*(g_lambda*KTA**2 + Alambda*(KTA+KTH)*g_H) &
                + (1.0-ft3)*(g_lambda*KTH**2 + Alambda*(KTA+KTH)*g_A)

!    KT2 = KT^ij KT_ij

     KT2 = g_A*KT2_A + g_B*KT2_B + g_H*KT2_H + two*g_C*KT2_C*r**2

  end if

! For the components of the physical extrinsic curvature
! we must undo the conformal rescaling and add the trace.

  KA = psi4*(KTA + third*A*trK)
  KB = psi4*(KTB + third*B*trK)
  KH = psi4*(KTH + third*H*trK)
  KC = psi4*(KTC + third*C*trK)

  if (angmom) then
     KC1 = psi4*(KTC1 + third*C1*trK)
     KC2 = psi4*(KTC2 + third*C2*trK)
  end if

! Trace of the full squared extrinsic curvature.

  K2 = KT2 + third*(trK + two*z4theta)**2


! **********************
! ***   BSSN Delta   ***
! **********************

! If we don't want to evolve the BSSN variable Delta, then
! here we calculate it in terms of derivatives of the metric.

  DelDef_r = - Dr_g_A - r*Dz_g_C - r*g_lambda - half*(g_A*Dr_hdet + r*g_C*Dz_hdet)*ihdet
  DelDef_z = - Dr_g_C*r - Dz_g_B - two*g_C - half*(r*g_C*Dr_hdet + g_B*Dz_hdet)*ihdet

  if (noDelta_r.or.(time==zero)) then
     Delta_r = DelDef_r
  end if

  if (noDelta_z.or.(time==zero)) then
     Delta_z = DelDef_z
  end if

  if (angmom) then
     DelDef_p = - Dr_g_C1*r - Dz_g_C2 - four*g_C1 - half*(r*g_C1*Dr_hdet + g_C2*Dz_hdet)*ihdet
     if (noDelta_p.or.(time==zero)) then
        Delta_p = DelDef_p
     end if
  end if

! Save initial values of Delta's (for shift conditions).

  if (time==zero) then
     Delta0_r = Delta_r
     Delta0_z = Delta_z
     if (angmom) then
        Delta0_p = Delta_p
     end if
  end if

! Derivatives of Delta_r.

  diffvar => Delta_r
  Dr_Delta_r = diff1r(-1)
  Dz_Delta_r = diff1z(+1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_Delta_r = diffadvr(-1)
        DAz_Delta_r = diffadvz(+1)
     else
        DAr_Delta_r = Dr_Delta_r
        DAz_Delta_r = Dz_Delta_r
     end if
  end if

! Derivatives of Delta_z.

  diffvar => Delta_z
  Dr_Delta_z = diff1r(+1)
  Dz_Delta_z = diff1z(-1)

  if (shift/="none") then
     if (shiftadvect) then
        DAr_Delta_z = diffadvr(+1)
        DAz_Delta_z = diffadvz(-1)
     else
        DAr_Delta_z = Dr_Delta_z
        DAz_Delta_z = Dz_Delta_z
     end if
  end if

! Derivatives of Delta_p.

  if (angmom) then

     diffvar => Delta_p
     Dr_Delta_p = diff1r(+1)
     Dz_Delta_p = diff1z(+1)

     if (shift/="none") then
        if (shiftadvect) then
           DAr_Delta_p = diffadvr(+1)
           DAz_Delta_p = diffadvz(+1)
        else
           DAr_Delta_p = Dr_Delta_p
           DAz_Delta_p = Dz_Delta_p
        end if
     end if

  end if

!                        r
! DD_Delta_rr = d ( Delta / r )
!               r

  auxarray = Delta_r/r

  diffvar => auxarray
  DD_Delta_rr = diff1r(+1)


! ********************
! ***   Z4 THETA   ***
! ********************

  if (formulation=="z4c") then

!    Derivatives of z4theta.

     diffvar => z4theta
     Dr_z4theta = diff1r(+1)
     Dz_z4theta = diff1z(+1)

     if (shift/="none") then
        if (shiftadvect) then
           DAr_z4theta = diffadvr(+1)
           DAz_z4theta = diffadvz(+1)
        else
           DAr_z4theta = Dr_z4theta
           DAz_z4theta = Dz_z4theta
        end if
     end if

  end if


! **********************************************
! ***   COVARIANT DERIVATIVES OF THE LAPSE   ***
! **********************************************

! Calculate second covariant derivative of lapse with both
! indices down and divided by corresponding r power.

! Covariant derivatives (conformal metric connection).

  if (angmom) then

     D2cov_alpha_A = Drr_alpha-Dr_alpha*r**4*g_C1*Dr_C1-(3.0*Dr_alpha&
           *g_C1*C1+Dz_alpha*g_C2*Dr_C1)*r**3-(Dr_alpha*g_C*Dr_C+3.0&
           *Dz_alpha*g_C2*C1)*r**2-(Dr_alpha*g_C*C-half*Dr_alpha*g_C&
           *Dz_A+half*Dz_alpha*g_C*Dr_A+Dz_alpha*g_B*Dr_C)*r-half&
           *Dr_alpha*g_A*Dr_A+half*Dz_alpha*g_B*Dz_A-Dz_alpha*g_B*C

     D2cov_alpha_B = Dzz_alpha-Dr_alpha*r**3*g_C1*Dz_C2-(Dz_alpha&
           *g_C2*Dz_C2+Dz_alpha*g_C*Dz_C)*r**2-(Dr_alpha*g_A*Dz_C&
           -half*Dz_alpha*g_C*Dr_B+half*Dr_alpha*g_C*Dz_B)*r+half&
           *Dr_alpha*g_A*Dr_B-half*Dz_alpha*g_B*Dz_B

     D2cov_alpha_H = -(-half*Dr_alpha*g_C*Dz_H-half*Dz_alpha*g_C&
           *Dr_H)*r+half*Dr_alpha*g_A*Dr_H+Dz_alpha*g_C*H+half&
           *Dz_alpha*g_B*Dz_H+Dr_alpha*g_A*H/r

     D2cov_alpha_C = Drz_alpha/r-half*Dr_alpha*g_C1*Dz_C1*r**3-(half&
           *Dz_alpha*g_C2*Dz_C1+half*Dr_alpha*g_C1*Dr_C2)*r**2&
           -(Dr_alpha*g_C1*C2+half*Dz_alpha*g_C2*Dr_C2)*r-half&
           *Dr_alpha*g_C*Dr_B-Dz_alpha*g_C2*C2-half*Dz_alpha*g_C*Dz_A&
           -(half*Dr_alpha*g_A*Dz_A+half*Dz_alpha*g_B*Dr_B)/r     

     D2cov_alpha_C1 = half*Dr_alpha*g_C*Dz_C1*r-half*Dr_alpha*g_C1&
           *Dr_H+half*Dz_alpha*g_B*Dz_C1-half*Dr_alpha*g_C*Dr_C2&
           -(Dr_alpha*g_C*C2+half*Dz_alpha*g_B*Dr_C2+Dr_alpha*g_C1*H&
           +half*Dz_alpha*g_C2*Dr_H)/r+C1*g_C*Dz_alpha

     D2cov_alpha_C2 = -half*Dz_alpha*g_C*Dz_C1*r**2-(half*Dr_alpha&
           *g_A*Dz_C1-half*Dz_alpha*g_C*Dr_C2+half*Dr_alpha*g_C1&
           *Dz_H)*r+half*Dr_alpha*g_A*Dr_C2-half*Dz_alpha*g_C2*Dz_H&
           +Dz_alpha*g_C*C2+Dr_alpha*g_A*C2/r
                  
     D2cov_alpha_lambda = -Dr_alpha*g_C1*Dr_C1*r**2 - (3.0*Dr_alpha&
           *g_C1*C1+Dz_alpha*g_C2*Dr_C1)*r-Dr_alpha*g_C*Dr_C-3.0&
           *Dz_alpha*g_C2*C1-(Dr_alpha*g_C*C-half*Dr_alpha*g_C*Dz_A&
           +half*Dz_alpha*g_C*Dr_A+Dz_alpha*g_B*Dr_C+half*Dr_alpha&
           *g_C*Dz_H+half*Dz_alpha*g_C*Dr_H)/r+DD_alphar/r+Dr_alpha/r&
           *(ft4*(r**2*C1*g_C1+C2*g_C2-g_lambda*H)+(1.0-ft4)*(r&
           **2*C1*g_C1+C*g_C+g_A*lambda))-half*Dr_alpha*g_A*(Dr_A&
           +Dr_H)/r**2+half*Dz_alpha*g_B*Dz_lambda+Dz_alpha*(g_C&
           *lambda+C1*g_C2)

  else

     D2cov_alpha_A = Drr_alpha - Dr_alpha*g_C*r**2*Dr_C - (Dz_alpha*g_B*Dr_C &
           + half*Dz_alpha*g_C*Dr_A + Dr_alpha*g_C*C - half*Dr_alpha*g_C*Dz_A)*r &
           - half*Dr_alpha*g_A*Dr_A - Dz_alpha*g_B*C + half*Dz_alpha*g_B*Dz_A

     D2cov_alpha_B = Dzz_alpha - Dz_alpha*g_C*r**2*Dz_C - (Dr_alpha*g_A*Dz_C &
           + half*Dr_alpha*g_C*Dz_B - half*Dz_alpha*g_C*Dr_B)*r &
           + half*Dr_alpha*g_A*Dr_B - half*Dz_alpha*g_B*Dz_B

     D2cov_alpha_H = -(-half*Dr_alpha*g_C*Dz_H - half*Dz_alpha*g_C*Dr_H)*r &
           + half*Dr_alpha*g_A*Dr_H + Dz_alpha*g_C*H &
           + half*Dz_alpha*g_B*Dz_H + Dr_alpha*g_A*H/r

     D2cov_alpha_C = Drz_alpha/r - half*Dr_alpha*g_C*Dr_B - half*Dz_alpha*g_C*Dz_A &
           - (half*Dr_alpha*g_A*Dz_A + half*Dz_alpha*g_B*Dr_B)/r
     
     D2cov_alpha_lambda = (DD_alphar)/r + Dr_alpha/r*(-ft4*g_lambda*H &
           + (1.0-ft4)*(g_C*C + g_A*lambda)) - Dr_alpha*g_C*Dr_C &
           - (Dz_alpha*g_B*Dr_C + Dr_alpha*g_C*C - half*Dr_alpha*g_C*Dz_A &
           + half*Dz_alpha*g_C*Dr_A + half*Dr_alpha*g_C*Dz_H &
           + half*Dz_alpha*g_C*Dr_H)/r + half*Dz_alpha*g_B*Dz_lambda &
           - half*Dr_alpha*g_A*(Dr_A + Dr_H)/r**2 + Dz_alpha*lambda*g_C

  end if

! Add the derivatives that complete the physical metric conection.

  auxarray = g_A*Dr_alpha*Dr_phi + g_B*Dz_alpha*Dz_phi &
           + r*g_C*(Dr_alpha*Dz_phi + Dz_alpha*Dr_phi)

  D2cov_alpha_A = D2cov_alpha_A - two*(two*Dr_alpha*Dr_phi - A*auxarray)
  D2cov_alpha_B = D2cov_alpha_B - two*(two*Dz_alpha*Dz_phi - B*auxarray)
  D2cov_alpha_H = D2cov_alpha_H + two*H*auxarray
  D2cov_alpha_C = D2cov_alpha_C - two*((Dr_alpha/r)*Dz_phi &
                + Dz_alpha*(Dr_phi/r) - C*auxarray)

  D2cov_alpha_lambda = D2cov_alpha_lambda - two*(two*(Dr_alpha/r)*(Dr_phi/r) &
                - lambda*auxarray)

  if (angmom) then
     D2cov_alpha_C1 = D2cov_alpha_C1 + two*C1*auxarray
     D2cov_alpha_C2 = D2cov_alpha_C2 + two*C2*auxarray
  end if

! Laplacian of lapse.

  if (angmom) then
     Lapla_alpha = g_A*D2cov_alpha_A + g_B*D2cov_alpha_B + g_H*D2cov_alpha_H &
          + two*r**2*(g_C*D2cov_alpha_C + r**2*g_C1*D2cov_alpha_C1 &
          + g_C2*D2cov_alpha_C2)
  else
     Lapla_alpha = g_A*D2cov_alpha_A + g_B*D2cov_alpha_B + g_H*D2cov_alpha_H &
          + two*r**2*g_C*D2cov_alpha_C
  end if

! Divide by conformal factor.

  Lapla_alpha = Lapla_alpha/psi4

! Below is the flat Laplacian (commented out).  It is useful for testing.

! Lapla_alpha = Drr_alpha + Dzz_alpha + Dr_alpha/r


! ************************
! ***   RICCI TENSOR   ***
! ************************

! The Ricci tensor has a contribution from the conformal factor and
! another coming from the conformal metric.
!
! R_ij = Rphi_ij + \hat R_ij
!
! Components of the Ricci tensor found with MAPLE.
!
! Terms coming from the derivatives of the conformal factor
!
!                                 2                                k
! Rphi  = - 2 D D phi - 2 gamma  D phi + 4 D phi D phi - 4 gamma  D phi D phi
!     ij       i j             ij           i     j             ij       k
!
! where D is the covariant conformal derivative.

! I) Terms coming from second covariant derivatives of phi.

  if (angmom) then

     RIC_A = - two*(Drr_phi-Dr_phi*r**4*g_C1*Dr_C1-(3.0*Dr_phi&
           *g_C1*C1+Dz_phi*g_C2*Dr_C1)*r**3-(Dr_phi*g_C*Dr_C+3.0&
           *Dz_phi*g_C2*C1)*r**2-(Dr_phi*g_C*C-half*Dr_phi*g_C&
           *Dz_A+half*Dz_phi*g_C*Dr_A+Dz_phi*g_B*Dr_C)*r-half&
           *Dr_phi*g_A*Dr_A+half*Dz_phi*g_B*Dz_A-Dz_phi*g_B*C)

     RIC_B = - two*(Dzz_phi-Dr_phi*r**3*g_C1*Dz_C2-(Dz_phi&
           *g_C2*Dz_C2+Dz_phi*g_C*Dz_C)*r**2-(Dr_phi*g_A*Dz_C&
           -half*Dz_phi*g_C*Dr_B+half*Dr_phi*g_C*Dz_B)*r+half&
           *Dr_phi*g_A*Dr_B-half*Dz_phi*g_B*Dz_B)

     RIC_H = - two*(-(-half*Dr_phi*g_C*Dz_H-half*Dz_phi*g_C&
           *Dr_H)*r+half*Dr_phi*g_A*Dr_H+Dz_phi*g_C*H+half&
           *Dz_phi*g_B*Dz_H+Dr_phi*g_A*H/r)

     RIC_C = - two*(Drz_phi/r-half*Dr_phi*g_C1*Dz_C1*r**3-(half&
           *Dz_phi*g_C2*Dz_C1+half*Dr_phi*g_C1*Dr_C2)*r**2&
           -(Dr_phi*g_C1*C2+half*Dz_phi*g_C2*Dr_C2)*r-half&
           *Dr_phi*g_C*Dr_B-Dz_phi*g_C2*C2-half*Dz_phi*g_C*Dz_A&
           -(half*Dr_phi*g_A*Dz_A+half*Dz_phi*g_B*Dr_B)/r )

     RIC_C1 = - two*(half*Dr_phi*g_C*Dz_C1*r-half*Dr_phi*g_C1&
           *Dr_H+half*Dz_phi*g_B*Dz_C1-half*Dr_phi*g_C*Dr_C2&
           -(Dr_phi*g_C*C2+half*Dz_phi*g_B*Dr_C2+Dr_phi*g_C1*H&
           +half*Dz_phi*g_C2*Dr_H)/r+C1*g_C*Dz_phi)

     RIC_C2 = - two*(-half*Dz_phi*g_C*Dz_C1*r**2-(half*Dr_phi&
           *g_A*Dz_C1-half*Dz_phi*g_C*Dr_C2+half*Dr_phi*g_C1&
           *Dz_H)*r+half*Dr_phi*g_A*Dr_C2-half*Dz_phi*g_C2*Dz_H&
           +Dz_phi*g_C*C2+Dr_phi*g_A*C2/r)

     RIC_lambda = - two*(-Dr_phi*g_C1*Dr_C1*r**2-(3.0*Dr_phi&
           *g_C1*C1+Dz_phi*g_C2*Dr_C1)*r-Dr_phi*g_C*Dr_C-3.0&
           *Dz_phi*g_C2*C1-(Dr_phi*g_C*C-half*Dr_phi*g_C*Dz_A&
           +half*Dz_phi*g_C*Dr_A+Dz_phi*g_B*Dr_C+half*Dr_phi&
           *g_C*Dz_H+half*Dz_phi*g_C*Dr_H)/r+DD_phir/r+Dr_phi/r&
           *(ft5*(r**2*C1*g_C1+C2*g_C2-g_lambda*H)+(1.0-ft5)*(r&
           **2*C1*g_C1+C*g_C+g_A*lambda))-half*Dr_phi*g_A*(Dr_A&
           +Dr_H)/r**2+half*Dz_phi*g_B*Dz_lambda+Dz_phi*(g_C&
           *lambda+C1*g_C2))
  else

     RIC_A = - two*(Drr_phi - Dr_phi*g_C*r**2*Dr_C - (Dz_phi*g_B*Dr_C &
           + half*Dz_phi*g_C*Dr_A + Dr_phi*g_C*C - half*Dr_phi*g_C*Dz_A)*r &
           - half*Dr_phi*g_A*Dr_A - Dz_phi*g_B*C + half*Dz_phi*g_B*Dz_A)

     RIC_B = - two*(Dzz_phi - Dz_phi*g_C*r**2*Dz_C - (Dr_phi*g_A*Dz_C &
           + half*Dr_phi*g_C*Dz_B - half*Dz_phi*g_C*Dr_B)*r &
           + half*Dr_phi*g_A*Dr_B - half*Dz_phi*g_B*Dz_B)

     RIC_H = - two*(-(-half*Dr_phi*g_C*Dz_H - half*Dz_phi*g_C*Dr_H)*r &
           + half*Dr_phi*g_A*Dr_H + Dz_phi*g_C*H &
           + half*Dz_phi*g_B*Dz_H + Dr_phi*g_A*H/r)

     RIC_C = - two*(Drz_phi/r - half*Dr_phi*g_C*Dr_B - half*Dz_phi*g_C*Dz_A &
           - (half*Dr_phi*g_A*Dz_A + half*Dz_phi*g_B*Dr_B)/r)
     
     RIC_lambda = - two*((DD_phir)/r + Dr_phi/r*(-ft5*g_lambda*H &
           + (1.0-ft5)*(g_C*C + g_A*lambda)) - Dr_phi*g_C*Dr_C &
           - (Dz_phi*g_B*Dr_C + Dr_phi*g_C*C - half*Dr_phi*g_C*Dz_A &
           + half*Dz_phi*g_C*Dr_A + half*Dr_phi*g_C*Dz_H &
           + half*Dz_phi*g_C*Dr_H)/r + half*Dz_phi*g_B*Dz_lambda &
           - half*Dr_phi*g_A*(Dr_A/r + Dr_H/r)/r + Dz_phi*lambda*g_C)

  end if

! II) Terms coming from Laplacian of phi.

! Since by now RIC_* just stores the 2nd covariant derivative (times -2)
! we take the trace and assign to auxarray to get the conformal laplacian
! (times -2).

  if (angmom) then
     auxarray = g_A*RIC_A + g_B*RIC_B + g_H*RIC_H &
              + two*r**2*(g_C*RIC_C + r**2*g_C1*RIC_C1 + g_C2*RIC_C2)
  else
     auxarray = g_A*RIC_A + g_B*RIC_B + g_H*RIC_H + two*r**2*g_C*RIC_C
  end if

  RIC_A = RIC_A + A*auxarray
  RIC_B = RIC_B + B*auxarray
  RIC_H = RIC_H + H*auxarray
  RIC_C = RIC_C + C*auxarray

  RIC_lambda = RIC_lambda + lambda*auxarray

  if (angmom) then
     RIC_C1 = RIC_C1 + C1*auxarray 
     RIC_C2 = RIC_C2 + C2*auxarray
  end if

! III) Quadratic derivatives of phi.

  RIC_A = RIC_A + four*(Dr_phi**2)
  RIC_B = RIC_B + four*(Dz_phi**2)
! RIC_H = RIC_H + four*(zero)/r**2
  RIC_C = RIC_C + four*(Dr_phi*Dz_phi)/r

  RIC_lambda = RIC_lambda + four*(Dr_phi/r)**2
! RIC_C1 = RIC_C2 + four*(zero)/r**2
! RIC_C2 = RIC_C2 + four*(zero)/r**2

! IV) Contracted quadratic derivatives (times -2 for convenience).

  auxarray = - two*(g_A*Dr_phi*Dr_phi + g_B*Dz_phi*Dz_phi &
           + two*r*g_C*Dr_phi*Dz_phi)

  RIC_A = RIC_A + two*A*auxarray
  RIC_B = RIC_B + two*B*auxarray
  RIC_H = RIC_H + two*H*auxarray
  RIC_C = RIC_C + two*C*auxarray

  RIC_lambda = RIC_lambda + two*lambda*auxarray

  if (angmom) then
     RIC_C1 = RIC_C1 + two*C1*auxarray 
     RIC_C2 = RIC_C2 + two*C2*auxarray
  end if

! V) Terms coming from the conformal Ricci Tensor.
! We calculate them on a separate subroutines.
! There are four contributions, calculated on
! separate subroutines each.

  call calc_conformalRicci1
  call calc_conformalRicci2
  call calc_conformalRicci3
  call calc_conformalRicci4

! Ricci scalar.

  if (angmom) then
     RSCAL = g_A*RIC_A + g_B*RIC_B + g_H*RIC_H &
           + two*r**2*(g_C*RIC_C + r**2*g_C1*RIC_C1 + g_C2*RIC_C2)
  else
     RSCAL = g_A*RIC_A + g_B*RIC_B + g_H*RIC_H + two*r**2*g_C*RIC_C
  end if

! Conformal rescaling of the Ricci scalar since it's calculated
! contracting with the full physical metric.

  RSCAL = RSCAL/psi4


! ****************************
! ***   DIVERGENCE OF KT   ***
! ****************************

! The divergence of KTA is calculated in a separate routine.

  call calc_DIVKTA


! **************************************
! ***   CONFORMAL SPHERICAL METRIC   ***
! **************************************

! Components of conformal metric in spherical coordinates
! (r,theta,phi).  The coordinate transformation is
! (remember theta is measured from the axis):
!
! r = rr sin(theta)
! z = rr cos(theta)
!
! Notice that phi is already one of our coordinates
! as we are using cylindrical coordinates (r,z,phi).

  grr = (r**2*A + z**2*B + two*r**2*z*C)/rr**2
  gtt = (z**2*A + r**2*B - two*r**2*z*C)
  grt = r*(z*(A-B) + (z**2 - r**2)*C)/rr
  gpp = r**2*H

  if (angmom) then
     grp = r**2*(r**2*C1 + z*C2)/rr
     gtp = r**3*(z*C1 - C2)
  else
     grp = 0.d0
     gtp = 0.d0
  end if

! Determinant of spherical metric.

  if (angmom) then
     hdetsph = grr*gtt*gpp + two*grt*grp*gtp &
             - grr*gtp**2 - gtt*grp**2 - gpp*grt**2
  else
     auxarray = (grr*gtt - grt**2)
     hdetsph  = gpp*auxarray
  end if

! Inverse conformal spherical metric.

  if (angmom) then
     ginvrr = (gtt*gpp - gtp**2)/hdetsph
     ginvtt = (grr*gpp - grp**2)/hdetsph
     ginvpp = (grr*gtt - grt**2)/hdetsph
     ginvrt = (grp*gtp - grt*gpp)/hdetsph
     ginvrp = (grt*gtp - grp*gtt)/hdetsph
     ginvtp = (grt*grp - gtp*grr)/hdetsph
  else
     ginvrr = gtt/auxarray
     ginvtt = grr/auxarray
     ginvpp = one/gpp
     ginvrt = - grt/auxarray
     ginvrp = zero
     ginvtp = zero
  end if

! Check if the inversion worked.  Should be commented out!
!
!  do i=0,Nr
!     do j=0,Nz
!        write(*,'(2I4,9F7.2)') i,j,grr(i,j)*ginvrr(i,j) + grt(i,j)*ginvrt(i,j) + grp(i,j)*ginvrp(i,j), &
!                                   grr(i,j)*ginvrt(i,j) + grt(i,j)*ginvtt(i,j) + grp(i,j)*ginvtp(i,j), &
!                                   grr(i,j)*ginvrp(i,j) + grt(i,j)*ginvtp(i,j) + grp(i,j)*ginvpp(i,j), &
!                                   grt(i,j)*ginvrr(i,j) + gtt(i,j)*ginvrt(i,j) + gtp(i,j)*ginvrp(i,j), &
!                                   grt(i,j)*ginvrt(i,j) + gtt(i,j)*ginvtt(i,j) + gtp(i,j)*ginvtp(i,j), &
!                                   grt(i,j)*ginvrp(i,j) + gtt(i,j)*ginvtp(i,j) + gtp(i,j)*ginvpp(i,j), &
!                                   grp(i,j)*ginvrr(i,j) + gtp(i,j)*ginvrt(i,j) + gpp(i,j)*ginvrp(i,j), &
!                                   grp(i,j)*ginvrt(i,j) + gtp(i,j)*ginvtt(i,j) + gpp(i,j)*ginvtp(i,j), &
!                                   grp(i,j)*ginvrp(i,j) + gtp(i,j)*ginvtp(i,j) + gpp(i,j)*ginvpp(i,j)
!     end do
!  end do


! ***************
! ***   END   ***
! ***************

  end subroutine auxiliary_geometry

