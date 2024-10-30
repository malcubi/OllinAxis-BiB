!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/matter/complex.f90,v 1.5 2021/03/05 23:02:31 malcubi Exp $ 

  subroutine sources_complex

! ********************************************
! ***   SOURCES FOR COMPLEX SCALAR FIELD   ***
! ********************************************

! This routine calculates the sources for a real
! scalar field.  The Klein-Gordon equation with
! an arbitrary potential V(phi) has the form:
!
! Box (phi)  -  (dV/dphi)  =  0
!
! The scalar field is evolved in first order form.
! Remember that we have defined the space and time
! derivative arrays as:
!                        m
! pi  =  ( d phi  -  beta d phi ) / alpha
!           t              m
!
! xi  =  d phi
!   i     i

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

! Source for scalar field.

  scomplex_phiR = alpha*complex_piR
  scomplex_phiI = alpha*complex_piI

  if (shift/="none") then
     scomplex_phiR = scomplex_phiR + beta_r*complex_xiR_r + beta_z*complex_xiR_z
     scomplex_phiI = scomplex_phiI + beta_r*complex_xiI_r + beta_z*complex_xiI_z
  end if

! Source for r derivatives.

  scomplex_xiR_r = alpha*Dr_complex_piR + complex_piR*Dr_alpha
  scomplex_xiI_r = alpha*Dr_complex_piI + complex_piI*Dr_alpha

  if (shift/="none") then
     scomplex_xiR_r = scomplex_xiR_r &
                    + beta_r*DAr_complex_xiR_r + complex_xiR_r*Dr_beta_r &
                    + beta_z*DAz_complex_xiR_r + complex_xiR_z*Dr_beta_z
     scomplex_xiI_r = scomplex_xiI_r &
                    + beta_r*DAr_complex_xiI_r + complex_xiI_r*Dr_beta_r &
                    + beta_z*DAz_complex_xiI_r + complex_xiI_z*Dr_beta_z
  end if

! Source for z derivatives.

  scomplex_xiR_z = alpha*Dz_complex_piR + complex_piR*Dz_alpha
  scomplex_xiI_z = alpha*Dz_complex_piI + complex_piI*Dz_alpha

  if (shift/="none") then
     scomplex_xiR_z = scomplex_xiR_z &
                    + beta_r*DAr_complex_xiR_z + complex_xiR_r*Dz_beta_r &
                    + beta_z*DAz_complex_xiR_z + complex_xiR_z*Dz_beta_z
     scomplex_xiI_z = scomplex_xiI_z &
                    + beta_r*DAr_complex_xiI_z + complex_xiI_r*Dz_beta_r &
                    + beta_z*DAz_complex_xiI_z + complex_xiI_z*Dz_beta_z
  end if

! Source for time derivative. There are many terms,
! so we build the sources step by step.

! Derivatives of xi (second derivatives of phi).

  if (complexmethod == "first") then

!    Terms with derivatives of xi.

     scomplex_piR = alpha*(g_A*Dr_complex_xiR_r + g_B*Dz_complex_xiR_z &
                  + r*g_C*(Dr_complex_xiR_z+Dz_complex_xiR_r))
     scomplex_piI = alpha*(g_A*Dr_complex_xiI_r + g_B*Dz_complex_xiI_z &
                  + r*g_C*(Dr_complex_xiI_z+Dz_complex_xiI_r))

!    Terms with derivatives of conformal factor.

     scomplex_piR = scomplex_piR + two*alpha &
                  *(Dr_phi*(g_A*complex_xiR_r + r*g_C*complex_xiR_z) &
                  + Dz_phi*(g_B*complex_xiR_z + r*g_C*complex_xiR_r))
     scomplex_piI = scomplex_piI + two*alpha &
                  *(Dr_phi*(g_A*complex_xiI_r + r*g_C*complex_xiI_z) &
                  + Dz_phi*(g_B*complex_xiI_z + r*g_C*complex_xiI_r))

!    Terms with derivatives of alpha.

     scomplex_piR = scomplex_piR &
                  + Dr_alpha*(g_A*complex_xiR_r + r*g_C*complex_xiR_z) &
                  + Dz_alpha*(g_B*complex_xiR_z + r*g_C*complex_xiR_r)
     scomplex_piI = scomplex_piI &
                  + Dr_alpha*(g_A*complex_xiI_r + r*g_C*complex_xiI_z) &
                  + Dz_alpha*(g_B*complex_xiI_z + r*g_C*complex_xiI_r)

!    Terms coming from derivatives of inverse conformal metric.

     scomplex_piR = scomplex_piR + alpha &
                  *(complex_xiR_r*(Dr_g_A + r*Dz_g_C) + complex_xiR_z*(Dz_g_B + r*Dr_g_C + g_C))
     scomplex_piI = scomplex_piI + alpha &
                  *(complex_xiI_r*(Dr_g_A + r*Dz_g_C) + complex_xiI_z*(Dz_g_B + r*Dr_g_C + g_C))

!    Terms coming from derivatives of the determinant of the
!    conformal metric g. Remember that det(g) = r**2*h.

     scomplex_piR = scomplex_piR + alpha*(g_A*complex_xiR_r/r + g_C*complex_xiR_z)
     scomplex_piI = scomplex_piI + alpha*(g_A*complex_xiI_r/r + g_C*complex_xiI_z)

     scomplex_piR = scomplex_piR + half*alpha*ihdet &
                  *(Dr_hdet*(g_A*complex_xiR_r + r*g_C*complex_xiR_z) &
                  + Dz_hdet*(g_B*complex_xiR_z + r*g_C*complex_xiR_r))
     scomplex_piI = scomplex_piI + half*alpha*ihdet &
                  *(Dr_hdet*(g_A*complex_xiI_r + r*g_C*complex_xiI_z) &
                  + Dz_hdet*(g_B*complex_xiI_z + r*g_C*complex_xiI_r))

  else

!    Terms with second derivatives of phi.

     scomplex_piR = alpha*(g_A*Drr_complex_phiR + g_B*Dzz_complex_phiR &
                  + two*r*g_C*Drz_complex_phiR)
     scomplex_piI = alpha*(g_A*Drr_complex_phiI + g_B*Dzz_complex_phiI &
                  + two*r*g_C*Drz_complex_phiI)

!    Terms with derivatives of conformal factor.

     scomplex_piR = scomplex_piR + two*alpha &
                  *(Dr_phi*(g_A*Dr_complex_phiR + r*g_C*Dz_complex_phiR) &
                  + Dz_phi*(g_B*Dz_complex_phiR + r*g_C*Dr_complex_phiR))
     scomplex_piI = scomplex_piI + two*alpha &
                  *(Dr_phi*(g_A*Dr_complex_phiI + r*g_C*Dz_complex_phiI) &
                  + Dz_phi*(g_B*Dr_complex_phiI + r*g_C*Dz_complex_phiI))

!    Terms with derivatives of alpha.

     scomplex_piR = scomplex_piR &
                  + Dr_alpha*(g_A*Dr_complex_phiR + r*g_C*Dz_complex_phiR) &
                  + Dz_alpha*(g_B*Dz_complex_phiR + r*g_C*Dr_complex_phiR)
     scomplex_piI = scomplex_piI &
                  + Dr_alpha*(g_A*Dr_complex_phiI + r*g_C*Dz_complex_phiI) &
                  + Dz_alpha*(g_B*Dz_complex_phiI + r*g_C*Dr_complex_phiI)

!    Terms coming from derivatives of inverse conformal metric.

     scomplex_piR = scomplex_piR + alpha &
                  *(Dr_complex_phiR*(Dr_g_A + r*Dz_g_C) + Dz_complex_phiR*(Dz_g_B + r*Dr_g_C + g_C))
     scomplex_piI = scomplex_piI + alpha &
                  *(Dr_complex_phiI*(Dr_g_A + r*Dz_g_C) + Dz_complex_phiI*(Dz_g_B + r*Dr_g_C + g_C))

!    Terms coming from derivatives of the determinant of the
!    conformal metric g. Remember that det(g) = r**2*h.

     scomplex_piR = scomplex_piR + alpha*(g_A*Dr_complex_phiR/r + g_C*Dz_complex_phiR)
     scomplex_piI = scomplex_piI + alpha*(g_A*Dr_complex_phiI/r + g_C*Dz_complex_phiI)

     scomplex_piR = scomplex_piR + half*alpha*ihdet &
                  *(Dr_hdet*(g_A*Dr_complex_phiR + r*g_C*Dz_complex_phiR) &
                  + Dz_hdet*(g_B*Dz_complex_phiR + r*g_C*Dr_complex_phiR))
     scomplex_piI = scomplex_piI + half*alpha*ihdet &
                  *(Dr_hdet*(g_A*Dr_complex_phiI + r*g_C*Dz_complex_phiI) &
                  + Dz_hdet*(g_B*Dz_complex_phiI + r*g_C*Dr_complex_phiI))

  end if

! Rescale by conformal factor since all the above terms
! were calculated with the conformal metric.

  scomplex_piR = scomplex_piR/psi4
  scomplex_piI = scomplex_piI/psi4

! Term proportional to trK.

  scomplex_piR = scomplex_piR + alpha*complex_piR*trK
  scomplex_piI = scomplex_piI + alpha*complex_piI*trK

! Potential term.

  if (complexpotential/="none") then
     scomplex_piR = scomplex_piR - alpha*complex_VPR
     scomplex_piI = scomplex_piI - alpha*complex_VPI
  end if

! Shift terms.

  if (shift/="none") then
     scomplex_piR = scomplex_piR + beta_r*DAr_complex_piR + beta_z*DAz_complex_piR
     scomplex_piI = scomplex_piI + beta_r*DAr_complex_piI + beta_z*DAz_complex_piI
  end if

! Dissipation.

  if (complexdiss/=zero) then

!    Real part.

     evolvevar => complex_phiR
     sourcevar => scomplex_phiR
     call dissipation(+1,+1,complexdiss)

     evolvevar => complex_piR
     sourcevar => scomplex_piR
     call dissipation(+1,+1,complexdiss)

     evolvevar => complex_xiR_r
     sourcevar => scomplex_xiR_r
     call dissipation(-1,+1,complexdiss)

     evolvevar => complex_xiR_z
     sourcevar => scomplex_xiR_z
     call dissipation(+1,-1,complexdiss)

!    Imaginary part.

     evolvevar => complex_phiI
     sourcevar => scomplex_phiI
     call dissipation(+1,+1,complexdiss)

     evolvevar => complex_piI
     sourcevar => scomplex_piI
     call dissipation(+1,+1,complexdiss)

     evolvevar => complex_xiI_r
     sourcevar => scomplex_xiI_r
     call dissipation(-1,+1,complexdiss)

     evolvevar => complex_xiI_z
     sourcevar => scomplex_xiI_z
     call dissipation(+1,-1,complexdiss)

  end if


! ********************************
! ***   RADIATIVE BOUNDARIES   ***
! ********************************

! Radiative boundaries are only applied to Pi.

  vl = one
  var0 = zero

  evolvevar => complex_piR
  sourcevar => scomplex_piR
  Dr_var => Dr_complex_piR
  Dz_var => Dz_complex_piR
  call radbound(vl,var0)

  evolvevar => complex_piI
  sourcevar => scomplex_piI
  Dr_var => Dr_complex_piI
  Dz_var => Dz_complex_piI
  call radbound(vl,var0)


! ***************
! ***   END   ***
! ***************

  end subroutine sources_complex
