!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/matter/scalar.f90,v 1.14 2021/03/05 23:01:56 malcubi Exp $ 

  subroutine sources_scalar

! ************************************
! ***   SOURCES FOR SCALAR FIELD   ***
! ************************************

! This routine calculates the sources for a real
! scalar field.  The Klein-Gordon equation with
! an arbitrary potential V(phi) has the form:
!
! Box (phi)  -  (dV/dphi)  =  0
!
! The above equation is evolved in first order form.
! Remember that we have defined the space and time
! derivative arrays as:
!
!                        m
! pi  =  ( d phi  -  beta d phi ) / alpha
!           t              m
!
! x   =  d phi
!  i      i
!
! The evolution equations in first order form are:
!
!                            m
! d phi  =  alpha pi  +  beta d phi
!  t                           m
!
!                                             m                m
! d x    =  alpha d pi  +  pi d alpha  +  beta d x  +  x d beta
!  t i             i           i                m i     m i
!
!                             m                    m
! d pi   =  alpha pi trK  +  d alpha d phi  +  beta d pi  -  alpha (dV/dphi)
!  t                                  m              m
!        +  alpha Lap(phi)
!
!
! with Lap(phi) the 3D Laplacian given by:
!
!                ij             k
! Lap(phi)  =   g  d d phi  -  G d phi
!                   i j           k
!
! Here g^ij is the physical inverse metric, and G^k
! are the physical contracted Christoffel symbols, which
! can be expressed in terms of the conformal metric and
! conformal factor psi as:
!
!  i    ~ i      ~ij
! G  =  G   -  2 g  d ln(psi)
!                    j
!
! Finally, for the contracted conformal Christoffel symbols
! we have:
!
! ~i        ~ im    im      ~
! G  =  - d g   -  g   d ln(g) / 2
!          m            m
!
! with g~ the determinant of the conformal metric.

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

  sscalar_phi = alpha*scalar_pi

  if (shift/="none") then
     sscalar_phi = sscalar_phi + beta_r*scalar_xi_r + beta_z*scalar_xi_z
  end if

! Source for xi_r.

  sscalar_xi_r = alpha*Dr_scalar_pi + scalar_pi*Dr_alpha

  if (shift/="none") then
     sscalar_xi_r = sscalar_xi_r &
                  + beta_r*DAr_scalar_xi_r + scalar_xi_r*Dr_beta_r &
                  + beta_z*DAz_scalar_xi_r + scalar_xi_z*Dr_beta_z 
  end if

! Source for xi_z.

  sscalar_xi_z = alpha*Dz_scalar_pi + scalar_pi*Dz_alpha

  if (shift/="none") then
     sscalar_xi_z = sscalar_xi_z &
                  + beta_r*DAr_scalar_xi_z + scalar_xi_r*Dz_beta_r &
                  + beta_z*DAz_scalar_xi_z + scalar_xi_z*Dz_beta_z 
  end if

! Source for time derivative. There are many terms,
! so we build the source step by step.
!
! Notice that there are two ways of doing this: using
! only derivatives of phi, or using the xi's.

  if (scalarmethod == "first") then

!    Terms with derivatives of xi.

     sscalar_pi = alpha*(g_A*Dr_scalar_xi_r + g_B*Dz_scalar_xi_z &
                + r*g_C*(Dr_scalar_xi_z + Dz_scalar_xi_r))

!    Terms with derivatives of conformal factor.

     sscalar_pi = sscalar_pi + two*alpha &
                *(Dr_phi*(g_A*scalar_xi_r + r*g_C*scalar_xi_z) &
                + Dz_phi*(g_B*scalar_xi_z + r*g_C*scalar_xi_r))

!    Terms with derivatives of alpha.

     sscalar_pi = sscalar_pi &
                + Dr_alpha*(g_A*scalar_xi_r + r*g_C*scalar_xi_z) &
                + Dz_alpha*(g_B*scalar_xi_z + r*g_C*scalar_xi_r)

!    Terms coming from derivatives of inverse conformal metric.

     sscalar_pi = sscalar_pi + alpha &
                *(scalar_xi_r*(Dr_g_A + r*Dz_g_C) + scalar_xi_z*(Dz_g_B + r*Dr_g_C + g_C))

!    Terms coming from derivatives of the determinant of the
!    conformal metric g. Remember that det(g) = r**2*h.

     sscalar_pi = sscalar_pi + alpha*(g_A*scalar_xi_r/r + g_C*scalar_xi_z)

     sscalar_pi = sscalar_pi + half*alpha*ihdet &
                *(Dr_hdet*(g_A*scalar_xi_r + r*g_C*scalar_xi_z) &
                + Dz_hdet*(g_B*scalar_xi_z + r*g_C*scalar_xi_r))

  else

!    Terms with second derivatives of phi.

     sscalar_pi = alpha*(g_A*Drr_scalar_phi + g_B*Dzz_scalar_phi &
                + two*r*g_C*Drz_scalar_phi)

!    Terms with derivatives of conformal factor.

     sscalar_pi = sscalar_pi + two*alpha &
                *(Dr_phi*(g_A*Dr_scalar_phi + r*g_C*Dz_scalar_phi) &
                + Dz_phi*(g_B*Dz_scalar_phi + r*g_C*Dr_scalar_phi))

!    Terms with derivatives of alpha.

     sscalar_pi = sscalar_pi &
                + Dr_alpha*(g_A*Dr_scalar_phi + r*g_C*Dz_scalar_phi) &
                + Dz_alpha*(g_B*Dz_scalar_phi + r*g_C*Dr_scalar_phi)

!    Terms coming from derivatives of inverse conformal metric.

     sscalar_pi = sscalar_pi + alpha &
                *(Dr_scalar_phi*(Dr_g_A + r*Dz_g_C) + Dz_scalar_phi*(Dz_g_B + r*Dr_g_C + g_C))

!    Terms coming from derivatives of the determinant of the
!    conformal metric g. Remember that det(g) = r**2*h.

     sscalar_pi = sscalar_pi + alpha*(g_A*Dr_scalar_phi/r + g_C*Dz_scalar_phi)

     sscalar_pi = sscalar_pi + half*alpha*ihdet &
                *(Dr_hdet*(g_A*Dr_scalar_phi + r*g_C*Dz_scalar_phi) &
                + Dz_hdet*(g_B*Dz_scalar_phi + r*g_C*Dr_scalar_phi))

  end if

! Rescale by conformal factor since all the above terms
! need to be divided by psi**4.

  sscalar_pi = sscalar_pi/psi4

! Term proportional to trK.

  sscalar_pi = sscalar_pi + alpha*scalar_pi*trK

! Potential term.

  if (scalarpotential/="none") then
     sscalar_pi = sscalar_pi - alpha*scalar_VP
  end if

! Shift terms.

  if (shift/="none") then
     sscalar_pi = sscalar_pi + beta_r*DAr_scalar_pi + beta_z*DAz_scalar_pi
  end if

! Dissipation.

  if (scalardiss/=zero) then

     evolvevar => scalar_phi
     sourcevar => sscalar_phi
     call dissipation(+1,+1,scalardiss)

     evolvevar => scalar_pi
     sourcevar => sscalar_pi
     call dissipation(+1,+1,scalardiss)

     evolvevar => scalar_xi_r
     sourcevar => sscalar_xi_r
     call dissipation(-1,+1,scalardiss)

     evolvevar => scalar_xi_z
     sourcevar => sscalar_xi_z
     call dissipation(+1,-1,scalardiss)

  end if


! ********************************
! ***   RADIATIVE BOUNDARIES   ***
! ********************************

! Radiative boundaries are only applied to Pi.

  vl = one
  var0 = zero

  evolvevar => scalar_pi
  sourcevar => sscalar_pi
  Dr_var => Dr_scalar_pi
  Dz_var => Dz_scalar_pi
  call radbound(vl,var0)


! ***************
! ***   END   ***
! ***************

  end subroutine sources_scalar

