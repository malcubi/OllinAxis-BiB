!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/idata_kerr.f90,v 1.7 2019/10/10 19:46:11 malcubi Exp $

  subroutine idata_kerr

! ****************************************
! ***   KERR BLACK HOLE INITIAL DATA   ***
! ****************************************

! Initial data for Kerr in quasi-isotropic coordinates.
! See: S.R. Brandt and E. Seidel, Phys.Rev.D 54, 1403.
!
! With some changes in notation, the metric can be written as:
!
!   2        4     2      2     2       2
! dl   =  psi  [ dr  +  dz  +  r  H dphi  ]
!
! where:
!
!    4         2   2
! psi  =  kappa / R
!
!                      4
! H    =  Sigma / kappa
!
! Here R=sqrt(rho**2+z**2) is the quasi-isotropic coordinate
! radius.  The quantity "kappa" is defined through:
!
!      2     2      2    2
! kappa  =  R   +  a  cos (theta)
!            BL
!
! where "a" is the spin parameter (angular momentum per unit mass),
! and R_BL is the standard Brill-Lindquist radius, related to
! the quasi-isotropic R radius through:
!
! R   =  R [ 1 + (M+a)/2R ] [ 1 + (M-a)/2R ]
!  BL
!
! Notice that in the paper by Seidel and Brandt the use "rho" instead
! of kappa, but I don't do that here to avoid confusion with the
! cylindrical coordinate rho (but remember that the code uses
! "r" instead of "rho", and "rr" instead of "R".
!
! Finally, the quantity Sigma is defined as:
!
!              2     2  2            2   2
! Sigma  =  ( R  +  a  )   -  Delta a sin (theta)
!              BL
!
! with Delta given by:
!
!          2               2
! Delta = R  -  2 M R  +  a
!          BL        BL
!
! For the Kerr metric the extrinsic curvature is not zero.  It turns
! out that in spherical coordinates the only non-zero components
! are the (R Phi) and (Theta Phi) components:
!
!                         2    2     2         2   2     2       2
! K           =  a M [ 2 R  ( R  +  a ) + kappa ( R  -  a ) ] sin (theta)
!  R phi                  BL   BL                  BL
!
!                          2                   4
!                sqrt(kappa / Sigma) / (r kappa )
!
!
! K                   3                                3
!  theta phi  =  - 2 a  M R  sqrt(Delta) cos(theta) sin (theta)
!                          BL
!
!                          2                 4
!                sqrt(kappa / Sigma) / (kappa )
!
! When transforming back to cylindrical coordinates, we find that
! only KC1 and KC2 are non-zero and are given by:
!
! KC1 = K_(rho,phi) / rho**3
! KC2 = K_(z,phi)   / rho**2
!
! Finally, KC1 and KC2 ae the physical components.  But notice that
! the physical extrinsic curvature is already traceless. To obtain
! the conformal traceless extrinsic curvature we then just need
! to rescale with the conformal factor.
!
! KTC1 = KC1/psi4
! KTC2 = KC2/psi4
!
! Notice that for this metric the outer horizon is located at:
!
! R_+ = sqrt(M**2 - a**2)/2
!
! NOTE:  In the code we use the parameters BH1mass for the mass "M",
! and BH1spin for the spin parameter "a". Notice that bot "M" and "a"
! have units of distance.

! Include modules.

  use param
  use arrays
  use derivatives

! Extra variables.

  implicit none

  real(8) zero,half,one,two

  real(8) rBL(1-ghost:Nrmax,1-ghost:Nzmax)     ! Boyer-Lindquist radial coordinate.
  real(8) del(1-ghost:Nrmax,1-ghost:Nzmax)     ! Delta function.
  real(8) sig(1-ghost:Nrmax,1-ghost:Nzmax)     ! Sigma function.
  real(8) ka2(1-ghost:Nrmax,1-ghost:Nzmax)     ! kappa**2 function.
  real(8) K_rp(1-ghost:Nrmax,1-ghost:Nzmax)    ! Extrinsic curvature spherical component (R PHI).
  real(8) K_tp(1-ghost:Nrmax,1-ghost:Nzmax)    ! Extrinsic curvature spherical component (THETA PHI).


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  half = 0.5d0
  one  = 1.d0
  two  = 2.d0


! ****************************
! ***   AUXILIARY ARRAYS   ***
! ****************************

! Notice that sin(theta)=rho/R and cos(theta)=z/R.

! Boyer-Lindquist radial coordinate defined in terms
! of the quasi-isotropic radius r as:
!
! R  =  R ( 1 + (M+a)/2R ) ( 1 + (M-a)/2R )
!  BL

  rBL = rr*(one + half*(BH1mass+BH1spin)/rr)*(one + half*(BH1mass-BH1spin)/rr)

! Delta function defined as:
!
!          2               2
! Delta = R  -  2 M R  +  a
!          BL        BL

  del = rBL**2 - two*BH1mass*rBL + BH1spin**2

! Sigma defined as:
!
!              2     2  2     2         2
! Sigma  =  ( R  +  a  )  -  a Delta sin (theta)
!              BL

  sig = (rBL**2 + BH1spin**2)**2 - del*(BH1spin*r/rr)**2

! kappa**2 function defined as:
!
!      2     2     2    2
! kappa  =  r  +  a  cos (theta)
!            BL

  ka2 = rBL**2 + (BH1spin*z/rr)**2


! *****************
! ***   LAPSE   ***
! *****************

! The lapse for the stationary metric is given through:
!
!        2                    2      2                2
! 1/alpha  =  1  +  2 M R  ( R   +  a ) / (Delta kappa )
!                        BL   BL
!
! But notice that we might not want to use it.

  if (ilapse=="quiso_kerr") then
     alpha = one/sqrt((one + two*BH1mass*rBL*(rBL**2 + BH1spin**2)/del/ka2))
  end if


! *****************
! ***   SHIFT   ***
! *****************

! And the shift is given by:
!
!     phi 
! beta    =   - 2 a M R  / Sigma
!                      BL

  if (ishift=="quiso_kerr") then
     beta_p = - two*BH1spin*BH1mass*rBL/sig
  end if


! ****************************
! ***   CONFORMAL FACTOR   ***
! ****************************

! Conformal factor psi4 given by:
!
!    4         2   2
! psi  =  kappa / R

  psi4 = ka2/rr**2

! Now find (psi,phi,chi).

  psi = sqrt(sqrt(psi4))
  phi = log(psi)
  chi = one/psi**2


! ****************************
! ***   CONFORMAL METRIC   ***
! ****************************

! The only non-trivial component of the
! conformal metric is H = Sigma/kappa**4.

  A = one
  B = one
  C = zero

  H = sig/ka2**2

  C1 = zero
  C2 = zero

! Regularization variable lambda = (A-H)/rho**2.

  lambda = (one-H)/r**2


! ***************************
! ***   DELTA VARIABLES   ***
! ***************************

! Since the metric component H is non-trivial, we
! find that Delta_r and Delta_z are also non-zero:
!
!      r
! Delta  =  rho lambda / H  -  (1/2) d   H / H
!                                     rho
!      z
! Delta  =  - (1/2) d H / H
!                    z

! Derivatives of H.

  diffvar => H
  Dr_H  = diff1r(+1)
  Dz_H  = diff1z(+1)

! Deltas.

  Delta_r = (one-H)/H/r - half*Dr_H/H
  Delta_z = - half*Dz_H/H


! *******************************
! ***   EXTRINSIC CURVATURE   ***
! *******************************

! For the Kerr metric the extrinsic curvature is not zero.  It turns
! out that in spherical coordinates the only non-zero components
! are the (R Phi) and (Theta Phi) components:
!
!                         2    2     2         2   2     2       2
! K           =  a M [ 2 R  ( R  +  a ) + kappa ( R  -  a ) ] sin (theta)
!  R phi                  BL   BL                  BL
!
!                          2                   4
!                sqrt(kappa / Sigma) / (r kappa )
!
!
! K                   3                                3
!  theta phi  =  - 2 a  M R  sqrt(Delta) cos(theta) sin (theta)
!                          BL
!
!                          2                 4
!                sqrt(kappa / Sigma) / (kappa )

  K_rp = BH1mass*BH1spin*(two*rBL**2*(rBL**2+BH1spin**2) + ka2*(rBL**2-BH1spin**2)) &
       *(r/rr)**2/ka2**2/abs(rr)*sqrt(ka2/sig)

  K_tp = - two*BH1spin**3*BH1mass*rBL*sqrt(del)*(z/rr)*(r/rr)**3/ka2**2*sqrt(ka2/sig)

! Transform back to cylindrical. Remember that:
!
! K_(rho,phi) = (rho/R) K_(R,phi) + (z/R**2)   K_(theta,phi)
! K_(z,phi)   = (z/R)   K_(R,phi) - (rho/R**2) K_(theta,phi)
!
! and also that:
!
! KC1 = K_(rho,phi) / rho**3
! KC2 = K_(z,phi)   / rho**2

  KC1 = (r*K_rp + z/rr*K_tp)/rr/r**3
  KC2 = (z*K_rp - r/rr*K_tp)/rr/r**2
  
! The above is the physical extrinsic curvature,
! which is already traceless.  To obtain the conformal
! traceless extrinsic curvature we then just need
! to rescale with the conformal factor.

  KTC1 = KC1/psi4
  KTC2 = KC2/psi4


! ***************
! ***   END   ***
! ***************

  end subroutine idata_kerr
