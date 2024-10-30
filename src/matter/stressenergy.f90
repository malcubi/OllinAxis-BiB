!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/matter/stressenergy.f90,v 1.11 2021/08/20 17:39:27 malcubi Exp $

  subroutine stressenergy

! *********************************************
! ***   EVALUATION OF STRESS-ENERGY TERMS   ***
! *********************************************

! This subroutine evaluates the stress-energy variables
! that appear in the Einstein equations.
!
! These quantities are defined in general as:
! 
!          mu  nu
! rho  =  n   n   T             (energy density)
!                  mu nu
!
!            mu  nu
! J    =  - n   P   T           (momentum density)
!  i             i   mu nu
!
!          mu  nu
! S    =  P   P    T            (stress tensor)
!  ij      i   j    mu nu
!
! with T_{mu,nu} the stress-energy tensor, n^mu the normal unit
! vector to the spatial hypersurfaces and P^mu_nu the projector
! operator onto the hypersurfaces.
!
! IMPORTANT:  For all types of matter, we always ADD to the values
! of the stress-energy variables.  This is because we might want
! to have more than one type of matter present in a given simulation.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  logical contains

  real(8) zero,half,one,two


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  half = 0.5d0
  one  = 1.d0
  two  = 2.d0


! *************************************************
! ***   MAKE SURE CONFORMAL FACTOR IS UPDATED   ***
! *************************************************

! This routine is typically called before updating the
! auxiliary geometry variables, so psi should not be
! assumed to have the correct value.

  if (chimethod) then
     phi = - half*log(chi)
  end if

  psi  = exp(phi)
  psi4 = psi**4


! ********************************
! ***   INITIALIZE TO VACUUM   ***
! ********************************

! Energy density.

  rho = zero

! Momentum density.

  J_r = zero
  J_z = zero

  if (angmom) then
     J_p = zero
  end if
  
! Stress tensor.

  S_A = zero
  S_B = zero
  S_H = zero
  S_C = zero

  S_lambda = zero

  if (angmom) then
     S_C1 = zero
     S_C2 = zero
  end if


! *********************************
! ***   COSMOLOGICAL CONSTANT   ***
! *********************************


! *****************************
! ***   REAL SCALAR FIELD   ***
! *****************************

! The stress-energy tensor for a scalar field "phi" with a
! self-interaction potential V has the form:
!
!                                     /  beta                     \
! T       =  d  phi d  phi  -  g      | d    phi d    phi  +  2 V | / 2
!  mu nu      mu     nu         mu nu \           beta            /
!
!
! Note: A scalar field in axisymmetry does not have angular dependency.
! From the above stress-energy tensor one finds for rho:
!
!                 2         rr  2         zz   2           rz
! rho  =  1/2 ( pi  +  gamma   X  +  gamma    X  +  2 gamma  X  X  )  +  V
!                               r              z              r  z
!
! where gamma is the physical spatial metric (so don't forget
! the psi^4 factors), pi := n^mu d_mu phi (with n^mu the unit normal
! vector to the spatial hypersurfaces), and X_i := d_i phi.
!
! For the momentum density J_i we find:
!
! J  =  - pi X
!  i          i
!
! (notice that J_phi=0).  This implies:
!
!  i         ij                 ij
! J  =  gamma   J  =  - pi gamma   X
!                j                  j
!
! Finally, for the stress tensor we find:
!
!                           2          mn
! S   =  X X  +  gamma  ( pi   -  gamma   X X   -  2 V) / 2
!  ij     i j         ij                   m n

  if (contains(mattertype,"scalar")) then

!    Energy density.

     rho = rho + half*(scalar_pi**2 + (g_A*scalar_xi_r**2 + g_B*scalar_xi_z**2 &
         + two*r*g_C*scalar_xi_r*scalar_xi_z)/psi4) + scalar_V

!    Momentum density (index up).

     J_r = J_r - (g_A*scalar_xi_r + r*g_C*scalar_xi_z)*scalar_pi/psi4
     J_z = J_z - (g_B*scalar_xi_z + r*g_C*scalar_xi_r)*scalar_pi/psi4

     if (angmom) then
        J_p = J_p - (r*g_C1*scalar_xi_r + g_C2*scalar_xi_z)*scalar_pi/psi4
     end if

!    Stress tensor.

     auxarray = half*(scalar_pi**2 &
              - (g_A*scalar_xi_r**2 + g_B*scalar_xi_z**2 &
              + two*r*g_C*scalar_xi_r*scalar_xi_z)/psi4) &
              - scalar_V

     S_A = S_A + scalar_xi_r**2 + A*psi4*auxarray
     S_B = S_B + scalar_xi_z**2 + B*psi4*auxarray
     S_C = S_C + scalar_xi_r*scalar_xi_z/r + C*psi4*auxarray
     S_H = S_H + H*psi4*auxarray

     if (angmom) then
        S_C1 = S_C1 + C1*psi4*auxarray
        S_C2 = S_C2 + C2*psi4*auxarray
     end if

!    S_lambda.

     if (.not.nolambda) then
        S_lambda = S_lambda + (scalar_xi_r/r)**2 + lambda*psi4*auxarray
     end if

  end if


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

! The stress-energy tensor for a scalar field "phi" with a
! self-interaction potential V has the form:
!
!                    *                  /  beta   *                 \
! T       =  (d  phi) d  phi  -  g      | d    phi d    phi  +  2 V | / 2
!  mu nu       mu      nu         mu nu \           beta            /
!
! But notice that in ths case the potential must be a function
! of the squqared norm of the scalar field.
!
! From this one gets the same values for (rho,J,S) as in the case
! of a real scalar field but twice (for real and imaginary parts).

  if (contains(mattertype,"complex")) then

!    Energy density.

     rho = rho + half*(complex_piR**2 + complex_piI**2 &
         +(g_A*(complex_xiR_r**2 + complex_xiI_r**2) &
         + g_B*(complex_xiR_z**2 + complex_xiI_z**2) &
         + two*r*g_C*(complex_xiR_r*complex_xiR_z + complex_xiI_r*complex_xiI_z))/psi4) &
         + complex_V

!    Momentum density (index up).

     J_r = J_r - ((g_A*complex_xiR_r + r*g_C*complex_xiR_z)*complex_piR &
               +  (g_A*complex_xiI_r + r*g_C*complex_xiI_z)*complex_piI)/psi4

     J_z = J_z - ((r*g_C*complex_xiR_r + g_B*complex_xiR_z)*complex_piR &
               +  (r*g_C*complex_xiI_r + g_B*complex_xiI_z)*complex_piI)/psi4

     if (angmom) then
        J_p = J_p - ((r*g_C1*complex_xiR_r + g_C2*complex_xiR_z)*complex_piR &
                  +  (r*g_C1*complex_xiI_r + g_C2*complex_xiI_z)*complex_piI)/psi4
     end if

!    Stress tensor.

     auxarray = half*(complex_piR**2 + complex_piI**2 &
              -(g_A*complex_xiR_r**2 + g_B*complex_xiR_z**2 &
              + two*r*g_C*complex_xiR_r*complex_xiR_z &
              + g_A*complex_xiI_r**2 + g_B*complex_xiI_z**2 &
              + two*r*g_C*complex_xiI_r*complex_xiI_z)/psi4) &
              - complex_V

     S_A = S_A + complex_xiR_r**2 + complex_xiI_r**2 + A*psi4*auxarray
     S_B = S_B + complex_xiR_z**2 + complex_xiI_z**2 + B*psi4*auxarray
     S_C = S_C + complex_xiR_r*complex_xiR_z/r + complex_xiI_r*complex_xiI_z/r + C*psi4*auxarray
     S_H = S_H + H*psi4*auxarray

     if (angmom) then
        S_C1 = S_C1 + C1*psi4*auxarray
        S_C2 = S_C2 + C2*psi4*auxarray
     end if

!    S_lambda.

     if (.not.nolambda) then
        S_lambda = S_lambda + (complex_xiR_r/r)**2 + (complex_xiI_r/r)**2 &
                 + lambda*psi4*auxarray
     end if

  end if


! **************************************
! ***   NO MORE MATTER FIELDS HERE   ***
! **************************************


! **********************************
! ***   TRACE OF STRESS TENSOR   ***
! **********************************

! Calculate the trace of stress tensor.

  if (angmom) then
     trS = g_A*S_A + g_B*S_B + g_H*S_H &
         + two*r**2*(g_C*S_C + r**2*g_C1*S_C1 + g_C2*S_C2)
  else
     trS = g_A*S_A + g_B*S_B + g_H*S_H + two*r**2*g_C*S_C
  end if

! Trace is taken with the physical metric, so we divide
! by psi4 to take this into account.

  trS = trS/psi4


! ***************
! ***   END   ***
! ***************

  end subroutine stressenergy
