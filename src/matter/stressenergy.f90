
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

  real(8) coef
  real(8) zero,half,one,two,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  half = 0.5d0
  one  = 1.d0
  two  = 2.d0

  smallpi = acos(-one)


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
!
! From the above stress-energy tensor one finds for rho:
!
!                 2        i
! rho  =  1/2 [ Pi  +  X  X  ]  +  V
!                       i
! 
!                 2         rr  2         zz   2           rz
!      =  1/2 [ Pi  +  gamma   X  +  gamma    X  +  2 gamma  X  X  ]  +  V
!                               r              z              r  z
!
! where gamma is the physical spatial metric (so don't forget
! the psi^4 factors), Pi := n^mu d_mu phi (with n^mu the unit normal
! vector to the spatial hypersurfaces), and X_i := d_i phi.
!
! For the momentum density J_i we find:
!
! J  =  - Pi X
!  i          i
!
! (notice that J_phi=0).  This implies:
!
!  i         ij                   ij
! J  =  gamma   J  =  - Pi ( gamma   X  )
!                j                    j
!
! Finally, for the stress tensor we find:
!
!                           2          mn
! S   =  X X  +  gamma  ( Pi   -  gamma   X X   -  2 V ) / 2
!  ij     i j         ij                   m n
!
!
! Notice in particular that:
!
!               2        i
! trS  =  ( 3 Pi  -  X  X  ) / 2  -  3 V
!                     i
!
! which implies:
!                     2
! rho + trS  =  2 ( Pi  -  V )
 
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
! But notice that in this case the potential must be a function
! of the squared norm of the scalar field.
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
              -(g_A*(complex_xiR_r**2 + complex_xiI_r**2) &
              + g_B*(complex_xiR_z**2 + complex_xiI_z**2) &
              + two*r*g_C*(complex_xiR_r*complex_xiR_z + complex_xiI_r*complex_xiI_z))/psi4) &
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

!    Corrections for rotating boson stars.

     if (angmom.and.(complex_l>0)) then

        auxarray = half*(g_H*complex_l**2*(complex_phiR**2 + complex_phiI**2)/r**2 &
                 + two*g_C1*complex_l*(complex_xiI_r*complex_phiR - complex_xiR_r*complex_phiI)*r &
                 + two*g_C2*complex_l*(complex_xiI_z*complex_phiR - complex_xiR_z*complex_phiI))

!       Energy density

        rho = rho + auxarray/psi4

!       Momentum density.

        J_r = J_r + r*g_C1*complex_l*(complex_piR*complex_phiI - complex_piI*complex_phiR)/psi4
        J_r = J_r +   g_C2*complex_l*(complex_piR*complex_phiI - complex_piI*complex_phiR)/psi4
        J_p = J_p +   g_H *complex_l*(complex_piR*complex_phiI - complex_piI*complex_phiR)/psi4/r**2

!       Stress tensor.

        S_A = S_A - A*auxarray
        S_B = S_B - B*auxarray
        S_C = S_C - C*auxarray
        S_H = S_H - H*auxarray + complex_l**2*(complex_phiR**2 + complex_phiI**2)/r**2

        S_C1 = S_C1 - C1*auxarray + complex_l*(complex_xiI_r*complex_phiR - complex_xiR_r*complex_phiI)/r**3
        S_C2 = S_C2 - C2*auxarray + complex_l*(complex_xiI_z*complex_phiR - complex_xiR_z*complex_phiI)/r**2

!       S_lambda.

        S_lambda = S_lambda - complex_l**2*(complex_phiR**2 + complex_phiI**2)/r**4 - lambda*auxarray

     end if

!    Boson density and current.  These are calculated here and not
!    in the analysis routine since for charged fields we need the
!    current for the source of the electric field.
!
!    For a complex scalar field there is a conserved Noether
!    current given by:
!
!                         *                       *
!    j    =  (i/2)  (  phi  d  phi  -  phi d   phi )
!     mu                     mu             mu
!
!         =  ( phiI  d  phiR  -  phiR  d  phiI )
!                     mu                mu
!
!    From this we can define the boson "density" as:
!
!                  mu
!    Bdens   =  - n   j   =  phiR*piI  -  phiI*piR
!                      mu
!
!    and the boson "flux" as:
!
!
!    Bflux   =  j
!         i      i

     complex_Bdens = complex_phiR*complex_piI - complex_phiI*complex_piR

     complex_Bflux_r = complex_phiI*complex_xiR_r - complex_phiR*complex_xiI_r   ! Index down.
     complex_Bflux_z = complex_phiI*complex_xiR_z - complex_phiR*complex_xiI_z   ! Index down.

  end if


! *************************
! ***   MAXWELL FIELD   ***
! *************************

  if (contains(mattertype,"electric")) then

     coef = 0.125d0/smallpi

     maxw_A2 = maxw_A_r * maxw_Au_r + maxw_A_z * maxw_Au_z
     maxw_E2 = maxw_E_r * maxw_Ed_r + maxw_E_z * maxw_Ed_z
     maxw_B2 = maxw_B_r * maxw_Bd_r + maxw_B_z * maxw_Bd_z

     if (angmom) then
        maxw_A2 = maxw_A2 + maxw_A_p * maxw_Au_p
        maxw_E2 = maxw_E2 + maxw_E_p * maxw_Ed_p
        maxw_B2 = maxw_B2 + maxw_B_p * maxw_Bd_p
     end if

!    Energy density.

     rho = rho + coef * (maxw_E2 + maxw_B2)

!    Momentum density (index up).
!    Notice that momentum density is zero without angular momentum.

     if (angmom) then
        J_r = J_r + 2.d0 * (1/psi2**3) * (1/hdet**0.5) * coef *             &
                           ( maxw_Ed_p*maxw_Bd_z - maxw_Ed_z*maxw_Bd_p )
        J_z = J_z + 2.d0 * (1/psi2**3) * (1/hdet**0.5) * coef *             &
                           ( maxw_Ed_r*maxw_Bd_p - maxw_Ed_p*maxw_Bd_r )
        J_p = J_p + 2.d0 * (1/psi2**3) * (1/hdet**0.5) * coef *             &
                           ( maxw_Ed_z*maxw_Bd_r - maxw_Ed_r*maxw_Bd_z )
     end if

!    Stress tensor.

     S_A = S_A + coef * ( A * psi4 * (maxw_E2 + maxw_B2)                 &
               - two * (maxw_Ed_r*maxw_Ed_r + maxw_Bd_r*maxw_Bd_r) )

     S_B = S_B + coef * ( B * psi4 * (maxw_E2 + maxw_B2)                 &
               - two * (maxw_Ed_z*maxw_Ed_z + maxw_Bd_z*maxw_Bd_z) )

     S_H = S_H + coef * ( H * psi4 * (maxw_E2 + maxw_B2) )

     S_C = S_C + coef * ( C * psi4 * (maxw_E2 + maxw_B2)                 &
               - (two/r) * (maxw_Ed_r*maxw_Ed_z + maxw_Bd_r*maxw_Bd_z) )

     if (angmom) then
        S_H  = S_H - coef * (two/r**2) * (maxw_Ed_p**2 + maxw_Bd_p**2)         ! Warning 1/r**2
        S_C1 = S_C1 + coef * (C1 * psi4 * (maxw_E2 + maxw_B2)                &
               - (two/r**3) * (maxw_Ed_r*maxw_Ed_p + maxw_Bd_r*maxw_Bd_p) )    ! Warning con 1/r**3
        S_C2 = S_C2 + coef * (C2 * psi4 * (maxw_E2 + maxw_B2)                &
               - (two/r**2) * (maxw_Ed_z*maxw_Ed_p + maxw_Bd_z*maxw_Bd_p) )    ! Warning con 1/r**2
     end if

!    S_lambda.

!    if (.not.nolambda) then
!       S_lambda = S_lambda + (scalar_xi_r/r)**2 + lambda*psi4*auxarray
!    end if

  end if


! *******************************
! ***   COMPLEX PROCA FIELD   ***
! *******************************

  if (contains(mattertype,"complexproca")) then

     coef = 0.125d0/smallpi

     proc_phi2 = proc_phiR**2 + proc_phiI**2

     proc_A2   = proc_AR_r * proc_ARu_r + proc_AI_r * proc_AIu_r   &
               + proc_AR_z * proc_ARu_z + proc_AI_z * proc_AIu_z
     proc_E2   = proc_ER_r * proc_ERd_r + proc_EI_r * proc_EId_r   &
               + proc_ER_z * proc_ERd_z + proc_EI_z * proc_EId_z
     proc_B2   = proc_BR_r * proc_BRd_r + proc_BI_r * proc_BId_r   &
               + proc_BR_z * proc_BRd_z + proc_BI_z * proc_BId_z

     if (angmom) then
        proc_A2   = proc_A2 + proc_AR_p * proc_ARu_p + proc_AI_p * proc_AIu_p
        proc_E2   = proc_E2 + proc_ER_p * proc_ERd_p + proc_EI_p * proc_EId_p
        proc_B2   = proc_B2 + proc_BR_p * proc_BRd_p + proc_BI_p * proc_BId_p
     end if

!    Energy density.

     rho = rho + coef * (proc_E2 + proc_B2 + proc_mass**2 * (proc_phi2 + proc_A2))

!    Momentum density (index up).

     J_r = J_r + 2.d0*coef*proc_mass**2 * (proc_ARu_r*proc_phiR + proc_AIu_r*proc_phiI)
     J_z = J_z + 2.d0*coef*proc_mass**2 * (proc_ARu_z*proc_phiR + proc_AIu_z*proc_phiI)

     if (angmom) then
        J_r = J_r + 2.d0 * (1/psi2**3) * (1/hdet**0.5) * coef *                 &
                           ( proc_ERd_p*proc_BRd_z + proc_EId_p*proc_BId_z      &
                           - proc_ERd_z*proc_BRd_p - proc_EId_z*proc_BId_p )
        J_z = J_z + 2.d0 * (1/psi2**3) * (1/hdet**0.5) * coef *                 &
                           ( proc_ERd_r*proc_BRd_p + proc_EId_r*proc_BId_p      &
                           - proc_ERd_p*proc_BRd_r - proc_EId_p*proc_BId_r )
        J_p = J_p + 2.d0 * (1/psi2**3) * (1/hdet**0.5) * coef *                 &
                           ( proc_ERd_z*proc_BRd_r + proc_EId_z*proc_BId_r      &
                           - proc_ERd_r*proc_BRd_z - proc_EId_r*proc_BId_z )    &
                  + 2.d0 * proc_mass**2 * coef *                                &
                           (proc_ARu_p*proc_phiR + proc_AIu_p*proc_phiI)
     end if

!    Stress tensor.

     S_A = S_A + coef * (A * psi4 * (proc_E2 + proc_B2)                      &
               - two * (proc_ERd_r*proc_ERd_r + proc_EId_r*proc_EId_r)       &
               - two * (proc_BRd_r*proc_BRd_r + proc_BId_r*proc_BId_r) )

     S_B = S_B + coef * (B * psi4 * (proc_E2 + proc_B2)                      &
               - two * (proc_ERd_z*proc_ERd_z + proc_EId_z*proc_EId_z)       &
               - two * (proc_BRd_z*proc_BRd_z + proc_BId_z*proc_BId_z) )

     S_H = S_H + coef * (H * psi4 * (proc_E2 + proc_B2) )

     S_C = S_C + coef * (C * psi4 * (proc_E2 + proc_B2)                          &
               - (two/r) * (proc_ERd_r*proc_ERd_z + proc_EId_r*proc_EId_z)       &
               - (two/r) * (proc_BRd_r*proc_BRd_z + proc_BId_r*proc_BId_z) )  ! Warning 1/r

     if (angmom) then
        S_H  = S_H - coef * (two/r**2) * (proc_ERd_p**2 + proc_EId_p**2             &  ! Warning 1/r**2
                                        + proc_BRd_p**2 + proc_BId_p**2)
        S_C1 = S_C1 + coef * (C1 * psi4 * (proc_E2 + proc_B2)                       &
               - (two/r**3) * (proc_ERd_r*proc_ERd_p + proc_EId_r*proc_EId_p)       &  ! Warning 1/r**3
               - (two/r**3) * (proc_BRd_r*proc_BRd_p + proc_BId_r*proc_BId_p) )        ! Warning 1/r**3
        S_C2 = S_C2 + coef * (C2 * psi4 * (proc_E2 + proc_B2)                       &
               - (two/r**2) * (proc_ERd_z*proc_ERd_p + proc_EId_z*proc_EId_p)       &  ! Warning 1/r**2
               - (two/r**2) * (proc_BRd_z*proc_BRd_p + proc_BId_z*proc_BId_p) )        ! Warning 1/r**2
     end if

!    We add mass terms to the stress tensor.

     S_A = S_A + coef * proc_mass**2 * ( - A*psi4 * (proc_A2 - proc_phi2)    &
                + two * ( proc_AR_r*proc_AR_r + proc_AI_r*proc_AI_r ) )

     S_B = S_B + coef * proc_mass**2 * ( - B*psi4 * (proc_A2 - proc_phi2)    &
                + two * ( proc_AR_z*proc_AR_z + proc_AI_z*proc_AI_z ) )

     S_H = S_H + coef * proc_mass**2 * ( - H*psi4 * (proc_A2 - proc_phi2) )

     S_C = S_C + coef * proc_mass**2 * ( - C*psi4 * (proc_A2 - proc_phi2)    &
                + (two/r) * ( proc_AR_r*proc_AR_z + proc_AI_r*proc_AI_z ) )   ! Warning 1/r

     if (angmom) then
        S_H  = S_H + coef * proc_mass**2 * (two/r**2)    &                    ! Warning 1/r**2
                          * (proc_AR_p*proc_AR_p + proc_AI_p*proc_AI_p)
        S_C1 = S_C1 + coef * proc_mass**2 * ( - C1*psi4 * (proc_A2 - proc_phi2)    &
                    + (two/r**3) * (proc_AR_r*proc_AR_p + proc_AI_r*proc_AI_p) )  ! Warning 1/r**3
        S_C2 = S_C2 + coef * proc_mass**2 * ( - C2*psi4 * (proc_A2 - proc_phi2)    &
                    + (two/r**2) * (proc_AR_z*proc_AR_p + proc_AI_z*proc_AI_p) )  ! Warning 1/r**2
     end if

!    S_lambda.

!    if (.not.nolambda) then
!       S_lambda = S_lambda + (scalar_xi_r/r)**2 + lambda*psi4*auxarray
!    end if

  end if


! **************************************
! ***   NO MORE MATTER FIELDS HERE   ***
! **************************************

! STOP:  DON'T ADD NEW MATTER FIELDS AFTER THIS POINT!


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

! The trace is taken with the physical metric,
! so we divide by psi4 to take this into account.

  trS = trS/psi4


! ***************
! ***   END   ***
! ***************

  end subroutine stressenergy
