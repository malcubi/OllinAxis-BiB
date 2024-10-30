!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/constraints.f90,v 1.12 2020/01/31 22:38:24 malcubi Exp $

  subroutine constraints(box,level)

! *************************************
! ***   EVALUATION OF CONSTRAINTS   ***
! *************************************

! This routine evaluates the hamiltonian, momentum and
! other constraints.
!
! Remember to check if the arrays have been allocated,
! since the constraints are declared as AUXILIARY and will
! only have memory allocated if one wants output for them.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer box,level
  integer i,j

  real(8) half,third,two,four,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  two   = 2.d0
  four  = 4.d0

  half  = 0.5d0
  third = 1.d0/3.d0

  smallpi = acos(-1.0)


! ***********************
! ***   HAMILTONIAN   ***
! ***********************

! The hamiltonian constraint in this case takes the form:
!
!                                   2
! ham  =  R  -  ( KT2 )  +  2/3  trK  -  16 pi rho
!
! Where R is the physical Ricci scalar.

  if (associated(grid(box,level)%ham)) then

     ham = RSCAL - KT2 + two*third*trK**2
  
     if (mattertype/="vacuum") then
        ham = ham - 16.d0*smallpi*rho
     end if

!    Boundaries. The calculation above is done all the way to the
!    boundary with one-sided differences, but this introduces
!    large artificial errors, so I fix it here.

     if (mod(rank,nprocr)==0) then
        ham(1-ghost,:) = ham(2-ghost,:)
     end if

     if (mod(rank+1,nprocr)==0) then
        ham(Nr,:) = ham(Nr-1,:)
     end if

     if (rank<nprocr) then
        ham(:,1-ghost) = ham(:,2-ghost)
     end if

     if (rank>=size-nprocr) then
        ham(:,Nz) = ham(:,Nz-1)
     end if

!    And fix the symmetries.

     if (ownaxis) then
        do i=1,ghost
           ham(1-i,:) = ham(i,:)
        end do
     end if

     if (eqsym.and.ownequator) then
        do j=1,ghost
           ham(:,1-j) = ham(:,j)
        end do
     end if

  end if


! ********************
! ***   MOMENTUM   ***
! ********************

! The momentum constraint in this case takes the form:
!
!    i        ji          i            ij                i   4
! mom  =  D KT  -  (2/3) d trK  +  6 KT  d phi  -  8 pi J psi
!          j                              j
!
! Notice that it has the indices up.

! r component.

  if (associated(grid(box,level)%mom_r)) then

     mom_r = DIV_KT_r + 6.d0*(KT_up_rr*Dr_phi + KT_up_rz*Dz_phi) &
           - two*third*(g_A*Dr_trK + r*g_C*Dz_trK)

     if (mattertype/="vacuum") then
        mom_r = mom_r - 8.d0*smallpi*J_r*psi4
     end if

  end if

! z component.

  if (associated(grid(box,level)%mom_z)) then

     mom_z = DIV_KT_z + 6.d0*(KT_up_rz*Dr_phi + KT_up_zz*Dz_phi) &
           - two*third*(r*g_C*Dr_trK + g_B*Dz_trK)

     if (mattertype/="vacuum") then
        mom_z = mom_z - 8.d0*smallpi*J_z*psi4
     end if

  end if

! Angular component.

  if (associated(grid(box,level)%mom_p).and.angmom) then

     mom_p = DIV_KT_p + 6.d0*(KT_up_rp*Dr_phi + KT_up_zp*Dz_phi) &
           - two*third*(r*g_C1*Dr_trK + g_C2*Dz_trK)

     if (mattertype/="vacuum") then
        mom_p = mom_p - 8.d0*smallpi*J_p*psi4
     end if

  end if


! ****************************
! ***   DELTA CONSTRAINT   ***
! ****************************

!       i           i        ij          ij              ^  i
! CDelta  :=   Delta - (- d g  -  (1/2) g  d ln hdet - Gamma  )
!                          j                j
!
! Here the ^ refers to the flat space value.

! r component.

  if (associated(grid(box,level)%CDelta_r)) then
     CDelta_r = Delta_r + (Dr_g_A + r*Dz_g_C + r*g_lambda &
              + half*(g_A*Dr_hdet + r*g_C*Dz_hdet)*ihdet)
  end if

! z component.

  if (associated(grid(box,level)%CDelta_z)) then
     CDelta_z = Delta_z + (Dr_g_C*r + Dz_g_B + two*g_C &
              + half*(r*g_C*Dr_hdet + g_B*Dz_hdet)*ihdet)
  end if

! p component.

  if (associated(grid(box,level)%CDelta_p).and.angmom) then
     CDelta_p = Delta_p + (Dr_g_C1*r + Dz_g_C2 + four*g_C1 &
              + half*(r*g_C1*Dr_hdet + g_C2*Dz_hdet)*ihdet)
  end if


! **************************************
! ***   REGULARIZATION CONSTRAINTS   ***
! **************************************

! Clambda:   r**2*lambda - (A - H)

  if (associated(grid(box,level)%Clambda)) then
     Clambda = r**2*lambda - (A-H)
  end if

! CAlambda:  r**2*Alambda - (KTA - KTH)

  if (associated(grid(box,level)%CAlambda)) then
     CAlambda = r**2*Alambda - (KTA-KTH)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine constraints

