
  subroutine sources_maxwell

! *****************************************
! ***   SOURCES FOR THE MAXWELL FIELD   ***
! *****************************************

! This routine calculates the sources for
! the Maxwell electromagnetic field.

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

! Aquí falta escribir


! Source for scalar potential.

  smaxw_phi = alpha * maxw_phi * trK                                          &
               - maxw_A_r*Dr_alpha - maxw_A_z*Dz_alpha                        &
               - alpha*maxw_Au_r * (8.d0*Dr_phi + (1.d0/hdet)*Dr_hdet)        &
               - alpha*maxw_Au_z * (8.d0*Dz_phi + (1.d0/hdet)*Dz_hdet)        &
               - (alpha/psi4) * ( r * ( maxw_A_r*Dz_g_C + maxw_A_z*Dr_g_C)    &
                     + maxw_A_r*Dr_g_A + maxw_A_z*g_C + maxw_A_z*Dz_g_B)      &
               - (alpha/psi4) * (g_A*Dr_maxw_A_r + r*g_C*Dz_maxw_A_r +        &
                                 g_B*Dz_maxw_A_z + r*g_C*Dr_maxw_A_z )

  if (shift/="none") then 
     smaxw_phi = smaxw_phi + zero   ! Pongo zero provisionalmente
  end if

! Source for vector potential.

  smaxw_A_r = - alpha * Dr_maxw_phi - maxw_phi * Dr_alpha      &
              - alpha * psi4 * (A*maxw_E_r + r*C*maxw_E_z)

  smaxw_A_z = - alpha * Dz_maxw_phi - maxw_phi * Dz_alpha      &
              - alpha * psi4 * (r*C*maxw_E_r + B*maxw_E_z)

  if (shift/="none") then 
     smaxw_A_r = smaxw_A_r + zero   ! Pongo cero provisionalmente
     smaxw_A_z = smaxw_A_z + zero   ! Pongo cero provisionalmente
  end if

! Source for electric field.

  smaxw_E_r = alpha*trK*maxw_E_r    ! En campo de Proca aqui debe sumarse un termino proporcional a la masa
  smaxw_E_z = alpha*trK*maxw_E_z

  if (shift/="none") then 
     smaxw_E_r = smaxw_E_r + zero   ! Pongo cero provisionalmente
     smaxw_E_z = smaxw_E_z + zero   ! Pongo cero provisionalmente
  end if

! Source for magnetic field.

  smaxw_B_r = alpha*trK*maxw_B_r    ! En campo de Proca aqui debe sumarse un termino proporcional a la masa
  smaxw_B_z = alpha*trK*maxw_B_z

  if (shift/="none") then 
     smaxw_B_r = smaxw_B_r + zero   ! Pongo cero provisionalmente
     smaxw_B_z = smaxw_B_z + zero   ! Pongo cero provisionalmente
  end if

! Dissipation.

  if (maxwelldiss/=zero) then

     evolvevar => maxw_phi
     sourcevar => smaxw_phi
     call dissipation(+1,+1,maxwelldiss)

     evolvevar => maxw_A_r
     sourcevar => smaxw_A_r
     call dissipation(-1,+1,maxwelldiss)

     evolvevar => maxw_A_z
     sourcevar => smaxw_A_z
     call dissipation(+1,-1,maxwelldiss)

     evolvevar => maxw_E_r
     sourcevar => smaxw_E_r
     call dissipation(-1,+1,maxwelldiss)

     evolvevar => maxw_E_z
     sourcevar => smaxw_E_z
     call dissipation(+1,-1,maxwelldiss)

     evolvevar => maxw_B_r
     sourcevar => smaxw_B_r
     call dissipation(-1,+1,maxwelldiss)

     evolvevar => maxw_B_z
     sourcevar => smaxw_B_z
     call dissipation(+1,-1,maxwelldiss)

  end if


! ********************************
! ***   RADIATIVE BOUNDARIES   ***
! ********************************

  vl = one
  var0 = zero

  evolvevar => maxw_phi
  sourcevar => smaxw_phi
  Dr_var => Dr_maxw_phi
  Dz_var => Dz_maxw_phi
  call radbound(vl,var0)

  evolvevar => maxw_A_r
  sourcevar => smaxw_A_r
  Dr_var => Dr_maxw_A_r
  Dz_var => Dz_maxw_A_r
  call radbound(vl,var0)

  evolvevar => maxw_A_z
  sourcevar => smaxw_A_z
  Dr_var => Dr_maxw_A_z
  Dz_var => Dz_maxw_A_z
  call radbound(vl,var0)

  evolvevar => maxw_E_r
  sourcevar => smaxw_E_r
  Dr_var => Dr_maxw_E_r
  Dz_var => Dz_maxw_E_r
  call radbound(vl,var0)

  evolvevar => maxw_E_z
  sourcevar => smaxw_E_z
  Dr_var => Dr_maxw_E_z
  Dz_var => Dz_maxw_E_z
  call radbound(vl,var0)
  
  evolvevar => maxw_B_r
  sourcevar => smaxw_B_r
  Dr_var => Dr_maxw_B_r
  Dz_var => Dz_maxw_B_r
  call radbound(vl,var0)

  evolvevar => maxw_B_z
  sourcevar => smaxw_B_z
  Dr_var => Dr_maxw_B_z
  Dz_var => Dz_maxw_B_z
  call radbound(vl,var0)


! ***************
! ***   END   ***
! ***************

  end subroutine sources_maxwell
