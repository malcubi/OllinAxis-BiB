
!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/sources_shift.f90,v 1.26 2021/03/01 21:06:59 malcubi Exp $

  subroutine sources_shift

! *****************************
! ***   SOURCES FOR SHIFT   ***
! *****************************

! This routine calculates the sources for the shift.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  real(8) sigma
  real(8) zero,two,third


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  two   = 2.d0
  third = 1.d0/3.d0

  if (bssnflavor=='eulerian') then
     sigma = 0.d0
  else if (bssnflavor=='lagrangian') then
     sigma = 1.d0
  end if


! ************************
! ***   STATIC SHIFT   ***
! ************************

  if ((shift=="zero").or.(shift=="static").or.(t(0,0)<shiftafter)) then

     sbeta_r = zero
     sbeta_z = zero

     if (angmom) then
        sbeta_p = zero
     end if


! ************************************************
! ***   GAMMADRIVER0: PARABOLIC GAMMA DRIVER   ***
! ************************************************

! This is the parabolic Gamma driver:
!
! d beta^i  =  csi d Delta^i
!  t                t

  else if (shift=="Gammadriver0") then

!    Shift source.

     sbeta_r = drivercsi*sDelta_r
     sbeta_z = drivercsi*sDelta_z

     if (angmom) then
        sbeta_p = drivercsi*sDelta_p
     end if

!    We are not evolving dtbeta.

     sdtbeta_r = zero
     sdtbeta_z = zero

     if (angmom) then
        sdtbeta_p = zero
     end if


! **************************************************************
! ***   GAMMADRIVER1:  FIRST ORDER HYPERBOLIC GAMMA DRIVER   ***
! **************************************************************

! This is the first order hyperbolic Gamma driver:
!
!       i        m      i             i
! d beta  =  beta d beta  +  csi Delta
!  t               m
!
! Where "csi" is the parameter that fixes the wave speed (csi=0.75
! corresponds to a speed of 1 for longitudinal modes in the
! asymptotic region).
!
! I have introduced one modification that helps avoid coordinate
! shocks:
!
!       i        m      i                     i
! d beta  =  beta d beta  +  csi fdriver Delta
!  t               m
!
! where "fdriver" is a dynamical function that evolves through:
!
!                  m
! d fdriver =  beta d fdriver  -  (2/3) fdriver DIV.beta
!  t                 m

  else if ((shift=="Gammadriver1").or.(shift=="Gammadrivershock1")) then

!    Source for shift.

     if (shift=="Gammadriver1") then

        sbeta_r = drivercsi*Delta_r + beta_r*DAr_beta_r + beta_z*DAz_beta_r
        sbeta_z = drivercsi*Delta_z + beta_r*DAr_beta_z + beta_z*DAz_beta_z

        if (driverD0) then
           sbeta_r = sbeta_r - drivercsi*Delta0_r
           sbeta_z = sbeta_z - drivercsi*Delta0_z
        end if

        if (angmom) then
           sbeta_p = drivercsi*Delta_p + beta_r*DAr_beta_p + beta_z*DAz_beta_p
           if (driverD0) then
              sbeta_p = sbeta_p - drivercsi*Delta0_p
           end if
        end if

     else if (shift=="Gammadrivershock1") then

        sfdriver = - two*third*fdriver*DIV_beta + beta_r*DAr_fdriver + beta_z*DAz_fdriver

        sbeta_r = drivercsi*fdriver*Delta_r + beta_r*Dr_beta_r + beta_z*Dz_beta_r
        sbeta_z = drivercsi*fdriver*Delta_z + beta_r*Dr_beta_z + beta_z*Dz_beta_z

        if (driverD0) then
           sbeta_r =sbeta_r - drivercsi*fdriver*Delta0_r
           sbeta_z =sbeta_z - drivercsi*fdriver*Delta0_z
        end if

        if (angmom) then
           sbeta_p = drivercsi*fdriver*Delta_p + beta_r*Dr_beta_p + beta_z*Dz_beta_p
           if (driverD0) then
              sbeta_p =sbeta_p - drivercsi*fdriver*Delta0_p
           end if
        end if

     end if

!    Damping terms.

     sbeta_r = sbeta_r - drivereta*beta_r
     sbeta_z = sbeta_z - drivereta*beta_z

     if (angmom) then
        sbeta_p = sbeta_p - drivereta*beta_p
     end if

!    We are not evolving dtbeta.

     sdtbeta_r = zero
     sdtbeta_z = zero

     if (angmom) then
        sdtbeta_p = zero
     end if


! **************************************************
! ***   GAMMADRIVER2:  HYPERBOLIC GAMMA DRIVER   ***
! **************************************************

! This is the standard second order hyperbolic Gamma driver used by
! everybody written in first order form by defining an extra variable
! "dtbeta" as the time derivative of beta:
!
!       i          m      i          i
! d beta    =  beta d beta  +  dtbeta
!  t                 m
!
!         i               i              i
! d dtbeta  =  csi d Delta  -  eta d beta
!  t                t               t
!
! As before, "csi" is the parameter that fixes the wave speed
! (csi=0.75 corresponds to a speed of 1 for longitudinal modes
! in the asymptotic region), and "eta" is a damping parameter.
! Notice that eta has dimensions of 1/distance, so that its value
! has to be changed depending on the scale of the problem.
!
! I have introduced one modification that helps avoid coordinate
! shocks:
!
!         i                        i                i
! d dtbeta  =  csi d (fdriver Delta )  -  eta d beta
!  t                t                          t
!
! where "fdriver" is a dynamical function that evolves through:
!
!                  m
! d fdriver =  beta d fdriver  -  (2/3) fdriver DIV.beta
!  t                 m

  else if ((shift=="Gammadriver2").or.(shift=="Gammadrivershock2")) then

!    Shift source.
!    The advection terms make the expression not tensorial,
!    but it seems to be much better than not having them.

     sbeta_r = dtbeta_r + beta_r*DAr_beta_r + beta_z*DAz_beta_r
     sbeta_z = dtbeta_z + beta_r*DAr_beta_z + beta_z*DAz_beta_z

     if (angmom) then
        sbeta_p = dtbeta_p + beta_r*DAr_beta_p + beta_z*DAz_beta_p
     end if

!    Source for dtbeta for standard Gamma driver.

     if (shift=="Gammadriver2") then

        sdtbeta_r = drivercsi*sDelta_r + beta_r*DAr_dtbeta_r + beta_z*DAz_dtbeta_r
        sdtbeta_z = drivercsi*sDelta_z + beta_r*DAr_dtbeta_z + beta_z*DAz_dtbeta_z

        if (angmom) then
           sdtbeta_p = drivercsi*sDelta_p + beta_r*DAr_dtbeta_p + beta_z*DAz_dtbeta_p
        end if

!    Source for dtbeta for shock avoiding Gamma driver.

     else if (shift=="Gammadrivershock2") then

        sfdriver = - two*third*fdriver*DIV_beta + beta_r*Dr_fdriver + beta_z*Dz_fdriver

        sdtbeta_r = drivercsi*(fdriver*sDelta_r + sfdriver*Delta_r) + beta_r*DAr_dtbeta_r + beta_z*DAz_dtbeta_r
        sdtbeta_z = drivercsi*(fdriver*sDelta_z + sfdriver*Delta_z) + beta_r*DAr_dtbeta_z + beta_z*DAz_dtbeta_z

        if (driverD0) then
           sdtbeta_r = sdtbeta_r - drivercsi*sfdriver*Delta0_r
           sdtbeta_z = sdtbeta_z - drivercsi*sfdriver*Delta0_z
        end if

        if (angmom) then
           sdtbeta_p = drivercsi*(fdriver*sDelta_p + sfdriver*Delta_p) + beta_r*DAr_dtbeta_p + beta_z*DAz_dtbeta_p
           if (driverD0) then
              sdtbeta_p = sdtbeta_p - drivercsi*sfdriver*Delta0_p
           end if
        end if

     end if

!    Damping term.

     sdtbeta_r = sdtbeta_r - drivereta*dtbeta_r
     sdtbeta_z = sdtbeta_z - drivereta*dtbeta_z

     if (angmom) then
        sdtbeta_p = sdtbeta_p - drivereta*dtbeta_p
     end if

!    Dissipation.

     if (geodiss/=0.d0) then

        evolvevar => dtbeta_r
        sourcevar => sdtbeta_r
        call dissipation(-1,+1,geodiss)

        evolvevar => dtbeta_z
        sourcevar => sdtbeta_z
        call dissipation(+1,-1,geodiss)

        if (angmom) then
           evolvevar => dtbeta_p
           sourcevar => sdtbeta_p
           call dissipation(+1,+1,geodiss)
        end if

     end if


! ************************
! ***   SANITY CHECK   ***
! ************************

  else

     if (rank==0) then
        print *
        print *, "Unkown shift condition."
        print *, 'Aborting! (subroutine sources_shift.f90)'
        print *
     end if

     call die

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_shift

