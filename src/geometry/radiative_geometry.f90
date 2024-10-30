!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/radiative_geometry.f90,v 1.16 2021/03/10 19:07:35 malcubi Exp $

  subroutine radiative_geometry

! ****************************************
! ***   RADIATIVE BOUNDARY CONDITION   ***
! ****************************************

! Radiative boundaries for geometry.
!
! Notice that in this routine we will always assume
! that the BSSN parameter eta is equal to 2.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  real(8) vl,va,vs,var0
  real(8) zero,one,two,third


! **********************************
! ***   DO WE WANT TO BE HERE?   ***
! **********************************

  if (spacetime=="background") return


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  one   = 1.d0
  two   = 2.d0
  third = 1.d0/3.d0

! Asymptotic value for all quantities.

  var0 = zero


! ***********************
! ***   EIGENSPEEDS   ***
! ***********************

! Find asymptotic eigenspeeds (vl,va,vs).
! Notice that here we assume that we are using the
! lagrangian formulation.  For the eulerian formulation
! the characteristic analysis becomes much more difficult.

  vl = 1.d0                          ! Speed of light far away.
  va = dsqrt(gauge_f)                ! Lapse eigenspeed far away.
  vs = dsqrt(4.d0*third*drivercsi)   ! Shift eigenspeed far away.


! ****************************
! ***   BOUNDARY FOR TRK   ***
! ****************************

! For Bona-Masso slicings, far away trK always behaves
! as a spherical wave with the slicing gauge speed.

  if (slicing/="maximal") then
     evolvevar => trK
     sourcevar => strK
     Dr_var => Dr_trK
     Dz_var => Dz_trK
     call radbound(va,var0)
  end if


! **********************************************
! ***   BOUNDARY FOR DIAGONAL KT COMPONENTS  ***
! **********************************************

! Component KTA (r,r).

  evolvevar => KTA
  sourcevar => sKTA
  Dr_var => Dr_KTA
  Dz_var => Dz_KTA

  call radbound(vl,var0)

! Component KTB (z,z).

  evolvevar => KTB
  sourcevar => sKTB
  Dr_var => Dr_KTB
  Dz_var => Dz_KTB
  call radbound(vl,var0)

! Component KTH (phi,phi).

  evolvevar => KTH
  sourcevar => sKTH
  Dr_var => Dr_KTH
  Dz_var => Dz_KTH
  call radbound(vl,var0)

! We need to correct the boundaries since they are
! coupled with the slicing gauge mode.

  if (mod(rank+1,nprocr)==0) then

     sKTA(Nr,:) = sKTA(Nr,:) + third*(one - vl/va)*strK(Nr,:)*two
     sKTB(Nr,:) = sKTB(Nr,:) - third*(one - vl/va)*strK(Nr,:)
     sKTH(Nr,:) = sKTH(Nr,:) - third*(one - vl/va)*strK(Nr,:)

  end if

  if ((.not.eqsym).and.(rank<nprocr)) then

     sKTA(:,1-ghost) = sKTA(:,1-ghost) - third*(one - vl/va)*strK(:,1-ghost)
     sKTB(:,1-ghost) = sKTB(:,1-ghost) + third*(one - vl/va)*strK(:,1-ghost)*two
     sKTH(:,1-ghost) = sKTH(:,1-ghost) - third*(one - vl/va)*strK(:,1-ghost)

  end if

  if (rank>=size-nprocr) then

     sKTA(:,Nz) = sKTA(:,Nz) - third*(one - vl/va)*strK(:,Nz)
     sKTB(:,Nz) = sKTB(:,Nz) + third*(one - vl/va)*strK(:,Nz)*two
     sKTH(:,Nz) = sKTH(:,Nz) - third*(one - vl/va)*strK(:,Nz)

  end if


! **************************************************
! ***   BOUNDARY FOR OFF-DIAGONAL KT COMPONENTS  ***
! **************************************************

! Component KTC (r,z).

  evolvevar => KTC
  sourcevar => sKTC
  Dr_var => Dr_KTC
  Dz_var => Dz_KTC
  call radbound(vl,var0)

! Components KTC1 (r,phi) and KTC2 (z,phi).

  if (angmom) then

     evolvevar => KTC1
     sourcevar => sKTC1
     Dr_var => Dr_KTC1
     Dz_var => Dz_KTC1
     call radbound(vl,var0)

     evolvevar => KTC2
     sourcevar => sKTC2
     Dr_var => Dr_KTC2
     Dz_var => Dz_KTC2
     call radbound(vl,var0)

  end if


! ********************************
! ***   BOUNDARY FOR ALAMBDA   ***
! ********************************

  if (.not.nolambda) then
     evolvevar => Alambda
     sourcevar => sAlambda
     Dr_var => Dr_Alambda
     Dz_var => Dz_Alambda
     call radbound(vl,var0)
  end if


! *************************************
! ***   BOUNDARY FOR DELTA/DTBETA   ***
! *************************************

! Radiative boundary for dtbeta.  At the moment,
! we only apply a boundary condition for the
! Gammadriver shift.

  if ((shift(1:11)=="Gammadriver").and.(drivercsi/=0.d0)) then

!    Gammadriver0.

     if (shift=="Gammadriver0") then

        sbeta_r(0,:) = 0.d0
        sbeta_z(0,:) = 0.d0

        sbeta_r(:,0) = 0.d0
        sbeta_z(:,0) = 0.d0

!    Gammadriver1.

     else if ((shift=="Gammadriver1").or.(shift=="Gammadrivershock1")) then

        evolvevar => Delta_r
        sourcevar => sDelta_r
        Dr_var => Dr_Delta_r
        Dz_var => Dz_Delta_r
        call radbound(vs,var0)

        evolvevar => Delta_z
        sourcevar => sDelta_z
        Dr_var => Dr_Delta_z
        Dz_var => Dz_Delta_z
        call radbound(vs,var0)

        if (angmom) then
           evolvevar => Delta_p
           sourcevar => sDelta_p
           Dr_var => Dr_Delta_p
           Dz_var => Dz_Delta_p
           call radbound(vs,var0)
        end if

!    Gammadriver2.

     else if ((shift=="Gammadriver2").or.(shift=="Gammadrivershock2")) then

        evolvevar => dtbeta_r
        sourcevar => sdtbeta_r
        Dr_var => Dr_dtbeta_r
        Dz_var => Dz_dtbeta_r
        call radbound(vs,var0)

        evolvevar => dtbeta_z
        sourcevar => sdtbeta_z
        Dr_var => Dr_dtbeta_z
        Dz_var => Dz_dtbeta_z
        call radbound(vs,var0)

        if (angmom) then
           evolvevar => dtbeta_p
           sourcevar => sdtbeta_p
           Dr_var => Dr_dtbeta_p
           Dz_var => Dz_dtbeta_p
           call radbound(vs,var0)
        end if

     end if

  end if


! **********************************
! ***   BOUNDARY FOR Z4C THETA   ***
! **********************************

! For the Z4c formulation we need to apply a
! boundary condition to z4theta.
!
! In order to impose the boundary condition we assume that
! far away z4theta behaves as an outgoing spherical wave:
!
! Theta  ~  g(r-t)/r
!
! I need to check this, I am not sure it is correct yet.

  if (formulation=="z4c") then
     evolvevar => z4theta
     sourcevar => sz4theta
     Dr_var => Dr_z4theta
     Dz_var => Dz_z4theta
     call radbound(vl,var0)

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine radiative_geometry

