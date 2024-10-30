!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/radbound.f90,v 1.4 2020/01/29 20:18:39 malcubi Exp $

  subroutine radbound(v,var0)

! ****************************************
! ***   RADIATIVE BOUNDARY CONDITION   ***
! ****************************************

! Radiative boundaries are applied directly to the sources.
!
! The radiative boundary conditions assume that at the
! boundaries the dynamical functions behave as:
!
! f = f0 + u(R-vt) / R
!
! where v is the wave speed, f0 is the asymptotic value of f,
! and R is the distance to the origin, which in our case
! is called "rr".
!
! On a Cartesian grid this can be shown to imply:
!
! d f  =  - v [  rho d   f  +  z d f + (f - f0) ] / R
!  t                  rho         z

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer i,j

  real(8) v,var0


! *************************************
! ***   RADIATIVE BOUNDARIES ON R   ***
! *************************************

! Only processors who own the boundary apply
! the boundary condition.

  if (mod(rank+1,nprocr)==0) then
     i = Nr
     sourcevar(i,:) = - v*(r(i,:)*Dr_var(i,:) + z(i,:)*Dz_var(i,:) &
                    + (evolvevar(i,:)-var0))/rr(i,:)
  end if


! *************************************
! ***   RADIATIVE BOUNDARIES ON Z   ***
! *************************************

! Only processors who own the boundary apply
! the boundary condition.

  if (rank>=size-nprocr) then
     j = Nz
     sourcevar(:,j) = - v*(r(:,j)*Dr_var(:,j) + z(:,j)*Dz_var(:,j) &
                    + (evolvevar(:,j)-var0))/rr(:,j)
  end if

  if ((.not.eqsym).and.(rank<nprocr)) then
     j = 1-ghost
     sourcevar(:,j) = - v*(r(:,j)*Dr_var(:,j) + z(:,j)*Dz_var(:,j) &
                    + (evolvevar(:,j)-var0))/rr(:,j)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine radbound







! NOTE: The following two routines correspond to the independent
! application of the previous boundary condition in each direction.

  subroutine radbound_r(v,var0)

! ****************************************
! ***   RADIATIVE BOUNDARY CONDITION   ***
! ***          rho direction           ***
! ****************************************

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer i

  real(8) v,var0


! *************************************
! ***   RADIATIVE BOUNDARIES ON R   ***
! *************************************

! Only processors who own the boundary apply
! the boundary condition.

  if (mod(rank+1,nprocr)==0) then
     i = Nr
     sourcevar(i,:) = - v*(r(i,:)*Dr_var(i,:) + z(i,:)*Dz_var(i,:) &
                    + (evolvevar(i,:)-var0))/rr(i,:)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine radbound_r








  subroutine radbound_z(v,var0)

! ****************************************
! ***   RADIATIVE BOUNDARY CONDITION   ***
! ***           z direction            ***
! ****************************************

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer j

  real(8) v,var0


! *************************************
! ***   RADIATIVE BOUNDARIES ON Z   ***
! *************************************

! Only processors who own the boundary apply
! the boundary condition.

  if (rank>=size-nprocr) then
     j = Nz
     sourcevar(:,j) = - v*(r(:,j)*Dr_var(:,j) + z(:,j)*Dz_var(:,j) &
                    + (evolvevar(:,j)-var0))/rr(:,j)
  end if

  if ((.not.eqsym).and.(rank<nprocr)) then
     j = 1-ghost
     sourcevar(:,j) = - v*(r(:,j)*Dr_var(:,j) + z(:,j)*Dz_var(:,j) &
                    + (evolvevar(:,j)-var0))/rr(:,j)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine radbound_z

