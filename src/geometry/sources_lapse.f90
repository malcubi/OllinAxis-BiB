!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/sources_lapse.f90,v 1.12 2021/03/09 00:37:59 malcubi Exp $

  subroutine sources_lapse

! *****************************
! ***   SOURCES FOR LAPSE   ***
! *****************************

! This routine calculates the sources for the lapse.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  real(8) va


! *****************************
! ***   STATIC OR MAXIMAL   ***
! *****************************

! Sources is set to zero for static or maximal slicing.

  if ((slicing=="static").or.(slicing=="maximal")) then

!    Set lapse source to 0.

     salpha = 0.d0


! **********************************
! ***   FIRST ORDER BONA-MASSO   ***
! **********************************

! Bona-Masso type slicings, including harmonic:
!
!                                            i
! d alpha  = - alpha**2 f(alpha) trK  +  beta  d alpha
!  t                                            i

  else if ((slicing=="harmonic").or.(slicing=="1+log").or. &
           (slicing=="shockavoid")) then

!    First the term with no shift:
!
!    d alpha = - alpha**2 f(alpha) trK
!     t                            
!
!    But remember that falpha = alpha**2*f(alpha).

     salpha = - falpha*trK

!    Shift terms.

     if (shift/="none") then
        salpha = salpha + beta_r*DAr_alpha + beta_z*DAz_alpha
     end if

!    David Hilditch adds this term following the evolution equations of the
!    generalized harmonic formulation. PhysRevD.96.104051

     if (gauge_eta/=0.d0) then
        salpha = salpha - gauge_eta*log(abs(alpha))
     end if

!    Dissipation.

     !if (geodiss/=0.d0) then
     !   evolvevar => alpha
     !   sourcevar => salpha
     !   call dissipation(+1,+1,geodiss)
     !end if

!    We don't use dtalpha.

     dtalpha = 0.d0


! ***********************************
! ***   SECOND ORDER BONA-MASSO   ***
! ***********************************

  else if ((slicing=="harmonic2").or.(slicing=="1+log2").or. &
           (slicing=="shockavoid2")) then

!    Source for alpha.

     salpha = dtalpha

!    Source for dtalpha.

     sdtalpha = - falpha*strK - gauge_eta*dtalpha

!    Shift terms.

     if (shift/="none") then
        salpha = salpha + beta_r*DAr_alpha + beta_z*DAz_alpha
        sdtalpha = sdtalpha + beta_r*DAr_dtalpha + beta_z*DAz_dtalpha
     end if

!    Dissipation.

     !if (geodiss/=0.d0) then
     !   evolvevar => dtalpha
     !   sourcevar => sdtalpha
     !   call dissipation(+1,+1,geodiss)
     !end if

!    Boundary.

     va = dsqrt(gauge_f)

     evolvevar => dtalpha
     sourcevar => sdtalpha
     Dr_var => Dr_dtalpha
     Dz_var => Dz_dtalpha
     call radbound(va,0.0d0)


! ************************
! ***   SANITY CHECK   ***
! ************************

  else

     print *
     print *, 'Unknown slicing condition.'
     print *, 'Aborting! (subroutine sources_lapse.f90)'
     print *

     call die

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_lapse

