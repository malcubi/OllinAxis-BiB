!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/detrace_KT.f90,v 1.2 2019/09/05 00:56:19 malcubi Exp $

  subroutine detrace_KT

! ********************************
! ***   SUBTRACT TRACE OF KT   ***
! ********************************

! Subtract trace of KT.  This routine is copied from
! the old code "OllinAxis" and was originally written
! by José Manuel Torres.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  real(8) third,two


! *******************
! ***   NUMBERS   ***
! *******************

  third = 1.d0/3.d0
  two   = 2.d0


! ****************************
! ***  FIND TRACE OF KTA   ***
! ****************************

! Calculate the trace of KT which is just numerical error
! (remember it should be traceless).

  if (.not.angmom) then
     auxarray = g_A*KTA + g_B*KTB + g_H*KTH + two*r**2*g_C*KTC
  else
     auxarray = g_A*KTA + g_B*KTB + g_H*KTH &
              + two*r**2*(g_C*KTC + r**2*g_C1*KTC1 + g_C2*KTC2)
  end if


! **************************
! ***   SUBTRACT TRACE   ***
! **************************

  KTA = KTA - third*A*auxarray
  KTB = KTB - third*B*auxarray
  KTH = KTH - third*H*auxarray
  KTC = KTC - third*C*auxarray

  if (angmom) then
     KTC1 = KTC1 - third*C1*auxarray
     KTC2 = KTC2 - third*C2*auxarray
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine detrace_KT

