
  subroutine syncall

! *****************************************
! ***   SYNCHRONIZE ACROSS PROCESSORS   ***
! *****************************************

! This routine should only be called for parallel runs.

  use param

! Extra variables.

  implicit none

! Synchronize geometric variables.

  call syncgeo

! Synchronize matter variables.
  
  if (mattertype/="vacuum") then
     call syncmatt
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine syncall



