
  subroutine auxiliary

! *******************************************************
! ***   AUXILARY QUANTITIES FOR GEOMETRY AND MATTER   ***
! *******************************************************

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  logical contains


! ********************
! ***   GEOMETRY   ***
! ********************

! Update auxiliary geometric variables.

  call auxiliary_geometry


! ******************
! ***   MATTER   ***
! ******************

! Update auxiliary matter variables.

  if (mattertype/="vacuum") then
     call auxiliary_matter
  end if


! ********************************
! ***   STRESS-ENERGY TENSOR   ***
! ********************************

  if (mattertype/="vacuum") then

!    Potentials for the different types of scalar fields (before stress-energy).

     if (contains(mattertype,"scalar").or.contains(mattertype,"complex")) then
        call potential
     end if

!    Calculate stress-energy variables.

     call stressenergy

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine auxiliary
