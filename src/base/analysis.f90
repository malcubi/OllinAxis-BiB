
  subroutine analysis(box,level)

! ************************************************
! ***   CALLS TO DIFFERENT ANALYSIS ROUTINES   ***
! ************************************************

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  integer box,level


! ***********************
! ***   CONSTRAINTS   ***
! ***********************

! Hamiltonian, momentum and other constraints.

  call constraints(box,level)


! *********************************
! ***   ANALYSIS FOR GEOMETRY   ***
! *********************************

  call analysis_geometry(box,level)


! *******************************
! ***   ANALYSIS FOR MATTER   ***
! *******************************

   if (mattertype/="vacuum") then
     call analysis_matter(box,level)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine analysis

