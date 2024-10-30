!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/analysis.f90,v 1.4 2019/11/15 18:08:48 malcubi Exp $

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


! ***************
! ***   END   ***
! ***************

  end subroutine analysis

