!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/matter/potential.f90,v 1.7 2019/11/14 20:30:17 malcubi Exp $

  subroutine potential

! *****************************************************
! ***   POTENTIAL FOR THE DIFFERENT SCALAR FIELDS   ***
! *****************************************************

! This routine calculates the potential and its derivative
! for the different types of scalar field.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  logical contains


! *****************************
! ***   REAL SCALAR FIELD   ***
! *****************************

  if (contains(mattertype,"scalar")) then

!    **************************
!    ***   ZERO POTENTIAL   ***
!    **************************

!    By default, set to zero.

     scalar_V  = 0.d0
     scalar_VP = 0.d0


!    ***************************
!    ***   PHI^2 POTENTIAL   ***
!    ***************************

!    A quadratic potential indicates a mass term
!    for the scalar field and has the form:
!
!    V  =  m^2 phi^2 / 2

     if (scalarpotential=="phi2") then

        scalar_V  = 0.5d0*scalar_mass**2*scalar_phi**2
        scalar_VP = scalar_mass**2*scalar_phi


!    ***************************
!    ***   PHI^4 POTENTIAL   ***
!    ***************************

!    A phi^4 potential indicates a self interaction of the scalar
!    field, its strenght is characterized by a dimensionless coupling
!    constant lambda and has the form:
!             
!    V  =  lambda phi^4 / 4

     else if (scalarpotential=="phi4") then

        scalar_V  = 0.5d0*scalar_mass**2*scalar_phi**2 &
                  + 0.25d0*scalar_lambda*scalar_phi**4
        scalar_VP = scalar_mass**2*scalar_phi &
                  + scalar_lambda*scalar_phi**3

     end if

  end if


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

  if (contains(mattertype,"complex")) then

!    **************************
!    ***   ZERO POTENTIAL   ***
!    **************************

!    By default, set to zero.

     complex_V = 0.d0

     complex_VPR = 0.d0
     complex_VPI = 0.d0


!    ***************************
!    ***   PHI^2 POTENTIAL   ***
!    ***************************

!    A quadratic potential indicates a mass term
!    for the complex field and has the form:
!
!    V  =  m^2 phi^2 / 2

     if (complexpotential=="phi2") then

        complex_V = 0.5d0*complex_mass**2*(complex_phiR**2 + complex_phiI**2)

        complex_VPR = complex_mass**2*complex_phiR
        complex_VPI = complex_mass**2*complex_phiI


!    ***************************
!    ***   PHI^4 POTENTIAL   ***
!    ***************************

!    A phi^4 potential indicates a self interaction of the scalar
!    field, its strenght is characterized by a dimensionless coupling
!    constant lambda and has the form:
!             
!    V  =  lambda phi^4 / 4

     else if (complexpotential=="phi4") then

        complex_V  = 0.5d0*complex_mass**2*(complex_phiR**2 + complex_phiI**2) &
                   + 0.25d0*complex_lambda*(complex_phiR**2 + complex_phiI**2)**2

        complex_VPR = complex_mass**2*complex_phiR &
                    + complex_lambda*complex_phiR*(complex_phiR**2 + complex_phiI**2)
        complex_VPI = complex_mass**2*complex_phiI &
                    + complex_lambda*complex_phiI*(complex_phiR**2 + complex_phiI**2)

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine potential

