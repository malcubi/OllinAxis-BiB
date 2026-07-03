
  subroutine idata_BosonstarCF

! *************************************************************
! ***   BOSON STAR INITIAL DATA IN CONFORMALLY FLAT GAUGE   ***
! *************************************************************

! Boson stars are solutions such that the spacetime is static
! and the complex scalar field has a harmonic dependence on time.
! This subroutine calculates initial data for a boson star
! in the conformally flat gauge using a shooting method.
!
! To obtain the initial data we assume that spacetime is
! static (K_ij=0), and also that the complex scalar field
! has the form:
!
! Phi(t,r) = phi(r) exp(i omega t)
!
! Notice that with this ansatz the stress-energy tensor
! is time-independent so the metric can be static.
!
! We also assume that the spatial metric is conformally flat.
!
! We then need to solve the Hamiltonian constraint for the
! conformal factor "psi", the Klein-Gordon equation for
! the scalar field "phi", and the maximal slicing equation
! for the lapse "alpha".
!
! The Hamiltonian constraint takes the form:
!
! __2                    5
! \/    psi  +  2 pi psi rho  =  0
!   flat
!
! where we use the flat Laplacian, and where rho is the energy
! density of the complex scalar field given by:
!
!               /          2         2       4                         2       \
! rho  =  (1/2) | [ (d phi) + (d phi) ] / psi  +  ( omega phi / alpha )  + 2 V |
!               \     r         z                                              /
!
! and V(phi) is the self-interaction potential.
!
! On the other hand, the Klein-Gordon equation becomes:
!
! Finally, the maximal slicing equation takes the form:
!
! __2
! \/ alpha  -  4 pi (rho + trS)  =  0
!
! where now:
!                                        2
! rho + trS  =  2 [ ( omega phi / alpha )  -  V ]

! Include modules.

  use param
  use arrays
  use procinfo
  use derivatives

! Extra variables.

  implicit none


! ***************
! ***   END   ***
! ***************

  end subroutine idata_BosonstarCF
