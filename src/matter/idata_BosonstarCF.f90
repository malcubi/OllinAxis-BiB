
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
! We also assume that the spatial metric is conformally flat:
!
!   2       4    2
! dl  =  psi   dl
!                flat
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
! __2         ij                                    2         
! \/ phi  +  g   d phi d ln(alpha)  +  (omega/alpha) phi  -  V'  = 0 
!                 i     j  
!
! Finally, the maximal slicing equation takes the form:
!
! __2
! \/ alpha  -  4 pi (rho + trS)  =  0
!
!
! where now:
!                                        2
! rho + trS  =  2 [ ( omega phi / alpha )  -  V ]
!
! Notice that the equation for psi involves the flat Laplacian,
! while the eqautions for (phi,alpha) involve the physical
! Laplacian.  For a conformally flat metric the two Laplacian
! operators for a scalar function F are related as:
!
! __2       __2         4       ij    
! \/  F  =  \/   F / psi  +  2 g   d F  d ln(psi)
!            flat                   i    j
!
!              __2          ~ij                      4
!        =  [  \/   F  +  2 g  d F  d ln(psi) ] / psi
!               flat            i    j
!
! where the tilde indicates the flat metric in cylindrical
! coordinates, so that we simply have:
!
! ~ij
! g  d F d G  =  d F d G  +  d F d G
!     i   j       r   r       z   z
!
! The equations for phi and alpha then become:
!
!
! __2          ~ij                                              4               2
! \/   phi  +  g   [ d ln(alpha)  +  2 d ln(psi) ] d phi  +  psi [ (omega/alpha) phi  -  V' ] = 0
!  flat               i                 i           j
!
! __2              ~ij
! \/   alpha  +  2 g  d ln(psi) d alpha  -  4 pi (rho + trS)  =  0
!  flat                i         j
!
!
! This routine takes as input parameter the value of the scalar
! field at the origin "boson_phi0", and solves the coupled system
! of elliptic equations.

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
