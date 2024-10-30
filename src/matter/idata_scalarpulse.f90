!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/matter/idata_scalarpulse.f90,v 1.22 2021/03/19 17:58:05 malcubi Exp $

  subroutine idata_scalarpulse

! *************************************
! ***   SCALAR PULSE INITIAL DATA   ***
! *************************************

! This initial data leaves the lapse, shift, conformal
! spatial metric and extrinsic curvature as in Minkowski
! (which have already been set up in "initial".)
!
! We then give a simple gaussian pulse in the scalar field
! and solve the Hamiltonian constraint for the conformal factor.
!
! Notice that the initial pulse is assumed to be time-symmetric,
! (i.e. pi=0), as otherwise we would have a non-zero momentum
! density and we would need to solve the coupled momentum and
! hamiltonian constraints.

! Include modules.

  use param
  use arrays
  use procinfo
  use derivatives

! Extra variables.

  implicit none

  integer box,level       ! Box number and level counters.

  real(8) zero,one,two,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  one  = 1.d0
  two  = 2.d0

  smallpi = acos(-one)


! *************************************************
! ***   SET UP GAUSSIAN PULSE IN SCALAR FIELD   ***
! *************************************************

! Loop over boxes and levels.

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Shell type perturbation. 

        if (scalarpulse_type=='shell') then

           if ((scalar_r0==zero).and.(scalar_z0==zero)) then
              scalar_phi = scalar_a0*exp(-(r**2/scalar_sr0**2 + z**2/scalar_sz0**2))
           else
              scalar_phi = scalar_a0 &
                         *(exp(-(rr-(r**2*scalar_r0 + z**2*scalar_z0)/rr**2)**2 &
                         *rr**2/(r**2*scalar_sr0**2 + z**2*scalar_sz0**2)) &
                         + exp(-(rr+(r**2*scalar_r0 + z**2*scalar_z0)/rr**2)**2 &
                         *rr**2/(r**2*scalar_sr0**2 + z**2*scalar_sz0**2)))
           end if

!       Torus type perturbation.

        else if (scalarpulse_type=='torus') then

           if ((scalar_r0==zero).and.(scalar_z0==zero)) then
              scalar_phi = scalar_a0*exp(-(r/scalar_sr0)**2)*exp(-(z/scalar_sz0)**2)
           else
              if (eqsym) then
                 scalar_phi = scalar_a0 &
                            *(exp(-((r-scalar_r0)/scalar_sr0)**2) + exp(-((r+scalar_r0)/scalar_sr0)**2)) &
                            *(exp(-((z-scalar_z0)/scalar_sz0)**2) + exp(-((z+scalar_z0)/scalar_sz0)**2))
              else
                 scalar_phi = scalar_a0*exp(-((z-scalar_z0)/scalar_sz0)**2) &
                            *(exp(-((r-scalar_r0)/scalar_sr0)**2) + exp(-((r+scalar_r0)/scalar_sr0)**2))
              end if
           end if

!       Baumgarte axisymmetric scalar field.

        else if (scalarpulse_type=='Baumgarte') then
           scalar_phi = scalar_a0*exp(-(r**2+(one-scalar_sz0)*z**2))
        end if

!       Spatial derivatives.

        diffvar => scalar_phi
        scalar_xi_r = diff1r(+1)
        scalar_xi_z = diff1z(+1)

!       Time derivative.

        scalar_pi = zero

!       Update potential.

        call potential

     end do
  end do


! *****************************
! ***   SOLVE HAMILTONIAN   ***
! *****************************

! Now we need to solve the Hamiltonian constraint for the
! conformal factor, which in this case takes the form:
!
! __2               5
! \/ psi +  2 pi psi rho  =  0
!
! with rho the energy density of the scalar field given by:
!
!             2      2          4
! rho  =  ( xi  +  xi  ) / 2 psi  +  V
!             r      z
!
! (assuming that the initial pulse is time-symmetric).
! Note that if V=0 the above equation becomes linear.


! ******************************************
! ***   ELLIPTIC EQUATION COEFFICIENTS   ***
! ******************************************

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       For a scalar pulse the source term is S0=0.

        ell_S0 = zero

!       Coefficient of linear term for elliptic equation.

        ell_S1 = smallpi*(scalar_xi_r**2 + scalar_xi_z**2)

!       Coefficient of psi^5.

        if (scalarpotential=="none") then
           ell_S5 = zero
        else
           ell_S5 = two*smallpi*scalar_V
        end if

     end do
  end do


! *******************************
! ***   WAVE ELLIPTIC SOLVER  ***
! *******************************

  if (ELL_solver=="wave") then

     if (rank==0) then
        print *, 'Calculating scalar pulse initial data with WaveElliptic solver'
        print *
     end if

     call wave_elliptic('flat','one')


! *******************************
! ***   SOR ELLIPTIC SOLVER   ***
! *******************************

  else if (ELL_solver=="sor") then

     if (rank==0) then
        print *, 'Calculating scalar pulse initial data with SOR elliptic solver'
        print *
     end if

     call sor_elliptic('flat','one')


! *************************
! ***   UNKNOWN SOLVER  ***
! *************************

  else

     if (rank==0) then
        print *
        print *, 'Unknow elliptic solver'
        print *, 'Aborting! (Subroutine idata_scalarpulse.f90)'
        print *
     end if

     call die

  end if


! **********************************
! ***   COPY SOLUTION ONTO PSI   ***
! **********************************

 do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Find psi,psi2,psi4.

        psi  = ell_u
        psi2 = psi**2
        psi4 = psi**4

!       Find phi,chi.

        phi = log(psi)
        chi = one/psi**dble(chipower)

     end do
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine idata_scalarpulse

