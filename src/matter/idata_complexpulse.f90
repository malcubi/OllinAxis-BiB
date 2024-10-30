!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/matter/idata_complexpulse.f90,v 1.15 2021/03/01 23:06:03 malcubi Exp $

  subroutine idata_complexpulse

! **************************************
! ***   COMPLEX PULSE INITIAL DATA   ***
! **************************************

! This initial data leaves the lapse, shift, conformal
! spatial metric and extrinsic curvature as in Minkowski
! (which have already been set up in "initial").
!
! We then give a simple gaussian pulse in the complex scalar field
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


! *********************************************************
! ***   SET UP GAUSSIAN PULSE IN COMPLEX SCALAR FIELD   ***
! *********************************************************

! Loop over boxes and levels.

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Shell type perturbation. 

        if (complexpulse_type=='shell') then

!          Real part.

           if ((complexR_r0==zero).and.(complexR_z0==zero)) then
              complex_phiR = complexR_a0*exp(-rr**4/(r**2*complexR_sr0**2 + z**2*complexR_sz0**2))
           else
              complex_phiR = complexR_a0 &
                      *(exp(-(rr-(r**2*complexR_r0 + z**2*complexR_z0)/rr**2)**2 &
                      *rr**2/(r**2*complexR_sr0**2 + z**2*complexR_sz0**2)) &
                      + exp(-(rr+(r**2*complexR_r0 + z**2*complexR_z0)/rr**2)**2 &
                      *rr**2/(r**2*complexR_sr0**2 + z**2*complexR_sz0**2)))
           end if

!          Imaginary part.

           if ((complexI_r0==zero).and.(complexI_z0==zero)) then
              complex_phiI = complexI_a0*exp(-rr**4/(r**2*complexI_sr0**2 + z**2*complexI_sz0**2))
           else
              complex_phiI = complexI_a0 &
                      *(exp(-(rr-(r**2*complexI_r0 + z**2*complexI_z0)/rr**2)**2 &
                      *rr**2/(r**2*complexI_sr0**2 + z**2*complexI_sz0**2)) &
                      + exp(-(rr+(r**2*complexI_r0 + z**2*complexI_z0)/rr**2)**2 &
                      *rr**2/(r**2*complexI_sr0**2 + z**2*complexI_sz0**2)))
           end if

!       Torus type perturbation.

        else if (complexpulse_type=='torus') then

!          Real part.

           if ((complexR_r0==zero).and.(complexR_z0==zero)) then
              complex_phiR = complexR_a0*exp(-(r/complexR_sr0)**2)*exp(-(z/complexR_sz0)**2)
           else
              if (eqsym) then
                 complex_phiR = complexR_a0 &
                         *(exp(-((r-complexR_r0)/complexR_sr0)**2) + exp(-((r+complexR_r0)/complexR_sr0)**2)) &
                         *(exp(-((z-complexR_z0)/complexR_sz0)**2) + exp(-((z+complexR_z0)/complexR_sz0)**2))
              else
                 complex_phiR = complexR_a0*exp(-((z-complexR_z0)/complexR_sz0)**2) &
                         *(exp(-((r-complexR_r0)/complexR_sr0)**2) + exp(-((r+complexR_r0)/complexR_sr0)**2))
              end if
           end if

!          Imaginary part.

           if ((complexI_r0==zero).and.(complexI_z0==zero)) then
              complex_phiI = complexI_a0*exp(-(r/complexI_sr0)**2)*exp(-(z/complexI_sz0)**2)
           else
              if (eqsym) then
                 complex_phiI = complexI_a0 &
                         *(exp(-((r-complexI_r0)/complexI_sr0)**2) + exp(-((r+complexI_r0)/complexI_sr0)**2)) &
                         *(exp(-((z-complexI_z0)/complexI_sz0)**2) + exp(-((z+complexI_z0)/complexI_sz0)**2))
              else
                 complex_phiI = complexI_a0*exp(-((z-complexI_z0)/complexI_sz0)**2) &
                         *(exp(-((r-complexI_r0)/complexI_sr0)**2) + exp(-((r+complexI_r0)/complexI_sr0)**2))
              end if
           end if

        end if

!       Spatial derivatives.

        diffvar => complex_phiR
        complex_xiR_r = diff1r(+1)
        complex_xiR_z = diff1z(+1)

        diffvar => complex_phiI
        complex_xiI_r = diff1r(+1)
        complex_xiI_z = diff1z(+1)

!       Time derivatives.

        complex_piR = zero
        complex_piI = zero

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
! with rho the energy density of the complex scalar field given by:
!
!                2         2          4
! rho  =  ( |xi |  +  |xi |  ) / 2 psi  +  V
!              r         z
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

!       For a complex pulse the source term is S0=0.

        ell_S0 = zero

!       Coefficient of linear term for elliptic equation.

        ell_S1 = smallpi*(complex_xiR_r**2 + complex_xiR_z**2 &
                        + complex_xiI_r**2 + complex_xiI_z**2)

!       Coefficient of psi^5.

        if (complexpotential=="none") then
           ell_S5 = zero
        else
           ell_S5 = two*smallpi*complex_V
        end if

     end do
  end do


! *******************************
! ***   WAVE ELLIPTIC SOLVER  ***
! *******************************

  if (ELL_solver=="wave") then

     if (rank==0) then
        print *, 'Calculating complex pulse initial data with WaveElliptic solver'
        print *
     end if

     call wave_elliptic('flat','one')


! *******************************
! ***   SOR ELLIPTIC SOLVER   ***
! *******************************

  else if (ELL_solver=="sor") then

     if (rank==0) then
        print *, 'Calculating complex pulse initial data with SOR elliptic solver'
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
        print *, 'Aborting! (Subroutine idata_complexpulse.f90)'
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

  end subroutine idata_complexpulse
