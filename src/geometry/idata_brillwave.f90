!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/idata_brillwave.f90,v 1.13 2022/10/24 22:00:34 malcubi Exp $

  subroutine idata_brillwave

! ***********************************
! ***   BRILL WAVE INITIAL DATA   ***
! ***********************************

! This routine sets up initial data for Brill wave spacetimes.
! In cylindrical coordinates the spatial metric has the form:
!                                                            
!       2       4               2    2     2                
!     dl  =  psi  [ exp(2q) ( dr + dz ) + r dphi ]          
!                                                            
! To obtain the initial data one solves the elliptic equation:
!                                                            
!       __2                                                  
!     ( \/  +  s ) psi  =  0                                
!                                                            
! where the flat laplacian is given by:                                                     
!                                                            
!      __2     2    2      2    2                             
!      \/  =  d / dr   +  d / dz  +  (1/r) d / dr         
!                                                            
! and
!                    2    2     2    2                          
!       s  =  1/4 ( d q/dr  +  d q/dz )                         
!
! Here the function q is given by:
!
!               2              2        2   
!       q = a0 r  exp[ - (r/sr) - (z/sz)  ]  (Holz-type) data.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer box,level           ! Box number and level counters.

  real(8) zero,half,one,two   ! Numbers.


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.0d0
  half = 0.5d0
  one  = 1.0d0
  two  = 2.0d0


! *********************************
! ***   BRILL SOURCE FUNCTION   ***
! *********************************

! Sanity check.

  if (ELL_solver=="none") then

     if (rank==0) then
        print *
        print *, 'For Brill wave data you need an elliptic solver!'
        print *, 'Aborting! (Subroutine idata_brillwave.f90)'
        print *
     end if

     call die

  end if

! Different Brill source functions.

  if (brilltype=="Holz") then

!    Loop over boxes and grids levels.

     do box=0,Nb
        do level=min(1,box),Nl(box)

!          Point to current grid.

           call currentgrid(box,level,grid(box,level))

!          Brill q function.

           brillq = brill_a0*r**2*exp(-(r/brill_sr0)**2-(z/brill_sz0)**2)

!          Coefficient of linear term for elliptic equation.

           ell_S1 = half*brillq*(one/r**2 &
                  + two*((r/brill_sr0**2)**2 + (z/brill_sz0**2)**2) &
                  - 5.d0/brill_sr0**2 - one/brill_sz0**2)

!          For Brill waves S0=S5=0.

           ell_S0 = zero
           ell_S5 = zero

        end do
     end do

! Unknown Brill source function.

  else

     if (rank==0) then
        print *
        print *, 'Unknow brilltype'
        print *, 'Aborting! (Subroutine idata_brillwave.f90)'
        print *
     end if

     call die

  end if

  
! ********************************
! ***   WAVE ELLIPTIC SOLVER   ***
! ********************************

  if (ELL_solver=="wave") then

     if (rank==0) then
        print *, 'Calculating Brill wave initial data with WaveElliptic solver ...'
     end if

     call wave_elliptic('flat','one')

     if (rank==0) then
        print *, 'Done'
        print *
     end if


! *******************************
! ***   SOR ELLIPTIC SOLVER   ***
! *******************************

  else if (ELL_solver=="sor") then

     if (rank==0) then
        print *, 'Calculating Brill wave initial data with SOR elliptic solver ...'
     end if

     call sor_elliptic('flat','one')

     if (rank==0) then
        print *, 'Done'
        print *
     end if


! *************************
! ***   UNKNOWN SOLVER  ***
! *************************

  else

     if (rank==0) then
        print *
        print *, 'Unknow elliptic solver'
        print *, 'Aborting! (Subroutine idata_brillwave.f90)'
        print *
     end if

     call die

  end if


! *******************************************************************
! ***   COPY SOLUTION ONTO PSI AND RECONSTRUCT CONFORMAL METRIC   ***
! *******************************************************************

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

!       Conformal metric coefficients.

        A = exp(two*brillq)
        B = A

!       Regularization variable.

        lambda = (A - H)/r**2

!       Delta functions. Easily calculated since
!       the initial metric is diagonal with A=B
!       and H=1.

        Delta_r = (one - one/A)/r
        Delta_z = zero

     end do
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine idata_brillwave

