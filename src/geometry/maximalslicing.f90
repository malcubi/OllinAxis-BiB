!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/geometry/maximalslicing.f90,v 1.7 2021/02/24 20:53:01 malcubi Exp $

  subroutine maximalslicing

! ***************************
! ***   MAXIMAL SLICING   ***
! ***************************

! Ths subroutine solves the maximal slicing condition:
!
! __2                        ij
! \/ alpha   =  alpha [ K   K    +   4 pi ( rho + S ) ]
!                        ij

! Include modules.

  use param
  use arrays
  use procinfo
  use derivatives

! Extra variables.

  implicit none

  integer box,level       ! Box number and level counters.

  real(8) zero,one,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  one  = 1.d0

  smallpi = acos(-one)


! ******************************************
! ***   ELLIPTIC EQUATION COEFFICIENTS   ***
! ******************************************

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       For maximal slicing S0=S5=0.

        ell_S0 = zero
        ell_S5 = zero

!       Coefficient of linear term for elliptic equation.
!       (Remember that it goes on the same side as the Laplacian,
!       so beware the sign).

        ell_S1 = - K2

!       Matter terms.

        if (mattertype/="vacuum") then
           ell_S1 = ell_S1 - 4.d0*smallpi*(rho + trS)
        end if

     end do
  end do


! *******************************
! ***   WAVE ELLIPTIC SOLVER  ***
! *******************************

  if (ELL_solver=="wave") then

     if (rank==0) then
        print *
        print *, 'Calculating maximal slicing with WaveElliptic solver'
     end if

     call wave_elliptic('physical','one')


! *******************************
! ***   SOR ELLIPTIC SOLVER   ***
! *******************************

  else if (ELL_solver=="sor") then

     if (rank==0) then
        print *
        print *, 'Calculating maximal slicing with SOR elliptic solver'
     end if

     call sor_elliptic('physical','one')


! *************************
! ***   UNKNOWN SOLVER  ***
! *************************

  else

     if (rank==0) then
        print *
        print *, 'Unknow elliptic solver'
        print *, 'Aborting! (Subroutine maximalslicing.f90)'
        print *
     end if

     call die

  end if


! ************************************
! ***   COPY SOLUTION ONTO ALPHA   ***
! ************************************

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Copy ell_u into alpha.

        alpha = ell_u

     end do
  end do


! ********************************
! ***   DERIVATIVES OF LAPSE   ***
! ********************************

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Derivatives of lapse.

        diffvar => alpha

        Dr_alpha = diff1r(+1)
        Dz_alpha = diff1z(+1)

        Drr_alpha = diff2r(+1)
        Dzz_alpha = diff2z(+1)
        Drz_alpha = diff2rz(+1,+1)

!       Regularization.
  
        auxarray = Dr_alpha/r

        diffvar => auxarray
        DD_alphar = diff1r(+1)

    end do
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine maximalslicing
