
  subroutine analysis_geometry(box,level)

! ********************************************************
! ***   CALCULATION OF ANALYSIS VARIABLES FOR OUTPUT   ***
! ********************************************************

! Include modules.

  use mpi
  use arrays
  use param
  use derivatives
  use procinfo

! Extra variables.

  implicit none

  logical flag

  integer i,j,k,l,m,n
  integer box,level

  real(8) interp
  real(8) rmax,zmax
  real(8) massr,massz,massd,aux
  real(8) third,half,one,two,smallpi

  real(8) integral
  real(8) Delr,Delz,Vr,Vz,Q1,Q2
  real(8) gu(1:3,1:3),Delta(1:3,1:3,1:3)


! *******************
! ***   NUMBERS   ***
! *******************

  third = 1.d0/3.d0
  half = 0.5d0

  one = 1.d0
  two = 2.d0
  smallpi = acos(-one)


! ************************************
! ***   PSEUDO SCWARZSCHILD MASS   ***
! ************************************

! The pseudo-Schwarzschild mass is obtained by noticing
! that for the Schwarzschild metric we have:
!
! g   =  1 / (1 - 2M/R)
!  RR
!
! with R the Schwarzschild or areal radius. This implies:
!
!                                             2         4
! M  =  R/2 ( 1 - 1/g  )  =  R/2 ( 1 - (dR/dr) / grr psi )
!                    RR
!
! with r the coordinate radius.  Notice also that in general:
!
!
!  d   =  ( rho d    +  z d  ) / r
!   r            rho       z
!
! But remember that the code uses r for "rho" and "rr" for r.
!
! We now take the areal radius as:
!
! R  =  psi**2 sqrt(gtt)
!
! where gtt is the angular conformal metric: gtt := g_(theta,theta).
! (The code calculates both grr and gtt in auxiliary_geometry.f90)
!
! Note that taking this value of R we are assuming spherical symmetry,
! which will only be true sufficiently far away.

! Find areal radius and save it in "auxarray".

  auxarray = psi2*sqrt(gtt)
  !auxarray = psi2*sqrt(abs(rr/r))*(gtt*gpp)**0.25d0 ! This is another possible way of defining R.

! Find derivatives of areal radius.

  diffvar => auxarray
  Dr_auxarray = diff1r(+1)
  Dz_auxarray = diff1z(+1)

! Calculate Schwarzschild mass.

  mass_sch = half*auxarray*(one - ((r*Dr_auxarray + z*Dz_auxarray)/rr)**2/(grr*psi4))

! Output Schwarzschild mass at boundaries.

  if ((box==0).and.(level==0)) then

     if (t(0,0)==0.d0) then

!       Single processor run.

        if (size==1) then

           massr = mass_sch(Nr, 0)
           massz = mass_sch( 0,Nz)
           massd = mass_sch(Nr,Nz)

!       Parallel run.

        else

           interpvar => mass_sch

           rmax = dble(Nrtotal-1)*dr0
           zmax = dble(Nztotal-1)*dz0

!          Mass along r direction.

           aux = interp(box,level,rmax,0.d0,flag)
           call MPI_Allreduce(aux,massr,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!          Mass along z direction.

           aux = interp(box,level,0.d0,zmax,flag)
           call MPI_Allreduce(aux,massz,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!          Mass along diagonal.

           aux = interp(box,level,rmax,zmax,flag)
           call MPI_Allreduce(aux,massd,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

        end if

!       Show results.

        if (rank==0) then
           print *
           write (*,'(A,ES13.6         )') ' Schwarzschild mass along r direction = ',massr
           write (*,'(A,ES13.6         )') ' Schwarzschild mass along z direction = ',massz
           write (*,'(A,ES13.6         )') ' Schwarzschild mass along diagonal    = ',massd
           print *
        end if

     end if

  end if


! ************************************
! ***   ADM MASS VOLUME INTEGRAL   ***
! ************************************

! ADM mass version 1.  For a momentarily  metric
! (with K_ij = 0) the ADM mass can be expressed as the
! following volume integral:
!
!                      / 
!    mass_ADM  =  2 pi | rho psi**5 r dr dz
!                      /
!
!                      /    i
!              +  1/8  | [ V  d  phi  +  Q1 - Q2 ] psi r sqrt(hdet) dr dz
!                      /       i
!
! where:
!
!     i     i     i
!    V  =  D  -  Z
!
! with:
!
!         i      mn      i
!    Delta   =  g   Delta
!                        mn
!
!     i                     im
!    Z       =  1/(2*hdet) g   d hdet
!                               m
!
!    Z       =  1/(2*hdet) d hdet
!     i                     i
!
! and:
!                m
!    Q1  =  Delta  Z
!                   m

!             mn     a        b
!    Q2  =  g   Delta    Delta
!                    mb       na
!
! and where the Delta^i_jk are just the difference between the
! conformal Christoffel symbols and the flat Christoffel symbols.
!
!         i           i       ^  i
!    Delta    =  Gamma   -  Gamma
!         jk          jk         jk
!
! Notice that the D^i above are just the Delta^i of the code.
!
! At the moment the volume integrals for the ADM mass
! are only done on the coarse grid.

  if ((box==0).and.(level==0)) then

!    Geometric contribution from V and Q1. These do not depend
!    on the angular momentum since neither d_i(phi), not Z_i
!    have angular components.

     do j=1-ghost,Nzmaxl(box)
        do i=1-ghost,Nrmaxl(box)

!          Find (Vr,Vz).

           Vr = - Dr_g_A(i,j) - r(i,j)*Dz_g_C(i,j) - r(i,j)*g_lambda(i,j) &
              - 1.d0/hdet(i,j)*(g_A(i,j)*Dr_hdet(i,j) + r(i,j)*g_C(i,j)*Dz_hdet(i,j))

           Vz = - Dz_g_B(i,j) - r(i,j)*Dr_g_C(i,j) - two*g_C(i,j) &
              - 1.d0/hdet(i,j)*(r(i,j)*g_C(i,j)*Dr_hdet(i,j) + g_B(i,j)*Dz_hdet(i,j))

!          Find (Delr,Delz) and Q2.

           Delr = - Dr_g_A(i,j) - r(i,j)*Dz_g_C(i,j) - r(i,j)*g_lambda(i,j) &
                - half/hdet(i,j)*(g_A(i,j)*Dr_hdet(i,j) + r(i,j)*g_C(i,j)*Dz_hdet(i,j))

           Delz = - Dz_g_B(i,j) - r(i,j)*Dr_g_C(i,j) - two*g_C(i,j) &
                - half/hdet(i,j)*(r(i,j)*g_C(i,j)*Dr_hdet(i,j) + g_B(i,j)*Dz_hdet(i,j))

           Q1 =  half/hdet(i,j)*(Delr*Dr_hdet(i,j) + Delz*Dz_hdet(i,j))

!          Add contributions to integrand.

           auxarray(i,j) = 0.125d0*(Vr*Dr_phi(i,j) + Vz*Dz_phi(i,j) + Q1)*psi(i,j)*r(i,j)*sqrt(hdet(i,j))

        end do
     end do

!    Geometric contributions from Q2.  This term
!    does depend on the angular momentum.

     if (angmom) then

!       Case with angular momentum.

        do j=1-ghost,Nzmaxl(box)
          do i=1-ghost,Nrmaxl(box)

!            Inverse metric.

             gu(1,1) = g_A(i,j)
             gu(2,2) = g_B(i,j)
             gu(3,3) = g_H(i,j)/r(i,j)**2
             gu(1,2) = r(i,j)*g_C(i,j)
             gu(1,3) = r(i,j)*g_C1(i,j)
             gu(2,3) = g_C2(i,j)

             gu(2,1) = gu(1,2)
             gu(3,1) = gu(1,3)
             gu(3,2) = gu(2,3)

!            Deltas. Most of them are equal to the Christoffels,
!            except Delta(1,3,3) and Delta(3,1,3).

             Delta(1,1,1) = chris_rrr(i,j)
             Delta(1,1,2) = chris_rrz(i,j)
             Delta(1,1,3) = chris_rrp(i,j)
             Delta(1,2,2) = chris_rzz(i,j)
             Delta(1,2,3) = chris_rzp(i,j)
             Delta(1,3,3) = chris_rpp(i,j) + r(i,j)

             Delta(2,1,1) = chris_zrr(i,j)
             Delta(2,1,2) = chris_zrz(i,j)
             Delta(2,1,3) = chris_zrp(i,j)
             Delta(2,2,2) = chris_zzz(i,j)
             Delta(2,2,3) = chris_zzp(i,j)
             Delta(2,3,3) = chris_zpp(i,j)

             Delta(3,1,1) = chris_prr(i,j)
             Delta(3,1,2) = chris_prz(i,j)
             Delta(3,1,3) = chris_prp(i,j) - one/r(i,j)
             Delta(3,2,2) = chris_pzz(i,j)
             Delta(3,2,3) = chris_pzp(i,j)
             Delta(3,3,3) = chris_ppp(i,j)

!            Symmetries.

             Delta(1,2,1) = Delta(1,1,2)
             Delta(2,2,1) = Delta(2,1,2)
             Delta(3,2,1) = Delta(3,1,2)

             Delta(1,3,1) = Delta(1,1,3)
             Delta(2,3,1) = Delta(2,1,3)
             Delta(3,3,1) = Delta(3,1,3)

             Delta(1,3,2) = Delta(1,2,3)
             Delta(2,3,2) = Delta(2,2,3)
             Delta(3,3,2) = Delta(3,2,3)

!            Now do the sum for Q2.

             Q2 = 0.d0

             do k=1,3
                do l=1,3
                   do m=1,3
                      do n=1,3
                         Q2 = Q2 + gu(m,n)*Delta(k,m,l)*Delta(l,n,k)
                      end do
                   end do
                end do
             end do

!            Add contributions to integrand.

             auxarray(i,j) = auxarray(i,j) - 0.125d0*Q2*psi(i,j)*r(i,j)*sqrt(hdet(i,j))

          end do
        end do

     else

!       No angular momentum.

        gu = 0.d0
        Delta = 0.d0

        do j=1-ghost,Nzmaxl(box)
          do i=1-ghost,Nrmaxl(box)

!            Inverse metric.

             gu(1,1) = g_A(i,j)
             gu(2,2) = g_B(i,j)
             gu(3,3) = g_H(i,j)/r(i,j)**2
             gu(1,2) = r(i,j)*g_C(i,j)

             gu(2,1) = gu(1,2)

!            Deltas. Most of them are equal to the Christoffels,
!            except Delta(1,3,3) and Delta(3,1,3).

             Delta(1,1,1) = chris_rrr(i,j)
             Delta(1,1,2) = chris_rrz(i,j)
             Delta(1,2,2) = chris_rzz(i,j)
             Delta(1,3,3) = chris_rpp(i,j) + r(i,j)

             Delta(2,1,1) = chris_zrr(i,j)
             Delta(2,1,2) = chris_zrz(i,j)
             Delta(2,2,2) = chris_zzz(i,j)
             Delta(2,3,3) = chris_zpp(i,j)

             Delta(3,1,3) = chris_prp(i,j) - one/r(i,j)
             Delta(3,2,3) = chris_pzp(i,j)

!            Symmetries.

             Delta(1,2,1) = Delta(1,1,2)
             Delta(2,2,1) = Delta(2,1,2)
             Delta(3,3,1) = Delta(3,1,3)
             Delta(3,3,2) = Delta(3,2,3)

!            Now do the sum for Q2.

             Q2 = 0.d0

             do k=1,3
                do l=1,3
                   do m=1,3
                      do n=1,3
                         Q2 = Q2 + gu(m,n)*Delta(k,m,l)*Delta(l,n,k)
                      end do
                   end do
                end do
             end do

!            Add contributions to integrand.

             auxarray(i,j) = auxarray(i,j) - 0.125d0*Q2*psi(i,j)*r(i,j)*sqrt(hdet(i,j))

          end do
        end do

     end if

!    Add matter contribution to integrand

     if (mattertype/="vacuum") then
        auxarray = auxarray + (2.d0*smallpi*r)*rho*psi**5
     end if

!    Integrate.

     mass_ADM_V = integral(0,0,auxarray)

     if (rank==0) then
        write (*,'(A,ES13.6)') ' Total integrated ADM mass   = ',mass_ADM_V
        print *
     end if

  end if


! *********************************
! ***   CHARACTERISTIC SPEEDS   ***
! *********************************

! Light speed.

  vl_rp = + abs(alpha)*sqrt(abs(g_A/psi4))
  vl_rm = - abs(alpha)*sqrt(abs(g_A/psi4))

  vl_zp = + abs(alpha)*sqrt(abs(g_B/psi4))
  vl_zm = - abs(alpha)*sqrt(abs(g_B/psi4))

! Slicing speeds.

  if (slicing=='maximal') then

     va_rp = vl_rp
     va_rm = vl_rm

     va_zp = vl_zp
     va_zm = vl_zm

  else if ((index(slicing,"harmonic")/=0).or. &
           (index(slicing,"1+log")/=0).or. &
           (index(slicing,"shockavoid")/=0)) then

     va_rp = + sqrt(abs(falpha*g_A/psi4))
     va_rm = - sqrt(abs(falpha*g_A/psi4))

     va_zp = + sqrt(abs(falpha*g_B/psi4))
     va_zm = - sqrt(abs(falpha*g_B/psi4))

  end if

! Shift speeds.

  if ((shift(1:11)=="Gammadriver").and.(drivercsi/=0.d0)) then

     vs_rp = + sqrt(abs(4.d0/3.d0*drivercsi*g_A/psi4))
     vs_rm = - sqrt(abs(4.d0/3.d0*drivercsi*g_A/psi4))

     vs_zp = + sqrt(abs(4.d0/3.d0*drivercsi*g_B/psi4))
     vs_zm = - sqrt(abs(4.d0/3.d0*drivercsi*g_B/psi4))

  end if

! Add shift contribution.

  if (shift/="none") then

     vl_rp = - beta_r + vl_rp
     vl_rm = - beta_r + vl_rm
     vl_zp = - beta_z + vl_zp
     vl_zm = - beta_z + vl_zm

     va_rp = - beta_r + va_rp
     va_rm = - beta_r + va_rm
     va_zp = - beta_z + va_zp
     va_zm = - beta_z + va_zm

     vs_rp = - beta_r + vs_rp
     vs_rm = - beta_r + vs_rm
     vs_zp = - beta_z + vs_zp
     vs_zm = - beta_z + vs_zm

  end if


! ****************************************
! ***   WEYL TENSOR AND WAVE EXTRACT   ***
! ****************************************

! Calculate Weyl tensor and curvature invariants.

  if (curvInv.or.wave_extract) then
     call weyl
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine analysis_geometry

