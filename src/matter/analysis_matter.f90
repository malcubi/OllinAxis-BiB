
  subroutine analysis_matter(box,level)

! ********************************************************
! ***   CALCULATION OF ANALYSIS VARIABLES FOR OUTPUT   ***
! ********************************************************

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical contains

  integer box,level

  real(8) rmax,zmax
  real(8) massr,massz
  real(8) complex_NB
  real(8) aux,bind


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

  if (contains(mattertype,"complex")) then

!    Norm of complex scalar field: sqrt(phiR**2 + phiI**2).

     if (associated(grid(box,level)%complex_phi_norm)) then
        complex_phi_norm = dsqrt(complex_phiR**2 + complex_phiI**2)
     end if

!    The integrated boson charge is:
!
!                   /             
!    complex_NB  =  |  complex_Bdens dV
!                   /
!
!    with dV the physical volumen element.
!
!    At the moment I only do the integral in the coarse
!    grid, and at t=0.

     if ((box==0).and.(level==0)) then

        if ((t(0,0)==0.d0).and.(idata=="bosonstar")) then

           call bosonintegral(complex_NB)

           if (rank==0) then
              write (*,'(A,ES12.5)') ' Total boson number NB = ',complex_NB
              print *
           end if

!          Output the total binding energy, but only
!          for certain type of initial data.

           if (size==1) then

              massr = mass_sch(Nr,0)
              massz = mass_sch(0,Nz)

           else

              rmax = r(Nr,0)
              zmax = z(0,Nz)

              if (ownequator.and.(rmax>dble(Nrtotal-1)*dr)) then
                 massr = mass_sch(Nr,0)
              else
                 massr = 0.d0
              end if

              if (ownaxis.and.(zmax>dble(Nztotal-1)*dz)) then
                 massz = mass_sch(0,Nz)
              else
                 massz = 0.d0
              end if

              aux = massr
              call MPI_ALLREDUCE(aux,massr,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

              aux = massz
              call MPI_ALLREDUCE(aux,massz,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

           end if

           if (rank==0) then
              bind = 0.5d0*(massr+massz) - complex_mass*complex_NB
              write(*,'(A,ES23.16)') ' Binding energy (M_sch - m*NB) = ',bind
              print *
           end if

        end if

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine analysis_matter

