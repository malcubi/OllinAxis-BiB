!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/cfl.f90,v 1.19 2022/10/24 21:57:59 malcubi Exp $

  subroutine cfl(step)

! ****************************
! ***   ADJUST TIME STEP   ***
! ****************************

! This subroutine modifies the time step in order
! to guarantee that the CFL condition is satisfied.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical firstcall         ! Is this the first call?

  integer i,j               ! Counters.
  integer step              ! Time step.
  integer box,level,proc    ! Box, level and processor counters.

  real(8) vrmax,vzmax       ! Local maximum speeds.
  real(8) aux               ! Auxiliary.

  real(8) vr_max(0:Nb,0:Nlmax),vz_max(0:Nb,0:Nlmax) ! Arrays with maximum speeds for different boxes and levels.

  character(20) filestatus

  data firstcall / .true. /


! **************************************
! ***   LOOP OVER BOXES AND LEVELS   ***
! **************************************

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))


!       *************************************************
!       ***   FIND CHARACTERISTIC SPEEDS IN R AND Z   ***
!       *************************************************

!       Light speed:
!
!       vl = alpha sqrt(g^ii / psi^4)

        vl_rp = sqrt(alpha*abs(g_A))/psi2
        vl_zp = sqrt(alpha*abs(g_B))/psi2

!       Slicing speeds.  For maximal slicing this
!       is just equal to the speed of light.  For
!       Bona-Masso type slicings it is:
!
!       va = +- alpha sqrt( f(alpha) g^ii / psi^4 )
!
!       But remember that the code defines:
!
!       falpha := alpha^2 f(alpha)

        if (slicing=='maximal') then

           va_rp = vl_rp
           va_zp = vl_zp

        else if ((index(slicing,"harmonic")/=0).or. &
                 (index(slicing,"1+log")/=0).or. &
                 (index(slicing,"shockavoid")/=0)) then

           va_rp = sqrt(abs(falpha*g_A))/psi2
           va_zp = sqrt(abs(falpha*g_B))/psi2

        end if

!       Shift speeds.  For Gammadriver type shifts we have:
!
!       vs = +- sqrt( (4/3) csi g^ii / psi^4)

        if ((shift(1:11)=="Gammadriver").and.(drivercsi/=0.d0)) then

           vs_rp = sqrt(abs(4.d0/3.d0*drivercsi*g_A/psi4))
           vs_zp = sqrt(abs(4.d0/3.d0*drivercsi*g_B/psi4))

        end if

!       Add shift contribution. With a non-zero shift,
!       one must add abs(beta^i) to ALL characteristic speeds.

        if (shift/="none") then

           vl_rp = abs(beta_r) + vl_rp
           vl_zp = abs(beta_z) + vl_zp

           va_rp = abs(beta_r) + va_rp
           va_zp = abs(beta_z) + va_zp

           vs_rp = abs(beta_r) + vs_rp
           vs_zp = abs(beta_z) + vs_zp

        end if


!       ******************************
!       ***   FIND MAXIMUM SPEED   ***
!       ******************************

!       Find maximum speeds in r and z.

        vrmax = 0.d0
        vzmax = 0.d0

        if (shift=="none") then

           do j=1-ghost,Nz
              do i=1-ghost,Nr

                 if (abs(vl_rp(i,j))>vrmax) vrmax=abs(vl_rp(i,j))
                 if (abs(vl_zp(i,j))>vzmax) vzmax=abs(vl_zp(i,j))

                 if (abs(va_rp(i,j))>vrmax) vrmax=abs(va_rp(i,j))
                 if (abs(va_zp(i,j))>vzmax) vzmax=abs(va_zp(i,j))

              end do
           end do

        else

           do j=1-ghost,Nz
              do i=1-ghost,Nr

                 if (abs(vl_rp(i,j))>vrmax) vrmax=abs(vl_rp(i,j))
                 if (abs(vl_zp(i,j))>vzmax) vzmax=abs(vl_zp(i,j))

                 if (abs(va_rp(i,j))>vrmax) vrmax=abs(va_rp(i,j))
                 if (abs(va_zp(i,j))>vzmax) vzmax=abs(va_zp(i,j))

                 if (abs(vs_rp(i,j))>vrmax) vrmax=abs(vs_rp(i,j))
                 if (abs(vs_rm(i,j))>vrmax) vrmax=abs(vs_rm(i,j))

              end do
           end do

        end if

!       For parallel runs find the global maximum.

        if (size>1) then

           call MPI_ALLREDUCE(vrmax,aux,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
           vrmax = aux

           call MPI_ALLREDUCE(vzmax,aux,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
           vzmax = aux

        end if

!       Save maximum speeds for current box and level.

        vr_max(box,level) = vrmax
        vz_max(box,level) = vzmax

     end do
  end do


! *******************************************************
! ***   FIND MAXIMUM SPEEDS ACROSS BOXES AND LEVELS   ***
! *******************************************************

  vrmax = 0.d0
  vzmax = 0.d0

  do box=0,Nb
     do level=min(1,box),Nl(box)
        if (abs(vr_max(box,level))>vrmax) vrmax=abs(vr_max(box,level))
        if (abs(vz_max(box,level))>vzmax) vzmax=abs(vz_max(box,level))
     end do
  end do


! ****************************
! ***   MODIFY TIME STEP   ***
! ****************************

! Calculate the next dt using the maximum value of the charactersitic
! speeds in r and z. But don't allow it to become too large.

  dt0 = dtfac*min(dr0,dz0,dr0/vrmax,dz0/vzmax)

  do level=0,Nlmax
     dtl(level) = dt0/2**level
  end do


! *********************************
! ***   SAVE TIME STEP TO FILE  ***
! *********************************

! On first call, replace file. Otherwise just append to it.

  if (firstcall) then
     firstcall = .false.
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Save value of dt0 to file.

  if (rank==0) then

     if (filestatus == 'replace') then
        open(1,file=trim(directory)//'/dt0.tl',form='formatted', &
             status=filestatus)
        write(1,*) '"dt0.tl'
     else
        open(1,file=trim(directory)//'/dt0.tl',form='formatted', &
                status=filestatus,position='append')
     end if
  
     write(1,"(2ES14.6)") t(0,0),dt0
     close(1)

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine cfl

