!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/setupboxes.f90,v 1.20 2021/02/26 19:25:51 malcubi Exp $

  subroutine setupboxes

! ***********************************
! ***   SET UP REFINEMENT BOXES   ***
! ***********************************

! This routine sets up location and size of each
! refinemnt box.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer box,level       ! Box and level counters.


! *********************************************
! ***   SET UP NUMBER OF BOXES AND LEVELS   ***
! *********************************************

! Sanity check.

  if ((Nb<0).or.(Nb>3)) then
     if (rank==0) then
        print *, 'The number of refinement boxes must be between 0 and 3'
        print *, 'Aborting! (subroutine setupboxes.f90)'
        print *
     end if
     call die
  end if

! Set number of levels for each refinement box.

  Nl(0) = Nl0
  Nl(1) = Nl1
  Nl(2) = Nl2
  Nl(3) = Nl3

  Nlmax = maxval(Nl)

! Output number of refinement boxes.

  if (Nb==0) then
     if (Nl(0)==0) then
        if (rank==0) then
           print *, 'No refinement levels, only base grid'
           print *
        end if
     else
        if (rank==0) then
           print *, 'One refinement box centered on the origin'
           write (*,'(A,I2)') ' Total grid levels (Nl+1) = ',Nl(0)+1
           print *
        end if
     end if
  else if (Nlmax>0) then
     if (rank==0) then
        write (*,'(A,I2)') ' Number of refinement boxes =',Nb
           do box=0,Nb
              write (*,'(A,I2,A,I2)') ' Box =',box,'  Nl =',Nl(box)
           end do
        print *
     end if
  else
     if (rank==0) then
        print *, 'If you choose Nb>0 you must also have Nl>0 for at least one box'
        print *, 'Aborting! (subroutine setupboxes.f90)'
        print *
     end if
     call die
  end if


! *************************************************************
! ***   ALLOCATE ARRAYS FOR GRID SPACING AND TIME STEPPING  ***
! *************************************************************

! Arrays for time and number of time steps.

  allocate(s(0:Nb,0:Nlmax),t(0:Nb,0:Nlmax))

! Arrays for time at previous time levels.ç

  allocate(t1(0:Nb,0:Nlmax),t2(0:Nb,0:Nlmax))

! Arrays for grid spacing and time stepping.

  allocate(drl(0:Nlmax),dzl(0:Nlmax),dtl(0:Nlmax))


! *****************************************************
! ***   GRID SPACING AND TIME STEP FOR EACH LEVEL   ***
! *****************************************************

! Sanity check.

  if ((dr/=0.d0).or.(dz/=0.d0).or.(dt/=0.d0)) then

     if (rank==0) then
        print *
        print *, 'The values of (dr,dz,dt)  should not be fixed in the parameter file'
        print *, 'Aborting! (subroutine setupboxes.f90)'
        print *
     end if

     call die

  end if

! Find grid spacing and time step for each level.

  dt0 = dtfac*min(dr0,dz0)

  do level=0,Nlmax
     dtl(level) = dt0/2**level
     drl(level) = dr0/2**level
     dzl(level) = dz0/2**level
  end do


! *********************************************
! ***   SET UP LOCATION AND SIZE OF BOXES   ***
! *********************************************

! If there is no equatorial symmetry, make sure that
! the total number of grid points for box 0, including
! ghost zones, is an odd number.  This guarantees that
! there is always a grid point at the equator.

  if (.not.eqsym) then
     if (mod(Nztotal+ghost,2)==0) then
        Nztotal = Nztotal + 1
     end if
  end if

! Set size for each refinement box.

  Nrbox(0) = Nrtotal
  Nrbox(1) = Nrbox1
  Nrbox(2) = Nrbox2
  Nrbox(3) = Nrbox3

  Nzbox(0) = Nztotal
  Nzbox(1) = Nzbox1
  Nzbox(2) = Nzbox2
  Nzbox(3) = Nzbox3

! Sanity check.

  do box=0,Nb
     if ((Nrbox(box)<2+ghost).or.(Nzbox(box)<2+ghost)) then
        if (rank==0) then
           print *, 'The values of (Nrtotal,Nztotal,Nrbox,Nzbox) should'
           print *, 'be greater than 2 plus the number of ghost zones'
           print *, 'Aborting! (subroutine setupboxes.f90)'
           print *
        end if
        call die
     end if
  end do

! Set center of box 0.

  rbox(0) = 0.d0
  zbox(0) = 0.d0

! Set center of refinement boxes.

  rbox(1) = rbox1
  rbox(2) = rbox2
  rbox(3) = rbox3

  zbox(1) = zbox1
  zbox(2) = zbox2
  zbox(3) = zbox3

! Now we need to adjust the position of these boxes
! so that the fine grid points will always coincide
! with those of box 0 if they where to overlap.

  do box=1,Nb

!    Since we stagger the axis, rbox(box) must always
!    be a half-integer multiple of drl(1).

     if (rbox(box)/=0.d0) then
        rbox(box) = (int(rbox(box)/drl(1)) + 0.5d0)*drl(1)
     end if

!    For equatorial symmetry we stagger the equator, so again
!    zbox(box) must always be a half-integer multiple of dzl(1).
!    On the other hand, if we don't have equatorial symmetry
!    there is always a frid point at the equator, so zbox(box)
!    must be an integer multiple of dzl(1).
 
     if (zbox(box)/=0.d0) then
        if (eqsym) then
           zbox(box) = (int(zbox(box)/dzl(1)) + 0.5d0)*dzl(1)
        else
           zbox(box) = int(zbox(box)/dzl(1))*dzl(1)
        end if
     end if

  end do

! Sanity check.

  do box=1,Nb

!    Check refinement boxes are contained within base grid in r direction.

     if (rbox(box)<0.d0) then
        if (rank==0) then
           print *, 'rbox must be larger than or equal to 0 for all refinement boxes'
           print *, 'Aborting! (subroutine setupboxes.f90)'
           print *
        end if
        call die
     else if ((rbox(box)>0.d0).and.(rbox(box)-0.5d0*Nrbox(box)*drl(1)<=0.d0)) then
        if (rank==0) then
           write(*,'(A,I2,A)') ' Refinement box',box,' extends to negative values of r'
           print *, 'Aborting! (subroutine setupboxes.f90)'
           print *
        end if
        call die
     else if (rbox(box)+0.5d0*Nrbox(box)*drl(1)>=(Nrbox(0)-ghost)*dr0) then
        if (rank==0) then
           write(*,'(A,I2,A)') ' Refinement box',box,' extends beyond base grid in r direction'
           print *, 'Aborting! (subroutine setupboxes.f90)'
           print *
        end if
        call die
     end if

!    Check refinement boxes are contained within base grid in z direction.

     if (eqsym) then

        if (zbox(box)<0.d0) then
           if (rank==0) then
              print *, 'For equatorial symmetry zbox must be larger than or equal to 0 for all refinement boxes'
              print *, 'Aborting! (subroutine setupboxes.f90)'
              print *
           end if
           call die
        else if ((zbox(box)>0.d0).and.(zbox(box)-0.5d0*Nzbox(box)*dzl(1)<=0.d0)) then
           if (rank==0) then
              write(*,'(A,I2,A)') ' Refinement box',box,' extends to negative values of z'
              print *, 'Aborting! (subroutine setupboxes.f90)'
              print *
           end if
           call die
        else if (zbox(box)+0.5d0*Nzbox(box)*dzl(1)>=(Nzbox(0)-ghost)*dz0) then
           if (rank==0) then
              write(*,'(A,I2,A)') ' Refinement box',box,' extends beyond base grid in z direction'
              print *, 'Aborting! (subroutine setupboxes.f90)'
              print *
           end if
           call die
        end if

     else

        if (zbox(box)-0.5d0*Nzbox(box)*dzl(1)<=-0.5d0*(Nzbox(0)-ghost)*dz0) then
           if (rank==0) then
              write(*,'(A,I2,A)') ' Refinement box',box,' extends beyond base grid in z direction'
              print *, 'Aborting! (subroutine setupboxes.f90)'
              print *
           end if
           call die
        else if (zbox(box)+0.5d0*Nzbox(box)*dzl(1)>=0.5d0*(Nzbox(0)-ghost)*dz0) then
           if (rank==0) then
              write(*,'(A,I2,A)') ' Refinement box',box,' extends beyond base grid in z direction'
              print *, 'Aborting! (subroutine setupboxes.f90)'
              print *
           end if
           call die
        end if

     end if

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine setupboxes

