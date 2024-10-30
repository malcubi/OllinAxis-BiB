!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/save1Ddata.f90,v 1.24 2021/03/28 00:44:46 malcubi Exp $

  subroutine save1Ddata

! ********************************
! ***   SAVE 1D DATA TO FILE   ***
! ********************************

! This routine saves 1D data, that is a 1D projection
! of the whole spatial arrays at a given time.
!
! Here, output is produced for the values in r-direction at the 
! equator (".rl") and in z-direction at the axis (".zl").

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical firstcall

  integer i,j,l
  integer box,level
  integer, allocatable, dimension (:) :: commas

  character(1)  aa
  character(20) filestatus

  data firstcall / .true. /

  save firstcall


! ***********************************
! ***   IS THIS THE FIRST CALL?   ***
! ***********************************

! On first call, replace file and figure
! out what needs output.

  if (firstcall) then

     firstcall = .false.

!    File status.

     filestatus = 'replace'

!    Find out length of string "outvars1D".

     l = len_trim(outvars1D)

!    Find out how many variables need output.

     nvars1D = 1

     do i=1,l
        aa = outvars1D(i:i)
        if (aa==",") then
           nvars1D = nvars1D + 1
        end if
     end do

!    Identify comma positions.

     allocate(commas(0:nvars1D))

     j = 0
     commas(0) = 0
     commas(nvars1D) = l+1

     do i=1,l
        aa = outvars1D(i:i)
        if (aa==",") then
           j = j+1
           commas(j)=i
        end if
     end do

!    Now read variable names, and eliminate spaces.

     allocate(outvars1Darray(1:nvars1D))

     do i=1,nvars1D
        outvars1Darray(i) = outvars1D(commas(i-1)+1:commas(i)-1)
        outvars1Darray(i) = trim(adjustl(outvars1Darray(i)))
     end do

!    Check if any name is repeated.

     if (rank==0) then
        do i=1,nvars1D
           do j=1,nvars1D
              if ((i.ne.j).and.(trim(outvars1Darray(i))==trim(outvars1Darray(j)))) then
                 print *, 'Error in parfile, array name repeated in outvars1D: ',outvars1Darray(i)
                 print *, 'Aborting! (subroutine save1Ddata)'
                 print *
                 call die
             end if
           end do
        end do
     end if

! Not first call.

  else

     filestatus = 'old'

  end if


! *********************
! ***   SAVE DATA   ***
! *********************

! Loop over boxes and grid levels.

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Save all variables that require output.

        do i=1,nvars1D
           call grabarray(trim(outvars1Darray(i)))
           call save1Dvariable(trim(outvars1Darray(i)),directory,box,level,outparallel,filestatus)
        end do

     end do
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine save1Ddata









  subroutine save1Dvariable(varname,outdir,box,level,type,filestatus)

! This subroutine outputs 1D files:
!
! 1) Files with extension .rl correspond to the equator.
!    For runs with equatorial symmetry the output is
!    in fact for the points with z=+dz/2.
!
! 2) Files with extension .zl correspond to the symmetry
!    axis (in fact points with r=+dr/2).
!
! 3) Files with extension .dl correspond to the diagonal
!    for runs with equatorial symmetry.  When there is
!    no equatorial symmetry two diagonals are output
!    with extensions (.d1l.d2l).
!
! Diagonal output is along the computational domain
! diagonals, and NOT along the (r=z,r=-z) lines.
!
! Notice that if the box is not a square in general there
! will not be grid points along the box diagonal, so we
! need to use interpolation for this output.

! Load modules.

  use mpi
  use param
  use arrays
  use procinfo

! Declare variables.

  implicit none

  logical flag

  integer i,j,p
  integer box,level
  integer i0,j0,imax,jmax,Naux,NN
  integer unit
  integer status(MPI_STATUS_SIZE)

  real(8) interp
  real(8) delta,rrmax,r0,z0,rr0
  real(8) aux1,aux2

  character(len=9) filen
  character(len=*) varname,outdir,type,filestatus

  real(8), dimension (1-ghost:Nrmax) :: varr,auxr
  real(8), dimension (1-ghost:Nzmax) :: varz,auxz


! *******************************
! ***   FILE NAME EXTENSION   ***
! *******************************

! Append box number and grid level to file name.

  write(filen,'(a,i1,a,i1)') '_b',box,'l',level

! For parallel runs that output one file per processor
! we also append the processor number.

  if ((size>1).and.(type=="fileperproc")) then
     if (size<10) then
        write(filen,'(a,i1,a,i1,a,i1)') '_b',box,'l',level,'p',rank
     else if (size<100) then
        write(filen,'(a,i1,a,i1,a,i2)') '_b',box,'l',level,'p',rank
     else
        write(filen,'(a,i1,a,i1,a,i3)') '_b',box,'l',level,'p',rank
     end if

  end if


! ************************************
! ***   SAVE VARIABLE ON EQUATOR   ***
! ************************************

  if (size==1) then

!    ********************************
!    ***   SINGLE PROCESSOR RUN   ***
!    ********************************

!    Save data only if the current box owns the equator.

     if (ownequator) then

        if (filestatus=='replace') then     
           open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.rl',form='formatted', &
           status='replace')
        else
           open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.rl',form='formatted', &
           status='old',position='append')
        end if

!       Write current time.

        write(1,"(A8,ES14.6)") '"Time = ',t(box,level)

!       Save data.

        do i=1-ghost,Nr
           write(1,"(2ES20.8E3)") r(i,eqz(box,0)),grabvar(i,eqz(box,0))
        end do

!       Leave blank space before next time.

        write(1,*)

!       Close file.

        close(1)

     end if


  else if (type=="singlefile") then

!    **************************************
!    ***   PARALLEL RUN - SINGLE FILE   ***
!    **************************************

!    Check if some processor in the current box owns the equator.

     do p=0,size-1
        if (eqz(box,p)/=-1) goto 10
     end do

     goto 20

!    Only processor 0 does output.

10   continue

     if (rank==0) then

!       Open file.

        if (filestatus=='replace') then     
           open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.rl',form='formatted', &
           status='replace')
        else
           open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.rl',form='formatted', &
           status='old',position='append')
        end if

!       Write current time.

        write(1,"(A8,ES14.6)") '"Time = ',t(box,level)

!       Save data from processor 0 if it owns the equator.

        if (eqz(box,0)>=1) then

           if (nprocr==1) then
              imax = Nrl(box,0)
           else
              imax = Nrl(box,0) - ghost
           end if

           do i=1-ghost,imax
              write(1,"(2ES20.8E3)") r(i,eqz(box,0)),grabvar(i,eqz(box,0))
           end do

        end if

!       Iterate over other processors.

        do p=1,size-1

!          Receive and save data from processors along equator
!          (be careful not to save ghost zones twice).

           if (eqz(box,p)>=1) then

              Naux = Nrmax+ghost
              call MPI_RECV(varr,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(auxr,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

              if (auxr(0)<0.d0) then
                 i0 = 1-ghost
              else
                 i0 = 1
              end if

              if (mod(p+1,nprocr)==0) then
                 imax = Nrl(box,p)
              else
                 imax = Nrl(box,p) - ghost
              end if

              do i=i0,imax
                 write(1,"(2ES20.8E3)") auxr(i),varr(i)
              end do

           end if

        end do

!       Leave blank space before next time.

        write(1,*)

!       Close file.

        close(1)

!    Other processors send data to processor 0.

     else

!       Only processors that own the equator send data.

        if (eqz(box,rank)>=1) then
           Naux = Nrmax+ghost
           varr = grabvar(:,eqz(box,rank))
           auxr = r(:,eqz(box,rank))
           call MPI_SEND(varr,Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(auxr,Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end if

     end if

     20 continue


  else

!    ********************************************
!    ***   PARALLEL RUN - ONE FILE PER PROC   ***
!    ********************************************

!    Each processor does its own output. First check it
!    the current processor owns the equator.

     if (ownequator) then

!       Open file.

        unit = rank+1

        if (filestatus=='replace') then     
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'.rl',form='formatted', &
           status='replace')
        else
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'.rl',form='formatted', &
           status='old',position='append')
        end if

!       Write current time.

        write(unit,"(A8,ES14.6)") '"Time = ',t(box,level)

!       Save data.

        do i=1-ghost,Nr
           write(unit,"(2ES20.8E3)") r(i,eqz(box,rank)),grabvar(i,eqz(box,rank))
        end do

!       Leave blank space before next time.

        write(unit,*)

!       Close file.

        close(unit)

     end if

  end if


! *********************************
! ***   SAVE VARIABLE ON AXIS   ***
! *********************************

  if (size==1) then

!    ********************************
!    ***   SINGLE PROCESSOR RUN   ***
!    ********************************

!    Save data only if the current box owns the axis.

     if (ownaxis) then

!       Open file.

        if (filestatus=='replace') then     
           open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.zl',form='formatted', &
           status='replace')
        else
           open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.zl',form='formatted', &
           status='old',position='append')
        end if

!       Write current time.

        write(1,"(A8,ES14.6)") '"Time = ',t(box,level)

!       Save data.

        do j=1-ghost,Nz
           write(1,"(2ES20.8E3)") z(1,j),grabvar(1,j)
        end do

!       Leave blank space before next time.

        write(1,*)

!       Close file.

        close(1)

     end if


  else if (type=="singlefile") then

!    *************************************
!    ***   PARALLEL RUN - SINGLE FILE  ***
!    *************************************

!    Check if some processor in the current box owns the axis.

     do p=0,size-1
        if (axis(box,p)/=-1) goto 30
     end do

     goto 40

!    Only processor 0 does output.

     30 continue

     if (rank==0) then

!       Open file.

        if (filestatus=='replace') then     
           open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.zl',form='formatted', &
          status='replace')
        else
           open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.zl',form='formatted', &
           status='old',position='append')
        end if

!       Write current time.

        write(1,"(A8,ES14.6)") '"Time = ',t(box,level)

!       Save data from processor 0.

        if (nprocz==1) then
           jmax = Nzl(box,0)
        else
           jmax = Nzl(box,0) - ghost
        end if

        do j=1-ghost,jmax
           write(1,"(2ES20.8E3)") z(1,j),grabvar(1,j)
        end do

!       Iterate over other processors.

        do p=1,size-1

!          Receive and save data from processors along axis
!          (be careful not to save ghost zones twice).

           if (mod(p,nprocr)==0) then

              Naux = Nzmax+ghost
              call MPI_RECV(varz,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(auxz,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

              if ((.not.eqsym).and.(p<nprocr)) then
                 j0 = 1-ghost
              else
                 j0 = 1
              end if

              if (p>=size-nprocr) then
                 jmax = Nzl(box,p)
              else
                 jmax = Nzl(box,p)-ghost
              end if

              do j=j0,jmax
                 write(1,"(2ES20.8E3)") auxz(j),varz(j)
              end do

           end if

        end do

!       Leave blank space before next time.

        write(1,*)

!       Close file.

        close(1)

!    Other processors send data to processor 0.

     else

        if (mod(rank,nprocr)==0) then
           Naux = Nzmax+ghost
           varz = grabvar(1,:)
           auxz = z(1,:)
           call MPI_SEND(varz,Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(auxz,Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end if

     end if

     40 continue


  else

!    ********************************************
!    ***   PARALLEL RUN - ONE FILE PER PROC   ***
!    ********************************************

!    Each processor does its own output. First check it
!    the current processor owns the axis.

     if (ownaxis) then

!       Open file.

        unit = rank+1

        if (filestatus=='replace') then     
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'.zl',form='formatted', &
           status='replace')
        else
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'.zl',form='formatted', &
           status='old',position='append')
        end if

!       Write current time.

        write(unit,"(A8,ES14.6)") '"Time = ',t(box,level)

!       Save data.

        do j=1-ghost,Nz
           write(unit,"(2ES20.8E3)") z(1,j),grabvar(1,j)
        end do

!       Leave blank space before next time.

        write(unit,*)

!       Close file.

        close(unit)

     end if

  end if


! **************************************
! ***   SAVE VARIABLE ON DIAGONALS   ***
! **************************************

! Only box 0 does output along diagonals. If there
! is no equatorial symmetry, we output two diagonals.
!
! Diagonal output is along the computational domain
! diagonals, and NOT along the line r=z.
!
! Notice that if the box is not a square in general there
! will not be grid points along the box diagonal, so we
! need to use interpolation for this output.

  if (box==0) then

     if (size==1) then

!       ********************************
!       ***   SINGLE PROCESSOR RUN   ***
!       ********************************

!       Equatorial symmetry.

        if (eqsym) then

!          Open file.

           if (filestatus=='replace') then     
              open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.dl',form='formatted', &
              status='replace')
           else
              open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.dl',form='formatted', &
              status='old',position='append')
           end if

!          Write current time.

           write(1,"(A8,ES14.6)") '"Time = ',t(box,level)

!          Save data.

           rrmax = sqrt(rmaxl(0,level)**2 + zmaxl(0,level)**2)
           delta = sqrt(drl(level)**2 + dzl(level)**2)

           NN = int(rrmax/delta + 0.5d0)
           delta = 1.d0/dble(NN)

           interpvar => grabvar

           do i=0,NN
              aux1 = (dble(i) - 0.5d0)*delta
              r0 = aux1*rmaxl(0,level)
              z0 = aux1*zmaxl(0,level)
              rr0 = aux1*rrmax
              aux2 = interp(box,level,r0,z0,flag)
              write(1,"(2ES20.8E3)") rr0,aux2
           end do

!          Leave blank space before next time.

           write(1,*)

!          Close file.

           close(1)

!       No equatorial symmetry.

        else

!          Open files.

           if (filestatus=='replace') then     
              open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.d1l',form='formatted', &
              status='replace')
              open(2,file=trim(outdir)//'/'//varname//trim(filen)//'.d2l',form='formatted', &
              status='replace')
           else
              open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.d1l',form='formatted', &
              status='old',position='append')
              open(2,file=trim(outdir)//'/'//varname//trim(filen)//'.d2l',form='formatted', &
              status='old',position='append')
           end if

!          Write current time.

           write(1,"(A8,ES14.6)") '"Time = ',t(box,level)
           write(2,"(A8,ES14.6)") '"Time = ',t(box,level)

!          Save data.

           rrmax = sqrt(rmaxl(0,level)**2 + zmaxl(0,level)**2)
           delta = sqrt(drl(level)**2 + dzl(level)**2)

           NN = int(rrmax/delta + 0.5d0)
           delta = 1.d0/dble(NN)

           interpvar => grabvar

           do i=0,NN
              aux1 = (dble(i)-0.5d0)*delta
              r0 = aux1*rmaxl(0,level)
              z0 = aux1*zmaxl(0,level)
              rr0 = aux1*rrmax
              aux2 = interp(box,level,r0,z0,flag)
              write(1,"(2ES20.8E3)") rr0,aux2
              z0 = - z0
              aux2 = interp(box,level,r0,z0,flag)
              write(2,"(2ES20.8E3)") rr0,aux2
           end do

!          Leave blank space before next time.

           write(1,*)
           write(2,*)

!          Close files.

           close(1)
           close(2)

        end if

     else

!       ************************
!       ***   PARALLEL RUN   ***
!       ************************

!       For the diagonal the output for parallel runs is
!       always a single file.  This is because the diagonal
!       may intercept many processors and sometimes with
!       just a couple of points, so it makes no sense to
!       output a separate file per processor.
!
!       Notice also that only processor 0 does output, but
!       all processors must call the interpolating routine.

!       Equatorial symmetry.

        if (eqsym) then

!          Open file and save current time.

           if (rank==0) then

              if (filestatus=='replace') then     
                 open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.dl',form='formatted', &
                 status='replace')
              else
                 open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.dl',form='formatted', &
                 status='old',position='append')
              end if

              write(1,"(A8,ES14.6)") '"Time = ',t(box,level)

           end if

!          Interpolate and save data.

           rrmax = sqrt(rmaxl(0,level)**2 + zmaxl(0,level)**2)
           delta = sqrt(drl(level)**2 + dzl(level)**2)

           NN = int(rrmax/delta + 0.5d0)
           delta = 1.d0/dble(NN)

           interpvar => grabvar

           do i=0,NN

              aux1 = (dble(i)-0.5d0)*delta
              r0 = aux1*rmaxl(0,level)
              z0 = aux1*zmaxl(0,level)
              rr0 = aux1*rrmax
              aux1 = interp(box,level,r0,z0,flag)
              call MPI_Allreduce(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

              if (rank==0) then
                 write(1,"(2ES20.8E3)") rr0,aux2
              end if

           end do

!          Leave blank space before next time and close file.

           if (rank==0) then
              write(1,*)
              close(1)
           end if

!       No equatorial symmetry.

        else

!          Open files and save current time.

           if (rank==0) then

              if (filestatus=='replace') then     
                 open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.d1l',form='formatted', &
                 status='replace')
                 open(2,file=trim(outdir)//'/'//varname//trim(filen)//'.d2l',form='formatted', &
                 status='replace')
              else
                 open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.d1l',form='formatted', &
                 status='old',position='append')
                 open(2,file=trim(outdir)//'/'//varname//trim(filen)//'.d2l',form='formatted', &
                 status='old',position='append')
              end if

              write(1,"(A8,ES14.6)") '"Time = ',t(box,level)
              write(2,"(A8,ES14.6)") '"Time = ',t(box,level)

           end if

!          Interpolate and save data.

           rrmax = sqrt(rmaxl(0,level)**2 + zmaxl(0,level)**2)
           delta = sqrt(drl(level)**2 + dzl(level)**2)

           NN = int(rrmax/delta + 0.5d0)
           delta = 1.d0/dble(NN)

           interpvar => grabvar

           do i=0,NN

              aux1 = (dble(i)-0.5d0)*delta
              r0 = aux1*rmaxl(0,level)
              z0 = aux1*zmaxl(0,level)
              rr0 = aux1*rrmax
              aux1 = interp(box,level,r0,z0,flag)
              call MPI_Allreduce(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

              if (rank==0) then
                 write(1,"(2ES20.8E3)") rr0,aux2
              end if

              aux1 = dble(i)/dble(NN)
              z0 = - z0
              aux1 = interp(box,level,r0,z0,flag)
              call MPI_Allreduce(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

              if (rank==0) then
                 write(2,"(2ES20.8E3)") rr0,aux2
              end if

           end do

!          Leave blank space before next time and close files.

           if (rank==0) then
              write(1,*)
              write(2,*)
              close(1)
              close(2)
           end if

        end if

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine save1Dvariable






