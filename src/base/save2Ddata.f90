!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/save2Ddata.f90,v 1.29 2022/10/24 21:58:46 malcubi Exp $

  subroutine save2Ddata

! ********************************
! ***   SAVE 2D DATA TO FILE   ***
! ********************************

! This routine saves "2D" data to files.

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

!    Find out length of string "outvars2D".

     l = len_trim(outvars2D)

!    Find out how many variables need output.

     nvars2D = 1

     do i=1,l
        aa = outvars2D(i:i)
        if (aa==",") then
           nvars2D = nvars2D + 1
        end if
     end do

!    Identify comma positions.

     allocate(commas(0:nvars2D))

     j = 0
     commas(0) = 0
     commas(nvars2D) = l+1

     do i=1,l
        aa = outvars2D(i:i)
        if (aa==",") then
           j = j+1
           commas(j)=i
        end if
     end do

!    Now read variable names, and eliminate spaces.

     allocate(outvars2Darray(1:nvars2D))

     do i=1,nvars2D
        outvars2Darray(i) = outvars2D(commas(i-1)+1:commas(i)-1)
        outvars2Darray(i) = trim(adjustl(outvars2Darray(i)))
     end do

!    Check if any name is repeated.

     if (rank==0) then
        do i=1,nvars2D
           do j=1,nvars2D
              if ((i.ne.j).and.(trim(outvars2Darray(i))==trim(outvars2Darray(j)))) then
                 print *, 'Error in parfile, array name repeated in outvars2D: ',outvars2Darray(i)
                 print *, 'Aborting! (subroutine save2Ddata)'
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

        do i=1,nvars2D
           call grabarray(trim(outvars2Darray(i)))
           call save2Dvariable(trim(outvars2Darray(i)),directory,box,level,outparallel,filestatus)
        end do

     end do
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine save2Ddata









  subroutine save2Dvariable(varname,outdir,box,level,type,filestatus)

! This subroutine outputs 2D files.
!
! At the moment I only do formatted output, and for each
! grid point I output a triplet of the form: (r,z,var).
!
! This is clearly not very efficient since:
!
!     1) For a structured grid we could deduce
!        the values of (r,z), or alternatively
!        just output one file for each, thus
!        making the output files smaller.
!
!     2) Formatted files are large, we could use
!        unformatted FORTRAN output, but then we
!        would need to make sure that the plotting
!        packages understand this.
!
! This all means that it might be a good idea to improve
! the output routine later.

! Load modules.

  use mpi
  use param
  use arrays
  use procinfo

! Declare variables.

  implicit none

  integer i,j,m,n,p
  integer box,level
  integer i0,j0,imax,jmax,Naux
  integer unit
  integer status(MPI_STATUS_SIZE)

  real(8) tiny

  real(8), dimension (1-ghost:Nrmax) :: varr,auxr,auxz

  character(len=9)  filen
  character(len=10) form
  character(len=*)  varname,outdir,type,filestatus


! *******************************
! ***   FILE NAME EXTENSION   ***
! *******************************

! Append box number, grid level and processor number to file name.

  if ((size==1).or.(type=="singlefile")) then
     write(filen,'(a2,i1,a1,i1)') '_b',box,'l',level
  else
     if (size<10) then
        write(filen,'(a,i1,a,i1,a,i1)') '_b',box,'l',level,'p',rank
     else if (size<100) then
        write(filen,'(a,i1,a,i1,a,i2)') '_b',box,'l',level,'p',rank
     else
        write(filen,'(a,i1,a,i1,a,i3)') '_b',box,'l',level,'p',rank
     end if
  end if


! *************************
! ***   OUTPUT FORMAT   ***
! *************************

! For checkpoint we need output with many significant figures,
! while for plotting we don't.  The easiest way to ensure this
! is saving 15 significant figures whenever the file is replaced,
! and saving fewer when it is only appended to.

  if (filestatus=='replace') then
     form = "(3ES23.15)"
  else
     form = "(3ES16.8)"
  end if


! *************************
! ***   SAVE VARIABLE   ***
! *************************

  tiny = 1.d-30

  if (size==1) then

!    *****************************
!    ***   ONE PROCESSOR RUN   ***
!    *****************************

!    Open file.

     if (filestatus=='replace') then
        open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.2D',form='formatted', &
             status=filestatus)
     else
        open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.2D',form='formatted', &
             status=filestatus,position='append')
     end if

!    Write current time.

     write(1,"(A8,ES23.15)") '#Time = ',t(box,level)

!    Write data.

     do j=1-ghost,Nz
        do i=1-ghost,Nr
     !do j=1,Nz
        !do i=1,Nr
           write(1,form) r(i,j),z(i,j),grabvar(i,j)+tiny
        end do
        write (1,*)
     end do

!    Leave blank line before next time.

     write(1,*)

!    Close file.

     close(1)


  else if (type=="singlefile") then

!    *************************************
!    ***   PARALLEL RUN - SINGLE FILE  ***
!    *************************************

!    If we have more than one processor and we want single
!    file output, then the logic is more complicated.
!    Basically, we output one full r line at a time and
!    slowly move up in z. The problem is that each processor
!    has many z lines, so at a given value of z we move from
!    one proc to the next along r. We then increment z and do
!    it again until we have covered all the values of z in the
!    given level of procs. Then we move up a level of procs 
!    and start again until we cover the whole grid.

!    Only processor 0 does output.

     if (rank==0) then

!       Open file.

        if (filestatus=='replace') then
           open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.2D',form='formatted', &
                status=filestatus)
        else
           open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.2D',form='formatted', &
                status=filestatus,position='append')
        end if

!       Write current time.

        write(1,"(A8,ES23.15)") '#Time = ',t(box,level)

!       Loop over processors in z direction.

        do m=0,nprocz-1

!          Loop over z.

           if (m==0) then
              j0 = 1-ghost
           else
              j0 = 1
           end if

           if (m==nprocz-1) then
              jmax = Nzl(box,m*nprocr)
           else
              jmax = Nzl(box,m*nprocr)-ghost
           end if

           do j=j0,jmax

!             Loop over processors in r direction.

              do n=0,nprocr-1

!                Find identity of processor that needs to output a line.

                 p = n + m*nprocr

!                Receive data and store it in (auxr,auxz,varr).

                 Naux = Nrmax + ghost

                 if (p==0) then
                    varr = grabvar(:,j)
                    auxr = r(:,j)
                    auxz = z(:,j)
                 else
                    call MPI_RECV(varr,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
                    call MPI_RECV(auxr,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
                    call MPI_RECV(auxz,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
                 end if

!                Save one line.

                 if (n==0) then
                    i0 = 1-ghost
                 else
                    i0 = 1
                 end if

                 if (n==nprocr-1) then
                    imax = Nrl(box,p)
                 else
                    imax = Nrl(box,p)-ghost
                 end if

                 do i=i0,imax
                    write(1,form) auxr(i),auxz(i),varr(i)+tiny
                 end do

              end do

!             Leave blank space.

              write(1,*)

           end do
        end do

!       Leave blank line before next time.

        write(1,*)

!       Close file.

        close(1)

!    Other processors send data to processor 0.

     else

        Naux = Nrmax + ghost

        if ((.not.eqsym).and.(rank<nprocr)) then
           j0 = 1-ghost
        else
           j0 = 1
        end if

        if (rank>=size-nprocr) then
           jmax = Nzl(box,rank)
        else
           jmax = Nzl(box,rank)-ghost
        end if

        do j=j0,jmax
           varr = grabvar(:,j)
           auxr = r(:,j)
           auxz = z(:,j)
           call MPI_SEND(varr,Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(auxr,Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(auxz,Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end do

     end if


  else

!    ********************************************
!    ***   PARALLEL RUN - ONE FILE PER PROC   ***
!    ********************************************

!    Each processor does its own output.

!    Open file.

     unit = rank+1

     if (filestatus=='replace') then
        open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'.2D',form='formatted', &
             status=filestatus)
     else
        open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'.2D',form='formatted', &
             status=filestatus,position='append')
     end if

!    Write current time.

     write(unit,"(A8,ES23.15)") '#Time = ',t(box,level)

!    Write data.

     do j=1-ghost,Nz
        do i=1-ghost,Nr
           write(unit,form) r(i,j),z(i,j),grabvar(i,j)+tiny
        end do
        write (unit,*)
     end do

!    Leave blank line before next time.

     write(unit,*)

!    Close file.

     close(unit)

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine save2Dvariable


