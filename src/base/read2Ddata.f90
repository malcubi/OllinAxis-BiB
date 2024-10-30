!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/read2Ddata.f90,v 1.8 2021/02/17 18:22:42 malcubi Exp $

  subroutine read2Ddata(varname,outdir,box,level,type)

! *************************
! ***   READ 2D FILES   ***
! *************************

! Ths subroutine reads 2D files.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical contains

  integer i,j             ! Counters.
  integer box,level       ! Box and level counters.
  integer unit            ! Output unit.

  real(8) rdata,zdata     ! Position data.
  real(8) tdata           ! Time data.

  real(8), dimension (1-ghost:Nrtotal,1-ghost:Nztotal) :: varg   ! Global array.

  character(1) first
  character(9) filen
  character(50) varname,outdir,type,line


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
! ***   READ VARIABLE   ***
! *************************

  if (size==1) then

!    *****************************
!    ***   ONE PROCESSOR RUN   ***
!    *****************************

!    Open file.

     open(1,file=trim(outdir)//'/'//trim(varname)//trim(filen)//'.2D',form='formatted',status='old')

!    Read first line and check if it is a comment.
!    If it is, then read the current time from it.

     read(1,'(A)') line

     line = adjustl(line)
     first = line(1:1)

     if (first=="#") then
        if (contains(line,"#Time")) then
           i = index(line,'=')
           line = trim(adjustl(line(i+1:len_trim(line))))
           read(line,*) tdata
           t = tdata
        end if
     else
        rewind(1)
     end if

!    Read data.

     do j=1-ghost,Nz
        do i=1-ghost,Nr
           read(1,*) rdata,zdata,grabvar(i,j)
        end do
     end do

!    Close file.

     close(1)


  else if (type=="singlefile") then

!    *************************************
!    ***   PARALLEL RUN - SINGLE FILE  ***
!    *************************************

!    Processor 0 reads the data and distributes it
!    to the other processors.

     if (rank==0) then

!       Open file.

        open(1,file=trim(outdir)//'/'//trim(varname)//trim(filen)//'.2D',form='formatted',status='old')

!       Read first line and check if it is a comment.
!       If it is then read the current time from it.

        read(1,'(A)') line

        line = adjustl(line)
        first = line(1:1)

        if (first=="#") then
           if (contains(line,"#Time")) then
              i = index(line,'=')
              line = trim(adjustl(line(i+1:len_trim(line))))
              read(line,*) tdata
              t = tdata
           end if
        else
           tdata = t(0,0)
           rewind(1)
        end if

!       Read data.

        do j=1-ghost,Nztotal
           do i=1-ghost,Nrtotal
              read(1,*) rdata,zdata,varg(i,j)
           end do
        end do

!       Close file.

        close(1)

     end if

!    Send time information to all other processors.

     call MPI_BCAST(tdata,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     !print *, rank,tdata

     t = tdata

!    Now distribute the data.

     call distribute(grabvar,varg)


  else

!    ********************************************
!    ***   PARALLEL RUN - ONE FILE PER PROC   ***
!    ********************************************

!    Each processor reads its own data.

!    Open file.

     unit = rank+1

     open(unit,file=trim(outdir)//'/'//trim(varname)//trim(filen)//'.2D',form='formatted',status='old')

!    Read first line and check if it is a comment.
!    If it is then read the current time from it.

     read(unit,'(A)') line

     line = adjustl(line)
     first = line(1:1)

     if (first=="#") then
        if (contains(line,"#Time")) then
           i = index(line,'=')
           line = trim(adjustl(line(i+1:len_trim(line))))
           read(line,*) tdata
           t = tdata
        end if
     else
        rewind(unit)
     end if

!    Read data.

     do j=1-ghost,Nz
        do i=1-ghost,Nr
           read(unit,*) rdata,zdata,grabvar(i,j)
        end do
     end do

!    Close file.

     close(unit)

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine read2Ddata


