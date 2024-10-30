!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/save0Ddata.f90,v 1.27 2021/03/28 00:05:16 malcubi Exp $

  subroutine save0Ddata

! ********************************
! ***   SAVE 0D DATA TO FILE   ***
! ********************************

! This routine saves "0D" data to files. By 0D data I mean reduced
! quantities as functions of time (minimum, maximum, norms, etc).

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

!    Find out length of string "outvars0D".

     l = len_trim(outvars0D)

!    Find out how many variables need output.

     nvars0D = 1

     do i=1,l
        aa = outvars0D(i:i)
        if (aa==",") then
           nvars0D = nvars0D + 1
        end if
     end do

!    Identify comma positions.

     allocate(commas(0:nvars0D))

     j = 0
     commas(0) = 0
     commas(nvars0D) = l+1

     do i=1,l
        aa = outvars0D(i:i)
        if (aa==",") then
           j = j+1
           commas(j)=i
        end if
     end do

!    Now read variable names, and eliminate spaces.

     allocate(outvars0Darray(1:nvars0D))

     do i=1,nvars0D
        outvars0Darray(i) = outvars0D(commas(i-1)+1:commas(i)-1)
        outvars0Darray(i) = trim(adjustl(outvars0Darray(i)))
     end do

!    Check if any name is repeated.

     if (rank==0) then
        do i=1,nvars0D
           do j=1,nvars0D
              if ((i.ne.j).and.(trim(outvars0Darray(i))==trim(outvars0Darray(j)))) then
                 print *, 'Error in parfile, array name repeated in outvars0D: ',outvars0Darray(i)
                 print *, 'Aborting! (subroutine save0Ddata)'
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

! Loop over boxes and grid levels. I loop over levels
! from fine to coarse. The reason for this is that I want
!to catch any possible NaN's on the finer grids first.

  do box=0,Nb
     do level=Nl(box),min(1,box),-1

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Save all variables that require output.

        do i=1,nvars0D
           call grabarray(trim(outvars0Darray(i)))
           call save0Dvariable(trim(outvars0Darray(i)),directory,box,level,s(0,0),t(0,0),filestatus)
        end do

     end do
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine save0Ddata








  subroutine save0Dvariable(varname,outdir,box,level,step,tt,filestatus)

! This subroutine calculates minimum, maximum, norm 1, norm 2,
! and total variation of a variable and outputs them.
!
! The quantities saved are:
!
! vmax:    Maximum value of the variable over the grid.
! vmin:    Minimum value of the variable over the grid.
!
! nm1:     Maximum absolute value.
! nm2:     Root mean square (rms).
!
! tvar:    Total variation.

! Load modules.

  use mpi
  use param
  use arrays
  use procinfo

! Declare variables.

  implicit none

  logical interpflag

  integer i,j
  integer step,box,level
  integer unit
  integer foundnan
  integer naux

  real(8) vmax,gmax,vmin,gmin,v0
  real(8) nm1,gnm1,nm2,gnm2
  real(8) NN,NNtot
  real(8) tt,interp
  real(8) aux

  character(len=1) comment
  character(len=9) filen
  character(len=*) varname,outdir,filestatus


! ************************
! ***   COMMENT TYPE   ***
! ************************

! Decide on comment type.

  if (commenttype=='xgraph') then
     comment = '"'
  else
     comment = '#'
  end if


! ****************************
! ***   FIND LOCAL NORMS   ***
! ****************************

! Set NaN catching flag to false.

  foundnan = 0

! Find minimum, maximum and norms.

  vmax = -1.d10
  vmin = +1.d10

  nm1 = 0.d0
  nm2 = 0.d0

  NN  = 0.d0

  do j=1,Nzl(box,rank)-ghost
     do i=1,Nrl(box,rank)-ghost

        aux = grabvar(i,j)

        if (isnan(aux)) then
           print *
           write(*,"(A,A,A,4I5,E18.8)") ' Found NaN in value of ',varname,' at (proc,box,level,step,time):',rank,box,level,step,tt
           write(*,"(A,2I5,2E18.8)") ' at point (i,j,r,z):',i,j,r(i,j),z(i,j)
           print *, 'Aborting! (subroutine save0Ddata.f90)'
           print *
           foundnan = 1
           goto 100
        end if

        NN = NN + 1.d0

        vmin = min(vmin,aux)
        vmax = max(vmax,aux)

        nm1 = max(nm1,abs(aux))
        nm2 = nm2 + aux**2

     end do
  end do

  100 continue

! If we found a NaN die.

  naux = foundnan; foundnan = 0
  call MPI_Allreduce(naux,foundnan,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierr)

  if (foundnan>0) then
     call die
  end if


! **************************************************
! ***   REDUCE QUANTITIES AMONG ALL PROCESSORS   ***
! **************************************************

! Global maximum.

  call MPI_Allreduce(vmax,gmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

! Global minimum.

  call MPI_Allreduce(vmin,gmin,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)

! Global nm1.

  call MPI_Allreduce(nm1,gnm1,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

! Global nm2.

  call MPI_Allreduce(NN,NNtot,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(nm2,gnm2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  gnm2 = sqrt(gnm2/NNtot)


! *****************************
! ***   SAVE GLOBAL NORMS   ***
! *****************************

  unit = 1

! Only processor 0 does output.

  if (rank==0) then

!    Append box number and grid level to file name.

     write(filen,'(a2,i1,a1,i1)') '_b',box,'l',level

!    Save maximum value.

     if (filestatus=='replace') then
        open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_max.tl',form='formatted', &
        status='replace')
        write(unit,*) comment//varname//'_max.tl'
     else
        open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_max.tl',form='formatted', &
        status='old',position='append')
     end if

     write(unit,"(2ES20.8E3)") t(box,level),gmax

     close(unit)

!    Save minimum value.

     if (filestatus=='replace') then
        open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_min.tl',form='formatted', &
        status='replace')
        write(unit,*) comment//varname//'_min.tl'
     else
        open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_min.tl',form='formatted', &
        status='old',position='append')
     end if

     write(unit,"(2ES20.8E3)") t(box,level),gmin

     close(unit)

!    Save norm 1.

     if (filestatus=='replace') then
        open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_nm1.tl',form='formatted', &
        status='replace')
        write(unit,*) comment//varname//'_nm1.tl'
     else
        open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_nm1.tl',form='formatted', &
        status='old',position='append')
     end if

     write(unit,"(2ES20.8E3)") t(box,level),gnm1

     close(unit)

!    Save norm 2.

     if (filestatus=='replace') then
        open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_nm2.tl',form='formatted', &
        status='replace')
        write(unit,*) comment//varname//'_nm2.tl'
     else
        open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_nm2.tl',form='formatted', &
        status='old',position='append')
     end if

     write(unit,"(2ES20.8E3)") t(box,level),gnm2

     close(unit)

  end if


! ********************************
! ***   SAVE VALUE AT ORIGIN   ***
! ********************************

! We only save origin data for box 0.

  if (box==0) then

     interpvar => grabvar

     aux = interp(box,level,0.d0,0.d0,interpflag)
     call MPI_ALLREDUCE(aux,v0,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!    Only processor 0 does output.

     unit = 1

     if (rank==0) then

        if (filestatus=='replace') then
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_origin.tl',form='formatted', &
           status='replace')
           write(unit,*) comment//varname//'_max.tl'
        else
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_origin.tl',form='formatted', &
           status='old',position='append')
        end if

        write(unit,"(2ES20.8E3)") t(box,level),v0

        close(unit)

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine save0Dvariable


