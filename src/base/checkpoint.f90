!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/checkpoint.f90,v 1.24 2021/02/17 18:22:20 malcubi Exp $

  subroutine checkpointsave

! **********************
! ***   CHECKPOINT   ***
! **********************

! This routine saves checkpoint files for a restart.

! Include modules.

  use param
  use arrays
  use procinfo
  use mpi

! Extra variables.

  implicit none

  logical firstcall

  integer i,j,k,l
  integer box,level
  integer nvars,unit
  integer, allocatable, dimension (:) :: commas

  character(1)   aa
  character(9)   filen
  character(20)  filestatus
  character(100) outdir,outtime

  character(50), allocatable, dimension (:) :: outvars

  data firstcall / .true. /

  save firstcall,nvars,outvars


! ***************************************
! ***   CREATE CHECKPOINT DIRECTORY   ***
! ***************************************

! Create a checkpoint directory and copy the parameter
! file inside it.

  write(outtime,"(F8.4)") t(0,0)
  outdir = trim(directory)//'/checkpoint'//'_t='//adjustl(trim(outtime))

  if (rank==0) then
     !call system('rm -rf '//outdir)
     call system('mkdir -p '//outdir)
     call system('cp '//trim(directory)//'/'//trim(parfile)//' '//outdir)
  end if

! I need to stop all processors here until they all reach
! the same point, as otherwise some processors might try
! to write checkpoint files before the directory is created.

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)


! ***********************************
! ***   IS THIS THE FIRST CALL?   ***
! ***********************************

! On first call figure out what needs output.

  if (firstcall) then

     firstcall = .false.

!    Find out length of string "checkvars".

     l = len_trim(checkvars)

!    Find out how many variables need output.

     nvars = 1

     do i=1,l
        aa = checkvars(i:i)
        if (aa==",") then
           nvars = nvars + 1
        end if
     end do

!    Identify comma positions.

     allocate(commas(0:nvars))

     j = 0
     commas(0) = 0
     commas(nvars) = l+1

     do i=1,l
        aa = checkvars(i:i)
        if (aa==",") then
           j = j+1
           commas(j)=i
        end if
     end do

!    Now read variable names, and eliminate spaces.

     allocate(outvars(1:nvars))

     do i=1,nvars
        outvars(i) = checkvars(commas(i-1)+1:commas(i)-1)
        outvars(i) = trim(adjustl(outvars(i)))
        !print *, 'Array: ',outvars(i)
     end do

  end if


! ***********************************
! ***   SAVE CRUCIAL PARAMETERS   ***
! ***********************************

! Some parameters are crucial for a restart:
! symmetries, grid structure, matter content,
! etcetera, and they need to be the same.
!
! It is not enough to re-read the old parameter
! file since many may have had their default
! values and that won't show up in the parfile.
! So I save them here to the file "checkparam".

  if (rank==0) then

!    Open file.

     open(1,file=trim(outdir)//'/checkparam',form='formatted',status='replace')

!    Save number of processors used in run.

     write(1,*) 'nproctot = ',nproctot

!    Save symmetry and angular momentum flags.

     write(1,*) 'eqsym = ',eqsym
     write(1,*) 'angmom = ',angmom

!    Save grid spacings.

     write(1,*) 'dr0 = ',dr0
     write(1,*) 'dz0 = ',dz0

!    Save grid size.

     write(1,*) 'Nrtotal = ',Nrtotal
     write(1,*) 'Nztotal = ',Nztotal

!    Save number of boxes, sizes, and number of levels.

     write(1,*) 'Nb = ',Nb

     write(1,*) 'Nl0 = ',Nl0
     write(1,*) 'Nl1 = ',Nl1
     write(1,*) 'Nl2 = ',Nl2
     write(1,*) 'Nl3 = ',Nl3

     write(1,*) 'Nrbox1 = ',Nrbox1
     write(1,*) 'Nzbox1 = ',Nzbox1
     write(1,*) 'Nrbox2 = ',Nrbox2
     write(1,*) 'Nzbox2 = ',Nzbox2
     write(1,*) 'Nrbox3 = ',Nrbox3
     write(1,*) 'Nzbox3 = ',Nzbox3

!    Save positions of refinement boxes.

     write(1,*) 'rbox1 = ',rbox1
     write(1,*) 'zbox1 = ',zbox1
     write(1,*) 'rbox2 = ',rbox2
     write(1,*) 'zbox2 = ',zbox2
     write(1,*) 'rbox3 = ',rbox3
     write(1,*) 'zbox3 = ',zbox3

!    Save matter content.

     write(1,*) 'mattertype = ',trim(mattertype)

!    Close file.

     close(1)

  end if


! *********************
! ***   SAVE DATA   ***
! *********************

! For a checkpoint we always replace the old files.

  filestatus = 'replace'

! Save name of parfile and list of variables
! that have output.

  if (rank==0) then

     open(1,file=trim(outdir)//'/checklist',form='formatted',status='replace')

     write(1,*) trim(parfile)

     do i=1,nvars
        write(1,*) trim(outvars(i))
     end do

     close(1)

  end if

! Loop over boxes and grid levels.

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Append box number, grid level and processor number to file name.

        if (size==1) then
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

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Save all variables that require output. Notice that
!       we call the subroutine "save2Dvariable" that is
!       found in the file "save2Ddata.f90".

        unit = rank + 1

        do i=1,nvars

!          Save 2D data.

           call grabarray(trim(outvars(i)))
           call save2Dvariable(trim(outvars(i)),outdir,box,level,'filepreproc',filestatus)

!          Don't save boundary data for coordinates (it does not exist).

           if ((outvars(i)=='r').or.(outvars(i)=='z').or.(outvars(i)=='rr')) then
              cycle
           end if

!          Lower r boundary.

           open(unit,file=trim(outdir)//'/'//trim(outvars(i))//'_bound_rL'//trim(filen),form='formatted',status='replace')

           write(unit,'(3ES23.15)') t(box,level),t1(box,level),t2(box,level)
           write(unit,*)

           do j=0,ghost-1
              do k=1-ghost,Nzmaxl(box)
                 write(unit,'(3ES23.15)') grabvar_bound_rL(j,k,0),grabvar_bound_rL(j,k,1),grabvar_bound_rL(j,k,2)
              end do
           end do

           close(unit)

!          Upper r boundary.

           open(unit,file=trim(outdir)//'/'//trim(outvars(i))//'_bound_rR'//trim(filen),form='formatted',status='replace')

           write(unit,'(3ES23.15)') t(box,level),t1(box,level),t2(box,level)
           write(unit,*)

           do j=0,ghost-1
              do k=1-ghost,Nzmaxl(box)
                 write(unit,'(3ES23.15)') grabvar_bound_rR(j,k,0),grabvar_bound_rR(j,k,1),grabvar_bound_rR(j,k,2)
              end do
           end do

           close(unit)

!          Lower z boundary.

           open(unit,file=trim(outdir)//'/'//trim(outvars(i))//'_bound_zL'//trim(filen),form='formatted',status='replace')

           write(unit,'(3ES23.15)') t(box,level),t1(box,level),t2(box,level)
           write(unit,*)

           do j=1-ghost,Nrmaxl(box)
              do k=0,ghost-1
                 write(unit,'(3ES23.15)') grabvar_bound_zL(j,k,0),grabvar_bound_zL(j,k,1),grabvar_bound_zL(j,k,2)
              end do
           end do

           close(unit)

!          Upper z boundary.

           open(unit,file=trim(outdir)//'/'//trim(outvars(i))//'_bound_zR'//trim(filen),form='formatted',status='replace')

           write(unit,'(3ES23.15)') t(box,level),t1(box,level),t2(box,level)
           write(unit,*)

           do j=1-ghost,Nrmaxl(box)
              do k=0,ghost-1
                 write(unit,'(3ES23.15)') grabvar_bound_zR(j,k,0),grabvar_bound_zR(j,k,1),grabvar_bound_zR(j,k,2)
              end do
           end do

           close(unit)

        end do

     end do
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine checkpointsave









  subroutine checkpointrestart

! *****************************
! ***   CHECKPOINT RESTART  ***
! *****************************

! This routine restarts from checkpoint files.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical flag
  logical eqsym_new,angmom_new

  integer i,j,k,nvars
  integer box,level
  integer nproctot_new,Nr_new,Nz_new
  integer Nb_new,Nl0_new,Nl1_new,Nl2_new,Nl3_new
  integer Nrbox1_new,Nzbox1_new,Nrbox2_new,Nzbox2_new,Nrbox3_new,Nzbox3_new
  integer unit
  integer error

  real(8) dr_new,dz_new
  real(8) rbox1_new,rbox2_new,rbox3_new
  real(8) zbox1_new,zbox2_new,zbox3_new

  character(9)    filen
  character(1000) matter_new
  character(1000) checkparam

  character(50), dimension (0:1000) :: outvars


! *****************
! ***   START   ***
! *****************

! Message to screen.

  if (rank==0) then
     print *, 'Starting from checkpoint!'
     print *
  end if

! Check if the checkpoint directory exists.

  inquire(FILE='./'//trim(checkpointfile),EXIST=flag)

  if (.not.flag) then
     if (rank==0) then
        print *, 'Checkpoint directory "',trim(checkpointfile),'" does not exist.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

! Warning.

  if (outparallel/="fileperproc") then
     if (rank==0) then
        print *, 'Checkpoint restart requires fileperproc type output.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if


! ************************
! ***   SANITY CHECK   ***
! ************************

! We need to make sure that the values of certain
! crucial parameters are the same as in the old run,
! or the new run will not work.
!
! Notice that many parameters can in fact change and it
! won't matter, but parameters related to symmetries,
! grid structure, matter content, etcetera, need to be
! the same.
!
! The first step it to save the values that those parameters
! have after reading the new parameter file.

  nproctot_new = nproctot

  eqsym_new  = eqsym
  angmom_new = angmom

  dr_new = dr0
  dz_new = dz0

  Nr_new = Nrtotal
  Nz_new = Nztotal

  Nb_new = Nb

  Nl0_new = Nl0
  Nl1_new = Nl1
  Nl2_new = Nl2
  Nl3_new = Nl3

  Nrbox1_new = Nrbox1
  Nzbox1_new = Nzbox1
  Nrbox2_new = Nrbox2
  Nzbox2_new = Nzbox2
  Nrbox3_new = Nrbox3
  Nzbox3_new = Nzbox3

  matter_new = mattertype

! Now parse the file "checkparam".

  checkparam = trim(checkpointfile)//'/checkparam'

  call parse(rank,trim(checkparam))

! Check number of processors.

  if (nproctot_new/=nproctot) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint we must use the same number of processors.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

! Check equatorial symmetry.

  if (eqsym_new.neqv.eqsym) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint the equatorial symmetry flag must be the same.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

! Check angular momentum flag.

  if (angmom_new.neqv.angmom) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint the angular momentum flag must be the same.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

! Check grid spacing.

  if ((dr0/=dr_new).or.(dz0/=dz_new)) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint the grid spacing must be the same.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

! Check grid size.

  if ((Nrtotal/=Nr_new).or.(Nztotal/=Nz_new)) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint the grid size must be the same.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

! Check number of boxes, sizes, and number of levels.

  if (Nb/=Nb_new) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint the number of refinement boxes must be the same.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

  if ((Nl0/=Nl0_new).or.(Nl1/=Nl1_new).or.(Nl2/=Nl2_new).or.(Nl3/=Nl3_new)) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint the number of levels on refinement boxes must be the same.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

  if ((Nrbox1/=Nrbox1_new).or.(Nzbox1/=Nzbox1_new).or. &
      (Nrbox2/=Nrbox2_new).or.(Nzbox2/=Nzbox2_new).or. &
      (Nrbox3/=Nrbox3_new).or.(Nzbox3/=Nzbox3_new)) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint the size of the refinement boxes must be the same.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

! Check matter content.

  if (matter_new/=mattertype) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint the matter type must be the same.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if


! ***************************
! ***   READ DATA FILES   ***
! ***************************

  unit = rank + 1

! Read names of variables in the checkpoint file.

  open(1,file=trim(checkpointfile)//'/checklist',form='formatted',status='old')

  i = 0

  do
     read(1,*,end=100) outvars(i)
     i=i+1
  end do

  100 nvars = i-1

  close(1)

! Loop over boxes and levels.

  do box=0,Nb
     do level=min(1,box),Nl(box)

!       Append box number and grid level to file name.

        if (size==1) then
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

!       Point to current grid.

        call currentgrid(box,level,grid(box,level))

!       Read arrays.

        do i=1,nvars

!          Read 2D data.

           call grabarray(trim(outvars(i)))
           call read2Ddata(outvars(i),checkpointfile,box,level,'fileperproc')

!          Don't read boundary data for coordinates.

           if ((outvars(i)=='r').or.(outvars(i)=='z').or.(outvars(i)=='rr')) then
              cycle
           end if

!          Lower r boundary.

           open(unit,file=trim(checkpointfile)//'/'//trim(outvars(i))//'_bound_rL'//trim(filen),form='formatted',status='old')

           read(unit,*) t(box,level),t1(box,level),t2(box,level)
           read(unit,*)

           do j=0,ghost-1
              do k=1-ghost,Nzmaxl(box)
                 read(unit,*) grabvar_bound_rL(j,k,0),grabvar_bound_rL(j,k,1),grabvar_bound_rL(j,k,2)
              end do
           end do
           
           close(unit)

!          Upper r boundary.

           open(unit,file=trim(checkpointfile)//'/'//trim(outvars(i))//'_bound_rR'//trim(filen),form='formatted',status='old')

           read(unit,*) t(box,level),t1(box,level),t2(box,level)
           read(unit,*)

           do j=0,ghost-1
              do k=1-ghost,Nzmaxl(box)
                 read(unit,*) grabvar_bound_rR(j,k,0),grabvar_bound_rR(j,k,1),grabvar_bound_rR(j,k,2)
              end do
           end do

           close(unit)

!          Lower z boundary.

           open(unit,file=trim(checkpointfile)//'/'//trim(outvars(i))//'_bound_zL'//trim(filen),form='formatted',status='old')

           read(unit,*) t(box,level),t1(box,level),t2(box,level)
           read(unit,*)

           do j=1-ghost,Nrmaxl(box)
              do k=0,ghost-1
                 read(unit,*) grabvar_bound_zL(j,k,0),grabvar_bound_zL(j,k,1),grabvar_bound_zL(j,k,2)
              end do
           end do

           close(unit)

!          Upper z boundary.

           open(unit,file=trim(checkpointfile)//'/'//trim(outvars(i))//'_bound_zR'//trim(filen),form='formatted',status='old')

           read(unit,*) t(box,level),t1(box,level),t2(box,level)
           read(unit,*)

           do j=1-ghost,Nrmaxl(box)
              do k=0,ghost-1
                 read(unit,*) grabvar_bound_zR(j,k,0),grabvar_bound_zR(j,k,1),grabvar_bound_zR(j,k,2)
              end do
           end do

           close(unit)

        end do

     end do
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine checkpointrestart

