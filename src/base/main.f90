!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/main.f90,v 1.32 2021/03/08 01:00:27 malcubi Exp $

! *****************************
! ***   PROGRAM OLLINAXIS   ***
! *****************************

  program main

! Load MPI module.  If we compile without MPI this will
! just load a dummy module.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical flag               ! Auxiliary flag.

  integer j,box,proc         ! Counters.

  real(8) aux                ! Auxiliary.

  character(60) message,string1,string2

! The variables introduced above are:
!
! parfile       Name of parameter file.
!
! message       Function to format messages to screen (see functions.f90)
! string*       Auxiliary strings.


! **************************
! ***   INITIALIZE MPI   ***
! **************************

! When compiling without MPI this call does nothing at all.

  call MPI_INIT(ierr)


! ********************************************************
! ***   FIND OUT NUMBER OF PROCESSORS AND LOCAL RANK   ***
! ********************************************************

! When compiling without MPI these calls just return size=1 and rank=0.

  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  nproctot = size


! ***************************
! ***   INITIAL MESSAGE   ***
! ***************************

  if (rank==0) then

     print *
     print *, '******************************************'
     print *, '******************************************'
     print *, '***                                    ***'
     print *, '***             OLLINAXIS              ***'
     print *, '***                                    ***'
     print *, '***   Evolving Einstein''s equations    ***'
     print *, '***         in axial symmetry          ***'
     print *, '***                                    ***'    
     print *, '******************************************'
     print *, '******************************************'
     print *

  end if


! *********************************
! ***   CALL PARAMETER PARSER   ***
! *********************************

! Get name of parameter file.

  call getarg(1,parfile)

  if (parfile==" ") then
     if (rank==0) then
        print *, 'Missing parfile name.'
        print *, 'Aborting! (subroutine main)'
        print *
     end if
     call die
  end if

! Call parser.

  call parse(rank,parfile)


! *************************
! ***   SANITY CHECKS   ***
! *************************

! Various sanity checks.  Might add more later.

! Ninfo must be positive.

  if (Ninfo<=0) then
     if (rank==0) then
        print *, 'Ninfo must be positive'
        print *, 'Aborting! (subroutine main)'
        print *
     end if
     call die
  end if

! Noutput must be non-negative.  Notice that 0 is allowed,
! but this is a special case for testing purposes only.
! When the corresponding parameter Noutput is equal to 0,
! we do output al ALL time steps for each box and level,
! and not only at the coarse time steps.

  if (Noutput0D<0) then
     if (rank==0) then
        print *, 'Noutput0D must be positive'
        print *, 'Aborting! (subroutine main)'
        print *
     end if
     call die
  end if

  if (Noutput1D<0) then
     if (rank==0) then
        print *, 'Noutput1D must be positive'
        print *, 'Aborting! (subroutine main)'
        print *
     end if
     call die
  end if

  if (Noutput2D<0) then
     if (rank==0) then
        print *, 'Noutput2D must be positive'
        print *, 'Aborting! (subroutine main)'
        print *
     end if
     call die
  end if

! Ncheckpoint must be positive.

  if (Ncheckpoint<=0) then
     if (rank==0) then
        print *, 'Ncheckpoint must be positive'
        print *, 'Aborting! (subroutine main)'
        print *
     end if
     call die
  end if


! **************************
! ***   OTHER MESSAGES   ***
! **************************

  if (rank==0) then

     print *, '******************************************'
     print *, '******************************************'
     print *, '***                                    ***'


!    Formulation message.

     string1 = "formulation:"
     string2 = trim(adjustl(formulation))
     print *, message(string1,string2)

!    Slicing message.

     print *, '***                                    ***'

     string1 = "slicing:"
     string2 = trim(adjustl(slicing))
     print *, message(string1,string2)

!    Shift message.

     string1 = "shift:"
     string2 = trim(adjustl(shift))
     print *, message(string1,string2)

!    Type of matter.

     print *, '***                                    ***'

     string1 = "matter type:"
     string2 = trim(adjustl(mattertype))
     print *, message(string1,string2)

!    Initial data message.

     print *, '***                                    ***'

     string1 = "initial data:"
     string2 = trim(adjustl(idata))

     print *, message(string1,string2)

!    Evolution method message.

     print *, '***                                    ***'

     string1 = "evolution method:"
     string2 = trim(adjustl(integrator))
     print *, message(string1,string2)

     string1 = "order of integration:"
     string2 = trim(adjustl(order))
     print *, message(string1,string2)

!    Boundary type message.

     string1 = "boundary condition:"
     string2 = trim(adjustl(boundtype))
     print *, message(string1,string2)

!    Output directory message.

     print *, '***                                    ***'

     string1 = "output directory:"
     string2 = " "
     print *, message(string1,string2)

     string1 = trim(adjustl(directory))
     string2 = " "
     print *, message(string1,string2)

!    End message.

     print *, '***                                    ***'   
     print *, '******************************************'
     print *, '******************************************'
     print *

  end if


! ***********************************
! ***   CREATE OUTPUT DIRECTORY   ***
! ***********************************

! Create output directory and copy the parameter
! file inside it.  But in case it already existed,
! remove it first.

  if (rank==0) then
     call system('rm -rf '//trim(directory))
     call system('mkdir -p '//trim(directory))
     call system('cp '//trim(parfile)//' '//trim(directory))
  end if


! ***********************
! ***   GHOST ZONES   ***
! ***********************

! For second order I use 2 ghost zones and for
! fourth order 3 ghost zones.  This might seem
! one more ghost zone than needed, but one must
! remember that for dissipation and one-sided
! derivatives we need an extra ghost zone.
!
! I allow the number of ghost zones to be set
! in the parameter file (though normally it
! shouldn't) for testing, but in that case it
! should be at least 2 for second order and 3
! for fourth order.

  if (ghost==0) then

!    If number of ghost zones was not set in parameter file.

     if (order=="two") then
        ghost = 2
     else if (order=="four") then
        ghost = 3
     end if

!    If number of ghost zones was set in parameter file check
!    that it is not smaller than it should be.

  else

     if (rank==0) then
        print *, 'Number of ghost zones fixed in parameter file'
     end if

     if ((order=="two").and.(ghost<2)) then
        ghost = 2
        if (rank==0) then
           print *, 'Number of ghost zones too small, changing it ...'
        end if
     else if ((order=="four").and.(ghost<3)) then
        ghost = 3
        if (rank==0) then
           print *, 'Number of ghost zones too small, changing it ...'
        end if
     end if

  end if

! Output number of ghost zones.

  if (rank==0) then
     write (*,'(A,I2)') ' Ghost zones  = ',ghost
     print *
  end if


! ****************************
! ***   REFINEMENT BOXES   ***
! ****************************

! Set up location and size of refinement boxes.

  call setupboxes


! ********************************
! ***   DOMAIN DECOMPOSITION   ***
! ********************************

! Do domain decomposition for parallel runs.

  call domaindecomp


! ***************************
! ***   ALLOCATE ARRAYS   ***
! ***************************

! Find out largest value of (Nrl,Nzl).  This is
! important since all arrays will have to be of the
! same size in order to make communications easier.

  Nrmaxl = 0
  Nzmaxl = 0

  do box=0,Nb
     do proc=0,size-1
        if (Nrl(box,proc)>Nrmaxl(box)) Nrmaxl(box)=Nrl(box,proc)
        if (Nzl(box,proc)>Nzmaxl(box)) Nzmaxl(box)=Nzl(box,proc)
     end do
  end do

! Allocate arrays and initialize to zero.

  call allocatearrays('on')


! ***********************************************************
! ***   FIGURE OUT POSITION OF AXIS, ORIGIN AND EQUATOR   ***
! ***********************************************************

! For several things, it is important to know the ownership
! of the axis, and also the ownership and z position
! of the origin and equator.

! Allocate arrays (origin,axis,eqz) and initialize to -1,
! which means that nobody owns them.

  allocate(origin(0:Nb,0:size-1))
  origin = -1

  allocate(axis(0:Nb,0:size-1))
  axis = -1

  allocate(eqz(0:Nb,0:size-1))
  eqz = -1

! Who owns the axis?
!
! 1) If rbox(box)=0, then all processors p such
!    that mod(proc,nprocr) = 0 own the axis.

  do box=0,Nb
     if (rbox(box)==0.d0) then
        do proc=0,size-1
           if (mod(proc,nprocr)==0) then
              axis(box,proc) = 1
           end if
        end do
     end if
  end do

! Who owns the origin and equator?
!
! 1) For equatorial symmetry this is easy, if zbox(box)=0
!    then processor 0 always owns the origin.
!
! 2) Also, for zbox(box)=0 processors with rank<nprocr
!    own the equator, and it corresponds to the value z=1.
!    Notice that this is in fact not precisely the equator,
!    since the grid is staggered, so there is no grid point
!    in the equator.

  if (eqsym) then

     do box=0,Nb
        if (zbox(box)==0.d0) then
           origin(box,0) = 1
           do proc=0,nprocr-1
              eqz(box,proc) = 1
           end do
        end if
     end do

! When there is no equatorial symmetry things are
! somewhat more complicated.  At the moment I
! assume that only boxes with zbox(box)=0 own the equator.
! The reason for doing this is that even
! if other refinement boxes intersect the 
! equator, they might not do so at all grid
! levels, and things can become quite confusing.

  else

!    Loop over boxes and processors.

     flag = .false.

     do box=0,Nb
        if (zbox(box)==0.d0) then

           do proc=0,size-1

!             Now loop over the points in the z direction
!             owned by the given processor.

              do j=1,Nzl(box,proc)-ghost

!                Find value of z at grid level 0.

                 aux = (dble(Nminl_z(box,proc)+j) - 0.5d0*dble(Nzbox(box)-ghost+1))*dz0

!                Now see if we have found the equator, or the first
!                grid point above it if the grid is staggered.
!                Also, only the first processor to own the equator
!                also owns the origin.

                 if ((aux==0.d0).or.(aux==0.5d0*dz0)) then
                    eqz(box,proc) = j
                    if (.not.flag) then
                       flag = .true.
                       origin(box,proc) = j
                    end if
                 end if

              end do
           end do

        end if
     end do

  end if


! *********************
! ***   EVOLUTION   ***
! *********************

  call evolve


! *****************************
! ***   DEALLOCATE ARRAYS   ***
! *****************************

! Deallocate grid arrays.

  call allocatearrays('off')

! Deallocate arrays with processor information.

  deallocate(Nrmaxl,Nminl_r,Nmaxl_r,Nrl)
  deallocate(Nzmaxl,Nminl_z,Nmaxl_z,Nzl)


! ************************
! ***   FINALIZE MPI   ***
! ************************

! When compiling without MPI this call does nothing at all.

  call MPI_FINALIZE(ierr)


! ***************
! ***   END   ***
! ***************

  if (rank==0) then
     print *,'------------------------------'
     print *
     print *, 'PROGRAM OLLINAXIS HAS FINISHED'
     print *
     print *, 'Have a nice day!'
     print *
     print *
     print *
  end if

  end program main


