
subroutine wavextract

  ! *****************************************
  ! ***   GRAVITATIONAL WAVE EXTRACTION   ***
  ! *****************************************

  ! Routine for gravitational wave extraction
  ! using the Newmann-PeNr0ose formalism and
  ! the multipolar decomposition of the
  ! Weyl scalar Psi4.

  ! Include modules.

  use mpi
  use param
  use arrays
  use procinfo
  use harmonix

  ! Extra variables.

  implicit none

  logical firstcall
  logical interpflag                        ! Interpolation flag.
  real(8), dimension(1:3) :: radii          ! Evaluation radii.
  real(8), dimension(1:3,1:3) :: Almr,Almi  ! l,m modes, real and imaginary parts.
  !	real(8), dimension(1:3) :: Almrp,Almip ! l,m modes, real and imaginary parts.
  integer box,level,bb,ll                   ! Box and level counters.
  integer Ntheta			                    ! Number of points.
  integer spinw,l,m		                    ! Spin, azimuthal and angular numbers.
  integer i,j,k			                    ! Counters.
  real(8) interp                            ! Interpolation function.
  real(8) theta,dth
  real(8) iWeylpsi4	                       ! Interpolated value of the Weyl scalar psi4.
  real(8) r0,z0 	                          ! Interpolated point
  real(8) half,smallpi
  real(8) aux                               ! Auxiliary.
  character(4)  filen
  character(20) filestatus

  data firstcall / .true. /


  ! ************************
  ! ***   SANITY CHECK   ***
  ! ************************

  ! Check that the extraction radii
  ! are greater than zero and smaller than rmax.

  if (rad1<0) then
     if (rank==0) then
        print *
        print *, 'Extraction radius 1 must be greater than zero.'
        print *, 'Aborting! (subroutine wave_extract)'
        print *
     end if
     call die
  end if

  if (rad2<0) then
     if (rank==0) then
        print *
        print *, 'Extraction radius 2 must be greater than zero.'
        print *, 'Aborting! (subroutine wave_extract)'
        print *
     end if
     call die

  end if

  if (rad3<0) then
     if (rank==0) then
        print *
        print *, 'Extraction radius 3 must be greater than zero.'
        print *, 'Aborting! (subroutine wave_extract)'
        print *
     end if
     call die

  end if

  if (rad1>min(Nrtotal,Nztotal)*min(dr0,dz0)) then
     if (rank==0) then
        print *
        print *, 'Extraction radius 1 must be smaller than rmax.'
        print *, 'Aborting! (subroutine wave_extract)'
        print *
     end if
     call die
  end if

  if (rad2>min(Nrtotal,Nztotal)*min(dr0,dz0)) then
     if (rank==0) then
        print *
        print *, 'Extraction radius 2 must be smaller than rmax.'
        print *, 'Aborting! (subroutine wave_extract)'
        print *
     end if
     call die
  end if

  if (rad3>min(Nrtotal,Nztotal)*min(dr0,dz0)) then
     if (rank==0) then
        print *
        print *, 'Extraction radius 3 must be smaller than rmax.'
        print *, 'Aborting! (subroutine wave_extract)'
        print *
     end if
     call die
  end if

  ! *******************
  ! ***   NUMBERS   ***
  ! *******************

  half = 0.5d0

  smallpi = acos(-1.d0)


  ! **********************
  ! ***   INITIALIZE   ***
  ! **********************

  ! Find total number of points and dth
  ! depending on the symmetry.  At the moment
  ! I force a fixed number of points in theta
  ! in order to have reasonable accuracy.

  if (eqsym) then
     Ntheta = 250
     dth = half*(smallpi)/dble(Ntheta)
  else
     Ntheta = 500
     dth = (smallpi)/dble(Ntheta)
  end if

  ! Fill the observation radii

  radii(1) = rad1
  radii(2) = rad2
  radii(3) = rad3

  ! Set to zero the multipole components.
  Almr = 0.0D0
  Almi = 0.0D0

  ! Weyl scalar psi4 has spin weight s=-2.
  ! In axial symmetry the only non-zero modes are the m=0.
  ! This loop calculates the l=2,4,6,m=0 modes.
  ! This may be confusiong, but with the loops arranged like this,
  ! we only need to interpolate one time per point
  ! in order to obtain the value of the 2,4,6 multipole modes.
  ! The integration method is fourth order accurate.

  ! Calculate weyl scalar Psi4.

  call weyl
  ! Set spin-weight and angular numbers.
  spinw=-2
  m=0

  do j=1,3

     if (radii(j)/= 0.d0) then

        do i=0,Ntheta

           theta = dble(i*dth)
           r0 = radii(j)*sin(theta)
           z0 = radii(j)*cos(theta)

           ! *****************************************************
           ! ***   FIND GRID BOX AND LEVEL FOR INTERPOLATION   ***
           ! *****************************************************

           ! We need to interpolate at the highest level
           ! available at the current location.

           box = 0
           level = 0

           !           if (Nlmax>0) then
           !              do bb=0,Nb
           !                 do ll=1,Nl(bb)
           !                    if ((r0>rminl(bb,ll)+drl(ll)).and.(r0<rmaxl(bb,ll)-drl(ll)).and. &
           !                         (z0>zminl(bb,ll)+dzl(ll)).and.(z0<zmaxl(bb,ll)-dzl(ll))) then
           !                       box = bb
           !                       level = ll
           !                    end if
           !                 end do
           !              end do
           !           end if


           interpvar => grid(box,level)%WEYL_psi4
           aux = interp(box,level,r0,z0,interpflag)
           call MPI_ALLREDUCE(aux,iweylpsi4,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
           do k=1,3
              ! Azimuthal number

              l = 2*k

              if (i == 0 .or. i == Ntheta) then
                 Almr(j,k) = Almr(j,k) + iweylpsi4 * sin(theta)*SWSHR(spinw,l,m,theta,0.0D0)
              else
                 ! Real part
                 if (mod(i+1,2) == 0) then
                    Almr(j,k) = Almr(j,k) + 4.0D0*iweylpsi4 * sin(theta)*SWSHR(spinw,l,m,theta,0.0D0)
                 else
                    Almr(j,k) = Almr(j,k) + 2.0D0*iweylpsi4 * sin(theta)*SWSHR(spinw,l,m,theta,0.0D0)
                 end if
              end if


           end do

        end do

     end if

  end do

  if (eqsym) then
     Almr = 4.0D0*smallpi*Almr*dth/3.0d0
  else
     Almr = 2.0D0*smallpi*Almr*dth/3.0d0
  end if

  ! ********************************
  ! ***   SAVE MULTIPOLAR MODES  ***
  ! ********************************

  ! Save horizon data.

  if (rank==0) then

     !    Determine file status and extension.

     if (firstcall) then
        firstcall = .false.
        filestatus = 'replace'
     else
        filestatus = 'old'
     end if

     write(filen,'(A)') '.tl'

     !    Open files.

     if (filestatus == 'replace') then
        open(1,file=trim(directory)//'/psi4_20'//trim(filen),form='formatted', &
             status=filestatus)
        open(2,file=trim(directory)//'/psi4_40'//trim(filen),form='formatted', &
             status=filestatus)
        open(3,file=trim(directory)//'/psi4_60'//trim(filen),form='formatted', &
             status=filestatus)
        !       Write the extraction radii.
        write(1,"(A,3ES20.8)") "extraction radii=",rad1,rad2,rad3
        write(2,"(A,3ES20.8)") "extraction radii=",rad1,rad2,rad3
        write(3,"(A,3ES20.8)") "extraction radii=",rad1,rad2,rad3

     else
        open(1,file=trim(directory)//'/psi4_20'//trim(filen),form='formatted', &
             status=filestatus,position='append')
        open(2,file=trim(directory)//'/psi4_40'//trim(filen),form='formatted', &
             status=filestatus,position='append')
        open(3,file=trim(directory)//'/psi4_60'//trim(filen),form='formatted', &
             status=filestatus,position='append')
     end if


     !    Save multipolar mode. The columns order are time, rad1,rad2,rad3

     write(1,"(4ES16.8)") t(0,0),Almr(1,1),Almr(2,1),Almr(3,1)
     write(2,"(4ES16.8)") t(0,0),Almr(1,2),Almr(2,2),Almr(2,2)
     write(3,"(4ES16.8)") t(0,0),Almr(1,3),Almr(3,3),Almr(3,3)

     !    Close files
     close(1)
     close(2)
     close(3)

  end if

end subroutine wavextract
