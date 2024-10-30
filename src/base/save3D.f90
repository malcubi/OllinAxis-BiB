
  subroutine save3D (varname,outdir)

! This subroutine transforms an axisymmetric 2D array to its cartesian 3D form.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

  implicit none

  integer i,j,k                      ! Counters
  integer box,level,bb,ll            ! Box and level counters.
  real(8),dimension(0:Nxx) :: xx      ! Grid in x direction
  real(8),dimension(0:Nyy) :: yy      ! Grid in y direction
  real(8),dimension(0:Nzz) :: zz      ! Grid in y direction
  real(8) r0,z0                         ! Radius
  real(8) tiny
  character(len=*)  varname,outdir
  character(len=9)  filen

  real(8) interp                     ! Interpolation function.
  logical interpflag                 ! Interpolation flag.
  real(8) interp_f                   ! Interpolated function.
  real(8) aux

  !real(8), dimension (1-ghost:Nrmax,1-ghost:Nzmax) :: array2D_to_3D

  tiny = 1.d-30

! Check if xmax,ymax is contained in the grid.

  if ((Nxx*dxx)**2+(Nyy*dyy)**2>rmaxl(0,0)**2) then
     if (rank==0) then
        print *
        print *, 'Cartesian grid is out of axisymmetric boundaries.'
        print *, 'Aborting! (subroutine save3D)'
        print *
     end if
     call die
  end if

  if ((Nzz*dz0)>zmaxl(0,0)) then
     if (rank==0) then
        print *
        print *, 'Cartesian grid is out of axisymmetric boundaries.'
        print *, 'Aborting! (subroutine save3D)'
        print *
     end if
     call die
  end if


! Fill cartesian arrays

  if (sgrid) then

     do i=0,Nxx
        xx(i) = (i+0.5d0)*dxx
     end do

     do i=0,Nyy
        yy(i) = (i+0.5d0)*dyy
     end do

     if (eqsym) then
        do i=0,Nzz
           zz(i) = (i+0.5d0)*dz0
        end do
     end if

  else

     do i=0,Nxx
        xx(i) = dble(i*dxx)
     end do

     do i=0,Nyy
        yy(i) = dble(i*dyy)
     end do

     if (eqsym) then
        do i=0,Nzz
           zz(i) = dble(i*dz0)
        end do
     end if

  end if




! Open file

  call grabarray(varname)
  
  open(1,file=trim(outdir)//'/'//varname//'.3D',form='formatted')
  !open(1,file='psi.3D',form='formatted')

!    Write data.

  do k = 0,Nzz   
     do j = 0,Nyy
        do i = 0,Nxx

           r0 = sqrt(xx(i)**2 + yy(j)**2)
           z0 = zz(k)

! *****************************************************
! ***   FIND GRID BOX AND LEVEL FOR INTERPOLATION   ***
! *****************************************************

! We need to interpolate at the highest level
! available at the current location.

           box = 0
           level = 0

           if (Nlmax>0) then
              do bb=0,Nb
                 do ll=1,Nl(bb)
                    if ((r0>rminl(bb,ll)+drl(ll)).and.(r0<rmaxl(bb,ll)-drl(ll)).and. &
                        (z0>zminl(bb,ll)+dzl(ll)).and.(z0<zmaxl(bb,ll)-dzl(ll))) then
                       box = bb
                       level = ll
                    end if
                 end do
              end do
           end if

! Interpolate array

          interpvar => grabvar!grid(box,level)%grabvar
          aux = interp(box,level,r0,z0,interpflag)
          call MPI_ALLREDUCE(aux,interp_f,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

          write(1,"(3ES23.15)") xx(i),yy(j),zz(k),interp_f+tiny

        end do
     end do
     !write(1,*)
  end do

  close(1)

  end subroutine save3D
