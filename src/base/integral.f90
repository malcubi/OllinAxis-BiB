
  real(8) function integral(box,level,var)

! ******************************************
! ***   CALCULATION OF VOLUME INTEGRAL   ***
! ******************************************

! This function integrates the array "var" on a
! specific box and level determined by the routine
! that calls it.
!
! At the moment the integral is done using the
! trapezium rule.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer i,j
  integer imax,jmax
  integer box,level

  real(8) rmax,zmax
  real(8) aux

  real(8) :: var(1-ghost:Nrmaxl(box),1-ghost:Nzmaxl(box))


! *********************
! ***   INTEGRATE   ***
! *********************

! Initialize.

  integral = 0.d0

  rmax = r(Nr,0)
  zmax = z(0,Nz)

! Find maximum values of (i,j).

  if ((nprocr==1).or.(rmax>dble(Nrtotal-1)*dr))  then
     imax = Nr
  else
     imax = Nr - ghost + 1
  end if

  if ((nprocz==1).or.(zmax>dble(Nztotal-1)*dz)) then
     jmax = Nz
  else
     jmax = Nz - ghost + 1
  end if

! Add contributions from different grid cells.

  do j=1,jmax-1
     do i=1,imax-1
        integral = integral + var(i,j) + var(i,j+1) + var(i+1,j) + var(i+1,j+1)
     end do
  end do

! If we own the axis we need to include the correct contribution
! from the grid cells that stagger the axis.

  if (ownaxis) then
     do j=1,jmax-1
        integral = integral + 0.5d0*(var(0,j) + var(0,j+1) + var(1,j) + var(1,j+1))
     end do
  end if

! If we have equatorial symmetry we need to add the correct
! contribution from the grid cells that stagger the equator.

  if (eqsym.and.ownequator) then

!    Equator.

     do i=1,imax-1
        integral = integral + 0.5d0*(var(i,0) + var(i+1,0) + var(i,1) + var(i+1,1))
     end do

!    Origin.

     if (ownorigin) then
        integral = integral + 0.25d0*(var(0,0) + var(1,0) + var(0,1) + var(1,1))
     end if

  end if

! Multiply by normalization factor.

  aux = 0.25d0*dr*dz
  integral = integral*aux

! Add the contributions from the different processors when needed.

  if (size>1) then
     aux = integral
     call MPI_ALLREDUCE(aux,integral,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if

!  For equatorial symmetry multiply with 2.

   if (eqsym) then
      integral = 2.d0*integral
   end if


! ***************
! ***   END   ***
! ***************

  end function integral
