!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/interp.f90,v 1.17 2021/02/24 22:34:06 malcubi Exp $

  real(8) function interp(box,level,r0,z0,flag)

! *************************
! ***   INTERPOLATION   ***
! *************************

! This function does linear or cubic interpolation
! on the gri¡d corresponding to box b and level l.
!
! Notice that the array that is interpolated is always
! the pointer array "interpvar", so one needs to make sure
! that this points to the variable one really wants to
! interpolate before calling this function.

! Include modules.

  use param
  use arrays
  use procinfo

! Declare variables
  
  implicit none

  logical flag            ! Flag to indicate ownership of interpolating point.
  logical forcelinear     ! Flag to force bilinear interpolation if we are too close to a boundary.

  integer i,j,m,n         ! Counters.
  integer box,level       ! Box number and level counters.
  integer i0,j0           ! Closest grid point with r<r0 and z<z0.

  real(8) r0,z0           ! Point at which we are interpolating.
  real(8) deltar,deltaz   ! Distance to grid point (i0,j0).
  real(8) rmin,rmax       ! Minimum and maximum r for curent grid.
  real(8) zmin,zmax       ! Minimum and maximum z for curent grid.
  real(8) small
  real(8) aux,auxr,auxz

  real(8) fa(-1:2,-1:2)   ! Small local array with values at nearest neighbors.


! *******************
! ***   NUMBERS   ***
! *******************

  small = 1.d-10


! ************************
! ***   SANITY CHECK   ***
! ************************

! First check that the interpolating point is inside box b and level l.

  rmin = rminl(box,level) - small
  rmax = rmaxl(box,level)

  zmin = zminl(box,level) - small
  zmax = zmaxl(box,level)

  if ((r0<rmin).or.(r0>rmax).or.(z0<zmin).or.(z0>zmax)) then
     if (rank==0) then
        print *
        print *, 'interp.f90: Interpolating point is outside the local grid:'
        write(*,'(A,2ES14.6)') ' Interpolating position =',r0,z0
        write(*,'(A,2ES14.6)') ' Grid edges in r = ',rmin,rmax
        write(*,'(A,2ES14.6)') ' Grid edges in z = ',zmin,zmax
        write(*,'(A,I0)')      ' Refinement box  = ',box
        write(*,'(A,I0)')      ' Grid level      = ',level
        print *
        print *, 'Aborting! (subroutine interp)'
        print *
     end if
     call die
  end if


! *******************************************************
! ***   CHECK IF THE LOCAL PROCESSOR OWNS THE POINT   ***
! *******************************************************

! For multiple processor runs check that you
! own the point to be interpolated.
!
! Notice that we need to make sure that only one
! processor owns the point, which might be a problem
! close to inter-processor boundaries where there is
! an overlap. When that happens here I choose the
! processor with lower (r,z) as owner of the point
! (but carefull with processor on the lower edge
! of the grid).

  if (size>1) then

!    Limits in r.

     if (mod(rank,nprocr)==0) then
        rmin = grid(box,level)%r(1-ghost,0) - small
     else
        rmin = grid(box,level)%r(1,0)
     end if

     if (mod(rank+1,nprocr)==0) then
        rmax = grid(box,level)%r(Nrl(box,rank),0)
     else
        rmax = grid(box,level)%r(Nrl(box,rank)-ghost+1,0)
     end if

!    Limits in z.

     if (rank<nprocr) then
        zmin = grid(box,level)%z(0,1-ghost) - small
     else
        zmin = grid(box,level)%z(0,1)
     end if

     if (rank>size-nprocr-1) then
        zmax = grid(box,level)%z(0,Nzl(box,rank))
     else
        zmax = grid(box,level)%z(0,Nzl(box,rank)-ghost+1)
     end if

!    If you don't own the grid point set the flag=.false.,
!    set interp=0, and return.  Otherwise set flag=.true.
!    and continue with the interpolation.

     if ((r0<=rmin).or.(r0>rmax).or.(z0<=zmin).or.(z0>zmax)) then
        flag = .false.
        interp = 0.d0
        return
     else
        flag = .true.
     end if

  end if


! *******************************************************
! ***   IDENTIFY GRID CELL WHERE (r0,z0) IS LOCATED   ***
! *******************************************************

! Notice that if we get to this point the current processor
! owns the point.

  i0 = int((r0-grid(box,level)%r(1-ghost,0))/drl(level)) + 1 - ghost
  deltar = r0 - grid(box,level)%r(i0,0)

  j0 = int((z0-grid(box,level)%z(0,1-ghost))/dzl(level)) + 1 - ghost
  deltaz = z0 - grid(box,level)%z(0,j0)

! Check that the point to be interpolated is indeed inside the cell.

  if ((deltar<-small).or.(deltaz<-small).or. &
      (deltar>drl(level)+small).or.(deltaz>dzl(level)+small)) then
     print *, 'This should never happen!'
     print *, 'Box = ',box,' Level = ',level
     print *, drl(level),deltar,dzl(level),deltaz
     print *, 'Aborting! (subroutine interp)'
     print *
  end if

! Check if we are in fact on top of a grid point, in which
! case we just copy the value and jump out.

  if ((deltar==0.d0).and.(deltaz==0.d0)) then
     interp = interpvar(i0,j0)
     return
  end if

! Now check if we are two close to a boundary for cubic
! interpolation.

  if (((rmaxl(box,level)-r0)<drl(level)).or.((zmaxl(box,level)-z0)<dzl(level)).or. &
      ((r0-rminl(box,level))<drl(level)).or.((z0-zminl(box,level))<dzl(level))) then
     forcelinear = .true.
  else
     forcelinear = .false.
  end if


! ********************************************************
! ***   COPY THE VALUES AT NEAREST GRID POINTS TO fa   ***
! ********************************************************

! Copy the values of the array to be interpolated on the
! nearest grid points to the small local array 'fa':
!
!    o        o        o        o   j0+2
!
!    o        o        o        o   j0+1
! 
!    o        o        o        o   j0
!
!    o        o        o        o   j0-1
!
!   i0-1     i0       i0+1     i0+2


  if (forcelinear) then
     do i=0,1
        do j=0,1
           fa(i,j) = interpvar(i0+i,j0+j)
        end do
     end do
  else
     do i=-1,2
        do j=-1,2
           fa(i,j) = interpvar(i0+i,j0+j)
        end do
     end do
  end if


! *******************************
! ***   LINEAR INTERPOLATON   ***
! *******************************

! Simple bilinear interpolation.

  if ((interporder=="bilinear").or.forcelinear) then

     auxr = deltar/drl(level)
     auxz = deltaz/dzl(level)

     interp = fa(0,0)*(1.d0-auxr)*(1.d0-auxz) + fa(1,0)*auxr*(1.d0-auxz) &
            + fa(0,1)*(1.d0-auxr)*auxz + fa(1,1)*auxr*auxz


! ******************************
! ***   CUBIC INTERPOLATON   ***
! ******************************

! Bicubic Lagrange polynomial interpolation.

  else if (interporder=="bicubic") then

     auxr = deltar/drl(level)
     auxz = deltaz/dzl(level)

     interp = 0.d0

     do i=-1,2
        do j=-1,2

           aux = 1.d0

           do m=-1,2
              if (m==i) cycle
              aux = aux*(auxr-dble(m))/dble(i-m)
           end do

           do n=-1,2
              if (n==j) cycle
              aux = aux*(auxz-dble(n))/dble(j-n)
           end do

           interp = interp + aux*fa(i,j)

        end do
     end do

  end if


! ***************
! ***   END   ***
! ***************

  end function interp




