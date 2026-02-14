
  module derivatives

  contains

! I use a Fortran module for the array-valued functions that
! calculate the first and second derivatives of arrays.
!
! It turns out that defining array-valued functions is not
! trivial, and putting them inside a module seems to be
! a way to solve the problem.



  function diff1r(symr)

! **********************************************************
! ***   CALCULATE FIRST DERIVATIVE WITH RESPECT TO RHO   ***
! **********************************************************

! Load modules.

  use param
  use arrays
  use procinfo

! Declare variables.

  implicit none

  integer i,symr

  real(8) idr,hidr

  real(8), dimension (1-ghost:Nrmax,1-ghost:Nzmax) :: diff1r


! *******************
! ***   NUMBERS   ***
! *******************

  idr  = 1.d0/dr
  hidr = 0.5d0*idr


! ************************
! ***   SECOND ORDER   ***
! ************************

  if (order=="two") then

!    Interior points: second order centered first derivative.

     do i=2-ghost,Nr-1
        diff1r(i,:) = hidr*(diffvar(i+1,:) - diffvar(i-1,:))
     end do

!    Point Nr:  2nd order fully one-sided.

     i = Nr
     diff1r(i,:) = + hidr*(3.d0*diffvar(i,:) - 4.d0*diffvar(i-1,:) + diffvar(i-2,:))


! ************************
! ***   FOURTH ORDER   ***
! ************************

  else if (order=="four") then

!    Interior points: fourth order centered first derivative.

     do i=1,Nr-2
        diff1r(i,:) = 0.25d0*idr*(8.d0*(diffvar(i+1,:) - diffvar(i-1,:)) &
                    - (diffvar(i+2,:) - diffvar(i-2,:)))/3.d0
     end do

!    Point Nr-1:  4th order semi one-sided (1 point to the right, 3 to the left).

     i = Nr-1
     diff1r(i,:) = idr*(3.d0*diffvar(i+1,:) + 10.d0*diffvar(i,:) - 18.d0*diffvar(i-1,:) &
              + 6.d0*diffvar(i-2,:) - diffvar(i-3,:))/12.d0

!    Point Nr:  4th order fully one-sided.

     i = Nr
     diff1r(i,:) = idr*(25.d0*diffvar(i,:) - 48.d0*diffvar(i-1,:) &
              + 36.d0*diffvar(i-2,:) - 16.d0*diffvar(i-3,:) + 3.d0*diffvar(i-4,:))/12.d0

  end if


! ******************************
! ***   SYMMETRIES AT AXIS   ***
! ******************************

  if (ownaxis) then
     do i=1,ghost
        diff1r(1-i,:) = - sign(1,symr)*diff1r(i,:)
     end do
  end if


! ***********************
! ***   SYNCHRONIZE   ***
! ***********************

! Here I only need to sync on r direction.

  if (nprocr>1) then
     call syncr(diff1r)
  end if


! ************************
! ***   END FUNCTION   ***
! ************************

  end function diff1r





  function diff1z(symz)

! ********************************************************
! ***   CALCULATE FIRST DERIVATIVE WITH RESPECT TO Z   ***
! ********************************************************

! Load modules.

  use param
  use arrays
  use procinfo

! Declare variables.

  implicit none

  integer j,symz

  real(8) idz,hidz

  real(8), dimension (1-ghost:Nrmax,1-ghost:Nzmax) :: diff1z


! *******************
! ***   NUMBERS   ***
! *******************

  idz  = 1.d0/dz
  hidz = 0.5d0*idz


! ************************
! ***   SECOND ORDER   ***
! ************************

  if (order=="two") then

!    Interior points: second order centered first derivative.

     do j=2-ghost,Nz-1
        diff1z(:,j) = hidz*(diffvar(:,j+1) - diffvar(:,j-1))
     end do

!    Boundary with one-sided differences.

     j = 1-ghost
     diff1z(:,j) = - hidz*(3.d0*diffvar(:,j) - 4.d0*diffvar(:,j+1) + diffvar(:,j+2))

     j = Nz
     diff1z(:,j) = + hidz*(3.d0*diffvar(:,j) - 4.d0*diffvar(:,j-1) + diffvar(:,j-2))


! ************************
! ***   FOURTH ORDER   ***
! ************************

  else if (order=="four") then

!    Interior points: fourth order centered first derivative.

     do j=3-ghost,Nz-2
        diff1z(:,j) = 0.25d0*idz*(8.d0*(diffvar(:,j+1) - diffvar(:,j-1)) &
                    - (diffvar(:,j+2) - diffvar(:,j-2)))/3.d0
     end do

!    Point 2-ghost: 4th order semi one-sided (1 point to the left, 3 to the right).

     j = 2-ghost
     diff1z(:,j) = - idz*(3.d0*diffvar(:,j-1) + 10.d0*diffvar(:,j) - 18.d0*diffvar(:,j+1) &
              + 6.d0*diffvar(:,j+2) - diffvar(:,j+3))/12.d0

!    Point 1-ghost: 4th order fully one-sided.

     j = 1-ghost
     diff1z(:,j) = - idz*(25.d0*diffvar(:,j) - 48.d0*diffvar(:,j+1) &
              + 36.d0*diffvar(:,j+2) - 16.d0*diffvar(:,j+3) + 3.d0*diffvar(:,j+4))/12.d0

!    Point Nz-1: 4th order semi one-sided (1 point to the right, 3 to the left).

     j = Nz-1
     diff1z(:,j) = + idz*(3.d0*diffvar(:,j+1) + 10.d0*diffvar(:,j) - 18.d0*diffvar(:,j-1) &
              + 6.d0*diffvar(:,j-2) - diffvar(:,j-3))/12.d0

!    Point Nz:  4th order fully one-sided.

     j = Nz
     diff1z(:,j) = + idz*(25.d0*diffvar(:,j) - 48.d0*diffvar(:,j-1) &
              + 36.d0*diffvar(:,j-2) - 16.d0*diffvar(:,j-3) + 3.d0*diffvar(:,j-4))/12.d0

  end if


! *********************************
! ***   SYMMETRIES ON EQUATOR   ***
! *********************************

  if (eqsym.and.ownequator) then
     do j=1,ghost
        diff1z(:,1-j) = - sign(1,symz)*diff1z(:,j)
     end do
  end if


! ***********************
! ***   SYNCHRONIZE   ***
! ***********************

! Here I only need to sync on z direction.

  if (nprocz>1) then
     call syncz(diff1z)
  end if


! ************************
! ***   END FUNCTION   ***
! ************************

  end function diff1z








  function diff2r(symr)

! **********************************************************
! ***   CALCULATE FIRST DERIVATIVE WITH RESPECT TO RHO   ***
! **********************************************************

! Load modules.

  use param
  use arrays
  use procinfo

! Declare variables.

  implicit none

  integer i,symr

  real(8) idr,idr2

  real(8), dimension (1-ghost:Nrmax,1-ghost:Nzmax) :: diff2r


! *******************
! ***   NUMBERS   ***
! *******************

  idr  = 1.d0/dr
  idr2 = idr**2


! ************************
! ***   SECOND ORDER   ***
! ************************

  if (order=="two") then

!    Interior points: second order centered second derivative.

     do i=2-ghost,Nr-1
        diff2r(i,:) = idr2*(diffvar(i+1,:) - 2.d0*diffvar(i,:) + diffvar(i-1,:))
     end do

!    Boundary with one-sided differences.

     i = 1-ghost
     diff2r(i,:) = idr2*(2.d0*diffvar(i,:) - 5.d0*diffvar(i+1,:) &
                 + 4.d0*diffvar(i+2,:) - diffvar(i+3,:))

     i = Nr
     diff2r(i,:) = idr2*(2.d0*diffvar(i,:) - 5.d0*diffvar(i-1,:) &
                 + 4.d0*diffvar(i-2,:) - diffvar(i-3,:))


! ************************
! ***   FOURTH ORDER   ***
! ************************

  else if (order=="four") then

!    Interior points: fourth order centered second derivative.

     do i=3-ghost,Nr-2
        diff2r(i,:) = - idr2*(30.d0*diffvar(i,:) &
                    - 16.d0*(diffvar(i+1,:) + diffvar(i-1,:)) &
                    +        diffvar(i+2,:) + diffvar(i-2,:))/12.d0
     end do

!    Point Nr-1:  4th order semi one-sided.

     i = Nr-1
     diff2r(i,:) = idr2*(10.d0*diffvar(i+1,:) - 15.d0*diffvar(i,:) - 4.d0*diffvar(i-1,:) &
                 + 14.d0*diffvar(i-2,:) - 6.d0*diffvar(i-3,:) + diffvar(i-4,:))/12.d0

!    Point Nr:  4th order fully one-sided.

     i = Nr
     diff2r(i,:) = idr2*(45.d0*diffvar(i,:) - 154.d0*diffvar(i-1,:) + 214.d0*diffvar(i-2,:) &
                 - 156.d0*diffvar(i-3,:) + 61.d0*diffvar(i-4,:) - 10.d0*diffvar(i-5,:))/12.d0

  end if


! ******************************
! ***   SYMMETRIES ON AXIS   ***
! ******************************

  if (ownaxis) then
     do i=1,ghost
        diff2r(1-i,:) = sign(1,symr)*diff2r(i,:)
     end do
  end if


! ***********************
! ***   SYNCHRONIZE   ***
! ***********************

! Here I only need to sync on r direction.

  if (nprocr>1) then
     call syncr(diff2r)
  end if


! ************************
! ***   END FUNCTION   ***
! ************************

  end function diff2r







  function diff2z(symz)

! ********************************************************
! ***   CALCULATE FIRST DERIVATIVE WITH RESPECT TO Z   ***
! ********************************************************

! Load modules.

  use param
  use arrays
  use procinfo

! Declare variables.

  implicit none

  integer j,symz

  real(8) idz,idz2

  real(8), dimension (1-ghost:Nrmax,1-ghost:Nzmax) :: diff2z


! *******************
! ***   NUMBERS   ***
! *******************

  idz  = 1.d0/dz
  idz2 = idz**2


! ************************
! ***   SECOND ORDER   ***
! ************************

  if (order=="two") then

!    Interior points: second order centered second derivative.

     do j=2-ghost,Nz-1
        diff2z(:,j) = idz2*(diffvar(:,j+1) - 2.d0*diffvar(:,j) + diffvar(:,j-1))
     end do

!    Boundary with one-sided diferences.

     j = 1-ghost
     diff2z(:,j) = idz2*(2.d0*diffvar(:,j) - 5.d0*diffvar(:,j+1) &
                 + 4.d0*diffvar(:,j+2) - diffvar(:,j+3))

     j = Nz
     diff2z(:,j) = idz2*(2.d0*diffvar(:,j) - 5.d0*diffvar(:,j-1) &
                 + 4.d0*diffvar(:,j-2) - diffvar(:,j-3))


! ************************
! ***   FOURTH ORDER   ***
! ************************

  else if (order=="four") then

!    Interior points: fourth order centered second derivative.

     do j=3-ghost,Nz-2
        diff2z(:,j) = - idz2*(30.d0*diffvar(:,j) &
                    - 16.d0*(diffvar(:,j+1) + diffvar(:,j-1)) &
                    +        diffvar(:,j+2) + diffvar(:,j-2))/12.d0
     end do

!    Point Nz-1:  4th order semi one-sided.

     j = 2-ghost
     diff2z(:,j) = idz2*(10.d0*diffvar(:,j-1) - 15.d0*diffvar(:,j) - 4.d0*diffvar(:,j+1) &
                 + 14.d0*diffvar(:,j+2) - 6.d0*diffvar(:,j+3) + diffvar(:,j+4))/12.d0

!    Point 1-ghost:  4th order fully one-sided.

     j = 1-ghost
     diff2z(:,j) = idz2*(45.d0*diffvar(:,j) - 154.d0*diffvar(:,j+1) + 214.d0*diffvar(:,j+2) &
                 - 156.d0*diffvar(:,j+3) + 61.d0*diffvar(:,j+4) - 10.d0*diffvar(:,j+5))/12.d0

!    Point Nz-1:  4th order semi one-sided.

     j = Nz-1
     diff2z(:,j) = idz2*(10.d0*diffvar(:,j+1) - 15.d0*diffvar(:,j) - 4.d0*diffvar(:,j-1) &
                 + 14.d0*diffvar(:,j-2) - 6.d0*diffvar(:,j-3) + diffvar(:,j-4))/12.d0

!    Point Nz:  4th order fully one-sided.

     j = Nz
     diff2z(:,j) = idz2*(45.d0*diffvar(:,j) - 154.d0*diffvar(:,j-1) + 214.d0*diffvar(:,j-2) &
                 - 156.d0*diffvar(:,j-3) + 61.d0*diffvar(:,j-4) - 10.d0*diffvar(:,j-5))/12.d0

  end if


! *********************************
! ***   SYMMETRIES ON EQUATOR   ***
! *********************************

  if (eqsym.and.ownequator) then
     do j=1,ghost
        diff2z(:,1-j) = sign(1,symz)*diff2z(:,j)
     end do
  end if


! ***********************
! ***   SYNCHRONIZE   ***
! ***********************

! Here I only need to sync on z direction.

  if (nprocz>1) then
     call syncz(diff2z)
  end if


! ************************
! ***   END FUNCTION   ***
! ************************

  end function diff2z







  function diff2rz(symr,symz)

! ********************************************************
! ***   CALCULATE FIRST DERIVATIVE WITH RESPECT TO Z   ***
! ********************************************************

! Load modules.

  use param
  use arrays
  use procinfo

! Declare variables.

  implicit none

  integer i,j,symr,symz

  real(8) idr,idz,idrz

  real(8), dimension (1-ghost:Nrmax,1-ghost:Nzmax) :: diff2rz


! *******************
! ***   NUMBERS   ***
! *******************

  idr  = 1.d0/dr
  idz  = 1.d0/dz

  idrz = idr*idz


! ************************
! ***   SECOND ORDER   ***
! ************************

  if (order=="two") then

!    Interior points: second order centered mixed derivative.

     do i=2-ghost,Nr-1
        do j=2-ghost,Nz-1
           diff2rz(i,j) = 0.25d0*idrz &
              *(diffvar(i+1,j+1) + diffvar(i-1,j-1) &
              - diffvar(i+1,j-1) - diffvar(i-1,j+1))
        end do
     end do

!    Boundaries on r.

     do j=2-ghost,Nz-1
        i = 1-ghost
        diff2rz(i,j) = - 0.25d0*idrz &
           *((3.d0*diffvar(i,j+1) - 4.d0*diffvar(i+1,j+1) + diffvar(i+2,j+1)) &
           - (3.d0*diffvar(i,j-1) - 4.d0*diffvar(i+1,j-1) + diffvar(i+2,j-1)))
        i = Nr
        diff2rz(i,j) = + 0.25d0*idrz &
           *((3.d0*diffvar(i,j+1) - 4.d0*diffvar(i-1,j+1) + diffvar(i-2,j+1)) &
           - (3.d0*diffvar(i,j-1) - 4.d0*diffvar(i-1,j-1) + diffvar(i-2,j-1)))
     end do

!    Boundaries on z.

     do i=2-ghost,Nr-1
        j = 1-ghost
        diff2rz(i,j) = - 0.25d0*idrz &
           *((3.d0*diffvar(i+1,j) - 4.d0*diffvar(i+1,j+1) + diffvar(i+1,j+2)) &
           - (3.d0*diffvar(i-1,j) - 4.d0*diffvar(i-1,j+1) + diffvar(i-1,j+2)))
        j = Nz
        diff2rz(i,j) = + 0.25d0*idrz &
           *((3.d0*diffvar(i+1,j) - 4.d0*diffvar(i+1,j-1) + diffvar(i+1,j-2)) &
           - (3.d0*diffvar(i-1,j) - 4.d0*diffvar(i-1,j-1) + diffvar(i-1,j-2)))
     end do

!    Corners.

     i = 1-ghost; j = 1-ghost
     diff2rz(i,j) = + 0.25d0*idrz*(9.d0*diffvar(i,j) &
                  +  16.d0*diffvar(i+1,j+1) + diffvar(i+2,j+2) &
                  - 12.d0*(diffvar(i  ,j+1) + diffvar(i+1,j  )) &
                  +  3.d0*(diffvar(i  ,j+2) + diffvar(i+2,j  )) &
                  -  4.d0*(diffvar(i+1,j+2) + diffvar(i+2,j+1)))

     i = 1-ghost; j = Nz
     diff2rz(i,j) = - 0.25d0*idrz*(9.d0*diffvar(i,j) &
                   +  16.d0*diffvar(i+1,j-1) + diffvar(i+2,j-2)&
                   - 12.d0*(diffvar(i+1,j  ) + diffvar(i  ,j-1)) &
                   +  3.d0*(diffvar(i+2,j  ) + diffvar(i  ,j-2)) &
                   -  4.d0*(diffvar(i+2,j-1) + diffvar(i+1,j-2)))

     i = Nr; j = 1-ghost
     diff2rz(i,j) = - 0.25d0*idrz*(9.d0*diffvar(i,j) &
                  +  16.d0*diffvar(i-1,j+1) + diffvar(i-2,j+2)&
                  - 12.d0*(diffvar(i  ,j+1) + diffvar(i-1,j  )) &
                  +  3.d0*(diffvar(i  ,j+2) + diffvar(i-2,j  )) &
                  -  4.d0*(diffvar(i-1,j+2) + diffvar(i-2,j+1)))

     i = Nr; j = Nz
     diff2rz(i,j) = + 0.25d0*idrz*(9.d0*diffvar(i,j) &
                  +  16.d0*diffvar(i-1,j-1) + diffvar(i-2,j-2) &
                  - 12.d0*(diffvar(i  ,j-1) + diffvar(i-1,j  )) &
                  +  3.d0*(diffvar(i  ,j-2) + diffvar(i-2,j  )) &
                  -  4.d0*(diffvar(i-1,j-2) + diffvar(i-2,j-1)))


! ************************
! ***   FOURTH ORDER   ***
! ************************

  else if (order=="four") then

!    Interior points: fourth order centered mixed derivative.

     do j=3-ghost,Nz-2
        do i=3-ghost,Nr-2
           diff2rz(i,j) = idrz*((diffvar(i+1,j+1) + diffvar(i-1,j-1) &
                               - diffvar(i+1,j-1) - diffvar(i-1,j+1))/3.d0 &
                              - (diffvar(i+2,j+2) + diffvar(i-2,j-2) &
                               - diffvar(i+2,j-2) - diffvar(i-2,j+2))/48.d0)
        end do
     end do

!    Boundaries on r (second order at the moment).

     do j=2-ghost,Nz-1
        i = 2-ghost
        diff2rz(i,j) = + 0.25d0*idrz &
           *(diffvar(i+1,j+1) + diffvar(i-1,j-1) &
           - diffvar(i+1,j-1) - diffvar(i-1,j+1))
        i = 1-ghost
        diff2rz(i,j) = - 0.25d0*idrz &
           *((3.d0*diffvar(i,j+1) - 4.d0*diffvar(i+1,j+1) + diffvar(i+2,j+1)) &
           - (3.d0*diffvar(i,j-1) - 4.d0*diffvar(i+1,j-1) + diffvar(i+2,j-1)))
        i = Nr-1
        diff2rz(i,j) = + 0.25d0*idrz &
           *(diffvar(i+1,j+1) + diffvar(i-1,j-1) &
           - diffvar(i+1,j-1) - diffvar(i-1,j+1))
        i = Nr
        diff2rz(i,j) = + 0.25d0*idrz &
           *((3.d0*diffvar(i,j+1) - 4.d0*diffvar(i-1,j+1) + diffvar(i-2,j+1)) &
           - (3.d0*diffvar(i,j-1) - 4.d0*diffvar(i-1,j-1) + diffvar(i-2,j-1)))
     end do

!    Boundaries on z (second order at the moment).

     do i=2-ghost,Nr-1
        j = 2-ghost
        diff2rz(i,j) = + 0.25d0*idrz &
           *(diffvar(i+1,j+1) + diffvar(i-1,j-1) &
           - diffvar(i-1,j+1) - diffvar(i+1,j-1))
        j = 1-ghost
        diff2rz(i,j) = - 0.25d0*idrz &
           *((3.d0*diffvar(i+1,j) - 4.d0*diffvar(i+1,j+1) + diffvar(i+1,j+2)) &
           - (3.d0*diffvar(i-1,j) - 4.d0*diffvar(i-1,j+1) + diffvar(i-1,j+2)))
        j = Nz-1
        diff2rz(i,j) = + 0.25d0*idrz &
           *(diffvar(i+1,j+1) + diffvar(i-1,j-1) &
           - diffvar(i-1,j+1) - diffvar(i+1,j-1))
        j = Nz
        diff2rz(i,j) = + 0.25d0*idrz &
           *((3.d0*diffvar(i+1,j) - 4.d0*diffvar(i+1,j-1) + diffvar(i+1,j-2)) &
           - (3.d0*diffvar(i-1,j) - 4.d0*diffvar(i-1,j-1) + diffvar(i-1,j-2)))
     end do

!    Corners (second order at the moment).

     i = 1-ghost; j = 1-ghost
     diff2rz(i,j) = + 0.25d0*idrz*(9.d0*diffvar(i,j) &
                  +  16.d0*diffvar(i+1,j+1) + diffvar(i+2,j+2) &
                  - 12.d0*(diffvar(i  ,j+1) + diffvar(i+1,j  )) &
                  +  3.d0*(diffvar(i  ,j+2) + diffvar(i+2,j  )) &
                  -  4.d0*(diffvar(i+1,j+2) + diffvar(i+2,j+1)))

     i = 1-ghost; j = Nz
     diff2rz(i,j) = - 0.25d0*idrz*(9.d0*diffvar(i,j) &
                   +  16.d0*diffvar(i+1,j-1) + diffvar(i+2,j-2)&
                   - 12.d0*(diffvar(i+1,j  ) + diffvar(i  ,j-1)) &
                   +  3.d0*(diffvar(i+2,j  ) + diffvar(i  ,j-2)) &
                   -  4.d0*(diffvar(i+2,j-1) + diffvar(i+1,j-2)))

     i = Nr; j = 1-ghost
     diff2rz(i,j) = - 0.25d0*idrz*(9.d0*diffvar(i,j) &
                  +  16.d0*diffvar(i-1,j+1) + diffvar(i-2,j+2)&
                  - 12.d0*(diffvar(i  ,j+1) + diffvar(i-1,j  )) &
                  +  3.d0*(diffvar(i  ,j+2) + diffvar(i-2,j  )) &
                  -  4.d0*(diffvar(i-1,j+2) + diffvar(i-2,j+1)))

     i = Nr; j = Nz
     diff2rz(i,j) = + 0.25d0*idrz*(9.d0*diffvar(i,j) &
                  +  16.d0*diffvar(i-1,j-1) + diffvar(i-2,j-2) &
                  - 12.d0*(diffvar(i  ,j-1) + diffvar(i-1,j  )) &
                  +  3.d0*(diffvar(i  ,j-2) + diffvar(i-2,j  )) &
                  -  4.d0*(diffvar(i-1,j-2) + diffvar(i-2,j-1)))

  end if


! ******************************************
! ***   SYMMETRIES ON AXIS AND EQUATOR   ***
! ******************************************

! Symmetries in z for equatorial symmetry.
! Must be done first to ensure that the corners
! end up with the correct values.

  if (eqsym.and.ownequator) then
     do j=1,ghost
        diff2rz(:,1-j) = - sign(1,symz)*diff2rz(:,j)
     end do
  end if

! Symmetries in r.

  if (ownaxis) then
     do i=1,ghost
        diff2rz(1-i,:) = - sign(1,symr)*diff2rz(i,:)
     end do
  end if


! ***********************
! ***   SYNCHRONIZE   ***
! ***********************

! Here we need to do a full sync.

  if (size>1) then
     call sync(diff2rz)
  end if


! ************************
! ***   END FUNCTION   ***
! ************************

  end function diff2rz






! **********************
! ***   END MODULE   ***
! **********************

  end module derivatives


