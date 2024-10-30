!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/derivadvect.f90,v 1.3 2019/10/04 15:49:36 malcubi Exp $

  module derivadvect

  contains

! I use a Fortran module for the array-valued functions that
! calculate the first and second derivatives of arrays.
!
! It turns out that defining array-valued functions is not
! trivial, and putting them inside a module seems to be
! a way to solve the problem.



  function diffadvr(symr)

! **********************************************************
! ***   CALCULATE FIRST DERIVATIVE WITH RESPECT TO RHO   ***
! **********************************************************

! Load modules.

  use param
  use arrays, only: diffvar,beta_r
  use procinfo

! Declare variables.

  implicit none

  integer i,j,symr

  real(8) idr,hidr

  real(8), dimension (1-ghost:Nrmax,1-ghost:Nzmax) :: diffadvr


! *******************
! ***   NUMBERS   ***
! *******************

  idr  = 1.d0/dr
  hidr = 0.5d0*idr


! ************************
! ***   SECOND ORDER   ***
! ************************

  if (order=="two") then

     do j=1-ghost,Nz

!       Interior points: second order one-sided derivative.

        do i=3-ghost,Nr-2

           if (beta_r(i,j)<=0.d0) then
              diffadvr(i,j) = + hidr*(3.d0*diffvar(i,j) - 4.d0*diffvar(i-1,j) + diffvar(i-2,j))
           else
              diffadvr(i,j) = - hidr*(3.d0*diffvar(i,j) - 4.d0*diffvar(i+1,j) + diffvar(i+2,j))
           end if

        end do

!       Point i=2-ghost

        i = 2-ghost

        if (beta_r(i,j)<=0.d0) then
           diffadvr(i,j) = hidr*(diffvar(i+1,j) - diffvar(i-1,j))
        else
           diffadvr(i,j) = - hidr*(3.d0*diffvar(i,j) - 4.d0*diffvar(i+1,j) + diffvar(i+2,j))
        end if

!       Point i=1-ghost

        i = 1-ghost

        diffadvr(i,j) = - hidr*(3.d0*diffvar(i,j) - 4.d0*diffvar(i+1,j) + diffvar(i+2,j))

!       Point i=Nr-1

        i = Nr-1

        if (beta_r(i,j)<=0.d0) then
           diffadvr(i,j) = hidr*(3.d0*diffvar(i,j) - 4.d0*diffvar(i-1,j) + diffvar(i-2,j))
        else
           diffadvr(i,j) = hidr*(diffvar(i+1,j) - diffvar(i-1,j))
        end if

!       Point i=Nr

        i = Nr

        diffadvr(i,j) = hidr*(3.d0*diffvar(i,j) - 4.d0*diffvar(i-1,j) + diffvar(i-2,j))

     end do


! ************************
! ***   FOURTH ORDER   ***
! ************************

  else if (order=="four") then

     do j=1-ghost,Nz

!       Interior points: fourth order semi-one-sided derivative.

        do i=4-ghost,Nr-3

           if (beta_r(i,j)<=0.d0) then
              diffadvr(i,j) = + 0.25d0*idr*(3.d0*diffvar(i+1,j) + 10.d0*diffvar(i,j) &
                   - 18.d0*diffvar(i-1,j) + 6.d0*diffvar(i-2,j) - diffvar(i-3,j))/3.d0
           else
              diffadvr(i,j) = - 0.25d0*idr*(3.d0*diffvar(i-1,j) + 10.d0*diffvar(i,j) &
                   - 18.d0*diffvar(i+1,j) + 6.d0*diffvar(i+2,j) - diffvar(i+3,j))/3.d0
           end if

        end do

!       Point i=3-ghost

        i = 3-ghost

        if (beta_r(i,j)<=0.d0) then
           diffadvr(i,j) = 0.25d0*idr*(8.d0*(diffvar(i+1,j) - diffvar(i-1,j)) &
                - (diffvar(i+2,j) - diffvar(i-2,j)))/3.d0
        else
           diffadvr(i,j) = - 0.25d0*idr*(3.d0*diffvar(i-1,j) + 10.d0*diffvar(i,j) &
                - 18.d0*diffvar(i+1,j) + 6.d0*diffvar(i+2,j) - diffvar(i+3,j))/3.d0
        end if

!       Point i=2-ghost

        i = 2-ghost

        diffadvr(i,j) = - 0.25d0*idr*(3.d0*diffvar(i-1,j) + 10.d0*diffvar(i,j) &
             - 18.d0*diffvar(i+1,j) + 6.d0*diffvar(i+2,j) - diffvar(i+3,j))/3.d0

!       Point i=1-ghost

        i = 1-ghost

        diffadvr(i,j) = - 0.25d0*idr*(25.d0*diffvar(i,j) - 48.d0*diffvar(i+1,j) &
             + 36.d0*diffvar(i+2,j) - 16.d0*diffvar(i+3,j) + 3.d0*diffvar(i+4,j))/3.d0

!       Point i=Nr-2

        i = Nr-2

        if (beta_r(i,j)<=0.d0) then
           diffadvr(i,j) = 0.25d0*idr*(3.d0*diffvar(i+1,j) + 10.d0*diffvar(i,j) &
                - 18.d0*diffvar(i-1,j) + 6.d0*diffvar(i-2,j) - diffvar(i-3,j))/3.d0
        else
           diffadvr(i,j) = 0.25d0*idr*(8.d0*(diffvar(i+1,j) - diffvar(i-1,j)) &
                - (diffvar(i+2,j) - diffvar(i-2,j)))/3.d0
        end if

!       Point i=Nr-1

        i = Nr-1

        diffadvr(i,j) = 0.25d0*idr*(3.d0*diffvar(i+1,j) + 10.d0*diffvar(i,j) &
             - 18.d0*diffvar(i-1,j) + 6.d0*diffvar(i-2,j) - diffvar(i-3,j))/3.d0

!       Point i=Nr

        i = Nr

        diffadvr(i,j) = 0.25d0*idr*(25.d0*diffvar(i,j) - 48.d0*diffvar(i-1,j) &
             + 36.d0*diffvar(i-2,j) - 16.d0*diffvar(i-3,j) + 3.d0*diffvar(i-4,j))/3.d0

     end do

  end if


! ******************************
! ***   SYMMETRIES AT AXIS   ***
! ******************************

  if (ownaxis) then
     do i=1,ghost
        diffadvr(1-i,:) = - symr*diffadvr(i,:)
     end do
  end if


! ***********************
! ***   SYNCHRONIZE   ***
! ***********************

! Here I only need to sync on r direction.

  if (nprocr>1) then
     call syncr(diffadvr)
  end if


! ************************
! ***   END FUNCTION   ***
! ************************

  end function diffadvr







  function diffadvz(symz)

! **********************************************************
! ***   CALCULATE FIRST DERIVATIVE WITH RESPECT TO RHO   ***
! **********************************************************

! Load modules.

  use param
  use arrays, only: diffvar,beta_z
  use procinfo

! Declare variables.

  implicit none

  integer i,j,symz

  real(8) idz,hidz

  real(8), dimension (1-ghost:Nrmax,1-ghost:Nzmax) :: diffadvz


! *******************
! ***   NUMBERS   ***
! *******************

  idz  = 1.d0/dz
  hidz = 0.5d0*idz


! ************************
! ***   SECOND ORDER   ***
! ************************

  if (order=="two") then

     do i=1-ghost,Nr

!       Interior points: second order one-sided derivative.

        do j=3-ghost,Nz-2

           if (beta_z(i,j)<=0.d0) then
              diffadvz(i,j) = + hidz*(3.d0*diffvar(i,j) - 4.d0*diffvar(i,j-1) + diffvar(i,j-2))
           else
              diffadvz(i,j) = - hidz*(3.d0*diffvar(i,j) - 4.d0*diffvar(i,j+1) + diffvar(i,j+2))
           end if

        end do

!       Point j=2-ghost

        j = 2-ghost

        if (beta_z(i,j)<=0.d0) then
           diffadvz(i,j) = hidz*(diffvar(i,j+1) - diffvar(i,j-1))
        else
           diffadvz(i,j) = - hidz*(3.d0*diffvar(i,j) - 4.d0*diffvar(i,j+1) + diffvar(i,j+2))
        end if

!       Point j=1-ghost

        j = 1-ghost

        diffadvz(i,j) = - hidz*(3.d0*diffvar(i,j) - 4.d0*diffvar(i,j+1) + diffvar(i,j+2))

!       Point j=Nz-1

        j = Nz-1

        if (beta_z(i,j)<=0.d0) then
           diffadvz(i,j) = hidz*(3.d0*diffvar(i,j) - 4.d0*diffvar(i,j-1) + diffvar(i,j-2))
        else
           diffadvz(i,j) = hidz*(diffvar(i,j+1) - diffvar(i,j-1))
        end if

!       Point j=Nz

        j = Nz

        diffadvz(i,j) = hidz*(3.d0*diffvar(i,j) - 4.d0*diffvar(i,j-1) + diffvar(i,j-2))

     end do


! ************************
! ***   FOURTH ORDER   ***
! ************************

  else if (order=="four") then

     do i=1-ghost,Nr

!       Interior points: fourth order semi-one-sided derivative.

        do j=4-ghost,Nz-3

           if (beta_z(i,j)<=0.d0) then
              diffadvz(i,j) = + 0.25d0*idz*(3.d0*diffvar(i,j+1) + 10.d0*diffvar(i,j) &
                   - 18.d0*diffvar(i,j-1) + 6.d0*diffvar(i,j-2) - diffvar(i,j-3))/3.d0
           else
              diffadvz(i,j) = - 0.25d0*idz*(3.d0*diffvar(i,j-1) + 10.d0*diffvar(i,j) &
                   - 18.d0*diffvar(i,j+1) + 6.d0*diffvar(i,j+2) - diffvar(i,j+3))/3.d0
           end if

        end do

!       Point j=3-ghost

        j = 3-ghost

        if (beta_z(i,j)<=0.d0) then
           diffadvz(i,j) = - 0.25d0*idz*(3.d0*diffvar(i,j-1) + 10.d0*diffvar(i,j) &
                - 18.d0*diffvar(i,j+1) + 6.d0*diffvar(i,j+2) - diffvar(i,j+3))/3.d0
        else
           diffadvz(i,j) = 0.25d0*idz*(8.d0*(diffvar(i,j+1) - diffvar(i,j-1)) &
                - (diffvar(i,j+2) - diffvar(i,j-2)))/3.d0
        end if

!       Point j=2-ghost

        j = 2-ghost

        diffadvz(i,j) = - 0.25d0*idz*(3.d0*diffvar(i,j-1) + 10.d0*diffvar(i,j) &
             - 18.d0*diffvar(i,j+1) + 6.d0*diffvar(i,j+2) - diffvar(i,j+3))/3.d0

!       Point j=1-ghost

        j = 1-ghost

        diffadvz(i,j) = - 0.25d0*idz*(25.d0*diffvar(i,j) - 48.d0*diffvar(i,j+1) &
             + 36.d0*diffvar(i,j+2) - 16.d0*diffvar(i,j+3) + 3.d0*diffvar(i,j+4))/3.d0

!       Point j=Nz-2

        j = Nz-2

        if (beta_z(i,j)<=0.d0) then
           diffadvz(i,j) = 0.25d0*idz*(3.d0*diffvar(i,j+1) + 10.d0*diffvar(i,j) &
                - 18.d0*diffvar(i,j-1) + 6.d0*diffvar(i,j-2) - diffvar(i,j-3))/3.d0
        else
           diffadvz(i,j) = 0.25d0*idz*(8.d0*(diffvar(i,j+1) - diffvar(i,j-1)) &
                - (diffvar(i,j+2) - diffvar(i,j-2)))/3.d0
        end if

!       Point j=Nz-1

        j = Nz-1

        diffadvz(i,j) = 0.25d0*idz*(3.d0*diffvar(i,j+1) + 10.d0*diffvar(i,j) &
             - 18.d0*diffvar(i,j-1) + 6.d0*diffvar(i,j-2) - diffvar(i,j-3))/3.d0

!       Point j=Nz

        j = Nz

        diffadvz(i,j) = 0.25d0*idz*(25.d0*diffvar(i,j) - 48.d0*diffvar(i,j-1) &
             + 36.d0*diffvar(i,j-2) - 16.d0*diffvar(i,j-3) + 3.d0*diffvar(i,j-4))/3.d0

     end do

  end if


! *********************************
! ***   SYMMETRIES ON EQUATOR   ***
! *********************************

  if (eqsym.and.ownequator) then
     do j=1,ghost
        diffadvz(:,1-j) = - symz*diffadvz(:,j)
     end do
  end if


! ***********************
! ***   SYNCHRONIZE   ***
! ***********************

! Here I only need to sync on z direction.

  if (nprocz>1) then
     call syncz(diffadvz)
  end if


! ************************
! ***   END FUNCTION   ***
! ************************

  end function diffadvz







! **********************
! ***   END MODULE   ***
! **********************

  end module derivadvect
