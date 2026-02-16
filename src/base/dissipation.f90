
  subroutine dissipation(symr,symz,diss)

! This routine adds Kreiss-Oliger dissipation to the sources.
!
! The idea behind Kreiss-Oliger dissipation is to add a
! dissipative term to the sources that will damp out
! the high frequency modes, which are the ones that
! typically go unstable.  The dissipative term has the
! form of an even difference operator that is two
! orders higher than the order of the numerical scheme.
!
! Why even?  Because only even difference operators are
! always dissipative (for small coefficients).  Odd
! difference operators dissipate modes propagating in a
! given direction while amplifying modes propagating in the
! oppostive direction.
!
! Why two orders above the order of the scheme? To ensure
! that the convergence order of the scheme is not spoiled.
! For example, for a second order scheme in 1 dimension "x"
! we add to the source of a function f a term of the form:
!
! df/dt  ->  df/dt  +  coeff ( Dx^3 d^4(f)/dx^4 )
!
! where Dx is the spatial interval. If we now approximate
! the fourth derivative as:
!
! d^4(f)/dx^4  =  D^4 f / Dx^4
!
! with D^4 the fourth order difference operator, we end up with:
!
! df/dt  ->  df/dt  +  (1/Dx) D^4(f)
!
! and when we do the finite differencing in time this becomes
! essentially:
!
! f(t+Dt) ~ f(t) + (Dt/Dx) D^4 f
!
! That is, the term we added to the source is actually third order,
! so it does not spoil the second order convergence. But at the
! finite difference level it amounts to adding a term that
! is always the same size as long at Dt/dx is constant.
!
! When we have a fourth order scheme we actually add D^6(f)/Dx
! to the sources.  In 2D we also add mixed terms with even
! difference operators that have the required order.
! So, for the second order scheme we add terms of the form D^4_r,
! D^4_z, and D^2_r D^2_z,  while for the fourth order scheme we add
! D^6_r, D^6_z and D^2_r D^4_z+D^2^z D^4_r.
!
! A final word about the signs:  It turns out that each extra D^2
! term in a given direction changes the sign from dissipation to
! amplification. The reason for this is easier to understand if
! we expand the solution in Fourier modes of the form f(x) = exp(ikx).
! Then each two derivatives bring a factor i^2=-1.
!
! This implies that a D^2 term must de ADDED, a D^4 term must be
! SUBTRACTED, and a D^6 term must be ADDED.
! For mixed terms, a term D^2_r D^2_z must be ADDED, while a term
! D^2_r D^4_z (or viceversa) must be SUBTRACTED.


! **************************************
! ***   ADD DISSIPATION TO SOURCES   ***
! **************************************

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  integer i,j
  integer symr,symz

  real(8) idt,diss


! *******************
! ***   NUMBERS   ***
! *******************

  idt = 1.d0/dt


! ************************
! ***   SANITY CHECK   ***
! ************************

  if ((abs(symr)/=1).or.(abs(symz)/=1)) then
     print *
     print *, 'dissipation: Symmetry parameters must be +-1'
     print *
     call die
  end if


! ************************
! ***   SECOND ORDER   ***
! ************************

! Here we add/subtract even fourth order operators.
!
! The fourth order difference in a given direction is:
!
! D^4 f  = 6 f  -  4 ( f    +  f   )  +  ( f    +  f   )
!      i      i         i+1     i-1         i+2     i-2
!
! The mixed operator has the form:
!
! D^2_r D^2_z f    =  4 f   -  2 ( f     +  f     +  f      +  f     )
!              i,j       i,j        i+1,j    i-1,j    i,j+1     i,j-1
!
!                  +  ( f        +  f        +  f        +  f       )
!                        i+1,j+1     i+1,j-1     i-1,j+1     i-1,j-1

  if (order=="two") then

!    Loop in r direction subtracting the operator:  D^4_r.

     do i=3-ghost,Nr-2
        sourcevar(i,:) = sourcevar(i,:) - diss*idt*(6.d0*evolvevar(i,:) &
                       - 4.d0*(evolvevar(i+1,:) + evolvevar(i-1,:)) &
                       +      (evolvevar(i+2,:) + evolvevar(i-2,:)))/6.d0
     end do

!    Loop in z direction subtracting the operator:  D^4_z.

     do j=3-ghost,Nz-2
        sourcevar(:,j) = sourcevar(:,j) - diss*idt*(6.d0*evolvevar(:,j) &
                       - 4.d0*(evolvevar(:,j+1) + evolvevar(:,j-1)) &
                       +      (evolvevar(:,j+2) + evolvevar(:,j-2)))/6.d0
     end do

!    Now loop in both directions adding the mixed operator: D^2_r D^2_z.

     !do j=3-ghost,Nz-2
     !   do i=3-ghost,Nr-2
     !      sourcevar(i,j) = sourcevar(i,j) + diss*idt*(4.d0*evolvevar(i,j) &
     !                     - 2.d0*(evolvevar(i+1,j  ) + evolvevar(i-1,j  ) &
     !                     +       evolvevar(i  ,j+1) + evolvevar(i  ,j-1)) &
     !                     +      (evolvevar(i+1,j+1) + evolvevar(i+1,j-1) &
     !                     +       evolvevar(i-1,j+1) + evolvevar(i-1,j-1)))/4.d0
     !   end do
     !end do


! ************************
! ***   FOURTH ORDER   ***
! ************************

! Here we add/subtract even sixth order operators.
!
! The sixth order difference in a given direction is:
!
! D^6 f  = - 20 f  +  15 ( f    +  f   )  -  6 ( f    +  f   )  +  ( f    +  f   )
!      i         i          i+1     i-1           i+2     i-2         i+3     i-3
!
! The mixed operator has the form:
!
! D^2_z D^4_r + D^2_r D^4_z  =  - 24 f     +  14 ( f     +  f     +  f      +  f     )
!                                     i,j           i+1,j    i-1,j    i,j+1     i,j-1
!
!                            -  8 ( f        +  f        +  f        +  f       )
!                                    i+1,j+1     i+1,j-1     i-1,j+1     i-1,j-1
!
!                            -  2 ( f        +  f        +  f        +  f       )
!                                    i+2,j       i-2,j       i,j+2       i,j-2
!
!                            +    ( f        +  f        +  f        +  f    
!                                    i+2,j+1     i+2,j-1     i-2,j+1     i-2,j-1
!
!                                +  f        +  f        +  f        +  f       )
!                                    i+1,j+2     i+1,j-2     i-1,j+2     i-1,j-2


  else if (order=="four") then

!    Loop in r direction adding the operator:  D^6_r.

     do i=4-ghost,Nr-3
        sourcevar(i,:) = sourcevar(i,:) - diss*idt*(20.d0*evolvevar(i,:) &
                       - 15.d0*(evolvevar(i+1,:) + evolvevar(i-1,:)) &
                       +  6.d0*(evolvevar(i+2,:) + evolvevar(i-2,:)) &
                       -       (evolvevar(i+3,:) + evolvevar(i-3,:)))/20.d0
     end do

!    Loop in z direction adding the operator:  D^6_z.

     do j=4-ghost,Nz-3
        sourcevar(:,j) = sourcevar(:,j) - diss*idt*(20.d0*evolvevar(:,j) &
                       - 15.d0*(evolvevar(:,j+1) + evolvevar(:,j-1)) &
                       +  6.d0*(evolvevar(:,j+2) + evolvevar(:,j-2)) &
                       -       (evolvevar(:,j+3) + evolvevar(:,j-3)))/20.d0
     end do

!    Now loop in both directions subtracting the mixed operator:  D^2_z D^4_ r + D^2_r D^4_z.

     !do j=4-ghost,Nz-3
     !   do i=4-ghost,Nr-3
     !      sourcevar(i,j) = sourcevar(i,j) - diss*idt*(-24.d0*evolvevar(i,j) &
     !                     + 14.d0*(evolvevar(i+1,j  ) + evolvevar(i-1,j  )  &
     !                     +        evolvevar(i  ,j+1) + evolvevar(i  ,j-1)) &
     !                     -  8.d0*(evolvevar(i+1,j+1) + evolvevar(i-1,j+1)  &
     !                     +        evolvevar(i+1,j-1) + evolvevar(i-1,j-1)) &
     !                     -  2.d0*(evolvevar(i+2,j  ) + evolvevar(i-2,j  )  &
     !                     +        evolvevar(i  ,j+2) + evolvevar(i  ,j-2)) &
     !                     +       (evolvevar(i+2,j+1) + evolvevar(i-2,j+1)  &
     !                     +        evolvevar(i+2,j-1) + evolvevar(i-2,j-1)  &
     !                     +        evolvevar(i+1,j+2) + evolvevar(i-1,j+2)  &
     !                     +        evolvevar(i+1,j-2) + evolvevar(i-1,j-2)))/24.d0
     !   end do
     !end do

  end if


! ******************************************
! ***   SYMMETRIES ON AXIS AND EQUATOR   ***
! ******************************************

! Symmetries in z for equatorial symmetry.
! Must be done first to ensure that the corners
! end up with the correct values.

  if (eqsym.and.ownequator) then
     do j=1,ghost
        sourcevar(:,1-j) = symz*sourcevar(:,j)
     end do
  end if

! Symmetries in r.

  if (ownaxis) then
     do i=1,ghost
        sourcevar(1-i,:) = symr*sourcevar(i,:)
     end do
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine dissipation
