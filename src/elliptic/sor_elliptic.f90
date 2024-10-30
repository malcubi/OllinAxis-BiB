!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/elliptic/sor_elliptic.f90,v 1.13 2021/02/24 22:35:01 malcubi Exp $

  subroutine sor_elliptic(type,init)

! *******************************
! ***   SOR ELLIPTIC SOLVER   ***
! *******************************

! This is a simple SOR elliptic solver for the Poisson equation:
!
!  __2                   5
!  \/ u  +  S1 u  +  S5 u  =  S0
!
! Here S0 is a source term, S1 is the coefficient of the 
! linear term, and S5 the coefficient of the u^5 term
! (this term often appears in the Hamiltonian constraint).
!
! Notice that at the moment this routine only works for
! a flat Laplacian.
!
! This routine is just a driver for the SOR routine from
! Numerical Recipes "sor" (below).  Notice that the SOR
! solver is only second order accurate!
!
! Also, in order to modify as little as possible the NR
! routine, the arrays (a,b,c,d,e,f,u) go from 1 to N+ghost,
! in contrast with the standard arrays from the code that go
! from 1-ghost to N.  We need to be careful with this when
! copying arrays.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical flag1,flag2        ! Interpolation flags.

  integer i,j,n,m            ! Counters.
  integer imax,jmax          ! Array size.
  integer box,level          ! Box number and level counters.

  real(8) rjac               ! Spectral radius of Jacobi iteration.
  real(8) r0,z0,interp       ! For interpolation.
  real(8) zero,half,one,two,smallpi
  real(8) aux1,aux2

  real(8), dimension(1:Nrmax+ghost,1:Nzmax+ghost) :: sor_a,sor_b,sor_c,sor_d,sor_e,sor_f,sor_u
  real(8), dimension(1:Nrmax+ghost,1:Nzmax+ghost) :: sor_g,sor_h,sor_i,sor_j

  character(*) type          ! Type of Laplacian (flat,conformal,physical).
  character(*) init          ! Initial guess.


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  half = 0.5d0
  one  = 1.d0
  two  = 2.d0

  smallpi = acos(-1.d0)

! Warning about second order.

  if (order=='four') then
     if (rank==0) then
        print *, 'WARNING: The SOR solver is only second order accurate!'
     end if
  end if


! *******************************
! ***   INITIALIZE SOLUTION   ***
! *******************************

! By default we initialize ell_u to 1, but we can
! overrdide this depending on the value of 'init'.

  do box=0,Nb
     do level=min(1,box),Nl(box)

        call currentgrid(box,level,grid(box,level))

        if (init=='zero') then
           ell_u = zero
        else if (init=='one') then
           ell_u = one
        else
           if (rank==0) then
              print *
              print *, 'Unknow initial data for sor_elliptic'
              print *, 'Setting to 1.'
              print *
           end if
           ell_u = one
        end if

     end do
  end do

! Initialize arrays for coefficients of elliptic equation.
! We set the source (f) and all off-diagonal coefficients
! (a,b,c,d,g,h,i,j) to zero, and we set the diagonal
! coefficient (e) to 1.

  sor_f = zero

  sor_e = one
  sor_a = zero
  sor_b = zero
  sor_c = zero
  sor_d = zero
  sor_g = zero
  sor_h = zero
  sor_i = zero
  sor_j = zero


! *************************************
! ***   FIRST SOLVE ON COARSE GRID  ***
! *************************************

! Jacobi spectral radius (for flat Laplacian).

  rjac = (dz/dr*cos(smallpi/dble(Nr+ghost)) &
       +  dr/dz*cos(smallpi/dble(Nz+ghost))) &
       / (dz/dr + dr/dz)

! Point to coarse grid.

  call currentgrid(0,0,grid(0,0))

! Loop over grid to find coefficients of
! elliptic equation.

  imax = Nr+ghost
  jmax = Nz+ghost

  do i=1-ghost,Nr
     do j=1-ghost,Nz

!       Initial guess.

        sor_u(i+ghost,j+ghost) = ell_u(i,j)

!       Coefficients of equation assuming a flat Laplacian.

        if (type=='flat') then

           sor_a(i+ghost,j+ghost) = (one/dr + half/r(i,j))/dr
           sor_b(i+ghost,j+ghost) = (one/dr - half/r(i,j))/dr

           sor_c(i+ghost,j+ghost) = one/dz**2
           sor_d(i+ghost,j+ghost) = one/dz**2

           sor_e(i+ghost,j+ghost) = - two*(one/dr**2 + one/dz**2) + ell_S1(i,j)

!       Physical Laplacian.

        else if ((type=='physical').or.(type=='conformal')) then

!          Terms that exist in flat Laplacian (with conformal metric).

           sor_a(i+ghost,j+ghost) = g_A(i,j)*(one/dr + half/r(i,j))/dr
           sor_b(i+ghost,j+ghost) = g_A(i,j)*(one/dr - half/r(i,j))/dr

           sor_c(i+ghost,j+ghost) = g_B(i,j)/dz**2
           sor_d(i+ghost,j+ghost) = g_B(i,j)/dz**2

           sor_e(i+ghost,j+ghost) = - two*(g_A(i,j)/dr**2 + g_B(i,j)/dz**2)

!          Mixed second derivative.

           sor_g(i+ghost,j+ghost) = + 0.5d0*r(i,j)*g_C(i,j)/(dr*dz)
           sor_h(i+ghost,j+ghost) = + 0.5d0*r(i,j)*g_C(i,j)/(dr*dz)
           sor_i(i+ghost,j+ghost) = - 0.5d0*r(i,j)*g_C(i,j)/(dr*dz)
           sor_j(i+ghost,j+ghost) = - 0.5d0*r(i,j)*g_C(i,j)/(dr*dz)

!          Derivatives of conformal metric.

           sor_a(i+ghost,j+ghost) = sor_a(i+ghost,j+ghost) &
                                  + half*(Dr_g_A(i,j) + r(i,j)*Dz_g_C(i,j))/dr
           sor_b(i+ghost,j+ghost) = sor_b(i+ghost,j+ghost) &
                                  - half*(Dr_g_A(i,j) + r(i,j)*Dz_g_C(i,j))/dr

           sor_c(i+ghost,j+ghost) = sor_c(i+ghost,j+ghost) &
                                  + half*(Dz_g_B(i,j) + r(i,j)*Dr_g_C(i,j) + g_C(i,j))/dz
           sor_d(i+ghost,j+ghost) = sor_d(i+ghost,j+ghost) &
                                  - half*(Dz_g_B(i,j) + r(i,j)*Dr_g_C(i,j) + g_C(i,j))/dz

!          Derivatives of determinant of conformal metric.

           sor_a(i+ghost,j+ghost) = sor_a(i+ghost,j+ghost) + 0.25d0*ihdet(i,j) &
                                  *(g_A(i,j)*Dr_hdet(i,j) + r(i,j)*g_C(i,j)*Dz_hdet(i,j))/dr
           sor_b(i+ghost,j+ghost) = sor_b(i+ghost,j+ghost) - 0.25d0*ihdet(i,j) &
                                  *(g_A(i,j)*Dr_hdet(i,j) + r(i,j)*g_C(i,j)*Dz_hdet(i,j))/dr

           sor_c(i+ghost,j+ghost) = sor_c(i+ghost,j+ghost) + 0.5d0*g_C(i,j)/dz + 0.25d0*ihdet(i,j) &
                                  *(g_B(i,j)*Dz_hdet(i,j) + r(i,j)*g_C(i,j)*Dr_hdet(i,j))/dz
           sor_d(i+ghost,j+ghost) = sor_d(i+ghost,j+ghost) - 0.5d0*g_C(i,j)/dz - 0.25d0*ihdet(i,j) &
                                  *(g_B(i,j)*Dz_hdet(i,j) + r(i,j)*g_C(i,j)*Dr_hdet(i,j))/dz

!          Derivatives of conformal factor.

           if (type/='conformal') then

              sor_a(i+ghost,j+ghost) = sor_a(i+ghost,j+ghost) &
                                     + (g_A(i,j)*Dr_phi(i,j) + r(i,j)*g_C(i,j)*Dz_phi(i,j))/dr
              sor_b(i+ghost,j+ghost) = sor_b(i+ghost,j+ghost) &
                                     - (g_A(i,j)*Dr_phi(i,j) + r(i,j)*g_C(i,j)*Dz_phi(i,j))/dr

              sor_c(i+ghost,j+ghost) = sor_c(i+ghost,j+ghost) &
                                     + (g_B(i,j)*Dz_phi(i,j) + r(i,j)*g_C(i,j)*Dr_phi(i,j))/dz
              sor_d(i+ghost,j+ghost) = sor_d(i+ghost,j+ghost) &
                                     - (g_B(i,j)*Dz_phi(i,j) + r(i,j)*g_C(i,j)*Dr_phi(i,j))/dz

           end if

!          Rescale by conformal factor since all the above terms
!          were calculated with the conformal metric.

           if (type/='conformal') then

              sor_a(i+ghost,j+ghost) = sor_a(i+ghost,j+ghost)/psi4(i,j)
              sor_b(i+ghost,j+ghost) = sor_b(i+ghost,j+ghost)/psi4(i,j)
              sor_c(i+ghost,j+ghost) = sor_c(i+ghost,j+ghost)/psi4(i,j)
              sor_d(i+ghost,j+ghost) = sor_d(i+ghost,j+ghost)/psi4(i,j)

              sor_e(i+ghost,j+ghost) = sor_e(i+ghost,j+ghost)/psi4(i,j)

              sor_g(i+ghost,j+ghost) = sor_g(i+ghost,j+ghost)/psi4(i,j)
              sor_h(i+ghost,j+ghost) = sor_h(i+ghost,j+ghost)/psi4(i,j)
              sor_i(i+ghost,j+ghost) = sor_i(i+ghost,j+ghost)/psi4(i,j)
              sor_j(i+ghost,j+ghost) = sor_j(i+ghost,j+ghost)/psi4(i,j)

           end if

!          Add ell_S1 to sor_e.

           sor_e(i+ghost,j+ghost) = sor_e(i+ghost,j+ghost) + ell_S1(i,j)

!       Unknown Laplacian type.

        else

           if (rank==0) then
              print *
              print *, 'Unknow Laplacian type'
              print *, 'Aborting! (subroutine sor_elliptic.f90)'
              print *
           end if

           call die

        end if

!       Source term.

        sor_f(i+ghost,j+ghost) = ell_S0(i,j)

     end do
  end do

! Call SOR.

  call sor(sor_a,sor_b,sor_c,sor_d,sor_e,sor_f,sor_g,sor_h,sor_i,sor_j,sor_u,imax,jmax,rjac,1.d0,'robin')

! Copy SOR solution into ell_u.

  do i=1-ghost,Nr
     do j=1-ghost,Nz
        ell_u(i,j) = sor_u(i+ghost,j+ghost)
     end do
  end do


! **********************
! ***   FINE GRIDS   ***
! **********************

! If we have refinement levels solve on finer grids.

  if (Nlmax>0) then

!    Iterate over boxes and levels higher than 0.

     do box=0,Nb
        do level=1,Nl(box)

!          If this level does not exist for this box cycle.

           if (Nl(box)<level) cycle

!          Point to current grid.

           call currentgrid(box,level,grid(box,level))

!          First we interpolate from coarser grid to find
!          initial guess and Dirichlet boundary conditions.         

!          Loop over ALL points in the current grid (and I do mean
!          all of them, not just the ones in the current processor).

           do m=0,Nrbox(box)+ghost-1
              do n=0,Nzbox(box)+ghost-1

!                Figure out (r,z) position for interpolation.

                 r0 = rminl(box,level) + dble(m)*drl(level)
                 z0 = zminl(box,level) + dble(n)*dzl(level)

!                Figure out to which grid point this (r0,z0) values
!                would correspond in the local processor at the
!                fine grid level.

                 i = nint((r0-r(1-ghost,0))/drl(level)) + 1 - ghost
                 j = nint((z0-z(0,1-ghost))/dzl(level)) + 1 - ghost

!                Notice now that the (i,j) values above might be outside
!                the range of the current processor.  This means that
!                we should not try to access this location as it belongs
!                to another processor.

                 if ((i>=1-ghost).and.(i<=Nrl(box,rank)).and.(j>=1-ghost).and.(j<=Nzl(box,rank))) then
                    flag1 = .true.
                 else
                    flag1 = .false.
                 end if

!                Interpolate variable from coarser grid level.

                 if (level==1) then
                    interpvar => grid(0,0)%ell_u
                    aux1 = interp(0,0,r0,z0,flag2)
                 else
                    interpvar => grid(box,level-1)%ell_u
                    aux1 = interp(box,level-1,r0,z0,flag2)
                 end if

                 call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!                But we only update the value if the location belongs to us.

                 if (flag1) then
                    ell_u(i,j) = aux2
                 end if

              end do
           end do

!          Now loop over current grid to find coefficients
!          of elliptic equation.

           imax = Nr+ghost
           jmax = Nz+ghost

           do i=1-ghost,Nr
              do j=1-ghost,Nz

!                Initial guess.

                 sor_u(i+ghost,j+ghost) = ell_u(i,j)

!                Coefficients of equation assuming a flat Laplacian.

                 if (type=='flat') then

                    sor_a(i+ghost,j+ghost) = (one/dr + half/r(i,j))/dr
                    sor_b(i+ghost,j+ghost) = (one/dr - half/r(i,j))/dr

                    sor_c(i+ghost,j+ghost) = one/dz**2
                    sor_d(i+ghost,j+ghost) = one/dz**2

                    sor_e(i+ghost,j+ghost) = - two*(one/dr**2 + one/dz**2) + ell_S1(i,j)

!                Physical Laplacian.

                 else if (type=='physical') then

!                   Terms that exist in flat Laplacian (with conformal metric).

                    sor_a(i+ghost,j+ghost) = g_A(i,j)*(one/dr + half/r(i,j))/dr
                    sor_b(i+ghost,j+ghost) = g_A(i,j)*(one/dr - half/r(i,j))/dr

                    sor_c(i+ghost,j+ghost) = g_B(i,j)/dz**2
                    sor_d(i+ghost,j+ghost) = g_B(i,j)/dz**2

                    sor_e(i+ghost,j+ghost) = - two*(g_A(i,j)/dr**2 + g_B(i,j)/dz**2)

!                   Mixed second derivative.

                    sor_g(i+ghost,j+ghost) = + 0.5d0*r(i,j)*g_C(i,j)/dr/dz
                    sor_h(i+ghost,j+ghost) = + 0.5d0*r(i,j)*g_C(i,j)/dr/dz
                    sor_i(i+ghost,j+ghost) = - 0.5d0*r(i,j)*g_C(i,j)/dr/dz
                    sor_j(i+ghost,j+ghost) = - 0.5d0*r(i,j)*g_C(i,j)/dr/dz

!                   Derivatives of conformal metric.

                    sor_a(i+ghost,j+ghost) = sor_a(i+ghost,j+ghost) &
                                           + half*(Dr_g_A(i,j) + r(i,j)*Dz_g_C(i,j))/dr
                    sor_b(i+ghost,j+ghost) = sor_b(i+ghost,j+ghost) &
                                           - half*(Dr_g_A(i,j) + r(i,j)*Dz_g_C(i,j))/dr

                    sor_c(i+ghost,j+ghost) = sor_c(i+ghost,j+ghost) &
                                           + half*(Dz_g_B(i,j) + r(i,j)*Dr_g_C(i,j) + g_C(i,j))/dz
                    sor_d(i+ghost,j+ghost) = sor_d(i+ghost,j+ghost) &
                                           - half*(Dz_g_B(i,j) + r(i,j)*Dr_g_C(i,j) + g_C(i,j))/dz

!                   Derivatives of determinant of conformal metric.

                    sor_a(i+ghost,j+ghost) = sor_a(i+ghost,j+ghost) + 0.25d0*ihdet(i,j) &
                                  *(g_A(i,j)*Dr_hdet(i,j) + r(i,j)*g_C(i,j)*Dz_hdet(i,j))/dr
                    sor_b(i+ghost,j+ghost) = sor_b(i+ghost,j+ghost) - 0.25d0*ihdet(i,j) &
                                  *(g_A(i,j)*Dr_hdet(i,j) + r(i,j)*g_C(i,j)*Dz_hdet(i,j))/dr

                    sor_c(i+ghost,j+ghost) = sor_c(i+ghost,j+ghost) + 0.5d0*g_C(i,j)/dz + 0.25d0*ihdet(i,j) &
                                  *(g_B(i,j)*Dz_hdet(i,j) + r(i,j)*g_C(i,j)*Dr_hdet(i,j))/dz
                    sor_d(i+ghost,j+ghost) = sor_d(i+ghost,j+ghost) - 0.5d0*g_C(i,j)/dz - 0.25d0*ihdet(i,j) &
                                  *(g_B(i,j)*Dz_hdet(i,j) + r(i,j)*g_C(i,j)*Dr_hdet(i,j))/dz

!                   Derivatives of conformal factor.

                    if (type/='conformal') then

                       sor_a(i+ghost,j+ghost) = sor_a(i+ghost,j+ghost) &
                                     + (g_A(i,j)*Dr_phi(i,j) + r(i,j)*g_C(i,j)*Dz_phi(i,j))/dr
                       sor_b(i+ghost,j+ghost) = sor_b(i+ghost,j+ghost) &
                                     - (g_A(i,j)*Dr_phi(i,j) + r(i,j)*g_C(i,j)*Dz_phi(i,j))/dr

                       sor_c(i+ghost,j+ghost) = sor_c(i+ghost,j+ghost) &
                                     + (g_B(i,j)*Dz_phi(i,j) + r(i,j)*g_C(i,j)*Dr_phi(i,j))/dz
                       sor_d(i+ghost,j+ghost) = sor_d(i+ghost,j+ghost) &
                                     - (g_B(i,j)*Dz_phi(i,j) + r(i,j)*g_C(i,j)*Dr_phi(i,j))/dz

                    end if

!                   Rescale by conformal factor since all the above terms
!                   were calculated with the conformal metric.

                    if (type/='conformal') then

                       sor_a(i+ghost,j+ghost) = sor_a(i+ghost,j+ghost)/psi4(i,j)
                       sor_b(i+ghost,j+ghost) = sor_b(i+ghost,j+ghost)/psi4(i,j)
                       sor_c(i+ghost,j+ghost) = sor_c(i+ghost,j+ghost)/psi4(i,j)
                       sor_d(i+ghost,j+ghost) = sor_d(i+ghost,j+ghost)/psi4(i,j)

                       sor_e(i+ghost,j+ghost) = sor_e(i+ghost,j+ghost)/psi4(i,j)

                       sor_g(i+ghost,j+ghost) = sor_g(i+ghost,j+ghost)/psi4(i,j)
                       sor_h(i+ghost,j+ghost) = sor_h(i+ghost,j+ghost)/psi4(i,j)
                       sor_i(i+ghost,j+ghost) = sor_i(i+ghost,j+ghost)/psi4(i,j)
                       sor_j(i+ghost,j+ghost) = sor_j(i+ghost,j+ghost)/psi4(i,j)

                    end if

!                   Add ell_S1 to sor_e.

                    sor_e(i+ghost,j+ghost) = sor_e(i+ghost,j+ghost) + ell_S1(i,j)

                 end if

!                Source term.

                 sor_f(i+ghost,j+ghost) = ell_S0(i,j)

              end do
           end do

!          Call SOR.

           call sor(sor_a,sor_b,sor_c,sor_d,sor_e,sor_f,sor_g,sor_h,sor_i,sor_j,sor_u,imax,jmax,rjac,1.d0,'dirichlet')

!          Copy SOR solution into ell_u.

           do i=1-ghost,Nr
              do j=1-ghost,Nz
                 ell_u(i,j) = sor_u(i+ghost,j+ghost)
              end do
           end do

        end do
     end do

  end if


! ****************************************************
! ***   RESTRICT FINE GRID DATA INTO COARSE GRID   ***
! ****************************************************

  if (Nlmax>0) then

     do box=0,Nb
        do level=1,Nl(box)

!          If this level does not exist for this box cycle.

           if (Nl(box)<level) cycle

!          Restriction for ell_u.

           finevar => grid(box,level)%ell_u

           if (level==1) then
              coarsevar => grid(0,level-1)%ell_u
           else
              coarsevar => grid(box,level-1)%ell_u
           end if

           call restrict(box,level,.false.)

        end do
     end do
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sor_elliptic









  subroutine sor(sor_a,sor_b,sor_c,sor_d,sor_e,sor_f,sor_g,sor_h,sor_i,sor_j,sor_u,imax,jmax,rjac,uinf,bound)

! **********************************************
! ***   SOR ROUTINE FROM NUMERICAL RECIPES   ***
! **********************************************

! Succesive overrelaxation with Chebyshev acceleration (from Numerical Recipes).
!
! The arrays (a,b,c,d,e,f) are input as coefficients of the equation, with
! dimensions (imax,jmax).  On input 'u' has the initial guess, and on output
! it has the solution. The quantity 'rjac' is the spectral radius of the
! Jacobi iteration.
!
! The finite difference equation is written as:
!
! a(i,j) u(i+1,l) + b(i,j) u(i-1,l) + c(i,j) u(i,j+1) + d(i,j) u(i,j-1) + e(i,j) u(i,j) = f(i,j)
!
! Notice that this implies that the solution is only second order accurate!
!
! The routine is taken directly from Numerical Recipes, and slightly modified
! to allow for Robin boundary conditions and to run in parallel.
!
! The routine also allows for a non-linear term of the type S5*u^5 but cheating.
! I write it as (S5*u^4)*u, and add the term on parenthesis to the linear term.
! This seems to work for the cases I have tried, but I can't guarantee it will
! always converge.
!
! Notice tha in order to modify as little as possible the original NR routine,
! the arrays (a,b,c,d,e,f,u) go from 1 to N+ghost, in contrast with the standard
! arrays from the code that go from 1-ghost to N.  We need to be careful with
! this when copying arrays.

! NOTE:  FOR THE MOMENT WORKING ON THE COARSE GRID ONLY (NO REFINEMENTS)!!

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer i,j,n
  integer isw,jsw,ipass
  integer imax,jmax

  real(8) rjac,uinf
  real(8) anorm,ganorm
  real(8) omega,resid
  real(8) cr,cz

  real(8), dimension(1:Nrmax+ghost,1:Nzmax+ghost) :: sor_a,sor_b,sor_c,sor_d,sor_e,sor_f,sor_u,sor_aux
  real(8), dimension(1:Nrmax+ghost,1:Nzmax+ghost) :: sor_g,sor_h,sor_i,sor_j

  character(*) bound


! ************************
! ***   SANITY CHECK   ***
! ************************

! Only Robin and Dirichlet boundary conditions allowd.
! For Dirichlet conditions, the boundary values should
! be set before entering this routine, as there are not
! touched here.

  if ((bound/='robin').and.(bound/='dirichlet')) then
     print *, 'Only Dirichlet and Robin boundary conditions allowed.'
     print *, 'Switching to Dirichlet (subroutine sor_elliptic.f90)'
     bound = 'dirichlet'
     print *
  end if


! *****************
! ***   START   ***
! *****************

! Initialize omega.

  omega = 1.d0

! Copy coefficient of linear term into sor_aux
! (to deal with non-linear term).

  sor_aux = sor_e

! Start iterations.

  do n=1,ELL_maxiter

     anorm = 0.d0
     isw = 1

!    For the nonlinear term S5*u^5, we first modify the coefficient
!    of the linear term to contain all non-linear terms.

     do i=2,imax-1
        do j=2,jmax-1
           sor_e(i,j) = sor_aux(i,j) + ell_S5(i-ghost,j-ghost)*sor_u(i,j)**4
        end do
     end do

!    Odd-even ordering for updating the solution.

     do ipass=1,2

        jsw = isw

        do i=2,imax-1

           do j=jsw+1,jmax-1,2
              resid = sor_a(i,j)*sor_u(i+1,j) + sor_b(i,j)*sor_u(i-1,j) &
                    + sor_c(i,j)*sor_u(i,j+1) + sor_d(i,j)*sor_u(i,j-1) &
                    + sor_g(i,j)*sor_u(i+1,j+1) + sor_h(i,j)*sor_u(i-1,j-1) &
                    + sor_i(i,j)*sor_u(i+1,j-1) + sor_j(i,j)*sor_u(i-1,j+1) &
                    + sor_e(i,j)*sor_u(i,j) - sor_f(i,j)
              anorm = anorm + abs(resid)
              sor_u(i,j) = sor_u(i,j) - omega*resid/sor_e(i,j)
           end do

           jsw = 3 - jsw

        end do

        isw = 3 - isw

!       Overrelaxation factor.

        if ((n.eq.1).and.(ipass.eq.1)) then
           omega = 1.d0/(1.d0 - 0.5d0*rjac**2)
        else
           omega = 1.d0/(1.d0 - 0.25d0*rjac**2*omega)
        end if

     end do

!    Reduce residual norm across processors.

     call MPI_Allreduce(anorm,ganorm,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!    Synchronize across processors.

     if (size>1) then
        call sync(sor_u)
     end if

!    Robin boundary conditions.

     if (bound=='robin') then

!       Robin condition on r boundary (second order).

        if (mod(rank+1,nprocr)==0) then
           cr = r(imax-ghost,1)/dr
           do j=2,jmax-1
              cz = z(1,j-ghost)/dz
              sor_u(imax,j) = (uinf + 0.5d0*cr*(4.d0*sor_u(imax-1,j) - sor_u(imax-2,j)) &
                            - 0.5d0*cz*(sor_u(imax,j+1) - sor_u(imax,j-1))) &
                            /(1.d0 + 1.5d0*cr)
           end do
        end if

!       Robin condition on upper z boundary (second order).

        if (rank>=size-nprocr) then
           cz = z(1,jmax-ghost)/dz
           do i=2,imax-1
              cr = r(i-ghost,1)/dr
              sor_u(i,jmax) = (uinf + 0.5d0*cz*(4.d0*sor_u(i,jmax-1) - sor_u(i,jmax-2)) &
                            - 0.5d0*cr*(sor_u(i+1,jmax) - sor_u(i-1,jmax))) &
                            /(1.d0 + 1.5d0*cz)
           end do
        end if

!       Robin condition on lower z boundary (second order).

        if ((.not.eqsym).and.(rank<nprocr)) then
           cz = z(1,1-ghost)/dz
           do i=2,imax-1
              cr = r(i-ghost,1)/dr
              sor_u(i,1) = (uinf - 0.5d0*cz*(4.d0*sor_u(i,2) - sor_u(i,3)) &
                         - 0.5d0*cr*(sor_u(i+1,1) - sor_u(i-1,1))) &
                         /(1.d0 - 1.5d0*cz)
           end do
        end if

!       Upper corner (second order).

        if ((mod(rank+1,nprocr)==0).and.(rank>=size-nprocr)) then
           cr = r(imax-ghost,jmax-ghost)/dr
           cz = z(imax-ghost,jmax-ghost)/dz
           sor_u(imax,jmax) = (uinf &
                            + 0.5d0*cr*(4.d0*sor_u(imax-1,jmax) - sor_u(imax-2,jmax)) &
                            + 0.5d0*cz*(4.d0*sor_u(imax,jmax-1) - sor_u(imax,jmax-2))) &
                            /(1.d0 + 1.5d0*(cr + cz))
        end if

!       Lower corner (second order).

        if ((.not.eqsym).and.(rank<nprocr).and.(mod(rank+1,nprocr)==0)) then
           cr = r(imax-ghost,1-ghost)/dr
           cz = z(imax-ghost,1-ghost)/dz
           sor_u(imax,1) = (uinf &
                         + 0.5d0*cr*(4.d0*sor_u(imax-1,1) - sor_u(imax-2,1)) &
                         - 0.5d0*cz*(4.d0*sor_u(imax  ,2) - sor_u(imax  ,3))) &
                         /(1.d0 + 1.5d0*(cr - cz))
        end if

     end if

!    Symmetries on axis and equator.

     if (ownaxis) then
        do i=1,ghost
           sor_u(i,:) = sor_u(2*ghost+1-i,:)
        end do
     end if

     if ((eqsym).and.(ownequator)) then
        do j=1,ghost
           sor_u(:,j) = sor_u(:,2*ghost+1-j)
        end do
     end if

!    If we reached the desired tolerance return.

     !print  *, rank,n,anorm,ganorm

     if (ganorm<ELL_epsilon) then
        if (rank==0) then
           write (*,'(A,i5,A)') ' SOR:   Solution converged after ',n,' iterations!'
        end if
        return
     end if

  end do

! If we get here we didn't converge.

  if (rank==0) then
     write (*,'(A,i6,A)') ' SOR:   Iterations did not converge after ',ELL_maxiter,' iterations.'
  end if

  !call die


! ***************
! ***   END   ***
! ***************

  end subroutine sor



