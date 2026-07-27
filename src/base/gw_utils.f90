! *******************************************************
! ***   GRAVITATIONAL WAVE EXTRACTION UTILITIES     ***
! *******************************************************

! Here I define utility modules and routines for gravitational
! wave extraction:
!
! Interpolate_BSSN_Point    Interpolates all BSSN metric, gauge, shift,
!                           and extrinsic curvature variables at a given 2D point.


MODULE gw_utils

! Include modules.

  USE param
  USE arrays
  USE procinfo
  USE mpi 

  IMPLICIT NONE

CONTAINS

! ***************************************************
! ***   SUBROUTINE INTERPOLATE_BSSN_POINT         ***
! ***************************************************

! This subroutine interpolates all BSSN spatial metric components,
! conformal factor, lapse, shift, and extrinsic curvature variables
! at an arbitrary extraction point (r0, z0) on a given refinement level.

  SUBROUTINE Interpolate_BSSN_Point(box, level, r0, z0, &
       val_A, val_B, val_H, val_C, val_phi, val_alpha, &
       val_KTA, val_KTB, val_KTH, val_KTC, val_trK, &
       val_betar, val_betaz, found_flag)

IMPLICIT NONE
    
! Input variables: Extraction point coordinates and grid info.

    INTEGER, INTENT(IN) :: box, level
    REAL(8), INTENT(IN) :: r0, z0

! Output variables: Metric, gauge, extrinsic curvature, and shift.

    REAL(8), INTENT(OUT) :: val_A, val_B, val_H, val_C, val_phi, val_alpha
    REAL(8), INTENT(OUT) :: val_KTA, val_KTB, val_KTH, val_KTC, val_trK
    REAL(8), INTENT(OUT) :: val_betar, val_betaz

    LOGICAL, INTENT(OUT) :: found_flag

! Local variables.

    INTEGER :: i0, j0, i, j, m, n
    
    REAL(8) :: deltar, deltaz, auxr, auxz, aux
    REAL(8) :: rmin, rmax, zmin, zmax, small
    LOGICAL :: forcelinear

! Initialize output variables.

    val_A = 0.d0; val_B = 0.d0; val_H = 0.d0
    val_C = 0.d0; val_phi = 0.d0; val_alpha = 0.d0
    val_KTA = 0.d0; val_KTB = 0.d0; val_KTH = 0.d0
    val_KTC = 0.d0; val_trK = 0.d0
    val_betar = 0.d0; val_betaz = 0.d0

    found_flag = .FALSE.
    small = 1.d-8 * min(deltar, deltaz)


! ***************************************************************
! ***   CHECK IF CURRENT PROCESSOR OWNS THE EXTRACTION POINT   ***
! ***************************************************************

    IF (size > 1) THEN

!      Find local processor boundaries taking ghost cells into account.

       IF (mod(rank,nprocr)==0) THEN
          rmin = grid(box,level)%r(1-ghost,0) - small
       ELSE
          rmin = grid(box,level)%r(1,0)
       END IF
       
       IF (mod(rank+1,nprocr)==0) THEN
          rmax = grid(box,level)%r(Nrl(box,rank),0)
       ELSE
          rmax = grid(box,level)%r(Nrl(box,rank)-ghost+1,0)
       END IF

       IF (rank < nprocr) THEN
          zmin = grid(box,level)%z(0,1-ghost) - small
       ELSE
          zmin = grid(box,level)%z(0,1)
       END IF
       
       IF (rank > size-nprocr-1) THEN
          zmax = grid(box,level)%z(0,Nzl(box,rank))
       ELSE
          zmax = grid(box,level)%z(0,Nzl(box,rank)-ghost+1)
       END IF

!      If point lies outside local domain, exit subroutine.

       IF ((r0 <= rmin).OR.(r0 > rmax).OR.(z0 <= zmin).OR.(z0 > zmax)) RETURN
    END IF
    
! Flag point as found on current MPI processor.

    found_flag = .TRUE.


! *********************************************
! ***   FIND GRID INDICES AND DISTANCES     ***
! *********************************************

    i0 = int((r0 - grid(box,level)%r(1-ghost,0))/drl(level)) + 1 - ghost
    j0 = int((z0 - grid(box,level)%z(0,1-ghost))/dzl(level)) + 1 - ghost
    
    deltar = r0 - grid(box,level)%r(i0,0)
    deltaz = z0 - grid(box,level)%z(0,j0)

! Determine interpolation order (force bilinear near boundaries).

    forcelinear = .FALSE.
    IF (interporder == "bilinear") forcelinear = .TRUE.
    IF (((rmaxl(box,level)-r0)<drl(level)).OR.((zmaxl(box,level)-z0)<dzl(level)).OR. &
        ((r0-rminl(box,level))<drl(level)).OR.((z0-zminl(box,level))<dzl(level))) THEN
       forcelinear = .TRUE.
    END IF

    auxr = deltar/drl(level)
    auxz = deltaz/dzl(level)


! *******************************************
! ***   INTERPOLATE ALL BSSN QUANTITIES   ***
! *******************************************

! Perform bilinear interpolation.

    IF (forcelinear) THEN
       DO i=0,1
          DO j=0,1
             aux = (1.d0-auxr)**(1-i) * auxr**i * (1.d0-auxz)**(1-j) * auxz**j
            
!            Metric variables.
             val_A     = val_A     + aux * grid(box,level)%A(i0+i, j0+j)
             val_B     = val_B     + aux * grid(box,level)%B(i0+i, j0+j)
             val_H     = val_H     + aux * grid(box,level)%H(i0+i, j0+j)
             val_C     = val_C     + aux * grid(box,level)%C(i0+i, j0+j)
             val_phi   = val_phi   + aux * grid(box,level)%phi(i0+i, j0+j)
             val_alpha = val_alpha + aux * grid(box,level)%alpha(i0+i, j0+j)
             
!            Extrinsic curvature variables.
             val_KTA   = val_KTA   + aux * grid(box,level)%KTA(i0+i, j0+j)
             val_KTB   = val_KTB   + aux * grid(box,level)%KTB(i0+i, j0+j)
             val_KTH   = val_KTH   + aux * grid(box,level)%KTH(i0+i, j0+j)
             val_KTC   = val_KTC   + aux * grid(box,level)%KTC(i0+i, j0+j)
             val_trK   = val_trK   + aux * grid(box,level)%trK(i0+i, j0+j)

!            Shift vector components (if active).
             IF (shift /= "none") THEN
                val_betar = val_betar + aux * grid(box,level)%beta_r(i0+i, j0+j)
                val_betaz = val_betaz + aux * grid(box,level)%beta_z(i0+i, j0+j)
             END IF

          END DO
       END DO
    ELSE

! Perform bicubic interpolation.

       DO i=-1,2
          DO j=-1,2
             aux = 1.d0
             DO m=-1,2
                IF (m==i) CYCLE
                aux = aux*(auxr-dble(m))/dble(i-m)
             END DO
             DO n=-1,2
                IF (n==j) CYCLE
                aux = aux*(auxz-dble(n))/dble(j-n)
             END DO
             
!            Metric variables.
             val_A     = val_A     + aux * grid(box,level)%A(i0+i, j0+j)
             val_B     = val_B     + aux * grid(box,level)%B(i0+i, j0+j)
             val_H     = val_H     + aux * grid(box,level)%H(i0+i, j0+j)
             val_C     = val_C     + aux * grid(box,level)%C(i0+i, j0+j)
             val_phi   = val_phi   + aux * grid(box,level)%phi(i0+i, j0+j)
             val_alpha = val_alpha + aux * grid(box,level)%alpha(i0+i, j0+j)

!            Extrinsic curvature variables.
             val_KTA   = val_KTA   + aux * grid(box,level)%KTA(i0+i, j0+j)
             val_KTB   = val_KTB   + aux * grid(box,level)%KTB(i0+i, j0+j)
             val_KTH   = val_KTH   + aux * grid(box,level)%KTH(i0+i, j0+j)
             val_KTC   = val_KTC   + aux * grid(box,level)%KTC(i0+i, j0+j)
             val_trK   = val_trK   + aux * grid(box,level)%trK(i0+i, j0+j)

!            Shift vector components (if active).
             IF (shift /= "none") THEN
                val_betar = val_betar + aux * grid(box,level)%beta_r(i0+i, j0+j)
                val_betaz = val_betaz + aux * grid(box,level)%beta_z(i0+i, j0+j)
             END IF

          END DO
       END DO
    END IF

! End.


END SUBROUTINE Interpolate_BSSN_Point

END MODULE gw_utils