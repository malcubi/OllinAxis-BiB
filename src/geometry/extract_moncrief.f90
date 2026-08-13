
  SUBROUTINE extract_moncrief

! *************************************************
! ***   MONCRIEF GRAVITATIONAL WAVE EXTRACTOR   ***
! *************************************************

! Routine for gravitational wave extraction using the Moncrief-Zerilli
! gauge-invariant perturbation formalism on 2D cylindrical grids.
!
! extract_moncrief    Calculates the polar metric gauge-invariant
!                     perturbation amplitude (Moncrief Q) at fixed radii
!                     and outputs multipolar modes to disk.
!
! Originally written by Carlos Joaquin.

  USE mpi
  USE param
  USE arrays
  USE procinfo, ONLY: ierr
  USE harmonix
  USE gw_utils

  IMPLICIT NONE

! Local variables.
  
  LOGICAL, SAVE :: firstcall = .TRUE.
  LOGICAL, SAVE :: capture_background = .TRUE.
  LOGICAL :: found_flag
  
  REAL(8), DIMENSION(1:3) :: radii
  
! Integration results at current time step.

  REAL(8), DIMENSION(1:3, 1:3) :: Moncrief_Q 
  
! Storage for initial static background at t=0.

  REAL(8), DIMENSION(1:3, 1:3), SAVE :: Moncrief_Background
  
  INTEGER :: box, level, bb, ll
  INTEGER :: i, j, k, l_mode
  INTEGER :: Ntheta

  REAL(8) :: dth, theta, smallpi, half  
  REAL(8) :: r0, z0, R_ext
  REAL(8) :: sin_th, cos_th

! Interpolated metric and gauge variables.

  REAL(8) :: loc_A, loc_B, loc_H, loc_C, loc_phi, loc_alp
  REAL(8) :: loc_KTA, loc_KTB, loc_KTH, loc_KTC, loc_trK, loc_betar, loc_betaz
  REAL(8), DIMENSION(13) :: send_buf, recv_buf
  REAL(8) :: val_A, val_B, val_H, val_C, val_phi, val_alp
  REAL(8) :: val_KTA, val_KTB, val_KTH, val_KTC, val_trK, val_betar, val_betaz
  
! Reconstructed physical quantities and source term.

  REAL(8) :: psi4_val, g_rho_rho, g_zz, g_rho_z
  REAL(8) :: g_th_th, g_phiph
  REAL(8) :: K_rr, K_zz, K_rz, K_thth, K_phiph
  REAL(8) :: Source_Term
  REAL(8) :: Legendre_P, simpson_w

  REAL(8) :: factor_l, fac_num, fac_den
  
! Output file handlers.

  CHARACTER(4)   :: filen
  CHARACTER(20)  :: filestatus
  CHARACTER(100) :: filename_out


! **********************************************
! ***   INITIAL SETUP AND PARAMETER CHECKS   ***
! **********************************************

! Exit if wave extraction is not requested.

  IF (.NOT. wave_extract) RETURN

  half = 0.5d0
  smallpi = ACOS(-1.d0)

! Set angular resolution based on equatorial symmetry.

  IF (eqsym) THEN
     Ntheta = 250
     dth = half*(smallpi)/dble(Ntheta)
  ELSE
     Ntheta = 500
     dth = (smallpi)/dble(Ntheta)
  END IF

! Store extraction radii.

  radii(1) = rad1
  radii(2) = rad2
  radii(3) = rad3
  
  Moncrief_Q = 0.0d0


! ***********************************************
! ***   MAIN ANGULAR AND RADIUS INTEGRATION   ***
! ***********************************************

  DO j=1,3

     R_ext = radii(j)

!    Validate extraction radius.

     IF (R_ext <= 0.d0) CYCLE
!    IF (R_ext > min(Nrtotal,Nztotal)*min(dr0,dz0)) CYCLE 

     DO i = 0, Ntheta

        theta = dble(i) * dth
        sin_th = SIN(theta)
        cos_th = COS(theta)

        r0 = R_ext * sin_th
        z0 = R_ext * cos_th

!       Find appropriate grid box and level for point (r0, z0).

        box = 0
        level = 0

        IF (Nlmax >= 0) THEN
           DO bb = 0, Nb
              DO ll = 0, Nl(bb)
                 IF ((r0 > rminl(bb,ll)+drl(ll)) .AND. &
                     (r0 < rmaxl(bb,ll)-drl(ll)) .AND. &
                     (z0 > zminl(bb,ll)+dzl(ll)) .AND. &
                     (z0 < zmaxl(bb,ll)-dzl(ll))) THEN
                    box = bb
                    level = ll
                 END IF
              END DO
           END DO
        END IF

!       Interpolate BSSN fields at point (r0, z0).

        CALL Interpolate_BSSN_Point(box, level, r0, z0, &
             loc_A, loc_B, loc_H, loc_C, loc_phi, loc_alp, &
             loc_KTA, loc_KTB, loc_KTH, loc_KTC, loc_trK, &
             loc_betar, loc_betaz, found_flag)

         send_buf = 0.d0

         IF (found_flag) THEN

           send_buf(1) = loc_A
           send_buf(2) = loc_B
           send_buf(3) = loc_H
           send_buf(4) = loc_C
           send_buf(5) = loc_phi
           send_buf(6) = loc_alp
           send_buf(7) = loc_KTA
           send_buf(8) = loc_KTB
           send_buf(9) = loc_KTH
           send_buf(10) = loc_KTC
           send_buf(11) = loc_trK
           send_buf(12) = loc_betar
           send_buf(13) = loc_betaz

        END IF

!       Reduce interpolated values across all processors.

        CALL MPI_ALLREDUCE(send_buf, recv_buf, 13, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        val_A   = recv_buf(1)
        val_B   = recv_buf(2)
        val_H   = recv_buf(3)
        val_C   = recv_buf(4)
        val_phi = recv_buf(5)
        val_alp = recv_buf(6)
        val_KTA = recv_buf(7)
        val_KTB = recv_buf(8)
        val_KTH = recv_buf(9)
        val_KTC = recv_buf(10)
        val_trK = recv_buf(11)

!       Reconstruct physical 3-metric components and regularize polar coordinate boundaries.

        psi4_val = EXP(4.0d0 * val_phi)

        g_rho_rho = psi4_val * val_A
        g_zz      = psi4_val * val_B
        g_rho_z   = psi4_val * r0 * val_C

!       Compute normalized phi-phi metric term using BSSN field H to avoid polar division by zero.

        g_phiph = psi4_val * (R_ext*sin_th)**2 * val_H

!       Tensor transformation for g_theta_theta component.

        g_th_th = (R_ext * cos_th)**2 * g_rho_rho &
                + (-R_ext * sin_th)**2 * g_zz &
                + 2.0d0 * (R_ext * cos_th) * (-R_ext * sin_th) * g_rho_z

!       Compute metric perturbation source term.

         K_thth  = R_ext**2*((psi4_val*(val_KTA + (1.d0/3.d0)*val_A*val_trK))*cos_th**2 &
                           + (psi4_val*(val_KTB + (1.d0/3.d0)*val_B*val_trK))*sin_th**2) 

        Source_Term = (psi4_val * val_H) - (g_th_th / R_ext**2)

!       Set Simpson's integration weights.

        IF (i==0 .OR. i==Ntheta) THEN
           simpson_w = 1.0d0
        ELSE IF (MOD(i,2) == 0) THEN
           simpson_w = 2.0d0
        ELSE
           simpson_w = 4.0d0
        END IF

!       Project source term onto spin-weighted spherical harmonics (s=-2).

        DO k=1,3

           l_mode = 2*k

!          Legendre_P = SWSHR(-2, l_mode, 0, theta, 0.0d0)
           Legendre_P = SWSHR(0, l_mode, 0, theta, 0.0d0)

           
           Moncrief_Q(j, k) = Moncrief_Q(j, k) + &
                simpson_w * Source_Term * Legendre_P * sin_th

        END DO

     END DO ! End Theta
  END DO ! End Radii

! Final angular integration normalization.

  DO k=1,3

      l_mode = 2*k

      fac_num  = FACT(l_mode - 2)
      fac_den  = FACT(l_mode + 2)
      factor_l = SQRT(fac_num / fac_den)

      IF (eqsym) THEN
        Moncrief_Q(:, k) = Moncrief_Q(:, k) * (4.0d0 * smallpi) * (dth / 3.0d0) * factor_l
      ELSE
        Moncrief_Q(:, k) = Moncrief_Q(:, k) * (2.0d0 * smallpi) * (dth / 3.0d0) * factor_l
      END IF

   END DO


! ***************************************************
! ***   BACKGROUND SUBTRACTION (OFFSET REMOVAL)   ***
! ***************************************************

  IF (capture_background) THEN

!    Store current value at t=0 as initial background.

     Moncrief_Background = Moncrief_Q
     capture_background = .FALSE.
     
!    Set initial perturbation signal to zero.

     Moncrief_Q = 0.0d0

     RETURN

  ELSE

!    Subtract initial static background to isolate physical radiation.

     Moncrief_Q = Moncrief_Q - Moncrief_Background

  END IF


! *********************************
! ***   SAVE MULTIPOLAR MODES   ***
! *********************************

  IF (rank == 0) THEN

!    Determine file status and extension.

     IF (firstcall) THEN
        firstcall = .FALSE.
        filestatus = 'replace'
     ELSE
        filestatus = 'old'
     END IF

     WRITE(filen, '(A)') '.tl'

!    Open files.

     IF (filestatus == 'replace') THEN

        OPEN(1, file=trim(directory)//'/moncrief_20'//trim(filen), form='formatted', &
             status=filestatus)
        OPEN(2, file=trim(directory)//'/moncrief_40'//trim(filen), form='formatted', &
             status=filestatus)
        OPEN(3, file=trim(directory)//'/moncrief_60'//trim(filen), form='formatted', &
             status=filestatus)

!       Write the extraction radii.

        WRITE(1, "(A,3ES20.8)") "extraction radii=", rad1, rad2, rad3
        WRITE(2, "(A,3ES20.8)") "extraction radii=", rad1, rad2, rad3
        WRITE(3, "(A,3ES20.8)") "extraction radii=", rad1, rad2, rad3

     ELSE

        OPEN(1, file=trim(directory)//'/moncrief_20'//trim(filen), form='formatted', &
             status=filestatus, position='append')
        OPEN(2, file=trim(directory)//'/moncrief_40'//trim(filen), form='formatted', &
             status=filestatus, position='append')
        OPEN(3, file=trim(directory)//'/moncrief_60'//trim(filen), form='formatted', &
             status=filestatus, position='append')

     END IF

!    Save multipolar mode. The columns order are time, rad1, rad2, rad3.

     WRITE(1, "(4ES16.8)") time, Moncrief_Q(1,1), Moncrief_Q(2,1), Moncrief_Q(3,1)
     WRITE(2, "(4ES16.8)") time, Moncrief_Q(1,2), Moncrief_Q(2,2), Moncrief_Q(3,2)
     WRITE(3, "(4ES16.8)") time, Moncrief_Q(1,3), Moncrief_Q(2,3), Moncrief_Q(3,3)

!    Close files.

     CLOSE(1)
     CLOSE(2)
     CLOSE(3)

  END IF


! ***************
! ***   END   ***
! ***************

  END SUBROUTINE extract_moncrief
