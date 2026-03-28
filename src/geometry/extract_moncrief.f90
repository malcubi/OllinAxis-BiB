SUBROUTINE extract_moncrief

  ! ******************************************************
  ! *** MONCRIEF-ZERILLI GRAVITATIONAL WAVE EXTRACTOR  ***
  ! ******************************************************
  
  USE mpi
  USE param
  USE arrays
  USE procinfo
  USE harmonix
  USE GW_EXTRACT_UTILS

  IMPLICIT NONE

  ! --- Variables Locales ---
  LOGICAL, SAVE :: firstcall = .TRUE.
  LOGICAL, SAVE :: capture_background = .TRUE.
  LOGICAL :: found_flag
  
  REAL(8), DIMENSION(1:3) :: radii
  
  ! Resultados de integracion actuales
  REAL(8), DIMENSION(1:3, 1:3) :: Moncrief_Q 
  
  ! Almacenamiento del Fondo (t=0) para restar
  REAL(8), DIMENSION(1:3, 1:3), SAVE :: Moncrief_Background
  
  INTEGER :: box, level, bb, ll
  INTEGER :: i, j, k, l_mode
  INTEGER :: Ntheta
  REAL(8) :: dth, theta, smallpi, half
  
  REAL(8) :: r0, z0, R_ext
  REAL(8) :: sin_th, cos_th
  
  ! Variables interpoladas
  REAL(8) :: loc_A, loc_B, loc_H, loc_C, loc_phi, loc_alp
  REAL(8), DIMENSION(6) :: send_buf, recv_buf
  REAL(8) :: val_A, val_B, val_H, val_C, val_phi, val_alp
  
  ! Variables Físicas
  REAL(8) :: psi4_val, g_rho_rho, g_zz, g_rho_z
  REAL(8) :: g_th_th, term_phi_phi_norm
  REAL(8) :: Source_Term
  REAL(8) :: Legendre_P, simpson_w
  
  ! Archivos
  CHARACTER(4)  :: filen
  CHARACTER(20) :: filestatus
  CHARACTER(100) :: filename_out

  ! **************************
  ! *** CONFIGURACIÓN    ***
  ! **************************

  ! Chequeo de seguridad: Si no se pidio extraer, salir.
  IF (.NOT. wave_extract) RETURN
  
  half = 0.5d0
  smallpi = ACOS(-1.d0)

  IF (eqsym) THEN
     Ntheta = 250
     dth = half*(smallpi)/dble(Ntheta)
  ELSE
     Ntheta = 500
     dth = (smallpi)/dble(Ntheta)
  END IF

  radii(1) = rad1
  radii(2) = rad2
  radii(3) = rad3
  
  Moncrief_Q = 0.0d0

  ! **************************
  ! *** BUCLE PRINCIPAL  ***
  ! **************************

  DO j = 1, 3
     R_ext = radii(j)
     ! Validar radio
     IF (R_ext <= 0.d0) CYCLE
     IF (R_ext > min(Nrtotal,Nztotal)*min(dr0,dz0)) CYCLE 

     DO i = 0, Ntheta
        
        theta = dble(i) * dth
        sin_th = SIN(theta)
        cos_th = COS(theta)
        
        r0 = R_ext * sin_th
        z0 = R_ext * cos_th

        ! 1. ENCONTRAR BOX Y LEVEL
        box = 0
        level = 0
        IF (Nlmax > 0) THEN
           DO bb = 0, Nb
              DO ll = 1, Nl(bb)
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

        ! 2. INTERPOLACIÓN (Usando tu modulo GW_EXTRACT_UTILS)
        CALL Interpolate_BSSN_Point(box, level, r0, z0, &
             loc_A, loc_B, loc_H, loc_C, loc_phi, loc_alp, found_flag)

        send_buf(1) = loc_A
        send_buf(2) = loc_B
        send_buf(3) = loc_H
        send_buf(4) = loc_C
        send_buf(5) = loc_phi
        send_buf(6) = loc_alp

        CALL MPI_ALLREDUCE(send_buf, recv_buf, 6, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        val_A   = recv_buf(1)
        val_B   = recv_buf(2)
        val_H   = recv_buf(3)
        val_C   = recv_buf(4)
        val_phi = recv_buf(5)
        val_alp = recv_buf(6)

        ! 3. RECONSTRUCCIÓN FÍSICA Y ESTABILIDAD EN POLOS
        ! Usamos psi4_val local para evitar conflicto de nombres
        psi4_val = EXP(4.0d0 * val_phi)
        
        g_rho_rho = psi4_val * val_A
        g_zz      = psi4_val * val_B
        g_rho_z   = psi4_val * r0 * val_C
        
        ! --- CORRECCIÓN CLAVE ---
        ! En lugar de g_phi_phi / sin^2, usamos la definicion BSSN H.
        ! g_phi_phi = psi^4 * rho^2 * H = psi^4 * R^2 * sin^2 * H
        ! Por lo tanto: g_phi_phi / sin^2 = psi^4 * R^2 * H
        ! Esto elimina la division por cero en los polos.
        
        term_phi_phi_norm = psi4_val * (R_ext**2) * val_H

        ! Transformacion Tensorial para g_theta_theta
        g_th_th = (R_ext * cos_th)**2 * g_rho_rho &
                + (-R_ext * sin_th)**2 * g_zz &
                + 2.0d0 * (R_ext * cos_th) * (-R_ext * sin_th) * g_rho_z

        ! 4. TERMINO FUENTE DE MONCRIEF
        ! Calculamos la deformación métrica.
        ! La resta analitica de Minkowski (R^2) ya no es estrictamente necesaria
        ! aqui porque restaremos el Background completo al final, pero ayuda
        ! a mantener los numeros pequeños durante la integral.
        
        Source_Term = (term_phi_phi_norm - g_th_th) / (R_ext**2)

        ! 5. INTEGRACIÓN (Spin-Weighted -2)
        ! Simpson's weights
        IF (i==0 .OR. i==Ntheta) THEN
           simpson_w = 1.0d0
        ELSE IF (MOD(i,2) == 0) THEN
           simpson_w = 2.0d0
        ELSE
           simpson_w = 4.0d0
        END IF
        
        DO k = 1, 3
           l_mode = 2*k
           ! Usamos armónicos con peso de espín -2 para filtrar monopolo
           Legendre_P = SWSHR(-2, l_mode, 0, theta, 0.0d0)
           
           Moncrief_Q(j, k) = Moncrief_Q(j, k) + &
                simpson_w * Source_Term * Legendre_P * sin_th
        END DO

     END DO ! Fin Theta
  END DO ! Fin Radios

  ! Normalización Final
  IF (eqsym) THEN
     Moncrief_Q = 4.0d0 * smallpi * Moncrief_Q * dth / 3.0d0
  ELSE
     Moncrief_Q = 2.0d0 * smallpi * Moncrief_Q * dth / 3.0d0
  END IF

  ! ---------------------------------------------------------
  ! SUSTRACCIÓN DE FONDO (BACKGROUND SUBTRACTION)
  ! ---------------------------------------------------------
  ! Esto elimina el offset inicial de la Onda de Brill y el ruido de gauge estático.
!! Corresponde a la versión donde no guarda datos a t = 0.  
!!  IF (capture_background) THEN
!!     Moncrief_Background = Moncrief_Q
!!     capture_background = .FALSE.
     
     ! En el primer paso, la "onda extraída" es cero por definición
!!     Moncrief_Q = 0.0d0 
!!  ELSE
     ! Pasos subsiguientes: Señal = Actual - Inicial
!!     Moncrief_Q = Moncrief_Q - Moncrief_Background
!!  END IF


  ! ***************************************
  ! *** SUSTRACCIÓN DE FONDO (OFFSET) ***
  ! ***************************************
  
  IF (capture_background) THEN
     ! Guardamos el valor actual (t=0) como fondo
     Moncrief_Background = Moncrief_Q
     capture_background = .FALSE.
     
     ! Como acabamos de capturar el fondo, la señal física en t=0 es cero.
     Moncrief_Q = 0.0d0
     
     ! Si estamos solo inicializando (t=0), salimos sin escribir a disco
     RETURN 
  ELSE
     ! Restamos el fondo guardado para obtener solo la perturbación
     Moncrief_Q = Moncrief_Q - Moncrief_Background
  END IF


  ! **************************
  ! *** GUARDAR DATOS    ***
  ! **************************

  IF (rank == 0) THEN
     IF (firstcall) THEN
        firstcall = .FALSE.
        filestatus = 'replace'
     ELSE
        filestatus = 'old'
     END IF
     
     DO k = 1, 3
        l_mode = 2*k
        WRITE(filename_out, "(A,I1,A)") trim(directory)//'/moncrief_l', l_mode, '.mq'

        IF (filestatus == 'replace') THEN
           OPEN(10+k, file=trim(filename_out), form='formatted', status='unknown')
           WRITE(10+k, "(A,3ES14.6)") "# Time vs Strain at R=", rad1, rad2, rad3
        ELSE
           OPEN(10+k, file=trim(filename_out), form='formatted', status='old', position='append')
        END IF

        WRITE(10+k, "(4ES16.8)") time, Moncrief_Q(1,k), Moncrief_Q(2,k), Moncrief_Q(3,k)
        CLOSE(10+k)
     END DO
  END IF

END SUBROUTINE extract_moncrief