!#include "ID_SF_utils.h"

  subroutine idata_procastar_ollinsphere

! ****************************************************
! ***   PROCA STAR INITIAL DATA FROM OLLINSPHERE   ***
! ****************************************************

! This subroutine import Proca Star initial data from Ollinsphere

  use param
  use arrays
  use procinfo
  use derivatives

  implicit none     
  
! Declaration of local variables

  logical :: contains

  integer :: box,level
  integer :: i,j,ir     ! Counters
  integer :: Nrtot,ios              

  real(8), allocatable :: metric_g(:), alpha_g(:)  ! radial metric and lapse global arrays.
  real(8), allocatable :: Aphi_g(:), Ar_g(:)       ! Scalar and vector potentials global arrays.
  real(8), allocatable :: Er_g(:)                  ! Electric field arrays.
  real(8), allocatable :: eAphi_g(:), eAr_g(:)     ! Scalar and vector potentials global arrays.
  real(8), allocatable :: eEr_g(:)                 ! Electric field arrays.

  real(8) :: col1,col2                             ! Auxiliar variables to read data.
  real(8) :: r1,r2,r3,r4                           ! Auxiliar r's for cubic interpolation.
  real(8) :: aux0,aux1,aux2,aux3,aux4              ! Auxiliar variables for cubic interpolation.
  real(8) :: radio,rcil
  real(8) :: half,smallpi

  real(8) :: dradio_dr, dradio_dp, dradio_dz       ! Jacobian components.
  real(8) :: dtheta_dr, dtheta_dp, dtheta_dz       ! Jacobian components.
  real(8) :: dphi_dr, dphi_dp, dphi_dz             ! Jacobian components.

  real(8) :: dr_dradio, dr_dtheta, dr_dphi         ! Jacobian components.
  real(8) :: dp_dradio, dp_dtheta, dp_dphi         ! Jacobian components.
  real(8) :: dz_dradio, dz_dtheta, dz_dphi         ! Jacobian components.

  character(len=128) :: dir   ! For the directory of initial data
  character(len=256) :: line


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0
  smallpi = acos(-1.d0)


! *********************
! ***   READ DATA   ***
! *********************

! Write the directory containing the initial data.

  dir = "/home/jorge/OllinSphere-BiB/exe/procastar/"

  if (contains(mattertype,"electric")) then
     dir = "/home/jorge/OllinSphere-BiB/exe/chargedprocastar/"
  end if

! Calculate Nrtot.

  Nrtot = - ghost

! row = 0.

  open(unit=10, file=trim(dir)//"alpha0.rl", status="old", action="read")

  do

     read(10, '(A)', iostat=ios) line

     if (ios /= 0) exit

     ! Stop if line is empty.

     if (len_trim(line) == 0) exit

     ! Ignore headers.

     if (line(1:1) /= "#") then
       Nrtot = Nrtot + 1
     end if

  end do

  close(10)

! Assignement of memory.

  allocate(metric_g(1-ghost:Nrtot))
  allocate(alpha_g(1-ghost:Nrtot))
  allocate(Aphi_g(1-ghost:Nrtot))
  allocate(Ar_g(1-ghost:Nrtot))
  allocate(Er_g(1-ghost:Nrtot))

! Read the file used to save the data.

  open(unit=10, file=trim(dir)//"cprocaPhi_R0.rl", status="old")     
     read(10, '(A)')               ! Skip header     
     do i = 1-ghost, Nrtot
        read(10, *) col1, col2
        Aphi_g(i) = col2
     end do
  close(10)
       
  open(unit=10, file=trim(dir)//"cprocaA_I0.rl", status="old")     
     read(10, '(A)')               ! Skip header
     do i = 1-ghost, Nrtot
        read(10, *) col1, col2
        Ar_g(i) = col2
     end do
  close(10)

  open(unit=10, file=trim(dir)//"cprocaE_R0.rl", status="old")     
     read(10, '(A)')               ! Skip header
     do i = 1-ghost, Nrtot
        read(10, *) col1, col2
        Er_g(i) = col2
     end do
  close(10)

  open(unit=10, file=trim(dir)//"alpha0.rl", status="old")     
     read(10, '(A)')               ! Skip header
     do i = 1-ghost, Nrtot
        read(10, *) col1, col2
        alpha_g(i) = col2
     end do
  close(10)

  open(unit=10, file=trim(dir)//"A0.rl", status="old")     
     read(10, '(A)')               ! Skip header
     do i = 1-ghost, Nrtot
        read(10, *) col1, col2
        metric_g(i) = col2
     end do
  close(10)


! **********************
! ***   START LOOP   ***
! **********************

  do box=0,Nb
     do level=min(1,box),Nl(box)

        call currentgrid(box,level,grid(box,level))

        do j = lbound(r,2), ubound(r,2)
           do i = lbound(r,1), ubound(r,1)

           radio = rr(i,j)
           rcil = r(i,j)

!          Cylindrycal coordinates: (r,p,z)
!          Spherical coordinates: (radio,theta,phi)

!          Jacobian components.

           dradio_dr = r(i,j)/radio
           dradio_dp = 0.d0
           dradio_dz = z(i,j)/radio

           dtheta_dr = z(i,j)/radio**2
           dtheta_dp = 0.d0
           dtheta_dz = -r(i,j)/radio**2

           dphi_dr = 0.d0
           dphi_dp = 1.d0
           dphi_dz = 0.d0

!          Jacobian components (inverses).

           dr_dradio = r(i,j)/radio
           dr_dtheta = z(i,j)
           dr_dphi   = 0.d0

           dp_dradio = 0.d0
           dp_dtheta = 0.d0
           dp_dphi   = 1.d0

           dz_dradio = z(i,j)/radio 
           dz_dtheta = -r(i,j)
           dz_dphi   = 0.d0        

!          Lagrange cubic interpolation.

           ir = int(radio/dr0 + 0.5d0) ! To make sure that r2<radio<r3

           r1 = (ir - 1.5) * dr0
           r2 = (ir - 0.5) * dr0
           r3 = (ir + 0.5) * dr0
           r4 = (ir + 1.5) * dr0

!          Proca scalar potential.

           aux1 = Aphi_g(ir-1)
           aux2 = Aphi_g(ir  )
           aux3 = Aphi_g(ir+1)
           aux4 = Aphi_g(ir+2)

           aux0 = aux1 * ((radio-r2)*(radio-r3)*(radio-r4))/((r1-r2)*(r1-r3)*(r1-r4)) + &
                  aux2 * ((radio-r1)*(radio-r3)*(radio-r4))/((r2-r1)*(r2-r3)*(r2-r4)) + &
                  aux3 * ((radio-r1)*(radio-r2)*(radio-r4))/((r3-r1)*(r3-r2)*(r3-r4)) + &
                  aux4 * ((radio-r1)*(radio-r2)*(radio-r3))/((r4-r1)*(r4-r2)*(r4-r3))

           proc_phiR(i,j) = aux0
           proc_phiI(i,j) = 0.d0

!          Proca vector potential.

           aux1 = Ar_g(ir-1)
           aux2 = Ar_g(ir  )
           aux3 = Ar_g(ir+1)
           aux4 = Ar_g(ir+2)

           aux0 = aux1 * ((radio-r2)*(radio-r3)*(radio-r4))/((r1-r2)*(r1-r3)*(r1-r4)) + &
                  aux2 * ((radio-r1)*(radio-r3)*(radio-r4))/((r2-r1)*(r2-r3)*(r2-r4)) + &
                  aux3 * ((radio-r1)*(radio-r2)*(radio-r4))/((r3-r1)*(r3-r2)*(r3-r4)) + &
                  aux4 * ((radio-r1)*(radio-r2)*(radio-r3))/((r4-r1)*(r4-r2)*(r4-r3))

           proc_AR_r(i,j) = 0.d0
           proc_AR_z(i,j) = 0.d0

           proc_AI_r(i,j) = dradio_dr * aux0
           proc_AI_z(i,j) = dradio_dz * aux0

!          Proca electric field.

           aux1 = Er_g(ir-1)
           aux2 = Er_g(ir  )
           aux3 = Er_g(ir+1)
           aux4 = Er_g(ir+2)

           aux0 = aux1 * ((radio-r2)*(radio-r3)*(radio-r4))/((r1-r2)*(r1-r3)*(r1-r4)) + &
                  aux2 * ((radio-r1)*(radio-r3)*(radio-r4))/((r2-r1)*(r2-r3)*(r2-r4)) + &
                  aux3 * ((radio-r1)*(radio-r2)*(radio-r4))/((r3-r1)*(r3-r2)*(r3-r4)) + &
                  aux4 * ((radio-r1)*(radio-r2)*(radio-r3))/((r4-r1)*(r4-r2)*(r4-r3))

           proc_ER_r(i,j) = dradio_dr * aux0
           proc_ER_z(i,j) = dradio_dz * aux0

           proc_EI_r(i,j) = 0.d0
           proc_EI_z(i,j) = 0.d0

!          Proca magnetic field.

           proc_BR_r = 0.0d0
           proc_BR_z = 0.0d0

           proc_BI_r = 0.0d0
           proc_BI_z = 0.0d0

!          If we consider initial data for charged Proca star

!          if (contains(mattertype,"electric")) then
            ! ...Pending
!          end if

!          Lapse alpha.

           aux1 = alpha_g(ir-1)
           aux2 = alpha_g(ir  )
           aux3 = alpha_g(ir+1)
           aux4 = alpha_g(ir+2)

           aux0 = aux1 * ((radio-r2)*(radio-r3)*(radio-r4))/((r1-r2)*(r1-r3)*(r1-r4)) + &
                  aux2 * ((radio-r1)*(radio-r3)*(radio-r4))/((r2-r1)*(r2-r3)*(r2-r4)) + &
                  aux3 * ((radio-r1)*(radio-r2)*(radio-r4))/((r3-r1)*(r3-r2)*(r3-r4)) + &
                  aux4 * ((radio-r1)*(radio-r2)*(radio-r3))/((r4-r1)*(r4-r2)*(r4-r3))

           alpha(i,j) = aux0

!          Metric.

           aux1 = metric_g(ir-1)
           aux2 = metric_g(ir  )
           aux3 = metric_g(ir+1)
           aux4 = metric_g(ir+2)

           aux0 = &
             aux1 * ((radio-r2)*(radio-r3)*(radio-r4))/((r1-r2)*(r1-r3)*(r1-r4)) + &
             aux2 * ((radio-r1)*(radio-r3)*(radio-r4))/((r2-r1)*(r2-r3)*(r2-r4)) + &
             aux3 * ((radio-r1)*(radio-r2)*(radio-r4))/((r3-r1)*(r3-r2)*(r3-r4)) + &
             aux4 * ((radio-r1)*(radio-r2)*(radio-r3))/((r4-r1)*(r4-r2)*(r4-r3))

           A(i,j) = dradio_dr**2 * aux0 + dtheta_dr**2 * radio**2

           B(i,j) = dradio_dz**2 * aux0 + dtheta_dz**2 * radio**2

           C(i,j) = (dradio_dr * dradio_dz * aux0 +    &
                     dtheta_dr * dtheta_dz * radio**2) * (1/rcil)

           H(i,j) = dphi_dp**2

!          Extrinsic curvature.

           KA(i,j) = 0.d0
           KB(i,j) = 0.d0
           KH(i,j) = 0.d0
           KC(i,j) = 0.d0

           end do
        end do

     psi = 1.d0
!    psi2 = psi**2
!    psi4 = psi2**2

!    phi = log(psi)
!    chi = one/psi**dble(chipower)


! ********************
! ***   END LOOP   ***
! ********************

     end do
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine idata_procastar_ollinsphere


