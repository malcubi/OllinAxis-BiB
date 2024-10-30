!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/procinfo.f90,v 1.10 2019/11/28 18:57:58 malcubi Exp $

! *********************************
! ***   PROCESSOR INFORMATION   ***
! *********************************

! size              Number of processors for the run.
! rank              Local processor number.
! ierr              MPI error code.
! errorcode         MPI error code.
!
! nprocr            Number of processors in r direction.
! nprocz            Number of processors in z direction.
!
! Nrl:              Array with local number of grid points in r direction.
! Nzl:              Array with local number of grid points in z direction.
!
! Nminl_r:          Array with local first grid point in r.
! Nmaxl_r:          Array with local last grid point in r.
! Nminl_z:          Array with local first grid point in z.
! Nmaxl_z:          Array with local last grid point in z.
!
! eqz:              Array with local position of equator.
! origin:           Array with local position of origin.

  module procinfo

  integer :: size
  integer :: rank
  integer :: ierr
  integer :: nprocr,nprocz

! These arrays have one index that identifies the refinement box.

  integer, allocatable, dimension(:) :: Nrmaxl,Nzmaxl  

! These arrays have two indices: one that identifies on which
! refinement box we are, and one for the processor number.

  integer, allocatable, dimension(:,:) :: Nrl,Nzl         ! Local number of points.
  integer, allocatable, dimension(:,:) :: Nminl_r,Nminl_z ! Local first grid point.
  integer, allocatable, dimension(:,:) :: Nmaxl_r,Nmaxl_z ! Local last grid point.

  integer, allocatable, dimension(:,:) :: eqz             ! Local position of equator.
  integer, allocatable, dimension(:,:) :: axis            ! Local position of axis.
  integer, allocatable, dimension(:,:) :: origin          ! Local position of origin.

  end module procinfo

