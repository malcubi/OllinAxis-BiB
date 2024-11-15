#!/usr/bin/env perl

#$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/prl/arrays.pl,v 1.51 2022/01/11 02:14:45 malcubi Exp $

# This perl script creates the files:
#
# arrays.f90
# accumulate.f90
# allocatearrays.f90
# currentgrid.f90
# grabarray.f90
# mytypes.f90
# saveold.f90
# simpleboundary.f90
# symmetries_r.f90
# symmetries_z.f90
# syncall.f90
# update.f90
#
# Plus the include files:
#
# bound_interp.inc
# restrict_interp.inc
# intersect_interp.inc

print "CREATING FILES FOR DEALING WITH ARRAYS ...\n\n";

# Open outfiles.

open(FILE_ARRAYS,">src/auto/arrays.f90") or die "Can't open grabarrays.f90: $!";
open(FILE_ACCUMULATE,">src/auto/accumulate.f90") or die "Can't open accumulate.f90: $!";
open(FILE_ALLOCATEARRAYS,">src/auto/allocatearrays.f90") or die "Can't open allocatearrays.f90: $!";
open(FILE_CURRENTGRID,">src/auto/currentgrid.f90") or die "Can't open currentgrid.f90: $!";
open(FILE_GRABARRAY,">src/auto/grabarray.f90") or die "Can't open grabarrays.f90: $!";
open(FILE_MYTYPES,">src/auto/mytypes.f90") or die "Can't open mytypes.f90: $!";
open(FILE_SAVEOLD,">src/auto/saveold.f90") or die "Can't open saveold.f90: $!";
open(FILE_SIMPLEBOUNDARY,">src/auto/simpleboundary.f90") or die "Can't open simpleboundary.f90: $!";
open(FILE_SYMMETRIES_R,">src/auto/symmetries_r.f90") or die "Can't open symmetries_r.f90: $!";
open(FILE_SYMMETRIES_Z,">src/auto/symmetries_z.f90") or die "Can't open symmetries_z.f90: $!";
open(FILE_SYNCALL,">src/auto/syncall.f90") or die "Can't open syncall.f90: $!";
open(FILE_UPDATE,">src/auto/update.f90") or die "Can't open update.f90: $!";

open(FILE_BOUNDINTERP,">src/auto/bound_interp.inc") or die "Can't open boundinterp.inc: $!";
open(FILE_RESTRICTINTERP,">src/auto/restrict_interp.inc") or die "Can't open restrict_copy.inc: $!";
open(FILE_INTERSECTINTERP,">src/auto/intersect_interp.inc") or die "Can't open restrict_copy.inc: $!";

# Write beginning of file arrays.f90

print FILE_ARRAYS "! Automatically generated file.  Do not edit!\n\n";
print FILE_ARRAYS "  module arrays\n\n";
print FILE_ARRAYS "  use mytypes\n\n";
print FILE_ARRAYS "  implicit none\n\n";
print FILE_ARRAYS "  integer, allocatable, dimension(:,:) :: s\n";
print FILE_ARRAYS "  real(8), allocatable, dimension(:,:) :: t\n\n";
print FILE_ARRAYS "  real(8), allocatable, dimension(:,:) :: t1\n";
print FILE_ARRAYS "  real(8), allocatable, dimension(:,:) :: t2\n\n";
print FILE_ARRAYS "  real(8), allocatable, dimension(:) :: drl\n";
print FILE_ARRAYS "  real(8), allocatable, dimension(:) :: dzl\n";
print FILE_ARRAYS "  real(8), allocatable, dimension(:) :: dtl\n\n";
print FILE_ARRAYS "  real(8), allocatable, dimension(:,:) :: rminl,rmaxl\n";
print FILE_ARRAYS "  real(8), allocatable, dimension(:,:) :: zminl,zmaxl\n\n";
print FILE_ARRAYS "  character(50), allocatable, dimension(:) :: outvars0Darray\n";
print FILE_ARRAYS "  character(50), allocatable, dimension(:) :: outvars1Darray\n";
print FILE_ARRAYS "  character(50), allocatable, dimension(:) :: outvars2Darray\n\n";
print FILE_ARRAYS "  type(gridfuncs), allocatable :: grid(:,:)\n\n";

# Write beginning of file accumulate.f90

print FILE_ACCUMULATE "! Automatically generated file.  Do not edit!\n\n";
print FILE_ACCUMULATE "  subroutine accumulate(k,niter,w)\n\n";
print FILE_ACCUMULATE "  use param\n";
print FILE_ACCUMULATE "  use arrays\n\n";
print FILE_ACCUMULATE "  implicit none\n\n";
print FILE_ACCUMULATE "  integer k,niter\n\n";
print FILE_ACCUMULATE "  real(8) w\n\n";

# Write beginning of file allocatearrays.f90

print FILE_ALLOCATEARRAYS "! Automatically generated file.  Do not edit!\n\n";
print FILE_ALLOCATEARRAYS "  subroutine allocatearrays(status)\n\n";
print FILE_ALLOCATEARRAYS "  use param\n";
print FILE_ALLOCATEARRAYS "  use arrays\n";
print FILE_ALLOCATEARRAYS "  use procinfo\n\n";
print FILE_ALLOCATEARRAYS "  implicit none\n";
print FILE_ALLOCATEARRAYS "  logical contains\n";
print FILE_ALLOCATEARRAYS "  integer box,level\n";
print FILE_ALLOCATEARRAYS "  character(len=*) status\n\n";
print FILE_ALLOCATEARRAYS "  if (trim(status)=='on') then\n";
print FILE_ALLOCATEARRAYS "     allocate(grid(0:Nb,0:Nlmax))\n";
print FILE_ALLOCATEARRAYS "  else\n";
print FILE_ALLOCATEARRAYS "     deallocate(grid)\n";
print FILE_ALLOCATEARRAYS "  end if\n\n";

# Write beginning of file currentgrid.f90

print FILE_CURRENTGRID "! Automatically generated file.  Do not edit!\n\n";
print FILE_CURRENTGRID "  subroutine currentgrid(box,level,localgrid)\n\n";
print FILE_CURRENTGRID "  use param\n";
print FILE_CURRENTGRID "  use arrays\n";
print FILE_CURRENTGRID "  use procinfo\n";
print FILE_CURRENTGRID "  use mytypes\n\n";
print FILE_CURRENTGRID "  implicit none\n\n";
print FILE_CURRENTGRID "  integer box,level\n\n";
print FILE_CURRENTGRID "  type(gridfuncs) :: localgrid\n\n";
print FILE_CURRENTGRID "  ownaxis = (axis(box,rank)/=-1)\n";
print FILE_CURRENTGRID "  ownequator = (eqz(box,rank)/=-1)\n";
print FILE_CURRENTGRID "  ownorigin = (origin(box,rank)/=-1)\n\n";
print FILE_CURRENTGRID "  Nr = Nrl(box,rank)\n";
print FILE_CURRENTGRID "  Nz = Nzl(box,rank)\n";
print FILE_CURRENTGRID "  Nrmax = Nrmaxl(box)\n";
print FILE_CURRENTGRID "  Nzmax = Nzmaxl(box)\n\n";
print FILE_CURRENTGRID "  time = t(box,level)\n\n";
print FILE_CURRENTGRID "  dt = dtl(level)\n";
print FILE_CURRENTGRID "  dr = drl(level)\n";
print FILE_CURRENTGRID "  dz = dzl(level)\n\n";

# Write beginning of file grabarray.f90

print FILE_GRABARRAY "! Automatically generated file.  Do not edit!\n\n";
print FILE_GRABARRAY "  subroutine grabarray(varname)\n\n";
print FILE_GRABARRAY "  use param\n";
print FILE_GRABARRAY "  use arrays\n";
print FILE_GRABARRAY "  use procinfo\n\n";
print FILE_GRABARRAY "  implicit none\n";
print FILE_GRABARRAY "  logical exists\n";
print FILE_GRABARRAY "  character(len=*) varname\n\n";
print FILE_GRABARRAY "  exists = .false.\n\n";

# Write beginning of file mytypes.f90

print FILE_MYTYPES "! Automatically generated file.  Do not edit!\n\n";
print FILE_MYTYPES "  module mytypes\n\n";
print FILE_MYTYPES "  type gridfuncs\n\n";

# Write beginning of file saveold.f90

print FILE_SAVEOLD "! Automatically generated file.  Do not edit!\n\n";
print FILE_SAVEOLD "  subroutine saveold\n\n";
print FILE_SAVEOLD "  use param\n";
print FILE_SAVEOLD "  use arrays\n\n";
print FILE_SAVEOLD "  implicit none\n\n";
print FILE_SAVEOLD "  integer i\n\n";

# Write beginning of file simpleboundary.f90

print FILE_SIMPLEBOUNDARY "! Automatically generated file.  Do not edit!\n\n";
print FILE_SIMPLEBOUNDARY "  subroutine simpleboundary\n\n";
print FILE_SIMPLEBOUNDARY "  use param\n";
print FILE_SIMPLEBOUNDARY "  use arrays\n\n";
print FILE_SIMPLEBOUNDARY "  implicit none\n\n";

# Write beginning of file symmetries_r.f90

print FILE_SYMMETRIES_R "! Automatically generated file.  Do not edit!\n\n";
print FILE_SYMMETRIES_R "  subroutine symmetries_r\n\n";
print FILE_SYMMETRIES_R "  use param\n";
print FILE_SYMMETRIES_R "  use arrays\n\n";
print FILE_SYMMETRIES_R "  implicit none\n\n";
print FILE_SYMMETRIES_R "  integer i\n\n";
print FILE_SYMMETRIES_R "  do i=1,ghost\n\n";

# Write beginning of file symmetries_z.f90

print FILE_SYMMETRIES_Z "! Automatically generated file.  Do not edit!\n\n";
print FILE_SYMMETRIES_Z "  subroutine symmetries_z\n\n";
print FILE_SYMMETRIES_Z "  use param\n";
print FILE_SYMMETRIES_Z "  use arrays\n\n";
print FILE_SYMMETRIES_Z "  implicit none\n\n";
print FILE_SYMMETRIES_Z "  integer j\n\n";
print FILE_SYMMETRIES_Z "  do j=1,ghost\n\n";

# Write beginning of file syncall.f90

print FILE_SYNCALL "! Automatically generated file.  Do not edit!\n\n";
print FILE_SYNCALL "  subroutine syncall\n\n";
print FILE_SYNCALL "  use param\n";
print FILE_SYNCALL "  use arrays\n\n";
print FILE_SYNCALL "  implicit none\n\n";

# Write beginning of file update.f90

print FILE_UPDATE "! Automatically generated file.  Do not edit!\n\n";
print FILE_UPDATE "  subroutine update(dtw)\n\n";
print FILE_UPDATE "  use param\n";
print FILE_UPDATE "  use arrays\n\n";
print FILE_UPDATE "  implicit none\n\n";
print FILE_UPDATE "  real(8) dtw \n\n";

# Write beginning of file bound_interp.inc

print FILE_BOUNDINTERP "! Automatically generated file.  Do not edit!\n\n";

# Write beginning of file restrict_interp.inc

print FILE_RESTRICTINTERP "! Automatically generated file.  Do not edit!\n\n";

# Write beginning of file intersect_interp.inc

print FILE_INTERSECTINTERP "! Automatically generated file.  Do not edit!\n\n";

# Open infile arrays.f90

open(INFILE,"src/base/arrays.config") or die "Can't open arrays.config: $!";

# Parse file arrays.f90 to identify declared arrays

while ($line=<INFILE>) {

   $nline = $nline+1;

#  Look only for lines declaring real or complex arrays, ignore all other lines.

   if (($line =~ /^\s*REAL/i)||($line =~ /^\s*COMPLEX/i)) {

#     Check that all keywords are present and grab array name (make sure to
#     ignore possible comment at the end).

      if ($line =~ /^\s+REAL/i) {
          $type = "real(8)";
          $zero = "0.d0";
          if ($line =~ /REAL\s*(\S+)\s*!\s*SYMMETRYR\s*=\s*(\S+)\s*,\s*SYMMETRYZ\s*=\s*(\S+)\s*,\s*INTENT\s*=\s*(\S+)\s*,\s*STORAGE\s*=\s*(.+)/i) {
             $var = $1;
             $symr = $2;
             $symz = $3;
             $intent = $4;
             $storage = $5;
          } else {
             die "arrays.pl: Bad syntax for REAL array assignment in line ",$nline," of file arrays.config\n\n";
          }
      }

#     Check if we have a 0D array.

      if ($line =~ /ZEROD/i) {
	  $zerod = "true";
      } else {
	  $zerod = "false";
      }

#     Check if the array has only one grid level.

      if ($line =~ /ONELEVEL/i) {
	  $onelevel = "true";
      } else {
	  $onelevel = "false";
      }

#     Check if the array does not require boundary conditions.

      if ($line =~ /NOBOUND/i) {
	  $nobound = "true";
      } else {
	  $nobound = "false";
      }

#     Check if the array should be checkpointed.

      if ($line =~ /CHECKPOINT/i) {
	  $checkpoint = "true";
      } else {
	  $checkpoint = "false";
      }

#     Check if the array does not require synchronization.

      if ($line =~ /NOSYNC/i) {
	  $nosync = "true";
      } else {
	  $nosync = "false";
      }

#     Check for commas.

      if ($var =~ /,/) {
          die "arrays.pl: Bad array assignment in line ",$nline," of file arrays.config\n",
              "           Only one array can be declared per line\n\n";
      }

#     Now write to FILE_ARRAYS code to declare arrays.

      if ($zerod eq "false") {

         if ($intent =~ /^POINTER$/i) {
            if ($onelevel eq "true") {
               print FILE_ARRAYS  "  ",$type,", pointer :: ",$var,"\n";
            } else {
               print FILE_ARRAYS  "  ",$type,", dimension (:,:), pointer :: ",$var,"\n";
            }
         } else {
            if ($onelevel eq "true") {
               print FILE_ARRAYS  "  ",$type,", pointer :: ",$var,"\n";
            } else {
               print FILE_ARRAYS  "  ",$type,", dimension (:,:), pointer :: ",$var,"\n";
            }
         }

         if (($intent =~ /^EVOLVE$/i) || ($intent =~ /^ELLIPTIC$/i)) {
            print FILE_ARRAYS  "  ",$type,", dimension (:,:), pointer :: s",$var,"\n";
            print FILE_ARRAYS  "  ",$type,", dimension (:,:), pointer :: ",$var,"_p\n";
            print FILE_ARRAYS  "  ",$type,", dimension (:,:), pointer :: ",$var,"_a\n";
            print FILE_ARRAYS  "  ",$type,", dimension (:,:), pointer :: ",$var,"_bound_rL(:,:,:)\n";
            print FILE_ARRAYS  "  ",$type,", dimension (:,:), pointer :: ",$var,"_bound_rR(:,:,:)\n";
            print FILE_ARRAYS  "  ",$type,", dimension (:,:), pointer :: ",$var,"_bound_zL(:,:,:)\n";
            print FILE_ARRAYS  "  ",$type,", dimension (:,:), pointer :: ",$var,"_bound_zR(:,:,:)\n";
         }

      } elsif ($zerod eq "true") {

         print FILE_ARRAYS  "  ",$type," :: ",$var,"\n";

         if ($intent =~ /^EVOLVE$/i) {
            print FILE_ARRAYS  "  ",$type," :: ",$var,"_p\n";
            print FILE_ARRAYS  "  ",$type," :: ",$var,"_a\n";
            print FILE_ARRAYS  "  ",$type," :: s",$var,"\n";
         }

      }

#     Now write to FILE_MYTYPES code to declare arrays.

      if ($zerod eq "false") {

         if ($intent !~ /^POINTER$/i) {
            if ($onelevel eq "true") {
               print FILE_MYTYPES  "     ",$type,", pointer :: ",$var,"\n";
            } else {
               print FILE_MYTYPES  "     ",$type,", dimension (:,:), pointer :: ",$var,"\n";
            }
         }

         if (($intent =~ /^EVOLVE$/i) || ($intent =~ /^ELLIPTIC$/i)) {
            print FILE_MYTYPES  "     ",$type,", dimension (:,:), pointer :: s",$var,"\n";
            print FILE_MYTYPES  "     ",$type,", dimension (:,:), pointer :: ",$var,"_p\n";
            print FILE_MYTYPES  "     ",$type,", dimension (:,:), pointer :: ",$var,"_a\n";
            print FILE_MYTYPES  "     ",$type,", dimension (:,:,:), pointer :: ",$var,"_bound_rL\n";
            print FILE_MYTYPES  "     ",$type,", dimension (:,:,:), pointer :: ",$var,"_bound_rR\n";
            print FILE_MYTYPES  "     ",$type,", dimension (:,:,:), pointer :: ",$var,"_bound_zL\n";
            print FILE_MYTYPES  "     ",$type,", dimension (:,:,:), pointer :: ",$var,"_bound_zR\n";
         }

      } elsif ($zerod eq "true") {

         print FILE_MYTYPES  "     ",$type," :: ",$var,"\n";

         if ($intent =~ /^EVOLVE$/i) {
            print FILE_MYTYPES  "     ",$type," :: ",$var,"_p\n";
            print FILE_MYTYPES  "     ",$type," :: ",$var,"_a\n";
            print FILE_MYTYPES  "     ",$type," :: s",$var,"\n";
         }

      }

#     Now write to FILE_CURRENTGRID code to point to arrays.

      if ($zerod eq "false") {

         if ($intent !~ /^POINTER$/i) {
            if ($onelevel eq "true") {
            } elsif (($intent !~ /^EVOLVE$/i) && ($intent !~ /^ELLIPTIC$/i)) {
               print FILE_CURRENTGRID  "  ",$var," => localgrid%",$var,"\n\n";
            } else {
               print FILE_CURRENTGRID  "  ",$var,"   => localgrid%",$var,"\n";
               print FILE_CURRENTGRID  "  s",$var,"  => localgrid%s",$var,"\n";
               print FILE_CURRENTGRID  "  ",$var,"_p => localgrid%",$var,"_p\n";
               print FILE_CURRENTGRID  "  ",$var,"_a => localgrid%",$var,"_a\n";
               print FILE_CURRENTGRID  "  ",$var,"_bound_rL => localgrid%",$var,"_bound_rL\n";
               print FILE_CURRENTGRID  "  ",$var,"_bound_rR => localgrid%",$var,"_bound_rR\n";
               print FILE_CURRENTGRID  "  ",$var,"_bound_zL => localgrid%",$var,"_bound_zL\n";
               print FILE_CURRENTGRID  "  ",$var,"_bound_zR => localgrid%",$var,"_bound_zR\n\n";
	    }
	 }

      }

#     Write to FILE_ACCUMULATE code to save old variables.

      if ($intent =~ /EVOLVE/i) {
         if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
            $cond = $1;
            print FILE_ACCUMULATE  "  if (",$cond,") then\n";
            print FILE_ACCUMULATE  "     if (k==1) then\n";
            print FILE_ACCUMULATE  "        ",$var,"_a = w*s",$var,"\n";
            print FILE_ACCUMULATE  "     else if (k<niter) then\n";
            print FILE_ACCUMULATE  "        ",$var,"_a = ",$var,"_a + w*s",$var,"\n";
            print FILE_ACCUMULATE  "     else\n";
            print FILE_ACCUMULATE  "        s",$var,"  = ",$var,"_a + w*s",$var,"\n";
            print FILE_ACCUMULATE  "     end if\n";
            print FILE_ACCUMULATE  "  end if\n\n";
         } else {
            print FILE_ACCUMULATE  "  if (k==1) then\n";
            print FILE_ACCUMULATE  "     ",$var,"_a = w*s",$var,"\n";
            print FILE_ACCUMULATE  "  else if (k<niter) then\n";
            print FILE_ACCUMULATE  "     ",$var,"_a = ",$var,"_a + w*s",$var,"\n";
            print FILE_ACCUMULATE  "  else\n";
            print FILE_ACCUMULATE  "     s",$var,"  = ",$var,"_a + w*s",$var,"\n";
            print FILE_ACCUMULATE  "  end if\n\n";
         }
      }

#     Write to FILE_ALLOCATEARRAYS code to allocate memory.

      if ($zerod eq "false") {

         if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
             $cond    = $1;
             $space   = "   ";
             $newline = "";
             print FILE_ALLOCATEARRAYS  "  if (",$cond,") then\n";
         } elsif ($storage !~ /^ALWAYS/i) {
             print $storage,"\n";
             die "arrays.pl: Bad STORAGE assignment in line ",$nline," of file arrays.config\n\n";
         } else {
             $space   = "";
             $newline = "\n";
         }

         if ($intent =~/^OUTPUT$/i) {
            print FILE_ALLOCATEARRAYS  $space,"  if (contains(outvars0D,\"",$var,"\").or. &\n";
            print FILE_ALLOCATEARRAYS  $space,"      contains(outvars1D,\"",$var,"\").or. &\n";
            print FILE_ALLOCATEARRAYS  $space,"      contains(outvars2D,\"",$var,"\")) then\n";
            print FILE_ALLOCATEARRAYS  $space,"     if (trim(status)=='on') then\n";
            print FILE_ALLOCATEARRAYS  $space,"     do box=0,Nb\n";
            print FILE_ALLOCATEARRAYS  $space,"        do level=min(1,box),Nl(box)\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"(1-ghost:Nrmaxl(box),1-ghost:Nzmaxl(box)))\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var," = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"        end do\n";
            print FILE_ALLOCATEARRAYS  $space,"     end do\n";
            print FILE_ALLOCATEARRAYS  $space,"     end if\n";
            print FILE_ALLOCATEARRAYS  $space,"  end if\n",$newline;
         } elsif ($intent =~ /^EVOLVE$/i) {
            print FILE_ALLOCATEARRAYS  $space,"  if (trim(status)=='on') then\n";
            print FILE_ALLOCATEARRAYS  $space,"     checkvars = trim(checkvars) // ',",$var,"'\n";
            print FILE_ALLOCATEARRAYS  $space,"     do box=0,Nb\n";
            print FILE_ALLOCATEARRAYS  $space,"        do level=min(1,box),Nl(box)\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"(1-ghost:Nrmaxl(box),1-ghost:Nzmaxl(box)))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%s",$var,"(1-ghost:Nrmaxl(box),1-ghost:Nzmaxl(box)))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"_p(1-ghost:Nrmaxl(box),1-ghost:Nzmaxl(box)))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"_a(1-ghost:Nrmaxl(box),1-ghost:Nzmaxl(box)))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"_bound_rL(0:ghost-1,1-ghost:Nzmaxl(box),0:2))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"_bound_rR(0:ghost-1,1-ghost:Nzmaxl(box),0:2))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"_bound_zL(1-ghost:Nrmaxl(box),0:ghost-1,0:2))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"_bound_zR(1-ghost:Nrmaxl(box),0:ghost-1,0:2))\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var," = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%s",$var," = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var,"_p = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var,"_a = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var,"_bound_rL = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var,"_bound_rR = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var,"_bound_zL = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var,"_bound_zR = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"        end do\n";
            print FILE_ALLOCATEARRAYS  $space,"     end do\n";
            print FILE_ALLOCATEARRAYS  $space,"  end if\n",$newline;
         } elsif ($intent =~ /^ELLIPTIC$/i) {
            print FILE_ALLOCATEARRAYS  $space,"  if (trim(status)=='on') then\n";
            print FILE_ALLOCATEARRAYS  $space,"     do box=0,Nb\n";
            print FILE_ALLOCATEARRAYS  $space,"        do level=min(1,box),Nl(box)\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"(1-ghost:Nrmaxl(box),1-ghost:Nzmaxl(box)))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%s",$var,"(1-ghost:Nrmaxl(box),1-ghost:Nzmaxl(box)))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"_p(1-ghost:Nrmaxl(box),1-ghost:Nzmaxl(box)))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"_a(1-ghost:Nrmaxl(box),1-ghost:Nzmaxl(box)))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"_bound_rL(0:ghost-1,1-ghost:Nzmaxl(box),0:2))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"_bound_rR(0:ghost-1,1-ghost:Nzmaxl(box),0:2))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"_bound_zL(1-ghost:Nrmaxl(box),0:ghost-1,0:2))\n";
            print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"_bound_zR(1-ghost:Nrmaxl(box),0:ghost-1,0:2))\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var," = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%s",$var," = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var,"_p = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var,"_a = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var,"_bound_rL = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var,"_bound_rR = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var,"_bound_zL = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var,"_bound_zR = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"        end do\n";
            print FILE_ALLOCATEARRAYS  $space,"     end do\n";
            print FILE_ALLOCATEARRAYS  $space,"  end if\n",$newline;
         } elsif ($intent =~ /^AUXILIARY$/i) {
            print FILE_ALLOCATEARRAYS  $space,"  if (trim(status)=='on') then\n";
            if ($checkpoint eq "true") {
               print FILE_ALLOCATEARRAYS  $space,"     checkvars = trim(checkvars) // ',",$var,"'\n";
            }
            if ($onelevel eq "true") {
            } else {
               print FILE_ALLOCATEARRAYS  $space,"     do box=0,Nb\n";
               print FILE_ALLOCATEARRAYS  $space,"        do level=min(1,box),Nl(box)\n";
               print FILE_ALLOCATEARRAYS  $space,"           allocate(grid(box,level)%",$var,"(1-ghost:Nrmaxl(box),1-ghost:Nzmaxl(box)))\n";
               print FILE_ALLOCATEARRAYS  $space,"           grid(box,level)%",$var," = ",$zero,"\n";
               print FILE_ALLOCATEARRAYS  $space,"        end do\n";
               print FILE_ALLOCATEARRAYS  $space,"     end do\n";
            }
            print FILE_ALLOCATEARRAYS  $space,"  end if\n",$newline;;
         } elsif ($intent =~ /^POINTER$/i) {
         } else {
            die "arrays.pl: Bad INTENT assignment in line ",$nline," of file arrays.config\n\n";
         }

         if ($storage =~ /^CONDITIONAL\s*\(.*\)/i) {
            print FILE_ALLOCATEARRAYS  "  else if (contains(outvars0D,\"",$var,"\").or. &\n";
            print FILE_ALLOCATEARRAYS  "           contains(outvars1D,\"",$var,"\").or. &\n";
            print FILE_ALLOCATEARRAYS  "           contains(outvars2D,\"",$var,"\")) then\n";
            print FILE_ALLOCATEARRAYS  "     if (rank==0) then\n";
            print FILE_ALLOCATEARRAYS  "        print *\n";
            print FILE_ALLOCATEARRAYS  "        print *, 'Error in parfile: array ",$var," has no storage,'\n";
            print FILE_ALLOCATEARRAYS  "        print *, 'so no output for it is possible.'\n";
            print FILE_ALLOCATEARRAYS  "        print *\n";
            print FILE_ALLOCATEARRAYS  "        print *, 'Aborting! (subroutine allocatearrays.f90)'\n";
            print FILE_ALLOCATEARRAYS  "        print *\n";
            print FILE_ALLOCATEARRAYS  "     end if\n";
            print FILE_ALLOCATEARRAYS  "     call die\n";
            print FILE_ALLOCATEARRAYS  "  end if\n\n";
         }

      }

#     Now write to FILE_GRABARRAY code to compare array name.

      if ($intent !~ /^POINTER$/i) {
         print FILE_GRABARRAY "  if (varname=='",$var,"') then\n";
         print FILE_GRABARRAY "     exists = .true.\n";
         print FILE_GRABARRAY "     grabvar => ",$var,"\n";
         if ($intent =~ /^EVOLVE$/i) {
            print FILE_GRABARRAY "     grabvar_bound_rL => ",$var,"_bound_rL\n";
            print FILE_GRABARRAY "     grabvar_bound_rR => ",$var,"_bound_rR\n";
            print FILE_GRABARRAY "     grabvar_bound_zL => ",$var,"_bound_zL\n";
            print FILE_GRABARRAY "     grabvar_bound_zR => ",$var,"_bound_zR\n";
         }
         print FILE_GRABARRAY "  end if\n\n";
      }

      if ($intent =~ /^EVOLVE$/i) {
         print FILE_GRABARRAY "  if (varname=='s",$var,"') then\n";
         print FILE_GRABARRAY "     exists = .true.\n";
         print FILE_GRABARRAY "     grabvar => s",$var,"\n";
         print FILE_GRABARRAY "  end if\n\n";
      }

#     Write to FILE_SAVEOLD code to save old variables.

      if ($intent =~ /EVOLVE/i) {
         if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
            $cond = $1;
            print FILE_SAVEOLD  "  if (",$cond,") then\n";
            print FILE_SAVEOLD  "     ",$var,"_p = ",$var,"\n";
            print FILE_SAVEOLD  "     do i=0,ghost-1\n";
            print FILE_SAVEOLD  "        ",$var,"_bound_rL(i,:,2) = ",$var,"_bound_rL(i,:,1)\n";
            print FILE_SAVEOLD  "        ",$var,"_bound_rL(i,:,1) = ",$var,"_bound_rL(i,:,0)\n";
            print FILE_SAVEOLD  "        ",$var,"_bound_rL(i,:,0) = ",$var,"(1-ghost+i,:)\n";
            print FILE_SAVEOLD  "        ",$var,"_bound_rR(i,:,2) = ",$var,"_bound_rR(i,:,1)\n";
            print FILE_SAVEOLD  "        ",$var,"_bound_rR(i,:,1) = ",$var,"_bound_rR(i,:,0)\n";
            print FILE_SAVEOLD  "        ",$var,"_bound_rR(i,:,0) = ",$var,"(Nr-i,:)\n";
            print FILE_SAVEOLD  "        ",$var,"_bound_zL(:,i,2) = ",$var,"_bound_zL(:,i,1)\n";
            print FILE_SAVEOLD  "        ",$var,"_bound_zL(:,i,1) = ",$var,"_bound_zL(:,i,0)\n";
            print FILE_SAVEOLD  "        ",$var,"_bound_zL(:,i,0) = ",$var,"(:,1-ghost+i)\n";
            print FILE_SAVEOLD  "        ",$var,"_bound_zR(:,i,2) = ",$var,"_bound_zR(:,i,1)\n";
            print FILE_SAVEOLD  "        ",$var,"_bound_zR(:,i,1) = ",$var,"_bound_zR(:,i,0)\n";
            print FILE_SAVEOLD  "        ",$var,"_bound_zR(:,i,0) = ",$var,"(:,Nz-i)\n";
            print FILE_SAVEOLD  "     end do\n";
            print FILE_SAVEOLD  "  end if\n\n";
         } else {
            print FILE_SAVEOLD  "  ",$var,"_p = ",$var,"\n";
            print FILE_SAVEOLD  "  do i=0,ghost-1\n";
            print FILE_SAVEOLD  "     ",$var,"_bound_rL(i,:,2) = ",$var,"_bound_rL(i,:,1)\n";
            print FILE_SAVEOLD  "     ",$var,"_bound_rL(i,:,1) = ",$var,"_bound_rL(i,:,0)\n";
            print FILE_SAVEOLD  "     ",$var,"_bound_rL(i,:,0) = ",$var,"(1-ghost+i,:)\n";
            print FILE_SAVEOLD  "     ",$var,"_bound_rR(i,:,2) = ",$var,"_bound_rR(i,:,1)\n";
            print FILE_SAVEOLD  "     ",$var,"_bound_rR(i,:,1) = ",$var,"_bound_rR(i,:,0)\n";
            print FILE_SAVEOLD  "     ",$var,"_bound_rR(i,:,0) = ",$var,"(Nr-i,:)\n";
            print FILE_SAVEOLD  "     ",$var,"_bound_zL(:,i,2) = ",$var,"_bound_zL(:,i,1)\n";
            print FILE_SAVEOLD  "     ",$var,"_bound_zL(:,i,1) = ",$var,"_bound_zL(:,i,0)\n";
            print FILE_SAVEOLD  "     ",$var,"_bound_zL(:,i,0) = ",$var,"(:,1-ghost+i)\n";
            print FILE_SAVEOLD  "     ",$var,"_bound_zR(:,i,2) = ",$var,"_bound_zR(:,i,1)\n";
            print FILE_SAVEOLD  "     ",$var,"_bound_zR(:,i,1) = ",$var,"_bound_zR(:,i,0)\n";
            print FILE_SAVEOLD  "     ",$var,"_bound_zR(:,i,0) = ",$var,"(:,Nz-i)\n";
            print FILE_SAVEOLD  "  end do\n\n";
         }
      }

#     Write to FILE_SIMPLEBOUNDARIES code to apply simple boundary conditions.

      if ($zerod eq "false") {

         if (($intent =~ /EVOLVE/i) && ($nobound eq "false")) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
                $cond = $1;
                print FILE_SIMPLEBOUNDARY  "  if (",$cond,") then\n";
	        print FILE_SIMPLEBOUNDARY  "     if (boundtype=='static') then\n";
	        print FILE_SIMPLEBOUNDARY  "        s",$var,"(Nr,:) = 0.d0\n";
	        print FILE_SIMPLEBOUNDARY  "        s",$var,"(:,Nz) = 0.d0\n";
	        print FILE_SIMPLEBOUNDARY  "        if (.not.eqsym) then\n";
	        print FILE_SIMPLEBOUNDARY  "           s",$var,"(:,1-ghost) = 0.d0\n";
	        print FILE_SIMPLEBOUNDARY  "        end if\n";
	        print FILE_SIMPLEBOUNDARY  "     else if (boundtype=='flat') then\n";
                print FILE_SIMPLEBOUNDARY  "        s",$var,"(Nr,:) = s",$var,"(Nr-1,:)\n";
                print FILE_SIMPLEBOUNDARY  "        s",$var,"(:,Nz) = s",$var,"(:,Nz-1)\n";
	        print FILE_SIMPLEBOUNDARY  "        if (.not.eqsym) then\n";
	        print FILE_SIMPLEBOUNDARY  "           s",$var,"(:,1-ghost) = s",$var,"(:,2-ghost)\n";
	        print FILE_SIMPLEBOUNDARY  "        end if\n";
                print FILE_SIMPLEBOUNDARY  "     end if\n";
                print FILE_SIMPLEBOUNDARY  "  end if\n";
            } else {
	        print FILE_SIMPLEBOUNDARY  "  if (boundtype=='static') then\n";
	        print FILE_SIMPLEBOUNDARY  "     s",$var,"(Nr,:) = 0.d0\n";
	        print FILE_SIMPLEBOUNDARY  "     s",$var,"(:,Nz) = 0.d0\n";
	        print FILE_SIMPLEBOUNDARY  "     if (.not.eqsym) then\n";
	        print FILE_SIMPLEBOUNDARY  "        s",$var,"(:,1-ghost) = 0.d0\n";
	        print FILE_SIMPLEBOUNDARY  "     end if\n";
	        print FILE_SIMPLEBOUNDARY  "  else if (boundtype=='flat') then\n";
	        print FILE_SIMPLEBOUNDARY  "     s",$var,"(Nr,:) = s",$var,"(Nr-1,:)\n";
                print FILE_SIMPLEBOUNDARY  "     s",$var,"(:,Nz) = s",$var,"(:,Nz-1)\n";
	        print FILE_SIMPLEBOUNDARY  "     if (.not.eqsym) then\n";
	        print FILE_SIMPLEBOUNDARY  "        s",$var,"(:,1-ghost) = s",$var,"(:,2-ghost)\n";
	        print FILE_SIMPLEBOUNDARY  "     end if\n";
                print FILE_SIMPLEBOUNDARY  "  end if\n";
            }
         }

      }

#     Write to FILE_SYMMETRIES_R code to apply symmetries at axis.

      if ($zerod eq "false") {

         if ($intent =~ /EVOLVE/i) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
               $cond = $1;
               print FILE_SYMMETRIES_R  "     if (",$cond,") then\n";
               if ($symr == "+1") {
                  print FILE_SYMMETRIES_R  "        ",$var,"(1-i,:) = + ",$var,"(i,:)\n";
               } elsif ($symr == "-1") {
                  print FILE_SYMMETRIES_R  "        ",$var,"(1-i,:) = - ",$var,"(i,:)\n";
               } else {
                  print FILE_SYMMETRIES_R  "        ",$var,"(1-i,:) = ",$symr,"*",$var,"(i,:)\n";
               }
               print FILE_SYMMETRIES_R  "     end if\n\n";
            } else {
               if ($symr == "+1") {
	          print FILE_SYMMETRIES_R  "     ",$var,"(1-i,:) = + ",$var,"(i,:)\n\n";
               } elsif ($symr == "-1") {
                  print FILE_SYMMETRIES_R  "     ",$var,"(1-i,:) = - ",$var,"(i,:)\n\n";
               } elsif ($symr != "0") {
                  print FILE_SYMMETRIES_R  "     ",$var,"(1-i,:) = ",$symr,"*",$var,"(i,:)\n\n";
               }
            }
         }

      }

#     Write to FILE_SYMMETRIES_Z code to apply symmetries at equator.

      if ($zerod eq "false") {

         if ($intent =~ /EVOLVE/i) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
               $cond = $1;
               print FILE_SYMMETRIES_Z  "     if (",$cond,") then\n";
               if ($symz == "+1") {
                  print FILE_SYMMETRIES_Z  "        ",$var,"(:,1-j) = + ",$var,"(:,j)\n";
               } elsif ($symz == "-1") {
                  print FILE_SYMMETRIES_Z  "        ",$var,"(:,1-j) = - ",$var,"(:,j)\n";
               } else {
                  print FILE_SYMMETRIES_Z  "        ",$var,"(:,1-j) = ",$symz,"*",$var,"(:,j)\n";
               }
               print FILE_SYMMETRIES_Z  "     end if\n\n";
            } else {
               if ($symz == "+1") {
	          print FILE_SYMMETRIES_Z  "     ",$var,"(:,1-j) = + ",$var,"(:,j)\n\n";
               } elsif ($symz == "-1") {
                  print FILE_SYMMETRIES_Z  "     ",$var,"(:,1-j) = - ",$var,"(:,j)\n\n";
               } elsif ($symr != "0") {
                  print FILE_SYMMETRIES_Z  "     ",$var,"(:,1-j) = ",$symz,"*",$var,"(:,j)\n\n";
               }
            }
         }

      }

#     Write to FILE_SYNCALL code to synchronize across processors.

      if ($zerod eq "false") {

         if ($intent =~ /EVOLVE/i) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
                $cond = $1;
                print FILE_SYNCALL  "  if (",$cond,") then\n";
	        print FILE_SYNCALL  "     call sync(",$var,")\n";
                print FILE_SYNCALL  "  end if\n\n";
            } else {
	        print FILE_SYNCALL  "  call sync(",$var,")\n\n";
            }
         }

      }

#     Write to FILE_UPDATE code to update variables.

      if ($intent =~ /EVOLVE/i) {
         if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
            $cond = $1;
            print FILE_UPDATE  "  if (",$cond,") then\n";
            print FILE_UPDATE  "     ",$var," = ",$var,"_p + dtw*s",$var,"\n";
            print FILE_UPDATE  "  end if\n\n";
         } else {
            print FILE_UPDATE  "  ",$var," = ",$var,"_p + dtw*s",$var,"\n\n";
         }
      }

#     Write to FILE_BOUNDINTERP code to interpolate variables at boundaries.

      if ($intent =~ /EVOLVE/i) {
      #if (($intent =~ /EVOLVE/i) && ($nobound eq "false")) {
         if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
            print FILE_BOUNDINTERP  "  if (",$cond,") then\n";
            print FILE_BOUNDINTERP  "     interpvar => grid(bbox,level-1)%",$var,"\n";
            print FILE_BOUNDINTERP  "     aux1 = interp(bbox,level-1,r0,z0,flag2)\n";
            print FILE_BOUNDINTERP  "     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)\n";
            print FILE_BOUNDINTERP  "     if (flag1) then\n";
            print FILE_BOUNDINTERP  "        if (border==1) then\n";
            print FILE_BOUNDINTERP  "           ",$var,"(i,j) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*",$var,"_p(i,j)\n";
            print FILE_BOUNDINTERP  "        else\n";
            print FILE_BOUNDINTERP  "           if (bd=='rL') then\n";
            print FILE_BOUNDINTERP  "              ",$var,"(i,j) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &\n";
            print FILE_BOUNDINTERP  "              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*",$var,"_bound_rL(m,j,0) &\n";
            print FILE_BOUNDINTERP  "              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*",$var,"_bound_rL(m,j,1)\n";
            print FILE_BOUNDINTERP  "           else if (bd=='rR') then\n";
            print FILE_BOUNDINTERP  "              ",$var,"(i,j) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &\n";
            print FILE_BOUNDINTERP  "              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*",$var,"_bound_rR(m,j,0) &\n";
            print FILE_BOUNDINTERP  "              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*",$var,"_bound_rR(m,j,1)\n";
            print FILE_BOUNDINTERP  "           else if (bd=='zL') then\n";
            print FILE_BOUNDINTERP  "              ",$var,"(i,j) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &\n";
            print FILE_BOUNDINTERP  "              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*",$var,"_bound_zL(i,n,0) &\n";
            print FILE_BOUNDINTERP  "              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*",$var,"_bound_zL(i,n,1)\n";
            print FILE_BOUNDINTERP  "           else if (bd=='zR') then\n";
            print FILE_BOUNDINTERP  "              ",$var,"(i,j) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &\n";
            print FILE_BOUNDINTERP  "              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*",$var,"_bound_zR(i,n,0) &\n";
            print FILE_BOUNDINTERP  "              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*",$var,"_bound_zR(i,n,1)\n";
            print FILE_BOUNDINTERP  "           end if\n";
            print FILE_BOUNDINTERP  "        end if\n";
            print FILE_BOUNDINTERP  "     end if\n";
            print FILE_BOUNDINTERP  "  end if\n\n";
         } else {
            print FILE_BOUNDINTERP  "  interpvar => grid(bbox,level-1)%",$var,"\n";
            print FILE_BOUNDINTERP  "  aux1 = interp(bbox,level-1,r0,z0,flag2)\n";
            print FILE_BOUNDINTERP  "  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)\n";
            print FILE_BOUNDINTERP  "  if (flag1) then\n";
            print FILE_BOUNDINTERP  "     if (border==1) then\n";
            print FILE_BOUNDINTERP  "        ",$var,"(i,j) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*",$var,"_p(i,j)\n";
            print FILE_BOUNDINTERP  "     else\n";
            print FILE_BOUNDINTERP  "        if (bd=='rL') then\n";
            print FILE_BOUNDINTERP  "           ",$var,"(i,j) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &\n";
            print FILE_BOUNDINTERP  "           + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*",$var,"_bound_rL(m,j,0) &\n";
            print FILE_BOUNDINTERP  "           + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*",$var,"_bound_rL(m,j,1)\n";
            print FILE_BOUNDINTERP  "        else if (bd=='rR') then\n";
            print FILE_BOUNDINTERP  "           ",$var,"(i,j) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &\n";
            print FILE_BOUNDINTERP  "           + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*",$var,"_bound_rR(m,j,0) &\n";
            print FILE_BOUNDINTERP  "           + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*",$var,"_bound_rR(m,j,1)\n";
            print FILE_BOUNDINTERP  "        else if (bd=='zL') then\n";
            print FILE_BOUNDINTERP  "           ",$var,"(i,j) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &\n";
            print FILE_BOUNDINTERP  "           + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*",$var,"_bound_zL(i,n,0) &\n";
            print FILE_BOUNDINTERP  "           + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*",$var,"_bound_zL(i,n,1)\n";
            print FILE_BOUNDINTERP  "        else if (bd=='zR') then\n";
            print FILE_BOUNDINTERP  "           ",$var,"(i,j) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &\n";
            print FILE_BOUNDINTERP  "           + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*",$var,"_bound_zR(i,n,0) &\n";
            print FILE_BOUNDINTERP  "           + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*",$var,"_bound_zR(i,n,1)\n";
            print FILE_BOUNDINTERP  "        end if\n";
            print FILE_BOUNDINTERP  "     end if\n";
            print FILE_BOUNDINTERP  "  end if\n\n";
         }
      }

#     Write to FILE_RESTRICTCOPY code to restrict variables to coarse grid.

      if ($intent =~ /EVOLVE/i) {
         if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
            $cond = $1;
            print FILE_RESTRICTINTERP  "  if (",$cond,") then\n";
            print FILE_RESTRICTINTERP  "     interpvar => grid(box,level)%",$var,"\n";
            print FILE_RESTRICTINTERP  "     aux1 = interp(box,level,r0,z0,flag2)\n";
            print FILE_RESTRICTINTERP  "     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)\n";
            print FILE_RESTRICTINTERP  "     if (flag1) then\n";
            print FILE_RESTRICTINTERP  "        grid(bbox,level-1)%",$var,"(i,j) = aux2\n";
            print FILE_RESTRICTINTERP  "     end if\n";
            print FILE_RESTRICTINTERP  "  end if\n\n";
         } else {
            print FILE_RESTRICTINTERP  "  interpvar => grid(box,level)%",$var,"\n";
            print FILE_RESTRICTINTERP  "  aux1 = interp(box,level,r0,z0,flag2)\n";
            print FILE_RESTRICTINTERP  "  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)\n";
            print FILE_RESTRICTINTERP  "  if (flag1) then\n";
            print FILE_RESTRICTINTERP  "     grid(bbox,level-1)%",$var,"(i,j) = aux2\n";
            print FILE_RESTRICTINTERP  "  end if\n\n";
         }
      }

#    Write to FILE_INTERSECTCOPY code to interpolate on box intersections.

     if (($intent =~ /EVOLVE/i || $intent =~ /AUXILIARY/) && ($nosync eq "false")) {
         if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
            $cond = $1;
            print FILE_INTERSECTINTERP  "  if (",$cond,") then\n";
            print FILE_INTERSECTINTERP  "     interpvar => grid(b2,level)%",$var,"\n";
            print FILE_INTERSECTINTERP  "     aux1 = interp(b2,level,r0,z0,flag2)\n";
            print FILE_INTERSECTINTERP  "     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)\n";
            print FILE_INTERSECTINTERP  "     if (flag1) then\n";
            print FILE_INTERSECTINTERP  "        grid(b1,level)%",$var,"(i,j) = aux2\n";
            print FILE_INTERSECTINTERP  "     end if\n";
            print FILE_INTERSECTINTERP  "  end if\n\n";
         } else {
            print FILE_INTERSECTINTERP  "  interpvar => grid(b2,level)%",$var,"\n";
            print FILE_INTERSECTINTERP  "  aux1 = interp(b2,level,r0,z0,flag2)\n";
            print FILE_INTERSECTINTERP  "  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)\n";
            print FILE_INTERSECTINTERP  "  if (flag1) then\n";
            print FILE_INTERSECTINTERP  "     grid(b1,level)%",$var,"(i,j) = aux2\n";
            print FILE_INTERSECTINTERP  "  end if\n\n";
         }
      }

#  Close two main conditional statements.

   } elsif (($line !~ /^\s*!.*/)&&($line !~ /^\s*$/)) { 
      die "arrays.pl: Bad syntax in line ",$nline," of file arrays.config\n\n";
   }
}

# Close INFILE.

close(INFILE);

# Write ending of file arrays.f90.

print FILE_ARRAYS "\n  end module arrays\n\n";

# Write ending of file accumulate.f90.

print FILE_ACCUMULATE "  end subroutine accumulate\n\n";

# Write ending of file allocatearrays.f90

print FILE_ALLOCATEARRAYS "  end subroutine allocatearrays\n\n";

# Write ending of file currentgrid.f90.

print FILE_CURRENTGRID "  end subroutine currentgrid\n\n";

# Write ending of file grabarray.f90.

print FILE_GRABARRAY "  if (.not.exists) then\n";
print FILE_GRABARRAY "     if (rank==0) then\n";
print FILE_GRABARRAY "        print *\n";
print FILE_GRABARRAY "        print *, 'Error in parfile, non-existent array: ',varname\n";
print FILE_GRABARRAY "        print *, 'Aborting! (subroutine grabarray.f90)'\n";
print FILE_GRABARRAY "        print *\n";
print FILE_GRABARRAY "     end if\n";
print FILE_GRABARRAY "     call die\n";
print FILE_GRABARRAY "  end if\n\n";
print FILE_GRABARRAY "  end subroutine grabarray\n\n";

# Write ending of file mytypes.f90.

print FILE_MYTYPES "\n";
print FILE_MYTYPES "  end type gridfuncs\n\n";
print FILE_MYTYPES "  end module mytypes\n\n";

# Write ending of file saveold.f90.

print FILE_SAVEOLD "  end subroutine saveold\n\n";

# Write ending of file simpleboundary.f90

print FILE_SIMPLEBOUNDARY "\n";
print FILE_SIMPLEBOUNDARY "  end subroutine simpleboundary\n\n";

# Write ending of file symmetries_r.f90.

print FILE_SYMMETRIES_R  "  end do\n\n";
print FILE_SYMMETRIES_R  "  end subroutine symmetries_r\n\n";

# Write ending of file symmetries_z.f90.

print FILE_SYMMETRIES_Z  "  end do\n\n";
print FILE_SYMMETRIES_Z  "  end subroutine symmetries_z\n\n";

# Write ending of file syncall.f90.

print FILE_SYNCALL  "  end subroutine syncall\n\n";

# Write ending of file update.f90.

print FILE_UPDATE "  end subroutine update\n\n";

# Close output files.

close(FILE_ARRAYS);
close(FILE_ACCUMULATE);
close(FILE_ALLOCATEARRAYS);
close(FILE_CURRENTGRID);
close(FILE_GRABARRAY);
close(FILE_MYTYPES);
close(FILE_SAVEOLD);
close(FILE_SIMPLEBOUNDARY);
close(FILE_SYMMETRIES_R);
close(FILE_SYMMETRIES_Z);
close(FILE_SYNCALL);
close(FILE_UPDATE);

close(FILE_BOUNDINTERP);
close(FILE_RESTRICTINTERP);
close(FILE_INTERSECTINTERP);


