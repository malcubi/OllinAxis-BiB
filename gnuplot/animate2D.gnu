
# Set terminal to qt (cross-platform interactive)

set terminal qt

# Find number of data sets using a perl call, and save it in N.
# The call below first lists find the first file with extension '2D'
# in the directory 'namedir':
#
# ls namedir/*2D | head -1
#
# It then pipes the result into a perl command using xargs to
# count the number of data sets.

N = `ls harmonic/*2D | head -1 | xargs perl -ne '$count += tr/#//; END{print "$count\n"}'`

# Set range in z. Best if you first look at the data in 1D
# to figure out the correct range.

set zrange [0.9:1.1]
set ticslevel 0

# Set viewing angle.

set view 60,300,1,1

# Initialize loop variable i.

i=0

# Load looper code.

load "loop2D.gnu"
