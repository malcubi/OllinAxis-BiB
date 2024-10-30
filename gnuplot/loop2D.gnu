
# Plot surface.  For a given variable, we read
# all files that have the extension '2D' and
# plot them on top of each other.

FILES = system("ls -1 harmonic/trK*2D")

splot for [data in FILES] data using 1:2:3 index i with lines

# Pause some time

pause 0.2

# Increment i.

i=i+1

# Check if we are done.

if (i<N) reread

