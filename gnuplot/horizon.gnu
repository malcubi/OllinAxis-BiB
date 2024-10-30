
# Set terminal to qt (cross-platform interactive)

set terminal qt

# Polar coordinates

set polar

# Set box range.  Change upper limit as desired.

set xrange [0:3]

set yrange [ 0:3]  # Equatorial symmetry.
#set yrange [-3:3]  # No equatorial symmetry.

# Plot horizon, but remember that our theta coordinate
# is measured from the axis, not the equator.

plot "ah_radius0.tl" using (pi/2-$1):2 with lines   # Equatorial symmetry.
#plot "ah_radius0.tl" using (pi-$1):2 with lines     # No equatorial symmetry.



