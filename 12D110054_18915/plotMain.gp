set size square
set xrange[0:40]
set yrange[0:40]
unset colorbox
unset key; unset tics; unset border

set multiplot layout 1, 2
set tmargin 2
set title "before"
unset key
plot "./initialProfile.dat" matrix with image
set title "after"
unset key
plot "./finalProfile.dat" matrix with image
unset multiplot
