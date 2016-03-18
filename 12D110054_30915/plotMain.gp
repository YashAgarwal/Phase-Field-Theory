set size square
set xrange[0:40]
set yrange[0:40]
unset colorbox
unset key; unset tics; unset border
set autoscale xfix
set autoscale yfix
set autoscale cbfix

set multiplot layout 1, 3
set tmargin 2

set title "Initial Profile"
unset key
plot "./finalProfile.dat" matrix with image

set title "HK stage I"
unset key
plot "./HKstage1.dat" matrix  with image

set title "HK Final Stage"
unset key
plot "./HKfinal.dat" matrix with image

unset multiplot
