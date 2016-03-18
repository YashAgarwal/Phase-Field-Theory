set xtics 0,25000,100000
set tics font "Times-Roman,8"
set multiplot layout 2, 2 title "Analysis Report" font "Times-Roman,14"
set tmargin 2
set title "Interface Position Vs Time"
unset key
plot "./interfacePosition.dat" using 2:1 with line
set title "Precipitate length Vs Time"
unset key
plot "./pptLength.dat" using 2:1 with line
set title "Interface Velocity Vs Time"
unset key
plot "./velocity.dat" using 2:1 with line
unset multiplot