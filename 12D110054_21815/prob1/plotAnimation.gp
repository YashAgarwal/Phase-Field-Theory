set yrange[0:1]
set xrange[0:512]
plot "./output/phi_0.dat" with lines
pause 0.5
plot "./output/phi_200.dat" with lines
pause 0.5
plot "./output/phi_400.dat" with lines
pause 0.5
plot "./output/phi_600.dat" with lines
pause 0.5
plot "./output/phi_800.dat" with lines
pause 0.5
plot "./output/phi_1000.dat" with lines
pause 0.5
plot "./output/phi_1200.dat" with lines
pause 0.5
plot "./output/phi_1400.dat" with lines
pause 0.5
plot "./output/phi_1600.dat" with lines
pause 0.5
plot "./output/phi_1800.dat" with lines
pause 0.5
plot "./output/phi_2000.dat" with lines
pause 0.5
set term png
set output "finalProfile.png"
replot
set term x11