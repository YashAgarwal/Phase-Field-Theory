set yrange[0:1]
set xrange[0:512]
plot "./output/c_0.dat" with lines
pause 0.5
plot "./output/c_100.dat" with lines
pause 0.5
plot "./output/c_200.dat" with lines
pause 0.5
plot "./output/c_300.dat" with lines
pause 0.5
plot "./output/c_400.dat" with lines
pause 0.5
plot "./output/c_500.dat" with lines
pause 0.5
plot "./output/c_600.dat" with lines
pause 0.5
plot "./output/c_700.dat" with lines
pause 0.5
plot "./output/c_800.dat" with lines
pause 0.5
plot "./output/c_900.dat" with lines
pause 0.5
plot "./output/c_1000.dat" with lines
pause 0.5
plot "./output/c_1100.dat" with lines
pause 0.5
plot "./output/c_1200.dat" with lines
pause 0.5
plot "./output/c_1300.dat" with lines
pause 0.5
plot "./output/c_1400.dat" with lines
pause 0.5
plot "./output/c_1500.dat" with lines
pause 0.5
plot "./output/c_1600.dat" with lines
pause 0.5
plot "./output/c_1700.dat" with lines
pause 0.5
plot "./output/c_1800.dat" with lines
pause 0.5
plot "./output/c_1900.dat" with lines
pause 0.5
plot "./output/c_2000.dat" with lines
pause 0.5
set term png
set output "finalProfile.png"
replot
set term x11