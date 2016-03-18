
set xrange[0:128]
set yrange[0:1]
plot "./output/c_0.dat" with lines
pause 1

plot "./output/c_0.dat" with lines
pause 0.3

plot "./output/c_200.dat" with lines
pause 0.3

plot "./output/c_400.dat" with lines
pause 0.3

plot "./output/c_600.dat" with lines
pause 0.3

plot "./output/c_800.dat" with lines
pause 0.3

plot "./output/c_1000.dat" with lines
pause 0.3

plot "./output/c_1200.dat" with lines
pause 0.3

plot "./output/c_1400.dat" with lines
pause 0.3

plot "./output/c_1600.dat" with lines
pause 0.3

plot "./output/c_1800.dat" with lines
pause 0.3

plot "./output/c_2000.dat" with lines
pause 0.3

plot "./output/c_2200.dat" with lines
pause 0.3

plot "./output/c_2400.dat" with lines
pause 0.3

plot "./output/c_2600.dat" with lines
pause 0.3

plot "./output/c_2800.dat" with lines
pause 0.3

plot "./output/c_3000.dat" with lines
pause 0.3

plot "./output/c_3200.dat" with lines
pause 0.3

plot "./output/c_3400.dat" with lines
pause 0.3

plot "./output/c_3600.dat" with lines
pause 0.3

plot "./output/c_3800.dat" with lines
pause 0.3

plot "./output/c_4000.dat" with lines
pause 0.3

plot "./output/c_4200.dat" with lines
pause 0.3

plot "./output/c_4400.dat" with lines
pause 0.3

plot "./output/c_4600.dat" with lines
pause 0.3

plot "./output/c_4800.dat" with lines
pause 0.3
set term png
set output "finalProfile.png"
replot
set term x11