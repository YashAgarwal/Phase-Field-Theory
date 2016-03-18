set size ratio 0.5
set yrange[-0.3:0.3]


plot "analysis_10.dat" using 2:1 with linespoints ls 1 title "@t=10" , \
	 "analysis_20.dat" using 2:1 with linespoints ls 2 title "@t=20", \
	 "analysis_50.dat" using 2:1 with linespoints ls 3 title "@t=50"
