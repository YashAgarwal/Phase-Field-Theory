set size ratio 0.75

set key outside box 
set  border

set title "fD vs t"
set xlabel "time"
set ylabel "fD"

plot "fD_1.dat"  with line  lw 3   title "band 1" , \
	 "fD_2.dat"  with line  lw 3   title "band 2", \
	 "fD_3.dat"  with line  lw 3   title "band 3", \
	 "fD_4.dat"  with line  lw 3   title "band 4", \
	 "fD_5.dat"  with line  lw 3   title "band 5", \
	 "fD_6.dat"  with line  lw 3   title "band 6", \
	 "SS_1.dat"  with line  lw 3   title "SS-1%", \
	 "SS_4.dat"  with line  lw 3   title "SS-4%"

