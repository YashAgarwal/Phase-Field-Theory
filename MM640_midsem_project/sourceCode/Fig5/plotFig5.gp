set size ratio 0.75

set key right bottom  box 
set  border


set title "R_1 vs t"
set xlabel "time"
set ylabel "R_1"

plot "R(t)_Ia.dat"  with line  title "Ia" , \
	 "R(t)_Ib.dat"  with line  title "Ib", \
	 "R(t)_Ic.dat"  with line  title "Ic"

