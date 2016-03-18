set cbrange[0:1]
set xrange[0:256]
set yrange[0:512]
set size ratio 2
set palette gray
unset colorbox
unset key; unset tics; unset border
set title "Final Profile"
set palette defined (0 'black',1 'white')
set cbrange[0:1]
plot "finalProfile_IIa_0.01.dat" matrix with image
