set size ratio 2
set xrange[0:256]
set yrange[0:512]
unset colorbox
unset key 
unset tics 
unset border

#set palette defined (0.49 'black',0.515 'white')
#set cbrange[0.495:0.51]

set palette defined (0 'black',1 'white')
set cbrange[0:1]


set title "Ia, delc=1%"

set multiplot layout 1, 3
set tmargin 2

set title "t=10"
unset key
#set palette defined (0.345 'black',0.605 'white')
#set cbrange[0.345:0.605]

plot "../mainData/output_Ia_0.04/c_10.dat" matrix with image

set title "t=50"
unset key
#set palette defined (0.15 'black',0.85 'white')
#set cbrange[0.15:0.85]

plot "../mainData/output_Ia_0.04/c_50.dat" matrix with image

set title "t=200"
unset key
#set palette defined (0 'black',1 'white')
#set cbrange[0:1]

plot "../mainData/output_Ia_0.04/c_200.dat" matrix with image

unset multiplot
