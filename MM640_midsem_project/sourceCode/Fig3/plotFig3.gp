set size ratio 0.75
set xrange[128:256]
set yrange[416:512]
unset colorbox
unset key; unset tics; unset border

set palette defined (0 'black',1 'white')

#set palette defined (0.49 'black',0.515 'white')
#set cbrange[0.495:0.51]

set multiplot layout 1, 3
set tmargin 2

set title "Ia, delc=4%"
unset key

plot "../mainData/output_Ia_0.04/c_100.dat" matrix with image

set title "Ia, delc=1%"
unset key
plot "../mainData/output_Ia_0.01/c_100.dat" matrix with image

set title "IIa, delc=1%"
unset key
plot "../mainData/output_IIa_0.01/c_100.dat" matrix with image

unset multiplot
