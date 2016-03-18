First set the values in the input file

/***input.dat****/
n_x
n_y
A
delta_x
delta_y
delta_y
T
T_write:after how many time steps do we print 
/****************/

first compile using  gcc test.c -lgsl -lgslcblas -lfftw3 -lm

then run ./a.out 

then run the gnuPlot script plotAnimation to see the concentration profile change with time 



