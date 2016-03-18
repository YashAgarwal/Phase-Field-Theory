This Code implements Allen Cahn equation for a particular input concentration

/*********How to run this code********/
First set all the parameter values in the file "input.dat"pe
the file input.dat has the following order

/**input.dat
n_x(number of nodes in the system)
delta_x(distance between consecutive nodes)
T(number of time steps)
delta_t(discrete time increment)
T_write(After how many time steps we print the profile in a file)
**/

then run the file run.sh by
./run.sh

Result
the final profile takes the shape of tanh stable interfaces as this is a metastable state for the system