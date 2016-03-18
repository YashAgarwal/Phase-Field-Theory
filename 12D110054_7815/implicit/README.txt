First Make the input file

////input.dat/////
N
delx
T
delt
alpha
Boundary-Condition-Left(1,2 or 3) Boundary-Condition-Right(1,2 or 3) 
//////////////////

then run ./make to make c_0.dat--initial profile

then run ./a.out-- the main program 

then run the octave script plotAnimation.oct to see the concentration profile change with time 



