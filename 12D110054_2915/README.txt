The input file is made as follows

***********input.dat****************
Case number(decides the values of A and beta)
n_x(number of nodes)
delta_x
T(number of time steps)
T_write(number of time steps after which the profile in saved)
delta_t
************************************

for the 3 cases use the following input parameters

***********CASE 1****************
1
1024
0.5
100000
1000
0.5

***********CASE 2****************
2
1024
0.5
100000
1000
0.5

***********CASE 3****************
3
1024
0.5
10000
1000
0.001


***********INSTRUCTIONS****************
After setting the input file run ./run.sh
and the output will be stored in the outputX.dat(X=case number)


***********RESULTS****************
for Interface Width
case 1 =X
case 2 =2X
case 3 =X/2

for Interfacial Energy
case 1 =X
case 2 =2X
case 3 =24X
