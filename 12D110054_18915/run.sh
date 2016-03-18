#! /bin/bash
gcc main.c -lm -lgsl -lgslcblas
./a.out
gnuplot --persist < plotMain.gp
rm a.out
