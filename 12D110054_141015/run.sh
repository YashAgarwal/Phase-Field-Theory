#! /bin/bash
gcc main.c -lm -lgsl -lgslcblas
./a.out
rm a.out
gnuplot --persist < plotAnimation.gp
