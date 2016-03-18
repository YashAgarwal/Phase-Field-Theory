#!/bin/bash
gcc main.c -lgsl -lgslcblas -lfftw3 -lm
./a.out
gnuplot --persist < plotAnimation.gp
