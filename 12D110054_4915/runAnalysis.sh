#!/bin/bash
gcc analysis.c -lgsl -lgslcblas -lfftw3 -lm
./a.out
gnuplot --persist < plotAnalysis.gp
