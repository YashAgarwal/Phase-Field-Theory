#!/bin/bash

#set your current working directory to this folder and then run this script

echo "This script will create the entire needed for the analysis"

echo "512" > ./input.dat
echo "1.0" >> ./input.dat
echo "256" >> ./input.dat
echo "1.0" >> ./input.dat
echo "0.42" >> ./input.dat
echo "200" >> ./input.dat
echo "0.04" >> ./input.dat
echo "1.1" >> ./input.dat
./run.sh

echo "20% complete"

echo "512" > ./input.dat
echo "1.0" >> ./input.dat
echo "256" >> ./input.dat
echo "1.0" >> ./input.dat
echo "0.42" >> ./input.dat
echo "200" >> ./input.dat
echo "0.01" >> ./input.dat
echo "1.1" >> ./input.dat
./run.sh

echo "40% complete"

echo "512" > ./input.dat
echo "1.0" >> ./input.dat
echo "256" >> ./input.dat
echo "1.0" >> ./input.dat
echo "0.42" >> ./input.dat
echo "100" >> ./input.dat
echo "0.04" >> ./input.dat
echo "1.2" >> ./input.dat
./run.sh

echo "60% complete"


echo "512" > ./input.dat
echo "1.0" >> ./input.dat
echo "256" >> ./input.dat
echo "1.0" >> ./input.dat
echo "0.42" >> ./input.dat
echo "100" >> ./input.dat
echo "0.04" >> ./input.dat
echo "1.3" >> ./input.dat
./run.sh

echo "80% complete"


echo "512" > ./input.dat
echo "1.0" >> ./input.dat
echo "256" >> ./input.dat
echo "1.0" >> ./input.dat
echo "0.42" >> ./input.dat
echo "100" >> ./input.dat
echo "0.01" >> ./input.dat
echo "2.1" >> ./input.dat
./run.sh

echo "100% complete"
