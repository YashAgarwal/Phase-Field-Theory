#!/bin/bash

#set your current working directory to this folder and then run this script

./clean.sh

cd ./data/withGb

./run.sh

cd ../withoutGb

echo "512" > ./input.dat
echo "1.0" >> ./input.dat
echo "256" >> ./input.dat
echo "1.0" >> ./input.dat
echo "0.42" >> ./input.dat
echo "200" >> ./input.dat
echo "0.04" >> ./input.dat
echo "1.1" >> ./input.dat
./run.sh

echo "512" > ./input.dat
echo "1.0" >> ./input.dat
echo "256" >> ./input.dat
echo "1.0" >> ./input.dat
echo "0.42" >> ./input.dat
echo "200" >> ./input.dat
echo "0.01" >> ./input.dat
echo "1.1" >> ./input.dat
./run.sh

