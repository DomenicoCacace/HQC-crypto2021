#!/bin/zsh

cmake -S .. -B ../build -DSECLVL=128 -DMODE="TIMING-PKE" -DCROSSCOMPILE=1 -DVERBOSE=1 -DMASKLVL=$1
make -C ../build flash
./serialRead.sh > pke &
while ps -p ${!} > /dev/null; do sleep 1; done;

cmake -S .. -B ../build -DSECLVL=128 -DMODE="TIMING-KEM" -DCROSSCOMPILE=1 -DVERBOSE=1 -DMASKLVL=$1
make -C ../build flash
./serialRead.sh > kem &
while ps -p ${!} > /dev/null; do sleep 1; done;

python parseTimingResults.py $(git rev-parse --short HEAD)
rm pke kem
