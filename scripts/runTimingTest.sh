#!/bin/zsh

cmake -S .. -B ../build -DSECLVL=128 -DMODE="TIMING-PKE" -DCROSSCOMPILE=1 -DVERBOSE=1 -DMASKLVL=$1
make -C ../build flash
./serialRead.sh > pke
cmake -S .. -B ../build -DSECLVL=128 -DMODE="TIMING-KEM" -DCROSSCOMPILE=1 -DVERBOSE=1 -DMASKLVL=$1
make -C ../build flash
./serialRead.sh > kem
python parseTimingResults.py $(git rev-parse --short HEAD)
rm pke kem
