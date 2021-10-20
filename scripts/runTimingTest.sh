#!/bin/zsh

cmake -S .. -B ../build -DSECLVL=128 -DMODE="TIMING-PKE" -DCROSSCOMPILE=1 -DVERBOSE=1
make -C ../build flash
./serialRead.sh > pke
cmake -S .. -B ../build -DSECLVL=128 -DMODE="TIMING-KEM" -DCROSSCOMPILE=1 -DVERBOSE=1
make -C ../build flash
./serialRead.sh > kem
python parseTimingResults.py $(git rev-parse --short HEAD)
rm pke kem

#cmake -S .. -B ../build -DSECLVL=192 -DMODE="TIMING" -DCROSSCOMPILE=1 -DVERBOSE=1
#make -C ../build flash
#./serialRead.sh

#cmake -S .. -B ../build -DSECLVL=256 -DMODE="TIMING" -DCROSSCOMPILE=1 -DVERBOSE=1
#make -C ../build flash
#./serialRead.sh
