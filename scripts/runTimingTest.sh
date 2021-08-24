#!/bin/zsh

cmake -S .. -B ../build -DSECLVL=128 -DMODE="TIMING" -DCROSSCOMPILE=1 -DVERBOSE=1
make -C ../build flash
./serialRead.sh > temp
python parseResults.py $(git rev-parse --short HEAD)
rm temp

#cmake -S .. -B ../build -DSECLVL=192 -DMODE="TIMING" -DCROSSCOMPILE=1 -DVERBOSE=1
#make -C ../build flash
#./serialRead.sh

#cmake -S .. -B ../build -DSECLVL=256 -DMODE="TIMING" -DCROSSCOMPILE=1 -DVERBOSE=1
#make -C ../build flash
#./serialRead.sh
