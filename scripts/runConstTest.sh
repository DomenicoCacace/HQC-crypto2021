#!/bin/zsh

cmake -S .. -B ../build -DSECLVL=128 -DMODE="CONST" -DCROSSCOMPILE=1 -DVERBOSE=1
make -C ../build flash
./serialRead.sh > temp
python parseConstResults.py $(git rev-parse --short HEAD)
rm temp

#cmake -S .. -B ../build -DSECLVL=192 -DMODE="CONST" -DCROSSCOMPILE=1 -DVERBOSE=1
#make -C ../build flash
#./serialRead.sh

#cmake -S .. -B ../build -DSECLVL=256 -DMODE="CONST" -DCROSSCOMPILE=1 -DVERBOSE=1
#make -C ../build flash
#./serialRead.sh
