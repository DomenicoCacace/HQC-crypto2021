#!/bin/zsh

cmake -S .. -B ../build -DSECLVL=128 -DMODE="CONST-KEM" -DCROSSCOMPILE=1 -DVERBOSE=1
make -C ../build flash
./serialRead.sh > kem
cmake -S .. -B ../build -DSECLVL=128 -DMODE="CONST-PKE" -DCROSSCOMPILE=1 -DVERBOSE=1
make -C ../build flash
./serialRead.sh > pke
python parseConstResults.py $(git rev-parse --short HEAD)
rm kem pke

#cmake -S .. -B ../build -DSECLVL=192 -DMODE="CONST" -DCROSSCOMPILE=1 -DVERBOSE=1
#make -C ../build flash
#./serialRead.sh

#cmake -S .. -B ../build -DSECLVL=256 -DMODE="CONST" -DCROSSCOMPILE=1 -DVERBOSE=1
#make -C ../build flash
#./serialRead.sh
