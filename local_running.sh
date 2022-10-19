#!/bin/bash 

make
export OMP_NUM_THREADS=16
./eb.out


