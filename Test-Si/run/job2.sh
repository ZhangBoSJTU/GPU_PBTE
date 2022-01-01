#!/bin/bash

cd ..
rm 1 -rf
set -e 
echo 2 > ./input/jobtype.txt
export LD_LIBRARY_PATH="../../external/spglib-1.7.4-intel16/lib:$LD_LIBRARY_PATH"
EXE=../../bin/release/GPUPBTE.x
SLURM_NPROCS=1
myproc=1
mkdir $myproc
echo $SLURM_NPROCS > $myproc/jobid.txt
echo $myproc >> $myproc/jobid.txt
cd $myproc
$EXE > out

