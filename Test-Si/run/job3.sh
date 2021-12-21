#!/bin/bash

set -e 
echo 3 > ../input/jobtype.txt
export LD_LIBRARY_PATH="../../external/spglib-1.7.4-intel16/lib:$LD_LIBRARY_PATH"
EXE=../../bin/release/GPUPBTE.x
$EXE > job3.out