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

# #!/bin/bash

# #SBATCH --job-name=Test-Si
# #SBATCH --partition=dgx2
# #SBATCH -N 1
# #SBATCH -n 1
# #SBATCH --ntasks-per-node 1
# #SBATCH --cpus-per-task=1
# #SBATCH --gres=gpu:1
# #SBATCH --output=%j.out
# #SBATCH --error=%j.err 

# cd ..
# rm -rf 1 2 3 4 5 6 7 8 9 10
# set -e 
# echo 2 > ./input/jobtype.txt
# module load intel-mkl
# module load intel-parallel-studio
# module load cuda/10.2.89-gcc-8.3.0
# module load nvhpc
# export LD_LIBRARY_PATH="../../external/spglib-1.7.4-intel16/lib:$LD_LIBRARY_PATH"
# EXE=../../bin/release/GPUPBTE.x
# SLURM_NPROCS=1
# myproc=1
# mkdir $myproc
# echo $SLURM_NPROCS > $myproc/jobid.txt
# echo $myproc >> $myproc/jobid.txt
# cd $myproc
# $EXE > out
# wait
# cd ..