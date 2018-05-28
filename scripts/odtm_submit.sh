#!/bin/ksh

#BSUB -x

#export I_MPI_FABRICS=shm:dapl
export OMP_NUM_THREADS=1

ulimit -c unlimited

mpirun -np $NPROCS ./odtm.exe
