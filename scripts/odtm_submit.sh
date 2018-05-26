#!/bin/ksh

#BSUB -q "cccr-res"
#BSUB -x
#BSUB -J ODTM
#BSUB -W 24:00
#BSUB -n 8
#BSUB -e stdlog.out
#BSUB -o stdlog.out

export I_MPI_FABRICS=shm:dapl
export OMP_NUM_THREADS=1

ulimit -c unlimited
set -xu

mpirun -np 8 ./odtm.exe

echo "model run exit status : " $?
