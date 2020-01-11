#!/bin/bash
#SBATCH -n 3                        # number of cores
#SBATCH -t 0-12:00                  # wall time (D-HH:MM)
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
module purge
module load intel-mpi/.2017.4
mpiifort -o test1 ParGaussV2.for 
export I_MPI_PMI_LIBRARY=/path/to/slurm/pmi/library/libpmi.so
mpirun -n 3 ./test1 
rm test1