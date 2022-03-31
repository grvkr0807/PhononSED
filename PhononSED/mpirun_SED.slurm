#!/bin/bash -l

# #SBATCH --job-name=LAMMPS_1ns
# #SBATCH --comment=LAMMPS_1ns

#SBATCH -t 100-0
#SBATCH -n 24
#SBATCH -w, --nodelist=node02
# ,node02,node03,node04,node06,node07,node08,node09,node10,node11
#SBATCH -e err_%j.out
module purge
module load mpi/openmpi-x86_64

export OMP_NUM_THREADS=1
cd $SLURM_SUBMIT_DIR
mpirun -n $SLURM_NTASKS -x OMP_NUM_THREADS -wdir $SLURM_SUBMIT_DIR ./PhononSED.x > SED_$SLURM_JOBID.out