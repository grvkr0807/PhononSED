#!/bin/bash

#SBATCH -t 28-0
#SBATCH -n 10
#SBATCH -w, --nodelist=node06
#SBATCH -e err
# # SBATCH --mail-user=gkumar1@umd.edu
# # SBATCH --mail-type=ALL
# # export OMP_NUM_THREADS=1

module purge
# export PATH=$PATH:/home/gkumar1/gulp-4.5/gulp-4.5/Src/
module load gulp/4.2.2/openmpi/1.10.3-3.el7

material='RDX';
super=1;
nk=1;
ktype='gamma';

# GULP simulation involving multiple unitcells and/or multiple kpoints(not shrink), specific naming for input files
gulp <input/${material}_N${super}_k${nk}.gin> output/${material}_N${super}_k${nk}.out
awk -f io_freq/awk_programs/get_frequencies_${ktype}.awk output/${material}_N${super}_k${nk}.out > output/freqs_${material}_N${super}_k${nk}.dat