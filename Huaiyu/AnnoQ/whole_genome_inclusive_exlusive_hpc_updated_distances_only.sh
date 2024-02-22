#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --mail-type=BEGIN,END,FAIL

module purge
module load gcc/11.3.0
module load python/3.9.12

python -m pip install pandas numpy --user

python whole_genome_inclusive_exclusive_hpc_updated_distances_only.py
