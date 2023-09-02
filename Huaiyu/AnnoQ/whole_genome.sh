#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=184GB
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL

module purge
module load gcc/8.3.0
module load python/3.7.6

python whole_genome.py
