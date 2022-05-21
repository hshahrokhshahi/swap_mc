#!/bin/bash
#SBATCH -J 'blank'
#SBATCH -t 301:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8Gb
#SBATCH --array=1
#SBATCH --output=haha.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.shahrokhshahi@gmail.com

##module load LAMMPS

time ./haha 0 ./ ./prev_run
