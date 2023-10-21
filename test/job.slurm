#!/bin/bash
#SBATCH --job-name=slurm-cmd
#SBATCH --output=slurm.log
#SBATCH --error=slurm.log
#SBATCH --ntasks=100

cd /raid3/nhoffman/kmatrix/test 
./run.sh test
