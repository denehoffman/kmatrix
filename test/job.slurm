#!/bin/bash
#SBATCH --job-name=slurm-cmd
#SBATCH --output=/raid3/nhoffman/kmatrix/test/slurm.log
#SBATCH --error=/raid3/nhoffman/kmatrix/test/slurm.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=6364324128@txt.att.net
#SBATCH --ntasks=100
#SBATCH --partition=blue

cd /raid3/nhoffman/kmatrix/test 
./run.sh test
