#!/bin/tcsh -f
#SBATCH --ntasks=100
#SBATCH --partition=blue
#SBATCH --output=/home/nhoffman/KMAT/run/logfile.log
#SBATCH --quiet
pwd; hostname; date; whoami
cd /home/nhoffman/KMAT/run
mpirun -np 100 /home/nhoffman/KMAT/bin/fitMPI -c kmat_f0_f2_nosym.cfg
echo DONE!; date
