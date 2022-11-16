#!/bin/tcsh -f
#SBATCH --ntasks=200
#SBATCH --partition=red
#SBATCH --output=/home/nhoffman/KMAT/run/logfile.log
#SBATCH --quiet
pwd; hostname; date; whoami
cd /home/nhoffman/KMAT/run
mpirun --mca btl_openib_allow_ib 1 /home/nhoffman/KMAT/bin/fitMPI -c kmat_all.cfg
echo DONE!; date
