#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=1:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

module load Python/3.6.6-foss-2018b
module load Anaconda2
source activate my_root

# Python Module OS installation: conda install -c jmcmurray os
# Python Module OS installation: conda install -c conda-forge glob2
# Python Module OS installation: conda install -c menpo pathlib

cd .
python LIMA.py