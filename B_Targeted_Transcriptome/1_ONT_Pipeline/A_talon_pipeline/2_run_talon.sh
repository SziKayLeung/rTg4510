#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=TALON.o
#SBATCH --error=TALON.e

# 11/10/2021: Run Minimap2 and TALON on ONT Targeted rTg4510 dataset

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
source $SC_ROOT/1_ONT_Pipeline/rTg4510_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh

##-------------------------------------------------------------------------

# 8) create_talon_db <talon_database_name>
create_talon_db $TALON_NAME

# 8c) Run TALON with a config file
# run_talon_and_quantify <name> <root_dir>
run_talon_and_quantify $NAME $WKD_ROOT

# talon_filter_quantify <name> <root_dir>
talon_filter_quantify $NAME $WKD_ROOT

# run_sqanti3 <sample> <mode=basic/full/nokallisto/lncrna> <input_dir> <dataset>
run_sqanti3 $NAME"_filtered_talon" basic $WKD_ROOT/6_talon/2_talon_full Unfiltered
run_sqanti3 $NAME"_filtered_talon" basic $WKD_ROOT/6_talon/3_talon_filter Filtered

