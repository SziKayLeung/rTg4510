#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=0:10:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=../../../bash_output/1a_prepare_samples.o
#SBATCH --error=../../../bash_output/1a_prepare_samples.e

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
FICLE_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/FICLE/
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/rTg4510_merged.config
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/01_source_function.sh
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous


##-------------------------------------------------------------------------

source activate sqanti2_py3

export dir=${CUPMERGE_DIR}
export samplename=all_iso_ont

cd ${dir}
mkdir -p 1_align 2_collapse 3_sqanti3


##-------------------------------------------------------------------------

# copy and replace filename in ONT transcript clean directory
echo "Replacing filenames in ONT transcript clean directory"
replace_filenames_with_csv.py --copy -i=${ONT_TCLEAN_DIR}/Batch2 -f=$META_ROOT/ONT_Batch2_Tclean_rename.csv -d=${ONT_TCLEAN_DIR}/AllBatch2Batch3
replace_filenames_with_csv.py --copy -i=${ONT_TCLEAN_DIR}/Batch3 -f=$META_ROOT/ONT_Batch3_Tclean_rename.csv -d=${ONT_TCLEAN_DIR}/AllBatch2Batch3
cd ${ONT_TCLEAN_DIR}/AllBatch2Batch3; rm *clean.sam* *clean.log*