#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job in hours:minutes:seconds
#PBS -A Research_Project-MRC148213 # research project to submit under.
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# 12/02/2020: chain whole transcriptome data, all samples using cupcake data from Isoseq3.2.1 
# 13/02/2020: run_fusion_finder on all individual samples and run_sqanti2_QC_fusioninput 

##*************************************** Define variables
sl693=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693
CHAIN=$sl693/WholeTranscriptome/Individual/Isoseq3.2.1/CHAIN
CHAIN_FUSION=$sl693/WholeTranscriptome/Individual/Isoseq3.2.1/CHAIN_FUSION
TOFU=$sl693/WholeTranscriptome/Individual/Isoseq3.2.1/TOFU
TOFU_FUSION=$sl693/WholeTranscriptome/Individual/Isoseq3.2.1/TOFU_FUSION
POLISHED=$sl693/WholeTranscriptome/Individual/Isoseq3.2.1/Isoseq3_WKD/CLUSTER
MAPPING=$sl693/WholeTranscriptome/Individual/Isoseq3.2.1/MAPPING
GENERAL_SCRIPT=$sl693/Scripts/general
SQANTI2=$sl693/WholeTranscriptome/Individual/Isoseq3.1.2/SQANTI2
SQANTI2_output_dir=$sl693/WholeTranscriptome/Individual/Isoseq3.2.1/SQANTI2

###*************************************** Process all whole transcriptome samples using tofu_cupcake chain.samples.py 
# Prepare files for chaining as according to Liz's chain.py 
SAMPLE_NAMES=(O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23)

cd $CHAIN
mkdir O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23
for i in ${SAMPLE_NAMES[@]}; do
    echo "Processing $i"
    cd $CHAIN/$i
    cp $TOFU/{*$i.collapsed.group.txt*,*$i.collapsed.filtered.gff*,*$i.collapsed.filtered.abundance.txt*,*$i.collapsed.filtered.rep.fq*} .
    echo "Copied over files"
    ls 
    rename $i Sample *
done

# Prepare configuration file 
cat >> Whole.config_file <<EOL
SAMPLE=O18;$CHAIN/O18
SAMPLE=K18;$CHAIN/K18
SAMPLE=S18;$CHAIN/S18
SAMPLE=L22;$CHAIN/L22
SAMPLE=Q20;$CHAIN/Q20
SAMPLE=K24;$CHAIN/K24
SAMPLE=Q21;$CHAIN/Q21
SAMPLE=K17;$CHAIN/K17
SAMPLE=M21;$CHAIN/M21
SAMPLE=O23;$CHAIN/O23
SAMPLE=S23;$CHAIN/S23
SAMPLE=K23;$CHAIN/K23

GROUP_FILENAME=Sample.collapsed.group.txt
GFF_FILENAME=Sample.collapsed.filtered.gff
COUNT_FILENAME=Sample.collapsed.filtered.abundance.txt
FASTQ_FILENAME=Sample.collapsed.filtered.rep.fq
EOL

module load Miniconda2
source activate sqanti2_py3
echo "Processing samples for chain_samples.py"
chain_samples.py $CHAIN/Whole.config_file count_fl 2> Whole.chain_samples.log
source deactivate 

###*************************************** Generate fusion.gtf using fusion script from Liz
source $GENERAL_SCRIPT/Post_IsoSeq/Post_Isoseq3_Functions.sh
# output in $TOFU_FUSION
run_fusion_finder O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23

###*************************************** SQANTI2 with fusion.gtf as input 
run_sqanti2_QC_fusioninput O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23