#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Targeted_Part2.o
#SBATCH --error=Targeted_Part2.e

# 14/12/2020: IsoSeq Analysis pipeline (Part2) for individual samples in TargetedSeq2): IsoSeq3, Mapping, TOFU

#************************************* DEFINE VARIABLES
# File directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/IsoSeq
Post_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq
#cd $Isoseq3_WKD; mkdir CCS LIMA REFINE CLUSTER
#cd $Post_Isoseq3_WKD; mkdir MAPPING TOFU CHAIN SQANTI
Targeted_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/

# For Pooled Targeted
### Important order of BAM files is the same as the sample names
SAMPLES_NAMES=(TargetedSeq1 Targeted_Seq_2 Targeted_Seq_3a Targeted_Seq_3b)
cat $Targeted_dir/Isoseq_Targeted_MouseRaw.txt
BAM_FILES=(`cat $Targeted_dir/Isoseq_Targeted_MouseRaw.txt | egrep -v "^\s*(#|$)"`)

# For Demultiplexing Samples
Pooled_Samples_Targeted1=(K19 K23 K21 K18 K20 K17)
Pooled_Samples_Targeted2=(S19 K24 L22 M21 O18 O23 O22 P19 T20)
Pooled_Samples_Targeted3a=(Q20A Q21A S18A S23A Q18A Q17A L18A Q23A T18A) # failed 3rd run
Pooled_Samples_Targeted3b=(Q20 Q21 S18 S23 Q18 Q17 L18 Q23 T18) # successful 3rd run
Barcoded_Targeted1_config_file=$Targeted_dir/Barcode_Configs/Isoseq_Mouse_Targeted1_barcode.config
Barcoded_Targeted2_config_file=$Targeted_dir/Barcode_Configs/Isoseq_Mouse_Targeted2_barcode.config
Barcoded_Targeted3_config_file=$Targeted_dir/Barcode_Configs/Isoseq_Mouse_Targeted3_barcode.config
ALL_SAMPLES_NAMES=(K19 K23 K21 K18 K20 K17 S19 K24 L22 M21 O18 O23 O22 P19 T20 Q20 Q21 S18 S23 Q18 Q17 L18 Q23 T18)

#************************************* Demultiplexing Targeted Samples: Run IsoSeq Analysis pipeline to Cupcake
# Isoseq3 workflow for individual samples in pooled dataset
#run_pipeline_to_tofu <Barcoded_config_file> <input_lima_sample>
source $FUNCTIONS/Isoseq3.2.2_Functions.sh
source $FUNCTIONS/Post_Isoseq3_Function.sh

tofu_re(){
    cd $4; mkdir $1
    cd $1
    source activate cupcake
    echo "Processing Sample $1 for TOFU"
    # Collapse
    collapse_isoforms_by_sam.py --input $2/$1.clustered.hq.fastq --fq -s $3/$1.clustered.hq.fastq.sorted.sam --dun-merge-5-shorter -o Sample &> $1.collapse.log

    # Create Abundance Script of full-length transcripts
    get_abundance_post_collapse.py Sample.collapsed $2/$1.clustered.cluster_report.csv 2> $1.abundance.log

    # Remove degraded isoforms (default setting)
    filter_away_subset.py Sample.collapsed 2> $1.filter.log

    source deactivate

    source activate sqanti2_py3
    # convert rep.fq to rep.fa for SQANTI2 input
    seqtk seq -a Sample.collapsed.filtered.rep.fq > Sample.collapsed.filtered.rep.fa
    echo "Processing Sample $1 for TOFU successful"
}

run_pipeline_to_tofu(){
  SAMPLES=$(echo "${@:3}")

  echo "Processing from $2: $SAMPLES"
  echo "Using: $1"

  for i in ${SAMPLES[@]}; do
    # run_targeted_REFINE $Input_Pooled_Sample $Input_config_file $Input_Lima_sample $Input_LIMA_directory $Output_directory
    #run_targeted_REFINE $i $1 $2 $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE
    #run_CLUSTER $i $Isoseq3_WKD/REFINE $Isoseq3_WKD/CLUSTER
    #convert_fa2fq $i".clustered.hq.fasta" $Isoseq3_WKD/CLUSTER
    #run_minimap2 $i $Isoseq3_WKD/CLUSTER mm10 $Post_Isoseq3_WKD/MAPPING
    #on_target_rate $Targeted_dir/Probes/FINAL_MOUSE.bed $Isoseq3_WKD/CLUSTER/$i".clustered.hq.fasta" $Post_Isoseq3_WKD/MAPPING/$i".clustered.hq.fastq.sam" $Post_Isoseq3_WKD/MAPPING/$i".fasta.sam.probe_hit.txt"
    tofu_re $i $Isoseq3_WKD/CLUSTER $Post_Isoseq3_WKD/MAPPING $Post_Isoseq3_WKD/TOFU
  done
}

run_pipeline_to_tofu $Barcoded_Targeted1_config_file TargetedSeq1 ${Pooled_Samples_Targeted1[@]}
run_pipeline_to_tofu $Barcoded_Targeted2_config_file Targeted_Seq_2 ${Pooled_Samples_Targeted2[@]}
run_pipeline_to_tofu $Barcoded_Targeted3_config_file Targeted_Seq_3a ${Pooled_Samples_Targeted3a[@]}
run_pipeline_to_tofu $Barcoded_Targeted3_config_file Targeted_Seq_3b ${Pooled_Samples_Targeted3b[@]}


#************************************* Chaining across all samples (use successful 3rd run)
cd $Post_Isoseq3_WKD/CHAIN
mkdir ${ALL_SAMPLES_NAMES[@]}

cat << EOF > $Post_Isoseq3_WKD/CHAIN/Chained_Configuration.txt
SAMPLE=K19;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/K19
SAMPLE=K23;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/K23
SAMPLE=K21;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/K21
SAMPLE=K18;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/K18
SAMPLE=K20;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/K20
SAMPLE=K17;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/K17
SAMPLE=S19;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/S19
SAMPLE=K24;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/K24
SAMPLE=L22;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/L22
SAMPLE=M21;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/M21
SAMPLE=O18;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/O18
SAMPLE=O23;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/O23
SAMPLE=O22;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/O22
SAMPLE=P19;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/P19
SAMPLE=T20;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/T20
SAMPLE=Q20;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/Q20
SAMPLE=Q21;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/Q21
SAMPLE=S18;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/S18
SAMPLE=S23;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/S23
SAMPLE=Q18;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/Q18
SAMPLE=Q17;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/Q17
SAMPLE=L18;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/L18
SAMPLE=Q23;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/Q23
SAMPLE=T18;/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq/CHAIN/T18

GROUP_FILENAME=Sample.collapsed.group.txt
GFF_FILENAME=Sample.collapsed.filtered.gff
COUNT_FILENAME=Sample.collapsed.filtered.abundance.txt
FASTQ_FILENAME=Sample.collapsed.filtered.rep.fq
EOF


echo "Processing samples for chain_samples.py"

module load Miniconda2
source activate cupcake
chain_samples.py $Post_Isoseq3_WKD/CHAIN/Chained_Configuration.txt count_fl --dun-merge-5-shorter 2> Chained_Configuration.log
