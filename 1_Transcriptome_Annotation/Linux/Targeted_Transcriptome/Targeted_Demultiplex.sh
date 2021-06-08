#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Targeted_Demultiplex.o
#SBATCH --error=Targeted_Demultiplex.e

# 15/12/2020: Merge all Samples and then demultiplex later for count matrix

############# Alternatively cat multiple refine files and therefore no need to cluster (which takes a lot of time)
#************************************* DEFINE GLOBAL VARIABLES
# File directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
DEMUX_SCRIPT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/1_Transcriptome_Annotation/Linux/Targeted_Transcriptome/
CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/cupcake/tofu/counting
REFINE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/IsoSeq/REFINE
All_Samples=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged
cd $All_Samples; mkdir CLUSTER MAPPING TOFU SQANTI

#************************************* All samples, WT, TG merged at refine
source $FUNCTIONS/Isoseq3.2.2_Functions.sh
source $FUNCTIONS/Post_Isoseq3_Function.sh
module load Miniconda2
module load R

#run_pipeline_to_chain <Sample_Name> <List.of.Samples...>
run_pipeline_to_tofu(){

  SAMPLES=$(echo "${@:2}")
  echo "Processing together: $SAMPLES"

  # merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
  # convert_fa2fq <file_name> <input_dir>
  # run_minimap2 <prefix_sample> <input_dir> <ERCC/mm10> <output_dir>
  # tofu <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir>
  # chain_prepare <ERCC/MOUSE>
  #merging_at_refine $REFINE $All_Samples/CLUSTER $1 ${SAMPLES[@]}
  #convert_fa2fq $1".clustered.hq.fasta" $All_Samples/CLUSTER
  #run_minimap2 $1 $All_Samples/CLUSTER mm10 $All_Samples/MAPPING
  #tofu $1 $All_Samples/CLUSTER $All_Samples/MAPPING $All_Samples/TOFU

  # demultiplex collapsed.read_stat.txt
  Rscript $DEMUX_SCRIPT/Demultiplex_Cupcake.R
}

run_sqanti2(){

    source activate sqanti2_py3
    SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
    CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence
    REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
    CAGE_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/CAGE

    echo "Processing Sample $1 for SQANTI2 QC"
    cd $2

    export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
    python $SQANTI2_DIR/sqanti_qc2.py -v

    if python $SQANTI2_DIR/sqanti_qc2.py -t 30 \
        --gtf $3 \
        $REFERENCE/gencode.vM22.annotation.gtf  \
        $REFERENCE/mm10.fa \
        --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed \
        --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt \
        --fl_count $4 \
        2> $1.sqanti.qc.log; then
        echo "Processing Sample $1 for SQANTI2 successful"
    else
        echo "Processing Sample $1 for SQANTI2 failed"
    fi

    mkdir Unfiltered_PNG
    mv *png* Unfiltered_PNG/

    if python $SQANTI2_DIR/sqanti_filter2.py \
        $1.collapsed.filtered_classification.txt \
        $1.collapsed.filtered_corrected.fasta  \
        $1.collapsed.filtered_corrected.gtf \
        -a 0.6 -c 3 2> $1.sqanti.filter.log; then
        echo "Processing Sample $1 for SQANTI2 filter successful"
    else
        echo "Processing Sample $1 for SQANTI2 filter failed"
    fi

    mkdir filtered_PNG
    mv *png* filtered_PNG/

    source deactivate

}

run_pipeline_to_tofu All_Targeted_Merged K17 K18 K19 K20 K21 K23 K24 L18 L22 M21 O18 O22 O23 P19 Q17 Q18 Q20 Q21 Q23 S18 S19 S23 T18 T20
# run_sqanti2_QC <sample_name> <working_directory> <collapsed.gtf> <count.txt>
run_sqanti2 All_Targeted_Merged $All_Samples/SQANTI $All_Samples/TOFU/All_Targeted_Merged.collapsed.filtered.gff $All_Samples/TOFU/All_Targeted_Merged.Demultipled_Abundance.txt
