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
#SBATCH --output=Individual_Part4.o
#SBATCH --error=Individual_Part4.e

# 14/12/2020: TAMA merge output from collapsed.gtf

#************************************* DEFINE GLOBAL VARIABLES
# Process tama collapse for both cap and no_cap to later run tama_go for degree of degradation
# tama_collapse_nocap $Sample_Name

TAMA=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/tama
TAMA_OUTPUT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual/TAMA
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
MAPPING=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual/MAPPING
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation

tama_collapse_nocap(){

		source activate sqanti2
		cd $3
    # -s = Sorted.sam file, -f = Genome Fasta File, -p = Output prefex, -d merge_duplicates, -x = no cap
    # python tama_collapse.py -s {sam} -f {genome} -p {prefix} -d merge_dup -x no_cap -a 100 -z 100 -sj sj_priority -log log_off (Kuo 13/05/2019)
    python $TAMA/tama_collapse.py -s $1.sorted.sam -f $REFERENCE/mm10.fa -p $1_nocap -d merge_dup -x no_cap -a 100 -z 100 -sj sj_priority &> $1.tama_nocap.log
    echo "Tama collapse of Sample $1 successful"
		source deactivate
}

run_tama_merge(){
    python $TAMA/tama_merge.py -f $2 -p $1_merged
}

SAMPLES_NAMES=(O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23)
for i in ${SAMPLES_NAMES[@]}; do tama_collapse_nocap $MAPPING/$i".clustered.hq.fastq" $TAMA_OUTPUT; done

### J.Humphrey
module load Miniconda2
source activate sqanti2_py3

TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples_OLD/mm10/TOFU
CLUSTER=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples_OLD/CLUSTER


cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples_OLD/mm10/CHAIN_DE
python /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/post_isoseq_cluster/demux_isoseq_with_genome.py --mapped_fafq $TOFU/All_Merged.collapsed.filtered.rep.fa --read_stat $TOFU/All_Merged.collapsed.read_stat.txt --cluster $CLUSTER/All_Merged.clustered.cluster_report.csv --output Individual --primer_names Samples.txt

TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Isoseq_Paper/WT8_ISOSEQ/IsoSeq3.1.2/TOFU
CLUSTER=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Isoseq_Paper/WT8_ISOSEQ/IsoSeq3.1.2/Isoseq3_WKD
CHAIN_DE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples_OLD/mm10/CHAIN_DE
python /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/post_isoseq_cluster/demux_isoseq_with_genome.py --mapped_fafq $TOFU/WT8IsoSeq.collapsed.filtered.rep.fa --read_stat $TOFU/WT8IsoSeq.collapsed.read_stat.txt --classify_csv $CLUSTER/WT8IsoSeq.flnc.report.csv --output WT8  --primer_names WT_Samples.txt

python /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/post_isoseq_cluster/demux_isoseq_with_genome.py --mapped_fafq $TOFU/WT8IsoSeq.collapsed.filtered.rep.fa --read_stat $TOFU/WT8IsoSeq.collapsed.read_stat.txt --classify_csv WT8IsoSeq_mod.flnc.report.csv --output WT8  --primer_names WT_Samples.txt

source $FUNCTIONS/Post_Isoseq3_Function.sh
run_sqanti3_QC WT8IsoSeq $TOFU/WT8IsoSeq.collapsed.filtered.gff WT8 mm10 NA NA $CHAIN_DE
