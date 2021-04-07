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


# 12/12/2020: IsoformSwitchAnalyze error from sqanti gtf
# trouble shoot by: 1) running SQANTI2 with and 2) without chaining, 3) running SQANTI3 without chaining

#************************************* DEFINE GLOBAL VARIABLES
# File directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
TROUBLESHOOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Testing/Linux/Troubleshoot
CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
CAGE_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/CAGE
SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2

cd $TROUBLESHOOT; mkdir SQANTI2_chain SQANTI2_nochain SQANTI3_nochain
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Testing/Linux
TOFU=$Isoseq3_WKD/TOFU
CHAIN=$Isoseq3_WKD/CHAIN

#************************************* SQANTI3 no chain (test on O18)
source $FUNCTIONS/Post_Isoseq3_Function.sh
# run_sqanti3_QC <output_prefix_sample> <input_fasta> <input_abundance_file> <coverage/genome=mm10_rnqaseq/mm10/hg38_gencode/hg38_chess/ERCC> <input_kallisto_file> <input_rnaseq_dir> <output_dir>
cd $TROUBLESHOOT/SQANTI3_nochain
run_sqanti3_QC O18 $TOFU/O18.collapsed.filtered.gff $TOFU/O18.collapsed.abundance.txt mm10 NA NA $TROUBLESHOOT/SQANTI3_nochain
# run_sqanti_Filter <input/output_dir> <classification.txt> <corrected.fasta> <corrected.gtf>
run_sqanti_Filter $TROUBLESHOOT/SQANTI3_nochain O18.collapsed.filtered_classification.txt O18.collapsed.filtered_corrected.fasta O18.collapsed.filtered_corrected.gtf

module load Miniconda2
source activate sqanti2_py3
#************************************* SQANTI2 no chain (test on O18)
cd $TROUBLESHOOT/SQANTI2_nochain
export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
python $SQANTI2_DIR/sqanti_qc2.py -t 30 \
            --gtf $TOFU/O18.collapsed.filtered.gff \
            $REFERENCE/gencode.vM22.annotation.gtf  \
            $REFERENCE/mm10.fa \
            --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed \
            --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt \
            --fl_count $TOFU/O18.collapsed.abundance.txt \
            &> O18.sqanti.qc.log

python $SQANTI2_DIR/sqanti_filter2.py O18.collapsed.filtered_classification.txt O18.collapsed.filtered_corrected.fasta O18.collapsed.filtered_corrected.gtf -a 0.6 -c 3


#************************************* SQANTI2 chain
cd $TROUBLESHOOT/SQANTI2_chain
python $SQANTI2_DIR/sqanti_qc2.py -t 30 \
            --gtf $CHAIN/all_samples.chained.gff \
            $REFERENCE/gencode.vM22.annotation.gtf  \
            $REFERENCE/mm10.fa \
            --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed \
            --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt \
            --fl_count $CHAIN/all_samples.chained_count.txt \
            &> O18.sqanti.qc.log

python $SQANTI2_DIR/sqanti_filter2.py all_samples.chained_classification.txt all_samples.chained_corrected.fasta all_samples.chained_corrected.gtf -a 0.6 -c 3

################ CONCLUSION: SQANTI2 chain works well on IsoformSwitchAnalyze (Allows input of gtf)
cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual_Samples_OLD/CHAIN/SQANTI2
export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
python $SQANTI2_DIR/sqanti_qc2.py -t 30 \
            --gtf ../all_samples.chained.gff \
            $REFERENCE/gencode.vM22.annotation.gtf  \
            $REFERENCE/mm10.fa \
            --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed \
            --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt \
            --fl_count ../all_samples.chained_count.txt \
            --skipORF \
            &> chained.sqanti.qc.log

python $SQANTI2_DIR/sqanti_filter2.py all_samples.chained_classification.txt all_samples.chained_corrected.fasta all_samples.chained_corrected.gtf -a 0.6 -c 3
