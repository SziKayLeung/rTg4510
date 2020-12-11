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
#SBATCH --output=Isoseq3_alls-%A_%a.o
#SBATCH --error=Isoseq3_all-%A_%a.e

################## chain individual samples
# Prepare files for chaining as according to Liz's chain.py

SAMPLE_NAMES=(O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23)
TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/TOFU
cd $CHAIN_INDIVIDUAL; mkdir ${SAMPLE_NAMES[@]}

# make sure not in conda environment when renaming
for sample in ${SAMPLE_NAMES[@]}; do
	echo "Processing $sample for mm10 Chaining"
	cd $CHAIN_INDIVIDUAL/$sample
	cp $TOFU/{*$sample.collapsed.group.txt*,*$sample.collapsed.filtered.gff*,*$sample.collapsed.filtered.abundance.txt*,*$sample.collapsed.filtered.rep.fq*} .
	echo "Copied over files"
	ls
	rename $sample Sample *
done


# Prepare configuration file
cat << EOF > $CHAIN_INDIVIDUAL/Chained_Configuration.txt
SAMPLE=O18;$CHAIN_INDIVIDUAL/O18
SAMPLE=K18;$CHAIN_INDIVIDUAL/K18
SAMPLE=S18;$CHAIN_INDIVIDUAL/S18
SAMPLE=L22;$CHAIN_INDIVIDUAL/L22
SAMPLE=Q20;$CHAIN_INDIVIDUAL/Q20
SAMPLE=K24;$CHAIN_INDIVIDUAL/K24
SAMPLE=Q21;$CHAIN_INDIVIDUAL/Q21
SAMPLE=K17;$CHAIN_INDIVIDUAL/K17
SAMPLE=M21;$CHAIN_INDIVIDUAL/M21
SAMPLE=O23;$CHAIN_INDIVIDUAL/O23
SAMPLE=S23;$CHAIN_INDIVIDUAL/S23
SAMPLE=K23;$CHAIN_INDIVIDUAL/K23

GROUP_FILENAME=Sample.collapsed.group.txt
GFF_FILENAME=Sample.collapsed.filtered.gff
COUNT_FILENAME=Sample.collapsed.filtered.abundance.txt
FASTQ_FILENAME=Sample.collapsed.filtered.rep.fq
EOF

echo "Processing samples for chain_samples.py"
CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/cupcake/tofu/counting
cd $CHAIN_INDIVIDUAL;python $CUPCAKE/chain_samples_py3.py $CHAIN_INDIVIDUAL/Chained_Configuration.txt count_fl --dun-merge-5-shorter 2> Chained_Configuration.log

SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence
export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
python $SQANTI2_DIR/sqanti_qc2.py --version
python $SQANTI2_DIR/sqanti_qc2.py -t 30 --gtf ../all_samples.chained.gff $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --expression ../all_samples.chained_count.txt

SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence

