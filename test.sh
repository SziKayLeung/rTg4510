CHAIN=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/CHAINED_ANALYSIS

All_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/All_Merged
TG_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/TG_Merged
WT_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/WT_Merged


SAMPLE_NAMES=(All_Merged TG_Merged WT_Merged)
SAMPLES_MM10_CHAIN_INPUT_DIR=($CHAIN/mm10/All_Merged $CHAIN/mm10/TG_Merged $CHAIN/mm10/WT_Merged) 
SAMPLES_ERCC_CHAIN_INPUT_DIR=($CHAIN/ERCC/All_Merged $CHAIN/ERCC/TG_Merged $CHAIN/ERCC/WT_Merged) 
SAMPLES_MM10_TOFU_DIR=($All_Isoseq3_WKD/TOFU/mm10 $TG_Isoseq3_WKD/TOFU/mm10 $WT_Isoseq3_WKD/TOFU/mm10)
SAMPLES_ERCC_TOFU_DIR=($All_Isoseq3_WKD/TOFU/ERCC $TG_Isoseq3_WKD/TOFU/ERCC $WT_Isoseq3_WKD/TOFU/ERCC)

# make sure not in conda environment when renaming
for i in {0..2}; do	
	sample=${SAMPLE_NAMES[$i]}
	echo "Processing $sample for mm10 Chaining"
	cd ${SAMPLES_MM10_CHAIN_INPUT_DIR[$i]}
	cp ${SAMPLES_MM10_TOFU_DIR[$i]}/{*$sample.collapsed.group.txt*,*$sample.collapsed.filtered.gff*,*$sample.collapsed.filtered.abundance.txt*,*$sample.collapsed.filtered.rep.fq*} .
	echo "Copied over files"
	ls
	rename $sample Sample *	
	
	echo "Processing $sample for ERCC Chaining"
	cd ${SAMPLES_ERCC_CHAIN_INPUT_DIR[$i]}
	cp ${SAMPLES_ERCC_TOFU_DIR[$i]}/{*$sample.collapsed.group.txt*,*$sample.collapsed.filtered.gff*,*$sample.collapsed.filtered.abundance.txt*,*$sample.collapsed.filtered.rep.fq*} .
	echo "Copied over files"
	ls
	rename $sample Sample *	
done

