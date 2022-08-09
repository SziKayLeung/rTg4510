conda create -n nanometh python=3.3
pip install megalodon
megalodon --help-long

wget https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.2.1_linux64.tar.gz
tar -xf ont-guppy-cpu_6.2.1_linux64.tar.gz

#https://community.nanoporetech.com/posts/guppy-v6-1-1-release
# https://labs.epi2me.io/notebooks/Modified_Base_Tutorial.html
# https://github.com/epi2me-labs/modbam2bed
export PATH=$PATH:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/ont-guppy-cpu/bin
ONT_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/ont/runs/minion/10074/20200807_1632_MC-110214_0_add313_506ffc5b/fast5_pass
SAMPLE=$ONT_DIR/barcode09
WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Mouse_Whole_Genome
GENOME_FASTA=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/mm10.fa
guppy_basecaller -i $SAMPLE -s $WKD -c dna_r9.4.1_450bps_modbases_5mc_cg_hac.cfg --bam_out --recursive --align_ref $GENOME_FASTA

source activate nanopore
cd $WKD/pass/; for i in *bam*; do samtools index $i; done
samtools index $WKD
source activate modbam2bed
cd $WKD
modbam2bed -e -m 5mC --cpg -t 4 $GENOME_FASTA $WKD/pass/*.bam > guppy.cpg.bam

cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Mouse_Whole_Genome/3_methylation
samtools view -h /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Mouse_Whole_Genome/2_basecalled/barcode12_old/pass/bam_runid_082328671fc89cf6e677dcf29cdcb679effd5880_0_0.bam > out.sam

pip install sorted-nearest
methplotlib -m sample.bedgraph -n out -w chr11:101,899,127-101,899,128

awk 'BEGIN{OFS="\t"} {$4=$11; NF=4; print}' old_guppy.cpg.bam > sample.bedgraph

# nanopolish
cd $WKD/3_methylation/nanopolish
sample=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Mouse_Whole_Genome/2_basecalled/barcode12_old/pass
export PATH=$PATH:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/nanopolish
nanopolish index -d /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Mouse_Whole_Genome/1_raw/barcode12/ $sample/fastq_runid_082328671fc89cf6e677dcf29cdcb679effd5880_0_0.fastq -s /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Mouse_Whole_Genome/2_basecalled/barcode12_old/sequencing_summary.txt

nanopolish call-methylation -t 16 -r $sample/fastq_runid_082328671fc89cf6e677dcf29cdcb679effd5880_0_0.fastq -b $sample/bam_runid_082328671fc89cf6e677dcf29cdcb679effd5880_0_0.bam -g /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/mm10.fa  > methylation_calls.tsv

/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/nanopolish/scripts/calculate_methylation_frequency.py methylation_calls.tsv > methylation_frequency.tsv

cat <(head -n1 methylation_frequency.tsv <(tail -n +2 methylation_frequency.tsv | sort -k2,2 -k3,3) | bgzip > methylation_frequency.tsv.gz)

source activate nanopore
methplotlib -m methylation_calls.tsv methylation_frequency.tsv \
            -n calls frequencies \
            -w chr10:108264784-108265842

            chr11	58978473	58978531


### slow5
export PATH=$PATH:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/slow5tools
slow5tools f2s add313_pass_barcode12_08232867_22.fast5 -o add313_pass_barcode12_08232867_22.slow5
guppy_basecaller -i add313_pass_barcode12_08232867_22.slow5 -s . -c dna_r9.4.1_450bps_modbases_5mc_cg_hac.cfg --bam_out --recursive --align_ref $GENOME_FASTA

barcode=(barcode09 barcode10 barcode11 barcode12)
dataset=(20200929_1506_2G_PAE62275_de0204be 20200808_1450_1F_PAF13790_4c803d35 20200929_1506_2G_PAE62275_de0204be)
RAW_FAST5_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/Aaron/PROMETHION/10074/CORTEX
RAW_SLOW5_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Mouse_Whole_Genome/1_raw
for d in ${dataset[@]}; do
  echo $d
  mkdir -p $RAW_SLOW5_DIR/$d
  #slow5tools f2s $RAW_FAST5_DIR/20200929_1506_2G_PAE62275_de0204be/fast5_fail/$d -d $RAW_SLOW5_DIR/$d  -p 8
  slow5tools f2s $RAW_FAST5_DIR/20200929_1506_2G_PAE62275_de0204be/fast5_pass/$d -d $RAW_SLOW5_DIR/$d  -p 8
done

for count in {0..3}; do
  bc=${barcode[count]}
  echo $bc

  for d in ${dataset[@]}; do
    echo $d
    mkdir -p $RAW_SLOW5_DIR/$bc $RAW_SLOW5_DIR/$bc/fail
    slow5tools f2s $RAW_FAST5_DIR/$d"/fast5_pass/"$bc -d $RAW_SLOW5_DIR/$d"/"$bc  -p 8
    slow5tools f2s $RAW_FAST5_DIR/$d"/fast5_fail/"$bc -d $RAW_SLOW5_DIR/$d"/"$bc"/fail" -p 8
    mv $RAW_SLOW5_DIR/$d/$bc/fail/* ..
  done
done
