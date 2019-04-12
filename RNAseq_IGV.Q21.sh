cd /mnt/data1/Szi/Q21
mkdir Q21_RNAseq_align
cd Q21_RNAseq_align

STAR --runThreadN 32 --genomeDir /mnt/data1/Szi/reference/STAR \
--readFilesIn /mnt/data1/isabel/RNA-seq/analysis_isabel_new/Tg4510_filtered/Q21/Q21_S6_R1_001.fastq.filtered \
/mnt/data1/isabel/RNA-seq/analysis_isabel_new/Tg4510_filtered/Q21/Q21_S6_R2_001.fastq.filtered

echo Q21 STAR Alignment Successful

samtools view -S -b Aligned.out.sam > Aligned.out.bam
samtools sort Aligned.out.bam Aligned.out.sorted
echo SAMTOOLS sort Successful
samtools index Aligned.out.sorted.bam
echo SAMTOOLS index Successful
samtools flagstat Aligned.out.sorted.bam > mappingstats.txt

mv Aligned.out.sorted.bam Q21_Tg4510_2mos_WT.sorted.bam
mv Aligned.out.sorted.bam.bai Q21_Tg4510_2mos_WT.sorted.bam.bai
rm Aligned.out.sam
