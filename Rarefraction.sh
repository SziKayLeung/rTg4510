#!/bin/sh

##############################################################################################################
# Date: 7th May 2019
# Rarefraction Curves after MatchAnnot 
#############################################################################################################

# determine path directory for input data
MATCHANNOT=/mnt/data1/Szi/softwares/MatchAnnot
REFERENCE=/mnt/data1/Szi/reference
SAM=/mnt/data1/Szi/IsoSeq_WholeTranscriptome/ToFU
CUPCAKE=/home/sLeung/cDNA_Cupcake/annotation

SAMPLES_NAMES=(L22 K18 S18 K17 O23)
source activate anaCogent
#############################################################################################################
# RareFaction using output from Matchannot 
cd /mnt/data1/Szi/IsoSeq_WholeTranscriptome/RareFaction/MatchAnnot
for sample in "${SAMPLES_NAMES[@]}"; do 
    echo "Processing $sample file after collapse for mapping"
    if time gmap -D $REFERENCE -d GRCm38.p4 -f samse -n 0 -t 12 -z sense_force        $SAM/$sample.collapsed.filtered.rep.fq  > $sample.fq.sam ;then
        echo "Mapped $sample successful"
    else
        echo "Mapped $sample failed"
   fi
   sort -k 3,3 -k 4,4n $sample.fq.sam > $sample.fq.sorted.sa
   python $MATCHANNOT/MA_works.py --gtf=$REFERENCE/gencode.vM20.primary_assembly.annotation.gtf
   $sample.fq.sorted.sam > $sample.matchannot.own.txt 
   python $CUPCAKE/parse_matchAnnot.py $SAM/$sample.collapsed.filtered.rep.fq $sample.matchannot.own.txt
   python $CUPCAKE/make_file_for_subsampling_from_collapsed.py -i $SAM/$sample.collapsed -o $sample.for_subsampling.txt -m1 $sample.matchannot.own.txt.parsed.txt 
   python $CUPCAKE/subsample.py --by pbgene --min_fl_count 2 --step 100 $sample.for_subsampling.txt  > $sample.by_pbgene.rarefaction.
   echo "Rarefaction of Sample $sample using MatchAnnot complete"
   python $CUPCAKE/subsample.py --by refgene --min_fl_count 2 --step 100 $sample.for_submsampling.txt  > $sample.by_refgene.rarefaction.txt
done