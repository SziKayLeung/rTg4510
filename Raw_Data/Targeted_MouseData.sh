#!/bin/bash

# Project_10040 - Iso-Seq Targeted Batch 1
# Project_10175 - Iso-Seq Targeted Batch 2 and Batch 3

# 22/10/2020: Cataloguing raw files and barocoded information on Targeted IsoSeq Mouse 


FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/targeted_isoseq/Raw_Data

#************************************* SUBREADS FILES
cd $FUNCTIONS
cat << EOF > $FUNCTIONS/Isoseq_Targeted_MouseRaw.txt 
#1.Targeted_Sequencing_1: K20, K18, K23, K21, K19, K17 (Project_10040)
/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/rawdata/Project_10040/02_raw_pb_reads/r54082_20191115_164308/2_B01/m54082_191116_131337.subreads.bam
#2.Targeted_Sequencing_2: S19, K24, L22, M21, O18, O23, O22, P19, T20 (Project_10175)
/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/rawdata/Project_10175/Batch2/r54082_20200730_135358/3_C01/m54082_200731_163617.subreads.bam
#3a.Targeted_Sequencing_3: Q20, Q21, S18, S23, Q18, Q17, L18, Q23, T18 (partial run, Project_10175)
/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/rawdata/Project_10175/Batch3/r54082_20200730_135358/4_D01/m54082_200801_130641.subreads.bam
#3b.Targeted_Sequencing_3: Q20, Q21, S18, S23, Q18, Q17, L18, Q23, T18 (complete run, Project_10175)
/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/rawdata/Project_10175/Batch3/r54082_20200807_141054/3_C06/m54082_200808_064617.subreads.bam
EOF

cat << EOF > $FUNCTIONS/Isoseq_Targeted_MouseRawDir.txt 
#1.Targeted_Sequencing_1: K20, K18, K23, K21, K19, K17 
/bioseqfs/PacBio/Sequel/3087_10040/r54082_20191115_164308/2_B01/m54082_191116_131337.subreadset.xml
/bioseqfs/PacBio/Sequel/3087_10040/r54082_20191115_164308/2_B01/m54082_191116_131337.subreads.bam
/bioseqfs/PacBio/Sequel/3087_10040/r54082_20191115_164308/2_B01/m54082_191116_131337.subreads.bam.pbi
#2.Targeted_Sequencing_2: S19, K24, L22, M21, O18, O23, O22, P19, T20 
/bioseqfs/PacBio/Sequel/10056_10175_10176/r54082_20200730_135358/3_C01/m54082_200731_163617.subreadset.xml
/bioseqfs/PacBio/Sequel/10056_10175_10176/r54082_20200730_135358/3_C01/m54082_200731_163617.subreads.bam
/bioseqfs/PacBio/Sequel/10056_10175_10176/r54082_20200730_135358/3_C01/m54082_200731_163617.subreads.bam.pbi
#3a.Targeted_Sequencing_3: Q20, Q21, S18, S23, Q18, Q17, L18, Q23, T18 
/bioseqfs/PacBio/Sequel/10056_10175_10176/r54082_20200730_135358/4_D01/m54082_200801_130641.subreadset.xml
/bioseqfs/PacBio/Sequel/10056_10175_10176/r54082_20200730_135358/4_D01/m54082_200801_130641.subreads.bam
/bioseqfs/PacBio/Sequel/10056_10175_10176/r54082_20200730_135358/4_D01/m54082_200801_130641.subreads.bam.pbi
#3b.Targeted_Sequencing_3: Q20, Q21, S18, S23, Q18, Q17, L18, Q23, T18 
/bioseqfs/PacBio/Sequel/10094_10175/r54082_20200807_141054/3_C06/m54082_200808_064617.subreadset.xml
/bioseqfs/PacBio/Sequel/10094_10175/r54082_20200807_141054/3_C06/m54082_200808_064617.subreads.bam
/bioseqfs/PacBio/Sequel/10094_10175/r54082_20200807_141054/3_C06/m54082_200808_064617.subreads.bam.pbi
EOF

#************************************* BARCODES 
cat << EOF > $FUNCTIONS/Isoseq_Mouse_Targeted1_barcode.config 
K19=fl.primer_5p--BC1001_3p.bam
K23=fl.primer_5p--BC1002_3p.bam
K21=fl.primer_5p--BC1003_3p.bam
K18=fl.primer_5p--BC1004_3p.bam
K20=fl.primer_5p--BC1005_3p.bam
K17=fl.primer_5p--BC1006_3p.bam
EOF

cat << EOF > $FUNCTIONS/Isoseq_Mouse_Targeted2_barcode.config 
S19=fl.primer_5p--BC1001_3p.bam
K24=fl.primer_5p--BC1002_3p.bam
L22=fl.primer_5p--BC1003_3p.bam
M21=fl.primer_5p--BC1004_3p.bam
O18=fl.primer_5p--BC1005_3p.bam
O23=fl.primer_5p--BC1006_3p.bam
O22=fl.primer_5p--BC1007_3p.bam
P19=fl.primer_5p--BC1008_3p.bam
T20=fl.primer_5p--BC1009_3p.bam
EOF

cat << EOF > $FUNCTIONS/Isoseq_Mouse_Targeted3_barcode.config 
Q20=fl.primer_5p--BC1001_3p.bam
Q21=fl.primer_5p--BC1002_3p.bam
S18=fl.primer_5p--BC1003_3p.bam
S23=fl.primer_5p--BC1004_3p.bam
Q18=fl.primer_5p--BC1005_3p.bam
Q17=fl.primer_5p--BC1006_3p.bam
L18=fl.primer_5p--BC1007_3p.bam
Q23=fl.primer_5p--BC1008_3p.bam
T18=fl.primer_5p--BC1009_3p.bam
EOF