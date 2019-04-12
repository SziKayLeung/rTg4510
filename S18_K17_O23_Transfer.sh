#$ -cwd -V -pe smp 32
<<<<<<< HEAD
# Date: 10th April 2019 
# Transfer S18, K17, O23 Files from Zeus to Isca

# S18_Tg4510_TG_2mos
scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/2_D11/m54082_190404_101400.subreads.bam sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data
=======

# Date: 12th April 2019 
# Transfer S18, K17, O23 Files from ISCA to ISCA personal scratch drive (to work from as path only recently established)

# L22_Tg4510_TG_8months
# scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190304_153723/3_C12/m54082_190306_083150.subreads.bam sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata 
scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190304_153723/3_C12/m54082_190306_083150.subreads.bam.pbi sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata 
# scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190304_153723/3_C12/m54082_190306_083150.subreadset.xml sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata 

# K18_Tg4510_WT_2months
# scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190304_153723/4_D12/m54082_190307_045507.subreads.bam sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata
scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190304_153723/4_D12/m54082_190307_045507.subreads.bam.pbi sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata
# scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190304_153723/4_D12/m54082_190307_045507.subreadset.xml sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata

# S18_Tg4510_TG_2mos
scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/2_D11/m54082_190404_101400.subreads.bam sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata
scp sl693@zeus.ex.ac.uk: /bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/2_D11/m54082_190404_101400.subreads.bam.pbi sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata
scp sl693@zeus.ex.ac.uk: /bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/2_D11/m54082_190404_101400.subreadset.xml sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata 
>>>>>>> 09754be2232a5ab291fd687d644a2a506386a54c

scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/2_D11/m54082_190404_101400.subreads.bam.pbi sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/2_D11/m54082_190404_101400.subreadset.xml sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

# K17_Tg4510_WT_2mos
<<<<<<< HEAD
scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/3_E11/m54082_190405_063832.subreads.bam sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/3_E11/m54082_190405_063832.subreads.bam.pbi sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data
=======
scp sl693@zeus.ex.ac.uk: /bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/3_E11/m54082_190405_063832.subreads.bam sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata
scp sl693@zeus.ex.ac.uk: /bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/3_E11/m54082_190405_063832.subreads.bam.pbi sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata
scp sl693@zeus.ex.ac.uk: /bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/3_E11/m54082_190405_063832.subreadset.xml sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata
>>>>>>> 09754be2232a5ab291fd687d644a2a506386a54c

scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/3_E11/m54082_190405_063832.subreadset.xml sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

# O23_Tg4510
<<<<<<< HEAD
scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190401_162910/1_A01/m54082_190401_165425.subreads.bam sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190401_162910/1_A01/m54082_190401_165425.subreads.bam.pbi sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190401_162910/1_A01/m54082_190401_165425.subreadset.xml sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data
=======
scp sl693@zeus.ex.ac.uk: /bioseqfs/PacBio/Sequel/2801/r54082_20190401_162910/1_A01/m54082_190401_165425.subreads.bam sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata
scp sl693@zeus.ex.ac.uk: /bioseqfs/PacBio/Sequel/2801/r54082_20190401_162910/1_A01/m54082_190401_165425.subreads.bam.pbi sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata
scp sl693@zeus.ex.ac.uk: /bioseqfs/PacBio/Sequel/2801/r54082_20190401_162910/1_A01/m54082_190401_165425.subreadset.xml sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/rawdata
>>>>>>> 09754be2232a5ab291fd687d644a2a506386a54c

echo transfer done