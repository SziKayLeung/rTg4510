#$ -cwd -V -pe smp 32
# Date: 10th April 2019 
# Transfer S18, K17, O23 Files from Zeus to Isca

# S18_Tg4510_TG_2mos
scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/2_D11/m54082_190404_101400.subreads.bam sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/2_D11/m54082_190404_101400.subreads.bam.pbi sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/2_D11/m54082_190404_101400.subreadset.xml sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

# K17_Tg4510_WT_2mos
scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/3_E11/m54082_190405_063832.subreads.bam sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/3_E11/m54082_190405_063832.subreads.bam.pbi sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/3_E11/m54082_190405_063832.subreadset.xml sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

# O23_Tg4510
scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190401_162910/1_A01/m54082_190401_165425.subreads.bam sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190401_162910/1_A01/m54082_190401_165425.subreads.bam.pbi sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190401_162910/1_A01/m54082_190401_165425.subreadset.xml sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/Tg4510_raw_data

echo transfer done