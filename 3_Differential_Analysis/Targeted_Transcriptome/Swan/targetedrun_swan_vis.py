#!/usr/bin/env python

import sys
import swan_vis as swan

# initialize a new SwanGraph
sg = swan.SwanGraph() 

# add the annotation gtf 
#annot_gtf = '/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/WholeIsoSeq_sqantitamafiltered.classification.gtf'
annot_gtf = '/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/gencode.vM22.annotation.gtf'
sg.add_annotation(annot_gtf)

# set variables 
targeted_gtf = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite.gtf"
#targeted_gtf = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3/"
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Vis"
#abundance_file = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Targeted_Transcriptome/Swan/swan_count_simple.tsv"
#novelinfo_dir="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Targeted_Transcriptome/Swan/"
#WT_2mos_gtf = targeted_gtf + "WT_2mos" + "_sqantitamafiltered.final.classification.gtf"
#TG_2mos_gtf = targeted_gtf + "TG_2mos" + "_sqantitamafiltered.final.classification.gtf"
#WT_8mos_gtf = targeted_gtf + "WT_8mos" + "_sqantitamafiltered.final.classification.gtf"
#TG_8mos_gtf = targeted_gtf + "TG_8mos" + "_sqantitamafiltered.final.classification.gtf"

# add datasets 
sg.add_dataset('Targeted', targeted_gtf)
#sg.add_dataset('WT 2 months', WT_2mos_gtf)
#sg.add_abundance(abundance_file, count_cols='WT_2mos', dataset_name='WT 2 months', tid_col='transcript_id')

#sg.add_dataset('WT 8 months', WT_8mos_gtf)
#sg.add_abundance(abundance_file, count_cols='WT_8mos', dataset_name='WT 8 months', tid_col='transcript_id')

#sg.add_dataset('TG 2 months', TG_2mos_gtf)
#sg.add_abundance(abundance_file, count_cols='TG_2mos', dataset_name='TG 2 months', tid_col='transcript_id')

#sg.add_dataset('TG 8 months', TG_8mos_gtf)
#sg.add_abundance(abundance_file, count_cols='TG_8mos', dataset_name='TG 8 months', tid_col='transcript_id')


#sg.save_graph('swan')
#sg = swan.SwanGraph('swan.p')

# SNCA plots for Thesis to show redundant transcripts
transcripts = ["PB.6948.1","PB.6948.115","PB.6948.15","PB.6948.160","PB.6948.35","PB.6948.81","PB.6948.7"]
for t in transcripts: 
  print("Generating figure:", t)
  sg.plot_transcript_path(t, indicate_novel=True)
  swan.save_fig(output_dir+'/SNCA/'+ t+'_Loop.png')
  sg.plot_transcript_path(t, browser=True)
  swan.save_fig(output_dir+'/SNCA/'+ t+'_Browser.png')



transcripts = ["PB.3915.175","PB.3915.120","PB.3915.5"]
for t in transcripts: 
  print("Generating figure:", t)
  sg.plot_transcript_path(t, indicate_novel=True)
  swan.save_fig(output_dir+ t+'_Loop.png')
  sg.plot_transcript_path(t, browser=True)
  swan.save_fig(output_dir+ t+'_Browser.png')
  

transcripts = ["PB.2634.299","PB.2634.365","PB.2634.2"]
for t in transcripts: 
  print("Generating figure:", t)
  sg.plot_transcript_path(t, indicate_novel=True)
  swan.save_fig(output_dir+'/'+ t+'_Loop.png')
  sg.plot_transcript_path(t, browser=True)
  swan.save_fig(output_dir+'/'+ t+'_Browser.png')
  

transcripts = ["PB.3915.1"]
for t in transcripts: 
  print("Generating figure:", t)
  sg.plot_transcript_path(t, indicate_novel=True)
  swan.save_fig(output_dir+'/'+ t+'_Loop.png')
  sg.plot_transcript_path(t, browser=True)
  swan.save_fig(output_dir+'/'+ t+'_Browser.png')





