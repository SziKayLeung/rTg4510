#!/usr/bin/env python

import sys
import swan_vis as swan
import pandas as pd

# initialize a new SwanGraph
sg = swan.SwanGraph() 

# add the annotation gtf 
wholetranscriptome_gtf = '/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/SQANTI_TAMA_FILTER/WholeIsoSeq_sqantitamafiltered.classification.gtf'
annot_gtf = '/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/gencode.vM22.annotation.gtf'
sg.add_annotation(annot_gtf)

# set variables 
whole_gtf = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/"
abundance_file="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/swan_count_simple.tsv"
output_dir="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Vis2"
novelinfo_dir="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Vis/InputFiles/"
WT_2mos_gtf = whole_gtf + "WT_2mos" + "_sqantitamafiltered.final.classification.gtf"
TG_2mos_gtf = whole_gtf + "TG_2mos" + "_sqantitamafiltered.final.classification.gtf"
WT_8mos_gtf = whole_gtf + "WT_8mos" + "_sqantitamafiltered.final.classification.gtf"
TG_8mos_gtf = whole_gtf + "TG_8mos" + "_sqantitamafiltered.final.classification.gtf"

# add datasets 
sg.add_dataset('WT 2 months', WT_2mos_gtf)
sg.add_abundance(abundance_file, count_cols='WT_2mos', dataset_name='WT 2 months', tid_col='transcript_id')

sg.add_dataset('WT 8 months', WT_8mos_gtf)
sg.add_abundance(abundance_file, count_cols='WT_8mos', dataset_name='WT 8 months', tid_col='transcript_id')

sg.add_dataset('TG 2 months', TG_2mos_gtf)
sg.add_abundance(abundance_file, count_cols='TG_2mos', dataset_name='TG 2 months', tid_col='transcript_id')

sg.add_dataset('TG 8 months', TG_8mos_gtf)
sg.add_abundance(abundance_file, count_cols='TG_8mos', dataset_name='TG 8 months', tid_col='transcript_id')

sg.save_graph('swan')
sg = swan.SwanGraph('swan.p')

# drop the existing, uninformative novelty category
sg.t_df.drop('novelty', axis=1, inplace=True)

# read novelty CSV with format transcript_id,novelty_category
novelty_df = pd.read_csv(novelinfo_dir + 'All_novelty.csv', names=['tid', 'novelty'])
sg.t_df = sg.t_df.merge(novelty_df, how='left', left_index=True, right_on='tid')


siggenes = ["Gfap","Cd34","Cd68","Osmr","C4b","Ubqln1","Gjb2","Adam23"]
for genes in siggenes: 
  print("Generating report:", genes)
  sg.gen_report(genes, prefix=output_dir + "/" + genes + "/" + genes, heatmap=True, novelty=True, indicate_novel=True)

# All Datasets
sg.add_dataset('Whole Transcriptome', WT_2mos_gtf)
sg.plot_transcript_path('TALONT000301953', indicate_novel=True)
