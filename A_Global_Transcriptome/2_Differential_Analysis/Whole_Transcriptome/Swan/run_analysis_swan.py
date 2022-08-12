#!/usr/bin/env python

import sys
import swan_vis as swan
import pandas as pd

# initialize a new SwanGraph
sg = swan.SwanGraph() 

# add the annotation gtf 
#annot_gtf = '/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/WholeIsoSeq_sqantitamafiltered.classification.gtf'
annot_gtf = '/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/gencode.vM22.annotation.gtf'
sg.add_annotation(annot_gtf)

config_df = pd.read_csv('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/swan_wholecompile.tsv', sep='\t')
config_df

sg.add_datasets('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/swan_wholecompile.tsv')
dataset_groups = [['K17_WT_2mos','Q21_WT_2mos','M21_WT_2mos'], ['K18_TG_2mos','S18_TG_2mos','O18_TG_2mos']]


# perform a differential gene expression 
# Wald test on the provided two lists of datasets
sg.de_gene_test(dataset_groups)
