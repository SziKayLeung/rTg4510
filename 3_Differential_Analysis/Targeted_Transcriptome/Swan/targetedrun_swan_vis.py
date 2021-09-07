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
targeted_gtf = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/COLLAPSE_FILTER/AllMouseTargeted_sqantisubset.final.classification.gtf"
output_dir="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Vis"

# add datasets 
sg.add_dataset('Targeted Mouse', targeted_gtf)
#sg.save_graph('swan')
#sg = swan.SwanGraph('swan.p')

sg.plot_transcript_path('PB.3742.1', indicate_novel=True)
swan.save_fig('PB.3742.1_Trem2.png')
#sg.plot_graph('Trem2', indicate_novel=True)
#swan.save_fig(output_dir+'/Trem2.png')

#sg.plot_graph('Trem2')

sg.plot_transcript_path('PB.5746.2', indicate_novel=True)
swan.save_fig('PB.5746.2_Abca1.png')

