#!/usr/bin/env python
# Szi Kay Leung

######## Run on command line ISCA
#module load Miniconda2
#source activate sqanti2_py3
#python sqanti_gtf_modify.py /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual_Samples/CHAIN/
#source deactivate
#############################

import pandas as pd
import csv
import sys

SQANTI_dir = sys.argv[1]
gtf_input_file = SQANTI_dir + "/all_samples.chained_classification.filtered_lite.gtf"
class_input_file = SQANTI_dir + "/all_samples.chained_classification.filtered_lite_classification.txt"


# aim: 1. read gtf file and classification file (for dictionary of isoform id and gene name)
#      2. replace the gene_id with the gene name
# input: gtf file, dictionary of pb_id: gene_name
# output: additional column with replaced dataframe
def replace_gene_id(gtf_file, class_file):
    # read Gtf file
    gtf_df = pd.read_csv(filepath_or_buffer = gtf_file, delim_whitespace=True, header=None,
        names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

    # read classification file
    class_df= pd.read_csv(filepath_or_buffer = class_file,sep='\t')

    # create a dictionary with isoform id and gene_file
    id_gene = dict(zip(class_df.isoform, class_df.associated_gene))

    # for formatting in output file
    def double_quote(word):
        return '"%s"' % word

    # create empty list after replacement
    new_attributes = []
    for info in gtf_df.attributes:
        # extracting the labels and names from attributes column
        gene_id_label = info.split(";")[0]
        print(gene_id_label)
        #gene_id = gene_id_label.split(" ")[1]
        # replace the gene_id with the gene_name from dictionary
        gene_name = [val for key, val in id_gene.items() if gene_id_label in key][0]
        print("Renaming", gene_id_label, "to", gene_name)
        new_attributes.append("gene_id " + double_quote(gene_name) + "; transcript_id" + gene_id_label +";")

    output = gtf_df.join(pd.DataFrame(new_attributes, columns = ["new_attributes"]))
    return(output)

def main():
    output = replace_gene_id(gtf_input_file,class_input_file)
    final = output[["seqid","source","type","start","end","score","strand","phase","new_attributes"]]
    final.to_csv(SQANTI_dir + 'all_samples.chained.rep_classification.filtered_lite_mod.gtf',sep='\t',index=False,
                 header=False, quoting=csv.QUOTE_NONE)

main()