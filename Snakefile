# Szi Kay Leung
# 30/10/2020: Adapted working script from https://github.com/nanoporetech/ont_tutorial_transcriptome.git

import os
import re
from os import path
import pandas as pd
from collections import OrderedDict

###*********** To run *************** 
#snakemake -j 1 --configfile <config.yaml>
#configfile: "config.yaml"
###**********************************  


# extract samples from the configfile
Samples = []
for i in range(len(config["Samples"])):
  dataSlice = config["Samples"][i]
  conditionSamples = list(list(dataSlice.items())[0][1].items())
  for j in range(len(conditionSamples)):
    sequenceFile = conditionSamples[j][1]
    sequenceFile = re.sub("RawData/","",sequenceFile) # this should be abstracted
    #print(sequenceFile)
    Samples.append(sequenceFile)

# Split the filenames into basename and extension
files = [filename.split('.', 1) for filename in Samples]
# Create a dictionary of mapping basename:extension
file_dict = {filename[0]: filename[1] if len(filename) == 2 else '' for filename in files}

rule all:
    input:
      expand("Analysis/samtools/{seqid}.flagstat", seqid=file_dict.keys()), 
      expand("Analysis/Salmon/{seqid}", seqid=file_dict.keys()), 
      config["annotation"]
      

rule build_minimap_index: ## build minimap2 index
    input:
        genome = config["transcriptome"]
    output:
        index = "Analysis/Minimap2/transcriptome_index.mmi"
    params:
        opts = config["minimap_index_opts"]
    threads: config["threads"]
    shell:"""
        minimap2 -t {threads} {params.opts} -I 1000G -d {output.index} {input.genome}
    """

rule map_reads: ## map reads using minimap2
    input:
       index = rules.build_minimap_index.output.index,
       fastq = lambda wc: "RawData/" + wc.seqid + "." + file_dict[wc.seqid]
    output:
       bam = "Analysis/Minimap2/{seqid}.bam",
       sbam = "Analysis/Minimap2/{seqid}.sorted.bam",
    params:
        opts = config["minimap2_opts"],
        msec = config["maximum_secondary"],
        psec = config["secondary_score_ratio"]
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} -ax map-ont -p {params.psec} -N {params.msec} {params.opts} {input.index} {input.fastq}\
    | samtools view -Sb > {output.bam};
    samtools sort -@ {threads} {output.bam} -o {output.sbam};
    samtools index {output.sbam};
    """

rule flagstat:
  input:
    "Analysis/Minimap2/{seqid}.bam"
  output:
    "Analysis/samtools/{seqid}.flagstat"
  shell:
    "samtools flagstat {input} > {output}"


rule count_reads:
    input:
        bam = "Analysis/Minimap2/{seqid}.sorted.bam",
        trs = config["transcriptome"],
    output:
        tsv = directory("Analysis/Salmon/{seqid}"),
    params:
        libtype = config["salmon_libtype"],
    threads: config["threads"]
    shell: """
        salmon quant --noErrorModel -p {threads} -t {input.trs} -l {params.libtype} -a {input.bam} -o {output.tsv}
    """
