# Isoseq3_Tg4510

This directory is a repository of scripts pertaining to Whole Transcriptome and Targeted Transcriptome Iso-Seq analysis of Tg4510 mice

## **Raw Data**
* Whole Transcriptome: List of location of rawdata on ISCA to feed into Iso-Seq Analysis pipeline (Linux commands)
* Targeted Transcriptome
    + List of location of rawdata on ISCA
    + Barcode Configuration Files: Samples with the relevant Pacbio barcode for demultiplexing
    + Probes: List of probes and locations for determining off-target rate


---
## 1. **Transcriptome_Annotation**:
Iso-Seq Analysis pipeline accurately annotates the whole and targeted transcriptome, and is performed either using Linux commands or via snakemake. Pipeline starts from IsoSeq3 tools (*CCS, lima, refine, cluster*), alignment to mouse transcriptome using *minimap2*, collapse and chaining of multiple samples with *Cupcake*, to transcriptome annotation with *SQANTI* using multiple input from RNA-Seq data, CAGE peaks and polyA motifs.

There are 3 main analyses:
1. **Whole Transcriptome - Individual Samples (n = 12 samples)** Global differential transcript expression and usage between WT and TG
    + Snakemake: Snakefile paired with config.yaml
    + Linux: Individual_Part1.sh for parallel analysis (IsoSeq3 tools, mapping and *Cupcake* collapse), Individual_Part2.sh for all chaining and sqanti

    <details>
      <summary>Output:</summary>

      + Sequel run quality: Number of polymerase reads, CCS reads, FL reads
      + Mapping read quality  
      + Rarefaction curves  
      + Read Lengths  
    </details>

<br>

2. **Pooled Samples (n = 3)**: Qualitative differences between pooled WT, TG and All.
    + Snakemake: Snakefile_All paired with config files for mouse transcriptome (config_all_mm10.yaml,config_all_WT.yaml, config_all_TG.yaml) and for ERCC alignment (config_all_ERCC.yaml,config_all_ERCC.yaml, config_all_ERCC.yaml)  
    + Linux: Batches_Part2.sh and Batches_Part2.sh for parallel analysis and chaining, *SQANTI* respectively

    <details>
      <summary>Output:</summary>

      + ERCC detection
      + Sequel run quality, read length differences between WT and TG
    </details>

<br>

3. **Targeted Transcriptome (n = 24)**: Differential transcript expression and usage between WT and TG on panelled genes
   + Snakemake: Snakefile_Targeted paired with config_targeted.yaml
   + Linux:

   <details>
     <summary>Output:</summary>

     + Off-target rate
     + Sample batches differences
   </details>

**RNASeq**: Short-read RNA-Seq data is used to 1) support junction sequence (splice site) from long-read Iso-Seq data from *STAR* output (SJ.bed) and 2) validate isoforms from RNA-Seq alignment using *Kallisto* - which is all integrated into *SQANTI*. All RNA-Seq data (n = 64) is used for junction support, whereas only the specific RNA-Seq data from sequenced sample is used for expression.  

---
