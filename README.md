# Isoseq3_Tg4510

This directory is a repository of scripts pertaining to Whole Transcriptome and Targeted Transcriptome Iso-Seq analysis of Tg4510 mice

## **Raw Data**
* Whole Transcriptome: List of location of rawdata on ISCA to feed into Iso-Seq Analysis pipeline (Linux commands)
* Targeted Transcriptome
    + List of location of rawdata on ISCA
    + Barcode Configuration Files: Samples with the relevant Pacbio barcode for demultiplexing
    + Probes: List of probes and locations for determining off-target rate

## **Datasets** 
There are 2 Iso-Seq datasets:
1. **Whole Transcriptome (n = 12 samples)**
2. **Targeted Transcriptome (n = 24 samples)** 

RNASeq: Short-read RNA-Seq data is used to 1) support junction sequence (splice site) from long-read Iso-Seq data from *STAR* output (SJ.bed) and 2) validate isoforms from RNA-Seq alignment using *Kallisto* - which is all integrated into *SQANTI*. 

---
## 1. **Transcriptome_Annotation**:
Iso-Seq Analysis pipeline accurately annotates the whole and targeted transcriptome. Pipeline starts from IsoSeq3 tools (*CCS, lima, refine, cluster*), alignment to mouse transcriptome using *minimap2*, collapse and chaining of multiple samples with *Cupcake*, to transcriptome annotation with *SQANTI* using multiple input from RNA-Seq data, CAGE peaks and polyA motifs. All analyses pertaining to the whole transcriptome annotation of the whole transcriptome can be found in the Whole_Transcriptome_Paper directory. The bioinformatic pipeline for analysing the targeted transcriptome dataset is the same as the whole transcriptome, with the exception of demultiplexing barcodes at *lima*. 

For both Iso-Seq datasets, all the respective samples were merged at *Iso-Seq3 Refine* to generate one big reference annotation. The abundance for each of the individual samples in whole transcriptome was then obtained from *cupcake*'s read_stat.txt, which tabulates all the FL transcripts, with unique id specific to the sample (sequel run id), associated to each isoform. Similarly, obtaining abundance for each of the barcoded samples in the targeted transcriptome ws obtained from *cupcake* with the addition of the refine report from the individual samples (before merging, but after demultiplexing) for differentiation of transcripts by the barcoded samples. 

---
## 2. **Transcriptome Characterisation**:


---
## 3. **Differential Analysis**:

Differential Analysis involves investigating changes at the gene and isoform level associated with rTg4510 pathology and progressively over age. Both datasets generated from the whole and targeted transcriptome approach were explored, and isoforms were visualised with Swan.

Different methods were trialled: IsoformSwitchAnalyzeR, DESeq2, pipeline from ONT and custom scripts (using non-parametic methods). In the end, we worked with results from tappAS as the tool was designed for time-series experiments, with the ability to investigate differential feature usage.

