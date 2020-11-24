# Szi Kay Leung
# 07/10/2019: Expression and Lengths of Target Genes from Gencode as Reference

#********************** Packages (install if not found)
list_of_packages <- c("stringr","dplyr","seqinr","CePa","ggplot2","ggrepel")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")

suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(seqinr))
suppressMessages(library(CePa))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
library(gridExtra)
library(ggsave)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("Rgraphviz")


#********************** Data input and Variables defined
output_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted"
ref_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019"
RSEM_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/RSEM/TEST"
## Input gencode genome gtf, required to extract relevant Transript_Ids from target genes
M22_gtf <- read.csv(print(paste0(ref_dir,"/gencode.vM22_gene_annotation_table.csv")))
head(M22_gtf)

## Input gencode transcript fa, required to extract length of transcripts from Transcript_Ids
# gencode.vM23.transcripts.fa downloaded from gencode (ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.transcripts.fa.gz)
M22_transcripts_fa <- read.fasta(print(paste0(ref_dir,"/gencode.vM23.transcripts.fa")))
M22_Transcripts <- data.frame(Fragments=names(M22_transcripts_fa), Seqs=unlist(getSequence(M22_transcripts_fa, as.string=T)))

## Input RNASeq RSEM results at isoform level, required for expression of transcripts (Mouse)
K17_RSEM <- read.table(print(paste0(RSEM_dir,"/K17_paired_end_quals.isoforms.results")))
colnames(K17_RSEM) <- as.character(unlist(K17_RSEM[1,]))
K17_RSEM <- K17_RSEM[-1, ]

## Specify name of targeted genes
targeted_genes <- c("Bin1","Trem2", "Cd33", "Vgf", "Fyn","Mapt","Trpa1", "Picalm", 
                    "Sorl1", "Abca7","Snca","Apoe","Abca1","App","Ank1","Clu", 
                    "Fus","Ptk2b","Rhbdf2","Tardbp")


#********************** Data Formating
## Remove extra words from column to only unique Transcript ID and Gene Names
M22_gtf$GeneName <- word(M22_gtf$GeneSymbol,c(3),sep = fixed (' '))
M22_Transcripts$Id <- word(M22_Transcripts$Fragments,c(1),sep = fixed ('|'))

# Checked for no duplicated Transcript_Ids 
which(duplicated(M22_Transcripts$TranscriptId) == "TRUE")


#********************** Collating and Merging relevant to target genes
## Extract from gencode genome gtf file, only rows containing target genes 
Targeted_gtf <- data.frame()
for(i in targeted_genes){
  df <- M22_gtf[which(M22_gtf$GeneName == i),]
  Targeted_gtf <- rbind(Targeted_gtf,df)
}
Targeted_gtf$Transcriptid <- word(Targeted_gtf$Transcriptid, c(3), sep = fixed (' '))

## Merge Targeted gtf file with transcript.fasta file to obtain corresponding fasta sequence (and subsequently length)
## Merge by Transcript Id
## all proteing coding?
Final_Transcripts <- merge(Targeted_gtf, M22_Transcripts, by.x = "Transcriptid", by.y = "Id")

# Reorder Final_Transcripts file and remove
colnames(Final_Transcripts)
Final_Transcripts <- Final_Transcripts[,c("Transcriptid","GeneName","Chromosome","Transcript_Length")]

# Tabulating and Plotting number of Transcripts per Gene
table(Final_Transcripts[,"GeneName"])   

p <- ggplot(Final_Transcripts, aes(GeneName)) + 
  geom_bar() + 
  labs(x = "Genes", y = "Number of Transcripts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))
  
print(p)


#********************** Transcript Length of Target genes 
## Calculate Length of Targeted transcripts by counting character of fasta sequence
Final_Transcripts$Seqs <- as.character(Final_Transcripts$Seqs)
Final_Transcripts$Transcript_Length <- nchar(Final_Transcripts$Seqs)

## Histogram of Transcript Lengths
hist(Final_Transcripts$Transcript_Length)

## Analysing number of transcripts with length > 4Kb (cut off point with no size selection)
# Number of Transcripts > 4000bp
dim(Final_Transcripts[Final_Transcripts$Transcript_Length > 4000,])[1]
# % of Transcripts > 4000bp
dim(Final_Transcripts[Final_Transcripts$Transcript_Length > 4000,])[1]/dim(Final_Transcripts)[1] * 100

# Tabulating and Plotting number of transcripts greater than 4000bp
table(Final_Transcripts[Final_Transcripts$Transcript_Length > 4000,"GeneName"])

Final_Transcripts$FourKb <- ifelse(Final_Transcripts$Transcript_Length < 4000,"Less","More")
Final_Transcripts$FourKb <- factor(Final_Transcripts$FourKb, levels = c("More","Less"))


p <- ggplot(Final_Transcripts, aes(GeneName)) + 
  geom_bar(aes(fill = FourKb) ) + 
  labs(x = "Genes", y = "Number of Transcripts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) 

print(p)

###************************************* Trasncript Expression of Target genes from RNASeq data (mouse)

## Check for no duplicated transcript
which(duplicated(K17_RSEM$transcript_id) == "TRUE")

## Merge Final_Transcript file with transcript length with Expression data, based on Transcript Ids
Expression_Transcripts <- merge(K17_RSEM, Final_Transcripts, by.x = "transcript_id", by.y = "Transcriptid", all.y = TRUE)

## Create function to plot Transcript Length and Expression of each Target Gene 
plot_per_gene <- function(gene){
  Expression_Transcripts$TPM <- as.numeric(as.character(Expression_Transcripts$TPM))
  Expression_Transcripts$IsoPct <- as.numeric(as.character(Expression_Transcripts$IsoPct))
  
  df <- Expression_Transcripts[Expression_Transcripts$GeneName == gene,]
  #print(head(df))
  
  ggplot(df, aes(x = Transcript_Length, y = log(TPM), label = transcript_id)) + 
    geom_point(aes(size = IsoPct),colour = "purple", alpha = 0.5) +
    geom_text(vjust = "inward", hjust = "inward", check_overlap = TRUE) +
    theme_classic() + 
    labs(title = print(paste(gene,":",length(df$transcript_id),"Transcripts"))) +
    xlim(0, 10000) + 
    ylim(-5,7) +
    geom_vline(xintercept = 4000, linetype="dotted", color = "red")
    
}

## Plot each gene and save 4 plots per page, by first saving plots to a list 
p = list()
counter = 1
for(i in targeted_genes){
  p[[counter]] <- plot_per_gene(i)
  counter = counter + 1 
}
ml <- marrangeGrob(p, nrow=2, ncol=3)
ggsave(print(paste0(output_dir,"/Target_Expression_Length.pdf")),ml,width = 30, height = 30, units = "cm",dpi = 300)

###************************************* Expression from GTEX (human)
#GTEX <- read.gct(print(paste0(ref_dir,"/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct")))
#head(GTEX)
