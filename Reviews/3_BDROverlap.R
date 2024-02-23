
library("ggplot2")
LOGEN_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
LOGEN="/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
SC_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510/Paper_Figures/"
source(paste0(SC_ROOT, "0_source_functions.R"))
source(paste0(SC_ROOT, "rTg4510_config.R"))
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Mouse_Isoseq/"

# note removed outliers in differential expression analysis BDR


## ---------- BDR input data -----------------

BDRONTclass <- read.table("/lustre/recovered/Research_Project-MRC148213/sl693/AD_BDR/D_ONT/5_cupcake/7_sqanti3/ontBDR_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt", sep = "\t", as.is = T)
BDRONT <- readRDS("/lustre/recovered/Research_Project-MRC148213/sl693/AD_BDR/01_figures_tables/Ont_DESeq2TranscriptLevel.RDS")
phenotype <- read.csv("/lustre/recovered/Research_Project-MRC148213/sl693/AD_BDR/0_metadata/B_ONT/Selected_ONTTargeted_BDR.csv", header = T)
TargetGene <- toupper(c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7",
                        "Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf"))

bdrGtf <- as.data.frame(rtracklayer::import("/lustre/recovered/Research_Project-MRC148213/sl693/AD_BDR/D_ONT/5_cupcake/7_sqanti3/ontBDR_collapsed.filtered_counts_filtered.gtf"))
refGtf <- as.data.frame(rtracklayer::import("/lustre/projects/Research_Project-MRC148213/lsl693/references/human/TargetGenes.gencode.v40.annotation.gtf"))
gtf$humanMerged <- rbind(bdrGtf [,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                         refGtf [,c("seqnames","strand","start","end","type","transcript_id","gene_id")])

BDRNormCountsTranscripts <- BDRONT$B2WaldBraak$norm_counts %>% filter(grepl("B2", sample))
lBDRNormCountsTranscripts <- BDRONT$B2WaldBraak$norm_counts_all %>% filter(grepl("B2", sample)) %>% select(sample,isoform,normalised_counts)
lBDRNormCountsTranscripts <- tidyr::spread(lBDRNormCountsTranscripts, sample, normalised_counts)

## ---------- IF across BDR genes -----------------

merged <- tabulateIF(BDRONTclass, "B2")
ggplot(merged, aes(x = associated_gene, y = as.numeric(perc), fill = forcats::fct_rev(structural_category))) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +
  #scale_color_manual(values = rep(NA, length(unique(minorgrouped$gene)))) + 
  labs(x = "Gene", y = "Isoform fraction (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name = "Isoform Classification", values = rev(c(alpha("#00BFC4",0.8),alpha("#00BFC4",0.3),
                                                                    alpha("#F8766D",0.8),alpha("#F8766D",0.3)))) +
  theme(legend.position = "None")


## ---------- Trem2 -----------------

pBdrTrem2 <- list(
  trans = BDR_plot(norm_counts=BDRNormCountsTranscripts, iso ="PB.81888.25", sampleExclude = "BBN00229416.1"),
  gene = BDR_plot(norm_counts=BDRNormCountsTranscripts, iso = NULL, gene ="PB.81888", sampleExclude = "BBN00229416.1"),
  IF = BDR_plot(norm_counts=lBDRNormCountsTranscripts, iso="PB.81888.25", IF = TRUE, Braak = TRUE), 
  track = ggTranPlots(gtf$humanMerged,BDRONTclass,isoList = c("ENST00000373113.8","PB.81888.25"),simple=TRUE, colour = c("black","red"))
)
plot_grid(plotlist = pBdrTrem2, labels = c("A","B","C","D"))


## ---------- Clu -----------------

pBDRClu <- list(
  trans1 = BDR_plot(norm_counts=BDRNormCountsTranscripts, iso="PB.92671.402", sampleExclude = NA),
  trans2 = BDR_plot(norm_counts=BDRNormCountsTranscripts, iso="PB.92671.402", Braak = TRUE),
  gene = BDR_plot(norm_counts=BDRNormCountsTranscripts, gene="PB.92671"),
  IF = BDR_plot(norm_counts=lBDRNormCountsTranscripts, iso="PB.92671.402", IF = TRUE, Braak = TRUE),
  track = ggTranPlots(gtf$humanMerged,BDRONTclass,isoList = c("ENST00000316403.15", "PB.92671.402"),simple=TRUE, colour = c("black","red"))
)
plot_grid(plot_grid(pBDRClu$trans1, pBDRClu$trans2, rel_widths = c(0.4,0.6), labels = c("i","ii")),
          pBDRClu$gene, pBDRClu$IF, pBDRClu$track, labels = c("A","B","C","D"), scale = 0.9)

CluNE <- paste0("PB.92671.",c("689","1693","3228","4691","3474","1988","1778","1875","1860","3549","4725"))
CluAF <- paste0("PB.92671.",c("2479","26717"))

pBDRCluTracks <- list(
  NE = ggTranPlots(gtf$humanMerged,BDRONTclass,isoList = c("PB.92671.402",CluNE ),simple=TRUE, colour = c("black",rep("red",12))),
  ggTranPlots(gtf$humanMerged,BDRONTclass,isoList = c("PB.92671.402",CluAF),simple=TRUE, colour = c("black",rep("red",12)))
)
plot_grid(plotlist = pBDRCluTracks, nrow = 2, rel_heights = c(0.8,0.2), labels = c("A","B"))



## Bin1
pBDRBin1 <- list(
  trans1 = BDR_plot(norm_counts=BDRNormCountsTranscripts, iso="PB.50706.13") + labs(title = "PB.50706.13"),
  trans2 = BDR_plot(norm_counts=BDRNormCountsTranscripts, iso="PB.50706.13", Braak = TRUE) + labs(title = "PB.50706.13"),
  trans3 = BDR_plot(norm_counts=BDRNormCountsTranscripts, iso="PB.50706.8") + labs(title = "PB.50706.8"),
  trans4 = BDR_plot(norm_counts=BDRNormCountsTranscripts, iso="PB.50706.8", Braak = TRUE) + labs(title = "PB.50706.8"),
  gene = BDR_plot(norm_counts=BDRNormCountsTranscripts, gene="PB.50706"),
  IF1 = BDR_plot(norm_counts=lBDRNormCountsTranscripts, iso="PB.50706.8", IF = TRUE, Braak = TRUE),
  IF2 = BDR_plot(norm_counts=lBDRNormCountsTranscripts, iso="PB.50706.8", IF = TRUE, Braak = FALSE),
  track = ggTranPlots(gtf$humanMerged,BDRONTclass,isoList = c("ENST00000316724.10","ENST00000409400.1","PB.50706.8", "PB.50706.13"),simple=TRUE, colour = c("black","black",rep("red",12)))
)

pdf("Bin1.pdf", width = 10, height = 20)
plot_grid(plot_grid(pBDRBin1$trans1, pBDRBin1$trans2, rel_widths = c(0.4,0.6), labels = c("i","ii")),
          plot_grid(pBDRBin1$trans3, pBDRBin1$trans4, rel_widths = c(0.4,0.6), labels = c("i","ii")),
          pBDRBin1$gene, 
          plot_grid(pBDRBin1$IF1, pBDRBin1$IF2, labels = c("i","ii")), 
          pBDRBin1$track, labels = c("A","B","C","D"), scale = 0.9,ncol=1)
dev.off()

BDR_plot(norm_counts=lBDRNormCountsTranscripts, iso="PB.50706.8", IF = TRUE)


## Final 
hAPOE = BDRONTclass[BDRONTclass$associated_gene == "APOE","isoform"]
APOE_p = ggTranPlots(inputgtf=gtf$humanMerged, classfiles=BDRONTclass,
                     isoList = hAPOE[hAPOE != "PB.45429.68"], gene = "APOE")
