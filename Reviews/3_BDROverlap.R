
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
lBDRNormCountsTranscripts <- BDRONT$B2WaldBraak$norm_counts_all %>% dplyr::filter(grepl("B2", sample)) %>% dplyr::select(sample,isoform,normalised_counts)
lBDRNormCountsTranscripts <- tidyr::spread(lBDRNormCountsTranscripts, sample, normalised_counts)

## ---------- IF across BDR genes -----------------

merged <- tabulateIF(BDRONTclass, "B2")
BDRIF <- ggplot(merged, aes(x = associated_gene, y = as.numeric(perc), fill = forcats::fct_rev(structural_category))) +
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

phenotype <- read.csv("/lustre/recovered/Research_Project-MRC148213/sl693/AD_BDR/0_metadata/B_ONT/Selected_ONTTargeted_BDR.csv", header = T)
Trem2TranscriptBDR <- BDRONT$B2WaldBraak$norm_counts %>% filter(isoform == "PB.81888.25") %>% filter(grepl("B2", sample)) %>% 
  mutate(sample = str_remove(sample, "B2.")) %>% 
  left_join(., phenotype, by = "sample") %>%
  filter(BraakTangle_numeric %in% c(0,1,2,5,6)) %>% 
  mutate(phenotype = ifelse(BraakTangle_numeric %in% c(0,1,2),"Control","AD")) %>%
  mutate(phenotype = factor(phenotype, levels = c("Control","AD"))) 
t.test(normalised_counts ~ phenotype, data = Trem2TranscriptBDR %>% filter(sample != "BBN00229416.1"))
summary(lm(normalised_counts ~ phenotype, data = Trem2TranscriptBDR %>% filter(sample != "BBN00229416.1")))

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

plot_grid(plot_grid(pBDRClu$trans1, pBDRClu$trans2),
          plot_transexp_overtime("Clu",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",rank=NULL,
                                 isoSpecific=c("PB.14646.139"),
                                 setorder=c("CONTROL","CASE")), nrow = 2, labels = c("A","B"))


CluTranscriptBDR <- BDRONT$B2WaldBraak$norm_counts %>% filter(isoform == "PB.92671.402") %>% filter(grepl("B2", sample)) %>% 
  mutate(sample = str_remove(sample, "B2.")) %>% 
  left_join(., phenotype, by = "sample") %>%
  filter(BraakTangle_numeric %in% c(0,1,2,5,6)) %>% 
  mutate(phenotype = ifelse(BraakTangle_numeric %in% c(0,1,2),"Control","AD")) %>%
  mutate(phenotype = factor(phenotype, levels = c("Control","AD"))) 
t.test(normalised_counts ~ phenotype, data = CluTranscriptBDR %>% filter(sample != "BBN00229416.1"))
summary(lm(normalised_counts ~ as.factor(BraakTangle_numeric), data = CluTranscriptBDR %>% filter(sample != "BBN00229416.1")))

## Bin1
pBDRBin1 <- list(
  trans1 = BDR_plot(norm_counts=BDRNormCountsTranscripts, iso="PB.50706.13") + labs(title = "LR.BIN1.13"),
  trans2 = BDR_plot(norm_counts=BDRNormCountsTranscripts, iso="PB.50706.13", Braak = TRUE) + labs(title = "LR.BIN1.13"),
  trans3 = BDR_plot(norm_counts=BDRNormCountsTranscripts, iso="PB.50706.8") + labs(title = "LR.BIN1.8"),
  trans4 = BDR_plot(norm_counts=BDRNormCountsTranscripts, iso="PB.50706.8", Braak = TRUE) + labs(title = "LR.BIN1.8"),
  gene = BDR_plot(norm_counts=BDRNormCountsTranscripts, gene="PB.50706"),
  IF1 = BDR_plot(norm_counts=lBDRNormCountsTranscripts, iso="PB.50706.8", IF = TRUE, Braak = TRUE),
  IF2 = BDR_plot(norm_counts=lBDRNormCountsTranscripts, iso="PB.50706.8", IF = TRUE, Braak = FALSE),
  track = ggTranPlots(gtf$humanMerged,BDRONTclass,isoList = c("ENST00000316724.10","ENST00000409400.1","PB.50706.8", "PB.50706.13"),simple=TRUE, colour = c("black","black",rep("red",12)))
)
Bin1TranscriptBDR <- BDRONT$B2WaldBraak$norm_counts %>% filter(isoform == "PB.50706.13") %>% filter(grepl("B2", sample)) %>% 
  mutate(sample = str_remove(sample, "B2.")) %>% 
  left_join(., phenotype, by = "sample") %>%
  filter(BraakTangle_numeric %in% c(0,1,2,5,6)) %>% 
  mutate(phenotype = ifelse(BraakTangle_numeric %in% c(0,1,2),"Control","AD")) %>%
  mutate(phenotype = factor(phenotype, levels = c("Control","AD"))) 
t.test(normalised_counts ~ phenotype, data = Bin1TranscriptBDR  %>% filter(sample != "BBN00229416.1"))
summary(lm(normalised_counts ~ as.factor(BraakTangle_numeric), data = Bin1TranscriptBDR  %>% filter(sample != "BBN00229416.1")))

pdf("Bin1.pdf", width = 10, height = 20)
plot_grid(plot_grid(pBDRBin1$trans1, pBDRBin1$trans2, rel_widths = c(0.4,0.6), labels = c("i","ii")),
          plot_grid(pBDRBin1$trans3, pBDRBin1$trans4, rel_widths = c(0.4,0.6), labels = c("i","ii")),
          pBDRBin1$gene, 
          plot_grid(pBDRBin1$IF1, pBDRBin1$IF2, labels = c("i","ii")), 
          pBDRBin1$track, labels = c("A","B","C","D"), scale = 0.9,ncol=1)
dev.off()

plot_grid(plot_grid(pBDRBin1$trans1, pBDRBin1$trans2, rel_widths = c(0.4,0.6), labels = c("i","ii")),
          plot_grid(pBDRBin1$trans3, pBDRBin1$trans4, rel_widths = c(0.4,0.6), labels = c("i","ii")))


plot_grid(plot_grid(pBDRBin1$trans1, pBDRBin1$trans2),
          plot_transexp_overtime("Bin1",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",rank=NULL,
                       isoSpecific=c("PB.22007.224"),
                       setorder=c("CONTROL","CASE")), nrow = 2, labels = c("A","B"))

plot_grid(plot_grid(pBDRBin1$trans3, pBDRBin1$trans4),
          plot_transexp_overtime("Bin1",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",rank=NULL,
                                 isoSpecific=c("PB.22007.99"),
                                 setorder=c("CONTROL","CASE")), nrow = 2, labels = c("A","B"))


## Final 
novelAPOE <- BDRONTclass[BDRONTclass$associated_gene == "APOE" & BDRONTclass$structural_category %in% c("NIC","NNC"),"isoform"]
APOEIso <- data.frame(
  Isoform = unlist(APOEIso <- list(
    Reference = c("ENST00000446996.5","ENST00000252486.9", "ENST00000434152.5", "ENST00000425718.1"),
    `FSM ISM` = BDRONTclass[BDRONTclass$associated_gene == "APOE" & BDRONTclass$structural_category %in% c("FSM","ISM"),"isoform"],
    `NIC NNC` = novelAPOE[novelAPOE != "PB.45429.68"]
  )),
  Category = rep(names(APOEIso), lengths(APOEIso))
)
APOEIso$colour <- c(rep(NA,length(APOEIso$Category)))


APOE_p = ggTranPlots(inputgtf=gtf$humanMerged, classfiles=BDRONTclass,
                     isoList = c(as.character(APOEIso$Isoform)), selfDf = APOEIso, gene = "APOE")


CLUIso <- data.frame(
  Isoform = unlist(CLUIso <- list(
    Reference = c("ENST00000316403.15","ENST00000405140.7","ENST00000522238.1","ENST00000523500.5"),
    `NE` = paste0("PB.92671.",c("689","1693","3228","3474","1988","1778","1860","3549","4725")),
    `AF` = paste0("PB.92671.",c("2479","4691","1875","26717"))
  )),
  Category = rep(names(CLUIso), lengths(CLUIso))
)
CLUIso$colour <- c(rep(NA,length(CLUIso$Category)))
CLU_p <- ggTranPlots(inputgtf=gtf$humanMerged, classfiles=BDRONTclass,
            isoList = c(as.character(CLUIso$Isoform)), selfDf = CLUIso, gene = "CLU") 

Trem2countsA <- BDR_plot(norm_counts=BDRNormCountsTranscripts, iso ="PB.81888.25", sampleExclude = "BBN00229416.1") + 
  scale_fill_manual(values = c(label_colour("WT"),label_colour("AD"))) + theme(legend.position = "None") +
  labs(title = "", subtitle = "ONT Transcript Expression")

Trem2countsB  <- BDR_plot(norm_counts=BDRNormCountsTranscripts, iso = NULL, gene ="PB.81888", sampleExclude = "BBN00229416.1") + 
  scale_fill_manual(values = c(label_colour("WT"),label_colour("AD"))) + theme(legend.position = "None") +
  labs(title = "", subtitle = "ONT Gene Expression")

BDRTrem2Iso <- data.frame(
  Isoform = unlist(BDRTrem2Iso <- list(
    `Mouse` = c("ENSMUST00000024791.14", "PB.20818.54"),
    `Human` = c("ENST00000373113.8", "PB.81888.25")
  )),
  Category = rep(names(BDRTrem2Iso), lengths(BDRTrem2Iso))
)
Trem2A <- ggTranPlots(inputgtf= gtf$targ_merged, classfiles=class.files$targ_filtered,
            isoList = c(as.character(BDRTrem2Iso$Isoform)), selfDf = BDRTrem2Iso, gene = "Trem2", squish=FALSE) 
TREM2B <- ggTranPlots(inputgtf= gtf$humanMerged, classfiles=BDRONTclass,isoList,
            isoList = c(as.character(BDRTrem2Iso$Isoform)), selfDf = BDRTrem2Iso, gene = "TREM2", squish=FALSE) 

pdf(paste0(output_dir,"/MainFiguresBDR.pdf"), width = 15, height = 13)
left <- plot_grid(BDRIF + mytheme + theme(legend.position = "top"),
                    plot_grid(plot_grid(Trem2A,TREM2B, ncol = 1, labels = c("D",NULL)),
                              plot_grid(Trem2countsA,Trem2countsB, nrow=1, labels = c("E","F")), nrow = 2, rel_heights = c(0.4,0.6)), ncol = 1, labels = c("A"))
right <- plot_grid(APOE_p, CLU_p, labels = c("B","C"), ncol = 1, rel_heights = c(0.6,0.4))
plot_grid(left, right, nrow = 1)
dev.off()


## ---- Clu Finalised -----  
BDRCluIso <- data.frame(
  Isoform = unlist(BDRCluIso <- list(
    `Mouse` = c("ENSMUST00000022616.13", "PB.14646.139"),
    `Human` = c("ENST00000316403.15", "PB.92671.402")
  )),
  Category = rep(names(BDRCluIso), lengths(BDRCluIso))
)

CluA <- ggTranPlots(inputgtf= gtf$targ_merged, classfiles=class.files$targ_filtered,
            isoList = c(as.character(BDRCluIso$Isoform)), selfDf = BDRCluIso, gene = "Clu", squish=FALSE) 
CluB <- ggTranPlots(inputgtf= gtf$humanMerged, classfiles=BDRONTclass,isoList,
                    isoList = c(as.character(BDRCluIso$Isoform)), selfDf = BDRCluIso, gene = "Clu", squish=FALSE) 
CluTranscriptB <- pBDRClu$trans2 + scale_fill_manual(values = c(label_colour("WT"),label_colour("AD"))) + theme(legend.position = "None")
CluTranscriptA <- plot_transexp_overtime("Clu",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="toprank",rank=1,isoSpecific=c("PB.14646.139"),setorder=c("CONTROL","CASE")) + scale_colour_manual(values = c(label_colour("WT"),label_colour("AD"))) + theme(legend.position = "None") + labs(title = NULL)

pdf("suppFigureBDRClu.pdf", width = 20, height = 10)
plot_grid(CluA,CluTranscriptA,CluB,CluTranscriptB, ncol = 2, labels = c("A","B","C","D"))
dev.off()

## ----- Bin1 Finalised ----- 
BDRBin1Iso <- data.frame(
  Isoform = unlist(BDRBin1Iso <- list(
    `Mouse` = c("ENSMUST00000234496.1","ENSMUST00000025239.8", "PB.22007.224","PB.22007.99"),
    `Human` = c("ENST00000316724.10", "ENST00000409400.1","PB.50706.13","PB.50706.8")
  )),
  Category = rep(names(BDRBin1Iso), lengths(BDRBin1Iso))
)

Bin1A <- ggTranPlots(inputgtf= gtf$targ_merged, classfiles=class.files$targ_filtered,
                     isoList = c(as.character(BDRBin1Iso$Isoform)), selfDf = BDRBin1Iso, gene = "Bin1", squish=TRUE) 
Bin1B <- ggTranPlots(inputgtf= gtf$humanMerged, classfiles=BDRONTclass,isoList,
                    isoList = c(as.character(BDRBin1Iso$Isoform)), selfDf = BDRBin1Iso, gene = "Bin1", squish=TRUE) 
Bin1TranscriptB <- pBDRBin1$trans2 + scale_fill_manual(values = c(label_colour("WT"),label_colour("AD"))) + theme(legend.position = "None")
Bin1TranscriptA <- plot_transexp_overtime("Bin1",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",isoSpecific=c("PB.22007.99"),setorder=c("CONTROL","CASE")) + scale_colour_manual(values = c(label_colour("WT"),label_colour("AD"))) + theme(legend.position = "None") 
Bin1TranscriptC <- pBDRBin1$trans4 + scale_fill_manual(values = c(label_colour("WT"),label_colour("AD"))) + theme(legend.position = "None")
Bin1TranscriptD <- plot_transexp_overtime("Bin1",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",isoSpecific=c("PB.22007.224"),setorder=c("CONTROL","CASE")) + scale_colour_manual(values = c(label_colour("WT"),label_colour("AD"))) + theme(legend.position = "None") 

pdf("suppFigureBDRBin1.pdf", width = 20, height = 10)
plot_grid(Bin1A,Bin1TranscriptD,Bin1TranscriptA,Bin1B,Bin1TranscriptC,Bin1TranscriptB, ncol = 3, labels = c("A","B","C","D","E","F"), 
          rel_widths = c(0.3,0.35,0.35))
dev.off()



                       