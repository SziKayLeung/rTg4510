
library("ggplot2")
LOGEN_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
SC_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510/Paper_Figures/"
source(paste0(SC_ROOT, "0_source_functions.R"))
source(paste0(SC_ROOT, "rTg4510_config.R"))
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Mouse_Isoseq/"

BDRONTclass <- read.table("/lustre/recovered/Research_Project-MRC148213/sl693/AD_BDR/D_ONT/5_cupcake/7_sqanti3/ontBDR_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt", sep = "\t", as.is = T)
BDRONT <- readRDS("/lustre/recovered/Research_Project-MRC148213/sl693/AD_BDR/01_figures_tables/Ont_DESeq2TranscriptLevel.RDS")
phenotype <- 
plotIFAll(Exp=Exp$targ_ont$normAll %>% select(-associated_gene),
          classf=class.files$targ_all,
          pheno=phenotype$targeted_rTg4510_ont,
          majorIso=row.names(TargetedDIU$ontDIUGeno$keptIso))

TargetGene <- toupper(c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7",
                        "Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf"))

tabulateIF <- function(gene, classf, countcol, genespecific=NULL){
  
  Counts <- classf %>% select(isoform,contains(countcol))
  rownames(Counts) <- Counts$isoform
  Counts <- Counts %>% select(-isoform)
  
  
  # Calculate the mean of normalised expression across all the samples per isoform
  meandf <- data.frame(meanvalues = apply(Counts,1,mean)) %>%
    rownames_to_column("isoform") %>% 
    # annotate isoforms with associated_gene and structural category
    left_join(., classf[,c("isoform","associated_gene","structural_category")], by = "isoform") 
  

  # Group meandf by associated_gene and calculate the sum of mean values for each group
  grouped <- aggregate(meandf$meanvalues, by=list(associated_gene=meandf$associated_gene), FUN=sum)
  
  # Calculate the proportion by merging back, and divide the meanvalues by the grouped values (x)
  IF <- meandf %>% 
    left_join(grouped, by = "associated_gene") %>%
    mutate(perc = meanvalues / x * 100) 
    
  return(IF)
  
}

tabIFOut <- lapply(TargetGene, function(x) tabulateIF(x, BDRONTclass, "B2"))
names(tabIFOut) <- TargetGene
tabIFOut  <- bind_rows(tabIFOut, .id = "associated_gene") 
ggplot(tabIFOut, aes(x = associated_gene, y = ""))

ggplot(merged, aes(x = associated_gene, y = as.numeric(perc), fill = forcats::fct_rev(structural_category))) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +
  #scale_color_manual(values = rep(NA, length(unique(minorgrouped$gene)))) + 
  labs(x = "Gene", y = "Isoform fraction (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name = "Isoform Classification", values = rev(c(alpha("#00BFC4",0.8),alpha("#00BFC4",0.3),
                                                                    alpha("#F8766D",0.8),alpha("#F8766D",0.3)))) +
  theme(legend.position = "None")


View(TREM2)

# human BDR dataset 
# ENST00000373113.8
phenotype <- read.csv("/lustre/recovered/Research_Project-MRC148213/sl693/AD_BDR/0_metadata/B_ONT/Selected_ONTTargeted_BDR.csv", header = T)
Trem2TranscriptBDR <- BDRONT$B2WaldBraak$norm_counts %>% filter(isoform == "PB.81888.25") %>% filter(grepl("B2", sample)) %>% 
  mutate(sample = str_remove(sample, "B2.")) %>% 
  left_join(., phenotype, by = "sample") %>%
  filter(BraakTangle_numeric %in% c(0,1,2,5,6)) %>% 
  mutate(phenotype = ifelse(BraakTangle_numeric %in% c(0,1,2),"Control","AD")) %>%
  mutate(phenotype = factor(phenotype, levels = c("Control","AD"))) 

t.test(normalised_counts ~ phenotype, data = Trem2TranscriptBDR %>% filter(sample != "BBN00229416.1"))

summary(lm(normalised_counts ~ phenotype, data = Trem2TranscriptBDR %>% filter(sample != "BBN00229416.1")))
ggplot(Trem2TranscriptBDR, aes(x = phenotype, y = normalised_counts)) + geom_boxplot() 

bdrGtf <- as.data.frame(rtracklayer::import("/lustre/recovered/Research_Project-MRC148213/sl693/AD_BDR/D_ONT/5_cupcake/7_sqanti3/ontBDR_collapsed.filtered_counts_filtered.gtf"))
refGtf <- as.data.frame(rtracklayer::import("/lustre/projects/Research_Project-MRC148213/lsl693/references/human/CLU.gencode.v40.annotation.gtf"))
gtf$humanMerged <- rbind(bdrGtf [,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                         refGtf [,c("seqnames","strand","start","end","type","transcript_id","gene_id")])

ggTranPlots(gtf$humanMerged,BDRONTclass,isoList = c("PB.92671.402","PB.92671.1702","ENST00000316403.15"),simple=TRUE, colour = c("blue","red", "red"))

CluTranscriptBDR <- BDRONT$B2WaldBraak$norm_counts %>% filter(isoform == "PB.92671.402") %>% filter(grepl("B2", sample)) %>% 
  mutate(sample = str_remove(sample, "B2.")) %>% 
  left_join(., phenotype, by = "sample") %>%
  filter(BraakTangle_numeric %in% c(0,1,2,5,6)) %>% 
  mutate(phenotype = ifelse(BraakTangle_numeric %in% c(0,1,2),"Control","AD")) %>%
  mutate(phenotype = factor(phenotype, levels = c("Control","AD"))) 
ggplot(CluTranscriptBDR, aes(x = phenotype, y = normalised_counts)) + geom_boxplot() 

BDRONT$B2WaldBraak$norm_counts %>% filter(isoform == "PB.92671.402") %>% filter(grepl("B2", sample)) %>% 
  mutate(sample = str_remove(sample, "B2.")) %>% 
  left_join(., phenotype, by = "sample") %>%
  mutate(BraakTangle_numeric = as.factor(BraakTangle_numeric)) %>%
  ggplot(., aes(x = BraakTangle_numeric, y = normalised_counts)) + geom_boxplot() 

CluNE <- paste0("PB.92671.",c("689","1693","3228","4691","3474","1988","1778","1875","1860","3549","4725"))
ggTranPlots(gtf$humanMerged,BDRONTclass,isoList = c("PB.92671.402",CluNE ),simple=TRUE, colour = c("blue",rep("red",12)))

CluAF <- paste0("PB.92671.",c("2479","26717"))

BDRONT$B2WaldBraak$norm_counts %>% filter(isoform == "PB.92671.26717") %>% filter(grepl("B2", sample)) %>% 
  mutate(sample = str_remove(sample, "B2.")) %>% 
  left_join(., phenotype, by = "sample") %>%
  mutate(BraakTangle_numeric = as.factor(BraakTangle_numeric)) %>%
  ggplot(., aes(x = BraakTangle_numeric, y = normalised_counts)) + geom_boxplot() 

BDRONTclass %>% select(isoform, contains("B2")) %>% 
  map_if(is.numeric, ~./sum(.) * 100) %>%
  as_data_frame() %>% 
  filter(isoform == "PB.92671.402") %>% 
  reshape2::melt(variable.name = "sample", value.name = "perc") %>% 
  mutate(sample = str_remove(sample, "B2.")) %>% 
  left_join(., phenotype[,c("sample","BraakTangle_numeric")], by = "sample") %>% 
  mutate(phenotype = ifelse(BraakTangle_numeric %in% c(0,1,2),"Control","AD")) %>%
  ggplot(., aes(x = as.factor(BraakTangle_numeric), y = perc)) + geom_boxplot()


## Bin1
ggTranPlots(gtf$humanMerged,BDRONTclass,isoList = c("PB.50706.8", "PB.50706.13"),simple=TRUE, colour = c("blue",rep("red",12)))
BDRONT$B2WaldBraak$norm_counts %>% filter(isoform == "PB.50706.8") %>% filter(grepl("B2", sample)) %>% 
  mutate(sample = str_remove(sample, "B2.")) %>% 
  left_join(., phenotype, by = "sample") %>%
  mutate(BraakTangle_numeric = as.factor(BraakTangle_numeric)) %>%
  ggplot(., aes(x = BraakTangle_numeric, y = normalised_counts)) + geom_boxplot() 

BDRONT$B2WaldBraak$norm_counts %>% filter(isoform == "PB.50706.8") %>% filter(grepl("B2", sample)) %>% 
  mutate(sample = str_remove(sample, "B2.")) %>% 
  left_join(., phenotype, by = "sample") %>%
  filter(BraakTangle_numeric %in% c(0,1,2,5,6)) %>% 
  mutate(BraakTangle_numeric = as.factor(BraakTangle_numeric)) %>%
  mutate(phenotype = ifelse(BraakTangle_numeric %in% c(0,1,2),"Control","AD")) %>%
  ggplot(., aes(x = phenotype, y = normalised_counts)) + geom_boxplot() 


BDRONTclass %>% select(isoform, contains("B2")) %>% 
  map_if(is.numeric, ~./sum(.) * 100) %>%
  as_data_frame() %>% 
  filter(isoform == "PB.50706.8") %>% 
  reshape2::melt(variable.name = "sample", value.name = "perc") %>% 
  mutate(sample = str_remove(sample, "B2.")) %>% 
  left_join(., phenotype[,c("sample","BraakTangle_numeric")], by = "sample") %>% 
  mutate(phenotype = ifelse(BraakTangle_numeric %in% c(0,1,2),"Control","AD")) %>%
  ggplot(., aes(x = as.factor(BraakTangle_numeric), y = perc)) + geom_boxplot()

dat <- BDRONTclass %>% select(isoform, contains("B2")) %>% 
  map_if(is.numeric, ~./sum(.) * 100) %>%
  as_data_frame() %>% 
  filter(isoform == "PB.50706.8") %>% 
  reshape2::melt(variable.name = "sample", value.name = "perc") %>% 
  mutate(sample = str_remove(sample, "B2.")) %>% 
  left_join(., phenotype, by = "sample") %>% 
  filter(BraakTangle_numeric %in% c(0,1,2,5,6)) %>%
  mutate(phenotype = ifelse(BraakTangle_numeric %in% c(0,1,2),"Control","AD")) 
ggplot(dat, aes(x =  phenotype, y = perc)) + geom_boxplot()
dat <- dat %>% filter(perc < 2)
summary(lm(perc ~ phenotype,  data = dat))

model = lm(perc ~ phenotype, data = dat, family = "binomial")
plot(model)

BDRONTclass[BDRONTclass$isoform == "PB.50706.8",]
BDRONTclass[BDRONTclass$isoform == "PB.50706.13",]
