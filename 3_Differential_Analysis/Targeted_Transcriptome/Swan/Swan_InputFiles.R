# Szi Kay Leung
# Generate Input files for Swan visualisation 

# SQANTI classification file and gtf 
class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt"
class.files <- read.table(class.names.files, as.is = T, sep = "\t", header = T)

class.names.gtf = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite.gtf"

phenotype <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/TargetedMouse_PhenotypeTAPPAS.txt", header = T) %>% mutate(group_name = ifelse(group == "CONTROL","WT","TG"), name = paste0(sample,"_",group_name,"_",time,"mos"))
phenotype_names <- reshape2::melt(phenotype[,c("sample","name")])

# other folders
swandir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Targeted_Transcriptome/Swan/"
sq_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3"

#### Abundance files for Swan 
swan_counts = class.files %>% select(starts_with("FL."))
# Match the column names to phenotype 
colnames(swan_counts) = phenotype_names$name[match(names(swan_counts),phenotype_names$sample)]
colnames(swan_counts) = substr(colnames(swan_counts),8,14)

# Sum the expression counts based on genotype and age
swan_count_simple = as.data.frame(sapply(unique(colnames(swan_counts)), function(x) rowSums(swan_counts[,grepl(x, colnames(swan_counts))])))
row.names(swan_count_simple) = class.files$isoform
swan_count_simple = swan_count_simple %>% rownames_to_column(., var = "transcript_id")
write.table(swan_count_simple,paste0(swandir,"/swan_count_simple.tsv"), quote = F, row.names = F, sep = "\t")


#### GTF files for each group 
# Subset gtf file
class.gtf <- read.table(class.names.gtf, sep = '\t')
class.gtf$isoform <- word(word(class.gtf $V9,c(4),sep = fixed(" ")),c(1),sep = ";")

# output gtf to sq_dir
subsetgtf_group <- function(var){
  # subset isoforms present by the group (genotype and age) i.e. count > 0
  iso_present = as.data.frame(swan_count_simple) %>% .[.[[var]] > 0, "transcript_id"]
  # subset gtf by the isoforms present 
  class_filter_gtf <- class.gtf  %>% filter(isoform %in% iso_present)
  write.table(class_filter_gtf[c(-10)],paste0(sq_dir,"/", var,"_sqantitamafiltered.classification.gtf"), quote = F, sep = "\t", col.names = F, row.names = F)
}
for(group in colnames(swan_count_simple)){subsetgtf_group(group)}

# on linux
#sq_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3
#cd $sq_dir 
#for gtf in WT_2mos TG_2mos WT_8mos TG_8mos;do 
#  sed 's/transcript_id \([^;]\+\)/transcript_id \"\1\"/g' $gtf"_sqantitamafiltered.classification.gtf" | sed 's/gene_id \([^;]\+\)/gene_id \"\1\"/g' | sed 's/gene_name \([^#;]\+\)/gene_name \"\1\"/g' | sed 's/ref_gene_id \([^;]\+\)/ref_gene_id \"\1\"/g' > $gtf"_sqantitamafiltered.final.classification.gtf"
#done 


#### Novelty Information 
write.csv(class.files[,c("isoform","structural_category")] %>% `colnames<-`(c("tid","novelty")) ,paste0(swandir,"/All_novelty.csv"), quote = F, row.names = T)
