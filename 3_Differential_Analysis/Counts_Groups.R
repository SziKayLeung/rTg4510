#### Swan visualisations
swan_counts = class.files %>% select(starts_with("FL."))
colnames(swan_counts) = substr(colnames(swan_counts),8,14)

# Sum the expression counts based on genotype and age
swan_count_simple = as.data.frame(sapply(unique(colnames(swan_counts)), function(x) rowSums(swan_counts[,grepl(x, colnames(swan_counts))])))
swan_count_simple = swan_count_simple %>% rownames_to_column(., var = "transcript_id")
write.table(swan_count_simple,paste0(output_helpfig_dir,"/swan_count_simple.tsv"), quote = F, row.names = F, sep = "\t")

# generate separate files for PBIDs present in group 
group_sqanti_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/GROUP/"
for(group in group_order){
  dat = swan_count_simple[swan_count_simple[[group]],"transcript_id"] %>% as.data.frame()
  write.table(dat,paste0(group_sqanti_dir,group,"_counts.txt"), row.names = F, quote = F, col.names = F)
}

# generate separate files for PBIDs present in each sample
sample_sqanti_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/SAMPLE/"
abundance =  class.files %>% select(starts_with("FL."))
for(col in colnames(abundance)){
  file = word(col,c(2), sep = fixed("."))
  cat("Processing",file,"\n")
  dat = rownames(abundance[abundance[[col]] != "0",]) %>% as.data.frame() 
  write.table(dat,paste0(sample_sqanti_dir,file,"_counts.txt"), row.names = F, quote = F, col.names = F)
  #nrow(abundance[abundance$FL.S23_WT_8mos != "0",])
  
  dat2 = abundance %>% select(col) %>% rownames_to_column(var = "transcript_id")
  write.table(dat2,paste0(sample_sqanti_dir,file,"_swan_counts.txt"), row.names = F, quote = F, col.names = T )
}

# Subset gtf file
class.gtf <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/WholeIsoSeq_sqantitamafiltered.final.classification.gtf", sep = '\t')
class.gtf$isoform <- word(word(class.gtf $V9,c(4),sep = fixed(" ")),c(1),sep = ";")

# output gtf to sq_dir
sq_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME"
subsetgtf_group <- function(var){
  # subset isoforms present by the group (genotype and age) i.e. count > 0
  iso_present = as.data.frame(swan_count_simple) %>% rownames_to_column(., var = "isoform") %>% .[.[[var]] > 0, "isoform"]
  # subset gtf by the isoforms present 
  class_filter_gtf <- class.gtf  %>% filter(isoform %in% iso_present)
  write.table(class_filter_gtf[c(-10)],paste0(sq_dir,"/", var,"_sqantitamafiltered.classification.gtf"), quote = F, sep = "\t", col.names = F, row.names = F)
}
for(group in colnames(swan_count_simple)){subsetgtf_group(group)}
# on linux
#sq_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME
#cd $sq_dir 
#for gtf in WT_2mos TG_2mos WT_8mos TG_8mos;do 
#  sed 's/transcript_id \([^;]\+\)/transcript_id \"\1\"/g' $gtf"_sqantitamafiltered.classification.gtf" | sed 's/gene_id \([^;]\+\)/gene_id \"\1\"/g' | sed 's/gene_name \([^#;]\+\)/gene_name \"\1\"/g' | sed 's/ref_gene_id \([^;]\+\)/ref_gene_id \"\1\"/g' > $gtf"_sqantitamafiltered.final.classification.gtf"
#done 

# novelty information 
swandir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Vis/InputFiles"
swangenes = c("Cd34","Cd68","Osmr","C4b","Ubqln1","Gjb2","Adam23")
for(gene in swangenes){
  dat <- class.files[class.files$associated_gene == gene, c("isoform","structural_category")]
  colnames(dat) <- c("tid","novelty")
  write.csv(dat,paste0(swandir,"/",gene,"_novelty.csv"), quote = F, row.names = F)
}
write.csv(class.files[,c("isoform","structural_category")] %>% `colnames<-`(c("tid","novelty")) ,paste0(swandir,"/All_novelty.csv"), quote = F, row.names = T)

