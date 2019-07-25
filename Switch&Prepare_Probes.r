
# Szi Kay Leung
# 25/07/2019: Created rscript for blast_output to switch coordinates if END > START 

# input user argument in bash: Rscript Switch$Prepare_Probes.r <path/blast_output_file> <path/output_directory>
args <- commandArgs(T)
tag <- tools::file_path_sans_ext(args[1])
blast_output_file <- args[1]
output_directory <- args[2]

blast <- read.table(blast_output_file)
# V9 == START coordinates, V10 == END Coordinates
print(paste("Number of probes with CORRECT orientation of probes in START and END:", dim(blast[which(blast$V9 > blast$V10),])[1]))
print(paste("Number of probes with WRONG orientation of probes in START and END:",   dim(blast[which(blast$V9 < blast$V10),])[1]))

# created new column of START and END 
blast$START <- ifelse(blast$V9 < blast$V10, blast$V9, blast$V10) # for START column: if START coordinates < END coordinates, keep START, otherwise switch END
blast$END <- ifelse(blast$V9 < blast$V10, blast$V10, blast$V9) # for END column: if START coordinates < END coordinates, keep END, otherwise switch START

# only retain chromosome, START, END, name of probe
final <- blast[,c("V2","START","END","V1")]
write.table(final,print(paste(output_directory,"final_blast.txt", sep = "/")), row.names = FALSE, quote = FALSE, col.names = FALSE)