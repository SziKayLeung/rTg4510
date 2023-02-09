## ---------- Script -----------------
##
## Purpose: input variables for differential analysis of Iso-Seq targeted mouse transcriptome datasets 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)


## ---------- packages -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))
suppressMessages(library("readxl"))
suppressMessages(library("cowplot"))


## ---------- Source function and config files -----------------

# source all general scripts related to long-read sequencing
LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"

source(paste0(LOGEN, "aesthetics_basics_plots/pthemes.R"))
sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"differential_analysis"), pattern="*.R", full = T), source,.GlobalEnv)


## ---------- Load tappAS files -----------------
loaded <- list(
  iso = input_tappasfiles(tappas_dirnames$iso),
  ont = input_tappasfiles(tappas_dirnames$ont)
)


## ---------- Annotate tappAS files -----------------
annotated <- list(
  iso = annotate_tappasfiles(input.class.files$iso,loaded$iso$input_norm,misc_input$iso),
  ont = annotate_tappasfiles(input.class.files$ont_unfil,loaded$ont$input_norm,misc_input$ont)
)


##### Differential analysis ############################# 
# Differential Gene 
pGeneExp <- list(iso = generate_diff_plots(TargetGene,annotated$iso,"IsoGene_Targeted","Iso-Seq Expression"),
                 ont = generate_diff_plots(TargetGene,annotated$ont,"ONTGene_Targeted","ONT Expression"))

# Differential Transcript 
pTransExp <- list(iso = generate_diff_plots(TargetGene,annotated$iso,"Transcript","Iso-Seq Expression"),
                  ont = generate_diff_plots(TargetGene,annotated$ont,"Transcript","ONT Expression"))

pTransExpTraj <- list(iso = generate_diff_plots(TargetGene,IsoExp,"Transcript Trajectory","Iso-Seq Expression",tappassigtrans$iso$TargetedIso_Transexp),
                      ont = generate_diff_plots(TargetGene,OntExp,"Transcript Trajectory","ONT Expression", tappassigtrans$ont$TargetedOnt_Transexp))


pdf(paste0(output_dir,"/Top30_ONT_DEA.pdf"), width = 10, height = 8)
for(i in 1:30){
  iso = tappassigtrans$ont$TargetedOnt_Transexp$isoform[i]
  print(plot_grid(plot_trans_exp_individual(iso, annotated$ont$Norm_transcounts),
            plot_trans_exp_individual_overtime(iso, annotated$ont$Norm_transcounts)))
}
dev.off()

dat <- data.frame()
for(i in 1:30){
  iso = tappassigtrans$ont$TargetedOnt_Transexp$isoform[i]
  dat[i,1] <- iso 
  
  dat[i,2] <- if(nrow(misc_input$ont_sq_reason %>% filter(filtered_isoform == iso)) == 0){
    "NA"
  }else{
    as.character(misc_input$ont_sq_reason[misc_input$ont_sq_reason$filtered_isoform == iso, "reason"])
  }
  
}

colnames(dat) <- c("isoform","filtered_reason")
dat
