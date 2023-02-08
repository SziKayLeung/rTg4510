`%!in%` = Negate(`%in%`)

segregate_tappasresults <- function(output_tappas, type){
  # tappas ouput with p values for each regression coefficient
  dat = output_tappas
  
  # apply Rsquared threshold at 0.5 
  dat = dat %>% filter(`p-value` < 0.05) %>% filter(`R-squared` > 0.5)
  
  if(type == "IsoSeq"){
    ### models
    # 1: casevscontrol yes, time no, timexcase no = genotype effect
    # 2. casevscontrol yes, time yes, timexcase no = genotype+age effect
    # 3. casevscontrol no, time yes, timexcase no = age effect
    # 4. casevscontrol no, time yes, timexcase yes = interaction effect
    # 5. casevscontrol no, time no, timexcase yes = interaction effect
    # 6. casevscontrol yes, time no, timexcase yes = interaction effect
    # 7. casevscontrol yes, time yes, timexcase yes = interaction effect
    ## yes = significant regression coefficient i.e != NA
    
    models = list(
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time) == "TRUE" & is.na(p.valor_TimexCASE) == "TRUE"), # model 1
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE"), # model 2
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "TRUE" & is.na(p.valor_Time) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE"), # model 3
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "TRUE" & is.na(p.valor_Time) == "FALSE" & is.na(p.valor_TimexCASE) == "FALSE"), # model 4
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "TRUE" & is.na(p.valor_Time) == "TRUE" & is.na(p.valor_TimexCASE) == "FALSE"), # model 5
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time) == "TRUE" & is.na(p.valor_TimexCASE) == "FALSE"), # model 6
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time) == "FALSE" & is.na(p.valor_TimexCASE) == "FALSE") # model 7
    )
    names(models) = c("Model 1 Genotype","Model 2 Genotype+Age","Model 3 Age","Model 4 Interaction", "Model 5 Interaction", "Model 6 Interaction","Model 7 Interaction")
    models = lapply(models, function(x) x %>% arrange(-desc(`p-value`)))
    
    
  }else if(type == "RNASeq"){
    # Model 1
    model1 = dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time) == "TRUE" & 
                              is.na(p.valor_TimexCASE) == "TRUE" & is.na(p.valor_Time2xCASE) == "TRUE" & is.na(p.valor_Time3xCASE) == "TRUE") 
    
    # Model 2
    model2a = dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE"
                             & is.na(p.valor_Time2xCASE) == "TRUE" & is.na(p.valor_Time3xCASE) == "TRUE") # model 2
    
    model2b = dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time2) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE"
                             & is.na(p.valor_Time2xCASE) == "TRUE" & is.na(p.valor_Time3xCASE) == "TRUE") # model 2
    
    model2c = dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time3) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE"
                             & is.na(p.valor_Time2xCASE) == "TRUE" & is.na(p.valor_Time3xCASE) == "TRUE") # model 2
    
    model2 = bind_rows(list(model2a,model2b,model2c)) %>% distinct()
    
    # Model 3
    model3a =  dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "TRUE" & is.na(p.valor_Time) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE", 
                              is.na(p.valor_Time2xCASE) == "TRUE", is.na(p.valor_Time3xCASE) == "TRUE")
    model3b =  dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "TRUE" & is.na(p.valor_Time2) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE",
                              is.na(p.valor_Time2xCASE) == "TRUE", is.na(p.valor_Time3xCASE) == "TRUE")
    model3c =  dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "TRUE" & is.na(p.valor_Time3) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE",
                              is.na(p.valor_Time2xCASE) == "TRUE", is.na(p.valor_Time3xCASE) == "TRUE")
    model3 = bind_rows(list(model3a,model3b,model3c)) %>% distinct()
    
    non_interaction = bind_rows(model1,model2,model3)
    interaction = dat[dat$...1 %!in% non_interaction$...1,]
    
    models = list(model1,model2,model3,interaction)
    names(models) = c("Model 1 Genotype","Model 2 Genotype+Age","Model 3 Age","Model 4 - 7 Interaction")
    models = lapply(models, function(x) x %>% arrange(-desc(`p-value`)))
    
  }else{
    print("IsoSeq or RNASeq required as argument")
  }
  
  p <- sapply(models, nrow) %>% reshape2::melt() %>% rownames_to_column(var = "model") %>% mutate(num = c(1:nrow(.))) %>%
    ggplot(., aes(x = reorder(model, -num), y = value)) + geom_bar(stat = "identity") +
    labs(y = "Number of Differentially Expressed Genes",x = "") + mytheme + coord_flip() 
  
  df = sapply(models, nrow) %>% reshape2::melt() %>% as.data.frame() %>% rownames_to_column(var = "Model") %>% mutate(value = as.numeric(as.character(value)))
  print(df)
  cat("Total of Number differentially expressed genes:", sum(df$value),"\n")
  cat("Number differentially expressed under Interaction effect:", sum(df[grepl("Interaction",df$Model),"value"]),"\n")
  cat("Number differentially expressed under Genotype effect:", sum(df[df$Model == "Model 1 Genotype","value"]),"\n")
  
  output = list(models,p)
  names(output) = c("models","p")
  return(output)
}