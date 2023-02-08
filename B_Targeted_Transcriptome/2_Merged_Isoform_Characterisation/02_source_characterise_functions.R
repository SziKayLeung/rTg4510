# plot label colour
label_colour <- function(genotype){
  if(genotype %in% c("WT","CONTROL")){colour = wes_palette("Royal1")[2]}else{
    if(genotype == "WT_2mos"){colour = alpha(wes_palette("Royal1")[2],0.5)}else{
      if(genotype %in% c("TG","CASE")){colour = wes_palette("Royal1")[1]}else{
        if(genotype == "TG_2mos"){colour = alpha(wes_palette("Royal1")[1],0.5)}else{
          if(genotype == "mouse"){colour = wes_palette("Royal1")[4]}else{
            if(genotype == "novel"){colour = wes_palette("Darjeeling1")[4]}else{
              if(genotype == "known"){colour = wes_palette("Darjeeling1")[5]}else{
              }}}}}}}
  return(colour)
}


plot_target_diff_col <- function(gene, type){
  if(type == "IsoSeq"){
    if(gene == "Trem2"){
      cset = c("PB.3742.1"=wes_palette("Darjeeling2")[3], "PB.3742.2"=wes_palette("Darjeeling2")[2],"PB.3742.3"=wes_palette("Darjeeling1")[1])
    }else if(gene == "Rhbdf2"){
      cset = c("PB.1691.13"=wes_palette("Darjeeling2")[1], "PB.1691.4"=wes_palette("Darjeeling2")[2])
    }else if(gene == "Ptk2b"){
      cset = c("PB.2637.336"=wes_palette("Darjeeling1")[1], "PB.2637.44"=wes_palette("Darjeeling2")[2],"PB.2637.334"=wes_palette("Darjeeling2")[3],
               "PB.2637.44"=wes_palette("Cavalcanti1")[2], "PB.2637.31"=wes_palette("IsleofDogs1")[2],
               "PB.2637.161"=wes_palette("Darjeeling1")[5], "PB.2637.159"=wes_palette("GrandBudapest2")[1]) 
    }else if(gene == "Apoe"){
      cset = c("PB.7333.32"=wes_palette("Darjeeling1")[1], "PB.7333.32"=wes_palette("Darjeeling2")[5],"PB.7333.39"=wes_palette("IsleofDogs1")[2],
               "PB.7333.25"=wes_palette("Cavalcanti1")[2], "PB.7333.85"=wes_palette("Darjeeling1")[2])
    }else if(gene == "Clu"){
      cset = c("PB.2634.2"=wes_palette("Darjeeling1")[1], "PB.2634.419"=wes_palette("Darjeeling2")[2],
               "PB.2634.256"=wes_palette("Darjeeling2")[3]) 
    }else if(gene == "App"){
      cset = c("PB.3388.103"=wes_palette("Darjeeling1")[1], "PB.3388.81"=wes_palette("Darjeeling2")[2],
               "PB.3388.1143"=wes_palette("Darjeeling2")[3],
               "PB.3388.2081"=wes_palette("Cavalcanti1")[2],"PB.3388.114"= wes_palette("Moonrise3")[2],
               "PB.3388.1628"=wes_palette("Rushmore1")[4],"PB.3388.1398"=wes_palette("Royal1")[3],
               "PB.3388.1653"=wes_palette("Darjeeling2")[5],"PB.3388.110"=wes_palette("Moonrise3")[3]) 
      
    }else if(gene == "Bin1"){
      cset = c("PB.3915.33"=wes_palette("Darjeeling1")[1], "PB.3915.5"=wes_palette("Darjeeling2")[2],
               "PB.3915.455"=wes_palette("Darjeeling2")[3], "PB.3915.1"=wes_palette("Cavalcanti1")[2], 
               "PB.3915.6"=wes_palette("IsleofDogs1")[2],
               "PB.3915.2"=wes_palette("Darjeeling1")[5]) 
    }else{
      cset = "NULL"
    }
    
  }else{
    if(gene == "Trem2"){
      cset = c("TALONT000740495"=wes_palette("IsleofDogs1")[1],
               "ENSMUST00000113237.3"=wes_palette("Cavalcanti1")[2],
               "ENSMUST00000024791.14"=wes_palette("Darjeeling1")[1])
    }else if(gene == "Abca1"){
      cset = c("TALONT000975469"=wes_palette("IsleofDogs1")[1], "TALONT000972642"=wes_palette("Cavalcanti1")[2],"ENSMUST00000030010.3"=wes_palette("Darjeeling1")[1])
    }else if(gene == "Abca7"){
      cset = c("TALONT000222205"=wes_palette("Darjeeling2")[3], "TALONT000226765"=wes_palette("Darjeeling2")[2],
               "TALONT000226746"=wes_palette("IsleofDogs1")[1], "TALONT000227148"=wes_palette("Cavalcanti1")[2],"TALONT000222245"=wes_palette("Darjeeling1")[1])
    }else if(gene == "Rhbdf2"){
      cset = c("TALONT000349463"=wes_palette("Cavalcanti1")[2], "ENSMUST00000103029.9"=wes_palette("Darjeeling2")[1])
    }else if(gene == "Rhbdf2"){
      cset = c("TALONT000349463"=wes_palette("Cavalcanti1")[2], "ENSMUST00000103029.9"=wes_palette("Darjeeling2")[1])
    }else if(gene == "Fyn"){
      cset = c("TALONT000199007"=wes_palette("Darjeeling1")[1], "ENSMUST00000099967.9"=wes_palette("Darjeeling2")[1], 
               "TALONT000196688"=wes_palette("IsleofDogs1")[1], "TALONT000199061"=wes_palette("Darjeeling2")[2])
    }else if(gene == "Sorl1"){
      cset = c("TALONT001406989"=wes_palette("Darjeeling2")[3], "TALONT001373799"=wes_palette("Darjeeling2")[2],"TALONT001403727"=wes_palette("Darjeeling1")[1],
               "ENSMUST00000060989.8"=wes_palette("IsleofDogs1")[1], "TALONT001373549"=wes_palette("Cavalcanti1")[2], "TALONT001403708"=wes_palette("IsleofDogs1")[2])
    }else if(gene == "Cd33"){
      cset = c("TALONT001237522"=wes_palette("Darjeeling2")[3], "TALONT001237572"=wes_palette("Darjeeling2")[2],"ENSMUST00000205503.1"=wes_palette("Darjeeling1")[1],
               "TALONT001237520"=wes_palette("IsleofDogs1")[1], "TALONT001237550"=wes_palette("Cavalcanti1")[2], "TALONT001237524"=wes_palette("IsleofDogs1")[2],
               "TALONT001237514"=wes_palette("Darjeeling1")[4], "TALONT001237510"=wes_palette("Darjeeling2")[4], "TALONT001237528"=wes_palette("Darjeeling2")[5])
    }else if(gene == "Fus"){
      cset = c("TALONT001283370"=wes_palette("Darjeeling2")[3], "TALONT001283437"=wes_palette("Darjeeling2")[2],"TALONT001283784"=wes_palette("Darjeeling1")[1],
               "ENSMUST00000174196.7"=wes_palette("Cavalcanti1")[2], "ENSMUST00000077609.11"=wes_palette("IsleofDogs1")[2],
               "ENSMUST00000106251.9"=wes_palette("Darjeeling1")[5])
    }else if(gene == "Snca"){
      cset = c("TALONT001103657"=wes_palette("Darjeeling2")[3], "TALONT001103777"=wes_palette("Darjeeling2")[2],"TALONT001104529"=wes_palette("Darjeeling1")[1],
               "ENSMUST00000163779.7"=wes_palette("Cavalcanti1")[2], "ENSMUST00000114268.4"=wes_palette("Darjeeling1")[5])
    }else if(gene == "Apoe"){
      cset = c("TALONT001164086"=wes_palette("Darjeeling2")[3], "ENSMUST00000173739.7"=wes_palette("Darjeeling2")[2],"TALONT001166657"=wes_palette("GrandBudapest2")[2],
               "TALONT001163706"=wes_palette("GrandBudapest2")[1], "ENSMUST00000174064.8"=wes_palette("Darjeeling1")[5])
    }else if(gene == "Ptk2b"){
      cset = c("ENSMUST00000089250.8"=wes_palette("Darjeeling1")[1], "TALONT000490475"=wes_palette("Darjeeling2")[2],
               "ENSMUST00000022622.13"=wes_palette("Darjeeling2")[3],
               "TALONT000490475"=wes_palette("Cavalcanti1")[2], "ENSMUST00000136216.7"=wes_palette("IsleofDogs1")[2],
               "TALONT000490424"=wes_palette("Darjeeling1")[5], "TALONT000492060"=wes_palette("GrandBudapest2")[2],
               "TALONT000492849"=wes_palette("GrandBudapest2")[4]) 
    }else if(gene == "Bin1"){
      cset = c("ENSMUST00000234496.1"=wes_palette("Darjeeling1")[1], "ENSMUST00000234857.1"=wes_palette("Darjeeling2")[2],
               "TALONT000785854"=wes_palette("Darjeeling2")[3],
               "ENSMUST00000025239.8"=wes_palette("Cavalcanti1")[2], 
               "TALONT000761829"=wes_palette("Darjeeling1")[5]) 
      
    }else if(gene == "Picalm"){
      cset = c("ENSMUST00000207225.1"=wes_palette("Darjeeling1")[1], "ENSMUST00000049537.8"=wes_palette("Darjeeling2")[2],
               "ENSMUST00000207225.1"=wes_palette("Darjeeling2")[3],
               "TALONT001255854"=wes_palette("Cavalcanti1")[2], 
               "TALONT001263109"=wes_palette("Darjeeling1")[5],"TALONT001263261"=wes_palette("GrandBudapest2")[4],
               "TALONT001260744"=wes_palette("GrandBudapest2")[2],"TALONT001260731"=wes_palette("IsleofDogs1")[2]) 
      
    }else if(gene == "Tardbp"){
      cset = c("ENSMUST00000045180.13"=wes_palette("Darjeeling1")[1], "ENSMUST00000147391.2"=wes_palette("Darjeeling2")[2],
               "ENSMUST00000084125.9"=wes_palette("Darjeeling2")[3],
               "TALONT00100759"=wes_palette("Cavalcanti1")[2], 
               "TALONT001007610"=wes_palette("Darjeeling1")[5],"TALONT001008280"=wes_palette("GrandBudapest2")[4],
               "TALONT001008011"=wes_palette("GrandBudapest2")[2]) 
      
      
    }else if(gene == "Clu"){
      cset = c("ENSMUST00000022616.13"=wes_palette("Darjeeling1")[1], "TALONT000483931"=wes_palette("Darjeeling2")[2],
               "TALONT000465283"=wes_palette("Darjeeling2")[3],
               "TALONT000489417"=wes_palette("Cavalcanti1")[2], 
               "TALONT000485940"=wes_palette("Darjeeling1")[5],"TALONT000440128"=wes_palette("GrandBudapest2")[4],
               "TALONT000440029"=wes_palette("IsleofDogs1")[2],"TALONT000440314"=wes_palette("Moonrise3")[2]) 
      
      
    }else if(gene == "App"){
      cset = c("ENSMUST00000227654.1"=wes_palette("Darjeeling1")[1], "ENSMUST00000005406.11"=wes_palette("Darjeeling2")[2],
               "TALONT000637770"=wes_palette("Darjeeling2")[3],
               "TALONT000706195"=wes_palette("Cavalcanti1")[2], 
               "TALONT000631124"=wes_palette("Darjeeling1")[5],"TALONT000632929"=wes_palette("GrandBudapest2")[4],
               "TALONT000633505"=wes_palette("IsleofDogs1")[2]) 
      
    }else if(gene == "Mapt"){
      cset = c("ENSMUST00000126820.1"=wes_palette("Darjeeling1")[1], "ENSMUST00000106992.9"=wes_palette("Darjeeling2")[2],
               "TALONT000338410"=wes_palette("Darjeeling2")[3],
               "TALONT000334577"=wes_palette("Cavalcanti1")[2], 
               "TALONT000333842"=wes_palette("Darjeeling1")[5],"TALONT000312346"=wes_palette("Moonrise3")[2]) 
      
    }else if(gene == "Vgf"){
      cset = c("ENSMUST00000041543.8"=wes_palette("Darjeeling1")[1], "TALONT001069675"=wes_palette("Darjeeling2")[2],
               "TALONT001069700"=wes_palette("Darjeeling2")[3],
               "TALONT00106972"=wes_palette("Cavalcanti1")[2], 
               "TALONT001069768"=wes_palette("Darjeeling1")[5],"TALONT001070640"=wes_palette("Moonrise3")[2],
               "TALONT001069788"=wes_palette("IsleofDogs1")[2],"TALONT001069765"=wes_palette("GrandBudapest2")[4]) 
      
    }else{
      cset = "NULL"
    }
  }
  return(cset)
}