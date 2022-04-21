#########################################
## Report Changes Regulatory Landscape
#########################################
##for now only adding the images
report_regulatoryLandscapeChanges<-function(reportUnit, patientResults){
  
  ###I need to know whether the considered gene is selected by means of a 
  ## LOF && GOF mechanism || LONG-RANGE--DIRECT EFFECT
  
  targetGene<-unlist(strsplit(x = reportUnit, split = "_", fixed = TRUE))[1]
  targetMech<-unlist(strsplit(x = reportUnit, split = "_", fixed = TRUE))[2]
  targetPhase<-unlist(strsplit(x = reportUnit, split = "_", fixed = TRUE))[3]
  
  #This way of approaching the GOF-LOF info, deprecated
  # pathoMechanism<-patientResults$resultsPerPhase[[targetPhase]][[targetPhase]][targetGene,]$Mechanism
  # pathoMechanism<-unlist(strsplit(x = pathoMechanism,
  #                                 split = ":",
  #                                 fixed = TRUE))[2]##GOF or LOF
  pathoMechanism<-targetMech
  
  ##To discriminate in dupl deletions, do something like, pathomechanism_Impact by long Range, or by direct effect
  ##As a gene duplication can trigger the disease due to neoTad and hence, dy a long-range mechanism
  #Deletion, truncation...important for the short explanation if not long-range effect
  
  targetRow_resultsPhase<-paste(targetGene,"--",targetMech, sep="")##In case interest on parsing resultsPhase
  pathomechanism_ImpactOverGene<-patientResults$resultsPerPhase[[targetPhase]][[targetPhase]][targetRow_resultsPhase,]$GeneImpact

  
  sv_type<-patientResults$patientInfo$TypeSV

  ##OPENING Section
  regLanChanges<-"<div class='regulatoryLandscapeChangesSection'>"
  
  ###################
  ## Text variables to be repeated
  text_gof1<-"To be a strong candidate for 'Gain of Function' through a long-range pathological mechanism, the gene should gain enhancer activity. To estimate enhancer activity we have considered the overall H3K27ac levels associated to the enhancers."
  
  text_lof1<-"To be a strong candidate for 'Loss of Function' through a long-range pathological mechanism, the gene should lose enhancer activity. To estimate enhancer activity we have considered the overall H3K27ac levels associated to the enhancers."
  
  
  ######################################################
  ##Adding images. 
  ##If long-range pathological mechanism
  ######################################################
  
  if((pathomechanism_ImpactOverGene ==  "LongRange") || 
     (pathomechanism_ImpactOverGene ==  "LongRange_geneDuplication") || 
     (pathomechanism_ImpactOverGene ==  "Direct_LongRange_geneDuplication")){
    
    ##LongRange: gene intact, 
    ##LongRange_geneDuplication means gene duplicated but pathological mechanism predicted by LongRange effects
    ##Direct_LongRange_geneDuplication, gene expressed duplicated & gaining a lot of enh difficult to predict main cause
    
    ##For fusion transcripts...I think we will use the state of the active promoter and that is all
    ##And not enh
    regLanChanges<-paste(regLanChanges,
                         "<div class='barplotEnh'>",
                         paste0("<img class='barplotEnhImage' src='graphicalSummaries/", reportUnit,"_", patientResults$job_UniCode, "_regulatoryLandscapeBarplotChanges.png", "'",
                                ">"),
                         "</div>",
                         sep="") 
  }else{
    #So direct impact over gene
    # Check here if need to improve this
    regLanChanges<-paste(regLanChanges,
                         "<div class='directImpactPlot'>",
                         #Quitando el image Id
                         paste0("<img class='directplotEnhImage' src='graphicalSummaries/", reportUnit,"_", patientResults$job_UniCode, "_regulatoryLandscapeBarplotChanges.png", "'",
                                ">"),
                         "</div>",
                         sep="")
    
  }
  
  ######################
  ## Adding Explanation
  ######################
  ##We need to know again the barplot information

  infoBarplot<-patientResults$resultsPerPhase[[targetPhase]][[targetPhase]][targetRow_resultsPhase,]
  
  ##Computing INITIAL numbers
  nEnh_Initial<-infoBarplot$avg_nEnh_InitialLeft + 
    infoBarplot$avg_nEnh_InitialRight
  
  acetilation_Initial<-infoBarplot$avg_AcetilationEnh_InitialLeft +
    infoBarplot$avg_AcetilationEnh_InitialRight
  
  ##Computing FINAL numbers
  if(pathoMechanism == "LOF"){
    ##For LOF we only consider cognate enh, since doing the balance gained-lost was losing some enhanceropathies related with SATB2 or SLC2A1
    nEnh_Final<-infoBarplot$avg_nEnh_KeptLeft + 
      infoBarplot$avg_nEnh_KeptRight 
    
    acetilation_Final<-infoBarplot$avg_AcetilationEnh_KeptLeft +
      infoBarplot$avg_AcetilationEnh_KeptRight   
    
  }else if(pathoMechanism == "GOF"){
    
    ##For GOF we consider what is kept and of course what is gained
    nEnh_Final<-infoBarplot$avg_nEnh_KeptLeft + 
      infoBarplot$avg_nEnh_KeptRight +
      infoBarplot$avg_nEnh_GainedFromTheOtherDomain
    
    acetilation_Final<-infoBarplot$avg_AcetilationEnh_KeptLeft +
      infoBarplot$avg_AcetilationEnh_KeptRight +
      infoBarplot$avg_AcetilationEnh_GainedFromTheOtherDomain  
    
  }
  
  #########################################################
  ##Only used on duplication explanations for now:
  
  ##nEnh_Gained tells you the number of NEW enh aquired thanks to the SV // NeoTad formation if Duplication
  nEnh_Gained<-infoBarplot$avg_nEnh_GainedFromTheOtherDomain
  
  ##Recall, kept only refers to intraTAD enhancers (we want to know with this, how many intra TAD enh duplicated, if any)
  ##If it is >0, it corresponds with the number of cognate duplicated enhancers for the target gene
  nEnh_cognate_Duplicated<-(infoBarplot$avg_nEnh_KeptLeft + infoBarplot$avg_nEnh_KeptRight) - 
                                (infoBarplot$avg_nEnh_InitialLeft + infoBarplot$avg_nEnh_InitialRight)
  
  ## To know whether Duplication is IntraTAD (if gene associated to both break), 
  ## if reg domain entirely duplicated (gene no association to breakpoint)
  ## Or if neoTAD happening, if gene associated just with one brekapoint (it means the other breakpoint falls in another domain)
  gene_Breakpoint<-infoBarplot$gene_Breakpoint 
  
  ##Modify regarding if it is for LOF or GOF
  #targetGene
  #pathoMechanism ##GOF or LOF
  #pathomechanism_ImpactOverGene ##long-range or direct

  textExplanation<-""
  
  if((pathoMechanism=="LOF")&&(pathomechanism_ImpactOverGene=="LongRange")){
    
    textExplanation<-paste(textExplanation,
                           "<b style='font-size:20px;'>Control vs Patient rearranged locus comparison </b>",
                           "<br><br>", 
                           targetGene, " sequence is not affected by the Structural Variant. However, there are potential long-range regulatory changes that could affect its expression.", 
                           "<br><br>",
                           targetGene," initially had ",
                           nEnh_Initial," enhancers  on its regulatory domain, ",
                           "and it lost ", round2(x = (nEnh_Initial - nEnh_Final),digits = 0),
                           " due to the ", sv_type, " (",round2(x = (abs(nEnh_Final-nEnh_Initial)/nEnh_Initial)*100, digits = 1) ,"% reduction). ",
                           "In addition, the total H3K27ac levels measured at its enhancers changed from ",
                           acetilation_Initial,
                           " to ",
                           acetilation_Final,
                           " RPGC (reads per genome coverage) ",
                           " (", round2(x = (abs(acetilation_Final-acetilation_Initial)/acetilation_Initial)*100, digits = 1),"% reduction).",
                           "<br><br>",
                           sep="")
    
    textExplanation<-paste(textExplanation,
                           text_lof1,
                           "<br><br>",
                           sep="")
  
  }else if((pathoMechanism=="GOF")&&(pathomechanism_ImpactOverGene=="LongRange")){
    
    #######################################
    ## Text for GOF && Long-range ##
    #######################################
    
    ##If initially the gene had no enhancers, saying that there were not, but not indicate an infinite gain of enh
    ##That makes no sense
    
    textExplanation<-paste(textExplanation,
                           "<b style='font-size:20px;'>Control vs Patient rearranged locus comparison </b>",
                           "<br><br>", 
                           targetGene, " sequence is not affected by the Structural Variant. However, there are potential long-range regulatory changes that could affect its expression.", 
                           "<br><br>",
                           targetGene,
                           " initially had ",
                           nEnh_Initial," enhancers  on its regulatory domain, ",
                           "and it gained ", round2(x = (nEnh_Final - nEnh_Initial),digits = 0),
                           " due to the ", sv_type,
                           sep = ""
    )
    
                           
    if(nEnh_Initial>0){
      ##So it makes sense to provide a ratio or percentage increase
      textExplanation<-paste(textExplanation,
                             " (", round2(x = (abs(nEnh_Final - nEnh_Initial)/nEnh_Initial)*100, digits = 1),
                             "% increase). In addition, the total H3K27ac levels measured at its enhancers changed from ",
                             acetilation_Initial,
                             " to ",
                             round2(x = acetilation_Final,digits = 1),
                             " RPGC (reads per genome coverage) (", round2(x = (abs(acetilation_Final-acetilation_Initial)/acetilation_Initial)*100, digits = 1),
                             "% increase). ",
                             "<br><br>",
                             sep="")
    }else{
      ##Skip ratios info, since if nEnh initial = 0, and there is a gain... the gain is infinite
      textExplanation<-paste(textExplanation,
                             "In addition, the total H3K27ac levels measured at its enhancers changed from ",
                             acetilation_Initial,
                             " to ",
                             round2(x = acetilation_Final, digits = 1),
                             " RPGC (reads per genome coverage).",
                             "<br><br>",
                             sep="")
    }
    
    # Not sure whether this is worth adding, putting this on user guide
    # #Ending sentence, if gain by TAD fusion, in deletion point it, as maybe it is something unexpected
    # if(sv_type == "Deletion"){
    #   textExplanation<-paste(textExplanation,
    #                          "The Deletion has caused the fusion of different TADs due to the removal of TAD boundaries.",
    #                          "<br><br>",
    #                          sep="")
    # }
    # else{
    #    textExplanation<-paste(textExplanation,
    #                           " ",
    #                           sep="")
    #  }
    

    
    textExplanation<-paste(textExplanation,
                           text_gof1,
                           "<br><br>",
                           sep="")
    
    
  }else if(pathomechanism_ImpactOverGene == "LongRange_geneDuplication"){
    ##---
    
    
    #######################################
    ## Text for GOF && Long-range ##
    #######################################
    
    ##If initially the gene had no enhancers, saying that there were not, but not indicate an infinite gain of enh
    ##That makes no sense
    
    textExplanation<-paste(textExplanation,
                           "<b style='font-size:20px;'>Control vs Patient rearranged locus comparison </b>",
                           "<br><br>", 
                           targetGene, " sequence is duplicated by the Structural Variant. In addition, there are potential long-range regulatory changes that could affect its expression.", 
                           "<br><br>", sep="")
    
    ##########################################
    ## Noticing if NeoTad, worth to mention
    ## We know neoTAD if nGained >0. Because gained come from different TAD and only way, by dupl mixing domains, hence NeoTAD
    
    if(gene_Breakpoint %in% c("Break1","Break2")){
      ##So the gene is associated just to one breakpoint, hence there is a regulatory domain merge
      ##Hence NeoTAD to be displayed
      textExplanation<-paste(textExplanation,
                             "Due to the Duplication, a new TAD is formed, where it is also now ",
                             targetGene,
                             ".",
                             sep="")    
      if(nEnh_Gained > 0 ){
        ##So at least a new enhancer gained from another domain.
        textExplanation<-paste(textExplanation,
                               " This new TAD contains ",
                               nEnh_Gained,
                               " enhancers from a different regulatory domain. Thus, not previously, in the one from ",
                              targetGene,
                               ".",
                               sep="")    
      }
      
      if(nEnh_cognate_Duplicated > 0 ){
        ##So at least a new enhancer gained from another domain.
        textExplanation<-paste(textExplanation,
                               " It also contains ",
                               nEnh_cognate_Duplicated,
                               " enhancers from ",
                               targetGene,
                               " cognate regulatory domain.",
                               sep="")
        
      }
      
      textExplanation<-paste(textExplanation,
                             " All the enhancers belonging to this new TAD have been duplicated. Hence, they appear both, at their normal domain, and also on the newly created.",
                             "<br><br>",
                             sep="")  
      
    }
    
    
    ####################################
    ##########Continuing as before +/-
    ####################################
    textExplanation<-paste(textExplanation,
          targetGene,
          " initially had ",
          nEnh_Initial," enhancers  on its regulatory domain.",
          " Now, there is a total of ", round2(x = (nEnh_Final),digits = 0),
          " enhancers around him, attending to the overall impact of the ", sv_type, ". ",
          sep="")       
    
    if(nEnh_Initial>0){
      ##So it makes sense to provide a ratio or percentage increase
      textExplanation<-paste(textExplanation,
                             "This implies an increase of ", round2(x =  (abs(nEnh_Final - nEnh_Initial)/nEnh_Initial)*100, digits = 1),
                             "%. In addition, the total H3K27ac levels measured at its enhancers changed from ",
                             acetilation_Initial,
                             " to ",
                             acetilation_Final,
                             " RPGC (reads per genome coverage) (",
                             round2(x = (abs(acetilation_Final-acetilation_Initial)/acetilation_Initial)*100, digits = 1),
                             "% increase). ",
                             "<br><br>",
                             sep="")
    }else{
      ##Skip ratios info, since if nEnh initial = 0, and there is a gain... the gain is infinite
      textExplanation<-paste(textExplanation,
                             "In addition, the total H3K27ac levels measured at its enhancers changed from ",
                             acetilation_Initial,
                             " to ",
                             acetilation_Final,
                             " RPGC (reads per genome coverage)",
                             "<br><br>",
                             sep="")
    }
    
    
    
    textExplanation<-paste(textExplanation,
                           text_gof1,
                           "<br><br>",
                           sep="")
    
    
  }else if (pathomechanism_ImpactOverGene =="Direct_LongRange_geneDuplication"){
    #---
    #For now same explanation that for LongRange_geneDuplication
    ##---
    
    
    ##---
    
    
    #######################################
    ## Text for GOF && Long-range ##
    #######################################
    
    ##If initially the gene had no enhancers, saying that there were not, but not indicate an infinite gain of enh
    ##That makes no sense
    
    textExplanation<-paste(textExplanation,
                           "<b style='font-size:20px;'>Control vs Patient rearranged locus comparison </b>",
                           "<br><br>", 
                           targetGene, " sequence is duplicated by the Structural Variant. In addition, there are potential long-range regulatory changes that could affect its expression.", 
                           "<br><br>", sep="")
    
    ##########################################
    ## Noticing if NeoTad, worth to mention
    ## We know neoTAD if nGained >0. Because gained come from different TAD and only way, by dupl mixing domains, hence NeoTAD
    
    if(gene_Breakpoint %in% c("Break1","Break2")){
      ##So the gene is associated just to one breakpoint, hence there is a regulatory domain merge
      ##Hence NeoTAD to be displayed
      textExplanation<-paste(textExplanation,
                             "Due to the Duplication, a new TAD is formed, where it is also now ",
                             targetGene,
                             ".",
                             sep="")    
      if(nEnh_Gained > 0 ){
        ##So at least a new enhancer gained from another domain.
        textExplanation<-paste(textExplanation,
                               " This new TAD contains ",
                               nEnh_Gained,
                               " enhancers from a different regulatory domain. Thus, not previously, in the one from ",
                               targetGene,
                               ".",
                               sep="")    
      }
      
      if(nEnh_cognate_Duplicated > 0 ){
        ##So at least a new enhancer gained from another domain.
        textExplanation<-paste(textExplanation,
                               " It also contains ",
                               nEnh_cognate_Duplicated,
                               " enhancers from ",
                               targetGene,
                               " cognate regulatory domain.",
                               sep="")
        
      }
      
      textExplanation<-paste(textExplanation,
                             " All the enhancers belonging to this new TAD have been duplicated. Hence, they appear both, at their normal domain, and also on the newly created.",
                             "<br><br>",
                             sep="")  
      
    }
    
    
    ####################################
    ##########Continuing as before +/-
    ####################################
    textExplanation<-paste(textExplanation,
                           targetGene,
                           " initially had ",
                           nEnh_Initial," enhancers  on its regulatory domain.",
                           " Now, there is a total of ", round2(x = (nEnh_Final),digits = 0),
                           " enhancers around him, attending to the overall impact of the ", sv_type,
                           sep="")       
    
    if(nEnh_Initial>0){
      ##So it makes sense to provide a ratio or percentage increase
      textExplanation<-paste(textExplanation,
                             " (", round2(x = (abs(nEnh_Final - nEnh_Initial)/nEnh_Initial)*100, digits = 1),
                             "% increase). In addition, the total H3K27ac levels measured at its enhancers changed from ",
                             acetilation_Initial,
                             " to ",
                             acetilation_Final,
                             " RPGC (reads per genome coverage) (", round2(x = (abs(acetilation_Final-acetilation_Initial)/acetilation_Initial)*100, digits = 1),
                             "% increase).",
                             "<br><br>",
                             sep="")
    }else{
      ##Skip ratios info, since if nEnh initial = 0, and there is a gain... the gain is infinite
      textExplanation<-paste(textExplanation,
                             ". In addition, the total H3K27ac levels measured at its enhancers changed from ",
                             acetilation_Initial,
                             " to ",
                             acetilation_Final,
                             " RPGC (reads per genome coverage).",
                             "<br><br>",
                             sep="")
    }
    
    
    
    textExplanation<-paste(textExplanation,
                           text_gof1,
                           "<br><br>",
                           sep="")
    

  }else if(pathomechanism_ImpactOverGene == "Direct_geneDuplication"){
    
    textExplanation<-paste(textExplanation,
                           "<b style='font-size:20px;'>Control vs Patient rearranged locus comparison </b>",
                           "<br><br>", 
                           targetGene, " sequence is duplicated by the Structural Variant.", 
                           "<br><br>", 
                           sep="")
    
    textExplanation<-paste(textExplanation,
                            "We are assuming the Structural Variant is happening in heterozygosis, since it is the most common scenario, hence only 1 allelle is affected. ", 
                           "<br><br>", sep="")
    
  }else if((pathomechanism_ImpactOverGene == "Direct_geneDeletion")){
    textExplanation<-paste(textExplanation,
                           "<b style='font-size:20px;'>Control vs Patient rearranged locus comparison </b>",
                           "<br><br>", 
                           targetGene, " sequence has been deleted. Thus, its expression will have been directly affected by the Structural Variant, without the need to consider the changes on the regulatory landscape (enhancers).",
                           "<br><br>",
                           sep="")
    
    textExplanation<-paste(textExplanation,
                           "We are assuming the Structural Variant is happening in heterozygosis, since it is the most common scenario, hence only 1 allelle is affected. ", 
                           "<br><br>",
                           sep="")
    
  }else if(pathomechanism_ImpactOverGene == "Direct_geneTruncation"){
    textExplanation<-paste(textExplanation,
                           "<b style='font-size:20px;'>Control vs Patient rearranged locus comparison </b>",
                           "<br><br>", 
                           targetGene, " sequence is truncated. Thus, its expression should be directly affected by the Structural Variant, without the need to consider the changes on the regulatory landscape (enhancers).",
                           "<br><br>",
                           sep="")
    
    textExplanation<-paste(textExplanation,
                           "We are assuming the Structural Variant is happening in heterozygosis, since it is the most common scenario, hence only 1 allelle is affected. ", 
                           "<br><br>", 
                           sep="")
    
    #Info number of truncated genes
    #For truncated genes, TAD map does not matter, gene position does not vary, hence regardless tad map or condition genes truncated.
    trunc_genes<-unique(patientResults$resultsPerPhase_secondaryInfo[[1]][[1]]$affectedGenes$brokenGenes)
    if(length(trunc_genes)>1){
      ##Indicate that a potential fusion transcript could have been created
      textExplanation<-paste(textExplanation,
                             "The Structural Variant is truncating ", length(trunc_genes)," genes: ", targetGene,
                             " and ", paste0(trunc_genes[trunc_genes!=targetGene],collapse=", ", sep=""),
                             ". Thus, a potential fusion protein could be generated. This should be evaluated in more detail if it is considered relevant.",
                             "<br><br>", 
                             sep="")
    }
    
  }
  
  ##Wrapping text explanation on the html
  
  regLanChanges<-paste(regLanChanges,
                       "<div class='explanationEnh'>",
                       paste("<p style='text-align:justify;'>",
                            textExplanation,
                             "</p>",
                            sep=""),
                       "</div>",
                       sep="")

  ##CLOSING Section
  regLanChanges<-paste(regLanChanges,
                       "</div>",
                       sep="")
  
  return(regLanChanges)
  
}