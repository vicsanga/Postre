##############################################
## FUNCTIONS to get LOF and GOF scores
## eval_lof_indirectEffect_score
## and eval_gof_indrectEffect_score
##############################################

eval_lof_indirectEffect_score<-function(matrixPhase,geneTransversalData, phenoScore, threshold_MaxExpresion, threshold_MinExpresion,
                                        minRatioEnhBalance, maxRatioEnhBalance){
  
  
  ##For a LOF if a gene not expressed makes no sense to consider anything else
  ##Hence
  ##IF gene poorly expressed, no need to consider anything else, for a LOF context gene being expressed is Essential
  if(matrixPhase[,"FPKM"]<threshold_MinExpresion){
    finalScore <- 0
    return(finalScore)
  }
  
  
  ##being consistently (highly) expressed gives points
  scoreExp<-0
  if(matrixPhase[,"FPKM"]>=threshold_MaxExpresion){
    scoreExp<-1
    
  }else if(matrixPhase[,"FPKM"]>threshold_MinExpresion){
    ###If expression between threshold_MinExpresion and threshold_MaxExpresion
    ##Let's normalize the value to a 0-1 scale
    expValue<-matrixPhase[,"FPKM"]
    ##Normalizing to 0-1 scale
    scoreExp<-(expValue-threshold_MinExpresion)/(threshold_MaxExpresion-threshold_MinExpresion)
    
  }
  ##Else (hence exp below threshold minExpresion) score exp kept to 0
  
  
  ################################################
  ## Score Polycomb
  
  scorePolycomb<-matrixPhase[,"polyComb_score"]
  
  if(scorePolycomb < 0){
    scorePolycomb <- 0 ##for those with -1. Non computed polyC score
  }
  
  
  ######################
  ## Gene Dosage Sensitivity Evaluation
  
  scoreDosageSensitive<-0
  
  ##En base a HI scores
  hi_scores<-c(geneTransversalData[,"nature_HI_score"], geneTransversalData[,"huang_HI_score"], geneTransversalData[,"clinGene_HI_score"])
  
  max_hi_score<-round(x = max(hi_scores), digits = 2)##Tomamos el valor máximo de los dos papers como referencia
  
  if(max_hi_score != -1){
    
    ##So in this situation at least for one source there is a HI score
    ##Computed in any of the two papers
    ##The score will reach max 1 if the gene is super HI
    
    if(max_hi_score >= 0.75){
      ##if moderatly high consider it, on the contrary, do not take into account
      scoreDosageSensitive<-scoreDosageSensitive + max_hi_score
    }
    
    
  }
  
  
  ######################################
  ##Enh balance
  ##Regarding NUMBER of enhancers
  ##Regarding ACETILATION of enhancers
  ######################################
  
  ############################
  ##Regarding NUMBER of enhancers
  scoreEnhBalance<-0
  ##ratio nInitial/nFinal (for LOF)
  nInitial<-matrixPhase[,"nEnhancers_initial"]
  
  nFinal<-matrixPhase[,"nEnhancers_kept"] ## + matrixPhase[,"nEnhancers_gained"] ##Quitamos lo que gana por los casos raros en que empata
  
  ## if...else Ladder ## if...else multiple conditions
  if ((nFinal > 0) && (nInitial > 0)) { #Scenario 1/4
    ##there were enh initially and there are enh left
    
    ratioEnhBalance<- nInitial/nFinal #LOF SPECIFIC ADJ
    
  } else if ((nFinal == 0) && (nInitial == 0)) { #Scenario 2/4
    
    ratioEnhBalance<-1
    
  } else if ((nFinal == 0) && (nInitial > 0)){ #Scenario 3/4
    ##there were some enhancers initially, but everything lost
    
    ratioEnhBalance<-maxRatioEnhBalance #LOF SPECIFIC ADJ
    
  } else if ((nFinal > 0) && (nInitial == 0 )){ #Scenario 4/4
    ##there were no enhancers initially, and there is a gain
    
    ratioEnhBalance<-minRatioEnhBalance #LOF SPECIFIC ADJ
    
  }
  
  ##converting ratioEnhBalance to 0-1 scoreEnhBalance
  if(ratioEnhBalance >= maxRatioEnhBalance){
    scoreEnhBalance <- 1
  }else if(ratioEnhBalance <= minRatioEnhBalance){
    scoreEnhBalance <- 0
  }else{
    ##Hence ratio enhBalance between max and min ratio
    ##scaling to 0-1
    scoreEnhBalance<-(ratioEnhBalance - minRatioEnhBalance)/
      (maxRatioEnhBalance - minRatioEnhBalance)
  }
  
  
  ############################
  ##Regarding ACETILATION of enhancers
  scoreAcetilationEnhBalance<-0
  ##ratio nInitial/nFinal (for LOF)
  nInitial<-matrixPhase[,"enhancers_acetilation_initial"]
  
  nFinal<-matrixPhase[,"enhancers_acetilation_kept"] ##+ matrixPhase[,"enhancers_acetilation_gained"] 
  
  ## if...else Ladder ## if...else multiple conditions
  if ((nFinal > 0) && (nInitial > 0)) { #Scenario 1/4
    ##there were enh initially and there are enh left
    
    ratioEnhBalance<- nInitial/nFinal #LOF SPECIFIC ADJ
    
  } else if ((nFinal == 0) && (nInitial == 0)) { #Scenario 2/4
    
    ratioEnhBalance<-1
    
  } else if ((nFinal == 0) && (nInitial > 0)){ #Scenario 3/4
    ##there were some enhancers initially, but everything lost
    
    ratioEnhBalance<-maxRatioEnhBalance #LOF SPECIFIC ADJ
    
  } else if ((nFinal > 0) && (nInitial == 0 )){ #Scenario 4/4
    ##there were no enhancers initially, and there is a gain
    
    ratioEnhBalance<-minRatioEnhBalance #LOF SPECIFIC ADJ
    
  }
  
  ##converting ratioEnhBalance to 0-1 scoreAcetilationEnhBalance
  if(ratioEnhBalance >= maxRatioEnhBalance){
    scoreAcetilationEnhBalance <- 1
  }else if(ratioEnhBalance <= minRatioEnhBalance){
    scoreAcetilationEnhBalance <- 0
  }else{
    ##Hence ratio enhBalance between max and min ratio
    ##scaling to 0-1
    scoreAcetilationEnhBalance<-(ratioEnhBalance - minRatioEnhBalance)/
      (maxRatioEnhBalance - minRatioEnhBalance)
  }
  
  
  ###########
  # finalScore<-(phenoScore+scoreDosageSensitive+scoreEnhBalance)/3
  
  ###After introducing scoreExp
  # finalScore<-(scoreExp+phenoScore+scoreDosageSensitive+scoreEnhBalance)/4
  
  #####Score Lof no puede pesar tanto como el phenoScore,because even if it is HI it does not mean it needs to be 
  ##causing the phenotype
  
  
  ############################################
  ## Assembling Score -- Specific for RunMode
  ############################################
  
  # gene_PreConditions<-(scoreExp + scoreDosageSensitive + scorePolycomb)/3
  if(max(c(scoreDosageSensitive,scorePolycomb),na.rm = TRUE)>=0.75){
    ##At least 1 higher or equal than 0.75, so considerable either polyC or dosageSensitivity
    gene_PreConditions<-(scoreExp + scoreDosageSensitive + scorePolycomb)/3
  }else{
    gene_PreConditions<-0
  }
  
  # gene_EnhLandscapeConditions<-(scoreEnhBalance + scoreAcetilationEnhBalance)/2
  ##Ignoring changes in enh number, considering only AcetilationEnhBalance
  gene_EnhLandscapeConditions<-scoreAcetilationEnhBalance
  
  finalScore<-(gene_PreConditions + gene_EnhLandscapeConditions + phenoScore)/3
  
  
  
  return(finalScore)
}

################################################
## evaluation Gain Of Function
eval_Gof_indirectEffect_score<-function(matrixPhase,geneTransversalData, phenoScore, threshold_MaxExpresion, threshold_MinExpresion,
                                        minRatioEnhBalance, maxRatioEnhBalance){
  ##########
  # Duda de si contar o No el phenoscore,
  # De momento no
  #por el tema de que nos movemos en un contexto de que ese gen es segun OMIM y MGI relevante
  ####
  
  ##If nEnhGained = 0. All the potential enhancer increase, come from intraTAD enh Duplication
  ##So in this case is important to consider the gene expression levels
  #Because if a gene not expressed with some enh on its TAD... if the enh duplicated... we do not expect much effect from this enh on the poorly exp gene
  #However if the gene is considerably expressed, probably it is regulated by this enh and its expression will get boosted
  
  
  ##If there is a gain in enhancers. As far as the gene is polycomb, we do not care about its expression
  ##If it is not exp, it can get activated. And if it is active its expression can be boosted if it gains more enh than it loses.
  nEnhGained<-matrixPhase[,"nEnhancers_gained"] 
  
  scoreExp<-0
  
  if(nEnhGained > 0){
    ##So new enhancers gained
    #Do not care about gene expression
    #An unactive gene can be activated, and an already active can be upregulated
    scoreExp<-1
    
  }else if(nEnhGained == 0){
    ##So no new or ectopic enhancers gained. Gene expression important
    ##If gene not expressed, the duplication of intraTAD enhancers...not considered relevant, not likely regulating it
    ##Hence return score 0
    if(matrixPhase[,"FPKM"]<threshold_MinExpresion){
      finalScore <- 0
      return(finalScore)
    }
    
    
    if(matrixPhase[,"FPKM"]>=threshold_MaxExpresion){
      scoreExp<-1
    }else if(matrixPhase[,"FPKM"]>threshold_MinExpresion){
      ###If expression between threshold_MinExpresion and threshold_MaxExpresion
      ##Let's normalize the value to a 0-1 scale
      expValue<-matrixPhase[,"FPKM"]
      ##Normalizing to 0-1 scale
      scoreExp<-(expValue-threshold_MinExpresion)/(threshold_MaxExpresion-threshold_MinExpresion)
      
    }
    ##Else (hence exp below threshold minExpresion) score exp kept to 0
    ##Gene not expressed, the duplication of intraTAD enhancers...not apparently relevant
  }
  
  ##Evaluacion de si polycomb o no in the promoter (hay genes no expresados sin polycomb)
  # scorePolycomb<-0
  scorePolycomb<-matrixPhase[,"polyComb_score"]
  
  if(scorePolycomb < 0){
    scorePolycomb <- 0 ##for those with -1. Non computed polyC score
  }
  
  ###########################################################################################################################
  ## Regarding geneDosage sensitivity, HI used  as proxy for gene dosage sensitive in general, for both DownReg and UpReg
  ###########################################################################################################################
  
  scoreDosageSensitive<-0
  
  ##En base a HI scores
  hi_scores<-c(geneTransversalData[,"nature_HI_score"], geneTransversalData[,"huang_HI_score"], geneTransversalData[,"clinGene_HI_score"])
  
  max_hi_score<-round(x = max(hi_scores), digits = 2)##Tomamos el valor máximo de los dos papers como referencia
  
  if(max_hi_score != -1){
    
    ##So in this situation at least for one source there is a HI score
    ##Computed in any of the two papers
    ##The score will reach max 1 if the gene is super HI
    
    if(max_hi_score >= 0.75){
      ##if moderatly high consider it, on the contrary, do not take into account
      scoreDosageSensitive<-scoreDosageSensitive + max_hi_score
    }
    
    
  }
  
  
  ########################
  ##Enh balance
  ##Regarding NUMBER of enhancers
  ##Regarding ACETILATION of enhancers
  ########################
  
  ############################
  ##Regarding NUMBER of enhancers
  scoreEnhBalance<-0
  ##ratio nFinal/nInitial (for GOF)
  nInitial<-matrixPhase[,"nEnhancers_initial"]
  
  nFinal<-matrixPhase[,"nEnhancers_kept"] + matrixPhase[,"nEnhancers_gained"] 
  
  ## if...else Ladder ## if...else multiple conditions
  if ((nFinal > 0) && (nInitial > 0)) { #Scenario 1/4
    ##there were enh initially and there are enh left
    
    ratioEnhBalance<- nFinal/nInitial ## GOF SPECIFIC ADJ
    
  } else if ((nFinal == 0) && (nInitial == 0)) { #Scenario 2/4
    
    ratioEnhBalance<-1
    
  } else if ((nFinal == 0) && (nInitial > 0)){ #Scenario 3/4
    ##there were some enhancers initially, but everything lost
    
    ratioEnhBalance<-minRatioEnhBalance #GOF SPECIFIC ADJ
    
  } else if ((nFinal > 0) && (nInitial == 0 )){ #Scenario 4/4
    ##there were no enhancers initially, and there is a gain
    
    ratioEnhBalance<-maxRatioEnhBalance #GOF SPECIFIC ADJ
    
  }
  
  ##converting ratioEnhBalance to 0-1 scoreEnhBalance
  if(ratioEnhBalance >= maxRatioEnhBalance){
    scoreEnhBalance <- 1
  }else if(ratioEnhBalance <= minRatioEnhBalance){
    scoreEnhBalance <- 0
  }else{
    ##Hence ratio enhBalance between max and min ratio
    ##scaling to 0-1
    scoreEnhBalance<-(ratioEnhBalance - minRatioEnhBalance)/
      (maxRatioEnhBalance - minRatioEnhBalance)
  }
  
  
  ############################
  ##Regarding ACETILATION of enhancers
  scoreAcetilationEnhBalance<-0
  ##ratio nInitial/nFinal (for LOF)
  nInitial<-matrixPhase[,"enhancers_acetilation_initial"]
  
  nFinal<-matrixPhase[,"enhancers_acetilation_kept"] + matrixPhase[,"enhancers_acetilation_gained"] 
  
  ###########
  ## if...else Ladder ## if...else multiple conditions
  if ((nFinal > 0) && (nInitial > 0)) { #Scenario 1/4
    ##there were enh initially and there are enh left
    
    ratioEnhBalance<- nFinal/nInitial ## GOF SPECIFIC ADJ
    
  } else if ((nFinal == 0) && (nInitial == 0)) { #Scenario 2/4
    
    ratioEnhBalance<-1
    
  } else if ((nFinal == 0) && (nInitial > 0)){ #Scenario 3/4
    ##there were some enhancers initially, but everything lost
    
    ratioEnhBalance<-minRatioEnhBalance #GOF SPECIFIC ADJ
    
  } else if ((nFinal > 0) && (nInitial == 0)){ #Scenario 4/4
    ##there were no enhancers initially, and there is a gain
    
    ratioEnhBalance<-maxRatioEnhBalance #GOF SPECIFIC ADJ
    
  }
  
  ##converting ratioEnhBalance to 0-1 scoreAcetilationEnhBalance
  if(ratioEnhBalance >= maxRatioEnhBalance){
    scoreAcetilationEnhBalance <- 1
  }else if(ratioEnhBalance <= minRatioEnhBalance){
    scoreAcetilationEnhBalance <- 0
  }else{
    ##Hence ratio enhBalance between max and min ratio
    ##scaling to 0-1
    scoreAcetilationEnhBalance<-(ratioEnhBalance - minRatioEnhBalance)/
      (maxRatioEnhBalance - minRatioEnhBalance)
  }
  
  ############################################
  ## Assembling Score -- Specific for RunMode
  ############################################
  
  # gene_PreConditions<-(scoreExp + scoreDosageSensitive + scorePolycomb)/3
  if(max(c(scoreDosageSensitive,scorePolycomb),na.rm = TRUE)>=0.75){
    ##At least 1 higher or equal than 0.75, so considerable either polyC or dosageSensitivity
    gene_PreConditions<-(scoreExp + scoreDosageSensitive + scorePolycomb)/3
  }else{
    gene_PreConditions<-0
  }
  
  
  # gene_EnhLandscapeConditions<-(scoreEnhBalance + scoreAcetilationEnhBalance)/2
  ##Ignoring changes in enh number, considering only AcetilationEnhBalance
  gene_EnhLandscapeConditions<-scoreAcetilationEnhBalance
  
  finalScore<-(gene_PreConditions + gene_EnhLandscapeConditions + phenoScore)/3
  
  return(finalScore)
}


########################################################
## evaluation Gain Of Function only by Gene Duplication
## active gene, expressed, duplicated (enh info ignored)
########################################################
eval_Gof_directEffect_score<-function(matrixPhase,geneTransversalData, phenoScore, threshold_MaxExpresion, threshold_MinExpresion,
                                      minRatioEnhBalance, maxRatioEnhBalance){
  ##########
  # Duda de si contar o No el phenoscore,
  # De momento no
  #por el tema de que nos movemos en un contexto de que ese gen es segun OMIM y MGI relevante
  ####
  ##being consistently expressed gives points
  ##Because here we are not considering enhancers!
  
  
  ##And the duplicaton of a non expressed gene can not be considered as pathogenic by GOF
  ##Hence
  if(matrixPhase[,"FPKM"]<threshold_MinExpresion){
    finalScore <- 0
    return(finalScore)
  }
  
  scoreExp<-0
  if(matrixPhase[,"FPKM"]>=threshold_MaxExpresion){
    scoreExp<-1
    
  }else if(matrixPhase[,"FPKM"]>threshold_MinExpresion){
    ###If expression between threshold_MinExpresion and threshold_MaxExpresion
    ##Let's normalize the value to a 0-1 scale
    expValue<-matrixPhase[,"FPKM"]
    ##Normalizing to 0-1 scale
    scoreExp<-(expValue-threshold_MinExpresion)/(threshold_MaxExpresion-threshold_MinExpresion)
    
  }
  ##Else (hence exp below threshold minExpresion) score exp kept to 0
  
  
  ##Evaluacion de si polycomb o no in the promoter
  ##Polycomb genes the ones with morphogenic skills, so relevant
  #scorePolycomb<-0
  scorePolycomb<-matrixPhase[,"polyComb_score"]
  
  if(scorePolycomb < 0){
    scorePolycomb <- 0 ##for those with -1. Non computed polyC score
  }
  
  ###########################################################################################################################
  ## Regarding geneDosage sensitivity, HI used  as proxy for gene dosage sensitive in general, for both DownReg and UpReg
  ###########################################################################################################################
  
  scoreDosageSensitive<-0
  
  ##En base a HI scores
  hi_scores<-c(geneTransversalData[,"nature_HI_score"], geneTransversalData[,"huang_HI_score"], geneTransversalData[,"clinGene_HI_score"])
  
  max_hi_score<-round(x = max(hi_scores), digits = 2)##Tomamos el valor máximo de los dos papers como referencia
  
  if(max_hi_score != -1){
    
    ##So in this situation at least for one source there is a HI score
    ##Computed in any of the two papers
    ##The score will reach max 1 if the gene is super HI
    
    if(max_hi_score >= 0.75){
      ##if moderatly high consider it, on the contrary, do not take into account
      scoreDosageSensitive<-scoreDosageSensitive + max_hi_score
    }
    
    
  }
  
  ############################################
  ## Assembling Score -- Specific for RunMode
  ############################################
  
  
  ###In this context I'm going to ignore the enh landscape changes
  # gene_PreConditions<-(scoreExp + scoreDosageSensitive + scorePolycomb)/3
  if(max(c(scoreDosageSensitive,scorePolycomb),na.rm = TRUE)>=0.75){
    ##At least 1 higher or equal than 0.75, so considerable either polyC or dosageSensitivity
    gene_PreConditions<-(scoreExp + scoreDosageSensitive + scorePolycomb)/3
  }else{
    gene_PreConditions<-0
  }
  
  #finalScore<-gene_PreConditions
  finalScore<-(gene_PreConditions + phenoScore)/2
  
  return(finalScore)
}


####################################################################################
## TO evaluate Direct EFFECTS (as gene truncation, or gene deletion) ###############
####################################################################################
eval_lof_directEffect_score<-function(matrixPhase,geneTransversalData, phenoScore, threshold_MaxExpresion, threshold_MinExpresion){
  
  # ##To evaluate the effect over either Deleted or Truncated genes
  
  ##For a LOF if a gene not expressed makes no sense to consider anything else
  ##Hence
  ##IF gene poorly expressed, no need to consider anything else, for a LOF context gene being expressed is Essential
  if(matrixPhase[,"FPKM"]<threshold_MinExpresion){
    finalScore <- 0
    return(finalScore)
  }
  
  
  ##being consistently (highly) expressed gives points
  scoreExp<-0
  if(matrixPhase[,"FPKM"]>=threshold_MaxExpresion){
    scoreExp<-1
    
  }else if(matrixPhase[,"FPKM"]>threshold_MinExpresion){
    ###If expression between threshold_MinExpresion and threshold_MaxExpresion
    ##Let's normalize the value to a 0-1 scale
    expValue<-matrixPhase[,"FPKM"]
    ##Normalizing to 0-1 scale
    scoreExp<-(expValue-threshold_MinExpresion)/(threshold_MaxExpresion-threshold_MinExpresion)
    
  }
  ##Else (hence exp below threshold minExpresion) score exp kept to 0
  
  
  ################################################
  ## Score Polycomb
  
  scorePolycomb<-matrixPhase[,"polyComb_score"]
  
  if(scorePolycomb < 0){
    scorePolycomb <- 0 ##for those with -1. Non computed polyC score
  }
  
  
  
  ######################
  ## Gene Dosage Sensitivity Evaluation
  
  scoreDosageSensitive<-0
  
  ##En base a HI scores
  hi_scores<-c(geneTransversalData[,"nature_HI_score"], geneTransversalData[,"huang_HI_score"], geneTransversalData[,"clinGene_HI_score"])
  
  max_hi_score<-round(x = max(hi_scores), digits = 2)##Tomamos el valor máximo de los dos papers como referencia
  
  if(max_hi_score != -1){
    
    ##So in this situation at least for one gene there is a HI score
    ##Computed in any of the two papers
    ##The score will reach max 1 if the gene is super HI
    
    if(max_hi_score >= 0.75){
      ##if moderatly high consider it, on the contrary, do not take into account
      scoreDosageSensitive<-scoreDosageSensitive + max_hi_score
    }
    
    
    
  }
  

  ############################################
  ## Assembling Score -- Specific for RunMode
  ############################################
  
  # gene_PreConditions<-(scoreExp + scoreDosageSensitive + scorePolycomb)/3
  if(max(c(scoreDosageSensitive,scorePolycomb),na.rm = TRUE)>=0.75){
    ##At least 1 higher or equal than 0.75, so considerable either polyC or dosageSensitivity
    gene_PreConditions<-(scoreExp + scoreDosageSensitive + scorePolycomb)/3
  }else{
    gene_PreConditions<-0
  }
  
  finalScore<-(gene_PreConditions + phenoScore)/2
  
  return(finalScore)
}





