###########################################################
## Auxiliar and common functions to avoid duplicated code
###########################################################

scoreExp_compute<-function(matrixPhase, threshold_MaxExpresion, threshold_MinExpresion){
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
  return(scoreExp)
}

scoreDS_LOF_compute<-function(geneTransversalData){
  ######################
  ## Gene Dosage Sensitivity Evaluation
  
  scoreDosageSensitive<-0
  
  ##En base a HI scores
  hi_scores<-c(geneTransversalData[,"nature_HI_score"], geneTransversalData[,"huang_HI_score"], geneTransversalData[,"clinGene_HI_score"], geneTransversalData[,"cell_HI_score"])
  
  max_hi_score<-round(x = max(hi_scores), digits = 2)
  
  ##The score will reach max 1 if the gene is super HI
  
  if(max_hi_score >= 0.85){
    ##if moderatly high consider it, on the contrary, do not take into account
    scoreDosageSensitive<-scoreDosageSensitive + max_hi_score
  }
  
  return(scoreDosageSensitive)
}

scoreDS_GOF_compute<-function(geneTransversalData){
  scoreDosageSensitive<-0
  
  ##En base a Triplosensitivity scores
  hi_scores<-geneTransversalData[,"cell_TriploSense_score"]
  
  max_hi_score<-round(x = max(hi_scores), digits = 2)
  
  ##The score will reach max 1 if the gene is super TS
  if(max_hi_score >= 0.85){
    ##if moderatly high consider it, on the contrary, do not take into account
    scoreDosageSensitive<-scoreDosageSensitive + max_hi_score
  }
    
  return(scoreDosageSensitive)
}


scoreEnhancer_LOF_compute<-function(matrixPhase){
  ######################################
  ## Enh balance
  ## Regarding ACETILATION of enhancers
  ######################################
  
  ############################
  ##Regarding ACETILATION of enhancers
  scoreAcetilationEnhBalance<-0
  ##ratio nInitial/nFinal (for LOF)
  nInitial<-matrixPhase[,"enhancers_acetilation_initial"]
  
  nFinal<-matrixPhase[,"enhancers_acetilation_kept"] ##+ matrixPhase[,"enhancers_acetilation_gained"] 
  
  if(nInitial > nFinal){
    ##Perdida de acetilacion
    scoreAcetilationEnhBalance<-1
    
  }
  
  return(scoreAcetilationEnhBalance)
}

scoreEnhancer_GOF_compute<-function(matrixPhase){
  ########################
  ##Enh balance
  ##Regarding ACETILATION of enhancers
  ########################
  
  ############################
  ##Regarding ACETILATION of enhancers
  scoreAcetilationEnhBalance<-0
  ##ratio nInitial/nFinal (for LOF)
  nInitial<-matrixPhase[,"enhancers_acetilation_initial"]
  
  nFinal<-matrixPhase[,"enhancers_acetilation_kept"] + matrixPhase[,"enhancers_acetilation_gained"] 
  
  
  if(nFinal > nInitial){
    ##Ganancia de acetilacion
    scoreAcetilationEnhBalance<-1
  } 
  
  return(scoreAcetilationEnhBalance)
  
}

scorePolycomb_compute<-function(matrixPhase){
  #ScorePolycomb
  scorePolycomb<-matrixPhase[,"polyComb_score"]
  if(scorePolycomb < 0){
    scorePolycomb <- 0 ##for those with -1. Non computed polyC score
  }
  return(scorePolycomb)
}

scoreGeneFeatures_standardMode<-function(scoreDosageSensitive, scorePolycomb){
  if(max(c(scoreDosageSensitive,scorePolycomb)) >= 0.85){
    ##At least 1 higher or equal than 0.85, so considerable either polyC or dosageSensitivity
    geneFeatures<-(scoreDosageSensitive + scorePolycomb)/2
  }else{
    geneFeatures<-0
  }
  return(geneFeatures)
}

scoreGeneFeatures_highSpecificityMode<-function(scoreDosageSensitive, scorePolycomb){
  if(all(c(scoreDosageSensitive,scorePolycomb) >= 0.85)){
    geneFeatures<-(scoreDosageSensitive + scorePolycomb)/2
  }else{
    geneFeatures<-0
  }
  return(geneFeatures)
}

