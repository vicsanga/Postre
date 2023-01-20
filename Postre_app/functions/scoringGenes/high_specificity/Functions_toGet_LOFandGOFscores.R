##############################################
## FUNCTIONS to get LOF and GOF scores
##############################################
source("functions/scoringGenes/auxFunctions/auxFunctionsToComputePathogenicScores.R", local = TRUE)

eval_lof_indirectEffect_score<-function(matrixPhase,geneTransversalData, phenoScore, threshold_MaxExpresion, threshold_MinExpresion,
                                        minRatioEnhBalance, maxRatioEnhBalance){
  
  ##ScoreExp
  scoreExp<-scoreExp_compute(matrixPhase = matrixPhase,
                             threshold_MaxExpresion = threshold_MaxExpresion, 
                             threshold_MinExpresion = threshold_MinExpresion)
  #ScoreEnhancer changes
  scoreEnh<-scoreEnhancer_LOF_compute(matrixPhase = matrixPhase)
  
  #ScorePolycomb
  scorePolycomb<-scorePolycomb_compute(matrixPhase = matrixPhase)
  
  #ScoreDosageSensitivity LOF
  scoreDosageSensitive<-scoreDS_LOF_compute(geneTransversalData = geneTransversalData)
  
  ##GeneFeatures STANDARD MODE
  geneFeatures<-scoreGeneFeatures_highSpecificityMode(scoreDosageSensitive = scoreDosageSensitive,
                                                      scorePolycomb = scorePolycomb)
  
  finalScore<-(geneFeatures + scoreEnh + phenoScore + scoreExp)/4
  
  return(list("finalScore"=finalScore,
              "geneEnhancerScore"=scoreEnh,
              "genePhenoScore"=phenoScore,
              "geneFeaturesScore"=geneFeatures,
              "dosageSensitivityScore"=scoreDosageSensitive,
              "polycombScore"=scorePolycomb,
              "geneExpressionScore"=scoreExp
  ))
  
}

################################################
## evaluation Gain Of Function
eval_Gof_indirectEffect_score<-function(matrixPhase,geneTransversalData, phenoScore, threshold_MaxExpresion, threshold_MinExpresion,
                                        minRatioEnhBalance, maxRatioEnhBalance){
  
  ##If nEnhGained = 0. All the potential enhancer increase, come from intraTAD enh Duplication
  ##So in this case is important to consider the gene expression levels
  #Because if a gene not expressed with some enh on its TAD... if the enh duplicated... we do not expect much effect from this enh on the poorly exp gene
  #However if the gene is considerably expressed, probably it is regulated by this enh and its expression will get boosted
  
  
  ##If there is a gain in enhancers. As far as the gene is polycomb, we do not care about its expression
  ##If it is not exp, it can get activated. And if it is active its expression can be boosted if it gains more enh than it loses.
  nEnhGained<-matrixPhase[,"nEnhancers_gained"] 
  
  scoreExp<-0
  
  caringAboutGeneExpresion<-TRUE
  
  if(nEnhGained > 0){
    ##So new enhancers gained
    #Do not care about gene expression
    #An unactive gene can be activated, and an already active can be upregulated
    scoreExp<-1
    caringAboutGeneExpresion<-FALSE
    
  }else if(nEnhGained == 0){
    ##So no new or ectopic enhancers gained. Gene expression important
    ##If gene not expressed, the duplication of intraTAD enhancers...not considered relevant, not likely regulating it
    ##ScoreExp
    scoreExp<-scoreExp_compute(matrixPhase = matrixPhase,
                               threshold_MaxExpresion = threshold_MaxExpresion, 
                               threshold_MinExpresion = threshold_MinExpresion)
    
  }
  
  #ScoreEnhancer changes
  scoreEnh<-scoreEnhancer_GOF_compute(matrixPhase = matrixPhase)
  
  #ScorePolycomb
  scorePolycomb<-scorePolycomb_compute(matrixPhase = matrixPhase)
  
  ##ScoreDosageSensitivity based on TS scores
  scoreDosageSensitive<-scoreDS_GOF_compute(geneTransversalData = geneTransversalData)
  
  
  ##GeneFeatures STANDARD MODE
  geneFeatures<-scoreGeneFeatures_highSpecificityMode(scoreDosageSensitive = scoreDosageSensitive,
                                                      scorePolycomb = scorePolycomb)
  
  
  if(caringAboutGeneExpresion==TRUE){
    finalScore<-(geneFeatures + scoreEnh + phenoScore + scoreExp)/4
    
    return(list("finalScore"=finalScore,
                "geneEnhancerScore"=scoreEnh,
                "genePhenoScore"=phenoScore,
                "geneFeaturesScore"=geneFeatures,
                "dosageSensitivityScore"=scoreDosageSensitive,
                "polycombScore"=scorePolycomb,
                "geneExpressionScore"=scoreExp
    ))
    
  }else if (caringAboutGeneExpresion==FALSE){
    finalScore<-(geneFeatures + scoreEnh + phenoScore)/3
    
    return(list("finalScore"=finalScore,
                "geneEnhancerScore"=scoreEnh,
                "genePhenoScore"=phenoScore,
                "geneFeaturesScore"=geneFeatures,
                "dosageSensitivityScore"=scoreDosageSensitive,
                "polycombScore"=scorePolycomb,
                "geneExpressionScore"=NA
    ))
  }
  
}


########################################################
## evaluation Gain Of Function only by Gene Duplication
## active gene, expressed, duplicated (enh info ignored)
########################################################
eval_Gof_directEffect_score<-function(matrixPhase,geneTransversalData, phenoScore, threshold_MaxExpresion, threshold_MinExpresion,
                                      minRatioEnhBalance, maxRatioEnhBalance){
  
  ##being consistently expressed gives points
  ##Because here we are not considering enhancers!
  ##And the duplicaton of a non expressed gene can not be considered as pathogenic by GOF
  
  ##ScoreExp
  scoreExp<-scoreExp_compute(matrixPhase = matrixPhase,
                             threshold_MaxExpresion = threshold_MaxExpresion, 
                             threshold_MinExpresion = threshold_MinExpresion)
  
  #ScorePolycomb
  scorePolycomb<-scorePolycomb_compute(matrixPhase = matrixPhase)
  
  ##ScoreDosageSensitivity based on TS scores
  scoreDosageSensitive<-scoreDS_GOF_compute(geneTransversalData = geneTransversalData)
  
  ##GeneFeatures STANDARD MODE
  geneFeatures<-scoreGeneFeatures_highSpecificityMode(scoreDosageSensitive = scoreDosageSensitive,
                                                      scorePolycomb = scorePolycomb)
  
  finalScore<-(geneFeatures + phenoScore + scoreExp)/3
  
  return(list("finalScore"=finalScore,
              "geneEnhancerScore"=NA,
              "genePhenoScore"=phenoScore,
              "geneFeaturesScore"=geneFeatures,
              "dosageSensitivityScore"=scoreDosageSensitive,
              "polycombScore"=scorePolycomb,
              "geneExpressionScore"=scoreExp
  ))
  
}


####################################################################################
## TO evaluate Direct EFFECTS (as gene truncation, or gene deletion) ###############
####################################################################################
eval_lof_directEffect_score<-function(matrixPhase,geneTransversalData, phenoScore, threshold_MaxExpresion, threshold_MinExpresion){
  
  # ##To evaluate the effect over either Deleted or Truncated genes
  
  ##ScoreExp
  scoreExp<-scoreExp_compute(matrixPhase = matrixPhase,
                             threshold_MaxExpresion = threshold_MaxExpresion, 
                             threshold_MinExpresion = threshold_MinExpresion)
  
  #ScorePolycomb
  scorePolycomb<-scorePolycomb_compute(matrixPhase = matrixPhase)
  
  #ScoreDosageSensitivity LOF
  scoreDosageSensitive<-scoreDS_LOF_compute(geneTransversalData = geneTransversalData)
  
  ##GeneFeatures STANDARD MODE
  geneFeatures<-scoreGeneFeatures_highSpecificityMode(scoreDosageSensitive = scoreDosageSensitive,
                                                      scorePolycomb = scorePolycomb)
  
  finalScore<-(geneFeatures + phenoScore + scoreExp)/3
  
  return(list("finalScore"=finalScore,
              "geneEnhancerScore"=NA,
              "genePhenoScore"=phenoScore,
              "geneFeaturesScore"=geneFeatures,
              "dosageSensitivityScore"=scoreDosageSensitive,
              "polycombScore"=scorePolycomb,
              "geneExpressionScore"=scoreExp
  ))
  
}





