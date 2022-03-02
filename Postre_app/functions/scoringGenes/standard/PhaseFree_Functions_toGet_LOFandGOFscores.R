##################################################################################
## Functions to compute LOF-GOF scores
## In the Phase-free scenario
## Phase-free scenario (neither Gene Expression nor Enhancer activity considered)
###################################################################################

##For these set of functions
## We do not care for expression or enhancer activity

####################################################################################
## TO evaluate Direct EFFECTS (as gene truncation, or gene deletion) ###############
####################################################################################
eval_lof_directEffect_phase_FREE_score<-function(info_gene, geneTransversalData, phenoScore){
  
  # ##To evaluate the effect over either Deleted or Truncated genes
  
  ##################################################################################
  #### Si no esta expresado, pero su phenoScore es alto y es HI, ha de puntuar alto
  ##En algun momento se causara dicha patologia a lo largo del desarrollo
  
  
  ################################################
  ## Score Polycomb
  
  scorePolycomb<-info_gene[,"polyComb_score"]
  
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
    
  }##Here as we do not know the gene expression, we do not pay attention to its expression on the current phase
  
  
  ############################################
  ## Assembling Score -- Specific for RunMode
  ############################################
  
  if(max(c(scoreDosageSensitive,scorePolycomb),na.rm = TRUE)>=0.75){
    ##At least 1 higher or equal than 0.75, so considerable either polyC or dosageSensitivity
    gene_PreConditions<-(scoreDosageSensitive + scorePolycomb)/2
    
  }else{
    gene_PreConditions<-0
  }
  
  
  finalScore<-(phenoScore + gene_PreConditions)/2 
  
  return(finalScore)
}

####################################################################################
## TO evaluate Long Range EFFECTS                                    ###############
####################################################################################

##For now it is going to be the same criterira as GeneBroken
##Since for whatever gene involved, if it appears on the analysis it is because its regulatory domain
##Has been altered
##Maybe in the future we can consider taking into account the % of the TAD that is disrupted
##But I don't find it that essential for now
##But well, if 99% of the TAD is mantained... it is likely not LOF long range effect
##Whereas for a broken gene, sure, there is an effect

eval_lof_indirectEffect_phase_FREE_score<-function(info_gene, geneTransversalData, phenoScore){
  
  ##For now not predicting for phase free
  ##If want the data, go to backup before 28 oct 2021
  
  return(0)
  
}

########################################################
## evaluation Gain Of Function only by Gene Duplication
########################################################
eval_Gof_directEffect_phase_FREE_score<-function(info_gene, geneTransversalData, phenoScore){
  
  ##We only have polycomb here......
  ##Because of course enh or expression info cannot be used
  ##Pretty far from being optiMAL
  
  ##Evaluacion de si polycomb o no in the promoter
  ##Polycomb genes the ones with morphogenic skills, so relevant
  # scorePolycomb<-0
  scorePolycomb<-info_gene[,"polyComb_score"]
  
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
  
  if(max(c(scoreDosageSensitive,scorePolycomb),na.rm = TRUE)>=0.75){
    ##At least 1 higher or equal than 0.75, so considerable either polyC or dosageSensitivity
    gene_PreConditions<-(scoreDosageSensitive + scorePolycomb)/2
    
  }else{
    gene_PreConditions<-0
  }
  
  finalScore<-(phenoScore + gene_PreConditions)/2 
  
  #finalScore<-gene_PreConditions
  
  return(finalScore)
}







