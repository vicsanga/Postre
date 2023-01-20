############################################
## Predicting genes involvment in disease
############################################

#################
##Function
rankingGenes<-function(genesData, phase, runMode, patientInfo){
  
  #########################################################################
  ##  Loading required scoring functions depending on the running Mode
  ## Este source es critico mantenerlo aqui
  #########################################################################
  
  if(runMode == "Standard"){
    source("functions/scoringGenes/standard/Functions_toGet_LOFandGOFscores.R",local = TRUE)
    source("functions/scoringGenes/standard/phenoScore_function.R",local = TRUE)
    
  }else if(runMode == "High-Specificity"){
    source("functions/scoringGenes/high_specificity/Functions_toGet_LOFandGOFscores.R",local = TRUE)
    source("functions/scoringGenes/high_specificity/phenoScore_function.R",local = TRUE)
    
  }
  
  # ##Thresholds de expression para LOF y GOFs
  # source("scripts_To_Load_Data/ExpressionThresholds_ForGofAndLof.R",
  #        local = TRUE)
  
  
  ##############################################
  ###Return matrix with scores per phase
  ##############################################
  
  matScores<-as.data.frame(matrix(data = NA, nrow = nrow(genesData), ncol = 5))
  rownames(matScores)<-rownames(genesData)
  colnames(matScores)<-c("GOF_score","LOF_score","Max_type","Max_Score","type")##score: 0-1 (1 high chances of implicated), ##type: GOF--LOF
  ##Max_Score to hold the max(GOF,LOF)
  ##Max type, to indicate if the most_likely scenario is GOF or LOF (if 0.5 --- 0.5 ) equally_likely
  
  #########################################
  ##Adding info about subscores for main scores, so that we can provide also that info
  
  #GenePhenoScore will be the same, but separate in case this changes in the future
  matScores$genePhenoScore_GOF<-NA
  matScores$genePhenoScore_LOF<-NA
  
  matScores$geneEnhancerScore_GOF<-NA
  matScores$geneEnhancerScore_LOF<-NA
  
  matScores$geneFeaturesScore_GOF<-NA
  matScores$geneFeaturesScore_LOF<-NA
  
  ##For subscores of the geneFeatures score
  matScores$dosageSensitivityScore_GOF<-NA
  matScores$dosageSensitivityScore_LOF<-NA
  
  matScores$polycombScore_GOF<-NA
  matScores$polycombScore_LOF<-NA
  
  matScores$geneExpressionScore_GOF<-NA
  matScores$geneExpressionScore_LOF<-NA
  

  #######################################################################
  ## We provide results for the different considered phases, in the end only 1 phase 
  #######################################################################
  res_scores<-list()
  res_scores[[phase]]<-matScores

  ##In addition, phaseFree
  #res_scores[["phaseFree"]]<-matScores
  
  ################################
  ###Let's fill the matrices
  ################################
  for(gene in rownames(genesData)){
    
    info_gene<-genesData[gene,]
    
    #########################################################
    ## Potential relatioship of the gene with the phenotypes
    ##Indep of the phase
    #########################################################
    
    geneTransversalData<-info_gene[,c("associatedPhenotypeIn_OMIM","associatedPhenotypeIn_MGI",
                                      ###"n_PhenotypeRelated_Through_OMIM_Human","n_PhenotypeRelated_Through_MGI_Mice",
                                      "mainPhenotype_Through_OMIM_Human","mainPhenotype_Through_MGI_Mice",
                                      "nature_HI_score",
                                      "huang_HI_score",
                                      "clinGene_HI_score",
                                      "cell_HI_score",
                                      "cell_TriploSense_score",
                                      "TAU_exp")]
    
    ###Pheno Score, used if required afterwards
    
    dataForPhenoScore<-info_gene[,c("associatedPhenotypeIn_OMIM","associatedPhenotypeIn_MGI",
                                    "mainPhenotype_Through_OMIM_Human","mainPhenotype_Through_MGI_Mice")]
    
    ## phenoScore different for GOF than LOF so computing both, then use the one required depending the prediction
    phenoScore_LOF<-phenoScore_fun(geneData = dataForPhenoScore, gof_case = FALSE, patientInfo)
    phenoScore_GOF<-phenoScore_fun(geneData = dataForPhenoScore, gof_case = TRUE, patientInfo)
    
    ###############################
    ###Deciphering Etiology
    ###############################
    
    ###DOING IT PER CONSIDERED PHASE
    
    ##Creating standarized matrix per phase. MATRIX PHASEs
    phaseMatrixes<-list()
    
    ##Recall phases vector is now just a one element character vector
    if(phase != "phaseFree"){
      
      phaseMatrixes[[phase]]<-info_gene[,c(
        "polyComb_score",
        paste0("FPKM_",phase),
        paste0("nEnhancers_initial_",phase),
        paste0("nEnhancers_kept_",phase),
        paste0("nEnhancers_gained_",phase),
        paste0("nEnhancers_maxAvailableInTheOtherDomain_",phase),
        paste0("enhancers_acetilation_initial_",phase),
        paste0("enhancers_acetilation_kept_",phase),
        paste0("enhancers_acetilation_gained_",phase),
        paste0("enhancers_maxAcetilationAvailableInTheOtherDomain_",phase))]
      
      ######################
      ##Standarize colnames
      #Since we are gonna deal with standarized matrixes per phase
      for(n in 1:length(phaseMatrixes)){
        colnames(phaseMatrixes[[n]])<-c("polyComb_score",
                                        "FPKM",
                                        "nEnhancers_initial","nEnhancers_kept","nEnhancers_gained",
                                        "nEnhancers_maxAvailableInTheOtherDomain",
                                        "enhancers_acetilation_initial","enhancers_acetilation_kept",
                                        "enhancers_acetilation_gained",
                                        "enhancers_maxAcetilationAvailableInTheOtherDomain")
      }
    }
    
    ############################################
    ## 1st check
    ##DIRECT OR INDIRECT EFFECT
    
    ##Recall the relevance of the hierarchy when associating genes with regulatory mechanism. Top-down assignment.
    
    if(info_gene[,"RegulatoryMechanism"]=="LongRange"){
      ##Gene intact, so possible INDIRECT EFFECT
      ####
      ##Meter una sentence para el gen,del tipo, no aparece duplicado ni deleccionado, secuencia intacta
      ##Por ello el mecanismo sera indirecto
      
      ####################################
      ## Checking per phase
      #####################################
      for(phase in names(res_scores)){
        if(phase != "phaseFree"){
          #Hence for phases where expression data and enhancer data is considered
          
          matrixPhase<-phaseMatrixes[[phase]]
          
          # 2nd check
          # Expressed or Not Expressed
          
          ###Evaluate if LOF or GOF is likely
          
          ## LOF evaluation: LONG-RANGE EFFECT
          lof_score_metadata<-eval_lof_indirectEffect_score(matrixPhase = matrixPhase, 
                                                            geneTransversalData = geneTransversalData,
                                                            phenoScore = phenoScore_LOF,
                                                            threshold_MaxExpresion = threshold_MaxExpresion,
                                                            threshold_MinExpresion = threshold_MinExpresion,
                                                            minRatioEnhBalance = minRatioEnhBalance,
                                                            maxRatioEnhBalance = maxRatioEnhBalance_LOF)
          lof_score<-lof_score_metadata$finalScore
          
          ## GOF evaluation: LONG-RANGE EFFECT
          gof_score_metadata<-eval_Gof_indirectEffect_score(matrixPhase = matrixPhase, 
                                                            geneTransversalData = geneTransversalData,
                                                            phenoScore = phenoScore_GOF,
                                                            threshold_MaxExpresion = threshold_MaxExpresion,
                                                            threshold_MinExpresion = threshold_MinExpresion,
                                                            minRatioEnhBalance = minRatioEnhBalance,
                                                            maxRatioEnhBalance = maxRatioEnhBalance_GOF)
          gof_score<-gof_score_metadata$finalScore
          
          ##Get max score && type
          Scores<-c(lof_score, gof_score)
          names(Scores)<-c("LOF_score","GOF_score")
          
          max_score<-max(Scores)## (unless both equal check, like 0,0 or 0.8 0.8)
          
          if(lof_score==gof_score){
            typeMax<-"equally_likely"
          }else{
            ##one score bigger than the other
            typeMax<-names(Scores)[which(Scores==max_score)]
          }
          
          
          ##if there is no expression data for the gene, we do not compute anything
          ##We put a -1 to place it at the bottom
          if(matrixPhase[,"FPKM"]==-1){
            max_score<-0
            typeMax<-"equally_likely"
            
            gof_score<-0
            lof_score<-0
            
          }
          
          ##if there multiple expression values, distant from them (eg 0.5 and 2 FPKM)
          ##We put a -2 to place it at the bottom
          if(matrixPhase[,"FPKM"]==-2){
            max_score<-0
            typeMax<-"equally_likely"
            
            gof_score<-0
            lof_score<-0
            
          }
          
          ##Add results
          # colnames(matScores)<-c("GOF_score","LOF_score","Max_type","Max_Score","type")
          res_scores[[phase]][gene,"GOF_score"]<-gof_score
          res_scores[[phase]][gene,"LOF_score"]<-lof_score
          
          res_scores[[phase]][gene,"Max_type"]<-typeMax       
          res_scores[[phase]][gene,"Max_Score"]<-max_score 
          
          ##type mechanism: longrange...
          res_scores[[phase]][gene,"type"]<-"LongRange"
          
          ##RECORDING METADATA/Additional info pathogenic score
          res_scores[[phase]][gene,"genePhenoScore_GOF"]<-gof_score_metadata$genePhenoScore
          res_scores[[phase]][gene,"genePhenoScore_LOF"]<-lof_score_metadata$genePhenoScore
          
          res_scores[[phase]][gene,"geneEnhancerScore_GOF"]<-gof_score_metadata$geneEnhancerScore
          res_scores[[phase]][gene,"geneEnhancerScore_LOF"]<-lof_score_metadata$geneEnhancerScore
          
          res_scores[[phase]][gene,"geneFeaturesScore_GOF"]<-gof_score_metadata$geneFeaturesScore
          res_scores[[phase]][gene,"geneFeaturesScore_LOF"]<-lof_score_metadata$geneFeaturesScore
          
          res_scores[[phase]][gene,"dosageSensitivityScore_GOF"]<-gof_score_metadata$dosageSensitivityScore
          res_scores[[phase]][gene,"dosageSensitivityScore_LOF"]<-lof_score_metadata$dosageSensitivityScore
          
          res_scores[[phase]][gene,"polycombScore_GOF"]<-gof_score_metadata$polycombScore
          res_scores[[phase]][gene,"polycombScore_LOF"]<-lof_score_metadata$polycombScore
          
          res_scores[[phase]][gene,"geneExpressionScore_GOF"]<-gof_score_metadata$geneExpressionScore
          res_scores[[phase]][gene,"geneExpressionScore_LOF"]<-lof_score_metadata$geneExpressionScore
          
        }
      }
      
    }else if(info_gene[,"RegulatoryMechanism"]=="Direct_geneTruncation"){
      
      ########################################
      ## CASE FOR GENES TRUNCATED ############
      ########################################
      ##In the future, this one will be adapted to handle fusion proteins
      
      
      ##So the gene is broken
      ##score GOF 0
      ##score LOF 
      
      ####################################
      ## Checking per phase
      #####################################
      for(phase in names(res_scores)){
        if(phase != "phaseFree"){
          # print(phase)
          matrixPhase<-phaseMatrixes[[phase]]
          
          ##score GOF 0
          ##score LOF 
          
          ## LOF evaluation: DIRECT EFFECT
          lof_score_metadata<-eval_lof_directEffect_score(matrixPhase = matrixPhase, 
                                                          geneTransversalData = geneTransversalData,
                                                          phenoScore = phenoScore_LOF,
                                                          threshold_MaxExpresion = threshold_MaxExpresion,
                                                          threshold_MinExpresion = threshold_MinExpresion)
          lof_score<-lof_score_metadata$finalScore
          
          ## GOF evaluation
          # gof_score<-0
          gof_score_metadata<-list("finalScore"=0,
                                   "geneEnhancerScore"=NA,
                                   "genePhenoScore"=NA,
                                   "geneFeaturesScore"=NA,
                                   "dosageSensitivityScore"=NA,
                                   "polycombScore"=NA,
                                   "geneExpressionScore"=NA
          )
          
          gof_score<-gof_score_metadata$finalScore
          
          
          ##Get max score && type
          Scores<-c(lof_score, gof_score)
          names(Scores)<-c("LOF_score","GOF_score")
          
          max_score<-max(Scores)## (unless both equal check, like 0,0 or 0.8 0.8)
          
          if(lof_score==gof_score){
            typeMax<-"equally_likely"
          }else{
            ##one score bigger than the other
            typeMax<-names(Scores)[which(Scores==max_score)]
          }
          
          ##if there is no expression data for the gene, we do not compute anything
          ##We put a -1 to place it at the bottom
          if(matrixPhase[,"FPKM"]==-1){
            max_score<-0
            typeMax<-"equally_likely"
            
            gof_score<-0
            lof_score<-0
            
          }
          
          ##if there multiple expression values, distant from them (eg 0.5 and 2 FPKM)
          ##We put a -2 to place it at the bottom
          if(matrixPhase[,"FPKM"]==-2){
            max_score<-0
            typeMax<-"equally_likely"
            
            gof_score<-0
            lof_score<-0
            
          }
          
          
          ##Add results
          # colnames(matScores)<-c("GOF_score","LOF_score","Max_type","Max_Score","type")
          res_scores[[phase]][gene,"GOF_score"]<-gof_score
          res_scores[[phase]][gene,"LOF_score"]<-lof_score
          
          res_scores[[phase]][gene,"Max_type"]<-typeMax       
          res_scores[[phase]][gene,"Max_Score"]<-max_score 
          
          ##type mechanism:
          res_scores[[phase]][gene,"type"]<-"Direct_geneTruncation"
          
          ##RECORDING METADATA/Additional info pathogenic score
          res_scores[[phase]][gene,"genePhenoScore_GOF"]<-gof_score_metadata$genePhenoScore
          res_scores[[phase]][gene,"genePhenoScore_LOF"]<-lof_score_metadata$genePhenoScore
          
          res_scores[[phase]][gene,"geneEnhancerScore_GOF"]<-gof_score_metadata$geneEnhancerScore
          res_scores[[phase]][gene,"geneEnhancerScore_LOF"]<-lof_score_metadata$geneEnhancerScore
          
          res_scores[[phase]][gene,"geneFeaturesScore_GOF"]<-gof_score_metadata$geneFeaturesScore
          res_scores[[phase]][gene,"geneFeaturesScore_LOF"]<-lof_score_metadata$geneFeaturesScore
          
          res_scores[[phase]][gene,"dosageSensitivityScore_GOF"]<-gof_score_metadata$dosageSensitivityScore
          res_scores[[phase]][gene,"dosageSensitivityScore_LOF"]<-lof_score_metadata$dosageSensitivityScore
          
          res_scores[[phase]][gene,"polycombScore_GOF"]<-gof_score_metadata$polycombScore
          res_scores[[phase]][gene,"polycombScore_LOF"]<-lof_score_metadata$polycombScore
          
          res_scores[[phase]][gene,"geneExpressionScore_GOF"]<-gof_score_metadata$geneExpressionScore
          res_scores[[phase]][gene,"geneExpressionScore_LOF"]<-lof_score_metadata$geneExpressionScore
          
        }
      }
      
    }else if(info_gene[,"RegulatoryMechanism"]=="Direct_geneDeletion"){
      
      ########################################
      ## CASE FOR GENES DELETED ############
      ########################################
      
      ##So the gene is deleted
      ##score GOF 0. This gene can not suffer a gain of function. Makes no sense
      ##score LOF 
      
      ####################################
      ## Checking per phase
      #####################################
      for(phase in names(res_scores)){
        if(phase != "phaseFree"){
          # print(phase)
          matrixPhase<-phaseMatrixes[[phase]]
          
          ##score GOF 
          ##score LOF 
          
          ## LOF evaluation: DIRECT EFFECT
          lof_score_metadata<-eval_lof_directEffect_score(matrixPhase = matrixPhase, 
                                                          geneTransversalData = geneTransversalData,
                                                          phenoScore = phenoScore_LOF,
                                                          threshold_MaxExpresion = threshold_MaxExpresion,
                                                          threshold_MinExpresion = threshold_MinExpresion)
          lof_score<-lof_score_metadata$finalScore
          
          ## GOF evaluation
          # gof_score<-0
          
          gof_score_metadata<-list("finalScore"=0,
                                   "geneEnhancerScore"=NA,
                                   "genePhenoScore"=NA,
                                   "geneFeaturesScore"=NA,
                                   "dosageSensitivityScore"=NA,
                                   "polycombScore"=NA,
                                   "geneExpressionScore"=NA
          )
          
          gof_score<-gof_score_metadata$finalScore
          
          
          
          ##Get max score && type
          Scores<-c(lof_score, gof_score)
          names(Scores)<-c("LOF_score","GOF_score")
          
          max_score<-max(Scores)## (unless both equal check, like 0,0 or 0.8 0.8)
          
          if(lof_score==gof_score){
            typeMax<-"equally_likely"
          }else{
            ##one score bigger than the other
            typeMax<-names(Scores)[which(Scores==max_score)]
          }
          
          ##if there is no expression data for the gene, we do not compute anything
          ##We put a -1 to place it at the bottom
          if(matrixPhase[,"FPKM"]==-1){
            max_score<-0
            typeMax<-"equally_likely"
            
            gof_score<-0
            lof_score<-0
            
          }
          
          ##if there multiple expression values, distant from them (eg 0.5 and 2 FPKM)
          ##We put a -2 to place it at the bottom
          if(matrixPhase[,"FPKM"]==-2){
            max_score<-0
            typeMax<-"equally_likely"
            
            gof_score<-0
            lof_score<-0
            
          }
          
          
          ##Add results
          # colnames(matScores)<-c("GOF_score","LOF_score","Max_type","Max_Score","type")
          res_scores[[phase]][gene,"GOF_score"]<-gof_score
          res_scores[[phase]][gene,"LOF_score"]<-lof_score
          
          res_scores[[phase]][gene,"Max_type"]<-typeMax       
          res_scores[[phase]][gene,"Max_Score"]<-max_score 
          
          ##type mechanism:
          res_scores[[phase]][gene,"type"]<-"Direct_geneDeletion"
          
          ##RECORDING METADATA/Additional info pathogenic score
          res_scores[[phase]][gene,"genePhenoScore_GOF"]<-gof_score_metadata$genePhenoScore
          res_scores[[phase]][gene,"genePhenoScore_LOF"]<-lof_score_metadata$genePhenoScore
          
          res_scores[[phase]][gene,"geneEnhancerScore_GOF"]<-gof_score_metadata$geneEnhancerScore
          res_scores[[phase]][gene,"geneEnhancerScore_LOF"]<-lof_score_metadata$geneEnhancerScore
          
          res_scores[[phase]][gene,"geneFeaturesScore_GOF"]<-gof_score_metadata$geneFeaturesScore
          res_scores[[phase]][gene,"geneFeaturesScore_LOF"]<-lof_score_metadata$geneFeaturesScore
          
          res_scores[[phase]][gene,"dosageSensitivityScore_GOF"]<-gof_score_metadata$dosageSensitivityScore
          res_scores[[phase]][gene,"dosageSensitivityScore_LOF"]<-lof_score_metadata$dosageSensitivityScore
          
          res_scores[[phase]][gene,"polycombScore_GOF"]<-gof_score_metadata$polycombScore
          res_scores[[phase]][gene,"polycombScore_LOF"]<-lof_score_metadata$polycombScore
          
          res_scores[[phase]][gene,"geneExpressionScore_GOF"]<-gof_score_metadata$geneExpressionScore
          res_scores[[phase]][gene,"geneExpressionScore_LOF"]<-lof_score_metadata$geneExpressionScore
          
        }
      }
      
    }else if(info_gene[,"RegulatoryMechanism"]=="Direct_LongRange_geneDuplication"){
      
      ########################################
      ## CASE FOR DUPLICATED GENES ###########
      ########################################
      
      ##So the gene is duplicated
      ##In this case, the most tricky one
      ##The gene can have an upregulation by enh adoption (neoTad) or enh duplication
      ##Or because it is currently expressed and its expression gets boosted by its duplication
      ##score GOF 
      ##score LOF 0.This gene can not suffer a loss of function. UNLESS it is broken, but in that case its regulatory mechanism is Direct_geneTruncation
      
      ####################################
      ## Checking per phase
      #####################################
      for(phase in names(res_scores)){
        if(phase != "phaseFree"){
          # print(phase)
          matrixPhase<-phaseMatrixes[[phase]]
          
          ##score GOF 
          ##score LOF 
          
          ## LOF evaluation:
          # lof_score<-0
          
          lof_score_metadata<-list("finalScore"=0,
                                   "geneEnhancerScore"=NA,
                                   "genePhenoScore"=NA,
                                   "geneFeaturesScore"=NA,
                                   "dosageSensitivityScore"=NA,
                                   "polycombScore"=NA,
                                   "geneExpressionScore"=NA
          )
          
          lof_score<-lof_score_metadata$finalScore
          
          
          ## GOF evaluation TWO OPTIONS
          ##By LongRange.Which takes into account poor expressed&enh adoption
          ##Lacking enh duplication...but attending to shadow enh...maybe not very convenient
          ##By directEffect. Which takes into account if it is considerably expressed and its duplicated, this one ignores enhancers
          
          ##Long range only computed for the genes whose TAD is disrupted (and not for those whose TADs is entirely duplicated)
          ##If gene intact but enh duplicated, it will go to the LongRange geneMech hence to other scoring category
          ##For those whose TAD is entirely duplicated we are not even considering where the enhancers are located
          
          if(info_gene[,"TypeDomainInitial"]=="TAD_disrupted"){
            ##It will take into account nEnh before and after, as with LongRange translocation, inversion, or deletion
            gof_score_metadata_longRange<-eval_Gof_indirectEffect_score(matrixPhase = matrixPhase, 
                                                                        geneTransversalData = geneTransversalData,
                                                                        phenoScore = phenoScore_GOF,
                                                                        threshold_MaxExpresion = threshold_MaxExpresion,
                                                                        threshold_MinExpresion = threshold_MinExpresion,
                                                                        minRatioEnhBalance = minRatioEnhBalance,
                                                                        maxRatioEnhBalance = maxRatioEnhBalance_GOF)
            
            gof_score_longRange<-gof_score_metadata_longRange$finalScore
            
          }else{
            ##We assign this score to 0, so that this can never be higher than the one by direct effect
            # gof_score_longRange<-0
            ##We assign this score to 0, so that this can never be higher than the one by direct effect
            gof_score_metadata_longRange<-list("finalScore"=0,
                                               "geneEnhancerScore"=NA,
                                               "genePhenoScore"=NA,
                                               "geneFeaturesScore"=NA,
                                               "dosageSensitivityScore"=NA,
                                               "polycombScore"=NA,
                                               "geneExpressionScore"=NA
            )
            gof_score_longRange<-gof_score_metadata_longRange$finalScore
            
          }
          
          ###SEGUIR POR AQUI, HE DE METER EL DE _METADATA Y LUEGO PILLAR EL VALOR DEL FINAL SCORE
          gof_score_metadata_directEffect<-eval_Gof_directEffect_score(matrixPhase = matrixPhase, 
                                                                       geneTransversalData = geneTransversalData,
                                                                       phenoScore = phenoScore_GOF,
                                                                       threshold_MaxExpresion = threshold_MaxExpresion,
                                                                       threshold_MinExpresion = threshold_MinExpresion,
                                                                       minRatioEnhBalance = minRatioEnhBalance,
                                                                       maxRatioEnhBalance = maxRatioEnhBalance_GOF)
            
          gof_score_directEffect<-gof_score_metadata_directEffect$finalScore
          
          ##Then we retrieve the maximum score
          ##Afterwards, depending on which one bigger longRange_GeneDuplication or Direct_GeneDuplication
          gof_score<-max(gof_score_longRange,gof_score_directEffect)
          
          ##Get max score && type
          Scores<-c(lof_score, gof_score)
          names(Scores)<-c("LOF_score","GOF_score")
          
          max_score<-max(Scores)## (unless both equal check, like 0,0 or 0.8 0.8)
          
          if(lof_score==gof_score){
            typeMax<-"equally_likely"
          }else{
            ##one score bigger than the other
            typeMax<-names(Scores)[which(Scores==max_score)]
          }
          
          ##if there is no expression data for the gene, we do not compute anything
          ##We put a -1 to place it at the bottom
          if(matrixPhase[,"FPKM"]==-1){
            max_score<-0
            typeMax<-"equally_likely"
            
            gof_score<-0
            lof_score<-0
            
          }
          
          ##if there multiple expression values, distant from them (eg 0.5 and 2 FPKM)
          ##We put a -2 to place it at the bottom
          if(matrixPhase[,"FPKM"]==-2){
            max_score<-0
            typeMax<-"equally_likely"
            
            gof_score<-0
            lof_score<-0
            
          }
          
          
          ##Add results
          # colnames(matScores)<-c("GOF_score","LOF_score","Max_type","Max_Score","type")
          res_scores[[phase]][gene,"GOF_score"]<-gof_score
          res_scores[[phase]][gene,"LOF_score"]<-lof_score
          
          res_scores[[phase]][gene,"Max_type"]<-typeMax       
          res_scores[[phase]][gene,"Max_Score"]<-max_score 
          
          
          ##type mechanism:
          ##Narrow depending if longrange>direct impact
          #gof_score<-max(gof_score_longRange,gof_score_directEffect)
          if(gof_score_longRange > gof_score_directEffect ){
            res_scores[[phase]][gene,"type"]<-"LongRange_geneDuplication" ##This means NeoTAD
            
            ##So, adding long-range info
            res_scores[[phase]][gene,"genePhenoScore_GOF"]<-gof_score_metadata_longRange$genePhenoScore

            res_scores[[phase]][gene,"geneEnhancerScore_GOF"]<-gof_score_metadata_longRange$geneEnhancerScore

            res_scores[[phase]][gene,"geneFeaturesScore_GOF"]<-gof_score_metadata_longRange$geneFeaturesScore

            res_scores[[phase]][gene,"dosageSensitivityScore_GOF"]<-gof_score_metadata_longRange$dosageSensitivityScore

            res_scores[[phase]][gene,"polycombScore_GOF"]<-gof_score_metadata_longRange$polycombScore

            res_scores[[phase]][gene,"geneExpressionScore_GOF"]<-gof_score_metadata_longRange$geneExpressionScore

          }else if(gof_score_directEffect > gof_score_longRange){
            res_scores[[phase]][gene,"type"]<-"Direct_geneDuplication"
            
            ##Adding direct-effect info
            res_scores[[phase]][gene,"genePhenoScore_GOF"]<-gof_score_metadata_directEffect$genePhenoScore
            
            res_scores[[phase]][gene,"geneEnhancerScore_GOF"]<-gof_score_metadata_directEffect$geneEnhancerScore
            
            res_scores[[phase]][gene,"geneFeaturesScore_GOF"]<-gof_score_metadata_directEffect$geneFeaturesScore
            
            res_scores[[phase]][gene,"dosageSensitivityScore_GOF"]<-gof_score_metadata_directEffect$dosageSensitivityScore
            
            res_scores[[phase]][gene,"polycombScore_GOF"]<-gof_score_metadata_directEffect$polycombScore
            
            res_scores[[phase]][gene,"geneExpressionScore_GOF"]<-gof_score_metadata_directEffect$geneExpressionScore
            
          }else if(gof_score_directEffect == gof_score_longRange){
            res_scores[[phase]][gene,"type"]<-"Direct_LongRange_geneDuplication"
            
            #Since long range also has a high score adding long-range info
            res_scores[[phase]][gene,"genePhenoScore_GOF"]<-gof_score_metadata_longRange$genePhenoScore
            
            res_scores[[phase]][gene,"geneEnhancerScore_GOF"]<-gof_score_metadata_longRange$geneEnhancerScore
            
            res_scores[[phase]][gene,"geneFeaturesScore_GOF"]<-gof_score_metadata_longRange$geneFeaturesScore
            
            res_scores[[phase]][gene,"dosageSensitivityScore_GOF"]<-gof_score_metadata_longRange$dosageSensitivityScore
            
            res_scores[[phase]][gene,"polycombScore_GOF"]<-gof_score_metadata_longRange$polycombScore
            
            res_scores[[phase]][gene,"geneExpressionScore_GOF"]<-gof_score_metadata_longRange$geneExpressionScore
          }
          
          ##Only 1 lofscore computed
          ##RECORDING METADATA/Additional info pathogenic score
          res_scores[[phase]][gene,"genePhenoScore_LOF"]<-lof_score_metadata$genePhenoScore
          
          res_scores[[phase]][gene,"geneEnhancerScore_LOF"]<-lof_score_metadata$geneEnhancerScore
          
          res_scores[[phase]][gene,"geneFeaturesScore_LOF"]<-lof_score_metadata$geneFeaturesScore
          
          res_scores[[phase]][gene,"dosageSensitivityScore_LOF"]<-lof_score_metadata$dosageSensitivityScore
          
          res_scores[[phase]][gene,"polycombScore_LOF"]<-lof_score_metadata$polycombScore
          
          res_scores[[phase]][gene,"geneExpressionScore_LOF"]<-lof_score_metadata$geneExpressionScore
          
        }
      }
      
    }else if(info_gene[,"RegulatoryMechanism"]=="Direct_uncertaintyRegion"){
      
      ###########################################
      ## FOR GENES IN THE UNCERTAINTY REGION ####
      ## Por si luego queremos modelarlo de forma diferente por separado
      ##################################################################
      
      
      ##As of now, model it as if it was a truncated gen. Since it is in an area that has been 
      ##affected by a rearrangement.
      
      ##So the gene is in the UNCERTAINTY REGION
      ##score GOF 0
      ##score LOF 
      
      ####################################
      ## Checking per phase
      #####################################
      for(phase in names(res_scores)){
        if(phase != "phaseFree"){
          matrixPhase<-phaseMatrixes[[phase]]
          
          ##score GOF 0
          ##score LOF 
          
          ##Lo tratamos como si estuviera broken
          ##como los direct effects
          ## LOF evaluation: DIRECT EFFECT
          lof_score_metadata<-eval_lof_directEffect_score(matrixPhase = matrixPhase, 
                                                          geneTransversalData = geneTransversalData,
                                                          phenoScore = phenoScore_LOF,
                                                          threshold_MaxExpresion = threshold_MaxExpresion,
                                                          threshold_MinExpresion = threshold_MinExpresion)
          lof_score<-lof_score_metadata$finalScore
          
          ## GOF evaluation
          # gof_score<-0
          
          gof_score_metadata<-list("finalScore"=0,
                                   "geneEnhancerScore"=NA,
                                   "genePhenoScore"=NA,
                                   "geneFeaturesScore"=NA,
                                   "dosageSensitivityScore"=NA,
                                   "polycombScore"=NA,
                                   "geneExpressionScore"=NA
          )
          
          gof_score<-gof_score_metadata$finalScore
          
          
          
          ##Get max score && type
          Scores<-c(lof_score, gof_score)
          names(Scores)<-c("LOF_score","GOF_score")
          
          max_score<-max(Scores)## (unless both equal check, like 0,0 or 0.8 0.8)
          
          if(lof_score==gof_score){
            typeMax<-"equally_likely"
          }else{
            ##one score bigger than the other
            typeMax<-names(Scores)[which(Scores==max_score)]
          }
          
          ##if there is no expression data for the gene, we do not compute anything
          ##We put a -1 to place it at the bottom
          if(matrixPhase[,"FPKM"]==-1){
            max_score<-0
            typeMax<-"equally_likely"
            
            gof_score<-0
            lof_score<-0
            
          }
          
          ##if there multiple expression values, distant from them (eg 0.5 and 2 FPKM)
          ##We put a -2 to place it at the bottom
          if(matrixPhase[,"FPKM"]==-2){
            max_score<-0
            typeMax<-"equally_likely"
            
            gof_score<-0
            lof_score<-0
            
          }
          
          ##Add results
          # colnames(matScores)<-c("GOF_score","LOF_score","Max_type","Max_Score","type")
          res_scores[[phase]][gene,"GOF_score"]<-gof_score
          res_scores[[phase]][gene,"LOF_score"]<-lof_score
          
          res_scores[[phase]][gene,"Max_type"]<-typeMax       
          res_scores[[phase]][gene,"Max_Score"]<-max_score 
          
          ##type mechanism:
          res_scores[[phase]][gene,"type"]<-"Direct_uncertaintyRegion"
          
          ##RECORDING METADATA/Additional info pathogenic score
          res_scores[[phase]][gene,"genePhenoScore_GOF"]<-gof_score_metadata$genePhenoScore
          res_scores[[phase]][gene,"genePhenoScore_LOF"]<-lof_score_metadata$genePhenoScore
          
          res_scores[[phase]][gene,"geneEnhancerScore_GOF"]<-gof_score_metadata$geneEnhancerScore
          res_scores[[phase]][gene,"geneEnhancerScore_LOF"]<-lof_score_metadata$geneEnhancerScore
          
          res_scores[[phase]][gene,"geneFeaturesScore_GOF"]<-gof_score_metadata$geneFeaturesScore
          res_scores[[phase]][gene,"geneFeaturesScore_LOF"]<-lof_score_metadata$geneFeaturesScore
          
          res_scores[[phase]][gene,"dosageSensitivityScore_GOF"]<-gof_score_metadata$dosageSensitivityScore
          res_scores[[phase]][gene,"dosageSensitivityScore_LOF"]<-lof_score_metadata$dosageSensitivityScore
          
          res_scores[[phase]][gene,"polycombScore_GOF"]<-gof_score_metadata$polycombScore
          res_scores[[phase]][gene,"polycombScore_LOF"]<-lof_score_metadata$polycombScore
          
          res_scores[[phase]][gene,"geneExpressionScore_GOF"]<-gof_score_metadata$geneExpressionScore
          res_scores[[phase]][gene,"geneExpressionScore_LOF"]<-lof_score_metadata$geneExpressionScore
          
        }
      }
    }
    
    ### Update regulatory mechanism column of MatrixesGeneEvaluation
    ### Specially for those cases where a priori not clear what could be happening
    #######################
    ## Initially if a gene is duplicated it can be upregulated by LongRange or direct effect.
    ## That is why it appear on Matrixes Gene Evaluation as: Direct_LongRange_geneDuplication
    ## Upon scoring the prediction can pivot to Direct_geneDuplication or LongRange_geneDuplication. When this occur. Let's modify
    ## the evaluating gene situation matrixes because downstream is used for the report this column. 
    
    ##We add the one chosen for the res_scores
    genesData[gene,"RegulatoryMechanism"]<-res_scores[[phase]][gene,"type"]
    
  }
  
  
  #######################################################################
  ###If length GOF long range > 1. Re rank by proximity with enhancers
  ##Pendent.
  
  resultsPrediction<-list("MatrixesGeneEvaluation"= genesData, ##updated
                          "ScoresResults" = res_scores)
  

  return(resultsPrediction)
  
}
