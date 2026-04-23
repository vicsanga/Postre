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

  ##In addition, CellTypeAgnostic
  #res_scores[["CellTypeAgnostic"]]<-matScores
  
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
    phenoScore_LOF<-phenoScore_fun(geneData = dataForPhenoScore, gof_case = FALSE, patientInfo = patientInfo)
    phenoScore_GOF<-phenoScore_fun(geneData = dataForPhenoScore, gof_case = TRUE, patientInfo = patientInfo)
    
    ###############################
    ###Deciphering Etiology
    ###############################
    
    ###DOING IT PER CONSIDERED PHASE
    
    ##Creating standarized matrix per phase. MATRIX PHASEs
    phaseMatrixes<-list()
    
    ##Recall phases vector is now just a one element character vector
    #The following matrixes were previously created only for phase != "CellTypeAgnostic"
    #But as of 29 april 2025, to keep logic, and leave door open for upgrading functionality
    #Generating the matrix for them too now. BUT, as of now these enhancers info will be just
    #EMPTY CONTENT, we have it since it comes from the info gene data
    #AND IT WILL SIMPLY BE IGNORED DOWNSTREAM FOR NOW
    ##IT IS ALL EQUAL TO 0 AND FPKM -1
    
    # if(phase != "CellTypeAgnostic"){}
    #Lo siguiente estaba antes anidado dentro del if previo
    
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
        if(phase != "CellTypeAgnostic"){
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
          
        }else{
          ## CellTypeAgnostic (CELL TYPE-AGNOSTIC SCENARIO)
          ## HERE, THE CODE IS LOCATED IN THE CONTEXT OF LONG-RANGE EVALUATION
          ##Current version as of april 2025, setting scores to 0 no prediction performed
          #But possibility to extend in the future.
          
          # Deprecated for now. If gene Intact, and no enh info... let's not assume long range effect without evidence
          # That would be kind of a prediction of a prediction (because long range you can not be sure neither). For now let's avoid that
          ##BUT I SHOULD ADD INFORMATION AS EMPTY/SET TO O OR STH, TO ENSURE THAT THE DOWNSTREAM PROCESSING 
          ##DOES NOT PRESENT ERRORS 
          ##TO DO THAT, GONNA USE PHASE FREE SCORING FUNCTIONS, EVEN IF NO LOGIC (EMPTY) INSIDE
          
          ##AL SER CellTypeAgnostic matrixPhase VA A ESTAR VACIA PORQUE NO HE ALMACENADO NADA
          #SE NECESITA DOWNSTREAM (no tenerla genera errores), PERO PARA LONG-RANGE VA todo a 0
          matrixPhase<-phaseMatrixes[[phase]]
          
          ###Evaluate if LOF or GOF is likely
          
          ## LOF evaluation: LONG-RANGE EFFECT
          lof_score_metadata<-eval_lof_indirectEffect_CellTypeAgnostic_score(matrixPhase = matrixPhase, 
                                                                      geneTransversalData = geneTransversalData,
                                                                      phenoScore = phenoScore_LOF)
          lof_score<-lof_score_metadata$finalScore
          
          ## GOF evaluation: LONG-RANGE EFFECT
          gof_score_metadata<-eval_Gof_indirectEffect_CellTypeAgnostic_score(matrixPhase = matrixPhase, 
                                                                      geneTransversalData = geneTransversalData,
                                                                      phenoScore = phenoScore_GOF)
          gof_score<-gof_score_metadata$finalScore
        }
        
        ## RECORDING RESULTS
        res_scores[[phase]][gene,"type"]<-"LongRange"
        ##Running shared code to minimize repetition (put after ifs, once properly propagated)
        source("functions/scoringGenes/auxFunctions/sharedCode.R",local = TRUE)
        
      }
      
    }else if(info_gene[,"RegulatoryMechanism"]=="Direct_geneTruncation"){
      
      ########################################
      ## CASE FOR GENES TRUNCATED ############
      ########################################
      ##So the gene is broken (the SV falls between its TSS and TTS)
      ##In the future, this one will be adapted to handle fusion proteins
      ##FURTHERMORE!!
      #It could still be a non-coding deletion in the middle of an intron, but in the current version assumed to be potentially truncating
      #E.g. could bee a deep intronic variant altering splicing (Leave as potential improvement for the future)
      ##score GOF 0
      ##score LOF 
      ##Also there could be an enhancer duplication in intronic region... that upregulate the gene with the enh in the intron...
      #IN THIS CATEGORY THERE IS DEFINITELY ROOM FOR IMPROVEMENT
      
      ####################################
      ## Checking per phase
      #####################################
      for(phase in names(res_scores)){
        if(phase != "CellTypeAgnostic"){
          
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
          
        }else{
          ## CellTypeAgnostic scenario
          
          
          matrixPhase<-phaseMatrixes[[phase]]
          
          ##score GOF 0
          ##score LOF 
          
          ## LOF evaluation: DIRECT EFFECT
          lof_score_metadata<-eval_lof_directEffect_CellTypeAgnostic_score(matrixPhase = matrixPhase, 
                                                                    geneTransversalData = geneTransversalData,
                                                                    phenoScore = phenoScore_LOF)
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
          
        }
        ## RECORDING RESULTS
        res_scores[[phase]][gene,"type"]<-"Direct_geneTruncation"
        ##Running shared code to minimize repetition (put after ifs, once properly propagated)
        source("functions/scoringGenes/auxFunctions/sharedCode.R",local = TRUE)
        
      }
      
    }else if(info_gene[,"RegulatoryMechanism"]=="Direct_geneDeletion"){
      
      ########################################
      ## CASE FOR GENES DELETED ############
      ########################################
      
      ##So the gene is deleted
      ##score GOF 0. Here the gene can not suffer a gain of function. NO DOUBT. Makes no sense
      ##score LOF 
      
      ####################################
      ## Checking per phase
      #####################################
      for(phase in names(res_scores)){
        if(phase != "CellTypeAgnostic"){
          
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
          
        }else{
          ##CellTypeAgnostic
          matrixPhase<-phaseMatrixes[[phase]]
          
          ## LOF evaluation: DIRECT EFFECT
          lof_score_metadata<-eval_lof_directEffect_CellTypeAgnostic_score(matrixPhase = matrixPhase, 
                                                                    geneTransversalData = geneTransversalData,
                                                                    phenoScore = phenoScore_LOF)
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
          
        }
        
        ## RECORDING RESULTS
        res_scores[[phase]][gene,"type"]<-"Direct_geneDeletion"
        ##Running shared code to minimize repetition (put after ifs, once properly propagated)
        source("functions/scoringGenes/auxFunctions/sharedCode.R",local = TRUE)
        
        
      }
      
    }else if(info_gene[,"RegulatoryMechanism"]=="Direct_LongRange_geneDuplication"){
      
      ##NOTA VICTOR 3 JULIO, aqui se me va a quedar mas codigo duplicado, veo mas complejo gestionarlo
      #Ya le dare una vuelta, de momento voy a simplemente desplegar la logica
      
      ########################################
      ## CASE FOR DUPLICATED GENES ###########
      ########################################
      
      ##So the gene is duplicated
      ##In this case, the most tricky one
      ##The gene can have an upregulation by enh adoption (neoTad) or enh duplication
      ##Or because it is currently expressed and its expression gets boosted by its duplication
      ##score GOF 
      ##score LOF 0.This gene can not suffer a loss of function. UNLESS it is broken, but in that case its regulatory mechanism is Direct_geneTruncation
      ##Although the latter can be better adressed (but at least a gene expression alteration will be assessed and fished)
      
      ####################################
      ## Checking per phase
      #####################################
      for(phase in names(res_scores)){
        if(phase != "CellTypeAgnostic"){
          
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
          ##By directEffect. Which takes into account if it is considerably expressed and its duplicated, this one ignores enhancers
          
          ##Long range only computed for the genes whose TAD is disrupted (and not for those whose TADs is entirely duplicated)
          ##If gene intact but enh duplicated, it will go to the LongRange geneMech hence to other scoring category
          ##For those whose TAD is entirely duplicated we are not even considering where the enhancers are located
          
          if(info_gene[,"TypeDomainInitial"]=="TAD_disrupted"){
            ###################################################
            ## The gene is duplicated and its TAD is disrupted
            ###################################################
            
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
            ############################################
            ##This means the WHOLE TAD IS DUPLICATED
            ############################################
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
          
          ####Computing coding effect
          gof_score_metadata_directEffect<-eval_Gof_directEffect_score(matrixPhase = matrixPhase, 
                                                                       geneTransversalData = geneTransversalData,
                                                                       phenoScore = phenoScore_GOF,
                                                                       threshold_MaxExpresion = threshold_MaxExpresion,
                                                                       threshold_MinExpresion = threshold_MinExpresion,
                                                                       minRatioEnhBalance = minRatioEnhBalance,
                                                                       maxRatioEnhBalance = maxRatioEnhBalance_GOF)
            
          gof_score_directEffect<-gof_score_metadata_directEffect$finalScore
          
        }else{
          ##PHASE FREE 
          
          matrixPhase<-phaseMatrixes[[phase]]
          
          #Assessing LOF & GOF
          
          ## LOF evaluation: #Set to 0 in gene duplication case as indicated above
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
          ##By directEffect. Which takes into account if it is considerably expressed and its duplicated, this one ignores enhancers
          
          ##Long range only computed for the genes whose TAD is disrupted (and not for those whose TADs is entirely duplicated)
          ##If gene intact but enh duplicated, it will go to the LongRange geneMech hence to other scoring category
          ##For those whose TAD is entirely duplicated we are not even considering where the enhancers are located
          
          ##NOTE 3 JULY 2025, as of now not considering enhancer adoption in this context (can be seen through the function)
          #But just running it to introduce the logic and set the path for future upgrades
          
          ##IN PHASE FREE STILL ONE TAD MAP HAS BEEN MATCHED to the phase. In case is used in a future upgrade
          ##That is why the following situation still has to be adressed
          
          if(info_gene[,"TypeDomainInitial"]=="TAD_disrupted"){
            ###################################################
            ## The gene is duplicated and its TAD is disrupted
            ###################################################
            
            ##It will take into account nEnh before and after, as with LongRange translocation, inversion, or deletion
            gof_score_metadata_longRange<-eval_Gof_indirectEffect_CellTypeAgnostic_score(matrixPhase = matrixPhase, 
                                                                                  geneTransversalData = geneTransversalData,
                                                                                  phenoScore = phenoScore_GOF)
            
            gof_score_longRange<-gof_score_metadata_longRange$finalScore
            
          }else{
            ############################################
            ##This means the WHOLE TAD IS DUPLICATED
            ############################################
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
          
          ####Computing coding effect
          gof_score_metadata_directEffect<-eval_Gof_directEffect_CellTypeAgnostic_score(matrixPhase = matrixPhase, 
                                                                                 geneTransversalData = geneTransversalData,
                                                                                 phenoScore = phenoScore_GOF)
          
          gof_score_directEffect<-gof_score_metadata_directEffect$finalScore
          
        }
        ##Downstream of this point, handling the scenario is a bit por complicated than in the rest of the cases
        #Reason why sharedCode_DupLongRange.R code used
        
        ## RECORDING RESULTS
        #res_scores[[phase]][gene,"type"]<-##HERE THIS IS NOT AS SIMPLE, 3 POSSIBLE LABELS ASSESSED (see shared code)
        ##Running shared code to minimize repetition (put after ifs, once properly propagated)
        #NOT SO sure about puting after ifs, given particularity present here that calls a slightly modified sharedCode
        source("functions/scoringGenes/auxFunctions/sharedCode_DupLongRange.R",local = TRUE)
        
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
        if(phase != "CellTypeAgnostic"){
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
          
        }else{
          ##CellTypeAgnostic scenario
          matrixPhase<-phaseMatrixes[[phase]]
          
          ##score GOF 0
          ##score LOF 
          
          ## LOF evaluation: DIRECT EFFECT
          lof_score_metadata<-eval_lof_directEffect_CellTypeAgnostic_score(matrixPhase = matrixPhase, 
                                                                    geneTransversalData = geneTransversalData,
                                                                    phenoScore = phenoScore_LOF)
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
          
        }
        ## RECORDING RESULTS
        res_scores[[phase]][gene,"type"]<-"Direct_uncertaintyRegion"
        ##Running shared code to minimize repetition (put after ifs, once properly propagated)
        source("functions/scoringGenes/auxFunctions/sharedCode.R",local = TRUE)
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
