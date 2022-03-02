#################################################
## function for prediction of patient etiology
## based on a single TAD map
################################################
patientPrediction_basedOnSingleTADmap<-function(patientInfo,##patient breakpoints && phenotype
                                                mapTads,##TADs map
                                                betweenTads, ##coordinates of the regions between TADs
                                                Master_GeneExpression,
                                                MasterEnh_map,
                                                phase,
                                                phasesVector,
                                                formatedPhenotype,
                                                runMode){
  ######################################
  ##Loading required functions
  ######################################
  source("scripts_To_Load_Data/cargarFunciones.R",local = TRUE)
  source("functions/RankingGenes_Deciphering_etiology.R", local = TRUE)

  ##Function required to load gene expression and enh data
  # source("functions/GeneExpression_AndEnhancerData_Loader.R")
  
  #####################################
  ##Loading required data
  #####################################
  source("scripts_To_Load_Data/cargaDatos.R",local = TRUE)
  

  #phasesVector<-enh_and_exp_data$phasesVector ##Contains the name of the different developmental stages, or phases considered
  #phasesVector<-phase##I think I need to do this

  
  #####################################
  ## Predicting
  #####################################
  
  #################
  ##Checking input
  ##And converting coordinates to hg19 if coordinates in hg38
  patientInfo<-checkPatientFeatures(patientInfo = patientInfo) ##If sth wrong an error will be rised and will stop the program execution
  
  ##################
  ##Define Locations Affected: Either TADs or TAD boundaries
  ##mensajes informativos quitarlos, dejarlos para el report no por consola ahora dando pc
  info_affectedRegions<-affectedRegions(dataPatient = patientInfo, tadsMap = mapTads, regionsBetweenTADs = betweenTads)
  
  ###################
  ##Retrieve Genes in the affected Locations
  info_affectedGenes<-affectedGenes(regionsAltered = info_affectedRegions, gtf = gtf_annotation, onlyProteinCoding = TRUE, patientInfo = patientInfo)
  
  ##############
  ##Retrieve Enhancers in the affected locations
  ##Number of them
  info_affectedEnhancers<-affectedEnhancers(regionsAltered = info_affectedRegions, enhMap = MasterEnh_map)
  ##uncertainty region, we only will have data for a breakpointUncertainty if there is this region, on the contrary there is not
  ##we know from which breakpoint comes for because it is in the name
  
  ########################
  ## Analysis continues if there is at least one gene in the Domains
  ## if no gene in any domain
  ## apparently this SV will not affect any gene, so no possible to predict effect
  
  genesInvolved<-unique(info_affectedGenes$genesPosition$geneSymbol)
  
  if(length(genesInvolved)==0){
    print("NO GENES INVOLVED IN THIS STRUCTURAL VARIATION")
  }else{
    ###Since at least a gene potentially affected
    ##Let's continue with the analysis
    
    
    ##############
    ##function to compute score matrices for affected genes
    ####Constructing rearrangement
    ####And evaluating the genes situation
    ##We introduce the type of SV, in this first case inversion
    matrixesGenesEvaluation<-evaluatingGenesSituation(genesInfo = info_affectedGenes,
                                                      domainsInfo = info_affectedRegions,
                                                      enhancersInfo = info_affectedEnhancers,
                                                      dataPatient = patientInfo,
                                                      humanPhenotype = humanBased_genePhenotype,
                                                      micePhenotype = miceBased_genePhenotype,
                                                      polyCombGenes = rts_allGenes,##genes with polycomb, now we have the table with RTS
                                                      enhMap = MasterEnh_map,
                                                      tau_exp_scores = tau_exp_scores,
                                                      poorExp = poorlyExpressedGenes,
                                                      natureLof = nature_lof_scores,##HI scores nature 2016
                                                      huangLof = huangScores,##Huang Hi scores 2010,
                                                      clinGen_hiInfo = clinGen_hiInfo, ##ClinGene HI info
                                                      geneExp = Master_GeneExpression,##fpkm for the genes for the different stages
                                                      mainPatientPhenotype = formatedPhenotype,
                                                      phasesVector = phase)##Modification, now we want to create the matrixesGeneEvaluationFocuse on Just one gene
    
    ## Up to here if we have a gene in the uncertainty region the balance of kept gain enh for it corresponds with NA
    ## In addition if enh in uncertainty, the are not considered for the kept or gained categories, they are treated
    ## as if they were lost
    
    
    ###########################################################
    ## RANKING GENES BASED ON THEIR RELATIONSHIP WITH DISEASE
    ############################################################

    ##Running prediction
    resultsRankingPrediction<-rankingGenes(genesData = matrixesGenesEvaluation,
                                    phasesVector = phase,
                                    runMode = runMode)
    
    ##Genes Scoring
    scoresDisease<-resultsRankingPrediction$ScoresResults
    
    ##Updated Matrixes Genes Evaluation regulatory mechanism column (for those cases not clear eg Direct_LongRange_GeneDuplication)
    matrixesGenesEvaluation<-resultsRankingPrediction$MatrixesGeneEvaluation
      
    
    ##Adding type domain info, to know if regulatory domain disrupted or left intact
    scoresDisease[[phase]]$"TypeDomainInitial"<-matrixesGenesEvaluation[rownames(scoresDisease[[phase]]),
                                                                        "TypeDomainInitial"]
    
  }
  
  #########################################################
  ## INFO TO RETURN
  #########################################################
  
  resultsData<-list()
  
  if(length(genesInvolved)==0){
    ##So, no genes involved in this SV
    resultsData<-NULL##asi podemos hacer rbinds
  }else{
    ##So at least a gene involved in this SV
    ########
    ## Record matrixes genes evaluation in the output
    resultsData$matrixesGenesEvaluation<-matrixesGenesEvaluation
    
    ###store genes scores sorted from higher to smaller
    for(phase in names(scoresDisease)){
      matrixPhase<-scoresDisease[[phase]]
      matrixPhase<-matrixPhase[order(matrixPhase$Max_Score,decreasing = TRUE),]
      resultsData[[phase]]<-matrixPhase
    }
  }
  
  
  #####################################################
  ## We need to track also affected Domain coordinates
  ## At least for graphical outputs
  resultsData$affectedRegions<-info_affectedRegions
  resultsData$affectedEnhancers<-info_affectedEnhancers
  resultsData$affectedGenes<-info_affectedGenes
  resultsData$phasesVector<-phasesVector
  resultsData$formatedPhenotype<-formatedPhenotype
  
  #For future, we need to name each FPKM column and enh source equally, to be able to match them
  # resultsData$phases<-unique(MasterEnh_map$source)##Then we use this to automatize the phase names for other functions
  #And regarding the phenotype chosen we will select one or the other file
  return(resultsData)
}