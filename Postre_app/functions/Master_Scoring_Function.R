############################################################################################
## Root Function designed  to compute patient prediction, based only on patient information
############################################################################################

master_scoring_function<-function(patientInfo,runMode, user_tadMapInfo, MultiDataList){
  
  ###For the Studied Phenotype. For each of the Developmental Stages or Phases. 
  ##We need to load TAD its corresponding TAD maps and run predictions
  
  #######################################################################################
  ## LOADING GENOMIC DATA specific for the phenotype of interest
  ## Expression Data, Enhancers Data, TADs data for the different developmental phases
  ########################################################################################
  # source("functions/GenomicData_Loader.R")
  
  genomic_data<-genomic_data_loader(patientPheno = patientInfo$Phenotype)

  Master_GeneExpression<-genomic_data$Master_GeneExpression #The one I need to activate when everything fine
  MasterEnh_map<-genomic_data$MasterEnh_map
  Master_RegulatoryDomains<-genomic_data$Master_RegulatoryDomains##List of TAD maps for the different developmental stages
  phasesVector<-genomic_data$phasesVector ##Contains the name of the different developmental stages, or phases considered
  formatedPhenotype<-genomic_data$formatedPhenotype
  genomeBrowser_links<-genomic_data$genomeBrowser_links
  
  ########################################################
  #Modifying here TAD info if user using its own TAD map
  ########################################################
  
  if(patientInfo$userTADmap=="yes"){

    ##Update TAD map info, modify it for the user one
    for(targetPhase in phasesVector){
      ##Erase previous info
      Master_RegulatoryDomains$TADs[[targetPhase]]<-NULL
      Master_RegulatoryDomains$betweenTADs[[targetPhase]]<-NULL
      
      ##Record new info
      Master_RegulatoryDomains$TADs[[targetPhase]][["TADmapUser"]]<-user_tadMapInfo$TAD_map
      Master_RegulatoryDomains$betweenTADs[[targetPhase]][["betweenTADs_TADmapUser"]]<-user_tadMapInfo$Between_TAD_map
    }
    
  }
  
  ##############
  ##Function predicting etiology based on a SINGLE tad map && patient info
  ###############
  # source("functions/PatientPrediction_basedOnSingleTADmap.R", 
  #        local = TRUE)
  # 
  # ##Try to put it here to avoid massive reloading
  # source("scripts_To_Load_Data/cargarFunciones.R",local = TRUE)
  # source("functions/RankingGenes_Deciphering_etiology.R", local = TRUE)

  
  ##############################################
  ##Carrying prediction per Developmental Stage
  ##############################################
  
  resultsPerPhase<-list()##To store MAIN results per phase
  resultsPerPhase_secondaryInfo<-list()##To store additional information such as evaluationMatrixes per TAD map per Phase
  for(phase in phasesVector){
    # print(phase) ##Line commented 11/09/2022
    resultsPerTADmap<-list()
    # browser()
    ##Only 1 tadmap per phase, so remove for loop
    nameMap<-names(Master_RegulatoryDomains$TADs[[phase]])[1]
    
    ##Here we are dealing with the regulatory domain maps associated to each developmental condition
    #Which can be the same for all stages (if no specific data available).
    #OR stage specific. The best, as regulatory domains can vary in time and space

    ##change mapTADs
    mapTads<-Master_RegulatoryDomains$TADs[[phase]][[nameMap]]
    
    ##change mapBetweenTads domains
    nameBetweenTadsMap<-paste0("betweenTADs_",nameMap)
    betweenTads<-Master_RegulatoryDomains$betweenTADs[[phase]][[nameBetweenTadsMap]]
    
    #####Computing prediction per Patient
    ##This function needs to get the expression data,enhancer data AND PHASE to focus on it for the prediction
    
    patientAnalysis<-patientPrediction_basedOnSingleTADmap(patientInfo = patientInfo,
                                                           mapTads = mapTads,
                                                           betweenTads = betweenTads,
                                                           Master_GeneExpression = Master_GeneExpression,
                                                           MasterEnh_map = MasterEnh_map,
                                                           phase=phase, ##IsOnly1Phase
                                                           phasesVector=phasesVector,##Not removed in case captured somewhere for some output
                                                           formatedPhenotype = formatedPhenotype,
                                                           runMode = runMode,
                                                           MultiDataList = MultiDataList)
    resultsPerTADmap[[nameMap]]<-patientAnalysis
    
    ####################################
    ## Integrating results per TAD map
    ####################################

    # source("functions/Integrating_TAD_Predictions.R",
    #        local = TRUE)
    ##Attending to the current phase, integrate the results of the different TAD maps
    ##In case there is more than one TAD map used
    resultsIntegratingTADs<-integratingTAD_predictions(resultsPerTADmap = resultsPerTADmap,
                                                       targetPhase = phase) 
    
    #####################################################
    ## Appending Phase Results to patientResults
    #####################################################
    ##Just to have the info, to have all the data we want
    resultsPerPhase_secondaryInfo[[phase]]<-resultsPerTADmap
    
    ###Main Results. Results per Phase
    resultsPerPhase[[phase]]<-resultsIntegratingTADs
    
  }
  
  ##Summary, of all affected genes, where are they located (used for graphical display)
  # source("functions/summary_positionalInfoGenes.R")
  allAffectedGenes_positionalInfo<-summary_positionalInfoGenes(resultsPerPhase = resultsPerPhase,
                                                               onlyProteinCoding = TRUE,
                                                               gtf_annotation = MultiDataList$gtf_annotation)
  
  patientResults<-list("patientInfo" = patientInfo,
                       "resultsPerPhase_secondaryInfo" = resultsPerPhase_secondaryInfo,##Aqui poner el results per secondary TADmap
                       "resultsPerPhase" = resultsPerPhase,
                       "allAffectedGenes_positionalInfo" = allAffectedGenes_positionalInfo,
                       "genomeBrowser_links" = genomeBrowser_links,
                       "MasterEnh_map"= MasterEnh_map)
  
  #######################################
  ## Getting masterSummaryResultsMatrix
  #######################################
  # source("functions/Master_Summary_Matrix_IntegratesMainPhaseResults.R",
  #        local = TRUE)
  
  patientResults$masterSummaryResultsMatrix<-master_Summary_resultsMatrix(patientResults = patientResults)
  
  #####################################
  ## Generating Run//Prediction unique code
  ## Will be appended to images. It will change each time input modified
  ## or alsoWithTime
  #Created to avoid caches issues reloading old images
  #For instance if pathogenic deletion neural crest late affecting TFAP2A, image name 
  #same as if pathogenic inversion neural crest affecting TFAP2A, 
  #The browser will load the oldest image... even if it does not exist at the folder... using from cache
  #So append id to image name containing exactly all the input submitted info
  #So if when re-submitting, anything changes on the input from the last submission
  #The image will be reloaded
  #####################################

  ##Build with time stamp to be sure 100% is unique because when allowing tad upload this problem appeared again
  #this will change each second
  job_UniCode<-Sys.time()
  job_UniCode<-gsub(pattern = ",",replacement = "", x=job_UniCode)
  job_UniCode<-gsub(pattern = "-",replacement = "", x=job_UniCode)
  job_UniCode<-gsub(pattern = " ",replacement = "", x=job_UniCode)
  job_UniCode<-gsub(pattern = ":",replacement = "", x=job_UniCode)
  
  ##Now to formatted time stamp adding type SV and formated Phenotype, with that more than enough
  job_UniCode<-paste0(job_UniCode,
                      "_",
                      patientInfo$TypeSV,
                      "_",
                      patientResults$resultsPerPhase_secondaryInfo[[1]][[1]]$formatedPhenotype,
                      collapse = "")

  ##Remove ","
  # job_UniCode<-gsub(pattern = ",",replacement = "", x=job_UniCode)
  
  #Adding jobUnicode to patient results
  patientResults$job_UniCode<-job_UniCode
  #######################################
  ## Error Control
  #######################################
  ##If script reached this point, no error rised
  ##Distinguish between normal situation, where some genes predicted, with respect to that one where
  ##No gene is associated with SV
  
  ##Check if it is character or not, prior to the comparison, to handle appropriate data types
  if(is.character(patientResults$masterSummaryResultsMatrix)){
    ##Hence we should have the following string
    if(patientResults$masterSummaryResultsMatrix == "NO GENES ASSOCIATED WITH SV"){
      ##To highlight the anomaly of NO genes associated with SV
      ##Happening for controls, non pathogenic SVs
      patientResults$Status<-"OK, but NO genes associated with SV"
    }else{
      stop("THIS CODE PART SHOULD BE UNREACHABLE")
    }
  }else{
    ##So it is a dataframe, and there are results
    patientResults$Status<-"OK" 
  }
  
  ##Return results
  return(patientResults)
  
}
