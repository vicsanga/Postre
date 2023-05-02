##########################################################
## POSTRE R function for integration in pipelines
## master R function for performing multiple SV analyses
##########################################################

################################################################################################
#USER HAS TO SET THIS!!!!
# If you have downloaded POSTRE github folder, uncompress and provide path to Postre_app folder, eg:
# setwd("~/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/Postre/Postre_app")
setwd("~/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/Postre/Postre_app")
################################################################################################


POSTRE_multiSV<-function(SVs, userTADmap=NULL, runMode="Standard", genePhenoMatch="yes"){
  
  ##############################################
  ##NO TOUCHING FROM HERE!!
  ##############################################
  consideredPheno<-c("head_neck",
                     "cardiovascular",
                     "limbs",
                     "neurodevelopmental")##As more phenos considered they will appear here
  
  ####################################
  ###Let's load required Functions
  ####################################
  # Load data into the default environment inside the function
  source("scripts_To_Load_Data/metaFunctionLoad.R",local = environment())
  
  ##Required object for Single Prediction
  ##Loading multidata object, to avoid multiple reloading
  load("data/MultiDataList.RData")
  
  
  ##Scores For Considering Relevant result
  ##minRequired scores for heatmap coloring and geneReport generation
  minScore<-0.8
  highScore<-0.9
  
  
  ###Meter como variables, todos los valores que  estan estaticos en el menu
  relevantChr<-c(paste("chr",1:22,sep = ""), "chrX")##chrY excluded not all data available for chrY
  
  #genePhenoMatch options: "yes"-"no"
  
  ##Setting column and row names to SVs file
  colnames(SVs)<-c("chr_Break1","coord_Break1","chr_Break2",
                         "coord_Break2","TypeSV","Phenotype","patientID") ##For internal usage, patientID, for frontEnd SV_ID. To clearly show that a patient can carry multiple SVs
  
  rownames(SVs)<-SVs$patientID
  
  ############################################################################################################################
  ## Ensuring data in proper character format, for the coordinates, needed in character 
  ## to avoid downstream internal errors when splitting by "," when range provided
  ##########################################################################################################################
  SVs$chr_Break1<-as.character(SVs$chr_Break1)
  SVs$coord_Break1<-as.character(SVs$coord_Break1)
  
  SVs$chr_Break2<-as.character(SVs$chr_Break2)
  SVs$coord_Break2<-as.character(SVs$coord_Break2)
  
  ###########################################
  ## Handling parameters & adding to table
  ###########################################
  SVs$runMode<-runMode
  SVs$genePhenoConsideration<-genePhenoMatch
  SVs$refGenome<-"hg19"
  SVs$userTADmap<-"no" ##Default option, used for UCSC warning when yes to upload TAD coord
  
  #For user TAD map, in principle is set to FALSE to track whether processing done
  userTadProcessed<-FALSE
  user_tadMapInfo<-list()##It will be kept empty unless user provides its own TADmap
  
  ##################################
  ## Performing predictions
  ##################################
  all_patientResults<-list()##To store results
  ###Tracker variables
  nPatient<-0  ##Tracking Status
  cohortTractablePhenos<-character() ##To track patient provided phenotypes, only sections related to patient phenotypes will be shown. No sense on showing cardiovascular, if no cardiovascular.
  ##Phenos considered when we also have data form them (eg if patient limb but yet no limb data also not section)
  
  userTadProcessed<-FALSE ##To track if own TAD has to be processed and doing it only once 
  start.time<-Sys.time()
  
  for(patient in rownames(SVs)){
    
    nPatient<-nPatient+1
    cat("n Structural Variant: ",nPatient,"\n")
    
    ########################################
    ## Selecting Patient Info
    patientInfo<-SVs[patient,]
    
    all_patientPheno<-unlist(strsplit(x = patientInfo$Phenotype, split = ",", fixed = TRUE))
    
    
    for(pheno in all_patientPheno){
      
      ##If pheno among the ones considered in the app, run prediction
      if(pheno %in% consideredPheno){
        
        ############################
        ## Carrying prediction #####
        ############################
        cohortTractablePhenos<-c(cohortTractablePhenos, pheno)##We track this phenotype, and we will provide information for it in the multi pat section
        
        patientResults<-list()
        
        patientResults<-tryCatch({
          
          ##Check if user has selected its own TAD map for the analysis
          ##When TRUE is because NO TAD map selected
          ##For between TAD map generation the reference genome is taken into consideration for chr sizes
          
          if(is.null(userTADmap) == FALSE){
            
            ##So, user selected TAD map
            ##Check that TADmap NOT already processed
            if(userTadProcessed == FALSE){
              
              ##So, not processed yet
              ##Do processing
              
              colnames(userTADmap)<-c("chr","start","end") ##For internal usage, patientID, for frontEnd SV_ID. To clearly show that a patient can carry multiple SVs
              
              ## Filter for chromosomes of interest
              ## Consider all chr among all patients
              patientBreakpoints<-unique(c(SVs$chr_Break1,
                                           SVs$chr_Break2))
              
              ##Thus to speed computations, skiping info of non required chromosomes
              userTADmap<-userTADmap[userTADmap$chr %in% patientBreakpoints,]
              
              ##Getting between TAD map
              userBetweenTADmap<-getBetweenTADmap(TADmap=userTADmap)
              
              ##Capturing processed information  
              user_tadMapInfo$TAD_map<-userTADmap
              user_tadMapInfo$Between_TAD_map<-userBetweenTADmap
              
              ##Recording that TAD map has already been processed
              userTadProcessed<-TRUE
            }
            
            ##To track wether user_tadMapInfo must be used or not in downstream functions
            patientInfo$userTADmap<-"yes"
            
            ##Updating patientInfo with userTADmap "yes", so that it is clarified that the user tad map has been used for the prediction
            SVs$userTADmap<-"yes" 
          }
          
          ##Working on each pheno separately, running prediction
          monoPheno_patientInfo<-patientInfo
          monoPheno_patientInfo$Phenotype<-pheno 
          
          ##If there is an error the following instruction will not be terminated
          ##In multiple screening we do not generate graphics,so only master_scoring_function used, and not the wrapper for graphics one
  
          patientResults<-suppressWarnings(
            master_scoring_function(patientInfo = monoPheno_patientInfo, runMode = runMode, user_tadMapInfo=user_tadMapInfo, MultiDataList = MultiDataList)
            )
          
          ##If there was no error patientResults$Status == "OK" or "OK, but NO genes associated with SV"
        },error = function(err){
          patientResults$Status<-"ERROR"
          return(patientResults)
          
        })
        ################################################
        ## Storaging Results per Patient & Phenotype
        ################################################
        ##To avoid the object size to be unnecessary huge, only mantain strictly necessary info
        ##We can maybe even simplify this more by using the matrix behind heatmap
        ##but let's see for now how it goes
        patientResults$resultsPerPhase_secondaryInfo<-NULL
        patientResults$genomeBrowser_links<-NULL
        patientResults$allAffectedGenes_positionalInfo<-NULL
        patientResults$MasterEnh_map<-NULL
        patientResults$resultsPerPhase<-NULL
        
        all_patientResults[[patient]][[pheno]]<-patientResults 
        
      }
    }
  }

  ###################################################
  ## Parsing Patient Results Once Prediction is done
  ###################################################
  cohortTractablePhenos<-unique(cohortTractablePhenos)##Phenos that have appeared on the patients, and that we can predict impact
  ##For the phenos that appear in the patients
  
  cohort_results<-cohortResults_Parser(minScore = minScore, all_patientResults = all_patientResults,
                                       consideredPheno = cohortTractablePhenos,
                                       discardRelevantByBrokenGene = FALSE,
                                       AllPatientsInfo = SVs)
  

  ##########################
  #Creating output object
  ##########################
  
  outpInfo<-list()
  outpInfo$pathogenicityPrediction_per_SV<-cohort_results$candidateGenesInfo
  #Removing some info
  outpInfo$pathogenicityPrediction_per_SV$TypeSV<-NULL
  outpInfo$pathogenicityPrediction_per_SV$N_LR_Mech<-NULL
  outpInfo$pathogenicityPrediction_per_SV$N_Direct_Mech<-NULL
  
  outpInfo$geneStats_recurrencyAndPathomech<-cohort_results$geneRecurrencyInfo
  outpInfo$errors<-cohort_results$error_General_Info
  return(outpInfo)
  
}

