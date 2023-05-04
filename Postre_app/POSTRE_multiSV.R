##########################################################################################
## Generation of POSTRE R function for performing analyses of multiple SV in an R script
##########################################################################################

### DO NOT MODIFY ANYTHING

# This script generates the function:
# POSTRE_multiSV(SVs, pathTo_Postre_app_Folder, userTADmap=NULL, runMode="Standard", genePhenoMatch="yes")
# Which allows to use POSTRE in an R script with just a function
# To use POSTRE_multiSV in an R script check instructions in POSTRE GitHub page https://github.com/vicsanga/Postre

########################################
## Generating POSTRE_multiSV function
########################################

POSTRE_multiSV<-function(SVs=NULL, pathTo_Postre_app_Folder=NULL, userTADmap=NULL, runMode="Standard", genePhenoMatch="Yes"){

if(is.null(SVs)){
    
cat('

The function presents the following signature:
POSTRE_multiSV(SVs=NULL, pathTo_Postre_app_Folder=NULL, userTADmap=NULL, runMode="Standard", genePhenoMatch="Yes", )

It is mandatory to provide the values for 2 arguments: SVs and pathTo_Postre_app_Folder


####################################################
Information and options for the  input parameters:
####################################################


####
SVs
####
Data frame with 7 columns containing SVs information (you can use as test file any of those provided in the GitHub testFiles folder)
 Column 1: Chromosome for the breakpoint 1
 Column 2: Genomic coordinates (hg19) for the breakpoint 1. When not base pair resolution, provide a comma separated range, e.g. 85092268,85092269.
 Column 3: Chromosome for the breakpoint 2
 Column 4: Genomic coordinates (hg19) for the breakpoint 2. When not base pair resolution, provide a comma separated range, e.g. 85092268,85092269.
 Column 5: Structural Variant Type. Current options: Inversion, Translocation, Deletion or Duplication.
 Column 6: Comma separated list of phenotypes associated with the structural variant. Current options are: head_neck, limbs, neurodevelopmental or cardiovascular. For instance: head_neck,neurodevelopmental,cardiovascular.
 Column 7: Structural variant unique identifier e.g. (Patient1_SV3)

 Note: For the case of structural variants happening strictly in one chromosome (deletions, inversions, duplications) the breakpoint 1 is the one associated with a smaller genomic coordinate, and the breakpoint 2 the one associated with a larger genomic coordinate. For translocations, it does not matter.


###########################
pathTo_Postre_app_Folder
###########################
After downloading POSTRE GitHub repository, uncompress it and provide the path to "Postre_app" folder:
 i.e. pathTo_Postre_app_Folder="/home/victor/Downloads/Postre-main/Postre_app"

##############
 userTADmap
############## 
This parameter can be NULL. If not NULL, it has to be a data frame with 3 columns (chr, start, end) with the TAD coordinates
 if a TAD map is provided, it will be taken as reference for the regulatory domains definition for all cell types
 if a TAD map is not provided (default, userTADmap=NULL) the ones considered in POSTRE database are taken as reference

#########
runMode
#########
2 options ("Standard" or "High-Specificity")

###############
genePhenoMatch
###############
2 options ("Yes","No")
 if "Yes" (default) it is required a known association of the candidate genes with the patient phenotype for being pathogenic candidates (in terms of gene-phenotype association requirements)
 if "No" it is enough for the candidate genes to be associated with any disease to be pathogenic candidates (in terms of gene-phenotype association requirements)


#####################################################
Information about output object
#####################################################

The object returned by POSTRE_multiSV is a list with 3 different elements:

#################################
pathogenicityPrediction_per_SV
#################################
This table provides a pathogenic prediction for each of the SVs and associated phenotypes analyzed. It presents different columns:
   
 + SV ID: Identifier of the SV as provided in the file uploaded to POSTRE.
 + Phenotype: Phenotype associated with the SV considered for pathogenic evaluation.
 + Pathogenic Score: Pathogenic score (0-1) computed for the SV-phenotype association. It corresponds with the maximum pathogenic score computed for all the candidate genes.
 + Pathogenic: It indicates whether the SV-phenotype association is predicted pathogenic (Yes) or not (No).
 + Causative genes: List of genes predicted as disease causative. This cell will be empty if pathogenicity is not predicted.
 + Candidate genes (Pathogenic Score): List of candidate genes (gene whose regulatory domain (TAD) or sequence (e.g gene deletion) is altered by a SV). A candidate gene is not necessarily involved in the disease etiology (i.e. candidate genes include both causative and non-causative genes). For each candidate gene, the maximum pathogenic score (0-1) computed along all cell types considered is provided in brackets.


###################################
geneStats_recurrencyAndPathomech
###################################
This table is an aggregation of the pathogenic predictions per gene, phenotype and pathogenic mechanism (coding, long-range). It presents different columns:
   
 + Gene: Name of the gene.
 + Phenotype: Phenotype associated with the gene where pathogenicity has been predicted.
 + N SVs: Number of SVs where pathogenicity has been predicted (either by coding or long-range pathogenic mechanisms).
 + N SVs Long-Range: Number of SVs where pathogenicity has been predicted through a long-range (enhancer mediated) pathogenic event.
 + N SVs Coding: Number of SVs where pathogenicity has been predicted through a coding mechanism (e.g. gene deletion).
 + SV IDs: Identifier of the SVs where pathogenicity has been predicted (either by coding or long-range pathogenic mechanisms).
 + Long-Range SV IDs: Identifier of the SVs where pathogenicity has been predicted through long-range (enhancer mediated) pathological mechanisms.
 + Coding SV IDs: Identifier of the SVs where pathogenicity has been predicted through coding (e.g. gene deletion) pathological mechanisms.


##########
errors
##########
Contains the ids of the SVs which rised an error (if they occur) during their interpretation.')
    
    
    
    
    return(NULL)
  }
  
  
  
  
  
  ##Getting current working directory to return to it after function execution
  user_workingDirectory<-getwd()
  
  #Setting as working directory Postre_app folder
  if(dir.exists(pathTo_Postre_app_Folder)){
    setwd(pathTo_Postre_app_Folder)
  }else{
    cat('Error: The path to Postre_app folder does not exist. Did you put it correctly? \n')
    return(NULL)
  }
  
  ##Checking genePhenoMatch input, and also ensuring lowercase characters for downstream analyses
  genePhenoMatch<-tolower(genePhenoMatch)
  if((genePhenoMatch != "yes") & (genePhenoMatch != "no")){
    cat('Error: Wrong input parameter, did you set genePhenoMatch to Yes or No? \n')
    return(NULL)
  }
  
  ##############################################
  consideredPheno<-c("head_neck",
                     "cardiovascular",
                     "limbs",
                     "neurodevelopmental")
  
  ####################################
  ###Let's load required Functions
  ####################################
  source("scripts_To_Load_Data/metaFunctionLoad.R", local = TRUE)
  
  ##Required object for Single Prediction
  ##Loading multidata object, to avoid multiple reloading
  load("data/MultiDataList.RData")
  
  ##Scores For Considering Relevant result
  ##minRequired scores for heatmap coloring and geneReport generation
  minScore<-0.8
  highScore<-0.9
  
  relevantChr<-c(paste("chr",1:22,sep = ""), "chrX")##chrY excluded not all data available for chrY
  
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
              
              colnames(userTADmap)<-c("chr","start","end") 
              
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
            
            ##To track whether user_tadMapInfo must be used or not in downstream functions
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
  
  ##Defining again as working directory the one that the user was using before executing the function
  setwd(user_workingDirectory)
  
  ##Returning output object
  return(outpInfo)
  
}

