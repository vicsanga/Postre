####################################################################################
## Genomic data loader
## Function to load the specific genomic data for the phenotype of interest
## Gene Expression, Enhancer Maps and TAD maps for the different developmental phases
## And also an object containing the names of the phases considered in each phenotype
#####################################################################################
# setwd("~/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/SV_app")

###########################################
##Info for multiple patient condition
##########################################
##Important object, pheno-phasesTable
##########################################

pheno_Phases<-list(
  ##Phases are obtained from their own phasesVector
  
  ##head_neck
  ##"head_neck"=c("NeuralCrestEarly", "NeuralCrestLate", "PalateCS20", "PalateCSMix")
  ##attempt removing phaseFree
  #"head_neck"=c("NeuralCrestEarly", "NeuralCrestLate", "PalateCS20", "PalateCSMix"),
  "head_neck"=c("NeuralCrestEarly", "NeuralCrestLate", "PalateCS20"),
  
  ##Cardiovascular
  "cardiovascular"=c("day5", "day7", "day15", "day80"),
  
  ##Limb
  "limbs"=c("EmbryonicLimb1", "EmbryonicLimb2"),
  
  ##Brain (behaviour, neurological, cognitive)
  ##"behaviour_neurological_cognitive"=c("PfcGw15", "PfcGw18")
  #Dudas que tenia, al final termino mas generico, nervous system
  # "nervous_system"=c("PfcGw15", "PfcGw18"),
  
  ##NEURODEVELOPMENTAL
  "neurodevelopmental"=c("PfcGw15", "PfcGw18"),
  
  ##VISION-EYE
  "vision_eye"=c("Retina","RPE")
  
  ##...
)



################
###FUNCTION
genomic_data_loader<-function(patientPheno){
  
  ##patientPheno value is given by App input, either through the single or multiple condition
  if((patientPheno == "Head & Neck") || (patientPheno == "head_neck")){
    
    ###########################################################
    ## Loading fpkms for the genes 
    #Prepared in:/home/victor/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/RNA-Seq/integrating_rnaseqData_forSoftware/fpkms_dataframesCreation/headNeck_fpkms/3_Master_HeadNeck_fpkms.R
    load(file = "data/specificData_PerPhenotype/head_neck/Master_GeneExpression_HeadNeck.RData")
    
    ###########################################################
    ##Load master enhancer map
    ##Created in: /home/victor/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/enhacers/software_MasterEnhMap_creation/HeadNeck/3_AddingAcetilation_CreatingDefinitiveMasterEnhMap.R
    load("data/specificData_PerPhenotype/head_neck/Master_EnhMap_HeadNeck.RData")
    
    ############################################################
    ## Load TAD maps data per developmental Stage
    ##Created in:~/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/preparingDataForSoftware/genomicData_PhenotypeSpecific/Head_Neck/TAD_maps/1_HeadNeck_Preparing_TAD_map_Objects.R
    load(file ="data/specificData_PerPhenotype/head_neck/Master_RegulatoryDomains_HeadNeck.RData")
    
    ################################################################################################################################
    ##Loading name of the phases, which will be exactly the same, to do the matching, between Expression AND Enhancers AND TAD maps
    phasesVector<-c("NeuralCrestEarly", "NeuralCrestLate", "PalateCS20")
    #####
    ## Correct phenotype category to be properly processed afterwards
    ##Hence the way the software handles it, without, spaces and so on
    formatedPhenotype<-"head_neck"
    
    
    ###############################
    ##Cargar UCSC_links table relating each dev stage with its corresponding ucsc session
    #colnames: Phenotype	Stage	SessionId	SessionLinkBaseName
    
    genomeBrowser_links<-read.delim(file="data/specificData_PerPhenotype/head_neck/head_neck_UCSC_links.tab",
                           sep="\t",
                           header = FALSE,
                           stringsAsFactors = FALSE)
    colnames(genomeBrowser_links)<-c("Phenotype", "Stage", "SessionId", "SessionLinkBaseName")
    
  }else if((patientPheno == "Cardiovascular") || (patientPheno == "cardiovascular")){
    ##Loading Data For Heart
    
    
    ###########################################################
    ## Loading fpkms for the genes 
    load(file = "data/specificData_PerPhenotype/heart/Master_GeneExpression_Heart.RData")
    
    ###########################################################
    ##Load master enhancer map
    load("data/specificData_PerPhenotype/heart/Master_EnhMap_Heart.RData")
    
    ############################################################
    ## Load TAD maps data per developmental Stage
    ##Created in:~/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/preparingDataForSoftware/genomicData_PhenotypeSpecific/Head_Neck/TAD_maps/1_HeadNeck_Preparing_TAD_map_Objects.R
    load(file ="data/specificData_PerPhenotype/heart/Master_RegulatoryDomains_Heart.RData")
    
    ################################################################################################################################
    ##Loading name of the phases, which will be exactly the same, to do the matching, between Expression AND Enhancers AND TAD maps
    #phasesVector<-c("day5", "day7", "day15", "day80")##
    phasesVector<-c("day5", "day7", "day15", "day80")##
    #####
    ## Correct phenotype category to be properly processed afterwards
    ##Hence the way the software handles it, without, spaces and so on
    formatedPhenotype<-"cardiovascular"
    
    
    ###############################
    ##Cargar UCSC_links table relating each dev stage with its corresponding ucsc session
    #colnames: Phenotype	Stage	SessionId	SessionLinkBaseName
    
    genomeBrowser_links<-read.delim(file="data/specificData_PerPhenotype/heart/cardiovascular_UCSC_links.tab",
                                    sep="\t",
                                    header = FALSE,
                                    stringsAsFactors = FALSE)
    colnames(genomeBrowser_links)<-c("Phenotype", "Stage", "SessionId", "SessionLinkBaseName")
    
    
  }else if((patientPheno == "Limbs") || (patientPheno == "limbs")){
    ##Loading Data For Limbs
    
    
    ###########################################################
    ## Loading fpkms for the genes 
    load(file = "data/specificData_PerPhenotype/limb/Master_GeneExpression_Limb.RData")
    
    ###########################################################
    ##Load master enhancer map
    load("data/specificData_PerPhenotype/limb/Master_EnhMap_Limb.RData")
    
    ############################################################
    ## Load TAD maps data per developmental Stage
    ##Created in:~/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/preparingDataForSoftware/genomicData_PhenotypeSpecific/Head_Neck/TAD_maps/1_HeadNeck_Preparing_TAD_map_Objects.R
    load(file ="data/specificData_PerPhenotype/limb/Master_RegulatoryDomains_Limb.RData" )
    
    ################################################################################################################################
    ##Loading name of the phases, which will be exactly the same, to do the matching, between Expression AND Enhancers AND TAD maps
    phasesVector<-c("EmbryonicLimb1", "EmbryonicLimb2")
    
    #####
    ## Correct phenotype category to be properly processed afterwards
    ##Hence the way the software handles it, without, spaces and so on
    formatedPhenotype<-"limbs" ##the way it is in the gene-pheno table
    
    
    ###############################
    ##Cargar UCSC_links table relating each dev stage with its corresponding ucsc session
    #colnames: Phenotype	Stage	SessionId	SessionLinkBaseName
    ##Acabar de crear ucsc sessions
    genomeBrowser_links<-read.delim(file="data/specificData_PerPhenotype/limb/limbs_UCSC_links.tab",
                                    sep="\t",
                                    header = FALSE,
                                    stringsAsFactors = FALSE)
    colnames(genomeBrowser_links)<-c("Phenotype", "Stage", "SessionId", "SessionLinkBaseName")
    
    
  }else if((patientPheno == "Neurodevelopmental") || (patientPheno == "neurodevelopmental")){
    ##Primer intento de brain-nervous system, en SV backup 5 oct 2021
    
    ###########################################################
    ## Loading fpkms for the genes 
    load(file = "data/specificData_PerPhenotype/neurodevelopmental/Master_GeneExpression_neurodevelopmental.RData")
    
    ###########################################################
    ##Load master enhancer map
    load("data/specificData_PerPhenotype/neurodevelopmental/Master_EnhMap_neurodevelopmental.RData")
    
    ############################################################
    ## Load TAD maps data per developmental Stage
    ##Created in:~/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/preparingDataForSoftware/genomicData_PhenotypeSpecific/...
    load(file ="data/specificData_PerPhenotype/neurodevelopmental/MRD_neurodevelopmental.RData")
    
    ################################################################################################################################
    ##Loading name of the phases, which will be exactly the same, to do the matching, between Expression AND Enhancers AND TAD maps
    phasesVector<-c("PfcGw15", "PfcGw18")
    
    #####
    ## Correct phenotype category to be properly processed afterwards
    ##Hence the way the software handles it, without, spaces and so on
    
    formatedPhenotype<-"neurodevelopmental"
    
    ###############################
    ##Cargar UCSC_links table relating each dev stage with its corresponding ucsc session
    #colnames: Phenotype	Stage	SessionId	SessionLinkBaseName
    genomeBrowser_links<-read.delim(file="data/specificData_PerPhenotype/neurodevelopmental/neurodevelopmental_UCSC_links.tab",
                                    sep="\t",
                                    header = FALSE,
                                    stringsAsFactors = FALSE)
    colnames(genomeBrowser_links)<-c("Phenotype", "Stage", "SessionId", "SessionLinkBaseName")
    
    
  }else if((patientPheno == "Vision-Eye") || (patientPheno == "vision_eye")){
    ##14 Agosto 2021, incorporacion fenotipo
    
    ###########################################################
    ## Loading fpkms for the genes 
    load(file = "data/specificData_PerPhenotype/vision_eye/Master_GeneExpression_Vision_Eye.RData")
    
    ###########################################################
    ##Load master enhancer map
    load("data/specificData_PerPhenotype/vision_eye/Master_EnhMap_VisionEye.RData")
    
    ############################################################
    ## Load TAD maps data per developmental Stage
    ##Created in:~/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/preparingDataForSoftware/genomicData_PhenotypeSpecific/...
    load(file ="data/specificData_PerPhenotype/vision_eye/Master_RegulatoryDomains_Vision_Eye.RData")
    
    ################################################################################################################################
    ##Loading name of the phases, which will be exactly the same, to do the matching, between Expression AND Enhancers AND TAD maps
    phasesVector<-c("Retina", "RPE")
    
    #####
    ## Correct phenotype category to be properly processed afterwards
    ##Hence the way the software handles it, without, spaces and so on
    
    formatedPhenotype<-"vision_eye" #Ha de ser igual que el de la columna de tablas gene-pheno relationships
    
    ###############################
    ##Cargar UCSC_links table relating each dev stage with its corresponding ucsc session
    #colnames: Phenotype	Stage	SessionId	SessionLinkBaseName
    genomeBrowser_links<-read.delim(file="data/specificData_PerPhenotype/vision_eye/vision_eye_UCSC_links.tab",
                                    sep="\t",
                                    header = FALSE,
                                    stringsAsFactors = FALSE)
    colnames(genomeBrowser_links)<-c("Phenotype", "Stage", "SessionId", "SessionLinkBaseName")
    
    
  }
  
  genomic_data<-list("Master_GeneExpression" = Master_GeneExpression,
                     "MasterEnh_map" = MasterEnh_map,
                     "Master_RegulatoryDomains" = Master_regulatoryDomains_list,
                     "phasesVector" = phasesVector,
                     "formatedPhenotype" = formatedPhenotype,
                     "genomeBrowser_links"=genomeBrowser_links)
  
  return(genomic_data)
}