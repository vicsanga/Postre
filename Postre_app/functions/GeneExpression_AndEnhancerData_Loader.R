###############################################################
## Loading Expression And Enhancer Data
## Regarding each of the possible phenotypes the files change
###############################################################

geneExp_Enh_loader<-function(patientPheno){
  
  ##patientPheno value is given by App input
  if(patientPheno=="Head & Neck"){
    
    ###########################################################
    ## Loading fpkms for the genes 
    #Prepared in:/home/victor/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/RNA-Seq/integrating_rnaseqData_forSoftware/fpkms_dataframesCreation/headNeck_fpkms/3_Master_HeadNeck_fpkms.R
    load(file = "data/specificData_PerPhenotype/head_neck/Master_GeneExpression_HeadNeck.RData")
    
    ###########################################################
    ##Load master enhancer map
    ##Created in: /home/victor/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/enhacers/analyzingAcetilationDistribution/1_RetrievingAcetilationLevelsPerEnhancer.R
    load("data/specificData_PerPhenotype/head_neck/Master_EnhMap_HeadNeck.RData")
    
    ############################################################
    ##Loading name of the phases, which will be exactly the same, to do the matching, between Expression AND Enhancers Data
    phasesVector<-c("Rada", "Prescot", "CS20", "CSMix")
  }
  
   #####
   ## Correct phenotype category to be properly processed afterwards
   ##Hence the way the software handles it, without, spaces and so on
   formatedPhenotype<-"head_neck"

   enh_exp_data<-list("Master_GeneExpression" = Master_GeneExpression,
                     "MasterEnh_map" = MasterEnh_map,
                     "phasesVector" = phasesVector,
                     "formatedPhenotype" = formatedPhenotype)
   
   return(enh_exp_data)
  
}
