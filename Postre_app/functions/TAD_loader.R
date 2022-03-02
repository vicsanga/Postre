###############################################################
## Loading Expression And Enhancer Data
## Regarding each of the possible phenotypes the files change
###############################################################

TAD_loader<-function(patientPheno){
  ##patientPheno value is given by App input
  if(patientPheno=="Head & Neck"){
    
    folderPath<-"data/softwareUsed_Maps/Head_Neck/"
    TADs_paths<-list.files(path = folderPath)
    
    TAD_maps<-list()
    
    for(TAD_file in TADs_paths){
      ##Tad Maps
      ## load all the elements in a folder
      # #print(TAD_file)
      wholePath<-paste0(folderPath,TAD_file)
      loadedTAD<-read.table(file = wholePath, header = F, sep="\t",stringsAsFactors = F)
      colnames(loadedTAD)<-c("chr","start","end")
      ## in case it is wrongly done, lets sort TADs regarding: 1st Chr, 2nd Start pos, 3rd End pos
      loadedTAD<-loadedTAD[order(loadedTAD$chr,loadedTAD$start,loadedTAD$end),]
      TAD_maps[[TAD_file]]<-loadedTAD
    }
    rm(loadedTAD)
    
    #########################
    ##Between TAD maps
    ##The coordinates of the regions between TADs
    load(file = "data/softwareUsed_Maps/BetweenTAD_Maps_Head_Neck.RData")
  }
  
  loadedData<-list("TAD_maps"=TAD_maps,
                   "between_TAD_maps"=between_TAD_maps)
  
  return(loadedData)
  
}
