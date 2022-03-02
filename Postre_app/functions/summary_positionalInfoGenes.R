##########################################################
## Retrieval Positional info All Affected Genes ##########
##########################################################
summary_positionalInfoGenes<-function(resultsPerPhase, onlyProteinCoding){
  
  ##We need gtf_annotation info to prepare the all genes coordinate object
  load("data/genesAnnotation.RData")
  
  ##Filtering gtf_annotation info
  if(onlyProteinCoding == TRUE){
    gtf_annotation<-subset(gtf_annotation, geneType=="protein_coding") 
  }

  allGenes<-character()##To track affected genes
  for(phase in names(resultsPerPhase)){
    allGenes<-c(allGenes, resultsPerPhase[[phase]][[phase]]$affected_gene)
  }
  allGenes<-unique(allGenes)##remove duplicated gene names, each appear twice, one for GOF one for LOF
  
  allAffectedGenesPosition<-gtf_annotation[allGenes,c("chr","start","end","TSS","strand","geneSymbol","geneType")]##object to store the infos
  allAffectedGenesPosition<-unique(allAffectedGenesPosition)
  
  return(allAffectedGenesPosition)
}

