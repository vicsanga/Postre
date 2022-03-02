###################################################################################
## To locate genes, just looking at their TSS
## Function to get the name of the genes located at the affected Regions 
## Atending to both the domains, and the breakpoint segments
###################################################################################
affectedGenes<-function(regionsAltered, gtf, onlyProteinCoding, patientInfo){
  ##regionsAltered for info_affectedRegions object
  domainsData<-regionsAltered$domainsAffected
  segmentsData<-regionsAltered$segmentsBreakP
  uncertaintyData<-regionsAltered$uncertaintyRegions
  
  ##Filtering gtf info
  if(onlyProteinCoding == TRUE){
    gtf<-subset(gtf, geneType=="protein_coding") 
  }

  #############################################
  #####Look for genes contained in the domains
  #####Just paying attention to TSS
  #####Doubt on whether paying attention to whole contained genes in future
  
  genesInDomain<-list()
  for(dominio in rownames(domainsData)){
    
    targetDomain_Chr<-domainsData[dominio,"chr"]
    targetDomain_Start<-domainsData[dominio,"start"]
    targetDomain_End<-domainsData[dominio,"end"]
    
    gtf$isChr<-gtf$chr==targetDomain_Chr
    gtf$isBelowUpperLimit<-targetDomain_End>=gtf$TSS
    gtf$isAboveLowerLimit<-targetDomain_Start<=gtf$TSS
    
    ##Is TAD holding the gene TSS?
    logicalInfoGenes<-gtf[,c("isChr","isBelowUpperLimit","isAboveLowerLimit")]
    isAGeneInside<-rowSums(logicalInfoGenes)==3##if ==3 fullfits all the requirements
    
    genesInDomain[[dominio]]<-gtf[isAGeneInside,"geneSymbol"]
  }
  
  ######################################################
  ## Look for genes contained in each breakpoint segment
  ## They do not necessarily need to be entirely inside (the whole gene)
  ## of the segment, if the region they cross is the domain border
  ## But ont the contrary, If they cross the breakpoint, in that case, they will be considered broken
  ## Then looking if any is broken, important
  genesInSegment<-list()
  for(segment in rownames(segmentsData)){
    
    targetSegment_Chr<-segmentsData[segment,"chr"]
    targetSegment_Start<-segmentsData[segment,"start"]
    targetSegment_End<-segmentsData[segment,"end"]
    
    gtf$isChr<-gtf$chr==targetSegment_Chr
    gtf$isBelowUpperLimit<-targetSegment_End>=gtf$TSS 
    gtf$isAboveLowerLimit<-targetSegment_Start<=gtf$TSS
    
    ##Is segment holding the gene TSS?
    logicalInfoGenes<-gtf[,c("isChr","isBelowUpperLimit","isAboveLowerLimit")]
    isAGeneInside<-rowSums(logicalInfoGenes)==3##if == 3 fullfits all the requirements
    
    genesInSegment[[segment]]<-gtf[isAGeneInside,"geneSymbol"]
  }
  
  #############################################################
  ## ¿Any gene has its coding sequence broken by the breakpoint?
  ## LeftSegment---BREAKPOINT---RigthSegment
  ## A gene is considered to be broken if the end of any of the left segments
  ## or the start of any of the right segments falls in the gene region (between its start or end)
  ## brokenGenes
  
  ##Broken genes considered in general (not discriminating by domain), if we have any, then handling in which it is
  brokenGenes<-character()
  
  ##For left segments we get the end
  for(segment in rownames(segmentsData)){
    
    targetSegment_Chr<-segmentsData[segment,"chr"]
    targetSegment_Start<-segmentsData[segment,"start"]
    targetSegment_End<-segmentsData[segment,"end"]
    
    if(grepl("left",segment, fixed = TRUE)){
      puntoRotura<-targetSegment_End
    }else{
      ##it means it is the right segment so we take the start
      puntoRotura<-targetSegment_Start
    }
    
    gtf$isChr<-gtf$chr==targetSegment_Chr
    gtf$isBreakBeforeEnd<-puntoRotura<=gtf$end ##gtf$end is the end of the gene, the 3' coord from the reference strand (evethough it can be the TSS, but it is the bigger pos)
    gtf$isBreakAfterStart<-puntoRotura>=gtf$start
    
    ##
    logicalInfoGenes<-gtf[,c("isChr","isBreakBeforeEnd","isBreakAfterStart")]
    isAGeneBroken<-rowSums(logicalInfoGenes)==3##if ==3 fullfits all the requirements
    brokenGenes<-c(brokenGenes, gtf$geneSymbol[isAGeneBroken])
  }
  brokenGenes<-unique(brokenGenes)##there can be repeated genes if we have a single coordinate or two super close for the breakpoint, in this case the broken gene will appear two times, so we do the unique
  
  
  ##############################################
  ## Looked for genes entirely contained in the 
  ## uncertainty regions. 
  ## If they are not entirely contained it is because
  ## they cross any of the breakpoint limits, 
  ## so they will be considered broken
  ###############################################
  
  ##Uncertainty genes considered in general (not discriminating by domain), if we have any, then handling in which it is
  genesInUncertainty<-character()
  
  for(breakp in names(uncertaintyData)){
    info_uncert<-uncertaintyData[[breakp]]
    if(is.null(info_uncert)==FALSE){
      ##so there is an uncertainty region
      ##Lets get the coordinates for this segment (uncertainty segment)
      ##If there is data we have a 1 row matrix with the region coord
      targetSegment_Chr<-info_uncert[1,"chr"]
      targetSegment_Start<-info_uncert[1,"start"]
      targetSegment_End<-info_uncert[1,"end"]
      
      gtf$isChr<-gtf$chr==targetSegment_Chr
      gtf$isBelowUpperLimit<-targetSegment_End>=gtf$end 
      gtf$isAboveLowerLimit<-targetSegment_Start<=gtf$start
      
      ##Is Uncertainty region holding the whole gene?
      ##If it is inside but not entirely, it will cross a segment limit
      ##then, it will be considered as broken
      logicalInfoGenes<-gtf[,c("isChr","isBelowUpperLimit","isAboveLowerLimit")]
      isAGeneInside<-rowSums(logicalInfoGenes)==3##if ==3 fullfits all the requirements
      
      genesInUncertainty<-c(genesInUncertainty, gtf[isAGeneInside,"geneSymbol"])
    }
  }
  
  genesInUncertainty<-unique(genesInUncertainty)
  
  
  ######################################################################################################
  ## If we are dealing with a Deletin or a Duplication
  ## That affects >1 TAD. If they are not adjacent
  ## All the genes that are in between are going to be directly affected, either deleted or duplicated
  ######################################################################################################
  
  ##Creating by default the objects, if it is an inversion or translocation, obj will not be used
  deletedGenes<-character()
  duplicatedGenes<-character()
  
  if((patientInfo$TypeSV=="Deletion") || (patientInfo$TypeSV=="Duplication")){
    
    if(patientInfo$TypeSV=="Deletion"){
      ##Look for Deleted genes
      ##All genes located between bkp1 and bkp2 deleted, does not matter whether they are on the affected domain or not
      ##If Affected Domain is not the same
      ##Defining Deletion area
      deletionArea_chr<-unique(segmentsData$chr)
      deletionArea_start<-segmentsData["breakpoint_1_right_segment","start"]
      deletionArea_end<-segmentsData["breakpoint_2_left_segment","end"]
      
      ##Check If affected Domain is the same,
      if(length(unique(domainsData$domainId))>1){
        ##So the affected Domain is not the same
        ##Hence looking also for deleted genes between domains
        gtf$isChr<-gtf$chr==deletionArea_chr
        gtf$isBelowUpperLimit<-deletionArea_end>=gtf$end 
        gtf$isAboveLowerLimit<-deletionArea_start<=gtf$start
        
        ##Catching genes
        logicalInfoGenes<-gtf[,c("isChr","isBelowUpperLimit","isAboveLowerLimit")]
        isAGeneInside<-rowSums(logicalInfoGenes)==3##if ==3 fullfits all the requirements
        
        deletedGenes<-unique(c(deletedGenes, gtf[isAGeneInside,"geneSymbol"]))
        
        
      }else if(length(unique(domainsData$domainId))==1){
        ##So the affected Domain is the same,
        ##Hence deleted genes only those in bkp1 right segment && bkp2 left segment
        ##BUT excluding those which are not shared, drawing intratad deletion to comprehend
        deletedGenes<-unique(intersect(genesInSegment$breakpoint_1_right_segment,genesInSegment$breakpoint_2_left_segment))
        
      } 
    }else {
      ##Working with duplications

      ##Look for duplicated genes
      ##All genes located between bkp1 and bkp2 duplicated, does not matter whether they are on the affected domain or not
      ##If Affected Domain is not the same
      ##Defining duplication area
      duplicationArea_chr<-unique(segmentsData$chr)
      duplicationArea_start<-segmentsData["breakpoint_1_right_segment","start"]
      duplicationArea_end<-segmentsData["breakpoint_2_left_segment","end"]
      
      ##Check If affected Domain is the same,
      if(length(unique(domainsData$domainId))>1){
        ##So the affected Domain is not the same
        ##Hence looking also for duplicated genes between domains
        gtf$isChr<-gtf$chr==duplicationArea_chr
        gtf$isBelowUpperLimit<-duplicationArea_end>=gtf$end 
        gtf$isAboveLowerLimit<-duplicationArea_start<=gtf$start
        
        ##Catching genes
        logicalInfoGenes<-gtf[,c("isChr","isBelowUpperLimit","isAboveLowerLimit")]
        isAGeneInside<-rowSums(logicalInfoGenes)==3##if ==3 fullfits all the requirements
        
        duplicatedGenes<-unique(c(duplicatedGenes, gtf[isAGeneInside,"geneSymbol"]))
        
      }else if(length(unique(domainsData$domainId))==1){
        ##So the affected Domain is the same,
        ##Hence duplicated genes only those in bkp1 right segment && bkp2 left segment
        ##BUT excluding those which are not shared, drawing intratad duplication to comprehend
        duplicatedGenes<-unique(intersect(genesInSegment$breakpoint_1_right_segment,genesInSegment$breakpoint_2_left_segment))
        
      } 
      
    }
  }
  
  
  ####Retrieve in an object the genes Position por all the possibly affected genes
  #allGENES<-unique(c(genesInDomain$domain_breakpoint_1, genesInDomain$domain_breakpoint_2,brokenGenes, genesInUncertainty))
  allGENES<-unique(c(genesInDomain$domain_breakpoint_1, genesInDomain$domain_breakpoint_2,brokenGenes, genesInUncertainty, deletedGenes, duplicatedGenes))
  
  genesPosition<-gtf[allGENES,c("chr","start","end","TSS","strand","geneSymbol","geneType")]
  
  
  ##rearreglar luego el objeto de salida
  ##If we retrieve the info for a gene, but it is not assigned to any domain, its because it is affected due to deletion or duplication ocurring between domains
  info_affectedGenes<-list("genesInDomain"=genesInDomain, 
                           "genesInSegment"=genesInSegment, 
                           "brokenGenes"=brokenGenes,
                           "deletedGenes"=deletedGenes,
                           "duplicatedGenes"=duplicatedGenes,
                           "genesInUncertainty"=genesInUncertainty,
                           "genesPosition"=genesPosition)
  return(info_affectedGenes)
  
}