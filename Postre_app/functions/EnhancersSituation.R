#####################################################
## Analyse enhancers scenario in the affected regions
#####################################################
#####################################################
## Function to get the number of ENHANCERS (DON'T KNOW IF SOURCE ALSO AND COORDS) located at the affected Regions 
## Atending to both the domains, and the breakpoint segments
#####################################################

###If an enhancer is broken we do not count it, We add it in the enhancers broken section.
###we put it as broken. Ya que esto creo que podria afectar a su funcion tb

affectedEnhancers<-function(regionsAltered,enhMap){
  ##regionsAltered for info_affectedRegions object
  domainsData<-regionsAltered$domainsAffected
  segmentsData<-regionsAltered$segmentsBreakP
  uncertaintyData<-regionsAltered$uncertaintyRegions
  
  #############################################
  #####Look for ENHANCERS contained in the domains
  ## Focus on ENHANCERS whole inside the domain (if Alvaro do not like this, put the eye on whichever side)
  
  enhancersInDomain<-list()
  #Get a matrix per Domain, with a row per enhancer, with its coordinates and source
  for(dominio in rownames(domainsData)){

    targetDomain_Chr<-domainsData[dominio,"chr"]
    targetDomain_Start<-domainsData[dominio,"start"]
    targetDomain_End<-domainsData[dominio,"end"]

    enhMap$isChr<-enhMap$chr==targetDomain_Chr
    enhMap$isBelowUpperLimit<-targetDomain_End>=enhMap$end ##enhMap$end is the end of the enhancer
    enhMap$isAboveLowerLimit<-targetDomain_Start<=enhMap$start

    ##Is TAD holding the whole enh?
    logicalInfoEnh<-enhMap[,c("isChr","isBelowUpperLimit","isAboveLowerLimit")]
    isAEnhInside<-rowSums(logicalInfoEnh)==3##if ==3 fullfits all the requirements
  
    if(sum(isAEnhInside)>=1){
      ##So we have at least an enhancer
      resMatrix<-enhMap[isAEnhInside,c("chr","start","end","source","acetilation")]

    }else{
      # resMatrix<-NULL ##Problem. If there are no, will not be added
      ##change by matrix with no rows so it will be added, by name
      
      resMatrix<-as.data.frame(matrix(data = NA,nrow = 0,##nrow =0 we want it EMPTY
                                      ncol = 5))
      colnames(resMatrix)<-c("chr","start","end","source","acetilation")
    }
    
    enhancersInDomain[[dominio]]<-resMatrix
  }
  
  
  ######################################################
  ## Look for enhancers ENTIRELY contained in each breakpoint segment
  enhancersInSegment<-list()
  for(segment in rownames(segmentsData)){
    
    targetSegment_Chr<-segmentsData[segment,"chr"]
    targetSegment_Start<-segmentsData[segment,"start"]
    targetSegment_End<-segmentsData[segment,"end"]
    
    enhMap$isChr<-enhMap$chr==targetSegment_Chr
    enhMap$isBelowUpperLimit<-targetSegment_End>=enhMap$end ##enhMap$end is the end of the enh
    enhMap$isAboveLowerLimit<-targetSegment_Start<=enhMap$start
    
    ##Is segment holding the whole enh?
    logicalInfoEnh<-enhMap[,c("isChr","isBelowUpperLimit","isAboveLowerLimit")]
    isAEnhInside<-rowSums(logicalInfoEnh)==3##if ==3 fullfits all the requirements
    
    if(sum(isAEnhInside)>=1){
      ##So we have at least an enhancer
      resMatrix<-enhMap[isAEnhInside,c("chr","start","end","source","acetilation")]

    }else{
      resMatrix<-NULL
    }
    
    enhancersInSegment[[segment]]<-resMatrix
  }
  
  ############################################################################
  ## Look for enhancers entirely contained in the 
  ## uncertainty regions. (breakpoint given by a region and not at single bp)
  ## If they are not entirely contained it is because
  ## they cross any of the breakpoint limits, 
  ## so they will be considered broken
  #############################################################################
  
  enhancersInUncertainty<-NULL##if there are we will keep them as a list of matrixes, one per breakpoint
  for(breakp in names(uncertaintyData)){
    info_uncert<-uncertaintyData[[breakp]]
    if(is.null(info_uncert)==FALSE){
      ##so there is an uncertainty region
      ##Lets get the coordinates for this segment (uncertainty segment)
      ##If there is data we have a 1 row matrix with the region coord
      targetSegment_Chr<-info_uncert[1,"chr"]
      targetSegment_Start<-info_uncert[1,"start"]
      targetSegment_End<-info_uncert[1,"end"]
      
      enhMap$isChr<-enhMap$chr==targetSegment_Chr
      enhMap$isBelowUpperLimit<-targetSegment_End>=enhMap$end ##enhMap$end is the end of the enh
      enhMap$isAboveLowerLimit<-targetSegment_Start<=enhMap$start
      
      ##Is segment holding the whole enh?
      logicalInfoEnh<-enhMap[,c("isChr","isBelowUpperLimit","isAboveLowerLimit")]
      isAEnhInside<-rowSums(logicalInfoEnh)==3##if ==3 fullfits all the requirements
      
      if(sum(isAEnhInside)>=1){
        ##So we have at least an enhancer
        resMatrix<-enhMap[isAEnhInside,c("chr","start","end","source","acetilation")]
        
      }else{
        ##Even there is uncertainty region, there is no enhancer inside
        resMatrix<-NULL
      }
    }else{
      ##So no uncertainty region
      resMatrix<-NULL
    }
    if(is.null(resMatrix)==FALSE){
      ##so there is an enhancer at least in the uncertainty region
      ##captured by the matrix
      enhancersInUncertainty[[breakp]]<-resMatrix
    }
  }
  
  #############################################################
  ## ¿Any ENHANCER has its sequence broken by the breakpoint?
  ## LeftSegment---BREAKPOINT---RigthSegment
  ## A gene is considered to be broken if the end of any of the left segments
  ## or the start of any of the right segments falls in the Enhancer region (between its start or end)
  ## brokenEnhancers
  
  ##Broken Enhancers considered per Domain
  ##If we then want to know which segment look for it
  temp_brokenEnhancers<-list()
  
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
    
    enhMap$isChr<-enhMap$chr==targetSegment_Chr
    enhMap$isBreakBeforeEnd<-puntoRotura<=enhMap$end ##enhMap$end is the end of the gene
    enhMap$isBreakAfterStart<-puntoRotura>=enhMap$start
    
    ##
    ##Is enhancer broken?
    logicalInfoEnh<-enhMap[,c("isChr","isBreakBeforeEnd","isBreakAfterStart")]
    isAEnhBroken<-rowSums(logicalInfoEnh)==3##if ==3 fullfits all the requirements
    
    if(sum(isAEnhBroken)>=1){
      ##So we have at least an enhancer
      resMatrix<-enhMap[isAEnhBroken,c("chr","start","end","source","acetilation")]

    }else{
      resMatrix<-NULL
    }
    
    temp_brokenEnhancers[[segment]]<-resMatrix
  }
  
  ########################################
  ##Get unique broken enhancers per domain
  #########################################
  ##there can be repeated broken enhancers if we have a single coordinate or two super close for the breakpoint, 
  ##in this case the broken enhancer will appear two times, so we do the unique
  
  domain1_BrokenEnh<-unique(rbind(temp_brokenEnhancers[["breakpoint_1_left_segment"]],
                          temp_brokenEnhancers[["breakpoint_1_right_segment"]]))
  
  domain2_BrokenEnh<-unique(rbind(temp_brokenEnhancers[["breakpoint_2_left_segment"]],
                            temp_brokenEnhancers[["breakpoint_2_right_segment"]]))
  
  brokenEnhancers<-list("domain_breakpoint_1"=domain1_BrokenEnh,
                        "domain_breakpoint_2"=domain2_BrokenEnh)
  
  ################
  ## Returning data
  info_affectedEnhancers<-list("enhancersInDomain"=enhancersInDomain,
                               "enhancersInSegment"=enhancersInSegment,
                               "brokenEnhancers"=brokenEnhancers,
                               "enhancersInUncertainty"=enhancersInUncertainty)
  
}