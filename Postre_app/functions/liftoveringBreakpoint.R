liftoveringBreakpoint<-function(chr, bpCoord){
  #Trick from https://genviz.org/module-01-intro/0001/06/02/liftoverTools/
  # "data/specificData_PerPhenotype/neurodevelopmental/MRD_neurodevelopmental.RData"
  

  ##Get breakpoint coordinates, may be 1 or 2: 1 if bp resolution, 2 if interval provided
  bp_allCoord<-unlist(strsplit(x=bpCoord, split = ",", fixed = TRUE))
  
  #Liftovering
  liftoveredInfo_chr<-character()
  liftoveredInfo_coords<-character()
  
  ## import the chain file
  #This can be done outside of here in a top environment to reduce the number of times this operation is run
  chainObject <- import.chain("data/hg38ToHg19.over.chain")
  
  for(genomicCoord in bp_allCoord){
    genomicCoord<-as.numeric(genomicCoord)
    
    # specify coordinates to liftover
    grObject <- GRanges(seqnames=chr, ranges=IRanges(start=genomicCoord, end=genomicCoord))  
    
    # ## import the chain file
    # chainObject <- import.chain("data/hg38ToHg19.over.chain")
    
    ## run liftOver
    results <- as.data.frame(liftOver(grObject, chainObject))
    
    ##Tracking liftoveredInfo
    liftoveredInfo_chr<-c(liftoveredInfo_chr, as.character(results$seqnames))
    
    if(length(unique(c(results$start, results$end)))>1){
      stop("ERROR, LIFTOVER PRODUCING MULTIPLE COORD FOR SINGLE COORD, UNCERTAIN HOW TO PROCEED")
    }
    
    liftoveredInfo_coords<-c(liftoveredInfo_coords, as.character(results$start)) ##Liftovering only one position so, start and end should be the same
  }

  # ##Results
  if(length(unique(liftoveredInfo_chr))>1){
    stop("ERROR, LIFTOVER PRODUCING MULTIPLE CHROMOSOMES, UNCERTAIN HOW TO PROCEED")
  }
  
  ##Formatted info
  
  if(length(liftoveredInfo_coords)>1){
    ##So multiple coords, uncertain region
    liftoveredInfo_coords<-paste0(liftoveredInfo_coords,
                                 sep="",
                                 collapse = ",")
  }
  
  return(list("chr"=as.character(unique(liftoveredInfo_chr)),
        "coord"=as.character(liftoveredInfo_coords))
        )
  
}