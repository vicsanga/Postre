#####################################################################################################
## Function to retrieve the Affected Regions. And the 4 segments to play with
## eg: ---TAD---Region1--breakpoint1--Region2---TAD     ---TAD---Region1--breakpoint2--Region2---TAD           
#####################################################################################################

##function to retrieve the domain where the breakpoint is located
##and in case that the breakpoint is given by a region insted of by a point that the whole region falls in the same domain

affectedRegions<-function(dataPatient, mapTads, regionsBetweenTADs){
  # print("AND FOR BREAKPOINTS OCCURRING IN DIFFERENT TADS")
  # print("------------------------------------------------")
  
  breakPoint_1_chr<-dataPatient$chr_Break1
  breakpoint_1_coord<-unique(as.integer(unlist(strsplit(dataPatient$coord_Break1,split = ",", fixed = TRUE))))
  
  breakPoint_2_chr<-dataPatient$chr_Break2
  breakpoint_2_coord<-unique(as.integer(unlist(strsplit(dataPatient$coord_Break2,split = ",", fixed = TRUE))))
  
  breakpData<-list("breakpoint_1"=list("chr"=breakPoint_1_chr, "coord"=breakpoint_1_coord),
                   "breakpoint_2"=list("chr"=breakPoint_2_chr, "coord"=breakpoint_2_coord))
  
  ##To retrieve the TADs
  ##In mind to have two options, TAD, betweenTAD
  
  # source("functions/CheckSameDomain.R", local = TRUE)
  
  ##function to check that both coordinates for a breakpoint fall in the same Domain, and report it
  ##if there is just one coordinate, the result will be TRUE, and sure a domain will be reported
  ##because if it has two coord, we have to check if both fall in the same Domain
  ##on the contrary we won't process that SV, for that TAD map
  ##Prediction being really complicated in that case, under our model of domains
  
  ##For breakpoint 1
  bp1_domainMapping<-checkSameDomain(breakpInfo = breakpData$breakpoint_1, 
                                     tadMap = mapTads,
                                     regionsBetweenTADs = regionsBetweenTADs)
  ##For breakpoint 2
  bp2_domainMapping<-checkSameDomain(breakpInfo = breakpData$breakpoint_2,
                                     tadMap = mapTads,
                                     regionsBetweenTADs = regionsBetweenTADs)
  
  ##To store the domains location per breakpoint
  ##In mind to have three options, TAD && between_TAD
  if(bp1_domainMapping$sameDomain==TRUE && bp2_domainMapping$sameDomain==TRUE){
    domainsAffected<-rbind(bp1_domainMapping$domainInfo, bp2_domainMapping$domainInfo)
    rownames(domainsAffected)<-c("domain_breakpoint_1","domain_breakpoint_2")
  }else{
    "PRINT THERE IS AT LEAST A BREAKPOINT WHICH FALLS BETWEEN DOMAINS"
    "DO NOTHING WITH THIS SV"
    "errBreakp"
    stop("The breakpoint falls between Domains Limits. Not clear how to proceed.")
  }

  ############################
  ##Let's deal with both breakpoints in the same Domain
  ##Well we do not have to remove this, right?? unless it is for an inversion
  ##Becuase if it is a duplication or a deletion...
  ##en matrix Ev domain Id, and if it is an inversion, check that is in different TADs
  ##Up to now apparently won't have an effect
  
  
  #############################
  ##Let's get segments coordinates
  segmentsBreakP<-as.data.frame(matrix(data = NA, nrow = 4, ncol = 3))
  colnames(segmentsBreakP)<-c("chr","start","end")##type of domain in domains info
  rownames(segmentsBreakP)<-c("breakpoint_1_left_segment","breakpoint_1_right_segment",
                              "breakpoint_2_left_segment","breakpoint_2_right_segment")
  
  #####filling the info
  #######################
  ##For breakpoint1
  ##chr
  segmentsBreakP[c("breakpoint_1_left_segment","breakpoint_1_right_segment"),"chr"]<-breakPoint_1_chr
  
  ##coord
  #Left
  segmentsBreakP["breakpoint_1_left_segment","start"]<-domainsAffected["domain_breakpoint_1","start"]
  segmentsBreakP["breakpoint_1_left_segment","end"]<-min(breakpoint_1_coord)
  ##Right
  segmentsBreakP["breakpoint_1_right_segment","end"]<-domainsAffected["domain_breakpoint_1","end"]
  segmentsBreakP["breakpoint_1_right_segment","start"]<-max(breakpoint_1_coord)
  
  ###################
  ##For breakpoint2
  ##chr
  segmentsBreakP[c("breakpoint_2_left_segment","breakpoint_2_right_segment"),"chr"]<-breakPoint_2_chr
  
  ##coord
  #Left
  segmentsBreakP["breakpoint_2_left_segment","start"]<-domainsAffected["domain_breakpoint_2","start"]
  segmentsBreakP["breakpoint_2_left_segment","end"]<-min(breakpoint_2_coord)
  ##Right
  segmentsBreakP["breakpoint_2_right_segment","end"]<-domainsAffected["domain_breakpoint_2","end"]
  segmentsBreakP["breakpoint_2_right_segment","start"]<-max(breakpoint_2_coord)

  #############################
  ##Let's get uncertainty regions
  ######UP TO HERE SEGUIR POR AQUI
  # source("functions/GetUncertainty.R",local = TRUE)
  
  bp1_uncertaintyRegion<-getUncertainty(breakpInfo = breakpData$breakpoint_1)##Return the coord or NULL if there isn't
  bp2_uncertaintyRegion<-getUncertainty(breakpInfo = breakpData$breakpoint_2)
  
  #rbind(bp1_uncertaintyRegion, bp2_uncertaintyRegion)
  
  ##Building uncertainty region object
  uncertaintyRegions<-list("breapkoint1_uncertainty" = bp1_uncertaintyRegion, "breakpoint2_uncertainty" = bp2_uncertaintyRegion)
  ##if there is uncertainty it is given in a 1 row matrix, with the coord of the region
  ##on the contrary NULL is returned
  
  
  ############################
  # Adding affected domains id
  ############################
  domainsAffected$domainId<-paste(domainsAffected$chr,"_",
                                 domainsAffected$start,"_",
                                 domainsAffected$end,
                                 sep="")
  
  #############################
  # Checking if SV is IntraTAD or InterTAD
  # Variable called SV landing
  # If both breakpoints fall in the same TAD
  
  if(domainsAffected["domain_breakpoint_1","domainId"] == domainsAffected["domain_breakpoint_2","domainId"]){
    ##Intra TAD SV
    SV_landing<-"IntraTAD"
    
  }else if(domainsAffected["domain_breakpoint_1","domainId"] != domainsAffected["domain_breakpoint_2","domainId"]){
    ##Inter TAD SV
    ##Because each breakpoint falls in a different Domain
    SV_landing<-"InterTAD"
  }
  
  
  #########################
  ## Objects to be returned
  #########################
  info_affectedRegions<-list("domainsAffected" = domainsAffected, "segmentsBreakP" = segmentsBreakP, 
                             "uncertaintyRegions" = uncertaintyRegions, "SV_landing" = SV_landing)
  ##add uncertainty regions in case there are, if there aren't it is NULL
  return(info_affectedRegions)
  
}








