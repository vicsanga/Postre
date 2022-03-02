##################################
##Check same domain
##Function to check that the breakP either given by one bp or by a region fall in the same Domain
##and report the Domain
##################################
checkSameDomain<-function(breakpInfo,
                          tadMap,
                          regionsBetweenTADs)## coordinates of the regions between TADs
{
  chr<-breakpInfo$chr
  coord<-breakpInfo$coord
  ##combinar el tadMap y el boundary map
  ##y sacar el dominio que sea y el tipo a la vez
  ##Si es un punto reportar en que dominio cae, juntar boundaryYtadMap y coger columna "typo" para saber que es
  tadMap$domainType<-"TAD"
  regionsBetweenTADs$domainType<-"between_TADs"
  
  domains<-rbind(tadMap,regionsBetweenTADs)
  
  ###sorting not necessary since we do not iterate over rows,
  ##its a simultaneous vectorial search
  ##but lets do it just to see the effect, of how intercalated are the TAD and spaces between TADs
  domains<-domains[order(domains$chr, domains$start),]
  rownames(domains)<-1:nrow(domains)
  
  ##Si es una region, mirar si hay un dominio que contenga la region entera y reportarlo, de lo contrario decir que no se puede modelizar la variacion
  ##estructural pues la resolucion de almenos uno de los breakpoints cubre diferentes dominios
  counterPos<-0
  for(position in coord){
    # print(position)
    counterPos<-counterPos+1
    
    domains$isChr<-domains$chr==chr
    domains$isBelowUpperLimit<-domains$end>=position
    domains$isAboveLowerLimit<-domains$start<=position     
    
    ###is domain holding the breakpoint
    logicalInfoDomains<-domains[,c("isChr","isBelowUpperLimit","isAboveLowerLimit")]
    isADomainValid<-rowSums(logicalInfoDomains)==3##if ==3 fullfits all the requirements
    
    ##retrieve the domain coordinates
    if(counterPos==1){
      domainInfo<-domains[isADomainValid,c("chr","start","end","domainType")]
    }else{
      ##already one region reported, so do rbind
      domainInfo<-rbind(domainInfo, domains[isADomainValid,c("chr","start","end","domainType")])
    }
  }
  
  ##If both coord for breakp(if it has two) in the same domain
  ##just one row will be left
  domainInfo<-unique(domainInfo)
  
  ###once the loop is done
  if(nrow(domainInfo)==1){
    ##it means the whole breakpoint
    ##if it is a unique bp, of course, or if its delimited by a region
    ##falls entirely inside of a Domain
    # print("whole breakpoint Falls in a Domain")
    ##if just a domain let's return the info for just the first row
    infoToReturn<-list("sameDomain"=TRUE, "domainInfo"=domainInfo)
  }else{
    # print("the breakpoint region falls between 2 Domains")
    infoToReturn<-list("sameDomain"=FALSE, "domainInfo"=domainInfo)
  }
}