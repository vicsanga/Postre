################################################################
## Get uncertainty regions of the breakpoints
## Uncertainty region: space between the limits of a breakpoint
################################################################
getUncertainty<-function(breakpInfo){
  chr<-breakpInfo$chr
  coord<-breakpInfo$coord
  
  if(length(coord)==1){
    ##there is no uncertainty
    ##bp mapped at single bp resolution
    return(NULL)
  }else{
    ##so there are two coordinates for the bp
    if((coord[2]-coord[1])<=1){
      ###so the coord are either the same one(should not be due to previous unique,having two coord I mean)
      ##or they are consecutive bp, so there is no uncertainty region
      return(NULL)
    }else{
      ##it means there is at least 1 bp between the coord, so there is this uncertainty region
      uncertMatrix<-as.data.frame(matrix(data = NA, nrow = 1, ncol = 3))
      colnames(uncertMatrix)<-c("chr","start","end")
      
      uncertMatrix[1,"chr"]<-chr
      uncertMatrix[1,"start"]<-coord[1]
      uncertMatrix[1,"end"]<-coord[2]
      
      return(uncertMatrix)
    }
  }
}

