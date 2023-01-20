######################
## Liftovering Genomic Coordinates
source(file = "functions/liftoveringBreakpoint.R",local = TRUE)

processing_GenomicCoordinates<-function(patientInfo){

  ##Currently the only option, liftover of hg38 to hg19, so, no need to additional checks on that
  
  ##Since each SV breakpoint may present 1 or 2 coordinates (if interval provided) each SV breakpoint must be
  ##Treated independently
  
  ##Processing breakp1
  breakp1_result<-liftoveringBreakpoint(chr = patientInfo$chr_Break1, 
                                        bpCoord = patientInfo$coord_Break1 
                                        )
  
  breakp2_result<-liftoveringBreakpoint(chr = patientInfo$chr_Break2, 
                                        bpCoord = patientInfo$coord_Break2 
                                        )  

  ##Reintroducing formatted coordinates 
  patientInfo$chr_Break1<-breakp1_result$chr
  patientInfo$coord_Break1<-breakp1_result$coord
  
  patientInfo$chr_Break2<-breakp2_result$chr
  patientInfo$coord_Break2<-breakp2_result$coord
  
  return(patientInfo)
}