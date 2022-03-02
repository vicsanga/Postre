######################################################
## Function developed to check patient input features
## And liftover coordinates from hg38 to hg19 if necessary
#####################################################
# https://master.bioconductor.org/packages/release/workflows/html/liftOver.html
# https://master.bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html
# library(liftOver)
# library(GenomicRanges)
# library(rtracklayer)

checkPatientFeatures<-function(patientInfo){
  ##To track if any error rised
  
  ###CHECK breakp1 pos menor que breakp2???
  ##dentro de pos para un breakp, a la izquierda de la coma el nro menor
  ##y a la derecha el mayor
  
  ##Y el breakpoint 1 menor que el breakpoint2
  # unique(as.integer(unlist(strsplit(dataPatient$coord_Break1,split = ",", fixed = TRUE))))
  # Get Breakp coord as numeric, with the left coord enough for the purpose of the order Testing (check if the smaller is the first one)
  break1_coord<-as.integer(unlist(strsplit(patientInfo$coord_Break1,split = ",", fixed = TRUE)))[1]
  break2_coord<-as.integer(unlist(strsplit(patientInfo$coord_Break2,split = ",", fixed = TRUE)))[1]  
  
  if(patientInfo$TypeSV %in% c("Inversion","Deletion","Duplication")){
    ##So SV happening on the same chromosome
    if(patientInfo$chr_Break1 != patientInfo$chr_Break2){
      # print("Error: for an inversion both breakpoints expected in the same chromosome")
      stop("For an inversion both breakpoints expected in the same chromosome")
    }
    
    if(break1_coord > break2_coord){
      ##error, because the break1 coord must be smaller, for SV happening on the same chr
      stop("The break1 coordinates must be smaller than break2 coordinates, for SV happening on the same chr. Check them")
    }
  }
  
  ################################################################
  ## If code reaches this point, no error rised
  ## Converting genome coordinates to hg19 if they come from hg38
  ################################################################
  # if(patientInfo$refGenome == "hg38"){
  #   ##We need to liftover the coordinates to hg19
  #   #browser()
  #   
  #   ##Get breakp1 start
  #   ##Get breakp1 end
  #   ##Get breakp1 chr
  #   chr_breakp1<-patientInfo$chr_Break1
  #   breakp1_start<-as.integer(unlist(strsplit(patientInfo$coord_Break1,split = ",", fixed = TRUE)))[1]
  #   breakp1_end<-as.integer(unlist(strsplit(patientInfo$coord_Break1,split = ",", fixed = TRUE)))[2]
  #   
  #   
  #   ##Get breakp2 start
  #   ##Get breakp2 end
  #   ##Get breakp2 chr
  #   chr_breakp2<-patientInfo$chr_Break2
  #   breakp2_start<-as.integer(unlist(strsplit(patientInfo$coord_Break2,split = ",", fixed = TRUE)))[1]
  #   breakp2_end<-as.integer(unlist(strsplit(patientInfo$coord_Break2,split = ",", fixed = TRUE)))[2]
  #   
  #   ##If breakpoint at single base pair resolution it means only start coord provided, this will imply that value for end is NA
  #   ##So repeat again the same position for the end
  #   if(is.na(breakp1_end)){
  #     breakp1_end <- breakp1_start
  #   }
  #   
  #   if(is.na(breakp2_end)){
  #     breakp2_end <- breakp2_start
  #   }
  #   
  #   ##Prepare data frame
  #   ##Test from tutorial
  #   # https://master.bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html
  #   ##Si es bp resolution, start y end el mismo
  #   df_coord_hg38<-data.frame(chr=c(chr_breakp1, chr_breakp2),
  #                             start=c(breakp1_start, breakp2_start),
  #                             end=c(breakp1_end, breakp2_end))
  #   
  #   ##Creating GRanges object
  #   granges_hg38<-makeGRangesFromDataFrame(df_coord_hg38, seqnames.field = "chr",  start.field = "start", end.field = "end")
  #   
  #   ##Getting hg38 to hg19 chain file object, provided by the UCSC
  #   # The transformation to hg19 coordinates is defined by a chain file provided by UCSC. rtracklayer::import.chain will bring the data into R.
  #   path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
  #   ch = import.chain(path)
  #   ##Option. ESTE OBJETO GUARDARLO EN LA INFO DEL APP, por si depende de internet para sacarlo...
  #   ##rename coordinates, specifying UCSC style
  #   seqlevelsStyle(granges_hg38) = "UCSC"
  #   
  #   ################################
  #   ##Getting Coordinates in hg19
  #   ################################
  #   granges_hg19<-liftOver(x=granges_hg38, chain = ch)
  #   ##Now convert granges to dataframe
  #   df_coord_hg19<-data.frame(granges_hg19)
  #   df_coord_hg19<-df_coord_hg19[,c("seqnames","start","end")]
  #   colnames(df_coord_hg19)<-c("chr","start","end")  
  # 
  #   ##If any coordinate cannot be converted rise a warning, or error
  #   if(nrow(df_coord_hg19)!=2){
  #     stop("Not possible to convert all coordiantes to hg19")
  #   }
  #   
  #   rownames(df_coord_hg19)<-c("breakp1","breakp2")
  #   
  #   #################################################################################
  #   ##Now, formatting patientInfo to hg19 coordinates according to liftover results
  #   #################################################################################
  #   ##Correcting refGenome
  #   patientInfo$refGenome<-"hg19"
  #   
  #   ##Adding chr data
  #   patientInfo$chr_Break1<-as.character(df_coord_hg19["breakp1","chr"])
  #   patientInfo$chr_Break2<-as.character(df_coord_hg19["breakp2","chr"])    
  #   
  #   ##Adding coord breakp1 data
  #   breakp1_hg19_coord_info<-c(df_coord_hg19["breakp1","start"],
  #                         df_coord_hg19["breakp1","end"])
  #   
  #   if(length(unique(breakp1_hg19_coord_info))==1){
  #     patientInfo$coord_Break1<-trimws(paste(unique(breakp1_hg19_coord_info),
  #                                     ",",
  #                                     sep="",
  #                                     collapse = ""))
  #     
  #   }else{
  #     ##So it will be equal to two coordinates
  #     patientInfo$coord_Break1<-trimws(paste(df_coord_hg19["breakp1","start"],
  #                                            ",",
  #                                            df_coord_hg19["breakp1","end"],
  #                                            sep="",
  #                                            collapse = ""))
  #   }
  #   
  #   
  #   ##Adding coord breakp2 data  
  #   breakp2_hg19_coord_info<-c(df_coord_hg19["breakp2","start"],
  #                              df_coord_hg19["breakp2","end"])
  #   
  #   if(length(unique(breakp2_hg19_coord_info))==1){
  #     patientInfo$coord_Break2<-trimws(paste(unique(breakp2_hg19_coord_info),
  #                                            ",",
  #                                            sep="",
  #                                            collapse = ""))
  #     
  #   }else{
  #     ##So it will be equal to two coordinates
  #     patientInfo$coord_Break2<-trimws(paste(df_coord_hg19["breakp2","start"],
  #                                            ",",
  #                                            df_coord_hg19["breakp2","end"],
  #                                            sep="",
  #                                            collapse = ""))
  #   }
  #   
  #   
  #   
  # }
  
  #######################
  ##Return patient info, same as it was, if hg19, or converted to hg19 if it was hg38
  #######################
  return(patientInfo)
  
}