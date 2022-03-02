############
##Secondary Function Used to avoid code repetition
addingInfoToMatrix<-function(gene, targetMatrix, phase, patientId){
  ######################################
  ##First adding to anyMechanism matrix
  if(gene %in% rownames(targetMatrix)){
    ##Adding 1 to its corresponding cell if it is not an NA. In that case assign to 1
    if(is.na(targetMatrix[gene,phase])){
      targetMatrix[gene,phase]<-1
    }else{
      ##Adding one to the already present counter
      targetMatrix[gene,phase] <- targetMatrix[gene,phase] + 1
    }
    
  }else{
    ##Creating gene entry, hence setting it to 1
    targetMatrix[gene,phase]<-1
  }
  
  ###########################
  ##Adding patient Id Info
  patientCellInfo<-targetMatrix[gene,"patients"]
  if(any(is.na(patientCellInfo))){
    ##Adding First Patient
    targetMatrix[gene,"patients"]<-patientId
  }else{
    ##Appending to already present Patients
    ##If patient not previously added, due to gene appearing in different phases

    if(!grepl(pattern = patientId, patientCellInfo)){
      #So info patient not previously added
      patientCellInfo<-c(patientCellInfo, patientId)
    }
    ##to add to the matrix we need to paste the elements into a single one
    ##to avoid errors
    targetMatrix[gene,"patients"]<-paste(unique(patientCellInfo), collapse = ",")
  }
  return(targetMatrix)
}