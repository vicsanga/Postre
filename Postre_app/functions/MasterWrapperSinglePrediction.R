############################################
## Master Wrapper Single Prediction 
## This creates all the required plots, htmls. etc.
############################################
####Loading Required Functions
###############################################
##To tho the main ranking//scoring
source(file = "functions/Master_Scoring_Function.R", local = TRUE)

## To generate outputs from prediction results
source(file = "functions/Gene_textReport.R",local = TRUE)
source(file = "functions/ucsc_view.R",local = TRUE)
source(file = "functions/heatmap_SummaryResults.R",local = TRUE)
source(file = "functions/graphicalSummary_generation.R",local = TRUE)
source(file = "functions/graphSummary_imgLoader.R",local = TRUE)
source(file = "functions/plots_RegulatoryEnvironmentChanges.R",local = TRUE)
source(file = "functions/plots_geneSusceptibility.R",local = TRUE)


masterWrapperSinglePrediction<-function(patientInfo,minScore, highScore, runMode, user_tadMapInfo, MultiDataList){
  
  ##Initial Prediction
  patientResults<-master_scoring_function(patientInfo = patientInfo, runMode = runMode, user_tadMapInfo = user_tadMapInfo, MultiDataList = MultiDataList)

  #############################################################
  ## Check whether there are genes associated with the SV
  ## If there are NOT, do nothing, afterwards will be provided html saying: "NO GENES FOUND ASSOCIATED"
  ## If there are, provide normal graphics output
  
  if(patientResults$Status == "OK, but NO genes associated with SV"){
    ##Just do nothing
  }else if(patientResults$Status == "OK"){
    
    ##Compute the normal and expected graphics and reports
    ##Generating graphical summary for the interesting genes
    #Interesting genes are tracked on the returned object
    patientResults<-graphicalSummary_generation(patientResults = patientResults,
                                                minPathogenicScore = minScore)
    
    ##Generating EnhancerLanscapeChanges Barplots
    plots_regulatoryEnvironmentChanges(patientResults = patientResults)
    
    ##Generating geneSusceptibility to GOF|LOF section plots
    plots_geneSusceptibility(patientResults = patientResults)
    
    ###################################
    ## CREATING MAIN RESULTS PAGE
    ###Creating Heatmap summary result
    
    patientResults$heatmapSummary<-heatmap_summaryResults(patientResults = patientResults,
                                                          minRequiredScore = minScore,
                                                          highScore = highScore )
    
    #######################################################################################
    ## Error Control, specifying it again "confirming"
    ## Since errors can arise on the additional steps
    #######################################################################################
    ##If script reached this point, no error rised
    patientResults$Status<-"OK"
    
  }
  
  ######################
  return(patientResults)
}