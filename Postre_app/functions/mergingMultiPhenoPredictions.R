################################################################
## Merging heatmaps...from multiple pheno selected prediction
################################################################
mergingMultiPhenoPredictions<-function(patientResults, selectedPheno, allPostreAvailablePheno){
  
  #Eliminate from selectedPheno all those pheno not currently handled by POSTRE
  #To avoid provide meaningless information for them
  selectedPheno<-selectedPheno[selectedPheno %in% allPostreAvailablePheno]
  masterHtml<-""
  
  ##Metemos el indice si resultados para mas de 1 fenotipo

  if(length(selectedPheno)>1){
    ## Enrich this html
    masterHtml<-paste(masterHtml,
                      "<h1 id='ResultsOverview'><b>Results Index</b></h1>",
                      sep = "",
                      collapse = "")
    
    ##Meter los indices para el fenotipo correspondiente
    contentIndex<-"<ul>"
    for(patientPheno in selectedPheno){
      contentIndex<-paste(contentIndex,
                          "<li><a href='#PhenoRes_",
                          gsub(" ", "", patientPheno, fixed = TRUE), ##Removing white spaces to avoid problems
                          "'>See Results for ",
                          patientPheno,
                          "</a></li>",
                          sep = "",
                          collapse = "")
    }
    contentIndex<-paste(contentIndex,
                        "</ul>",
                        sep = "",
                        collapse = "")
    
    ##Me faltara meter estos, al masterHtml variable 
    masterHtml<-paste(masterHtml,
                      contentIndex,
                      sep = "",
                      collapse = "")
  }
  
  ##Attach heathampHtml for each of the individual phenotypes
  for(patientPheno in selectedPheno){

    masterHtml<-paste("<div class='wrapperMainSingleResults'>",#Div, to control/make disappear the output when re-submitting
                      masterHtml,
                      ##Adding a div to the heatmap to link from the index
                      "<div id='PhenoRes_",
                      gsub(" ", "", patientPheno, fixed = TRUE), ##Removing white spaces to avoid problems
                      "'>",
                      patientResults[[patientPheno]]$heatmapSummary,
                      "</div>",
                      "</div>",
                      sep = "",
                      collapse = "")
    
  }

  ##Assigning merged html to output variable
  patientResults$heatmapSummary<-masterHtml
  
  ##If script reached this point, no error rised
  patientResults$Status<-"OK"
  
  return(patientResults)
  
}