#######################################################
## Report Gene Susceptibility to Condition [LOF|GOF]
#######################################################
report_geneSusceptibilityToCondition<-function(reportUnit, patientResults){

  ###I need to know whether the considered gene is selected by means of a 
  ## LOF && GOF mechanism || LONG-RANGE--DIRECT EFFECT

  targetGene<-unlist(strsplit(x = reportUnit, split = "_", fixed = TRUE))[1]
  targetMech<-unlist(strsplit(x = reportUnit, split = "_", fixed = TRUE))[2]
  targetPhase<-unlist(strsplit(x = reportUnit, split = "_", fixed = TRUE))[3]
  
  #This way of approaching the GOF-LOF info, deprecated
  # pathoMechanism<-patientResults$resultsPerPhase[[targetPhase]][[targetPhase]][targetGene,]$Mechanism
  # pathoMechanism<-unlist(strsplit(x = pathoMechanism,
  #                                 split = ":",
  #                                 fixed = TRUE))[2]##GOF or LOF
  pathoMechanism<-targetMech
  
  ##To discriminate in dupl deletions, do something like, pathomechanism_Impact by long Range, or by direct effect
  ##As a gene duplication can trigger the disease due to neoTad and hence, dy a long-range mechanism
  targetRow_resultsPhase<-paste(targetGene,"--",targetMech, sep="")##In case interest on parsing resultsPhase
  pathomechanism_ImpactOverGene<-patientResults$resultsPerPhase[[targetPhase]][[targetPhase]][targetRow_resultsPhase,]$GeneImpact
  
  ######################
  ##OPENING Section
  geneSuscept<-"<div class='geneSusceptibilityToConditionSection'>"
  
  ####################################
  ## Adding Images//InfoGraphic Part
  ####################################
  
  ## TO DO
  ##Here it will be a custom selection regarding if it is lof or gof, or the phase, and gene conditions
  ##if expressed but not HI... maybe we don't load the HI part?
  
  ##get file names in graphical summaries
  fileNames<-list.files(path = "www/graphicalSummaries/")
  
  ##If working on phaseFree, not considered Expression Data
  if(targetPhase != "phaseFree"){
    ##Check if expPlot exists before adding
    expPlotName<-paste(reportUnit,"_", patientResults$job_UniCode, "_expressionLevelPlot.png", sep="")
    if(expPlotName %in% fileNames){
      geneSuscept<-paste(geneSuscept,
                         "<div class='expressionLevel'>",
                         paste0("<img class='geneSusceptImgs' src='graphicalSummaries/",expPlotName , "'",
                                ">"),
                         "</div>",
                         sep="")  
    }
  }
  
  ######################################################################################
  ###Adding Secondary images, Gene Dosage and polybcomb score plots
  ######################################################################################
  ###############################
  ##Ading Dosage Sensitive Plot
  ##Check if HI plot exists before adding
  
  ##Check if there is hi score plot for the target gene
  hiPlotName<-paste(targetGene,"_", patientResults$job_UniCode, "_HI_plot.png", sep="")
  if(hiPlotName %in% fileNames){
    
    if(targetPhase != "phaseFree"){
      geneSuscept<-paste(geneSuscept,
                         "<div class='dosageSensitivePlot'>",
                         paste0("<img class='hiPlotImg geneSusceptImgs' src='graphicalSummaries/", hiPlotName, "'",
                                ">"),
                         "</div>",
                         sep="") 
    }else{
      ##So PHASE FREE
      ##Here we are not adding the Expression Level plot, expression not considered, so relocate the plot to the left
      ##To avoid leaving a white space at the left part
      geneSuscept<-paste(geneSuscept,
                         "<div class='dosageSensitivePlot_whenPrimaryPlotNotPresent'>",
                         paste0("<img class='hiPlotImg geneSusceptImgs' src='graphicalSummaries/", hiPlotName, "'",
                                ">"),
                         "</div>",
                         sep="")
    }
  }
  
  ################################  
  ##Adding Polycomb GOF image
  ##id='polycGene_Image'
  
  ##Check if there is polycomb score plot for the target gene
  polycPlotName<-paste(targetGene,"_", patientResults$job_UniCode, "_polyCscore_plot.png", sep="")
  
  if(polycPlotName %in% fileNames){
    
    geneSuscept<-paste(geneSuscept,
                       "<div class='polyCombPlot'>",
                       paste0("<img class='polyCplotImg geneSusceptImgs' src='graphicalSummaries/", polycPlotName, "'",
                              ">"),
                       "</div>",
                       sep="") 
  }
  
  
  
  
  #####################################
  ## Adding Explanation 
  #####################################
  ##Modify regarding if it is for LOF or GOF
  #targetGene
  #pathoMechanism ##GOF or LOF
  #pathomechanism_ImpactOverGene ##long-range or direct
  
  textExplanation<-"<br><br>"
  
  if((pathoMechanism=="LOF") && (targetPhase!="phaseFree")){
    ##Be sure the gene is on the TAD map we are assessing, so screen them
    
    for(nPhase in 1:length(patientResults$resultsPerPhase_secondaryInfo[[targetPhase]])){
      
      evaluationMatrix<-patientResults$resultsPerPhase_secondaryInfo[[targetPhase]][[nPhase]]$matrixesGenesEvaluation
      
      if(targetGene %in% rownames(evaluationMatrix)){
        ##So the gene appears in the currently used TAD map
        ##Hence we have the info we need from it
        targetExpression<-evaluationMatrix[targetGene,
                                           paste("FPKM_",
                                                 targetPhase,
                                                 sep="")]
        
        targetHI_score<-evaluationMatrix[targetGene,
                                         c("nature_HI_score",
                                           "huang_HI_score","clinGene_HI_score")] ##-1 when no score assigned
        
        targetHI_score<-max(targetHI_score)##if it is -1, it is because no HI score assigned
        
        
        polycScore<-evaluationMatrix[targetGene,"polyComb_score"]
        
        break
      }
    }
    
    
    textExplanation<-paste(textExplanation,
                           # "<b style='font-size:20px;'>Wild Type (WT) vs Patient rearranged locus comparison </b>",
                           # "<br><br>", 
                           targetGene,
                           " presents an expression of ",
                           round2(x = targetExpression,
                                  digits = 1),
                           " fpkm at the ", targetPhase,
                           " stage. ",
                           ##"<br><br>",
                           sep="")
    
    #########################################################################################
    ##Adding Haploinsufficiency score info (if there is)
    ##If there is not we will talk about the gene expression pattern if tissue specific
    ##########################################################################################
    
    if(targetHI_score>-1){
      
      ##So an HI score provided
      textExplanation<-paste(textExplanation,
                             "The Dosage Sensitivity Score (DS) obtained for ",
                             targetGene,
                             " is : ",
                             round2(targetHI_score, digits=2),
                             ". This score must be interpreted on a 0-1 scale (1 means that for the candidate gene, producing more or less protein than normal is likely pathogenic, 0 means that deviations from the normal amount of generated protein are not detrimental).",
                             # "<br><br>",
                             sep="") 
    }
    
    ###############################
    ##Polycomb Score Information
    ###############################
    textExplanation<-paste(textExplanation,
                           " The Polycomb Score obtained for ",
                           targetGene,
                           " is: ",
                           round2(polycScore,digits = 2),
                           ". This score must be interpreted on a 0-1 scale (1 means strong evidence for being a H3K27me3-polycomb based regulatory gene, 0 implies the opposite).",
                           "<br><br>",
                           sep="") 
    
    
    
    textExplanation<-paste(textExplanation,
                           "To be a strong candidate for 'Loss of Function', at a particular stage, the candidate gene must be expressed. ",
                           "In addition, it is important to know the dosage sensitivity for the target gene. Because, having more or less protein than normal can be pathogenic for some but not all genes, Haploinsufficiency metrics are currently used as proxy for this aim. ",
                           "Furthermore, we also find relevant to study whether the candidate gene is a polycomb one (genes covered by broad domains of H3K27me3 when inactive and with tissue specific expression patterns), since these genes tend to regulate cellular identity and altering their expression can lead to strong phenotypic alterations. Check our manuscript to know more about how these metrics have been obtained. ",
                           "<br><br>",
                           sep="")

    
    
  }else if((pathoMechanism=="GOF") && (targetPhase!="phaseFree")){
    ##Be sure the gene is on the TAD map we are assessing, so screen them
    
    for(nPhase in 1:length(patientResults$resultsPerPhase_secondaryInfo[[targetPhase]])){
      
      evaluationMatrix<-patientResults$resultsPerPhase_secondaryInfo[[targetPhase]][[nPhase]]$matrixesGenesEvaluation
      
      if(targetGene %in% rownames(evaluationMatrix)){
        ##So the gene appears in the currently used TAD map
        ##Hence we have the info we need from it
        targetExpression<-evaluationMatrix[targetGene,
                                           paste("FPKM_",
                                                 targetPhase,
                                                 sep="")]
        
        targetHI_score<-evaluationMatrix[targetGene,
                                         c("nature_HI_score",
                                           "huang_HI_score","clinGene_HI_score")] ##-1 when no score assigned
        
        targetHI_score<-max(targetHI_score)##if it is -1, it is because no HI score assigned
        
        polycScore<-evaluationMatrix[targetGene,"polyComb_score"]
        
        break
      }
    }
    
    
    
    textExplanation<-paste(textExplanation,
                           # "<b style='font-size:20px;'>Wild Type (WT) vs Patient rearranged locus comparison </b>",
                           # "<br><br>", 
                           targetGene,
                           " presents an expression of ",
                           round2(x = targetExpression,
                                  digits = 1),
                           " fpkm at the ", targetPhase,
                           " stage. ",
                           ##"<br><br>",
                           sep="")
    
    #########################################################################################
    ##Adding Haploinsufficiency score info (if there is)
    ##If there is not we will talk about the gene expression pattern if tissue specific
    ##########################################################################################
    
    if(targetHI_score>-1){
      
      ##So an HI score provided
      textExplanation<-paste(textExplanation,
                             "The Dosage Sensitivity Score (DS) obtained for ",
                             targetGene,
                             " is : ",
                             round2(targetHI_score, digits=2),
                             ". This score must be interpreted on a 0-1 scale (1 means that for the candidate gene, producing more or less protein than normal is likely pathogenic, 0 means that deviations from the normal amount of generated protein are not detrimental).",
                             # "<br><br>",
                             sep="") 
    }
    
    
    ###############################
    ##Polycomb Score Information
    ###############################
    textExplanation<-paste(textExplanation,
                           " The Polycomb Score obtained for ",
                           targetGene,
                           " is: ",
                           round2(polycScore,digits = 2),
                           ". This score must be interpreted on a 0-1 scale (1 means strong evidence for being a H3K27me3-polycomb based regulatory gene, 0 implies the opposite).",
                           "<br><br>",
                           sep="") 
    
    
    textExplanation<-paste(textExplanation,
                           "To be a strong candidate for 'Gain of Function', we consider important being a polycomb gene (genes covered by broad domains of H3K27me3 when inactive and with tissue specific expression patterns), since these genes tend to regulate cellular identity and altering their expression can lead to strong phenotypic alterations. ",
                           "It is also relevant to know the dosage sensitivity for the target gene. Because, having more or less protein than normal can be pathogenic for some but not all genes, Haploinsufficiency metrics are currently used as proxy for this aim. Check our manuscript to know more about how these metrics have been obtained.",
                           "<br><br>",
                           sep="")
    
    
  }else if((pathoMechanism=="LOF") && (targetPhase=="phaseFree")){
    ####EXPLANATION FOR LOF AND PHASE FREE
    ##So gene truncated and deleted, even though without considering its expression
    ##If it is associated with the phenotype, whenever it is expressed it will have an impact
    
    ##Be sure the gene is on the TAD map we are assessing, so screen them
    
    for(nPhase in 1:length(patientResults$resultsPerPhase_secondaryInfo[[targetPhase]])){
      
      evaluationMatrix<-patientResults$resultsPerPhase_secondaryInfo[[targetPhase]][[nPhase]]$matrixesGenesEvaluation
      
      if(targetGene %in% rownames(evaluationMatrix)){
        ##So the gene appears in the currently used TAD map
        ##Hence we have the info we need from it
        
        targetHI_score<-evaluationMatrix[targetGene,
                                         c("nature_HI_score",
                                           "huang_HI_score","clinGene_HI_score")] ##-1 when no score assigned
        
        targetHI_score<-max(targetHI_score)##if it is -1, it is because no HI score assigned
        
        polycScore<-evaluationMatrix[targetGene,"polyComb_score"]
        
        break
      }
    }
    
    
    ##Adding Haploinsufficiency score info (if there is)
    ##If there is not we will talk about the gene expression pattern if tissue specific
    if(targetHI_score>-1){
      
      ##So an HI score provided
      textExplanation<-paste(textExplanation,
                             "The Dosage Sensitivity Score (DS) obtained for ",
                             targetGene,
                             " is : ",
                             round2(targetHI_score, digits=2),
                             ". This score must be interpreted on a 0-1 scale (1 means that for the candidate gene, producing more or less protein than normal is likely pathogenic, 0 means that deviations from the normal amount of generated protein are not detrimental).",
                             # "<br><br>",
                             sep="") 
    }
    
    ###############################
    ##Polycomb Score Information
    ###############################
    textExplanation<-paste(textExplanation,
                           " The Polycomb Score obtained for ",
                           targetGene,
                           " is: ",
                           round2(polycScore,digits = 2),
                           ". This score must be interpreted on a 0-1 scale (1 means strong evidence for being a H3K27me3-polycomb based regulatory gene, 0 implies the opposite).",
                           "<br><br>",
                           sep="") 
    
    textExplanation<-paste(textExplanation,
                           "To be a strong candidate for 'Loss of Function',",
                           "it is relevant to know the dosage sensitivity for the target gene, Haploinsufficiency metrics are used as proxy for this aim.",
                           "Furthermore, we also find relevant to study whether the candidate gene is a polycomb one (genes covered by broad domains of H3K27me3 when inactive and with tissue specific expression patterns), since these genes tend to regulate cellular identity and altering their expression can lead to strong phenotypic alterations. ",
                           "<br><br>",
                           sep="")
    
    
    
  }else if((pathoMechanism=="GOF") && (targetPhase=="phaseFree")){
    
    ##Be sure the gene is on the TAD map we are assessing, so screen them
    
    for(nPhase in 1:length(patientResults$resultsPerPhase_secondaryInfo[[targetPhase]])){
      
      evaluationMatrix<-patientResults$resultsPerPhase_secondaryInfo[[targetPhase]][[nPhase]]$matrixesGenesEvaluation
      
      if(targetGene %in% rownames(evaluationMatrix)){
        ##So the gene appears in the currently used TAD map
        ##Hence we have the info we need from it
        
        targetHI_score<-evaluationMatrix[targetGene,
                                         c("nature_HI_score",
                                           "huang_HI_score","clinGene_HI_score")] ##-1 when no score assigned
        
        targetHI_score<-max(targetHI_score)##if it is -1, it is because no HI score assigned
        
        polycScore<-evaluationMatrix[targetGene,"polyComb_score"]
        
        break
      }
    }
    
    #########################################################################################
    ##Adding Haploinsufficiency score info (if there is)
    ##If there is not we will talk about the gene expression pattern if tissue specific
    ##########################################################################################
    
    if(targetHI_score>-1){
      
      ##So an HI score provided
      textExplanation<-paste(textExplanation,
                             "The Dosage Sensitivity Score (DS) obtained for ",
                             targetGene,
                             " is : ",
                             round2(targetHI_score, digits=2),
                             ". This score must be interpreted on a 0-1 scale (1 means that for the candidate gene, producing more or less protein than normal is likely pathogenic, 0 means that deviations from the normal amount of generated protein are not detrimental).",
                             # "<br><br>",
                             sep="") 
    }
    
    ###############################
    ##Polycomb Score Information
    ###############################
    textExplanation<-paste(textExplanation,
                           " The Polycomb Score obtained for ",
                           targetGene,
                           " is: ",
                           round2(polycScore,digits = 2),
                           ". This score must be interpreted on a 0-1 scale (1 means strong evidence for being a polycomb regulatory gene, 0 implies the opposite)",
                           "<br><br>",
                           sep="") 
    
    
    textExplanation<-paste(textExplanation,
                           "To be a strong candidate for 'Gain of Function',",
                           "it is relevant to know the dosage sensitivity for the target gene, Haploinsufficiency metrics are used as proxy for this aim.",
                           "Furthermore, we also find relevant to study whether the candidate gene is a polycomb one (genes covered by broad domains of H3K27me3 when inactive and with tissue specific expression patterns), since these genes tend to regulate cellular identity and altering their expression can lead to strong phenotypic alterations. ",
                           "<br><br>",
                           sep="")
    
  }
  
  ##Wrapping text explanation on the html
  
  geneSuscept<-paste(geneSuscept,
                     "<div class='explanationGeneSuscept'>",
                     paste("<p style='text-align:justify;'>",
                           textExplanation,
                           "</p>",
                           sep=""),
                     "</div>",
                     sep="")
  
  ######################
  ##CLOSING Section
  geneSuscept<-paste(geneSuscept,
                     "</div>",
                     sep="")
  
  return(geneSuscept)
  
}