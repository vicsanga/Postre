#######################################
## PLOTS for Gene Susceptibility
#######################################

##required data, thresholds expression
source(file = "scripts_To_Load_Data/ExpressionThresholds_ForGofAndLof.R",
       local = TRUE)

source("functions/expressionLevelPlot.R", local = TRUE)
source("functions/haploinsufficiencyScorePlot.R", local = TRUE)
source("functions/polycombScorePlot.R", local = TRUE)

##Always gene-expression barplot
plots_geneSusceptibility<-function(patientResults){
  ##We are going to store the plots in the same folder than the graphical abstracts www/graphicalSummaries/
  
  ##HI plot is done just Once per gene, since it is the same, so get all relevant genes upon parsing and then parse again doing the HI score plot
  alreadyHI_plotted<-character() ##To track if plot already done
  
  ##Same applies for polycomb score
  already_PolyCscore_plotted<-character()
  
  for(targetElement in patientResults$genesConditions_ToReport){
    ##targetElement<-"TFAP2A_LOF_NeuralCrestLate"
    
    ###I need to know whether the considered gene is selected by means of a 
    ## LOF && GOF mechanism || LONG-RANGE--DIRECT EFFECT
    
    targetGene<-unlist(strsplit(x = targetElement, split = "_", fixed = TRUE))[1]
    targetMech<-unlist(strsplit(x = targetElement, split = "_", fixed = TRUE))[2]
    targetPhase<-unlist(strsplit(x = targetElement, split = "_", fixed = TRUE))[3]
    
    targetRow_resultsPhase<-paste(targetGene,"--",targetMech, sep="")##In case interest on parsing resultsPhase
    
    #This way of approaching the GOF-LOF info, deprecated
    # pathoMechanism<-patientResults$resultsPerPhase[[targetPhase]][[targetPhase]][targetGene,]$Mechanism
    # pathoMechanism<-unlist(strsplit(x = pathoMechanism,
    #                                 split = ":",
    #                                 fixed = TRUE))[2]##GOF or LOF
    
    pathoMechanism<-targetMech
    
    ##To discriminate in dupl deletions, do something like, pathomechanism_Impact by long Range, or by direct effect
    ##As a gene duplication can trigger the disease due to neoTad and hence, by a long-range mechanism
    pathomechanism_ImpactOverGene<-patientResults$resultsPerPhase[[targetPhase]][[targetPhase]][targetRow_resultsPhase,]$GeneImpact ##GOF or LOF
    
    ###################################
    ###DOING expressionLevelPlot, this will be done twice in the case that same gene good candidate for GOF an LOF
    ##As a matter of  simplicity redoing plot for GOF and LOF if happens in same phase, 
    ##If affects performance notoriously, avoid redoing plot for gof and lof, and doing just once
    
    if(targetPhase!="phaseFree"){
      ##So we have an expression value
      targetExpression<-patientResults$resultsPerPhase_secondaryInfo[[targetPhase]][[1]]$matrixesGenesEvaluation[targetGene,
                                                                                                                 paste("FPKM_",
                                                                                                                       targetPhase,
                                                                                                                       sep="")]
      outpPath<-"www/graphicalSummaries/"
      ##outpPath<-"graphicalSummaries/"
      fullOutpPath<-paste0(outpPath,targetGene,"_",targetMech,"_",targetPhase,"_", patientResults$job_UniCode,"_expressionLevelPlot.png")
      
      ####################################
      ##png with maximum resolution 300dpi
      # png(filename = fullOutpPath, width = 12, height = 8, units = "in", res = 300 )
      png(filename = fullOutpPath, width = 480, height = 480, res = 300, pointsize = 2 )
      ##Jugar con dimensiones, y con tamanyo de labels que asi como esta no se van a ver en el html
      
      
      ##Adjust image, to exclude spaces outside canvas (drawing area)
      par(mar = c(0,0,0,0))
      
      expressionLevelPlot(targetExpression = targetExpression,
                          targetGene = targetGene,
                          targetPhase = targetPhase,
                          maxExpression = threshold_MaxExpresion,
                          minExpression = threshold_MinExpresion)
      
      dev.off()##Saving Graph
    }
    
    ########################################################################
    ###DOING Haploinsufficiency Plot (SecondaryPlot Gene Susceptibility)
    ## Plot for Dosage sensitivity study
    ########################################################################
    
    if(!(targetGene %in% alreadyHI_plotted)){
      ##Doing HI score plot
      ##If a gene do not have HI score, think about how to do plot
      ##Be sure the gene is on the TAD map we are assessing, so screen them
      
      for(nPhase in 1:length(patientResults$resultsPerPhase_secondaryInfo[[targetPhase]])){
        
        evaluationMatrix<-patientResults$resultsPerPhase_secondaryInfo[[targetPhase]][[nPhase]]$matrixesGenesEvaluation
        
        if(targetGene %in% rownames(evaluationMatrix)){
          ##So the gene appears in the currently used TAD map
          ##Hence we have the info we need from it
          targetHI_score<-evaluationMatrix[targetGene,
                                           c("nature_HI_score",
                                             "huang_HI_score","clinGene_HI_score")] ##-1 when no score assigned
          break
        }
      }
      
      targetHI_score<-max(targetHI_score)##if it is -1, it is because no HI score assigned
      
      ##Plotting Haploinsufficiency score info (if there is)
      if(targetHI_score>-1){ ##So as one hi score !=-1 there is one assigned
        
        # ##Plotting HI score
        # source("functions/haploinsufficiencyScorePlot.R", local = TRUE)
        
        outpPath<-"www/graphicalSummaries/"
        ##outpPath<-"graphicalSummaries/"
        fullOutpPath<-paste0(outpPath,targetGene,"_", patientResults$job_UniCode, "_HI_plot.png")
        
        ####################################
        ##png with maximum resolution 300dpi
        # png(filename = fullOutpPath, width = 12, height = 8, units = "in", res = 300 )
        png(filename = fullOutpPath, width = 480, height = 480, res = 300, pointsize = 2 )
        ##Jugar con dimensiones, y con tamanyo de labels que asi como esta no se van a ver en el html
        
        
        ##Adjust image, to exclude spaces outside canvas (drawing area)
        par(mar = c(0,0,0,0))
        
        haploinsufficiencyScorePlot(targetGene = targetGene,
                            targetHI_score = targetHI_score)
        
        dev.off()##Saving Graph
        
        
      }
      
      #############################
      ##add to already HI plotted
      alreadyHI_plotted<-c(alreadyHI_plotted, targetGene)
      
    }
    
    
    ##############################################################################
    ##       DOING POLYCOMB SCORE Plot (SecondaryPlot Gene Susceptibility)
    ##############################################################################
    
    if(!(targetGene %in% already_PolyCscore_plotted)){
      ##Doing Polycomb score plot
      ##If a gene do not have score, think about how to do plot
      ##Be sure the gene is on the TAD map we are assessing, so screen them
      
      for(nPhase in 1:length(patientResults$resultsPerPhase_secondaryInfo[[targetPhase]])){
        
        evaluationMatrix<-patientResults$resultsPerPhase_secondaryInfo[[targetPhase]][[nPhase]]$matrixesGenesEvaluation
        
        if(targetGene %in% rownames(evaluationMatrix)){
          ##So the gene appears in the currently used TAD map
          ##Hence we have the info we need from it
          polyCscore<-evaluationMatrix[targetGene,"polyComb_score"]##-1 when no score assigned
          break
        }
      }
      
      ##if it is -1, it is because no HI score assigned
      ##Plotting Polycomb score info (if there is)
      if(polyCscore>-1){ ##So as one polyCscore !=-1 there is one assigned
        
        # ##Plotting polyCscore
        # source("functions/polycombScorePlot.R", local = TRUE)
        
        outpPath<-"www/graphicalSummaries/"
        ##outpPath<-"graphicalSummaries/"
        fullOutpPath<-paste0(outpPath,targetGene,"_", patientResults$job_UniCode, "_polyCscore_plot.png")
        
        ####################################
        ##png with maximum resolution 300dpi
        # png(filename = fullOutpPath, width = 12, height = 8, units = "in", res = 300 )
        png(filename = fullOutpPath, width = 480, height = 480, res = 300, pointsize = 2 )
        ##Jugar con dimensiones, y con tamanyo de labels que asi como esta no se van a ver en el html
        
        
        ##Adjust image, to exclude spaces outside canvas (drawing area)
        par(mar = c(0,0,0,0))
        
        polycombScorePlot(targetGene = targetGene,
                       polyCscore = polyCscore)
        
        dev.off()##Saving Graph
        
        
      }
      
      #############################
      ##add to already HI plotted
      already_PolyCscore_plotted<-c(already_PolyCscore_plotted, targetGene)
      
    }
  }
}




