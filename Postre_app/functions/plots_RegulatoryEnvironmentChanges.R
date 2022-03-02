###################################################################
## Function to make the plots regarding enhancer Landscape Changes
## It is only done if long-range and not PhaseFree conditions
## If direct impact. And no enh info... unnecessary
## Another image will be displayed in those conditions
####################################################################
source("functions/sv_barplotting.R")
source("functions/makingGraphs/regulatoryImpact_directEffects.R")
plots_regulatoryEnvironmentChanges<-function(patientResults){

  ##We are going to store the plots in the same folder than the graphical abstracts www/graphicalSummaries/
  # allGenes<-patientResults$allAffectedGenes_positionalInfo$geneSymbol
  # allPhases<-names(patientResults$resultsPerPhase)

  for(targetElement in patientResults$genesConditions_ToReport){

    targetGene<-unlist(strsplit(x = targetElement, split = "_", fixed = TRUE))[1]
    targetMech<-unlist(strsplit(x = targetElement, split = "_", fixed = TRUE))[2]
    targetPhase<-unlist(strsplit(x = targetElement, split = "_", fixed = TRUE))[3]
    
    targetRow_resultsPhase<-paste(targetGene,"--",targetMech, sep="")
    geneImpact<-patientResults$resultsPerPhase[[targetPhase]][[targetPhase]][targetRow_resultsPhase,"GeneImpact"]

    ##GOF or LOF
    #This way of approaching the GOF-LOF info, deprecated
    # pathoMechanism<-patientResults$resultsPerPhase[[targetPhase]][[targetPhase]][targetGene,]$Mechanism
    # pathoMechanism<-unlist(strsplit(x = pathoMechanism,
    #                                 split = ":",
    #                                 fixed = TRUE))[2]##GOF or LOF
    
    pathoMechanism<-targetMech
    
    ##Probably not interested on doing for phaseFree at least for now
    ##But we could make some total addup and predict how many enh is missing the gene at some point in development
    ##FilterOut PhaseFree
    
    #######################
    ## Getting numbers
    infoBarplot<-patientResults$resultsPerPhase[[targetPhase]][[targetPhase]][targetRow_resultsPhase,]
    
    ##Computing INITIAL numbers
    nEnh_Initial<-infoBarplot$avg_nEnh_InitialLeft + 
      infoBarplot$avg_nEnh_InitialRight
    
    acetilation_Initial<-infoBarplot$avg_AcetilationEnh_InitialLeft +
      infoBarplot$avg_AcetilationEnh_InitialRight
    
    ##Computing FINAL numbers
    if(pathoMechanism == "LOF"){
      ##For LOF we only consider cognate enh, since doing the balance gained-lost was losing some enhanceropathies related with SATB2 or SLC2A1
      nEnh_Final<-infoBarplot$avg_nEnh_KeptLeft + 
        infoBarplot$avg_nEnh_KeptRight 
      
      acetilation_Final<-infoBarplot$avg_AcetilationEnh_KeptLeft +
        infoBarplot$avg_AcetilationEnh_KeptRight   
      
    }else if(pathoMechanism == "GOF"){
      
      ##For GOF we consider what is kept and of course what is gained
      nEnh_Final<-infoBarplot$avg_nEnh_KeptLeft + 
        infoBarplot$avg_nEnh_KeptRight +
        infoBarplot$avg_nEnh_GainedFromTheOtherDomain
      
      acetilation_Final<-infoBarplot$avg_AcetilationEnh_KeptLeft +
        infoBarplot$avg_AcetilationEnh_KeptRight +
        infoBarplot$avg_AcetilationEnh_GainedFromTheOtherDomain  
      
    }

    
    outpPath<-"www/graphicalSummaries/"
    ##outpPath<-"graphicalSummaries/"
    fullOutpPath<-paste0(outpPath,targetGene,"_",pathoMechanism,"_",targetPhase,
                         "_", patientResults$job_UniCode,"_regulatoryLandscapeBarplotChanges.png")
    
    ##For Direct Impacts it will will be another graphical display of the SV impact
    
    ####################################
    ##png with maximum resolution 300dpi
    # png(filename = fullOutpPath, width = 12, height = 8, units = "in", res = 300 )
    
    
    
    if(((geneImpact == "LongRange") ||
      (geneImpact == "LongRange_geneDuplication") ||
      (geneImpact == "Direct_LongRange_geneDuplication")) &&
      (targetPhase != "phaseFree") ){
      png(filename = fullOutpPath, width = 12, height = 8, units = "in", res = 300 )
      ## It is only done (for now) if long-range and not PhaseFree conditions
      ## Becuase in phaseFree no enh info. And if not LongRange enh Info not considered
      ##LongRange: gene intact, 
      ##LongRange_geneDuplication means gene duplicated but pathological mechanism predicted by LongRange effects
      ##Direct_LongRange_geneDuplication, gene expressed duplicated & gaining a lot of enh difficult to predict main cause
      #################
      ##Barplotting
      #same function changing axis, and title
      
      #source("functions/sv_barplotting.R")##devolver el p1 o p2 object de ggplot2??
      
      # outpPath<-"www/graphicalSummaries/"
      # ##outpPath<-"graphicalSummaries/"
      # fullOutpPath<-paste0(outpPath,targetGene,"_",targetPhase,
      #                      "_", patientResults$job_UniCode,"_regulatoryLandscapeBarplotChanges.png")
      # 
      # ####################################
      # ##png with maximum resolution 300dpi
      # png(filename = fullOutpPath, width = 12, height = 8, units = "in", res = 300 )
      ##Jugar con dimensiones, y con tamanyo de labels que asi como esta no se van a ver en el html
      
      par(mfrow=c(1,2))
      sv_barplotting(barplotValues = c(c(nEnh_Initial, nEnh_Final)), 
                     tagsBarplotValues = c("Control","Patient"),
                     title= "Enh. Number Changes")
      
      sv_barplotting(barplotValues = c(acetilation_Initial, acetilation_Final), 
                     tagsBarplotValues = c("Control","Patient"),
                     title= "Enh. H3K27ac Changes")
      
      dev.off()##Saving Graph
      
    }else if(geneImpact == "Direct_geneDeletion"){
      png(filename = fullOutpPath, width = 6, height = 2, res = 300, units="in") ##el del barplot de HI score
      ##Adjust image, to exclude spaces outside canvas (drawing area)
      par(mar = c(0,0,0,0))
      
      plot(x=0:20, y=0:20, #type = "n",
           ylim = c(5,12),
           xlim = c(-5,40),
           xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '') ##Upon completion, REMOVE AXIS
    
      ##Ploting WT-Control, reference situation
      wt_section_RegulatoryImpact(targetGene = targetGene)
      
      #Ploting SV, rearranged part
      deletion_regulatoryImpact(targetGene = targetGene)
      
      dev.off()##Saving Graph
    }else if(geneImpact == "Direct_geneTruncation"){
      png(filename = fullOutpPath, width = 6, height = 2, res = 300, units="in") ##el del barplot de HI score
      ##Adjust image, to exclude spaces outside canvas (drawing area)
      par(mar = c(0,0,0,0))
      
      plot(x=0:20, y=0:20, #type = "n",
           ylim = c(5,12),
           xlim = c(-5,40),
           xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '') ##Upon completion, REMOVE AXIS
      
      ##Ploting WT-Control, reference situation
      wt_section_RegulatoryImpact(targetGene = targetGene)
      
      #Ploting SV, rearranged part
      truncation_regulatoryImpact(targetGene = targetGene)
      
      
      dev.off()##Saving Graph
    }else if(geneImpact == "Direct_geneDuplication"){
      ##browser()
      png(filename = fullOutpPath, width = 6, height = 2, res = 300, units="in") ##el del barplot de HI score
      ##Adjust image, to exclude spaces outside canvas (drawing area)
      par(mar = c(0,0,0,0))
      
      plot(x=0:20, y=0:20, #type = "n",
           ylim = c(5,12),
           xlim = c(-5,40),
           xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '') ##Upon completion, REMOVE AXIS
      
      ##Ploting WT-Control, reference situation
      wt_section_RegulatoryImpact(targetGene = targetGene)
      
      #Ploting SV, rearranged part
      duplication_regulatoryImpact(targetGene = targetGene)
      
      
      dev.off()##Saving Graph
    }
    # dev.off()##Saving Graph
  }
}