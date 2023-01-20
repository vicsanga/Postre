#####################################################
## Creating Master Summary Matrix
## Which integrates the main results from each phase 
#####################################################

master_Summary_resultsMatrix<-function(patientResults){

  phases<-names(patientResults$resultsPerPhase)
  ##Integrated results matrix
  genesInvolved<-rownames(patientResults$allAffectedGenes_positionalInfo)
  
  ##Before Continuing if nGenesInvolved == 0, so no genes involved for any context for the phenotype
  ##Skip computations and notify it, since an HTML will be generated explaining this
  
  if(length(genesInvolved)==0){
    
    return("NO GENES ASSOCIATED WITH SV")
    
  }else{
    ##So there are genes involved for at least a context, hence it makes sense to generate the heatmap and provide information
    Matrix_integratedResults<-as.data.frame(matrix(data = NA, nrow = length(genesInvolved), ncol = 1 +  length(phases)))
    rownames(Matrix_integratedResults)<-genesInvolved
    colnames(Matrix_integratedResults)<-c("affected_gene",phases)
    ##affected genes column
    Matrix_integratedResults$affected_gene<-rownames(Matrix_integratedResults)
    
    ##Add max score column, to retrieve the maxScore for the gene along the phases
    ##We will use this value to rank the matrix in decreasing order
    ##Once the matrix is sorted, we will remove this column
    Matrix_integratedResults$maxScore<-NA
    
    
    ##In the phase columns we are going to store some metadata:
    ## Phase Avg Score; Type of Mechanism; GeneImpact , Even if the gene is duplicated, so there is an impact on it, the effect can be due to a long range mechanism
    ##We record it for both, LOF and GOF score. Since if both relevant, both highlighted
    ## {PhaseAvgScore};{LOF|GOF|MIX};{long-range//direct}
    ## eg 0.95;LOF;long-range__0.5;GOF;long-range
    
    ##Gene Impact Column
    ##Direct|LongRange
    
    Matrix_integratedResults$GeneImpact<-NA
    
    
    ######Let's fill the dataframe
    for (gene in genesInvolved) {
      
      scoresTracking<-numeric() ##To track avg scores along phases for a gene. We will get all, both GOF and LOF
      trackingGenesImpact<-character()##will be always the same except for dup that can be longRange or Direct impact
      
      for(targetPhase in phases){
        ##Track average score, to get max, so that afterwords we rank the matrix by this value
        resultsTargetPhase<-patientResults$resultsPerPhase[[targetPhase]][[targetPhase]]
        
        ##As with phase the TAD map changes, if a gene does not appear, do not try to get info (will rise error)
        if(gene %in% resultsTargetPhase$affected_gene){
          ##So the gene is involved in the current phase
          ##We need to get info for both, GOF and LOF scores
          
          cellData<-character()##To store gene info
          cellDataTracking<-list() ##to store LOF and GOF
          for(targetScore in c("LOF","GOF")){
            
            targetGeneScoreName<-paste(gene, "--", targetScore, sep="")
            
            targetPhase_results<-resultsTargetPhase[targetGeneScoreName,]
            
            #Get score, and introduce it in the metadata info
            currentScore<-round2(targetPhase_results$Average_Score,digits = 2)
            scoresTracking<-c(scoresTracking, currentScore)
            
            trackingGenesImpact<-c(trackingGenesImpact, targetPhase_results$GeneImpact)
            
            ## Getting Type of Mechanism
            mech_split<-unlist(strsplit(x = targetPhase_results$Mechanism,
                                        split = ":"))
            
            if(mech_split[1]=="MIX"){
              ##Hence in some TAD map, GOF, in other LOF... 
              ##Hence we summarize that with mix
              stop("ERROR: THIS SHOULD NOT BE REACHABLE, REGARDING LAST IMPLEMENTATION IN NOV 2021")
              currentMechanism<-mech_split[1]
            }else{
              ##Hence, either ALL:LOF or GOF mechanism in all TAD maps
              currentMechanism<-mech_split[2]
            }
            ##Get geneImpact
            geneImpact<-targetPhase_results$GeneImpact
            
            cellData<-paste(currentScore,
                            ";",
                            currentMechanism,
                            ";",
                            geneImpact,
                            sep = "",
                            collapse = "")
            
            cellDataTracking[[targetScore]]<-cellData
            
          }
          ##Maybe removing pathomech, when dup, at some point, but think how to handle it to avoid downstream problems
          Matrix_integratedResults[gene,targetPhase]<-paste(cellDataTracking$LOF, "--", cellDataTracking$GOF, sep = "")
          
        }
        
      }
      
      #############################################
      ##Adding Overall geneImpact to Last Column 
      #############################################
      
      trackingGenesImpact<-unique(trackingGenesImpact)
      
      if((length(trackingGenesImpact) == 1)){
        ##So for all the phases, same gene impact.
        Matrix_integratedResults[gene,"GeneImpact"]<-trackingGenesImpact
        
      }else if((length(trackingGenesImpact) > 1) && (patientResults$patientInfo$TypeSV == "Duplication")){
        ##So more than one definition of gene impact
        ##it is expected to occur for gene duplications if we predict sth by LongRange eg "LongRange_geneDuplication" or sth direct "Direct_geneDuplication"
        ##For the overall geneImpact summary column we put the broader concept. But for the stage specific analysis, we work with the more specific condition, for the reports
        
        Matrix_integratedResults[gene,"GeneImpact"]<-"Direct_LongRange_geneDuplication"
        
      }else{
        ##multiple impacts on none duplications... Not expected
        stop("Multiple Gene impacts in a not Duplication. Not expected. Check what is wrong.")
      }
      
      ###Once all phases have been screened
      ##Now, add max score, to sort them afterwards per maxScore. Whatever gene at least will appear in one phase, so, at least, there is one couple of scores (LOF-GOF) to get the max
      Matrix_integratedResults[gene,"maxScore"]<-max(scoresTracking)
      
    }
    
    ##Ordering per MaxScore column
    Matrix_integratedResults<-Matrix_integratedResults[order(Matrix_integratedResults$maxScore, decreasing = TRUE),]
    
    ##MaxScores smaller than 0,  set to 0
    Matrix_integratedResults$maxScore[Matrix_integratedResults$maxScore<0]<-0
    
    # #Once sorted per MaxScore, remove MaxScore column
    # Matrix_integratedResults$maxScore<-NULL
    # I need it to filter for genes to do the report
    return(Matrix_integratedResults) 
  }
  
}




