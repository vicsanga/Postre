###################################################################
## Function to integrate predictions after Single TAD predictions
###################################################################
integratingTAD_predictions<-function(resultsPerTADmap, targetPhase){
  
  #####Retrieving list of involved genes, all of them will appear in the Final Result
  genesInvolved<-character()
  #getting also the landing of SV tendency, if is IntraTAD in all or InterTAD or Uncertain(if for some intraTAD and for others inter)
  SV_landing_tendency<-character()
  
  for(tadMap in names(resultsPerTADmap)){
    genesAffectedInMap<-rownames(resultsPerTADmap[[tadMap]]$matrixesGenesEvaluation)
    genesInvolved<-unique(c(genesInvolved,genesAffectedInMap))
    
    SV_landing_tendency<-c(SV_landing_tendency, resultsPerTADmap[[tadMap]]$affectedRegions$SV_landing)  
    
  }
  
  ##Object to return results
  list_integratedResults<-list()
  
  ##IF AT THIS POINT, NO GENES ASSOCIATED WITH THE SV SKIP THE REST OF THE SCRIPT
  ##LEAVE EMPTY RESULTS TABLE
  if(length(genesInvolved)==0){
    ##So, no genes affected
    ##Recording results per phase  
    list_integratedResults[[targetPhase]]<-NULL
  }else{
    ##So, at least one gene affected
    ##Averaged//Integrated results matrix (we need a list, to store one per phase)
    IntegratedResultsMatrixColumns<-c("Mechanism","Average_Score","GeneImpact",
                                      "N_TAD_maps_AppearInvolved",
                                      "minScore","maxScore",
                                      "avg_nEnh_InitialLeft",
                                      "avg_nEnh_InitialRight",
                                      "avg_nEnh_InTheOtherDomain",
                                      ##Rearranged scenario
                                      "avg_nEnh_KeptLeft",
                                      "avg_nEnh_KeptRight",
                                      "avg_nEnh_GainedFromTheOtherDomain",
                                      
                                      "gene_Breakpoint",
                                      
                                      ###We want to track also enh acetilation levels
                                      ##WT situation
                                      "avg_AcetilationEnh_InitialLeft",
                                      "avg_AcetilationEnh_InitialRight",
                                      "avg_AcetilationEnh_InTheOtherDomain",
                                      ##Rearranged scenario
                                      "avg_AcetilationEnh_KeptLeft",
                                      "avg_AcetilationEnh_KeptRight",
                                      "avg_AcetilationEnh_GainedFromTheOtherDomain",
                                      
                                      #Tracking whether SV is intraTAD, interTAD, Unclear (whether for some map is IntraTAD, for some other InterTAD,can happen if TAD smaller)
                                      #SV_landing: "intraTAD", "interTAD", "Unclear"
                                      "SV_landing"
    )
    Matrix_integratedResults<-as.data.frame(matrix(data = NA, nrow = length(genesInvolved)*2, ##*2 because we track GOF and LOF separately
                                                   ncol = length(IntegratedResultsMatrixColumns) ))
    colnames(Matrix_integratedResults)<-IntegratedResultsMatrixColumns
    
    ##I want to track the results independently for GOF and for LOF
    ##To append LOF and GOF adding "--"LOF and --GOF since no gene name presents "--" but some do present "-" or "_"
    ##So for posterior splits, if we do split(geneName, by="-") if the gene already presents "-" the result will be a mess 
    targetRowNames<-sort(c(
      paste(genesInvolved,"--LOF", sep=""),
      paste(genesInvolved,"--GOF", sep="")))
    
    rownames(Matrix_integratedResults)<-targetRowNames
    
    ##MaxScore por si se lleva en un mapa un score de 1 trackearlo
    ##Mech can be ALL: GOF, ALL: LOF, MIX: GOF-LOF
    ##Impact over gene, it can be, long-range, direct-geneTruncated...
    
    
    ###Let's fill the matrix
    
    ##Matrix to hold results per phase
    targetMatrix<-Matrix_integratedResults
    
    ##Let's obtain gene info along TAD maps
    for(geneCond in targetRowNames){
      ##geneCond<-"TFAP2A_LOF"
      gene<-unlist(strsplit(x = geneCond, split="--", fixed = TRUE))[1]
      targetMech<-unlist(strsplit(x = geneCond, split="--", fixed = TRUE))[2] ##Aggregating GOF or LOF information
      maxScores<-numeric()
      nTADmapsAffected<-0
      geneImpacts<-character()##longrange or directGeneTruncated ...
      mechanisms<-character()##GOF vs LOF, or both
      
      geneBreakp<-character()##To track gene breakpoint...
      ##It can be NA if is located in the domains between directly affected TADs. 
      ##For instance a TAD entirely deleted between two partially deleted TADs.
      ##Lo llamaremos como "No_AssociatedBreakpoint"
      ##Also there can be multiple association
      ##It can occur a tad map slightly bigger so a gene not previously associated to a breakpoint (hence its breakpoint was NA) 
      ##now it is...
      ##Or lets imagine two consecutive TADs one small deletion, depending if tad maps change or no maybe we associate first a
      ##gene in the border of a tad to a breakpoin, and with another tad map we associate it to another breakpoint
      ##Hence we can have >1 "breakpoint" associated to a gene.If that occur, lets assign "notClear_BreakpointAssociation"
      ##lo llamaremos como:NotClear_BreakpointAssociation
      
      #To track WT situation
      nEnh_initial_LEFT<-numeric()##To track number enh to the left side (smaller positions) of the gene
      nEnh_initial_RIGHT<-numeric()
      nEnh_initial_OtherDOMAIN<-numeric()##Number of enhancers available in the other domain with whom the rearrangement is produced
      
      ###To track rearranged situation
      nEnh_kept_LEFT<-numeric()##To track number enh to the left side (smaller positions) of the gene
      nEnh_kept_RIGHT<-numeric()
      nEnh_gained_OtherDOMAIN<-numeric()##Number of enhancers available in the other domain with whom the rearrangement is produced
      
      ###Same for enhancer acetilation levels
      #To track WT situation
      acetilationEnh_initial_LEFT<-numeric()
      acetilationEnh_initial_RIGHT<-numeric()
      acetilationEnh_initial_OtherDOMAIN<-numeric()
      
      #To track rearranged situation
      acetilationEnh_kept_LEFT<-numeric()
      acetilationEnh_kept_RIGHT<-numeric()
      acetilationEnh_gained_OtherDOMAIN<-numeric()
      
      for(tadMap in names(resultsPerTADmap)){
        ##tadMap<-"H1-ESC_Dixon2015-raw_TADs.txt"
        resultsPhase<-resultsPerTADmap[[tadMap]][[targetPhase]]
        
        ##getting MatrixesGeneEvaluation info
        geneEvaluationMatrix_phase<-resultsPerTADmap[[tadMap]]$matrixesGenesEvaluation
        
        ##subset for the gene
        ##if it is found affected in this TADs map
        if(gene %in% rownames(resultsPhase)){
          resultsPhase<-resultsPhase[gene,]
          
          ##nTads it Appears involved
          nTADmapsAffected<-nTADmapsAffected+1
          
          geneImpacts<-unique(c(geneImpacts, resultsPhase$type))
          
          ###Working for each mechanism, per separte, GOF and LOF, so forget about Max_Types. We record information per separate
          ##So forget about Maxs, which imply merging info GOF-LOF and we decided to keep them per separate
          # mechanisms<-unique(c(mechanisms,resultsPhase$Max_type))
          mechanisms<-unique(c(mechanisms,targetMech))
          
          geneBreakp<-unique(c(geneBreakp, geneEvaluationMatrix_phase[gene,"gene_Breakpoint"]))
          
          ##It can occur a tad map slightly bigger so a gene not previously associated to a breakpoint (hence its breakpoint was NA) 
          ##now it is...
          ##Or lets imagine two consecutive TADs one small deletion, depending if tad maps change or no maybe we associate first a
          ##gene in the border of a tad to a breakpoint, and with another tad map we associate it to another breakpoint
          ##Hence we can have >1 "breakpoint" associated to a gene.If that occur, lets assign "notClear_BreakpointAssociation"
          
          
          ##not unique here, we want to do the average afterwards
          ##Again, forgetting about maxScores, we want to work with both LOF and GOF scores per separate
          # maxScores<-c(maxScores, resultsPhase$Max_Score)
          
          maxScores<-c(maxScores, resultsPhase[[paste(targetMech,"_score", sep="")]])
          
          if(targetPhase != "phaseFree"){
            ##Because for this phase we do not take into account enh info as we do not have
            ##########################################
            ## Gett matrixesGeneEvaluation for plots
            #nEnhInitialLeft id column
            idCol<-paste0("nEnh_ToTheLeft_initial_", targetPhase)
            nEnh_initial_LEFT<-c(nEnh_initial_LEFT, 
                                 geneEvaluationMatrix_phase[gene,idCol])
            
            #nEnhInitialRight id column
            idCol<-paste0("nEnh_ToTheRight_initial_", targetPhase)
            nEnh_initial_RIGHT<-c(nEnh_initial_RIGHT, 
                                  geneEvaluationMatrix_phase[gene,idCol])
            
            #nEnhInTheOtherDomain
            idCol<-paste0("nEnhancers_maxAvailableInTheOtherDomain_", targetPhase)
            nEnh_initial_OtherDOMAIN<-c(nEnh_initial_OtherDOMAIN, 
                                        geneEvaluationMatrix_phase[gene,idCol])
            
            #nEnhKeptLeft id column
            idCol<-paste0("nEnh_ToTheLeft_kept_", targetPhase)
            nEnh_kept_LEFT<-c(nEnh_kept_LEFT, 
                              geneEvaluationMatrix_phase[gene,idCol])
            
            #nEnhkeptRight id column
            idCol<-paste0("nEnh_ToTheRight_kept_", targetPhase)
            nEnh_kept_RIGHT<-c(nEnh_kept_RIGHT, 
                               geneEvaluationMatrix_phase[gene,idCol])
            
            #nEnhInTheOtherDomain
            idCol<-paste0("nEnhancers_gained_", targetPhase)
            nEnh_gained_OtherDOMAIN<-c(nEnh_gained_OtherDOMAIN, 
                                       geneEvaluationMatrix_phase[gene,idCol])
            
            ###########################################
            ## For enhancer acetilation levels
            ###########################################
            ###Info enhancer Acetilation levels
            idCol<-paste0("enhancers_acetilation_ToTheLeft_initial_", targetPhase)
            acetilationEnh_initial_LEFT<-c(acetilationEnh_initial_LEFT, 
                                           geneEvaluationMatrix_phase[gene,idCol])
            
            #acetilation EnhInitialRight id column
            idCol<-paste0("enhancers_acetilation_ToTheRight_initial_", targetPhase)
            acetilationEnh_initial_RIGHT<-c(acetilationEnh_initial_RIGHT, 
                                            geneEvaluationMatrix_phase[gene,idCol])
            
            #acetilation EnhInTheOtherDomain
            idCol<-paste0("enhancers_maxAcetilationAvailableInTheOtherDomain_", targetPhase)
            acetilationEnh_initial_OtherDOMAIN<-c(acetilationEnh_initial_OtherDOMAIN, 
                                                  geneEvaluationMatrix_phase[gene,idCol])
            
            #acetilation EnhKeptLeft id column
            idCol<-paste0("enhancers_acetilation_ToTheLeft_kept_", targetPhase)
            acetilationEnh_kept_LEFT<-c(acetilationEnh_kept_LEFT, 
                                        geneEvaluationMatrix_phase[gene,idCol])
            
            #acetilation EnhkeptRight id column
            idCol<-paste0("enhancers_acetilation_ToTheRight_kept_", targetPhase)
            acetilationEnh_kept_RIGHT<-c(acetilationEnh_kept_RIGHT, 
                                         geneEvaluationMatrix_phase[gene,idCol])
            
            #acetilation EnhInTheOtherDomain
            idCol<-paste0("enhancers_acetilation_gained_", targetPhase)
            acetilationEnh_gained_OtherDOMAIN<-c(acetilationEnh_gained_OtherDOMAIN, 
                                                 geneEvaluationMatrix_phase[gene,idCol])
          }
        }
        
        
      }
      
      ######################################################################################
      ##Type of SV (inter TAD) or intra TAD is independent on the gene, only depend on TADs
      ######################################################################################
      ##resultsPerTADmap 
      
      
      ###UP to here, we have the info to fill the matrix for this gene for this phase
      ##Average, over the total number of TAD maps considered
      
      averageMaxScores<-sum(maxScores)/length(names(resultsPerTADmap))
      
      ##WT situation
      average_nEnh_Initial_LEFT<-round2(sum(nEnh_initial_LEFT)/length(names(resultsPerTADmap))
                                        ,digits = 0)
      average_nEnh_Initial_RIGHT<-round2(sum(nEnh_initial_RIGHT)/length(names(resultsPerTADmap)),
                                         digits = 0)
      average_nEnh_Initial_OtherDOMAIN<-round2(sum(nEnh_initial_OtherDOMAIN)/length(names(resultsPerTADmap)),
                                               digits = 0)
      
      #Rearranged scenario
      average_nEnh_Kept_LEFT<-round2(sum(nEnh_kept_LEFT)/length(names(resultsPerTADmap))
                                     ,digits = 0)
      average_nEnh_Kept_RIGHT<-round2(sum(nEnh_kept_RIGHT)/length(names(resultsPerTADmap)),
                                      digits = 0)
      average_nEnh_Gained_OtherDOMAIN<-round2(sum(nEnh_gained_OtherDOMAIN)/length(names(resultsPerTADmap)),
                                              digits = 0)
      
      ##########################################
      ##For enhancers acetilation levels
      ##########################################
      ##WT situation
      average_acetilationEnh_Initial_LEFT<-round2(sum(acetilationEnh_initial_LEFT)/length(names(resultsPerTADmap))
                                                  ,digits = 0)
      average_acetilationEnh_Initial_RIGHT<-round2(sum(acetilationEnh_initial_RIGHT)/length(names(resultsPerTADmap)),
                                                   digits = 0)
      average_acetilationEnh_Initial_OtherDOMAIN<-round2(sum(acetilationEnh_initial_OtherDOMAIN)/length(names(resultsPerTADmap)),
                                                         digits = 0)
      
      #Rearranged scenario
      average_acetilationEnh_Kept_LEFT<-round2(sum(acetilationEnh_kept_LEFT)/length(names(resultsPerTADmap))
                                               ,digits = 0)
      average_acetilationEnh_Kept_RIGHT<-round2(sum(acetilationEnh_kept_RIGHT)/length(names(resultsPerTADmap)),
                                                digits = 0)
      average_acetilationEnh_Gained_OtherDOMAIN<-round2(sum(acetilationEnh_gained_OtherDOMAIN)/length(names(resultsPerTADmap)),
                                                        digits = 0)
      
      #######################
      ##Filling the matrix
      #######################
      ##Here not gene, but target row
      targetMatrix[geneCond,"Average_Score"]<-averageMaxScores
      targetMatrix[geneCond,"minScore"]<-min(maxScores)
      targetMatrix[geneCond,"maxScore"]<-max(maxScores)
      
      ##Average WT condition
      targetMatrix[geneCond,"avg_nEnh_InitialLeft"]<-average_nEnh_Initial_LEFT
      targetMatrix[geneCond,"avg_nEnh_InitialRight"]<-average_nEnh_Initial_RIGHT
      targetMatrix[geneCond,"avg_nEnh_InTheOtherDomain"]<-average_nEnh_Initial_OtherDOMAIN
      
      ##Average rearranged scenario
      targetMatrix[geneCond,"avg_nEnh_KeptLeft"]<-average_nEnh_Kept_LEFT
      targetMatrix[geneCond,"avg_nEnh_KeptRight"]<-average_nEnh_Kept_RIGHT
      targetMatrix[geneCond,"avg_nEnh_GainedFromTheOtherDomain"]<-average_nEnh_Gained_OtherDOMAIN
      
      ###############################
      ##For enh acetilation levels
      ##Average WT condition
      targetMatrix[geneCond,"avg_AcetilationEnh_InitialLeft"]<-average_acetilationEnh_Initial_LEFT
      targetMatrix[geneCond,"avg_AcetilationEnh_InitialRight"]<-average_acetilationEnh_Initial_RIGHT
      targetMatrix[geneCond,"avg_AcetilationEnh_InTheOtherDomain"]<-average_acetilationEnh_Initial_OtherDOMAIN
      
      ##Average rearranged scenario
      targetMatrix[geneCond,"avg_AcetilationEnh_KeptLeft"]<-average_acetilationEnh_Kept_LEFT
      targetMatrix[geneCond,"avg_AcetilationEnh_KeptRight"]<-average_acetilationEnh_Kept_RIGHT
      targetMatrix[geneCond,"avg_AcetilationEnh_GainedFromTheOtherDomain"]<-average_acetilationEnh_Gained_OtherDOMAIN
      
      ######################################################
      ##Keept track of nTAD maps the gene appears involved
      targetMatrix[geneCond,"N_TAD_maps_AppearInvolved"]<-nTADmapsAffected
      
      ##Gene associated breakpoint
      ##For +info on the options check code coments above
      geneBreakp<-unique(geneBreakp)
      if(length(geneBreakp)>1){
        targetMatrix[geneCond,"gene_Breakpoint"]<-"NotClear_BreakpointAssociation"
      }else if(is.na(geneBreakp)){
        targetMatrix[geneCond,"gene_Breakpoint"]<-"No_AssociatedBreakpoint"
      }else{
        targetMatrix[geneCond,"gene_Breakpoint"]<-geneBreakp
      }
      
      ##if it is or not good for validation will be for all the same
      ##since if it is long range always good, if broken, depends on expression not on TAD map
      
      if(length(geneImpacts)==1){
        targetMatrix[geneCond,"GeneImpact"]<-geneImpacts 
        
      }else if((length(geneImpacts)==2) && ("Direct_LongRange_geneDuplication" %in% geneImpacts)){
        ##To deal with a bug
        ##gene impacts was presenting: "LongRange_geneDuplication"        "Direct_LongRange_geneDuplication"
        ##Para un gen GCNT2
        ##Si tenemos los dos esos, por temas de predicciones asociada a TADs variability
        ##Poner Direct_LongRange_geneDuplication y a funcionar (que engloba a ambos)
        targetMatrix[geneCond,"GeneImpact"]<-"Direct_LongRange_geneDuplication" 
      }
      
      if(length(mechanisms)>1 || mechanisms=="equally_likely"){
        ##I think this is not longer possible, the way we are handling the patients
        stop("ERROR: COULD SHOULD NOT ENTER HERE, BECAUSE WE TREAT PER SEPARATE GOF AND LOF. The aggregated handling is deprecated")
        ##it implies either in some TAD maps for a gene equally likely GOF or LOF
        ##or in a map more likely to be GOF and in another to be LOF
        ##due to enhancers playing around scenario
        ##I guess for this cases scores will tend to be low, but in any case, track this
        targetMatrix[geneCond, "Mechanism"]<-"MIX:GOF-LOF"
      }else{
        ##quitamos el _score
        mechanisms<-unlist(strsplit(mechanisms,split = "_", fixed = TRUE))[1]
        targetMatrix[geneCond, "Mechanism"]<-paste0("ALL:",mechanisms)
      }
      
      ##Last sort the matrix in decreasing order of average score
      targetMatrix<-targetMatrix[order(targetMatrix$Average_Score,decreasing = TRUE),]
    }  
    
    ##Adding a column (first one)
    ##With the name of the studied gene
    geneOrderMatrix<-unlist(
      lapply(X = rownames(targetMatrix), FUN = function(x){
        return(
          unlist(strsplit(x=x, split = "--", fixed = TRUE))[1]
        )
      }))
    
    targetMatrix<-cbind(affected_gene=geneOrderMatrix,targetMatrix)
    #Convert affected gene name to character, not factor (as it is introduced)
    targetMatrix$affected_gene<-as.character(targetMatrix$affected_gene)
    
    ##Add one, with info, affectedGene + mech??
    
    ##########################
    ##Adding SV_landing info
    
    if(all(SV_landing_tendency=="IntraTAD")){
      targetMatrix$SV_landing<-"IntraTAD"
      
    }else if(all(SV_landing_tendency=="InterTAD")){
      targetMatrix$SV_landing<-"InterTAD"
      
    }else{
      ##TAD dynamics, for some intraTAD for other InterTAD... unclear
      targetMatrix$SV_landing<-"Uncertain"
    }
    
    ##Recording results per phase  
    list_integratedResults[[targetPhase]]<-targetMatrix
    
  }
  
  
  ######################################
  ## Return merged per TAD map results
  ######################################
  
  return(list_integratedResults)
  
}


###Igual luego decir, bueno si buscamos y no hay high scores, mirar a ver si a nivel individual hay. Ya que el paciente esta enfermo,
###aunque bueno esto no es tan robusto como si sale en todos los mapas de TADs