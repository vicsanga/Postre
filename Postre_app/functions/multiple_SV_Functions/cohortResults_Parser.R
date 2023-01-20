##########################
##Parsing Cohort Results
##########################
##minscore, minimum required score to highlight as pathogenic

cohortResults_Parser<-function(minScore,all_patientResults, consideredPheno, discardRelevantByBrokenGene, AllPatientsInfo){

  ##Load required auxiliar function
  source(file = "functions/multiple_SV_Functions/addingInfoToMatrix.R",
         local = TRUE)
  
  ##Load list, associated phases per Pheno
  source(file = "functions/GenomicData_Loader.R", local = TRUE)
  
  ##To store results
  perPheno_cohortResults<-list()
  
  ###To track for which patient//SV we get an error. In general, do not matter at which point
  sv_withErrorRised<-character()
  #To track in same way no gene found
  sv_with_NO_geneFound<-character()
  
  ##To track number of candidate genes per SV and SV pathogenic score
  ##We can also estimate the %of predictedSV in general afterwards from this object
  # candidateGenesColumns<-c("SVid","Phenotype","PathogenicScore","Ncandidates","Ncausative","Candidates","Causative","TypeSV","N_LR_Mech","N_Direct_Mech")
  candidateGenesColumns<-c("SVid","Phenotype","PathogenicScore","Pathogenic","Ncausative","Causative","Ncandidates","Candidates","TypeSV","N_LR_Mech","N_Direct_Mech")
  ##"N_LR_Mech","N_Direct_Mech" to store the number of times that LR or Direct predicted, eg for 2 genes Direct and 1 LR
  ##Used for internal analyses
  candidateGenesInfo<-as.data.frame(matrix(data = NA, nrow = 0, ncol = length(candidateGenesColumns)))
  colnames(candidateGenesInfo)<-candidateGenesColumns

  
  ##Aggregate info
  geneRecurrencyDfColumns<-c("Gene","Phenotype","N SVs", "N Long-Range", "N Coding","SVid","Long-Range","Coding")
  geneRecurrencyDf<-as.data.frame(matrix(data = NA, nrow = 0, ncol = length(geneRecurrencyDfColumns)))
  colnames(geneRecurrencyDf)<-geneRecurrencyDfColumns
  
  ###############################################################
  ## For each of the considered phenotypes the analysis are run
  ###############################################################
  
  for(target_pheno in consideredPheno){
    ##maybe error rised only at specific phenos
    sv_phenoSpecific_withErrorRised<-character()
    #To track in same way No gene found
    sv_phenoSpecific_with_NO_geneFound<-character()
    
    # print(target_pheno)
    #target_pheno<-"head_neck"
    ##get considered phases
    ################################
    ##Results Matrices preparation
    ################################
    namesPhases<-pheno_Phases[[target_pheno]]##phases associated to phenotype of interest
    nPhases<-length(namesPhases)
    
    features_summaryMatrix<-c(namesPhases,"patients")
    
    ##Building results matrix
    backBone_summaryMatrix<-as.data.frame(matrix(data=NA, nrow=0, ncol=length(features_summaryMatrix)))
    colnames(backBone_summaryMatrix)<-features_summaryMatrix
    
    ##Object to store all results for a phenotype
    cohort_results<-list("anyMechanism"=backBone_summaryMatrix,
                         "DirectEffectLOF"=backBone_summaryMatrix,
                         "DirectEffectGOF"=backBone_summaryMatrix,
                         "LongRangeLOF"=backBone_summaryMatrix,
                         "LongRangeGOF"=backBone_summaryMatrix)
    
    ######################################
    ##Matrix to store All patients info
    featuresConsideredInfoPatient<-colnames(AllPatientsInfo)##Here you have patients from many different types of phenotypes, but colnames the same
    df_allPatientsInfo<-as.data.frame(matrix(data=NA, nrow=0, ncol=length(featuresConsideredInfoPatient)))
    colnames(df_allPatientsInfo)<-featuresConsideredInfoPatient
    
    ###Adding 4 additional columns to the anyMechanism//Overview table
    cohort_results$anyMechanism<-cbind(cohort_results$anyMechanism,
                                       "DirectEffectLOF"=numeric())
    cohort_results$anyMechanism<-cbind(cohort_results$anyMechanism,
                                       "LongRangeLOF"=numeric())
    cohort_results$anyMechanism<-cbind(cohort_results$anyMechanism,
                                       "DirectEffectGOF"=numeric())
    cohort_results$anyMechanism<-cbind(cohort_results$anyMechanism,
                                       "LongRangeGOF"=numeric())
    
    for(patient in names(all_patientResults)){
      ##Name of patient is SV id. if >1 SV think on how handling
      
      ##If patient presents the currently considered phenotype
      
      patientResults<-all_patientResults[[patient]]
      
      if(target_pheno %in% names(patientResults)){
        patientResults<-patientResults[[target_pheno]]
        #So patient falls in the current phenotype
        ##Add to all patients info dataframe
        df_allPatientsInfo<-rbind(df_allPatientsInfo,
                                  AllPatientsInfo[patient,])
        
        if(patientResults$Status=="OK"){
          
          
          ##Tracking to which category we have already added the mechanism, for the AnyMechanism Statistics//Summary
          ##So if a gene-mech already considered for a patient, not reconsider
          mechOverview<-character() ## DirectEffectLof||DirectEffectGOF||LongRangeLOF||LongRangeGOF
          genePatient_trackingMechOverview<-character() ## DirectEffectLof_Gene||DirectEffectGOF_Gene||LongRangeLOF_Gene||LongRangeGOF_Gene
          
          ###Getting Candidate Genes
          masterSummaryResMat<-patientResults$masterSummaryResultsMatrix
          candidateGenes<-unique(masterSummaryResMat$affected_gene)
          max_PathoScore<-max(masterSummaryResMat$maxScore, na.rm = TRUE)
          
          
          masterSummaryResMat$combGeneScore<-paste0(masterSummaryResMat$affected_gene,
                                                    "(",
                                                    masterSummaryResMat$maxScore,
                                                    ")")
          
          candidateGenesAndPathoScore<-paste0(masterSummaryResMat$combGeneScore,sep="", collapse = ",")

          
          ##Adding info to candidateGenesInfo Matrix
          targetRow<-nrow(candidateGenesInfo)+1
          
          candidateGenesInfo[targetRow,"SVid"]<-patient
          candidateGenesInfo[targetRow,"Phenotype"]<-target_pheno
          candidateGenesInfo[targetRow,"Ncandidates"]<-length(candidateGenes)
          candidateGenesInfo[targetRow,"Candidates"]<-candidateGenesAndPathoScore
          candidateGenesInfo[targetRow,"PathogenicScore"]<-max_PathoScore
          candidateGenesInfo[targetRow,"TypeSV"]<-AllPatientsInfo[patient,"TypeSV"]
          
          if(max_PathoScore>=0.8){
            candidateGenesInfo[targetRow,"Pathogenic"]<-"Yes"  
          }else{
            candidateGenesInfo[targetRow,"Pathogenic"]<-"No"  
          }
          
          
          filt_pathog<-subset(patientResults$masterSummaryResultsMatrix, maxScore>=minScore)
          causativeGenes<-unique(filt_pathog$affected_gene)
          
          ##Adding info causative (if there are not this will be 0 and empty)
          candidateGenesInfo[targetRow,"Ncausative"]<-length(causativeGenes)
          candidateGenesInfo[targetRow,"Causative"]<-paste0(causativeGenes, sep="",collapse = "," )
          
          ##By default if no pathog mech predicted
          N_LR_Mech<-0
          N_Direct_Mech<-0
          ##CHECK N ROWS WITH AN IF
          if(nrow(filt_pathog)>=1){
            ## So, at least 1 disease causative gene predicted
            ##Getting causative genes and pathomech
            causativeGenes<-filt_pathog$affected_gene
            
            ##Getting pathological mechanisms
            pathologicalMechanisms<-filt_pathog$GeneImpact
            ##First part of pathomech terms is either Direct or LongRange
            pathologicalMechanisms<-unlist(lapply(pathologicalMechanisms,FUN = function(x){
              return(unlist(strsplit(x=x, split = "_", fixed = TRUE))[1])
            }))
            
            ##Add causative genes, check this after or tomorrow
            N_LR_Mech<-sum(pathologicalMechanisms == "LongRange")
            N_Direct_Mech<-sum(pathologicalMechanisms == "Direct")
            
          }
          
          ##Adding data
          candidateGenesInfo[targetRow,"N_LR_Mech"]<-N_LR_Mech
          candidateGenesInfo[targetRow,"N_Direct_Mech"]<-N_Direct_Mech
          
          
          ##########################################################
          ## Parsing results matrix
          ######################################################
          ## Adding info to gene-mech table
          
          if(nrow(filt_pathog)>=1){
            ## So, at least 1 disease causative gene predicted
            ##Getting causative genes and pathomech
            causativeGenes<-filt_pathog$affected_gene
            
            for(DCgene in causativeGenes){
              targetInfo_filt_pathog<-filt_pathog[DCgene,]
              
              ##Getting pathological mechanisms
              t_pathoMech<-targetInfo_filt_pathog$GeneImpact
              
              
              ##Check case "Direct_LongRange_geneDuplication", aqui sumaremos en los dos, coding and long.range pq o not sure o ambos predichos
              
              if(t_pathoMech != "Direct_LongRange_geneDuplication"){
                
                ##Okey, so, only Direct or LongRange will be
                ##Extraemos esa info
                ##First part of pathomech terms is either Direct or LongRange
                t_pathoMech<-unlist(strsplit(x=t_pathoMech, split = "_", fixed = TRUE))[1]
              }else{
                ##t_pathomech is Direct_LongRange_geneDuplication
                ##Desgranar, si patog predicha solo via longRange...
                
                # "LongRange_geneDuplication", "Direct_geneDuplication", "Direct_LongRange_geneDuplication"
                
                # getPathoMechs
                
                process_filt_pathog<-filt_pathog
                process_filt_pathog$affected_gene<-NULL
                process_filt_pathog$maxScore<-NULL
                process_filt_pathog$GeneImpact<-NULL
                ##Solo quedan las fases
                trackPathoMech<-character()
                for(tphase in colnames(process_filt_pathog)){
                  
                  tphaseInfo<-process_filt_pathog[[tphase]]
                  for(tMech in 1:2){
                    
                    tInfo<-unlist(strsplit(tphaseInfo, split = "--", fixed = TRUE))[tMech]
                    tInfo_score<-as.numeric(unlist(strsplit(tInfo, split = ";", fixed = TRUE))[1])
                    tInfo_mech<-unlist(strsplit(tInfo, split = ";", fixed = TRUE))[3]
                    
                    if(!is.na(tInfo_score)){
                      if(tInfo_score >= minScore){
                        
                        trackPathoMech<-c(trackPathoMech, tInfo_mech)
                        trackPathoMech<-unique(trackPathoMech)
                      }
                    }
                  }
                }
                
                
                
                # if(length(trackPathoMech) > 1)){
                #   ##Two mechs predicted, so. t_mech se queda como esta
                #   # t_pathoMech <- "Direct_LongRange_geneDuplication"                  
                # }

                 if (length(trackPathoMech) == 1){
                  ##Only 1 mech predicted, if it is only direct or only long-range only put direct or long range in output table
                   if(trackPathoMech=="LongRange_geneDuplication"){
                     t_pathoMech<-"LongRange"
                   }else{
                     t_pathoMech<-"Direct"
                   }
                }
                
              }
              
              rowId<-paste0(DCgene,"_",target_pheno)
              if(rowId %in% rownames(geneRecurrencyDf)){
                ##Adding info, over existing, so increasing numbers
                geneRecurrencyDf[rowId,"N SVs"]<-geneRecurrencyDf[rowId,"N SVs"] + 1
                geneRecurrencyDf[rowId,"SVid"]<-paste0(c(geneRecurrencyDf[rowId,"SVid"],patient), collapse=",")
                
                if(t_pathoMech=="LongRange"){
                  geneRecurrencyDf[rowId,"Long-Range"]<-paste0(c(na.omit(geneRecurrencyDf[rowId,"Long-Range"]),patient), collapse=",")
                  geneRecurrencyDf[rowId,"N Long-Range"]<-geneRecurrencyDf[rowId,"N Long-Range"]+1
                }else if(t_pathoMech=="Direct"){
                  geneRecurrencyDf[rowId,"Coding"]<-paste0(c(na.omit(geneRecurrencyDf[rowId,"Coding"]),patient), collapse=",")
                  geneRecurrencyDf[rowId,"N Coding"]<-geneRecurrencyDf[rowId,"N Coding"]+1
                  
                }else if(t_pathoMech=="Direct_LongRange_geneDuplication"){
                  geneRecurrencyDf[rowId,"Long-Range"]<-paste0(c(na.omit(geneRecurrencyDf[rowId,"Long-Range"]),patient), collapse=",")
                  geneRecurrencyDf[rowId,"Coding"]<-paste0(c(na.omit(geneRecurrencyDf[rowId,"Coding"]),patient), collapse=",")
                  
                  geneRecurrencyDf[rowId,"N Long-Range"]<-geneRecurrencyDf[rowId,"N Long-Range"]+1
                  geneRecurrencyDf[rowId,"N Coding"]<-geneRecurrencyDf[rowId,"N Coding"]+1
                }
              
              }else{
                # browser()
                ##Adding de novo
                ##Adding info to candidateGenesInfo Matrix
                targetRow<-nrow(geneRecurrencyDf)+1
                geneRecurrencyDf[targetRow,]<-NA
                rownames(geneRecurrencyDf)[targetRow]<-rowId
                
                geneRecurrencyDf[rowId,"Gene"]<-DCgene                
                geneRecurrencyDf[rowId,"Phenotype"]<-target_pheno
                geneRecurrencyDf[rowId,"N SVs"]<-1
                geneRecurrencyDf[rowId,"SVid"]<-patient
                
                if(t_pathoMech=="LongRange"){
                  geneRecurrencyDf[rowId,"Long-Range"]<-patient
                  geneRecurrencyDf[rowId,"Coding"]<-NA
                  geneRecurrencyDf[rowId,"N Long-Range"]<-1
                  geneRecurrencyDf[rowId,"N Coding"]<-0
                  
                }else if(t_pathoMech=="Direct"){
                  geneRecurrencyDf[rowId,"Coding"]<-patient
                  geneRecurrencyDf[rowId,"Long-Range"]<-NA
                  geneRecurrencyDf[rowId,"N Long-Range"]<-0
                  geneRecurrencyDf[rowId,"N Coding"]<-1
                  
                }else if(t_pathoMech=="Direct_LongRange_geneDuplication"){
                  geneRecurrencyDf[rowId,"Long-Range"]<-patient
                  geneRecurrencyDf[rowId,"Coding"]<-patient
                  geneRecurrencyDf[rowId,"N Long-Range"]<-1
                  geneRecurrencyDf[rowId,"N Coding"]<-1
                }
                
              }
            }
          }
          
          
          for(phase in namesPhases){
            # print(phase)
            # print(patient)
            ##Per phase 
            ##phaseResults<-patientResults$resultsPerPhase[[phase]][[phase]]
            phaseResults<-patientResults$masterSummaryResultsMatrix[,c("affected_gene",phase)]
            
            ##Rows with a cell with a NA need to be removed, they come from genes not considered on the current phase
            ##Due to TAD changes dependent on stage
            ##remove them to avoid downstream problems
            rowsToMantain<-!is.na(phaseResults[,phase])
            phaseResults<-phaseResults[rowsToMantain,]
            
            ##Extract scores for each gene
            ##Exctract GOF and LOF, per separate
            phaseResults$LOF_Score<-as.numeric(unlist(lapply(X = phaseResults[[phase]],
                                                             FUN = function(x){
                                                               lof_info<-unlist(strsplit(x=x, split = "--", fixed = TRUE))[1]
                                                               lof_score<-unlist(strsplit(x=lof_info, split = ";", fixed = TRUE))[1]
                                                               return(lof_score)
                                                             })))
            
            phaseResults$GOF_Score<-as.numeric(unlist(lapply(X = phaseResults[[phase]],
                                                             FUN = function(x){
                                                               gof_info<-unlist(strsplit(x=x, split = "--", fixed = TRUE))[2]
                                                               gof_score<-unlist(strsplit(x=gof_info, split = ";", fixed = TRUE))[1]
                                                               return(gof_score)
                                                             })))
            
            ##Mechanism : LOF, GOF
            ##Deprecated since now we have info for both, GOF and LOF, so info for both considered
            # phaseResults$Mechanism<-unlist(lapply(X = phaseResults[[phase]],
            #                                                      FUN = function(x){unlist(strsplit(x=x, split = ";", fixed = TRUE))[2]}))
            
            ##GeneImpact: LONG_RANGE...
            ##Will be the same, but just as a matter of time doing duplicated
            phaseResults$GeneImpact_LOF<-unlist(lapply(X = phaseResults[[phase]],
                                                       FUN = function(x){
                                                         lof_info<-unlist(strsplit(x=x, split = "--", fixed = TRUE))[1]
                                                         lof_impact<-unlist(strsplit(x=lof_info, split = ";", fixed = TRUE))[3]
                                                         return(lof_impact)
                                                       }))
            
            phaseResults$GeneImpact_GOF<-unlist(lapply(X = phaseResults[[phase]],
                                                       FUN = function(x){
                                                         gof_info<-unlist(strsplit(x=x, split = "--", fixed = TRUE))[2]
                                                         gof_impact<-unlist(strsplit(x=gof_info, split = ";", fixed = TRUE))[3]
                                                         return(gof_impact)
                                                         }))
            
            ##Subset results for relevant minScore
            ##Mantain only those with at least 1 score (for gof or lof > minScore)
            phaseResults<-subset(phaseResults, LOF_Score >= minScore | GOF_Score >= minScore )
          
            ##If true, if among the genes there is one directly broken or deleted With a high score
            ##Skip the patient (this will be used for when looking for only long-range  mechanisms)
            if(discardRelevantByBrokenGene==TRUE){
              if(("Direct_geneTruncation" %in% c(phaseResults$GeneImpact_LOF,phaseResults$GeneImpact_GOF)) ||
                 ("Direct_geneDeletion" %in% c(phaseResults$GeneImpact_LOF,phaseResults$GeneImpact_GOF))){
                ##We empty the matrix so no relevant info to keep digging
                phaseResults<-phaseResults[0,]
              }
            }
            
            #############################################################################
            ##Take into account that both backbone matrix and phase results can be empty
            
            if(nrow(phaseResults)>=1){
              ##So at least info for one gene
              for(gene in phaseResults$affected_gene){
                
                ##Tengo que iterar what follows
                gene_scores<-c(phaseResults[gene,"LOF_Score"], phaseResults[gene,"GOF_Score"])
                names(gene_scores)<-c("LOF","GOF")
                
                relevantMechanisms<-names(gene_scores)[gene_scores>=minScore]
                
                for(gene_mechanism in relevantMechanisms){
                  
                  gene_impact<-unlist(strsplit(phaseResults[gene,paste0("GeneImpact_", gene_mechanism)],split = "_", fixed = TRUE))[1]##Hence, LongRange -> longrange pathomechanisms; direct -> direct impact gene deletion, duplication, truncation...
                  
                  
                  ##############################################
                  ##Adding info to parsing result Matrices #####
                  ##############################################
                  
                  ######################################
                  ##First adding to anyMechanism matrix
                  cohort_results$anyMechanism<-addingInfoToMatrix(gene=gene, 
                                                                  targetMatrix = cohort_results$anyMechanism, ## RELEVANT PARAMETER TO CHANGE
                                                                  phase = phase, 
                                                                  patientId = patient)
                  
                  ###############################################################################
                  ##Now adding to the corresponding dataframe associated to the pathomechanism
                  
                  ##Long-range LOF
                  if((gene_impact == "LongRange") && (gene_mechanism == "LOF")){
                    cohort_results$LongRangeLOF<-addingInfoToMatrix(gene=gene, 
                                                                    targetMatrix = cohort_results$LongRangeLOF, ## RELEVANT PARAMETER TO CHANGE
                                                                    phase = phase, 
                                                                    patientId = patient)
                    ##track mechanism
                    #If this mechanism not already considered for this patient&gene, add it to the anyMech overview table
                    currentTargetMech<-"LongRangeLOF"#variable to re-use the following code
                    current_genePatient_TargetMech<-paste(currentTargetMech,
                                                          "_",
                                                          gene,
                                                          sep="")
                    
                    if(!(current_genePatient_TargetMech %in% genePatient_trackingMechOverview)){
                      
                      ##Add +1 to anyMech table. It already has an entrance, as it has been already added most of the info
                      if(is.na(cohort_results$anyMechanism[gene,currentTargetMech])){
                        cohort_results$anyMechanism[gene,currentTargetMech]<-1
                      }else{
                        ##so already one patient with this mechanism for the same gene
                        cohort_results$anyMechanism[gene,currentTargetMech]<-cohort_results$anyMechanism[gene,currentTargetMech] + 1
                      }
                      
                      #######################################################################
                      ##Track it on the mechOverview, to not count it >1 for the same patient
                      genePatient_trackingMechOverview<-c(genePatient_trackingMechOverview, 
                                                          current_genePatient_TargetMech)
                    }
                    
                  }
                  
                  ##Long-range GOF
                  if((gene_impact == "LongRange") && (gene_mechanism == "GOF")){
                    cohort_results$LongRangeGOF<-addingInfoToMatrix(gene=gene, 
                                                                    targetMatrix = cohort_results$LongRangeGOF, ## RELEVANT PARAMETER TO CHANGE
                                                                    phase = phase, 
                                                                    patientId = patient)
                    ##track mechanism
                    #If this mechanism not already considered for this patient&gene, add it to the anyMech overview table
                    currentTargetMech<-"LongRangeGOF"#variable to re-use the following code
                    current_genePatient_TargetMech<-paste(currentTargetMech,
                                                          "_",
                                                          gene,
                                                          sep="")
                    
                    if(!(current_genePatient_TargetMech %in% genePatient_trackingMechOverview)){
                      
                      ##Add +1 to anyMech table. It already has an entrance, as it has been already added most of the info
                      if(is.na(cohort_results$anyMechanism[gene,currentTargetMech])){
                        cohort_results$anyMechanism[gene,currentTargetMech]<-1
                      }else{
                        ##so already one patient with this mechanism for the same gense
                        cohort_results$anyMechanism[gene,currentTargetMech]<-cohort_results$anyMechanism[gene,currentTargetMech] + 1
                      }
                      
                      #######################################################################
                      ##Track it on the mechOverview, to not count it >1 for the same patient
                      genePatient_trackingMechOverview<-c(genePatient_trackingMechOverview, 
                                                          current_genePatient_TargetMech)
                    }
                    
                    
                  }
                  
                  ##Direct LOF
                  if((gene_impact == "Direct") && (gene_mechanism == "LOF")){
                    cohort_results$DirectEffectLOF<-addingInfoToMatrix(gene=gene, 
                                                                       targetMatrix = cohort_results$DirectEffectLOF, ## RELEVANT PARAMETER TO CHANGE
                                                                       phase = phase, 
                                                                       patientId = patient)
                    ##track mechanism
                    #If this mechanism not already considered for this patient&gene, add it to the anyMech overview table
                    currentTargetMech<-"DirectEffectLOF"#variable to re-use the following code
                    current_genePatient_TargetMech<-paste(currentTargetMech,
                                                          "_",
                                                          gene,
                                                          sep="")
                    
                    if(!(current_genePatient_TargetMech %in% genePatient_trackingMechOverview)){
                      
                      ##Add +1 to anyMech table. It already has an entrance, as it has been already added most of the info
                      if(is.na(cohort_results$anyMechanism[gene,currentTargetMech])){
                        cohort_results$anyMechanism[gene,currentTargetMech]<-1
                      }else{
                        ##so already one patient with this mechanism for the same gense
                        cohort_results$anyMechanism[gene,currentTargetMech]<-cohort_results$anyMechanism[gene,currentTargetMech] + 1
                      }
                      
                      #######################################################################
                      ##Track it on the mechOverview, to not count it >1 for the same patient
                      genePatient_trackingMechOverview<-c(genePatient_trackingMechOverview, 
                                                          current_genePatient_TargetMech)
                    }
                    
                    
                  }
                  ##Direct GOF
                  if((gene_impact == "Direct") && (gene_mechanism == "GOF")){
                    cohort_results$DirectEffectGOF<-addingInfoToMatrix(gene=gene, 
                                                                       targetMatrix = cohort_results$DirectEffectGOF, ## RELEVANT PARAMETER TO CHANGE
                                                                       phase = phase, 
                                                                       patientId = patient)
                    ##track mechanism
                    #If this mechanism not already considered for this patient&gene, add it to the anyMech overview table
                    currentTargetMech<-"DirectEffectGOF"#variable to re-use the following code
                    current_genePatient_TargetMech<-paste(currentTargetMech,
                                                          "_",
                                                          gene,
                                                          sep="")
                    
                    if(!(current_genePatient_TargetMech %in% genePatient_trackingMechOverview)){
                      
                      ##Add +1 to anyMech table. It already has an entrance, as it has been already added most of the info
                      if(is.na(cohort_results$anyMechanism[gene,currentTargetMech])){
                        cohort_results$anyMechanism[gene,currentTargetMech]<-1
                      }else{
                        ##so already one patient with this mechanism for the same gense
                        cohort_results$anyMechanism[gene,currentTargetMech]<-cohort_results$anyMechanism[gene,currentTargetMech] + 1
                      }
                      
                      #######################################################################
                      ##Track it on the mechOverview, to not count it >1 for the same patient
                      genePatient_trackingMechOverview<-c(genePatient_trackingMechOverview, 
                                                          current_genePatient_TargetMech)
                    }
                  }
                }
              }
            }
          }
          
          ##Substitute NAs per 0s
          ##And add affected_gene column
          for(dataName in names(cohort_results)) {
            #print(dataName)
            matrix_dataName<-cohort_results[[dataName]]
            
            
            matrix_dataName[is.na(matrix_dataName)] <- 0
            matrix_dataName$affected_gene<-rownames(matrix_dataName)
            
            ##add nPatients column
            matrix_dataName$Num_Patients <- as.numeric(lapply(X = strsplit(x = as.character(matrix_dataName$patients),
                                                                         split = ",",fixed = TRUE),
                                                            FUN = length))
            
            ###Resorting columns
            # firstColumns<-c("affected_gene","Num_Patients","patients")
            firstColumns<-c("affected_gene","Num_Patients")## patients last column
            lastColumns<-c("patients")
            notFirstColumns<-colnames(matrix_dataName)[!(colnames(matrix_dataName) %in% c(firstColumns,lastColumns))]
            matrix_dataName<-matrix_dataName[,c(firstColumns,notFirstColumns,lastColumns)]
            
            ##Sort by nPatients column decreasing, and then by gene name increasing (alphabetical order)
            matrix_dataName<-matrix_dataName[order(-matrix_dataName$Num_Patients, matrix_dataName$affected_gene),]
            
            cohort_results[[dataName]]<-matrix_dataName
            
          } 
        }else if(patientResults$Status=="ERROR"){
          sv_withErrorRised<-c(sv_withErrorRised,patient)##Error in General, does not matter at which pheno
          sv_phenoSpecific_withErrorRised<-c(sv_phenoSpecific_withErrorRised, patient)##Errors specific of a Phenotype
          
        }else if(patientResults$Status=="OK, but NO genes associated with SV"){
          ##No Gene Associated with the SV
          ##DO NOTHING, NO information to add to any matrix
          ##Track info, info is power ;)
          sv_with_NO_geneFound<-c(sv_with_NO_geneFound,patient)##Not found in General, does not matter at which pheno
          sv_phenoSpecific_with_NO_geneFound<-c(sv_phenoSpecific_with_NO_geneFound, patient)##Not found, specific of a Phenotype
          
          ##Adding info to candidateGenesInfo Matrix
          targetRow<-nrow(candidateGenesInfo)+1
          
          candidateGenesInfo[targetRow,"SVid"]<-patient
          candidateGenesInfo[targetRow,"Phenotype"]<-target_pheno
          candidateGenesInfo[targetRow,"Ncandidates"]<-0
          candidateGenesInfo[targetRow,"Candidates"]<-""
          candidateGenesInfo[targetRow,"PathogenicScore"]<-0
          candidateGenesInfo[targetRow,"TypeSV"]<-AllPatientsInfo[patient,"TypeSV"]
          
          candidateGenesInfo[targetRow,"Pathogenic"]<-"No"  

          
        }
      }else{
        ##The current SV do not related with currently studied phenotype so, we do not take it into consideration
        #do nothing
      }
    }  
    ##For any mechanism remove phases column. Only interested in overall statistics
    cohort_results$anyMechanism[,namesPhases]<-NULL
    
    ##Once all SV associated to a pheno analyzed, they are introduced here
    perPheno_cohortResults[[target_pheno]]<-cohort_results
    ##We add also the all patient info associated to this phenotype
    ##resorting patients info
    firstCol<-"patientID"
    notFirstCol<-colnames(df_allPatientsInfo)[colnames(df_allPatientsInfo)!=firstCol]
    df_allPatientsInfo<-df_allPatientsInfo[,c(firstCol,notFirstCol)]
    perPheno_cohortResults[[target_pheno]][["patientsInfo"]]<-df_allPatientsInfo
    ##Error Info
    perPheno_cohortResults[[target_pheno]][["errorInfo"]]<-unique(sv_phenoSpecific_withErrorRised)##error associated to the phenotype
    
    ##No Gene Found Info
    perPheno_cohortResults[[target_pheno]][["noGeneFoundInfo"]]<-unique(sv_phenoSpecific_with_NO_geneFound)## associated to the phenotype
  }
  
  ##Adding error info
  perPheno_cohortResults[["error_General_Info"]]<-unique(sv_withErrorRised)##errors happening wherever, at any phenotype
  
  ##Adding No Gene Found Info
  perPheno_cohortResults[["noGeneFound_General_Info"]]<-unique(sv_with_NO_geneFound)##errors happening wherever, at any phenotype
  
  
  ##Adding Info NCandidate genes per SV
  #Let's rename the Candidates column, and also the causative genes one
  colnames(candidateGenesInfo)[which(colnames(candidateGenesInfo)=="Candidates")]<-"Candidate genes (Pathogenic Score)"
  colnames(candidateGenesInfo)[which(colnames(candidateGenesInfo)=="Causative")]<-"Causative genes"
  
  perPheno_cohortResults[["candidateGenesInfo"]]<-candidateGenesInfo
  
  ##Adding Info Recurrency Affected genes in different SVs
  # browser()
  ##Ordering per Recurrency column
  geneRecurrencyDf<-geneRecurrencyDf[order(geneRecurrencyDf$`N SVs`, decreasing = TRUE),]
  #Eliminar NAs
  geneRecurrencyDf[is.na(geneRecurrencyDf)] <- ""
  
  perPheno_cohortResults[["geneRecurrencyInfo"]]<-geneRecurrencyDf
  
  ##Returning patient parsed info
  return(perPheno_cohortResults)
}
