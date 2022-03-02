##########################
##Parsing Cohort Results
##########################
# setwd("~/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/SV_app")
##Filter for patients with no broken gene considered as harmful
##So only possible explanation long-range mechanism
##So if a patient has a gene broken, with a score > minScore. If we set discardRelevantBroken=TRUE
##The information for that patient will be skipped


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
  
  ###############################################################
  ## For each of the considered phenotypes the analysis are run
  ###############################################################
  for(target_pheno in consideredPheno){
    ##maybe error rised only at specific phenos
    sv_phenoSpecific_withErrorRised<-character()
    #To track in same way No gene found
    sv_phenoSpecific_with_NO_geneFound<-character()
    
    #print(target_pheno)
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
          
          for(phase in namesPhases){
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
  
  ##Returning patient parsed info
  return(perPheno_cohortResults)
}
