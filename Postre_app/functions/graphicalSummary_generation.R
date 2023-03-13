##############################
## Trying Graphical Summary ##
##############################

################
## IMPORTANT
##############

#Upon image adjustment creation, adjust padding to 0 to better fit in the html page

# Required libraries
# library(plotrix)##For enhancers, ellipse shape representation
# library(shape)##For curve arrows representation
# library(diagram)##For curve arrows representation

##Required additional functions
source("functions/aux_functions_forGraphicalSummary.R", local = TRUE) ##For WT TAD plots
source("functions/makingGraphs/aux_fun_Plots_Rearranged_TADs.R", local = TRUE) ##For rearranged TAD plots

##By default, only computing report for genes whose score > threshold, if not it takes too much time and to many plots

graphicalSummary_generation<-function(patientResults, minPathogenicScore){
  ##min pathogenic scores filters out genes non relevant (below threshold) so their graphical summary it is not generated
  
  ##To track genes & conditions for which report will be generated
  patientResults$genesConditions_ToReport<-character()
  
  # ################
  # # WATCHOUT
  # #For now, each time this function is called it is removed the current content of graphSummary folder
  # #But for future take into account if this is problematic, when multiple users simultaneously connected
  # #For now is not rising problems on shiny server
  # 
  # ##check if any image present in the folder, if there are, delete them
  # files_paths<-list.files(path = 'www/graphicalSummaries/')
  # if(length(files_paths) !=0){
  #   ##Hence there is at least an image to delete
  #   system('rm www/graphicalSummaries/*')  
  # }
  
  allGenes<-patientResults$allAffectedGenes_positionalInfo$geneSymbol
  allPhases<-names(patientResults$resultsPerPhase)
  
  # gene<-"TFAP2A"
  # phase<-"NeuralCrestLate"
  # minPathogenicScore<-0.8
  
  ####################
  ## Patient SV type
  ## Depending on the Type of SV a different type of plot is generated
  
  sv_type<-patientResults$patientInfo$TypeSV
  
  ##Breakpoints options, get them, used to get coord of opposite domain
  ##Take into account we can also have the option: "No_AssociatedBreakpoint"
  ##We will only work and get coordinates of opposite domain, when current gene associated with one break
  #breakpOptions<-unique(patientResults$resultsPerPhase[[1]][[1]]$gene_Breakpoint)##break1, break2
  binary_breakpOptions<-c("Break1", "Break2")
  
  
  
  # if((sv_type=="Inversion") || (sv_type=="Translocation")){
  for(gene in allGenes){
    
    ###############################################
    ##Getting data stage specific and doing plots
    ###############################################
    for(phase in allPhases){
      
      resultsPhase<-patientResults$resultsPerPhase[[phase]][[phase]]
      
      if(gene %in% resultsPhase$affected_gene){
        ##If gene score is above interesting, continue with all the calculations
        ##On the contrary, skip it
        phaseIntegratedResults<-resultsPhase
  
        for(targetMech in c("LOF","GOF")){
          ##Screen independently LOF and GOF results
          
          targetGeneMech<-paste(gene,"--",targetMech, sep = "")
          
          ##Phase results are provided for both LOF and GOF, so narrowing gene-Mechanism
          ##And all info retrieved from phaseIntegratedResults, gene-Mechanism!
          phase_score<-phaseIntegratedResults[targetGeneMech,"Average_Score"]
          gene_breakpoint<-phaseIntegratedResults[targetGeneMech,"gene_Breakpoint"]
          gene_mechanism<-phaseIntegratedResults[targetGeneMech,"GeneImpact"]
          SV_landing<-phaseIntegratedResults[targetGeneMech,"SV_landing"]
          
          gene_TSS<-patientResults$allAffectedGenes_positionalInfo[gene,"TSS"]
          chr_gene<-patientResults$allAffectedGenes_positionalInfo[gene,"chr"]
          
          # Deprecated since we now scan results for both
          # pathoMechanism<-resultsPhase[gene,"Mechanism"]
          # pathoMechanism<-unlist(strsplit(x = pathoMechanism,
          #                                 split = ":",
          #                                 fixed = TRUE))[2]##GOF or LOF or GOF-LOF
          pathoMechanism<-targetMech
          
          if(phase_score>=minPathogenicScore){
            ###############################################################################
            ###Okey, the gene is worth to plot, so continue with the plotting calculations
            ##And also track it to afterwards generate the gene associated report
            
            ##Add here the gene relevant tracking info
            ##If pathoMechanism is inconclusive (hence "MIX:GOF-LOF" we do not want to explain it because really uncertain, eventhough score > minScore)
            ##We also want to track directly affected genes, eventhough not significant score, proably a user will want to check them
            
            ##Track this gene, as a report will be generated for it
            geneConditionId<-paste(gene,
                                   "_",
                                   targetMech,
                                   "_",
                                   phase,
                                   sep="")
            
            patientResults$genesConditions_ToReport<-c(patientResults$genesConditions_ToReport, 
                                                       geneConditionId)
            
            
            
            ###Important here to know, whether the gene regulatory environment has been affected or not
            ##We know if has not been affected (so TAD entirely dup, or deleted if:gene_breakpoint=="No_AssociatedBreakpoint")
            ##If gene does not have an associated breakpoint, makes no sense to talk about an opposite domain
            ##And same applies if the gene is associated to both breakpoints, because the SV is intraTAD
            
            ##Diff options available:
            ## (there always could be more that I'm not picturing right now)
            ##If gene associated to just one break, break1 or break2, there is a TAD reshuffling so it makes sense to think about opposite domain (in this case breakpoint painted inside the TAD)
            ##If gene associated to none, is because TAD entirely Deleted or Duplicated (breakpoints painted outside the TAD)
            ##If gene associated with both (Break1&2) is because intraTAD SV (in this case bp painted inside of TAD both before or after gene)
            
            if(gene_breakpoint %in% binary_breakpOptions){
              ##Retrieve info too of opposite domain
              ##Here apply this if gene break is break1 or break2
              oppositeDomain_breakpoint<-binary_breakpOptions[binary_breakpOptions != gene_breakpoint]
              
              gene_breakpointCoord<-patientResults$patientInfo[[paste0("coord_",gene_breakpoint)]]
              oppositeDomain_breakpointCoord<-patientResults$patientInfo[[paste0("coord_",oppositeDomain_breakpoint)]]
              
              ##I'm gonna use just the start of the breakpoint as reference to determine where to place it
              gene_breakpointCoord_start<-as.integer(unlist(strsplit(x = gene_breakpointCoord,
                                                                     split = ",", fixed = TRUE))[1])
              
              oppositeDomain_breakpointCoord_start<-as.integer(unlist(strsplit(x = oppositeDomain_breakpointCoord,
                                                                               split = ",", fixed = TRUE))[1])
              # chr_oppositeDomain
              chromosomes<-unique(patientResults$allAffectedGenes_positionalInfo$chr)
              if(length(chromosomes)==1){
                chr_oppositeDomain<-chr_gene
              }else{
                #So it is a translocation, hence >1 chr involved
                chr_oppositeDomain<-chromosomes[chromosomes != chr_gene]
              }
              
              ############################################################################
              ## Defining where do the WT tad with the gene is painted
              ## And the relative position of the breakpoint with respect to the gene TSS
              ############################################################################
              
              ##First, relative position of gene breakpoint with respect to the gene TSS
              if(gene_breakpointCoord_start < gene_TSS){
                geneBreakP_Position_respectToTSS <-"beforeTSS"
              }else{
                geneBreakP_Position_respectToTSS <-"afterTSS"
              }
              
              ##Defining where do the WT tad with the gene is painted
              
              ## If the gene breakpoint is larger than the opposite breakpoint, it goes to the right side
              ## and the other way round
              ## For SV happening in the same chr it is the way that should be. 
              ## And for translocations it is probably the easiest way to understand if we consider that more than one gene can jump as significant
              ## And if that happens it can be messy that two gens from same TAD are painted on different TAD primary positions,
              ## Based on where do they have their breakpoints in their TAD
              ##So here we respect genome positioning
              
              if(gene_breakpointCoord_start>oppositeDomain_breakpointCoord_start){
                situation<-"primaryTAD_Dextral"##Painted at the right side of the window
              }else{
                situation<-"primaryTAD_Sinistral"##Painted at the left side of the window 
              }
              
            }else if(gene_breakpoint == "No_AssociatedBreakpoint"){
              ########################################################
              ##For entirely duplicated or deleted regulatory domains
              ########################################################
              geneBreakP_Position_respectToTSS<-"outOf_RegulatoryDomain"
              situation<-"primaryTAD_Central"
              
            }else if(gene_breakpoint == "Break1&2"){
              #######################
              ##For intraTAD SV
              #######################
              
              ## Regarding geneBreakP_Position_respectToTSS
              ## The only situation important is for long-range. Which we need to know if we paint it to the right or to the left
              ## Because if deleted... lines painted based on that, around the gene figure
              ## And same goes for Truncated
              ## And same goes for Duplicated
              
              if(gene_mechanism=="LongRange"){
                ##So for the case... just work with Breakpoint1 coord,
                ##Since in long-range IntraTAD cases both will be either to the left or to the right
                gene_breakpointCoord<-patientResults$patientInfo[[paste0("coord_","Break1")]]
                
                ##I'm gonna use just the start of the breakpoint as reference to determine where to place it
                gene_breakpointCoord_start<-as.integer(unlist(strsplit(x = gene_breakpointCoord,
                                                                       split = ",", fixed = TRUE))[1])
                if(gene_breakpointCoord_start < gene_TSS){
                  geneBreakP_Position_respectToTSS <-"beforeTSS"
                }else{
                  geneBreakP_Position_respectToTSS <-"afterTSS"
                }
                
              }else{
                ##Just generating this variable to avoid missing variable errors
                ##On function call eventhough won't be used
                geneBreakP_Position_respectToTSS<-"Unnecessary"
              }
              
              ##Since both breakpoints on its regulatory Domain... we only paint one reg domain
              situation<-"primaryTAD_Central"
              
            }else if(gene_breakpoint == "NotClear_BreakpointAssociation"){
              ##Not plot for now in this scenario, but doing sth to avoid app crash
              ##This can be that due to TADs variability association of genes with brekapoints not easy
              geneBreakP_Position_respectToTSS<-"outOf_RegulatoryDomain"
              situation<-"primaryTAD_Central"
            }
            
            ###########################################
            ##+ Data required for Plot 
            ##Averages of the TAD maps taken as the reference values
            
            ###Info cognate enhancers
            nEnh_initial_left<-phaseIntegratedResults[targetGeneMech,"avg_nEnh_InitialLeft"]
            nEnh_initial_right<-phaseIntegratedResults[targetGeneMech,"avg_nEnh_InitialRight"]
            
            nEnh_kept_left<-phaseIntegratedResults[targetGeneMech,"avg_nEnh_KeptLeft"]
            nEnh_kept_right<-phaseIntegratedResults[targetGeneMech,"avg_nEnh_KeptRight"]
            
            ##Info ectopic enhancers
            if(pathoMechanism == "LOF"){
              ##for LOF we do not consider them, only cognate enh, hence set to 0
              ##Just setting to 0
              nEnh_other_domain<-0
              nEnh_gained<-0
              
            }else if(pathoMechanism == "GOF"){
              ##For GOF we of course consider ectopic enh, such as for enh adoption
              nEnh_other_domain<-phaseIntegratedResults[targetGeneMech,"avg_nEnh_InTheOtherDomain"]
              nEnh_gained<-phaseIntegratedResults[targetGeneMech,"avg_nEnh_GainedFromTheOtherDomain"] 
            }
            
            ##The nEnk kept...gained, will be NA on some cases as gene deletions
            ##Since we do not care about enh, so for this cases, if we have NA...assign the value to 0
            if(is.na(nEnh_kept_left)==TRUE){
              nEnh_kept_left<-0
            }
            
            if(is.na(nEnh_kept_right)==TRUE){
              nEnh_kept_right<-0
            }
            
            if(is.na(nEnh_gained)==TRUE){
              nEnh_gained<-0
            }
            
            
            ####################################################
            ## Wee need to define the type of line to be painted with respect to the enhancers and the gene
            ## Is the gene broken? Captures All the enhancers?
            
            ## OPTIONS FOR WHERE TO LOCATE THE BREAKPOINT LINE ON THE GENE TAD
            ## Options for gene_breakp_line_type: 
            ## If the gene is broken:
            ##        gene_breakp_line_type<-"center"
            ##
            ## If the breakpoint is after the TSS (bigger position):
            ## We know this if geneBreakP_Position_respectToTSS<-"afterTSS"
            ## So we need to look at the number of enhancers kept on the right:
            ##        gene_breakp_line_type<-"afterTSS_removing_all"; If NO enh Kept
            ##        gene_breakp_line_type<-"afterTSS_removing_some"; If some enh kept
            ##        gene_breakp_line_type<-"afterTSS_removing_none"; If ALL enh kept
            ##
            ## If the breakpoint is before the TSS (smaller position):
            ## We know this if geneBreakP_Position_respectToTSS<-"beforeTSS"
            ## So we need to look at the number of enhancers kept on the left:
            ##        gene_breakp_line_type<-"beforeTSS_removing_all"; If NO enh Kept
            ##        gene_breakp_line_type<-"beforeTSS_removing_some"; If some enh kept
            ##        gene_breakp_line_type<-"beforeTSS_removing_none"; If ALL enh kept
            ##
            ## OPTIONS FOR WHERE TO LOCATE THE BREAKPOINT LINE ON THE OTHER AFFECTED DOMAIN
            ## So we need to look at nEnhGained (does not matter by which side, we relocate them on the broken geneTad side)
            ##  otherDomain_breakp_line_type<-"brings_all"; If all enh from the other TAD put into proximity of the gene
            ##  otherDomain_breakp_line_type<-"brings_some"; If some enh from the other TAD put into proximity of the gene
            ##  otherDomain_breakp_line_type<-"brings_none"; If none enh from the other TAD put into proximity of the gene
            
            ###########
            ##Defining the breakpoints line location: gene_breakp_line_type (to locate the breakpoints, on both the gene TAD and the other Domain)
            
            ########################
            ##Breakpoint Gene TAD
            ##First to check, is the gene sequence directly affected?
            ##Try like this putting here gene Deletion too
            ##Distinguish with Duplication, because for duplication kept number can be bigger than the initial
            #Because cognate enh duplicated (and gained, only considered referring to new enh from other domains)
            
            ###For Duplications missing explanation of the logic followed but similar to previous described. But with particularities of duplications
            
            ##For long range, exclude from first if duplications due to their particularities
            ##They are handled on their own if
            if((gene_mechanism == "LongRange") && (sv_type != "Duplication")){
              
              ##Check in which situation of the described above we are
              if(geneBreakP_Position_respectToTSS=="afterTSS"){
                ##It implies the breakpoint is after the TSS (bigger position):
                ##Check n enh kept on the right
                if(nEnh_kept_right==0){
                  ##for now, does not matter wheter it had initially or not
                  gene_breakp_line_type<-"afterTSS_removing_all"
                  
                }else if(nEnh_initial_right>nEnh_kept_right){
                  gene_breakp_line_type<-"afterTSS_removing_some"
                  
                }else if(nEnh_initial_right==nEnh_kept_right){
                  gene_breakp_line_type<-"afterTSS_removing_none"
                }
                
              }else if(geneBreakP_Position_respectToTSS=="beforeTSS"){
                ##It implies the breakpoint is before the TSS (smaller position):
                ##Check nEnh kept on the left
                if(nEnh_kept_left==0){
                  ##for now, does not matter wheter it had initially or not
                  gene_breakp_line_type<-"beforeTSS_removing_all"
                  
                }else if(nEnh_initial_left>nEnh_kept_left){
                  gene_breakp_line_type<-"beforeTSS_removing_some"
                  
                }else if(nEnh_initial_left==nEnh_kept_left){
                  gene_breakp_line_type<-"beforeTSS_removing_none"
                }
                
              }
              
            }else if(gene_mechanism == "Direct_geneTruncation"){
              ###############################
              # For DIRECT GENE TRUNCATIONS
              ###############################
              
              gene_breakp_line_type<-"center"
              
              
              ###For Now:
              ## Set n enhnacers initial and in the other TAD to 0
              ## If gene pathogenic effect by Direct Gene Truncation we do not care about enhancers
              ##So not painting them in the plot (achieved by enh = 0)
              ##But if in future we predict fusion transcripts... should be painted if that is the prediction
              nEnh_initial_left<-0
              nEnh_initial_right<-0
              nEnh_other_domain<-0
              
              nEnh_kept_left<-0
              nEnh_kept_right<-0
              nEnh_gained<-0
              
              
            }else if(gene_mechanism == "Direct_geneDeletion"){
              ############################
              # For DIRECT GENE DELETIONS
              ############################
              
              ##Here we need to distinguish between InterTAD situation (two TADs painted)
              ##Vs IntraTAD... just one TAD painted
              
              if(SV_landing == "IntraTAD"){
                ##So geneDeletion IntraTAD... both break must be surrounding him
                ##I think we may need to use this again for Duplications
                gene_breakp_line_type<-"surroundingGene"
                
              }else if(SV_landing == "InterTAD"){
                ##So multiple TADs affected
                gene_breakp_line_type<- geneBreakP_Position_respectToTSS ##To know where to locate the breakpoint, that is all
                
                ##We paint the line, simply surrounding the gene
                ##Check in which situation of the described above we are
                if(geneBreakP_Position_respectToTSS == "afterTSS"){
                  ##It implies the breakpoint is after the TSS (bigger position):
                  ##for now, does not matter wheter it had initially or not
                  gene_breakp_line_type<-"afterTSS_removing_all"
                  
                }else if(geneBreakP_Position_respectToTSS == "beforeTSS"){
                  ##It implies the breakpoint is before the TSS (smaller position):
                  ##for now, does not matter wheter it had initially or not
                  gene_breakp_line_type<-"beforeTSS_removing_all"
                  
                }else if(geneBreakP_Position_respectToTSS == "outOf_RegulatoryDomain"){
                  ##the breakpoints painted outside of the regulatory domain
                  gene_breakp_line_type<-"outOf_RegulatoryDomain"
                }
                
                ## Set n enhnacers initial and in the other TAD to 0
                ## If gene deleted we do not care about enhancers
                ##So not painting them in the plot (achieved by enh = 0)
                nEnh_initial_left<-0
                nEnh_initial_right<-0
                nEnh_other_domain<-0
                
                nEnh_kept_left<-0
                nEnh_kept_right<-0
                nEnh_gained<-0 
              }
              
              ## Set n enhnacers initial and in the other TAD to 0
              ## If gene pathogenic effect by Direct Gene Deletion we do not care about enhancers
              ##So not painting them in the plot (achieved by enh = 0)
              nEnh_initial_left<-0
              nEnh_initial_right<-0
              nEnh_other_domain<-0
              
              nEnh_kept_left<-0
              nEnh_kept_right<-0
              nEnh_gained<-0
              
              
            }else if(sv_type == "Duplication"){
              ################################################################################
              ##Duplications have their own kind of particularities so let's treat them a part
              ################################################################################
              
              ##Focusing on duplic that enh can be triggering the disease
              if((gene_mechanism == "LongRange") || ##In this case the gene is not duplicated
                 (gene_mechanism == "LongRange_geneDuplication") ||
                 (gene_mechanism =="Direct_LongRange_geneDuplication")){
                
                ##################################################################
                ## IMPORTANT PARTICULARITIE. "Kept" tracks cognate enhancers 
                ## So if nkept > intial is because some of the initial duplicated
                
                if(geneBreakP_Position_respectToTSS=="afterTSS"){
                  ##It implies the breakpoint is after the TSS (bigger position):
                  ##Check n enh kept on the right
                  if(nEnh_kept_right == nEnh_initial_right){
                    ##for now, does not matter wheter it had initially or not
                    gene_breakp_line_type<-"afterTSS_duplic_none"##we will put the line between gene and enh
                    
                  }else if(nEnh_kept_right == 2*nEnh_initial_right ){
                    gene_breakp_line_type<-"afterTSS_duplic_all" 
                    
                  }else if(nEnh_kept_right > nEnh_initial_right){
                    gene_breakp_line_type<-"afterTSS_duplic_some" ##we will put the line between enh
                  }
                  
                }else if(geneBreakP_Position_respectToTSS=="beforeTSS"){
                  ##It implies the breakpoint is before the TSS (smaller position):
                  ##Check nEnh kept on the left
                  if(nEnh_kept_left == nEnh_initial_left ){
                    ##for now, does not matter wheter it had initially or not
                    gene_breakp_line_type<-"beforeTSS_duplic_none"##we will put the line between gene and enh
                    
                  }else if(nEnh_kept_left == 2*nEnh_initial_left){
                    gene_breakp_line_type<-"beforeTSS_duplic_all"##all enhancer on that side duplicated
                    
                  }else if(nEnh_kept_left > nEnh_initial_left){
                    gene_breakp_line_type<-"beforeTSS_duplic_some" ##we will put the line between enh
                  }
                  
                }
                
              } else if(gene_mechanism == "Direct_geneDuplication"){
                ############################
                # For DIRECT GENE DUPLICATIONS
                ############################
                
                ##Here we need to distinguish between InterTAD situation (two TADs painted)
                ##Vs IntraTAD... just one TAD painted
                
                if(SV_landing == "IntraTAD"){
                  ##So geneDuplication, pathogenic IntraTAD... both break must be surrounding him
                  gene_breakp_line_type<-"surroundingGene"
                  
                }else if(SV_landing == "InterTAD"){
                  ##So multiple TADs affected
                  gene_breakp_line_type<- geneBreakP_Position_respectToTSS ##To know where to locate the breakpoint, that is all
                  
                  ##We paint the line, simply surrounding the gene
                  ##Check in which situation of the described above we are
                  if(geneBreakP_Position_respectToTSS == "afterTSS"){
                    ##It implies the breakpoint is after the TSS (bigger position):
                    ##for now, does not matter wheter it had initially or not
                    gene_breakp_line_type<-"afterTSS_duplic_none"
                    
                  }else if(geneBreakP_Position_respectToTSS == "beforeTSS"){
                    ##It implies the breakpoint is before the TSS (smaller position):
                    ##for now, does not matter wheter it had initially or not
                    gene_breakp_line_type<-"beforeTSS_duplic_none"
                    
                  }else if(geneBreakP_Position_respectToTSS == "outOf_RegulatoryDomain"){
                    ##the breakpoints painted outside of the regulatory domain
                    gene_breakp_line_type<-"outOf_RegulatoryDomain"
                  }
                }
                
                ## Set n enhnacers initial and in the other TAD to 0
                ## If gene pathogenic effect by Direct Gene Duplication we do not care about enhancers
                ##So not painting them in the plot (achieved by enh = 0)
                nEnh_initial_left<-0
                nEnh_initial_right<-0
                nEnh_other_domain<-0
                
                nEnh_kept_left<-0
                nEnh_kept_right<-0
                nEnh_gained<-0 
                
              }
              
            }
            
            
            ########################################
            ##  The other Domain Breakpoint Location 
            ########################################
            ##Only used when >1 domain painted
            if(nEnh_gained==0){
              ##For now does not matter whether there were initially or not
              otherDomain_breakp_line_type<-"brings_none" 
            }else if(nEnh_gained==nEnh_other_domain){
              ##So all the enhancers from the other domain have been won
              ##AND the  nEnh gained is different than 0
              otherDomain_breakp_line_type<-"brings_all"
            }else if(nEnh_gained<nEnh_other_domain){
              otherDomain_breakp_line_type<-"brings_some"
            }
            ##Depending on primary tad sinistral or dextral, the none and the all will be place at one side or the other of the TAD
            
            ####################
            ##STARTING WITH PLOT
            ####################
            
            ##We are going to do plots & reports for thins with "high" or "pathogenic scores"
            ##But also for direct impacts because for sure people will want to take a look at the directly affected genes
            ##So track it on the list of Not, main candidates
            
            if(gene_breakpoint != "NotClear_BreakpointAssociation"){
              ##So we generate graphical Abstract 
              ##We will see what to do with NotClear_BreakpointAssociation, but for now skip to avoid crash
              ###########
              #Image Output Path
              
              outpPath<-"www/graphicalSummaries/"
              ##outpPath<-"graphicalSummaries/"
              fullOutpPath<-paste0(outpPath, gene,"_",targetMech, "_", phase,"_", patientResults$job_UniCode, ".png")
              
              #############################
              ## Loading auxiliar functions for plotting
              ##source("functions/aux_functions_forGraphicalSummary.R")
              
              ####################################
              ##png with maximum resolution 300dpi
              png(filename = fullOutpPath, width = 12, height = 8, units = "in", res = 300 )
              ##x11()
              ##########################################
              ###Creating canvas for plotting
              yAxisLim<-c(-50,30)
              
              ##tagEnhancersLabel to avoid enhancershancers when neoTAD painting, to avoid overlap
              
              tagEnhancersLabel<-"enhancers"
              
              ##Defining xAxis Width
              if(((gene_mechanism == "LongRange_geneDuplication") || (gene_mechanism == "Direct_LongRange_geneDuplication") ) && (SV_landing == "InterTAD")){
                ##xAxiss wider if NEO TAD because a third tad will be painted
                ##Usar un lienzo mas grande, due to NEO-TAD situation
                xAxisLim<-c(-15,55)
                tagEnhancersLabel<-"Enh."
                
              }else{
                xAxisLim<-c(0,40)
              }
              
              
              yPos_chr_WT<-13 ##variable to hold the position of the chr text  in the WT situation in the Y axis
              
              ##Adjust image, to exclude spaces outside canvas (drawing area)
              par(mar = c(0,0,0,0))
              
              #When I want to se the axis, uncomment the 3 following lines
              # plot(x=1:2, y=1:2, type="n",
              #       ylim=yAxisLim,
              #       xlim=xAxisLim)
              
              ##OPTION REMOVING AXIS
              #UPON COMPLETION USING THIS
              plot(x=0:20, y=0:20, type = "n",
                   ylim = yAxisLim,
                   xlim = xAxisLim,
                   xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '') ##Upon completion, REMOVE AXIS
              
              ##############################################################################
              ##############################################################################
              
              texSizeChr<-2
              # enhNumbersSize<-1##Some variables for text sizes
              
              #Adding plot header
              text(x=20, y=30, label=paste0("Graphical summary of ", gene, " regulatory domain in ", phase), cex = 1.7)
              
              text(x=20,y=26, label="Simplification",cex=1.7)
              
              #WT allelle label
              text(x=20, y=20, label="Control Scenario", cex = 1.5)
              
              
              ##Check whether the current Gene gene regulatory domain is disrupted or not (it can occur for deletions or duplications)
              ##If the TAD gene is entirely deleted or duplicated, between the affected ones
              ##We will know that this is occurring if the gene_breakpoint takes the value<-"No_AssociatedBreakpoint"
              ##Still pendent to take into account if value: "NotClear_BreakpointAssociation"
              
              ##Again important to check whether gene regulatory domain altered or not
              #Info derived from gene breakpoint
              ##gene_regulatoryDomain_altered##TRUE OR FALSE
              ##situation "primaryTAD_Central" means we only paint one tad in the middle cause it is not disrupted
              
              ##Coord for tads over X axis (Depending on the situation either two TADs or just one will be painted)
              #Let's try make them wider two cm per side to ensure everything fits on the re-arrangement plots
              tad_XCoord_OnLeftSide<-c(3,17,10)#c(5,15,10) ##c(3,17,10)
              tad_XCoord_OnRightSide<-c(tad_XCoord_OnLeftSide[1]+10+11,
                                        tad_XCoord_OnLeftSide[2]+10+11,
                                        tad_XCoord_OnLeftSide[3]+10+11)
              tad_XCoord_OnCenter<-c(13,27,20)
              
              #Over Y axis, WT line
              tad_YCoord_WildTypeLine<-c(0,0,15)
              
              
              
              if(situation %in% c("primaryTAD_Sinistral", "primaryTAD_Dextral")){
                ##So we paint two TADs cause gene breakpoint disrupts its regulatory domain, and re-shuffles it with another
                
                if(situation =="primaryTAD_Sinistral"){
                  ########################
                  ## PAINTING GENE-TAD ###
                  ## on the left side 
                  ########################
                  
                  text(x=5, y=yPos_chr_WT, label=chr_gene, cex = texSizeChr)
                  
                  ##Creating TAD to represente gene WT scenario
                  
                  tad_X_cord<-tad_XCoord_OnLeftSide ##c(5,15,10)
                  info_drawingGENE_TAD<-paintGene_WT_TAD(tad_X_cord = tad_X_cord,
                                                         tad_Y_cord = tad_YCoord_WildTypeLine,
                                                         nEnh_initial_left = nEnh_initial_left, 
                                                         nEnh_initial_right = nEnh_initial_right,
                                                         gene = gene,
                                                         gene_breakp_line_type = gene_breakp_line_type,
                                                         situation = situation,
                                                         patientResults = patientResults, tagEnhancersLabel = tagEnhancersLabel)
                  
                  geneCenter<-info_drawingGENE_TAD$geneCenter
                  
                  # ##Adding interTAD dotted lines
                  # paintInterTAD_WT_lines(tad_X_cord = c(5,15,10))
                  
                  ################################################
                  ## PAINTING SECONDARY//OTHER AFFECTED DOMAIN ###
                  ## on the right side 
                  ################################################
                  text(x=26, y=yPos_chr_WT, label=chr_oppositeDomain, cex = texSizeChr)
                  ##Creating TAD to represente gene WT scenario
                  tad_X_cord<-tad_XCoord_OnRightSide
                  tad_Y_cord<-tad_YCoord_WildTypeLine
                  ##Info returned used to represent the rearrangement
                  
                  info_drawingSecondaryTAD<-paint_Enhancer_WT_Secondary_TAD(tad_X_cord = tad_X_cord,
                                                                            tad_Y_cord = tad_YCoord_WildTypeLine,
                                                                            nEnh_other_domain = nEnh_other_domain,
                                                                            geneCenter = geneCenter,
                                                                            otherDomain_breakp_line_type = otherDomain_breakp_line_type,
                                                                            situation = situation,
                                                                            geneBreakP_Position_respectToTSS = geneBreakP_Position_respectToTSS,
                                                                            patientResults = patientResults, tagEnhancersLabel = tagEnhancersLabel)
                  
                }else if(situation == "primaryTAD_Dextral"){
                  ########################
                  ## PAINTING GENE-TAD ###
                  ## on the Right Side side 
                  ########################
                  
                  text(x=26, y=yPos_chr_WT, label=chr_gene, cex = texSizeChr)
                  
                  ##Creating TAD to represente gene WT scenario
                  tad_X_cord<-tad_XCoord_OnRightSide##c(5+10+11,15+10+11,10+10+11)
                  ##Info returned used to represent the rearrangement
                  info_drawingGENE_TAD<-paintGene_WT_TAD(tad_X_cord = tad_X_cord,
                                                         tad_Y_cord = tad_YCoord_WildTypeLine,
                                                         nEnh_initial_left = nEnh_initial_left, 
                                                         nEnh_initial_right = nEnh_initial_right,
                                                         gene = gene,
                                                         gene_breakp_line_type = gene_breakp_line_type,
                                                         situation = situation,
                                                         patientResults = patientResults, tagEnhancersLabel = tagEnhancersLabel)
                  
                  geneCenter<-info_drawingGENE_TAD$geneCenter
                  
                  # ##Adding interTAD dotted lines
                  # paintInterTAD_WT_lines(tad_X_cord = c(5,15,10))
                  
                  ################################################
                  ## PAINTING SECONDARY//OTHER AFFECTED DOMAIN ###
                  ## on the LEFT side 
                  ################################################
                  text(x=5, y=yPos_chr_WT, label=chr_oppositeDomain, cex = texSizeChr)
                  ##Creating TAD to represente gene WT scenario
                  tad_X_cord<-tad_XCoord_OnLeftSide
                  tad_Y_cord<-tad_YCoord_WildTypeLine
                  
                  ##Info returned used to represent the rearrangement
                  
                  info_drawingSecondaryTAD<-paint_Enhancer_WT_Secondary_TAD(tad_X_cord = tad_X_cord,
                                                                            tad_Y_cord = tad_YCoord_WildTypeLine,
                                                                            nEnh_other_domain = nEnh_other_domain,
                                                                            geneCenter = geneCenter,
                                                                            otherDomain_breakp_line_type = otherDomain_breakp_line_type,
                                                                            situation = situation,
                                                                            geneBreakP_Position_respectToTSS = geneBreakP_Position_respectToTSS,
                                                                            patientResults = patientResults, tagEnhancersLabel = tagEnhancersLabel)
                }
                
              }else if(situation == "primaryTAD_Central"){
                ##So either the breakpoint falls in a different TAD, and this TAD is either duplicated or deleted entirely
                ##Or both breakpoint fall inside of the same TAD
                ##In any case, breakpoints do not affect its boundaries
                
                ##So paint it in the center
                ##With breakpoints flanking theTAD upon a certain distance
                ##source("functions/aux_functions_forGraphicalSummary.R")
                
                ########################
                ## PAINTING GENE-TAD ###
                ## on the Center (Only one TAD painted)
                ########################
                
                text(x=tad_XCoord_OnCenter[1]+2, y=yPos_chr_WT, label=chr_gene, cex = texSizeChr)
                
                ##Creating TAD to represente gene WT scenario
                
                tad_X_cord<-tad_XCoord_OnCenter
                
                info_drawingGENE_TAD<-paintGene_WT_TAD(tad_X_cord = tad_X_cord,
                                                       tad_Y_cord = tad_YCoord_WildTypeLine,
                                                       nEnh_initial_left = nEnh_initial_left, 
                                                       nEnh_initial_right = nEnh_initial_right,
                                                       gene = gene,
                                                       gene_breakp_line_type = gene_breakp_line_type,
                                                       situation = situation,
                                                       patientResults = patientResults, tagEnhancersLabel = tagEnhancersLabel)
                
                geneCenter<-info_drawingGENE_TAD$geneCenter
                
              }
              
              
              ###########################
              ## Adding Label SV type
              ###########################
              ##source("functions/aux_functions_forGraphicalSummary.R")
              paint_SV_labels(sv_type = sv_type, xAxisLim = xAxisLim, yAxisPos = -12.5)
              
              ##############################
              ###Paint interTad arrows
              ##############################
              
              ##Avoid if is translocation
              ##For now remove, too noisy
              # if(sv_type=="Inversion"){
              #   paintInterTAD_arrows() 
              # }
              
              ############################################################
              ## Painting Patient Re-arranged Scenario ###################
              ## Here is when plot changes, regarding SV
              ############################################################
              ##source("functions/makingGraphs/aux_fun_Plots_Rearranged_TADs.R")
              
              #Adding SV rearrangment header
              #WT allelle label
              text(x=20, y=20-40, label="Patient Scenario", cex = 1.5)
              #text(x=20,y=18-30, label="Simplification",cex=1.2)
              #text(x=20, y=16-30, label=paste0("Focusing on gene ", gene, " in phase ", phase), cex = 1.2)
              
              #Height at which TADs are displayed
              tad_YCoord_Rearrangements<-tad_YCoord_WildTypeLine-40
              
              #############################################################################################################
              ##Here, we need to differentiate depending on the SV
              ##Also if gene direct effect, for LOF (as truncation or deletion), do not paint the whole rearrangement
              #############################################################################################################
              
              if(gene_mechanism == "Direct_geneTruncation"){
                ##Plot GeneTruncation. Only one plot in the middle no TAD painted
                paint_Gene_Truncation(gene = gene,
                                      xAxisLim = xAxisLim,
                                      tad_X_cord = tad_XCoord_OnCenter,
                                      tad_YCoord_Rearrangements = tad_YCoord_Rearrangements)
                
              }else if(gene_mechanism == "Direct_geneDeletion"){
                ##Plot GeneDeletion. Only one plot in the middle no TAD painted
                paint_Gene_Deletion(gene = gene,
                                    xAxisLim = xAxisLim,
                                    tad_X_cord = tad_XCoord_OnCenter,
                                    tad_YCoord_Rearrangements = tad_YCoord_Rearrangements)
                
              }else if(gene_mechanism == "LongRange"){
                
                if((sv_type=="Inversion") || (sv_type=="Translocation")){
                  ##It implies a reshuffling of TADs, so no DNA gained or lost
                  
                  #Paint Gene SV re-arranged TAD
                  infoDrawing_gene_sv_tad<-paintGene_SV_TAD(nEnh_initial_left = nEnh_initial_left, 
                                                            nEnh_initial_right = nEnh_initial_right,
                                                            gene = gene,
                                                            gene_breakp_line_type = gene_breakp_line_type,
                                                            situation = situation,
                                                            nEnh_kept_left = nEnh_kept_left, 
                                                            nEnh_kept_right = nEnh_kept_right,
                                                            nEnh_gained = nEnh_gained,
                                                            nEnh_other_domain = nEnh_other_domain,
                                                            info_drawingGENE_TAD = info_drawingGENE_TAD, 
                                                            info_drawingSecondaryTAD = info_drawingSecondaryTAD,
                                                            tad_XCoord_OnLeftSide = tad_XCoord_OnLeftSide, 
                                                            tad_XCoord_OnRightSide = tad_XCoord_OnRightSide,
                                                            tad_YCoord_Rearrangements = tad_YCoord_Rearrangements,
                                                            geneBreakP_Position_respectToTSS = geneBreakP_Position_respectToTSS, tagEnhancersLabel = tagEnhancersLabel)
                  
                  #Paint SecondaryTAD SV re-arranged
                  #The one where the gene of interest is not located
                  paintSecondaryDomain_SV_TAD(nEnh_initial_left = nEnh_initial_left, 
                                              nEnh_initial_right = nEnh_initial_right,
                                              gene = gene,
                                              gene_breakp_line_type = gene_breakp_line_type,
                                              situation = situation,
                                              nEnh_other_domain = nEnh_other_domain,
                                              nEnh_kept_left = nEnh_kept_left, 
                                              nEnh_kept_right = nEnh_kept_right,
                                              nEnh_gained = nEnh_gained,
                                              info_drawingGENE_TAD = info_drawingGENE_TAD, 
                                              info_drawingSecondaryTAD = info_drawingSecondaryTAD,
                                              tad_XCoord_OnLeftSide = tad_XCoord_OnLeftSide, 
                                              tad_XCoord_OnRightSide = tad_XCoord_OnRightSide,
                                              tad_YCoord_Rearrangements = tad_YCoord_Rearrangements,
                                              infoDrawing_gene_sv_tad = infoDrawing_gene_sv_tad,  tagEnhancersLabel = tagEnhancersLabel)
                  
                }else if(sv_type=="Deletion"){
                  ##For Deletions, by long-range... so TAD disrupted... in any case intraTAD or betweenTADs only one TAD painted
                  paintGene_SV_Deletion_TAD(nEnh_initial_left = nEnh_initial_left, 
                                            nEnh_initial_right = nEnh_initial_right,
                                            gene = gene,
                                            gene_breakp_line_type = gene_breakp_line_type,
                                            situation = situation,
                                            nEnh_kept_left = nEnh_kept_left, 
                                            nEnh_kept_right = nEnh_kept_right,
                                            nEnh_gained = nEnh_gained,
                                            nEnh_other_domain = nEnh_other_domain,
                                            info_drawingGENE_TAD = info_drawingGENE_TAD, 
                                            info_drawingSecondaryTAD = info_drawingSecondaryTAD,
                                            tad_X_cord = tad_XCoord_OnCenter ,
                                            tad_YCoord_Rearrangements = tad_YCoord_Rearrangements,
                                            geneBreakP_Position_respectToTSS = geneBreakP_Position_respectToTSS, tagEnhancersLabel = tagEnhancersLabel)
                  
                }else if(sv_type == "Duplication"){
                  
                  ##In this context there is no gene duplication
                  ###If we are facing a LongRange mechanism, and the gene is not duplicated, as it is happening in this case
                  ###It must be intraTAD, because on the contrary, if a gene Not duplicated, and in a InterTAD SV, the gene regulatory 
                  ##Domain remains as it was.
                  
                  ##If it is not intra-TAD, raise error, it can be an interTAD where breakpoint falls in enh...hence counted as lost
                  ##So LOF happening by means of a Duplication... but rare. So if this is the case rise error to check what is going on
                  
                  
                  ##Check if SV is intraTAD (primaryTAD_central or phaseIntegratedResults$SV_landing IntraTAD)
                  
                  ##To be a NeoTAD process the gene needs to be duplicated
                  ##Check that THIS is INTRA-TAD  
                  
                  if(SV_landing == "IntraTAD"){
                    
                    #okey pintar
                    ##Pintar el TAD con los enh, si hay mas que habian en un lado pintar 6 u 8 si se doblan todos
                    ##Si hay menos... tal vez por duplication falling strictly in a enh... pintar menos.. Buff for now not considered
                    paintGene_SV_Duplication_OnlyLongRange_Intra_TAD(nEnh_initial_left = nEnh_initial_left, 
                                                                     nEnh_initial_right = nEnh_initial_right,
                                                                     gene = gene,
                                                                     gene_breakp_line_type = gene_breakp_line_type,
                                                                     situation = situation,
                                                                     nEnh_kept_left = nEnh_kept_left, 
                                                                     nEnh_kept_right = nEnh_kept_right,
                                                                     nEnh_gained = nEnh_gained,
                                                                     nEnh_other_domain = nEnh_other_domain,
                                                                     info_drawingGENE_TAD = info_drawingGENE_TAD, 
                                                                     info_drawingSecondaryTAD = info_drawingSecondaryTAD,
                                                                     tad_X_cord = tad_XCoord_OnCenter ,
                                                                     tad_YCoord_Rearrangements = tad_YCoord_Rearrangements,
                                                                     geneBreakP_Position_respectToTSS = geneBreakP_Position_respectToTSS,
                                                                     patientResults = patientResults, tagEnhancersLabel = tagEnhancersLabel)
                    
                  }else{
                    ##Error, no concibo esta situacion ahora mismo
                    stop("I don't contemplate how this can happen sensibly right now")
                  }
                  
                  
                  
                }
              }else if((gene_mechanism == "LongRange_geneDuplication") || (gene_mechanism == "Direct_LongRange_geneDuplication")){
                ##Here it either can occur, NeoTAD if SV interTAD or be intraTAD, but for sure gene Duplicated
                ##So two genes need to be painted
                
                ## NEO-TAD processes. If interTAD
                
                ##Filter for interTAD, if it is intraTAD i do not reasonably understand how this can be pathogenic
                ##IF not considered as Direct or Direct_LonRange_geneDuplication
                ##Because this would imply that gene poorly expressed, by cognate enh duplication pathogenic
                ##And we are trying to prevent this from being a good candidate
                
                
                
                if(SV_landing == "InterTAD" ){
                  
                  ######################
                  ## NEO-TAD process
                  ######################
                  ##Left and Right TADs are painted same as WT at the sides
                  ##So reuse the same script styling
                  ##Situation can not be primary tad central
                  
                  ######################
                  ## Painting borderTADs
                  tad_Y_cord<-tad_YCoord_Rearrangements
                  
                  tad_XCoord_MoreLefted_OnLeftSide<-tad_XCoord_OnLeftSide - 10
                  tad_XCoord_MoreRighted_OnRightSide<-tad_XCoord_OnRightSide + 10
                  
                  
                  ##Modificar TAD X coord tirarlos mas a  los lados
                  
                  if(situation =="primaryTAD_Sinistral"){
                    ########################
                    ## PAINTING GENE-TAD ###
                    ## on the left side 
                    ########################
                    
                    #text(x=5, y=yPos_chr_WT, label=chr_gene, cex = texSizeChr)
                    
                    ##Creating TAD to represente gene WT scenario
                    
                    tad_X_cord<-tad_XCoord_MoreLefted_OnLeftSide
                    info_drawingGENE_TAD<-paintGene_WT_TAD(tad_X_cord = tad_X_cord,
                                                           tad_Y_cord = tad_YCoord_Rearrangements,
                                                           nEnh_initial_left = nEnh_initial_left, 
                                                           nEnh_initial_right = nEnh_initial_right,
                                                           gene = gene,
                                                           gene_breakp_line_type = gene_breakp_line_type,
                                                           situation = situation,
                                                           patientResults = patientResults, tagEnhancersLabel = tagEnhancersLabel)
                    
                    geneCenter<-info_drawingGENE_TAD$geneCenter
                    
                    
                    ################################################
                    ## PAINTING SECONDARY//OTHER AFFECTED DOMAIN ###
                    ## on the right side 
                    ################################################
                    #text(x=26, y=yPos_chr_WT, label=chr_oppositeDomain, cex = texSizeChr)
                    ##Creating TAD to represente gene WT scenario
                    tad_X_cord<-tad_XCoord_MoreRighted_OnRightSide
                    
                    ##Info returned used to represent the rearrangement
                    
                    info_drawingSecondaryTAD<-paint_Enhancer_WT_Secondary_TAD(tad_X_cord = tad_X_cord,
                                                                              tad_Y_cord = tad_YCoord_Rearrangements,
                                                                              nEnh_other_domain = nEnh_other_domain,
                                                                              geneCenter = geneCenter,
                                                                              otherDomain_breakp_line_type = otherDomain_breakp_line_type,
                                                                              situation = situation,
                                                                              geneBreakP_Position_respectToTSS = geneBreakP_Position_respectToTSS,
                                                                              patientResults = patientResults, tagEnhancersLabel = tagEnhancersLabel)
                    
                  }else if(situation == "primaryTAD_Dextral"){
                    ########################
                    ## PAINTING GENE-TAD ###
                    ## on the Right Side side 
                    ########################
                    
                    #text(x=26, y=yPos_chr_WT, label=chr_gene, cex = texSizeChr)
                    
                    ##Creating TAD to represente gene WT scenario
                    tad_X_cord<-tad_XCoord_MoreRighted_OnRightSide
                    
                    ##Info returned used to represent the rearrangement
                    info_drawingGENE_TAD<-paintGene_WT_TAD(tad_X_cord = tad_X_cord,
                                                           tad_Y_cord = tad_YCoord_Rearrangements,
                                                           nEnh_initial_left = nEnh_initial_left, 
                                                           nEnh_initial_right = nEnh_initial_right,
                                                           gene = gene,
                                                           gene_breakp_line_type = gene_breakp_line_type,
                                                           situation = situation,
                                                           patientResults = patientResults, tagEnhancersLabel = tagEnhancersLabel)
                    
                    geneCenter<-info_drawingGENE_TAD$geneCenter
                    
                    
                    ################################################
                    ## PAINTING SECONDARY//OTHER AFFECTED DOMAIN ###
                    ## on the LEFT side 
                    ################################################
                    #text(x=5, y=yPos_chr_WT, label=chr_oppositeDomain, cex = texSizeChr)
                    ##Creating TAD to represente gene WT scenario
                    tad_X_cord<-tad_XCoord_MoreLefted_OnLeftSide
                    
                    
                    ##Info returned used to represent the rearrangement
                    
                    info_drawingSecondaryTAD<-paint_Enhancer_WT_Secondary_TAD(tad_X_cord = tad_X_cord,
                                                                              tad_Y_cord = tad_YCoord_Rearrangements,
                                                                              nEnh_other_domain = nEnh_other_domain,
                                                                              geneCenter = geneCenter,
                                                                              otherDomain_breakp_line_type = otherDomain_breakp_line_type,
                                                                              situation = situation,
                                                                              geneBreakP_Position_respectToTSS = geneBreakP_Position_respectToTSS,
                                                                              patientResults = patientResults, tagEnhancersLabel = tagEnhancersLabel)
                  }
                  
                  ########################
                  ##Painting the NEO TAD 
                  ########################
                  paintGene_SV_Duplication_NEO_TAD(nEnh_initial_left = nEnh_initial_left, 
                                                   nEnh_initial_right = nEnh_initial_right,
                                                   gene = gene,
                                                   gene_breakp_line_type = gene_breakp_line_type,
                                                   situation = situation,
                                                   nEnh_kept_left = nEnh_kept_left, 
                                                   nEnh_kept_right = nEnh_kept_right,
                                                   nEnh_gained = nEnh_gained,
                                                   nEnh_other_domain = nEnh_other_domain,
                                                   info_drawingGENE_TAD = info_drawingGENE_TAD, 
                                                   info_drawingSecondaryTAD = info_drawingSecondaryTAD,
                                                   tad_X_cord = tad_XCoord_OnCenter ,
                                                   tad_YCoord_Rearrangements = tad_YCoord_Rearrangements,
                                                   geneBreakP_Position_respectToTSS = geneBreakP_Position_respectToTSS, tagEnhancersLabel = tagEnhancersLabel)
                  
                  
                  
                  
                }else {
                  #It would be intraTAD dup, where pathogenic mech by enh dup and not by gene dup
                  #...rare, attending to what I have stated before
                  ##Error, no concibo esta situacion ahora mismo
                  stop("I don't contemplate how this can happen sensibly right now")
                }
                
                
                
              }else if(gene_mechanism == "Direct_geneDuplication"){
                ##We paint two times only the gene, and enough...
                
                paint_Gene_Duplicated(gene = gene,
                                      xAxisLim = xAxisLim,
                                      tad_X_cord = tad_XCoord_OnCenter,
                                      tad_YCoord_Rearrangements = tad_YCoord_Rearrangements)
                
                
              }
              
              #######################################
              ## CLOSING && STORING PLOT ############
              #######################################
              
              dev.off()##Saving Graph
            }
          }
        }
      }
    }
  }
  
  return(patientResults)##Contains info about genes filtered by minScore, relevant to do report
}
