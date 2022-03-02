#################################################
## Function to evaluate the genes situation
## Regarding the different considered variables
## Report it as a matrix

###Function to calculate the nr of enh kept and gained
source(file = "functions/CalculatingNrEnhBeforeAndAfter.R",local = TRUE)


##We can add more variables if we want in the future
##################################################
evaluatingGenesSituation<-function(genesInfo,domainsInfo,enhancersInfo,dataPatient, humanPhenotype, micePhenotype,
                                   polyCombGenes,##genes with H3K27me3 
                                   enhMap,
                                   tau_exp_scores,
                                   poorExp,##poorly expressed genes ids,
                                   natureLof,##natureLofScores //HI score for genes
                                   huangLof,##huang Hi scores
                                   clinGen_hiInfo,##clinGene HI info
                                   geneExp,##genesExpresion
                                   mainPatientPhenotype,##Patient considered phenotype to carry on analysis
                                   phasesVector)##vector holding the name of the different considered and studied phases 
  {
  typeSV<-dataPatient$TypeSV
  
  source("scripts_To_Load_Data/ExpressionThresholds_ForGofAndLof.R",
         local = TRUE)
  
  PoorlyExpThresh<-threshold_MinExpresion
  ConsiderablyExpThresh<-threshold_MaxExpresion 
  
  ######################################################
  ##Preparing Result Matrix
  ####la result matrix tiene una fila por gene afectado
  ##Y una columna por variable a estudio
  affectedGenes<-unique(genesInfo$genesPosition$geneSymbol)
  brokenGenes<-unique(genesInfo$brokenGenes)
  genesInUncertaintyRegions<-unique(genesInfo$genesInUncertainty)
  #In case there are
  deletedGenes<-unique(genesInfo$deletedGenes)
  duplicatedGenes<-unique(genesInfo$duplicatedGenes)
  
  ##For now we only compute enh statistics for genes that are in a TAD which is partially disrupted
  ## For the genes that are located in a TAD entirely deleted or duplicated we do not want enhancer info
  ## If it is expressed... expression boost, if it is not... duplicating the TAD not apparent consecuences
  genesInDisruptedTAD<-unique(c(genesInfo$genesInDomain$domain_breakpoint_1, genesInfo$genesInDomain$domain_breakpoint_2))
  
  ##Studied Features to appear in evaluation matrices
  studiedFeatures<-c("associatedPhenotypeIn_OMIM",
                     "associatedPhenotypeIn_MGI",
                     # "n_PhenotypeRelated_Through_OMIM_Human",##integer (nPhenotypes Matching)
                     # "n_PhenotypeRelated_Through_MGI_Mice",##integer (nPhenotypes Matching)
                     
                     "mainPhenotype_Through_OMIM_Human",##TRUE,FALSE #MainPhenotype is the one selected by the user (the one used to load enh, exp and TAD data)
                     "mainPhenotype_Through_MGI_Mice",##TRUE,FALSE
                     
                    
                     #####################################################
                     ##Gain of Function, Loss of Function candidates
                     ##Nature LOF scores "the probability of being loss-of-function intolerant (intolerant of both heterozygous and homozygous lof variants)"
                     "nature_HI_score",
                     "huang_HI_score",
                     "clinGene_HI_score",
                     
                     ##TAU_exp score
                     ##For consistently expressed genes
                     ##If not consistently the TAU==-1
                     ##if expressed at least 5fpkm in any early stage tissue, TAU (to avoid 0,0,1)
                     ##So if number dif -1 because consistently expressed in at least 1 condition
                     ##Then if close to 0 housekeeping, && if close to one tissue specific
                     ##inside of this intra-fetal system tissues (tissue specific at developmental stages would be)
                     "TAU_exp",
                     
                     ###Polycomb promoter
                     "polyComb_score",
                     
                     "RegulatoryMechanism",## Direct_geneTruncation, Direct_geneDeletion, direct_geneDuplication(for genes active), longRange_GeneDuplication... This column updated when predictions done | indirect_geneIntact indirect_geneDuplication
                     "TypeDomainInitial",##TAD or between_TAD
                     "TypeDomainMixedWith",##TAD or between_TAD (if it is something in between, not defined if merging a TAD with an interTAD)

                     "InitialDomain_ID",##To track if is the same or not for the different breakpoints. The id of the inital domain
                     
                     ##track enhancers of each type
                     ##GAINED is restricted to enh coming from gene-different Regulatory Domain (so in duplication (only tricky case) enh duplicated in the same tad not considered as gained but as "double" kept)
                     paste0("nEnhancers_initial_",phasesVector),
                     paste0("nEnhancers_kept_",phasesVector),
                     paste0("nEnhancers_gained_",phasesVector),
                     paste0("nEnhancers_maxAvailableInTheOtherDomain_",phasesVector),
                     
                     ###To retrieve fpkm for the different stages
                     ##AQUI MEJOR FPKM Y PHASES
                     paste0("FPKM_",phasesVector),
                     
                     ####To track the initial acetilation levels on the enh, how much acetilation is lost and gained on the SV
                     paste0("enhancers_acetilation_initial_",phasesVector),
                     paste0("enhancers_acetilation_kept_",phasesVector),
                     paste0("enhancers_acetilation_gained_",phasesVector),
                     paste0("enhancers_maxAcetilationAvailableInTheOtherDomain_",phasesVector),
                     
                     
                     ##To track the number of enhancers toward the upper and lower TAD border with respect to each gene
                     ##Useful for graphical displays
                     paste0("nEnh_ToTheLeft_initial_",phasesVector),
                     paste0("nEnh_ToTheRight_initial_",phasesVector),
                     paste0("nEnh_ToTheLeft_kept_",phasesVector),
                     paste0("nEnh_ToTheRight_kept_",phasesVector),
                     # paste0("nEnh_ToTheLeft_gained_",phasesVector),
                     # paste0("nEnh_ToTheRight_gained_",phasesVector),
                     
                     ##To track the number of enhancers toward the upper and lower TAD border with respect to each gene
                     ##Useful for graphical displays
                     paste0("enhancers_acetilation_ToTheLeft_initial_",phasesVector),
                     paste0("enhancers_acetilation_ToTheRight_initial_",phasesVector),
                     paste0("enhancers_acetilation_ToTheLeft_kept_",phasesVector),
                     paste0("enhancers_acetilation_ToTheRight_kept_",phasesVector),
                     
                     ###Gene breakpoint
                     "gene_Breakpoint"##To track whether gene is associated with breakpoint 1 or 2
                     )
  
  resultMatrix<-as.data.frame(matrix(data = NA, nrow =length(affectedGenes), 
                                     ncol= 1 + length(studiedFeatures)))##+1 for the affected gene column
  
  colnames(resultMatrix)<-c("affected_gene",studiedFeatures)
  rownames(resultMatrix)<-affectedGenes
  
  ##ADD affected_gene name column (easier to handle afterwards gene name for other analysis)
  resultMatrix$affected_gene<-affectedGenes
  
  ##########################################################
  ## Filling the results Matrix
  
  for (gene in affectedGenes) {
    
    ####################################
    ## Gene associated breakpoint
    if((gene %in% genesInfo$genesInDomain$domain_breakpoint_1)
       && (gene %in% genesInfo$genesInDomain$domain_breakpoint_2)){
      ##So the gene is associated to both breakpoints (can occur if intraTAD SV)
      resultMatrix[gene,"gene_Breakpoint"]<-"Break1&2"
    }else if(gene %in% genesInfo$genesInDomain$domain_breakpoint_1){
      resultMatrix[gene,"gene_Breakpoint"]<-"Break1"
    }else if (gene %in% genesInfo$genesInDomain$domain_breakpoint_2){
      resultMatrix[gene,"gene_Breakpoint"]<-"Break2"
    }
    
    #################################################################
    ##Filling main phenotype columns && n phenotypes patient related
    ##HUMAN DATA
    if(gene %in% rownames(humanPhenotype)){
      ##let's indicate it has an associated phenotype in OMIM
      resultMatrix[gene,"associatedPhenotypeIn_OMIM"]<-TRUE
      
      ##Main phenotype related column
      resultMatrix[gene,"mainPhenotype_Through_OMIM_Human"]<-humanPhenotype[gene,mainPatientPhenotype]==1    
      
      # ##All phentoypes patient
      # resultMatrix[gene,"n_PhenotypeRelated_Through_OMIM_Human"]<-sum(humanPhenotype[gene,patientPhenotypes])
      
    }else{
      ##NOT IN OMIM, by default we mark this as FALSE, but we could also put not in OMIM
      ##Altought they can be, but not in HPO classified, so that is the reason why we do not get them
      resultMatrix[gene,"associatedPhenotypeIn_OMIM"]<-FALSE
      
      resultMatrix[gene,"mainPhenotype_Through_OMIM_Human"]<-FALSE 
      
      # resultMatrix[gene,"n_PhenotypeRelated_Through_OMIM_Human"]<-0
    }
    
    ##MICE DATA
    if(gene %in% rownames(micePhenotype)){
      ##let's indicate it has an associated phenotype in MGI
      resultMatrix[gene,"associatedPhenotypeIn_MGI"]<-TRUE
      
      ##Main phenotype related column
      resultMatrix[gene,"mainPhenotype_Through_MGI_Mice"]<-micePhenotype[gene,mainPatientPhenotype]==1   
      
      # ##All phentoypes patient
      # resultMatrix[gene,"n_PhenotypeRelated_Through_MGI_Mice"]<-sum(micePhenotype[gene,patientPhenotypes])
    }else{
      ##NOT IN MGI, by default we mark this as FALSE, but we could also put not in OMIM
      ##Altought they can be, but not in MPO classified, so that is the reason why we do not get them
      resultMatrix[gene,"associatedPhenotypeIn_MGI"]<-FALSE
      
      resultMatrix[gene,"mainPhenotype_Through_MGI_Mice"]<-FALSE 
      
      # resultMatrix[gene,"n_PhenotypeRelated_Through_MGI_Mice"]<-0
    }
    
    
    #########################
    ## Type of Regulatory mechanism: Direct_disruptingGeneSequence Indirect_notDisruptingGeneSequence(Long Range enter here, unless duplications)
    #########################
    ##Defining also hierarchy to handle different scenarios
    #Hence if a gene is broken, it will not be considered as deleted or duplicated
    #THIS IS IMPORTANT
    if(gene %in% brokenGenes){
      ##if it is crossing the breakpoint segment limits
      ##we consider it as broken
      resultMatrix[gene,"RegulatoryMechanism"]<-"Direct_geneTruncation"
      
    }else if(gene %in% genesInUncertaintyRegions){
        ##if it is entirely in the uncertainty region of the breakpoint
        resultMatrix[gene,"RegulatoryMechanism"]<-"Direct_uncertaintyRegion"
        
    }else if(gene %in% deletedGenes){
      resultMatrix[gene,"RegulatoryMechanism"]<-"Direct_geneDeletion"
      
    }else if(gene %in% duplicatedGenes){
      ##A gene duplication can trigger an upregulation due to duplication of active gene
      ##or by neoTad formation which is a long range mechanism. That is why we here put Direct_LongRange
      ##or by enhancers duplications(but due to shadow enh... not very likely)
      ##Once scoring done this cell will be updated if longRange_geneDuplication more likely than Direct_geneDuplication... or if both equally likely
      resultMatrix[gene,"RegulatoryMechanism"]<-"Direct_LongRange_geneDuplication"
      
    }else{
        ##Esto es el final else
        ##LONG RANGE EFFECT
        resultMatrix[gene,"RegulatoryMechanism"]<-"LongRange" 
    }
    
    ##########################
    ## GENE FPKMs Retrieval
    ##########################
    if((gene %in% Master_GeneExpression$gene_name) && (phasesVector !="phaseFree")){
      ##Because we do not have expression data for phaseFree
      
      fpkms_list<-list()
      
      for(phase in phasesVector){
        phaseExpression<-Master_GeneExpression[[phase]][Master_GeneExpression$gene_name==gene]
        
        if(length(phaseExpression)==1){
          ##So it is not repeated the gene name
          ##Hence, just one expression value
          resultMatrix[gene,paste0("FPKM_",phase)]<-phaseExpression
        }else{
          ##So the same gene name has >1 expression value
          ##If all above or below thresholds, do the average. If they cross thresholds set as -2. Do not handle for the moment
          if((all(phaseExpression<=PoorlyExpThresh)) || (all(phaseExpression>=ConsiderablyExpThresh))){
            resultMatrix[gene,paste0("FPKM_",phase)]<-mean(phaseExpression)
          }else{
            resultMatrix[gene,paste0("FPKM_",phase)]<-(-2)##Disparity expression values
          }
        }
      }
    }else{
      ###assign -1 if expression not measured
      ##pa que no meta un NA y luego lo haga como factor... y convierta algun nro a factor...
      resultMatrix[gene,c(paste0("FPKM_",phasesVector))]<-(-1)
    }

    #####################################
    ##Filling Polycomb score column
    #####################################
    
    ##Get max in case 1 gene symbol repeated, to avoid bugs there
    polycScore<-unlist(polyCombGenes[polyCombGenes$`#Gene`==gene,"PolyC_score"])
    
    if(is.null(polycScore)){
      ##So not computed
      polycScore <- -1
    }else{
      ##Get max value, in case some gene symbol repeated, to avoid bugs here
      polycScore <- max(polycScore)
    }
    
    resultMatrix[gene,"polyComb_score"]<-polycScore
      

    ####################################
    ## Filling column HI scores
    ####################################
    
    if(gene %in% rownames(natureLof)){
      resultMatrix[gene,"nature_HI_score"]<-natureLof[gene,"pLI"]
    }else{
      ##We put -1, so it means, no reported score
      ##It is not the same an NA (-1 the way we define them) that a 0, que significa que se ha estudiado y la probabilidad es baja
      resultMatrix[gene,"nature_HI_score"]<-(-1)
    }
    
    if(gene %in% rownames(huangLof)){
      resultMatrix[gene,"huang_HI_score"]<-huangLof[gene,"pHI"]
    }else{
      ##We put -1, so it means, no reported score
      ##It is not the same an NA (-1 the way we define them) that a 0, que significa que se ha estudiado y la probabilidad es baja
      resultMatrix[gene,"huang_HI_score"]<-(-1)
    }
    
    ##ClinGene Hi info
    if(gene %in% rownames(clinGen_hiInfo)){
      resultMatrix[gene,"clinGene_HI_score"]<-clinGen_hiInfo[gene,"clinGene_01_score"]
    }else{
      ##We put -1 same criteria as for other HI datasets
      resultMatrix[gene,"clinGene_HI_score"]<-(-1)
    }
    
    #######################################
    ## Adding TAU_exp score
    #######################################
    if(gene %in% tau_exp_scores$gene_name){
      
      tauExp<-tau_exp_scores$tau_exp_score[tau_exp_scores$gene_name==gene]
      
      if(length(tauExp)==1){
        ##just one value
        resultMatrix[gene,"TAU_exp"]<-tauExp##Gene Expression Not Evaluated (We won't have FPKMs neither)
      }else{
        ##Multiple values TAU
        ##Do fucking average
        resultMatrix[gene,"TAU_exp"]<-mean(tauExp)
      }
      
    }else{
      tauExp<-(-2)
      resultMatrix[gene,"TAU_exp"]<-(-2)##Gene Expression Not Evaluated (We won't have FPKMs neither)
    }
    
    ##########################################################################################################
    ## Getting enhancer statistics. N enh gained. N enh lost etc.
    ## For now, we only want to do this for genes that are located in disrupted TADs,
    ## and hence can suffer long range rewiring
    ## So if there is a whole TAD deleted between two tads that are partially deleted
    ## For the genes that are located in a TAD entirely deleted or duplicated we do not want enhancer info
    ## If it is expressed... expression boost, if it is not... duplicating the TAD not apparent consecuences
    ###########################################################################################################
    if(gene %in% genesInDisruptedTAD){
      #######################################################################################
      ## Computing nEnh initially for each gene to the right and to the left part of the TAD
      ## the nEnh kept and gained is performed lateron
      ## Required for graphical representation plots
      ######################################################################################
      
      ###With respect to each gene, as it can vary, number of enhancer to the left and number of enhancers to the right 
      ##Initially
      if(gene %in% genesInfo$genesInDomain$domain_breakpoint_1){
        ##So the gene is in the breakpoint 1 Domain
        currentGeneDomain<-"domain_breakpoint_1"
      }else{
        ##So the gene is in the breakpoint 2 Domain
        currentGeneDomain<-"domain_breakpoint_2"
      }
      
      matrix_enhInDomain<-enhancersInfo$enhancersInDomain[[currentGeneDomain]]
      geneTSS<-info_affectedGenes$genesPosition[gene,"TSS"]
      
      for(sourceEnh in phasesVector){
        
        subs_enh<-subset(matrix_enhInDomain,source==sourceEnh)
        
        nEnhInitiallyToThe_Right<-sum(subs_enh$start > geneTSS) #take as reference the enh start
        nEnhInitiallyToThe_Left<-sum(subs_enh$start < geneTSS)
        
        acetilationEnhInitiallyToThe_Right<-sum(subs_enh$acetilation[subs_enh$start > geneTSS])
        acetilationEnhInitiallyToThe_Left<-sum(subs_enh$acetilation[subs_enh$start < geneTSS])

        ###########################################
        ####Adding numbers to the results matrix
        
        ###For the Right Side
        ###Colname in the results matrix
        colId<-paste0("nEnh_ToTheRight_initial_",sourceEnh)
        ##Add to all the genes in the domain the n enh for this stage
        resultMatrix[gene,colId]<-nEnhInitiallyToThe_Right
        
        ###For the Left Side
        ###Colname in the results matrix
        colId<-paste0("nEnh_ToTheLeft_initial_",sourceEnh)
        ##Add to all the genes in the domain the n enh for this stage
        resultMatrix[gene,colId]<-nEnhInitiallyToThe_Left
        
        ######For the acetilation enh levels
        ###For the Right Side
        ###Colname in the results matrix
        colId<-paste0("enhancers_acetilation_ToTheRight_initial_",sourceEnh)
        ##Add to all the genes in the domain the n enh for this stage
        resultMatrix[gene,colId]<-acetilationEnhInitiallyToThe_Right
        
        ###For the Left Side
        ###Colname in the results matrix
        colId<-paste0("enhancers_acetilation_ToTheLeft_initial_",sourceEnh)
        ##Add to all the genes in the domain the n enh for this stage
        resultMatrix[gene,colId]<-acetilationEnhInitiallyToThe_Left
      }
    }
    
  }
  
  ##Again the same in terms of relevant genes for enh calculations
  ##########################################################################################################
  ## Getting enhancer statistics. N enh gained. N enh lost etc.
  ## For now, we only want to do this for genes that are located in disrupted TADs,
  ## and hence can suffer long range rewiring
  ## So if there is a whole TAD deleted between two tads that are partially deleted
  ## For the genes that are located in a TAD entirely deleted or duplicated we do not want enhancer info
  ## If it is expressed... expression boost, if it is not... duplicating the TAD not apparent consecuences
  ###########################################################################################################
  
  ##We are constantly providing for the calculation of Nenh gained and Nenh lost the targetGenes. So the calculations only will appear
  ##For the genes we consider, hence we exclude those in TADs where there is no disruption
  
  
  #########################################
  ####N enhancers initially for each gene IN THEIR TAD
  #########################################
  ## INITIAL
  ## N enhancers initially
  ##Contamos enh en cada dominio cuantos enh hay, y los sumamos a los genes de ese dominio como Nro inicial
  ##Los contamos por stage (por complementariedad)
  for(domain in names(enhancersInfo$enhancersInDomain)){
    
    # print(domain)
    
    ###Enhancers in the domain
    matrix_enhInDomain<-enhancersInfo$enhancersInDomain[[domain]]
    # phasesVector<-unique(matrix_enhInDomain$source)
    phasesVector<-phasesVector
    
    ##Genes in the Domain
    
    genesInDomain<-genesInfo$genesInDomain[[domain]]
    
    for(sourceEnh in phasesVector){
      # print(sourceEnh)
      subs_enh<-subset(matrix_enhInDomain,source==sourceEnh)
      
      n_enh<-nrow(subs_enh)#number of enhancers 
      totalAcetilation<-sum(subs_enh$acetilation)#total acetilation of the enhancers
      
      ###Colname in the results matrix
      colId<-paste0("nEnhancers_initial_",sourceEnh)
      ##Add to all the genes in the domain the n enh for this stage
      resultMatrix[genesInDomain,colId]<-n_enh
      
      ##Add acetilation levels info
      colAcetilLevelId<-paste0("enhancers_acetilation_initial_",sourceEnh)
      resultMatrix[genesInDomain,colAcetilLevelId]<-totalAcetilation
      
    }
  }
  
  #########################################
  ####N enhancers max available for each gene
  #### IN THE OTHER DOMAIN
  #########################################
  ##"nEnhancers_maxAvailableInTheOtherDomain_"
  for(domain in names(enhancersInfo$enhancersInDomain)){
    
    pos_name_oppositeDomain<-which(names(genesInfo$genesInDomain)!=domain)
    
    # genesInOppositeDomain<-genesInfo$genesInDomain[[pos_name_oppositeDomain]]
    
    ##Genes in Domain
    genesInDomain<-genesInfo$genesInDomain[[domain]]
    
    ###Enhancers in the Opposite domain
    matrix_enhInOppositeDomain<-enhancersInfo$enhancersInDomain[[pos_name_oppositeDomain]]
    
    for(sourceEnh in phasesVector){
      # print(sourceEnh)
      subs_enh<-subset(matrix_enhInOppositeDomain,source==sourceEnh)
      
      n_enh<-nrow(subs_enh)#number of enhancers 
      totalAcetilation<-sum(subs_enh$acetilation)#total acetilation of the enhancers, in the OPPOSITE DOMAIN in this case
      
      ###Colname in the results matrix
      colId<-paste0("nEnhancers_maxAvailableInTheOtherDomain_",sourceEnh)
      ##Add to all the genes in the domain the n enh for this stage
      resultMatrix[genesInDomain,colId]<-n_enh
      
      ##Add acetilation levels info
      colAcetilLevelId<-paste0("enhancers_maxAcetilationAvailableInTheOtherDomain_",sourceEnh)
      resultMatrix[genesInDomain,colAcetilLevelId]<-totalAcetilation
      
    }
  }
  
  #############################################
  ## Type domain initial TAD or between_TAD
  #############################################
  ## INITIAL
  for(domain in rownames(info_affectedRegions$domainsAffected)){
    
    ##print(domain)##domain_breakpoint_1 || domain_breakpoint_2
    
    
    domainID<-info_affectedRegions$domainsAffected[domain,"domainId"]
    domainType<-info_affectedRegions$domainsAffected[domain,"domainType"]
    
    ###retrieve genes in that domain
    ##Genes in Domain
    genesInDomain<-genesInfo$genesInDomain[[domain]]
    
    ####Add in the matrix for the genes in the domain
    ##its domain type, and the identifier of the domain
    resultMatrix[genesInDomain,"TypeDomainInitial"]<-paste(domainType,
                                                           "_disrupted",
                                                           sep="")##To emphasize that the regulatory environment of the gene is disrupted
    
    resultMatrix[genesInDomain,"InitialDomain_ID"]<-domainID
    
  }
  
  ############################################################################################################################################################################
  ## Now, for those genes that have an NA on their TypeDomainInitial
  ## It means that their regulatory domain (e.g. TAD) or between TAD if they are between TADs is intact
  ## And they are associated to the SV because their TAD (or regulatory domain) has been entirely duplicated or deleted (if a TAD is completely inverted or translocated, not expected pathogenicity)
  ############################################################################################################################################################################
  resultMatrix$TypeDomainInitial[is.na(resultMatrix$TypeDomainInitial)]<-"regulatoryDomain_intact"
  
  
  ##Type domain initial...not sure if necessary
  
  ##Not used for now, I think, so just comment it
  # ## MIXED WITH
  # initialDomains<-unique(domainsInfo$domainsAffected$domainType)
  # 
  # if(length(initialDomains)==1){
  #   ##then either both bp are in TADs or both bp are in between_TADs
  #   ##so the type of domain with which you are gonna mix is the same one
  #   ##is.na to avoid putting this info for the genes in domains located in the middle
  #   ##by the way this info is not being used for now for anything
  #   resultMatrix$TypeDomainMixedWith[!is.na(resultMatrix$TypeDomainInitial)]<-initialDomains
  # }else{
  #   ##it means we are working with the 2 posible domains: TAD && between_TADs
  #   ##so the rows with type of domain A (TAD for instance), after the rearrangement will be mixed with the B type
  #   Adomain<-initialDomains[1]
  #   Bdomain<-initialDomains[2]
  #   
  #   resultMatrix[resultMatrix$TypeDomainInitial==Adomain, "TypeDomainMixedWith"]<-Bdomain
  #   resultMatrix[resultMatrix$TypeDomainInitial==Bdomain, "TypeDomainMixedWith"]<-Adomain    
  #   
  # }
  
  ##So just remove the column
  resultMatrix$TypeDomainMixedWith<-NULL
  
  ###############################
  ## n Enh Kept and Gained ####
  ###############################
  
  ##Model For inversions
  #N enhancers AFTER rearrangement
  ######1st for both breakpoints falling in TADs
  if(typeSV=="Inversion"){
    ##We  need to differentiate between intra TAD or between TADs // combining TADs
    ##Because modelling differs
    if(length(unique(domainsInfo$domainsAffected$domainId))>1){
      ##So there is more than one affected TAD, hence the SV is inter-TADs or inter-regulatory domains
      #####################################################
      ##after the inversion, lefts segments go together
      ##and right segments go together
      
      ##For genes(if there are) in breakp1 leftSegment
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_1_left_segment,
                                   keptEnh = enhancersInfo$enhancersInSegment$breakpoint_1_left_segment,
                                   gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_2_left_segment,
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)
      
      ##For genes in breakp2 leftSegment
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_2_left_segment,
                                   keptEnh = enhancersInfo$enhancersInSegment$breakpoint_2_left_segment,
                                   gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_1_left_segment,
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)
      
      ##For genes in breakp1 rightSegment
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_1_right_segment,
                                   keptEnh = enhancersInfo$enhancersInSegment$breakpoint_1_right_segment,
                                   gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_2_right_segment,
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)
      
      ##For genes in breakp2 rightSegment
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_2_right_segment,
                                   keptEnh = enhancersInfo$enhancersInSegment$breakpoint_2_right_segment,
                                   gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_1_right_segment,
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)
      
    }else if(length(unique(domainsInfo$domainsAffected$domainId))==1){
      ##So, it is an inversion strictyly occurring intra TAD, so unless gene or enh broken, no relevant change
      ##For now, that's our current view
      
      ##With only one calculation enough
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInDomain$domain_breakpoint_1,##domain is the same, so does not matter which we choose
                                   keptEnh = unique(rbind(enhancersInfo$enhancersInSegment$breakpoint_1_left_segment,
                                                   enhancersInfo$enhancersInSegment$breakpoint_1_right_segment,
                                                   enhancersInfo$enhancersInSegment$breakpoint_2_left_segment,
                                                   enhancersInfo$enhancersInSegment$breakpoint_2_right_segment)),#Here we do not take all into account all enh in the domain
                                                     ##because maybe some enh is broken
                                                     #And by looking at the enh per segment, if a enh is broken it will not be counted
                                                     #Because we associate an enh to a segment if it entirely falls inside
                                                     #And hence it will not be considered
                                   gainedEnh = NULL,
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)
    }


  }else if(typeSV=="Deletion"){
    
    ##We  need to differentiate between intra TAD or between TADs // combining TADs
    ##Because modelling differs 
    
    if(length(unique(domainsInfo$domainsAffected$domainId))>1){
      ##So there is more than one affected TAD, hence the SV is inter-TADs or inter-regulatory domains
      #####################################################
      
      ##Genes that are in bkp1 left segment, lose what is on bkp1 right, and gain what is in bkp2 right
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_1_left_segment,
                                   keptEnh = enhancersInfo$enhancersInSegment$breakpoint_1_left_segment,
                                   gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_2_right_segment,
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)
      
      ##For genes in breakp2 rightSegment
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_2_right_segment,
                                   keptEnh = enhancersInfo$enhancersInSegment$breakpoint_2_right_segment,
                                   gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_1_left_segment,
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)
      
      ##Genes that are in bkp1 right segm or bkp2 left segm are deleted, so directly affected
      ##Hence we do not consider enhancer situation for them
      
    }else if(length(unique(domainsInfo$domainsAffected$domainId))==1){
      ##So, it is an deletion strictyly occurring intra TAD
      ##In this context,and attending to the data we have, there can not be a GOF, so gained enh NULL. (unless some silencers deleted... but we don't have them annotated)
      
      ############################################
      ##For the long-range candidates. LEFT side
      ##On the one hand Genes that are to the LEFT of breakpoint1 
      
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_1_left_segment,
                                   keptEnh =rbind(enhancersInfo$enhancersInSegment$breakpoint_1_left_segment,
                                                  enhancersInfo$enhancersInSegment$breakpoint_2_right_segment),
                                   gainedEnh = NULL,
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)
      
      
      ############################################
      ##For the long-range candidates. RIGHT side
      ##On the one hand Genes that are to the RIGHT of breakpoint2 
      
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_2_right_segment,
                                   keptEnh =rbind(enhancersInfo$enhancersInSegment$breakpoint_2_right_segment,
                                                  enhancersInfo$enhancersInSegment$breakpoint_1_left_segment),
                                   gainedEnh = NULL,
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)
      
      ###############################################################
      ##Genes that are in betweenn are deleted
      ##So directly affected
      ##Hence we do not consider enhancer situation for them
      ## So there will be NAs on the cells associated to this values
    }
    
  }else if(typeSV=="Duplication"){
    
    ##We  need to differentiate between intra TAD or between TADs // combining TADs
    ##Because modelling differs 
    
    if(length(unique(domainsInfo$domainsAffected$domainId))>1){
      ##So there is more than one affected TAD, hence the SV is inter-TADs or inter-regulatory domains
      #####################################################
      
      ##Genes that are in bkp1 left segment, stay in an identical condition than in the WT 
      ##Unless an enhancer is broken by the breakpoint
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_1_left_segment,
                                   keptEnh = rbind(enhancersInfo$enhancersInSegment$breakpoint_1_left_segment,
                                                   enhancersInfo$enhancersInSegment$breakpoint_1_right_segment),
                                   gainedEnh = NULL,
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)
      
      ##The same applies for genes in breakp2 rightSegment
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_2_right_segment,
                                   keptEnh = rbind(enhancersInfo$enhancersInSegment$breakpoint_2_left_segment,
                                                   enhancersInfo$enhancersInSegment$breakpoint_2_right_segment),
                                   gainedEnh = NULL,
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)
      
      ##Genes that are in bkp1 right segm are, on the one hand duplicated
      ##But in addition, they gain interactions with what is located on the left side of the breakpoint 2
      ##So besides the direct impact there is a long range mechanism going on
      ##They also maintain the interactions with all its cognate enh (unless any broken)
      ##because there is a copy of its tad remaining exactly the same
      ##And the enhancers that are located on their segment are also duplicated
      
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_1_right_segment,
                                   keptEnh = rbind(enhancersInfo$enhancersInSegment$breakpoint_1_left_segment,
                                                   enhancersInfo$enhancersInSegment$breakpoint_1_right_segment,
                                                   enhancersInfo$enhancersInSegment$breakpoint_1_right_segment),
                                   gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_2_left_segment,##Gained enh reserved For NEW ENHANCERS (from different RegDomain)
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)
      
      ##Genes that are in bkp2 left segm are, on the one hand duplicated
      ##But in addition, they gain interactions with what is located on the right side of the breakpoint 1
      ##So besides the direct impact there is a long range mechanism going on
      ##They also maintain the interactions with all its cognate enh (unless any broken)
      ##because there is a copy of its tad remaining exactly the same
      ##And the enhancers that are located on their segment are also duplicated
      
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_2_left_segment,
                                   keptEnh = rbind(enhancersInfo$enhancersInSegment$breakpoint_2_left_segment,
                                                   enhancersInfo$enhancersInSegment$breakpoint_2_left_segment,
                                                   enhancersInfo$enhancersInSegment$breakpoint_2_right_segment),
                                   gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_1_right_segment, ##Gained enh reserved For NEW ENHANCERS ((from different RegDomain))
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)
  
      
    }else if(length(unique(domainsInfo$domainsAffected$domainId))==1){
      ##So, it is a duplication strictly occurring intra TAD
      ##In this context,and attending to the data we have, there can only be a GOF. (unless some silencers are duplicated... but we don't have them annotated)

      ###All genes keep all the enhancers (except if any enh is broken by a breakpoint, that is why we select them by segments which entirely contain them)
      ##Gained enhancers:Those which are in the intersection between bkp1 right segment and bkp2 left segment (duplication area)
      
      ##The same applies for the enhancers
      ##Getting duplicated enhancers
      duplicatedEnh<-rbind(enhancersInfo$enhancersInSegment$breakpoint_1_right_segment,
                           enhancersInfo$enhancersInSegment$breakpoint_2_left_segment)
      
      ##duplicated() returns true at the first duplicated value, and only for it, not for the match.
      #Perfect, we only want to capture one instance of the repeated ones in the duplicatedEnh dframe.Which are the enh in the duplication area.
      duplicatedEnh<-duplicatedEnh[duplicated(duplicatedEnh),]
      
      resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInDomain$domain_breakpoint_1,##domain is the same, so does not matter which we choose
                                   keptEnh = rbind(enhancersInfo$enhancersInSegment$breakpoint_1_left_segment,
                                                   duplicatedEnh,##enh located in the duplication area, counted only once
                                                   duplicatedEnh,
                                                   enhancersInfo$enhancersInSegment$breakpoint_2_right_segment),
                                   gainedEnh = NULL,##Gained enh reserved For NEW ENHANCERS ((from different RegDomain))
                                   matrixRes = resultMatrix,
                                   enhOrigin = phasesVector,
                                   genesInfo = genesInfo)

    }
    
  }else if(typeSV=="Translocation"){
    
    ###Loading centromere locations for hg19   
    ####load centromere positions
    ##R object generated in script: /home/victor/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/CentromeresPosition/processingCentromereData.R
    load(file = "data/hg19_centromerePositions.RData")
    
    ##We consider left coords as the smaller of the breakpoints coord, and right coord the bigger
    
    ##For breakpoint 1 Information
    bp1_leftCoord<-info_affectedRegions$segmentsBreakP["breakpoint_1_left_segment","end"]
    bp1_rightCoord<-info_affectedRegions$segmentsBreakP["breakpoint_1_right_segment","start"]
    bp1_chr<-info_affectedRegions$domainsAffected["domain_breakpoint_1","chr"]
    bp1_centromerCoord<-centromerePos[bp1_chr,"centromerCenter"]##for this chromosome where the bp (breakpoint) occurs, where is the centromere
    
    ##For breakpoint 2 Information
    bp2_leftCoord<-info_affectedRegions$segmentsBreakP["breakpoint_2_left_segment","end"]
    bp2_rightCoord<-info_affectedRegions$segmentsBreakP["breakpoint_2_right_segment","start"]
    bp2_chr<-info_affectedRegions$domainsAffected["domain_breakpoint_2","chr"]
    bp2_centromerCoord<-centromerePos[bp2_chr,"centromerCenter"]##for this chromosome where the bp occurs, where is the centromere
    
    ######If centromere in between put Fuck Centromere in between
    ##Fer un log file... whatever to track this shits
    if((bp1_centromerCoord<bp1_rightCoord && bp1_centromerCoord>bp1_leftCoord) || (bp2_centromerCoord<bp2_rightCoord && bp2_centromerCoord>bp2_leftCoord)){
      ##it means that at least one centromere is between Der breakpoints
      print(enf)
      print("WARNING, centromere between coordinates for a Breakpoint")
      ##Anyadir por si esta el breakpoint en si dentro del centromero??? Aunque dudo que se pueda mapear un breakpoint si cae en el centromero
      
    }else{
      #okey, both breakpoints coordinates located  at one side of the centromere center
      #after breaking dna one will keep the centromere and the other not
      #lets make the rearrangements with them
      ##Based on the assumption that the rearrangemets are done so that the parts connected with centromeres connect with the fragments that do not have it
      
      infoCarrying<-character()##Carry means carries centromer
      if(bp1_leftCoord > bp1_centromerCoord){
        ##It means left coord is after centromer, so this fragment is the one retaining the centromere
        infoCarryingBp1<-c("bp1_left_Carry")##Info carrying breakpoint 1
      }else{
        ##It means bp1_right segment carries the centromer
        infoCarryingBp1<-c("bp1_left_NOT_Carry")##Info carrying breakpoint 1 
      }
      
      
      if(bp2_leftCoord > bp2_centromerCoord){
        ##It means left coord is after centromer, so this fragment is the one retaining the centromere
        infoCarryingBp2<-c("bp2_left_Carry")##Info carrying breakpoint 2
      }else{
        ##It means bp2_right segment carries the centromer
        infoCarryingBp2<-c("bp2_left_NOT_Carry")##Info carrying breakpoint 2
      }
      
      infoCarrying<-c(infoCarryingBp1,infoCarryingBp2)
      
      #########################################################################
      ## TWO cases. Left segments carry centromer or do not carry
      ## This implies that right segments have or do not have the centromers
      ## So they have to cross with them to preserve centromeric equilibrium
      #########################################################################
      
      if(all(infoCarrying==c("bp1_left_Carry","bp2_left_Carry")) || all(infoCarrying==c("bp1_left_NOT_Carry","bp2_left_NOT_Carry"))){
        
        ##So left segments of bp1 are crossed with right segments of bp2
        ##and the other way round
        
        ##For genes(if there are) in breakp1 leftSegment
        ##They gain enh from right segment of the other breakpoint
        resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_1_left_segment,
                                     keptEnh = enhancersInfo$enhancersInSegment$breakpoint_1_left_segment,
                                     gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_2_right_segment,
                                     matrixRes = resultMatrix,
                                     enhOrigin = phasesVector,
                                     genesInfo = genesInfo)
        
        ##For genes in breakp2 leftSegment
        resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_2_left_segment,
                                     keptEnh = enhancersInfo$enhancersInSegment$breakpoint_2_left_segment,
                                     gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_1_right_segment,
                                     matrixRes = resultMatrix,
                                     enhOrigin = phasesVector,
                                     genesInfo = genesInfo)
        
        ##For genes in breakp1 rightSegment
        resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_1_right_segment,
                                     keptEnh = enhancersInfo$enhancersInSegment$breakpoint_1_right_segment,
                                     gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_2_left_segment,
                                     matrixRes = resultMatrix,
                                     enhOrigin = phasesVector,
                                     genesInfo = genesInfo)
        
        ##For genes in breakp2 rightSegment
        resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_2_right_segment,
                                     keptEnh = enhancersInfo$enhancersInSegment$breakpoint_2_right_segment,
                                     gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_1_left_segment,
                                     matrixRes = resultMatrix,
                                     enhOrigin = phasesVector,
                                     genesInfo = genesInfo)
      }
      
      ############################################################
      ## CASE 3 / 4 && 4/4
      ############################################################
      
      if(all(infoCarrying==c("bp1_left_Carry","bp2_left_NOT_Carry")) || all(infoCarrying==c("bp1_left_NOT_Carry","bp2_left_Carry"))){
        
        ##So left segments of bp1 are crossed with left segments of bp2
        ##and the other way round
        
        ##For genes(if there are) in breakp1 leftSegment
        resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_1_left_segment,
                                     keptEnh = enhancersInfo$enhancersInSegment$breakpoint_1_left_segment,
                                     gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_2_left_segment,
                                     matrixRes = resultMatrix,
                                     enhOrigin = phasesVector,
                                     genesInfo = genesInfo)
        
        ##For genes in breakp2 leftSegment
        resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_2_left_segment,
                                     keptEnh = enhancersInfo$enhancersInSegment$breakpoint_2_left_segment,
                                     gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_1_left_segment,
                                     matrixRes = resultMatrix,
                                     enhOrigin = phasesVector,
                                     genesInfo = genesInfo)
        
        ##For genes in breakp1 rightSegment
        resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_1_right_segment,
                                     keptEnh = enhancersInfo$enhancersInSegment$breakpoint_1_right_segment,
                                     gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_2_right_segment,
                                     matrixRes = resultMatrix,
                                     enhOrigin = phasesVector,
                                     genesInfo = genesInfo)
        
        ##For genes in breakp2 rightSegment
        resultMatrix<-calculatingEnh(targetGenes = genesInfo$genesInSegment$breakpoint_2_right_segment,
                                     keptEnh = enhancersInfo$enhancersInSegment$breakpoint_2_right_segment,
                                     gainedEnh = enhancersInfo$enhancersInSegment$breakpoint_1_right_segment,
                                     matrixRes = resultMatrix,
                                     enhOrigin = phasesVector,
                                     genesInfo = genesInfo)
      }
    }
  }

  ######################################
  ## Return matrix with the information
  ######################################
  return(resultMatrix)
}