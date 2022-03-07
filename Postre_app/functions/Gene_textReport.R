##############################################
## Function to create The gene report ########
##############################################

##We need the table of the genes associated with the phenotypes in both human and MICE
gene_textReport<-function(patientResults, minPatogenicScore, mainPhenotype, targetGene){
  ##We need the results for just one condition, for all the genes.
  ##Min patogenicScore, defines the minimum patogenic score to show gene report

  ########################################
  ## Loading gene-phenotype relationships
  ########################################
  
  ##Gene-phenotype relationships based in Human
  load(file ="data/gene_fenotipe_basedOn_hpoAndOmim.RData")
  humanBased_genePhenotype<-gene_fenotipe_table
  rm(gene_fenotipe_table)
  
  ##Human orthologues,Based on Mice genes and Mice Phenotypes
  ##10471 genes
  load(file ="data/gene_fenotipe_basedOn_mgi_and_MP.RData")
  miceBased_genePhenotype<-gene_fenotipe_table
  rm(gene_fenotipe_table)
  
  ###############################################################################
  ## Loading info linking genes with databases entries to provide webpages links
  ###############################################################################
  
  #####################################
  ### based on human databases
  
  ##from HPO web, including orphanet data
  ##This dframe also contains for each gene-report the HPOs associated
  ##So when working in a specific phenotype use to determine relevant reports (those where main pheno is present)
  load(file="data/gene_phenotypeInfo_hpoWeb.RData")
  
  ### lets load info linking OMIM ids with CUIs (used to build Medgen links)
  load(file = "data/medgene_Info.RData")
  
  
  ##Load  also hpo terms associated per phenotype
  ##Used to filter reports for those where main phenotype is involved
  ##I've seen this to be an issue with human data (Orphanet for instance and sox9 returning entries related with sexual stuff, not head neck)
  load(file = "data/hpos_perFenotipe_HUMAN.RData")
  
  #####################################
  ### based on mice databases
  ### lets load info linking GENE with MGI entry
  load(file = "data/table_humanGene_mouseGene_MGI_id.RData")
  
  ##Apply Any Filtering of minimum threshold  
  targetMatrix<-patientResults$masterSummaryResultsMatrix
  
  #genesToWriteReport<-rownames(targetMatrix)[targetMatrix$maxScore>=minPatogenicScore]
  genesToWriteReport<-targetGene
  
  gene_reports<-list()
  
  # ## initial part for geneText_html object
  ## DEPRECATED UNNECESSARY
  # headerPart<-##"<html> 
  # "<head>
  # <meta name='viewport' content='width=device-width, initial-scale=1'>
  # <style>
  # </style>
  # </head>
  # <body>"
  
  headerPart<-""
  
  
  ##Let's remove all line breaks, introduced from internet copypaste java & html code
  headerPart<-gsub("[\r\n]", "", headerPart)
  
  for (gene in genesToWriteReport) {

    geneText<-character()
    
    #geneText<-paste(... = "<h1>",gene,"</h1>", collapse = "")
    
    ################################
    ##Phenotypic Part
    ################################
    
    ###Gene previous association with patient disease
    # geneText<-paste(geneText, "<h2>Gene previous association with patient phenotype</h2>", collapse = "")
    
    ##Start:Content of the section start
    geneText<-paste(geneText, "<div class='geneHumanPhenotype_content'>", collapse = "")
    
    #########
    # HUMANS
    #########
    ##Gathering links info
    
    
    ##ORPHANET & OMIM From HPO web and orphanet
    gene_webInfo<-subset(gene_phenotype, entrez_gene_symbol == gene)
    
    ##Variable used if true to filter downstream link results
    geneAssociatedMainPhenotype<-FALSE
    
    if(nrow(gene_webInfo)>=1 ##so there is a WEB entry at least
       ){
      ##Hence, there is gene-phenotype information
      
      ##But maybe, that OMIM or MEDGEN case, was not phenotypically categorized by HPO
      if(gene %in% rownames(humanBased_genePhenotype)){
        #So the gene is in the human-HPO phenotype categorized matrix

        if(humanBased_genePhenotype[gene,mainPhenotype]==1){
          geneText<-paste(geneText, "<p style='text-align:justify;'>",
                          targetGene,
                          " has been previously associated with the patient phenotype category in humans.</p>", 
                          collapse = "")
          
          ###Ensure we only mantain entries associated with the main phenotype
          ##To avoid providing entries associated with different phenotypes such as sex reversal for SOX9 in orphanet where nothing mentioned about craniofacial
          geneAssociatedMainPhenotype<-TRUE ##So an additional filtering of the links will be applied to provide only those with Main Pheno HPO terms associated
          
          
        }else{
          geneText<-paste(geneText, "<p>Eventhough ", 
                          targetGene,
                          " has not been linked with the patient phenotype category, it has been associated with other pathologies.</p>", collapse = "")
          
        }
      }else{
        ##Indicate that the gene presents an entry either in MedGen or in OMIM
        ##But we have not been able to get its standarized phenotype information
      
        ##Introduce space lines between elements if it gets two close
        geneText<-paste(geneText, "<p></p>", collapse = "")
        geneText<-paste(geneText, "<p>We are sorry to inform you that altough ",
                        targetGene,
                        " presents a link to a disease, we could not get any standarized phenotype information. So, we could not take it into consideration. Feel free to explore it.</p>",
                        "<p>If you are sure it is associated to a phenotype tell us, and we will take it into account.</p>",
                        collapse = "")
      }
      
      
      #######################################
      ## Now, Gathering and Adding WEB LINKS
      #######################################
      ##If gene associated with main phenotype, filter results to avoid providing reports with the gene involved but unrelated with main pheno
      hpos_MainPheno<-hpos_PerFenotipe[[mainPhenotype]]
      
      if(geneAssociatedMainPhenotype == TRUE){
        ##Filter gene webInfo attending to HPOs related with main phenotype
        gene_webInfo<-gene_webInfo[gene_webInfo$HPO_Term_ID %in% hpos_MainPheno,]
      }
      
      ##############################################################################################################################
      ## Go on with mim and orpha terms left, filtered or unfilitered if no relevant gene-pheno found, hence all results returned
      
      ##two possible sources mim2gene orphadata
      gene_webInfo_OMIM<-subset(gene_webInfo, G_D_source =="mim2gene")$disease_ID_for_link
      gene_webInfo_orphanet<-subset(gene_webInfo, G_D_source =="orphadata")$disease_ID_for_link
      ##We need to remove OMIM and ORPHA prefixes (everything before the :)
      gene_webInfo_OMIM<-gsub(pattern = "^.*:", replacement = "", x=gene_webInfo_OMIM)
      gene_webInfo_orphanet<-gsub(pattern = "^.*:", replacement = "", x=gene_webInfo_orphanet)
      
      #######Gathering entries per webpage
      #OMIM
      omim_entries<-unique(gene_webInfo_OMIM)
      
      #MEDGENE (we gather the CUI identifier to build the link to medgene)
      medgene_entries<-medgene_info[omim_entries,"X.OMIM_CUI"]
      medgene_entries<-unique(medgene_entries)
      
      ##ORPHANET
      orphanet_entries<-unique(gene_webInfo_orphanet)
      
      webEntries<-c(omim_entries, medgene_entries, orphanet_entries)
      
      ##############################
      ####Adding OMIM INFO
      if(any(omim_entries!="NULL")){
        ##"NULL" is the way absence of info is annotated on the table
        ##It occurs for MedGene ids, in case it happens in the future for MIM ids, let's keep it
        geneText<-paste(geneText, "<p>+ info in OMIM:</p>",collapse = "")
        
        ####Create html link
        all_urls<-paste("https://www.omim.org/entry/",omim_entries, sep = "")
        nEntry<-0
        for(entry in all_urls){
          nEntry<-nEntry+1
          geneText<-paste(geneText, "<p class='linkEntry'><a href='",entry,"' target='_blank'>","OMIM Entry:",nEntry,"</a></p>", 
                          collapse = "")
        }
      }
      
      ##############################
      ####Adding MedGen INFO
      ##"NULL" is the way absence of info is annotated on the table
      if(any(medgene_entries!="NULL")){
        ##Introduce space lines between elements if it gets two close
        geneText<-paste(geneText, "<p></p>", collapse = "")
        geneText<-paste(geneText, "<p>+ info in MedGen:</p>",collapse = "")
        
        ####Create html link
        all_urls<-paste("https://www.ncbi.nlm.nih.gov/medgen/",unique(medgene_entries), sep = "")
        
        nEntry<-0
        for(entry in all_urls){
          nEntry<-nEntry+1
          geneText<-paste(geneText, "<p class='linkEntry'><a href='",entry,"' target='_blank'>","MedGen Entry:",nEntry,"</a></p>", 
                          collapse = "")
        } 
      }
      
      ##############################
      ####Adding Orphanet INFO
      if(any(orphanet_entries!="NULL")){
        ##Introduce space lines between elements if it gets two close
        geneText<-paste(geneText, "<p></p>", collapse = "")
        geneText<-paste(geneText, "<p>+ info in Orphanet:</p>",collapse = "")
        
        ####Create html link
        all_urls<-paste("https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=EN&Expert=",unique(orphanet_entries), sep = "")
        
        nEntry<-0
        for(entry in all_urls){
          nEntry<-nEntry+1
          geneText<-paste(geneText, "<p class='linkEntry'><a href='",entry,"' target='_blank'>","Orphanet Entry:",nEntry,"</a></p>", 
                          collapse = "")
          
        } 
        
      }
      
    }else{
      geneText<-paste(geneText, "<p>There is NO phenotypic information.</p>", collapse = "")
    }
    
    ##Closing humanPhenotype section
    geneText<-paste(geneText,
                    "</div>", 
                    collapse = "")
    
    # ##Introduce space lines between elements if it gets two close
    # geneText<-paste(geneText, "<p></p>", collapse = "")
    
    
    #########
    # MICE
    #########

    geneText<-paste(geneText, "<div class='geneMicePhenotype_content'>", collapse = "")
    
    gene_webInfo<-subset(table_humanGene_mouseGene_MGI_id,human_symbol==gene)
    if(nrow(gene_webInfo)>=1 ##so there is an MGI id
    ){
      ##Hence, there is information from MGI
      
      ##But maybe that MGI case was not phenotypically categorized by MP or no phentoype found
      if(gene %in% rownames(miceBased_genePhenotype)){
        #So the gene is in the mice-phenotype categorized matrix
        ##It can have or not, an association with phenotype
        if(miceBased_genePhenotype[gene,mainPhenotype]==1){
          geneText<-paste(geneText, "<p style='text-align:justify;'>The mouse homologous gene for ",
                          targetGene,
                          " has been associated in mice with a phenotype category similar to the one reported in the patient.</p>", 
                          collapse = "")
          
        }else{
          if(rowSums(miceBased_genePhenotype[gene,])>=1){
            #So the gene is associated with at least a phenotype (but not the main one)
            geneText<-paste(geneText, "<p>Eventhough the mouse homologous gene for ", 
                            targetGene,
                            " has not been linked with a phenotype category similar to the one reported in the patient, it has been associated with other pathologies.</p>", collapse = "")
            
          }
        }
      }else{
        ##Indicate that the gene presents an entry in MGI
        ##But we have not been able to get its standarized phenotype information (if there is)
        
        ##Introduce space lines between elements if it gets two close
        geneText<-paste(geneText, "<p></p>", collapse = "")
        geneText<-paste(geneText, "<p>We are sorry to inform you that altough the mouse homologous gene for ",
                        targetGene,
                        " presents a link to an MGI entry, we could not get any standarized phenotype information. So, we could not take it into consideration. Feel free to explore it.</p>",
                        "<p>If you are sure it is associated to a phenotype tell us, and we will take it into account.</p>",
                        collapse = "")
      }
      
      ##############################
      ####Adding MGI INFO Web Links
      ##############################
      
      if(any(gene_webInfo$MGI_id!="NULL")){
        
        geneText<-paste(geneText, "<p>+ info in MGI:</p>",collapse = "")
        
        ####Create html link
        all_urls<-paste("http://www.informatics.jax.org/marker/",unique(gene_webInfo$MGI_id), sep = "")
        nEntry<-0
        for(entry in all_urls){
          nEntry<-nEntry+1
          geneText<-paste(geneText, "<p class='linkEntry'><a href='",entry,"' target='_blank'>","MGI Entry:",nEntry,"</a></p>", 
                          collapse = "")
        }
      }
      
    }else{
      geneText<-paste(geneText, "<p>There is NO phenotypic information.</p>", collapse = "")
    }
    
    ##Closing micePhenotype section
    geneText<-paste(geneText,
                    "</div>", 
                    collapse = "")
    
    #######################################
    ##Adding links to human and mice logos
    info_humanLogo<-"<div class='logoHuman_location'><img class='genePhenoImgs' src='personLogo.png' alt='human logo' style='width:100%;height:100%;'></div>"
    info_miceLogo<-"<div class='logoMice_location'><img class='genePhenoImgs' src='miceLogo.png' alt='mice logo' style='width:100%;height:100%;'></div>"
    
    ##Embedding geneText in a div to apply a grid css configuration
    
    geneText<-paste("<div class='genePhenotype_section'>",
                    geneText, 
                    info_humanLogo,
                    info_miceLogo,
                    "</div>", 
                    collapse = "")
    
    # ##Introduce space lines between elements if it gets two close
    # geneText<-paste(geneText, "<p></p>", collapse = "")
    
    
    ##################################
    ##Store gene html in list
    ##################################
    gene_reports[[gene]]<-geneText
    
  }  
  
  ###Here we have a list with an html code per gene explaining what is going on with it
  ###It is good to have the object in this way because in the future if we want to format the code by blocks
  ###Here, we have it by blocks per gene
  ##Esto lo tengo que integrar https://www.w3schools.com/howto/tryit.asp?filename=tryhow_js_collapsible
  geneMasterHTML<-paste(unlist(gene_reports),sep = "",collapse = "")

  ##AddingCollapsableList INfo
  # endPart<-"</body>"##Deprecated. More info to be added
  endPart<-""
  ##Let's remove all line breaks, introduced from internet copypaste java & html code
  endPart<-gsub("[\r\n]", "", endPart)
  
  geneMasterHTML<-paste(headerPart, geneMasterHTML,endPart,sep = "",collapse = "")

  #####################
  ## return geneMasterHTML
  return(geneMasterHTML)
  
}
