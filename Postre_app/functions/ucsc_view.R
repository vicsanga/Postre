#########################################################################
## Function to create an html with different links
## For the user to visualize in UCSC the differently
## affected domains, with colored regions to easily visualize them
#########################################################################


ucsc_view<-function(patientResults, browserSessionId, targetGene, devStage){
  ##Table holding genome browser links
  genomeBrowser_links<-patientResults$genomeBrowser_links
  
  ##Zoom Out Distance From Domain
  ###How many surrounding bp we want to see in bp
  zoomOutDist<-500000
  
  borderHighlight_extension<-500 ##Extension of the TAD border highlight in bp, if put at single bp probably difficult to see in large zoom out
  
  breakpointHighlight_extension<-500 
  
  ##Check for how many domains affected,
  ##If only 1, only 1 link required. If 2 two links
  ##Get the coordinates of the affected TADs regarding all TAD maps used
  # names(patientResults$resultsPerTADmap)
  
  ##################
  ## To start with the approach let's paint with respect to the first map, for the devStage of interest
  affectedRegions<-patientResults$resultsPerPhase_secondaryInfo[[devStage]][[1]]$affectedRegions  
  
  ##Unique, if everyting happens in the same TAD, only 1 link provided
  ##But in that case we need to provide and paint 2 breakpoints in the same TAD
  
  domainsCoord<-unique(affectedRegions$domainsAffected)
  breakpCoord<-unique(affectedRegions$segmentsBreakP)
  
  ###########################
  ## Starting HTML ##########
  shiny_html_ucsc<-"<p style='font-size:25px'><b>UCSC links to Affected Domains in the Control condition</b></p>"
  # shiny_html_ucsc<-paste(shiny_html_ucsc,
  #                        "<p>For each affected Domain, it is depicted in blue the regulatory domain limits, and in red the breakpoint location.</p>",
  #                        sep = "",
  #                        collapse = "")
  
  ucsc_domain_tracks<-list()
  
  #######################################################################
  ##Aqui es donde deberia diferenciar entre Del|Dup, vs Inv|Transloc 
  ##Si es inv or translocation pintamos 2 tads si 2 tads affected
  ##Ya que tb lo de en medio de los TADs no ha sido disrupted
  ##Si es deletion o duplication pintamos 1 link solo, aunque se afecten multiples TADs
  #Ya que todo lo que este en medio de los breakp, alterado
  
  sv_type<-patientResults$patientInfo$TypeSV
  
  
  if((sv_type == "Inversion") || (sv_type == "Translocation")){
    
    ##Here we paint the two affected domains per separate, hence in >1 link (if >1 TAD affected)
    
    for(nDomain in 1:nrow(domainsCoord)){
      ##beginningLink
      ##Here take the UCSC base name link, dependent on phenotype and Phase
      ## eg "https://genome.ucsc.edu/s/sgvic/SV_RadaLab"
      linkBaseName<-genomeBrowser_links$SessionLinkBaseName[genomeBrowser_links$SessionId==browserSessionId]
      
      linkString<-paste(linkBaseName,
                        "?db=hg19&position=",
                        sep="")   ##This is the require beggining, we are directing to 
      
      ##To display, extend TAD coordinates +- 10Kb
      start_browser<-domainsCoord[nDomain,"start"] - zoomOutDist
      end_browser<-domainsCoord[nDomain,"end"] + zoomOutDist
      chr_browser<-domainsCoord[nDomain,"chr"]
      ##Formatting chr coordinates
      
      formated_browserView<-paste(chr_browser,":",
                                  start_browser,"-", 
                                  end_browser, sep = "")
      
      linkString<-paste(linkString,formated_browserView,sep = "" ,collapse = "")
      
      ##We do not want user cookies to be loaded that can interfere with our display
      linkString<-paste(linkString,"&ignoreCookie=1",sep = "" ,collapse = "")
      
      ##############################################
      ##Adding Higlight 
      ##############################################
      
      browser_highlight<-"&highlight="
      
      ###############################
      ## ADDING BREAKPOINT HIGHLIGHT
      ###############################
      
      ##If we have just 1 row in domainsCoord, it implies
      ##That everyting occurs in the same TAD,
      ##Hence in that case, we still have two breakpoints
      ##So add the highlight for both breakpoints in the same TAD
      
      if(nrow(domainsCoord)>1){
        
        ##So we have the breakpoints in 2 different tracks, hence we will create 2 tracks in the output html
        
        ##Let's retrieve breakpoint ID to get subsequently its coord
        breakpoint_ID<-unlist(strsplit(x = rownames(domainsCoord)[nDomain], split = "domain_",fixed = TRUE))[2]
        
        ##Now let's get both breakpointLimits
        breakpoint_start<-breakpCoord[paste(breakpoint_ID,"_left_segment",collapse = "",sep = "")
                                      ,"end"]
        
        breakpoint_end<-breakpCoord[paste(breakpoint_ID,"_right_segment",collapse = "",sep = "")
                                    ,"start"]
        
        ##########
        ## Extending coord to highlight a bigger region and facilitate visualization
        breakp_HighlightStart<-breakpoint_start - breakpointHighlight_extension
        breakp_HighlightEnd<-breakpoint_end + breakpointHighlight_extension
        
        breakP_region<-paste("hg19.",
                             chr_browser,"%3A",
                             breakp_HighlightStart,"-",
                             breakp_HighlightEnd,
                             collapse = "",
                             sep = "")
        
        ##RegionColor For the Regions ##BLUE
        colorBreakp<-paste("%23",##%23 is essential is the conversion somehow to url syntax of #
                           "ff3300",
                           sep = "",
                           collapse = "")
        
        ##Insterting breakp highlight on track
        browser_highlight<-paste(browser_highlight,
                                 "%7C",
                                 breakP_region,
                                 colorBreakp,
                                 sep="",
                                 collapse = "")
      }else{
        if(nrow(domainsCoord)==1){
          ##It implies that everything occur in the same TAD , so we need to put both breakpoints highlights in the same track
          
          for(breakpoint_ID in c("breakpoint_1","breakpoint_2")){
            ##Now let's get both breakpointLimits
            breakpoint_start<-breakpCoord[paste(breakpoint_ID,"_left_segment",collapse = "",sep = "")
                                          ,"end"]
            
            breakpoint_end<-breakpCoord[paste(breakpoint_ID,"_right_segment",collapse = "",sep = "")
                                        ,"start"]
            
            ##########
            ## Extending coord to highlight a bigger region and facilitate visualization
            breakp_HighlightStart<-breakpoint_start - breakpointHighlight_extension
            breakp_HighlightEnd<-breakpoint_end + breakpointHighlight_extension
            
            breakP_region<-paste("hg19.",
                                 chr_browser,"%3A",
                                 breakp_HighlightStart,"-",
                                 breakp_HighlightEnd,
                                 collapse = "",
                                 sep = "")
            
            ##RegionColor For the Regions ##BLUE
            colorBreakp<-paste("%23",##%23 is essential is the conversion somehow to url syntax of #
                               "ff3300",
                               sep = "",
                               collapse = "")
            
            ##Insterting breakp highlight on track
            browser_highlight<-paste(browser_highlight,
                                     "%7C",
                                     breakP_region,
                                     colorBreakp,
                                     sep="",
                                     collapse = "")
            
          }
        }
      }

      #############################################################
      ## Highligthing All enhancers in the affected TAD
      #############################################################
      
      masterEnh<-patientResults$MasterEnh_map
      masterEnh<-subset(masterEnh, source==devStage)
      #subset for genome region
      masterEnh$isChr<-masterEnh$chr==chr_browser
      
      ##Before displaying all enhancers in browser displayed region. Now filtering only for TAD where disease gene is present.
      # masterEnh$isBeforeEnd<-masterEnh$end<=end_browser
      # masterEnh$isAfterStart<-masterEnh$start>=start_browser
      # start_browser<-domainsCoord[nDomain,"start"] - zoomOutDist
      # end_browser<-domainsCoord[nDomain,"end"] + zoomOutDist
      
      ##Usar midpoint en lugar de enhancer start-end?? We are considering whole enh to be inside TADs on prediction scripts, so continue as it is.
      masterEnh$isBeforeEnd<-masterEnh$end<=domainsCoord[nDomain,"end"]
      masterEnh$isAfterStart<-masterEnh$start>=domainsCoord[nDomain,"start"]
      
      
      
      targetEnh<-rowSums(masterEnh[,c("isChr","isBeforeEnd","isAfterStart")])==3
      ##filter for targetEnh
      masterEnh<-masterEnh[targetEnh,]
      
      if(nrow(masterEnh)>0){
        #so at least one enh to highlight
        for(nEnh in 1:nrow(masterEnh)){
          #print(nEnh)
          #enhChr we know, is the chr_browser
          enh_start<-masterEnh[nEnh,"start"]
          enh_end<-masterEnh[nEnh,"end"]
          
          
          ##Add ucsc info to highlight
          enh_region<-paste("hg19.",
                            chr_browser,"%3A",
                            enh_start,"-",
                            enh_end,
                            collapse = "",
                            sep = "")
          
          ##RegionColor For the Regions #b3ffb3 33cc33 85e085
          colorEnh<-paste("%23",##%23 is essential is the conversion somehow to url syntax of #
                          "85e085",
                          sep = "",
                          collapse = "")
          
          ##Insterting breakp highlight on track
          browser_highlight<-paste(browser_highlight,
                                   "%7C",
                                   enh_region,
                                   colorEnh,
                                   sep="",
                                   collapse = "")
        }
      }
      
      #############################################################
      ## Highligthing The Target Gene if its domain is being linked
      #############################################################
      genesInfo<-patientResults$allAffectedGenes_positionalInfo
      genesInfo<-genesInfo[targetGene,]
      gene_chr<-genesInfo$chr
      gene_start<-genesInfo$start
      gene_end<-genesInfo$end
      gene_tss<-genesInfo$TSS
      
      ##################################################################
      ##if currently the genome browser is covering the target gene
      ##then highlight it
      ##################################################################
      coveringCurrentlyTargetGene<-FALSE##To track and customize UCSC link when showing multiple links, to highlight the one of the targetGene
      if((gene_chr==chr_browser) && (gene_tss>=start_browser) && (gene_tss<=end_browser)){
        ##highlight the gene body, as it is included in the currently created link area
        ##To track and customize UCSC link when showing multiple links, to highlight the one of the targetGene
        coveringCurrentlyTargetGene<-TRUE
        ##Add ucsc info to highlight
        gene_region<-paste("hg19.",
                           chr_browser,"%3A",
                           gene_start,"-",
                           gene_end,
                           collapse = "",
                           sep = "")
        
        ##RegionColor For the Regions #3399ff
        colorGene<-paste("%23",##%23 is essential is the conversion somehow to url syntax of #
                         "3399ff",
                         sep = "",
                         collapse = "")
        
        ##Insterting breakp highlight on track
        browser_highlight<-paste(browser_highlight,
                                 "%7C",
                                 gene_region,
                                 colorGene,
                                 sep="",
                                 collapse = "")
      }
      
      
      ########################################################################
      ### IF WE ARE DEALING WITH TRUNCATION
      ## So, the gene partially altered
      ## PAINT AGAIN THE BREAKPOINT TO AVOID BEING HIDDEN BY THE GENE BLUEISH
      ########################################################################
      
      if(patientResults$masterSummaryResultsMatrix[targetGene,"GeneImpact"] == "Direct_geneTruncation"){
        ###############################
        ## ADDING BREAKPOINT HIGHLIGHT
        ###############################
        
        ##If we have just 1 row in domainsCoord, it implies
        ##That everyting occurs in the same TAD,
        ##Hence in that case, we still have two breakpoints
        ##So add the highlight for both breakpoints in the same TAD
        
        if(nrow(domainsCoord)>1){
          
          ##So we have the breakpoints in 2 different tracks, hence we will create 2 tracks in the output html
          
          ##Let's retrieve breakpoint ID to get subsequently its coord
          breakpoint_ID<-unlist(strsplit(x = rownames(domainsCoord)[nDomain], split = "domain_",fixed = TRUE))[2]
          
          ##Now let's get both breakpointLimits
          breakpoint_start<-breakpCoord[paste(breakpoint_ID,"_left_segment",collapse = "",sep = "")
                                        ,"end"]
          
          breakpoint_end<-breakpCoord[paste(breakpoint_ID,"_right_segment",collapse = "",sep = "")
                                      ,"start"]
          
          ##########
          ## Extending coord to highlight a bigger region and facilitate visualization
          breakp_HighlightStart<-breakpoint_start - breakpointHighlight_extension
          breakp_HighlightEnd<-breakpoint_end + breakpointHighlight_extension
          
          breakP_region<-paste("hg19.",
                               chr_browser,"%3A",
                               breakp_HighlightStart,"-",
                               breakp_HighlightEnd,
                               collapse = "",
                               sep = "")
          
          ##RegionColor For the Regions ##BLUE
          colorBreakp<-paste("%23",##%23 is essential is the conversion somehow to url syntax of #
                             "ff3300",
                             sep = "",
                             collapse = "")
          
          ##Insterting breakp highlight on track
          browser_highlight<-paste(browser_highlight,
                                   "%7C",
                                   breakP_region,
                                   colorBreakp,
                                   sep="",
                                   collapse = "")
        }else{
          if(nrow(domainsCoord)==1){
            ##It implies that everything occur in the same TAD , so we need to put both breakpoints highlights in the same track
            
            for(breakpoint_ID in c("breakpoint_1","breakpoint_2")){
              ##Now let's get both breakpointLimits
              breakpoint_start<-breakpCoord[paste(breakpoint_ID,"_left_segment",collapse = "",sep = "")
                                            ,"end"]
              
              breakpoint_end<-breakpCoord[paste(breakpoint_ID,"_right_segment",collapse = "",sep = "")
                                          ,"start"]
              
              ##########
              ## Extending coord to highlight a bigger region and facilitate visualization
              breakp_HighlightStart<-breakpoint_start - breakpointHighlight_extension
              breakp_HighlightEnd<-breakpoint_end + breakpointHighlight_extension
              
              breakP_region<-paste("hg19.",
                                   chr_browser,"%3A",
                                   breakp_HighlightStart,"-",
                                   breakp_HighlightEnd,
                                   collapse = "",
                                   sep = "")
              
              ##RegionColor For the Regions ##BLUE
              colorBreakp<-paste("%23",##%23 is essential is the conversion somehow to url syntax of #
                                 "ff3300",
                                 sep = "",
                                 collapse = "")
              
              ##Insterting breakp highlight on track
              browser_highlight<-paste(browser_highlight,
                                       "%7C",
                                       breakP_region,
                                       colorBreakp,
                                       sep="",
                                       collapse = "")
              
            }
          }
        }
      }
      
      ##########################################
      ## Adding highligted regions to UCSC link
      ##########################################
      
      linkString<-paste(linkString,
                        browser_highlight,
                        sep = "",
                        collapse = "")
      
      ##Affected Domain Nr
      NrDomain<-paste("domain_",
                      nDomain,
                      sep = "",
                      collapse = "")
      
      ucsc_domain_tracks[NrDomain]<-linkString
      
      ################################################
      ##Let's add the linkString to the shiny html
      ################################################
      
      ##If currently handling the target gene, add info to link
      
      if(coveringCurrentlyTargetGene==TRUE){
        shiny_html_ucsc<-paste(shiny_html_ucsc,
                               "<p></p><p style='font-size:22px'><a href='",linkString, "' target='_blank'>", "<b>Visualize</b>" ," Affected Regulatory Domain: ",nDomain,
                               " (",
                               targetGene,
                               ")",
                               "</a></p>", 
                               sep = "",collapse = "")
      }else{
        ##As it was before, just indicating the number. There is only 1 link, hence no need to differentiate them
        shiny_html_ucsc<-paste(shiny_html_ucsc,
                               "<p></p><p style='font-size:22px'><a href='",linkString, "' target='_blank'>", "<b>Visualize</b>" ," Affected Regulatory Domain: ",nDomain,
                               "</a></p>", 
                               sep = "",collapse = "")
      }

    }
    
  }else if((sv_type == "Deletion") || (sv_type == "Duplication")){
    ##Here we paint everything for now in Just one link, as everything between breakpoints gets affected
    nDomain<-1
    
    ##beginningLink
    ##Here take the UCSC base name link, dependent on phenotype and Phase
    ## eg "https://genome.ucsc.edu/s/sgvic/SV_RadaLab"
    linkBaseName<-genomeBrowser_links$SessionLinkBaseName[genomeBrowser_links$SessionId==browserSessionId]
    
    linkString<-paste(linkBaseName,
                      "?db=hg19&position=",
                      sep="")   ##This is the require beggining, we are directing to 
    
    ##To display, extend TAD coordinates +- 10Kb
    ##start_browser, left most TAD coord
    ##end_browser, right most TAD coord
    chr_browser<-unique(domainsCoord$chr) ##Del or Dup, so only 1 chromosome, happening in the same one
    start_browser<-min(domainsCoord$start) - zoomOutDist
    end_browser<-max(domainsCoord$end) + zoomOutDist
    
    ##Formatting chr coordinates
    
    formated_browserView<-paste(chr_browser,":",
                                start_browser,"-", 
                                end_browser, sep = "")
    
    linkString<-paste(linkString,formated_browserView,sep = "" ,collapse = "")
    
    ##We do not want user cookies to be loaded that can interfere with our display
    linkString<-paste(linkString,"&ignoreCookie=1",sep = "" ,collapse = "")
    
    ##############################################
    ##Adding Higlights
    ##############################################
    browser_highlight<-"&highlight="
    
    ##############################################################################
    ##For TAD limits +- 100bp
    ##Aqui si que pintamos los breakp de los dos tads Disrupted (si >1 disrupted)
    ##############################################################################
    #For now not adding tad limits because we have the bed coord below
    #If we need to paint an average tad coord... maybe doing it by rethinking the track
    #Or adding beed coord simultaneously, or painting regarding the coord of the TAD providing a high score
    #Internally we take into account the average
    #But not painting TAD coord in bed + grey vertical lines because too much, emborrona mucho
    
    # for(nDomainAffected in 1:nrow(domainsCoord)){
    #   print(nDomainAffected)
    #   
    #   ####For each TAD limit, the colouring coords are TADlimit+-100bp
    #   start_TAD<-domainsCoord[nDomainAffected,"start"]
    #   extension_start<-start_TAD - borderHighlight_extension
    #   
    #   end_TAD<-domainsCoord[nDomainAffected,"end"]
    #   extension_end<-end_TAD + borderHighlight_extension
    #   
    #   #startRegion
    #   startRegion<-paste("hg19.",
    #                      chr_browser,"%3A",
    #                      extension_start,"-",
    #                      start_TAD,
    #                      collapse = "",
    #                      sep = "")
    #   #EndRegion
    #   endRegion<-paste("hg19.",
    #                    chr_browser,"%3A",
    #                    end_TAD,"-",
    #                    extension_end,
    #                    collapse = "",
    #                    sep = "")
    #   
    #   ##RegionColor For the Regions ##gray
    #   colorBorders<-paste("%23",##%23 is essential is the conversion somehow to url syntax of #
    #                       "877878",
    #                       sep = "",
    #                       collapse = "")##colorCode  # el del breakpoint "#e60000"
    #   
    #   browser_highlight<-paste(browser_highlight,
    #                            startRegion,
    #                            colorBorders,
    #                            "%7C",
    #                            endRegion,
    #                            colorBorders,
    #                            sep="",
    #                            collapse = "")
    # }
    
    
    ###############################
    ## ADDING BREAKPOINT HIGHLIGHT
    ###############################
    
    ##Let's retrieve breakpoint ID to get subsequently its coord
    # breakpoint_ID<-unlist(strsplit(x = rownames(domainsCoord)[nDomain], split = "domain_",fixed = TRUE))[2]
    breakpoint_ID<-"breakpoint_1&2"
    
    ##Now let's get both breakpointLimits
    breakpoint_start<-min(as.numeric(unlist(strsplit(x = patientResults$patientInfo$coord_Break1,
                              split = ",", fixed = TRUE))))
    
    breakpoint_end<-max(as.numeric(unlist(strsplit(x = patientResults$patientInfo$coord_Break2,
                                                   split = ",", fixed = TRUE))))
    
    ##########
    ## Extending coord to highlight a bigger region and facilitate visualization
    ## I do not think this is necessary for Del Dup
    ########## 
    breakp_HighlightStart<-breakpoint_start #- breakpointHighlight_extension
    breakp_HighlightEnd<-breakpoint_end #+ breakpointHighlight_extension
    
    breakP_region<-paste("hg19.",
                         chr_browser,"%3A",
                         breakp_HighlightStart,"-",
                         breakp_HighlightEnd,
                         collapse = "",
                         sep = "")
    
    ##RegionColor For the Regions ##BLUE
    colorBreakp<-paste("%23",##%23 is essential is the conversion somehow to url syntax of #
                       "ff3300",
                       sep = "",
                       collapse = "")
    
    ##Insterting breakp highlight on track
    browser_highlight<-paste(browser_highlight,
                             "%7C",
                             breakP_region,
                             colorBreakp,
                             sep="",
                             collapse = "")
    
    #############################################################
    ## Highligthing All enhancers in the affected genome region
    #############################################################
    
    masterEnh<-patientResults$MasterEnh_map
    masterEnh<-subset(masterEnh, source==devStage)
    #subset for genome region
    masterEnh$isChr<-masterEnh$chr==chr_browser
    masterEnh$isBeforeEnd<-masterEnh$end<=end_browser
    masterEnh$isAfterStart<-masterEnh$start>=start_browser
    
    targetEnh<-rowSums(masterEnh[,c("isChr","isBeforeEnd","isAfterStart")])==3
    ##filter for targetEnh
    masterEnh<-masterEnh[targetEnh,]
    
    if(nrow(masterEnh)>0){
      #so at least one enh to highlight
      for(nEnh in 1:nrow(masterEnh)){
        #print(nEnh)
        #enhChr we know, is the chr_browser
        enh_start<-masterEnh[nEnh,"start"]
        enh_end<-masterEnh[nEnh,"end"]
        
        
        ##Add ucsc info to highlight
        enh_region<-paste("hg19.",
                          chr_browser,"%3A",
                          enh_start,"-",
                          enh_end,
                          collapse = "",
                          sep = "")
        
        ##RegionColor For the Regions #b3ffb3 33cc33 85e085
        colorEnh<-paste("%23",##%23 is essential is the conversion somehow to url syntax of #
                        "85e085",
                        sep = "",
                        collapse = "")
        
        ##Insterting breakp highlight on track
        browser_highlight<-paste(browser_highlight,
                                 "%7C",
                                 enh_region,
                                 colorEnh,
                                 sep="",
                                 collapse = "")
      }
    }
    
    
    #############################################################
    ## Highligthing The Target Gene if its domain is being linked
    #############################################################
    genesInfo<-patientResults$allAffectedGenes_positionalInfo
    genesInfo<-genesInfo[targetGene,]
    gene_chr<-genesInfo$chr
    gene_start<-genesInfo$start
    gene_end<-genesInfo$end
    gene_tss<-genesInfo$TSS
    
    ##if currently the genome browser is covering the target gene
    ##then highlight it
    if((gene_chr==chr_browser) && (gene_tss>=start_browser) && (gene_tss<=end_browser)){
      ##highlight the gene body, as it is included in the currently created link area
      
      ##Add ucsc info to highlight
      gene_region<-paste("hg19.",
                         chr_browser,"%3A",
                         gene_start,"-",
                         gene_end,
                         collapse = "",
                         sep = "")
      
      ##RegionColor For the Regions #3399ff
      colorGene<-paste("%23",##%23 is essential is the conversion somehow to url syntax of #
                       "3399ff",
                       sep = "",
                       collapse = "")
      
      ##Insterting breakp highlight on track
      browser_highlight<-paste(browser_highlight,
                               "%7C",
                               gene_region,
                               colorGene,
                               sep="",
                               collapse = "")
    }
    
    ########################################################################
    ### IF WE ARE DEALING WITH TRUNCATION
    ## So, the gene partially altered
    ## PAINT AGAIN THE BREAKPOINT TO AVOID BEING HIDDEN BY THE GENE BLUEISH
    ########################################################################
    
    if(patientResults$masterSummaryResultsMatrix[targetGene,"GeneImpact"] == "Direct_geneTruncation"){
      ###############################
      ## ADDING BREAKPOINT HIGHLIGHT
      ###############################
      
      ##If we have just 1 row in domainsCoord, it implies
      ##That everyting occurs in the same TAD,
      ##Hence in that case, we still have two breakpoints
      ##So add the highlight for both breakpoints in the same TAD
      
      if(nrow(domainsCoord)>1){
        
        ##So we have the breakpoints in 2 different tracks, hence we will create 2 tracks in the output html
        
        ##Let's retrieve breakpoint ID to get subsequently its coord
        breakpoint_ID<-unlist(strsplit(x = rownames(domainsCoord)[nDomain], split = "domain_",fixed = TRUE))[2]
        
        ##Now let's get both breakpointLimits
        breakpoint_start<-breakpCoord[paste(breakpoint_ID,"_left_segment",collapse = "",sep = "")
                                      ,"end"]
        
        breakpoint_end<-breakpCoord[paste(breakpoint_ID,"_right_segment",collapse = "",sep = "")
                                    ,"start"]
        
        ##########
        ## Extending coord to highlight a bigger region and facilitate visualization
        breakp_HighlightStart<-breakpoint_start - breakpointHighlight_extension
        breakp_HighlightEnd<-breakpoint_end + breakpointHighlight_extension
        
        breakP_region<-paste("hg19.",
                             chr_browser,"%3A",
                             breakp_HighlightStart,"-",
                             breakp_HighlightEnd,
                             collapse = "",
                             sep = "")
        
        ##RegionColor For the Regions ##BLUE
        colorBreakp<-paste("%23",##%23 is essential is the conversion somehow to url syntax of #
                           "ff3300",
                           sep = "",
                           collapse = "")
        
        ##Insterting breakp highlight on track
        browser_highlight<-paste(browser_highlight,
                                 "%7C",
                                 breakP_region,
                                 colorBreakp,
                                 sep="",
                                 collapse = "")
      }else{
        if(nrow(domainsCoord)==1){
          ##It implies that everything occur in the same TAD , so we need to put both breakpoints highlights in the same track
          
          for(breakpoint_ID in c("breakpoint_1","breakpoint_2")){
            ##Now let's get both breakpointLimits
            breakpoint_start<-breakpCoord[paste(breakpoint_ID,"_left_segment",collapse = "",sep = "")
                                          ,"end"]
            
            breakpoint_end<-breakpCoord[paste(breakpoint_ID,"_right_segment",collapse = "",sep = "")
                                        ,"start"]
            
            ##########
            ## Extending coord to highlight a bigger region and facilitate visualization
            breakp_HighlightStart<-breakpoint_start - breakpointHighlight_extension
            breakp_HighlightEnd<-breakpoint_end + breakpointHighlight_extension
            
            breakP_region<-paste("hg19.",
                                 chr_browser,"%3A",
                                 breakp_HighlightStart,"-",
                                 breakp_HighlightEnd,
                                 collapse = "",
                                 sep = "")
            
            ##RegionColor For the Regions ##BLUE
            colorBreakp<-paste("%23",##%23 is essential is the conversion somehow to url syntax of #
                               "ff3300",
                               sep = "",
                               collapse = "")
            
            ##Insterting breakp highlight on track
            browser_highlight<-paste(browser_highlight,
                                     "%7C",
                                     breakP_region,
                                     colorBreakp,
                                     sep="",
                                     collapse = "")
            
          }
        }
      }
    }
    
    ##########################################
    ## Adding highligted regions to UCSC link
    ##########################################
    
    linkString<-paste(linkString,
                      browser_highlight,
                      sep = "",
                      collapse = "")
    
    ##Affected Domain Nr
    NrDomain<-paste("domain_",
                    nDomain,
                    sep = "",
                    collapse = "")
    
    ucsc_domain_tracks[NrDomain]<-linkString
    
    #############
    ##Let's add the linkString to the shiny html
    shiny_html_ucsc<-paste(shiny_html_ucsc,
                           "<p></p><p style='font-size:22px'><a href='",linkString, "' target='_blank'>", "<b>Visualize</b>" ," Affected Regulatory Domain: ",nDomain,"</a></p>", 
                           sep = "",collapse = "")
    
  }
  
  
  
  ################################################################
  ###Lets add the UCSC logo on the  left part of the string
  ###and the ucsc links info on the right
  ##to get info about the css check: "MainInterfaceStyling.html"
  ################################################################
  info_imageLogo<-"<div class='logoUCSC_location'><img src='UCSC_logo.png' alt='UCSC logo' style='width:100%;'></div>"

  shiny_html_ucsc<-paste("<div class='linksUCSC_location'>",
                         shiny_html_ucsc,
                         "</div>",
                         sep=""
                         )
  
  shiny_html_ucsc<-paste("<div class='infoLogoLink'>",
                         info_imageLogo,
                         shiny_html_ucsc,
                         "</div>",
                         sep="")
  
  
  ##############################################################
  ### Let's add the color code Legend
  ##############################################################

  info_colorCode<-paste("<div class='colorCodeSection'>
	  <div class='cs_colorRegDom'></div>
	  <div class='cs_colorBreak'></div>
	  <div class='cs_colorGene'></div>
	  <div class='cs_colorEnh'></div>
	  <div class='cs_textRegDom'><p>TADs - Regulatory Domains</p></div>
    <div class='cs_textBreak'><p>",
    if(patientResults$patientInfo$TypeSV=="Deletion"){
        "Deleted Region"
    }else if(patientResults$patientInfo$TypeSV=="Duplication"){
        "Duplicated Region"
    }else{
        ##Hence Inversion or Translocation
        ##Here the reddish region represents solely the breakpoints
        "Breakpoint"
    },
    "</p></div>
	  <div class='cs_textGene'><p>",
    targetGene,               
    "</p></div>
	  <div class='cs_textEnh'><p>Enhancers</p></div>
  </div>",
                        sep="")
  
  ##removing line breaks
  info_colorCode<-gsub("[\r\n]", "", info_colorCode)
  
  ##legend Title
  info_legendTitle<-"<div class='legendHeader'>
  <p style='font-weight:bold;font-size:22.5px'>Legend for colors used in Browser</p>
  </div>"
  
  info_legendTitle<-gsub("[\r\n]", "", info_legendTitle)
  
  ###Assembling info legend color code
  macroInfo_colorCode<-""
  macroInfo_colorCode<-paste("<div class='wrapperColorLegend'>",
                             info_legendTitle,
                             info_colorCode,
                             "</div>",
                             sep="")
  
  #adding colorCode info to html
  shiny_html_ucsc<-paste(shiny_html_ucsc,
                         macroInfo_colorCode,
                         sep="")
  
  
  ########
  #Wrapping all info for the ucsc links section
  ########
  shiny_html_ucsc<-paste("<div class='ucsc_infoSection'>",
                         shiny_html_ucsc,
                         "</div>")
  
  return(shiny_html_ucsc)
}







