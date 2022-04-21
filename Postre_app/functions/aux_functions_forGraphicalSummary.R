######################################################################################
## Auxiliar Functions developed for Graphical Summary plotting
## Maybe this functions are only valid for Inversions And Translocations Between TADs
######################################################################################

##NOTE: To place text labels, take as reference the center point of interest,
##since the label is centered over that coordinate

paintGene_WT_TAD<-function(tad_X_cord, tad_Y_cord, nEnh_initial_left, nEnh_initial_right, gene, gene_breakp_line_type, situation, patientResults){
  ########################################################
  ## Function to paint GENE WT TAD Initial representation
  ########################################################
  
  
  enh_y_positions<-rep.int(x=tad_Y_cord[1], times = 4)
  
  ## Let's mantain color codes by TAD position
  ## If TAD is on the left blue or whatever we choos 1, if it is on the right orange or whatever we choose 2
  ## Regardless whether it is the gene or the secondary domain TAD
  ## In order not to confound the user if > 1 gene is affected
  ##For central tads... all of them paint in kind of orange... not to confound...consider for future
  if(situation == "primaryTAD_Central"){
    colorTAD<-"#9999ff"
  }else if(tad_X_cord[1]<20){
    colorTAD<-"#acdfeb" ##originally blueish
  }else if(tad_X_cord[1]>20){
    colorTAD<-"#e5edb2" ##oranginally orangeish
  }
  
  
  ##Painting TAD
  polygon(tad_X_cord, tad_Y_cord, col = colorTAD, border = "#ffffff")
  
  ##Defining some gene features required to plot it
  #And for other calculations of breakpoint and enhs locations
  x_space<-c(tad_X_cord[1],tad_X_cord[2])
  
  x_space_start<-x_space[1]
  x_space_end<-x_space[2]
  
  #gene position at the center -0.5 (since we provide the left coord over which the gene is painted)
  genePos<-(x_space_start + (x_space_end-x_space_start)/2)-0.5##o coger el centro del TAD con TAD3 y restarle 0.5 ya que mide 1
  
  #genePos<-tad_X_cord[1]+4.5
  geneSize<-1
  geneCenter<-genePos + geneSize/2
  gene_X_cord<-c(genePos, genePos+geneSize, genePos+geneSize, genePos)
  gene_Y_cord<-c(tad_Y_cord[1]-1,tad_Y_cord[1]-1,tad_Y_cord[1]+1,tad_Y_cord[1]+1)
  
  
  #########################################
  ## Defining coordinates breakpoint line
  ## And additional relevant variables
  #########################################
  enhNumbersSize<-0.8
  
  distanceBreakpFromGene<-0.3 ##Distance breakpoint from the gene, when the breakpoint is between the gene & enh
  
  distanceEnhFromGene<-1.5 ##Distance between the enhancers and the gene
  
  distanceBreakpFromEnh<-0.3 ##Distance breakpoint from enhancers when the breakpoint located between enh and TAD border
  
  sizeEnhRegion<-1##Touching this will not really affect enh size, because created regardless of this
  
  distance_Yaxis_geneLabel<-2
  distance_Yaxis_EnhGene<-2
  distance_Yaxis_EnhLabel_EnhNumber<-1.3
  
  ##############
  ## I think I do not need the primary TAD info situation
  
  ##Important to have in mind the TAD info situation if it is tad central... we paint two breakpoints there
  
  if(situation %in% c("primaryTAD_Dextral","primaryTAD_Sinistral")){
    ###################################
    ##TWO TADS ARE GOING TO BE PAINTED, so in the current TAD representation
    ## only ONE BREAKPOINT will be displayed
    ###################################
    
    ##combining with duplic but attending to the logic behind removing or duplicating to place the lines
    ##(backup prior to this merges from 26March)
    
    if(gene_breakp_line_type=="center"){
      ##It implies the gene has been broken so the breakpoint located over gene
      #No need to consider enhancer changes to refine breakpoint location
      breakPos<-geneCenter
      
    }else if((gene_breakp_line_type=="afterTSS_removing_all") ||
             (gene_breakp_line_type=="afterTSS_duplic_none")){
      ## The breakpoint is located right after the end of the gene
      
      breakPos<-gene_X_cord[2] + distanceBreakpFromGene
      
    }else if((gene_breakp_line_type=="afterTSS_removing_some") ||
             (gene_breakp_line_type=="afterTSS_duplic_some")){
      
      enhXpos<-genePos + distanceEnhFromGene##Aqui va el enh distance from gene ##Obtained from painting enhancers right side code part. Enh cluster Start
      breakPos<-enhXpos+0.5
      
    }else if((gene_breakp_line_type=="afterTSS_removing_none") ||
             (gene_breakp_line_type=="afterTSS_duplic_all")){
      
      ## Breakpoint before TAD end. Between enh cluster end TAD end
      ##just after enhancers
      breakPos<-genePos+distanceEnhFromGene+sizeEnhRegion+distanceBreakpFromEnh
      
    }else if((gene_breakp_line_type=="beforeTSS_removing_all") ||
             (gene_breakp_line_type=="beforeTSS_duplic_none")){
      ## The breakpoint is located right before the begining of the gene
      breakPos<-genePos-distanceBreakpFromGene
      
    }else if((gene_breakp_line_type=="beforeTSS_removing_some") ||
             (gene_breakp_line_type=="beforeTSS_duplic_some")){
      
      enhXpos<-genePos-distanceEnhFromGene##Obtained from painting enhancers left side code part. Enh cluster Start
      breakPos<-enhXpos+0.5
      
    }else if((gene_breakp_line_type=="beforeTSS_removing_none") ||
             (gene_breakp_line_type=="beforeTSS_duplic_all")){
      ## Breakpoint after TAD start
      #breakPos<-tad_X_cord[1]+1.5 ##1.5 Units after TAD start
      ##just before enhancers
      breakPos<-genePos-distanceEnhFromGene-distanceBreakpFromEnh
      
    }
    
  }else if(situation == "primaryTAD_Central"){
    ###########################################################################
    ##   ONE TAD IS GOING TO BE PAINTED, so in the current TAD representation
    ##   TWO BREAKPOINTS will be displayed
    ##   It is due to the fact that both ,or none, breakpoint directly touch the TAD
    ## This is going to occur for Deletions And Duplications
    ###########################################################################
    
    if(gene_breakp_line_type=="outOf_RegulatoryDomain"){
      ##Breakpos will be painted out of the regulatory domain... and we will paint two
      breakPos<-c(tad_X_cord[1]-2, tad_X_cord[2]+2)
      
    }else if(gene_breakp_line_type=="center"){
      #######################
      ##So GENE TRUNCATED.
      #######################
      
      ###Check if both breakpoints intragenic
      ##or one intragenic and the other outside of the gene
      
      ##I'm just going to take the start position of the breakpoints for this purpose
      #Chr is going to be the same, since this scenario where we are now is occurring only for deletions or duplications
      #(scenario of painting just one TAD)
      bp1_start<-as.numeric(unlist(strsplit(x = patientResults$patientInfo$coord_Break1,
                          split = ",", fixed = TRUE))[1])
      
      bp2_start<-as.numeric(unlist(strsplit(x = patientResults$patientInfo$coord_Break2,
                          split = ",", fixed = TRUE))[1])
      
      gene_start<-patientResults$allAffectedGenes_positionalInfo[gene,"start"]
      
      gene_end<-patientResults$allAffectedGenes_positionalInfo[gene,"end"]
      
      ##Recall, bp1 smaller than bp2 for SV happening on the same chr (all except Translocation)
      if((bp1_start >= gene_start) && (bp2_start <= gene_end)){
        
        ##So both breakpoints painted inside of the gene
        breakPos<-c(geneCenter - 0.3,
                    geneCenter + 0.3)
        
      }else if((bp1_start >= gene_start) && (bp2_start > gene_end)){
        ##So one breakpoint intraGene and the other beyond
        breakPos<-c(geneCenter,
                    genePos + geneSize + 0.5)
        
      }else if((bp1_start < gene_start) && (bp2_start <= gene_end)){
        ##So one breakpoint before and the other intragene
        breakPos<-c(genePos - 0.5,
                    geneCenter)
      }
      
      ##It implies the gene has been broken so the breakpoint located over gene
      #No need to consider enhancer changes to refine breakpoint location
      #breakPos<-geneCenter
      
    }else if(gene_breakp_line_type=="surroundingGene"){
      ##One breakpoint before the gene
      ##One breakpoint after the gene
      breakPos<-c(genePos - 0.5,
                  genePos + geneSize + 0.5)
      
    }else if((gene_breakp_line_type=="afterTSS_removing_all") || ##so all enh deleted 
             (gene_breakp_line_type=="afterTSS_duplic_all"))##so all enh duplicated  
                ##So we want lines surrounding the enh cluster
                ##combining with duplic but attending to the logic behind removing or duplicating to place the lines
                ##(backup prior to this merges from 26March)
      {
      ## The breakpoint is located right after the end of the gene
      
      breakPos<-c(gene_X_cord[2] + distanceBreakpFromGene, ##Before Enh
                  genePos + distanceEnhFromGene + sizeEnhRegion + distanceBreakpFromEnh)##After Enh
      
    }else if((gene_breakp_line_type=="afterTSS_removing_some") ||
             (gene_breakp_line_type=="afterTSS_duplic_some")){
      
      enhXpos<-genePos + distanceEnhFromGene
      breakPos<-c(enhXpos + sizeEnhRegion/2,
                  enhXpos + sizeEnhRegion + distanceBreakpFromEnh)##One breakpoint inside of the cluster, the other beyond
                  
    }else if((gene_breakp_line_type=="afterTSS_removing_none") ||
             (gene_breakp_line_type=="afterTSS_duplic_none")){
      
      ##just after enhancers
      positionBreak<-genePos+distanceEnhFromGene+sizeEnhRegion+distanceBreakpFromEnh
      breakPos<-c(positionBreak,
                  positionBreak+1)##after enhcluster and beyond
      
    }else if((gene_breakp_line_type=="beforeTSS_removing_all") ||
             (gene_breakp_line_type=="beforeTSS_duplic_all")){
      # The breakpoint is located right before the begining of the gene
      # genePos-distanceEnhFromGene = start coordinate of the enh cluster
      breakPos<-c(genePos - distanceBreakpFromGene, ##between enh and cluster
                  genePos - distanceEnhFromGene - distanceBreakpFromEnh)
      
    }else if((gene_breakp_line_type=="beforeTSS_removing_some") ||
             (gene_breakp_line_type=="beforeTSS_duplic_some")){
      
      enhXpos<-genePos-distanceEnhFromGene##Obtained from painting enhancers left side code part. Enh cluster Start
      
      breakPos<-c(enhXpos-distanceBreakpFromEnh,
                  enhXpos+0.5)##One breakpoint before the cluster, and the other inside
      
      
    }else if((gene_breakp_line_type=="beforeTSS_removing_none") || 
             (gene_breakp_line_type=="beforeTSS_duplic_none")){
      ## Breakpoint after TAD start
      #breakPos<-tad_X_cord[1]+1.5 ##1.5 Units after TAD start
      ##just before enhancers
      positionBreak<-genePos-distanceEnhFromGene-distanceBreakpFromEnh
      
      breakPos<-c(positionBreak-1,
                  positionBreak)
      
    }
      
  }
  
  #########################
  ##Painting breakpoint/s
  #########################
  nIterations<-0##labels only painted on the last iteration
  for(targetBreakpoint in breakPos){
    #print(targetBreakpoint)
    nIterations<-nIterations + 1
    
    ##Defining Breakpoint position
    breakp_x_coord<-c(targetBreakpoint, targetBreakpoint)##Tad central point for gene broken
    breakp_y_coord<-c(tad_Y_cord[1]-7, tad_Y_cord[3]+3)
    
    ############################
    ##Painting Breakpoint LINES
    breakpointColor<-"brown3"
    lines(x=breakp_x_coord,
          y=breakp_y_coord,
          col=breakpointColor,
          lwd=2,
          lty=3)
    
    #################################
    ## Painting breakpoint Labels
    if(nIterations==length(breakPos)){
      ##So we are on the last iteration, let's paint the LABELS
      ############################
      ##Painting breakpoint labels
      #############################
      #If they are >1 & super close, paint just "Breakpoints", in the middle of the breakpoints space
      if(length(breakPos)>1){
        if(abs(breakPos[2]-breakPos[1])<5){
          
          ##Using min() because smaller coord can be on any pos
          breakLabelPosition<-(min(breakPos)+abs(breakPos[2]-breakPos[1])/2)  ##Midpoint breakpoint lines
          
          #Breakpoint label
          ##Boxed labels: xpad,ypad arguments, used to expand the box beyond the text
          boxed.labels(x=breakLabelPosition,  y=breakp_y_coord[1] - 1,labels = "Breakpoints", cex=0.8, border = NA, bg ="white", 
                       xpad=1,
                       ypad=2, ##To allow the text to breathe
                       col=breakpointColor
          )
          
        }else{
          ##Paint both labels
          for(targetBreakpoint in breakPos){
            ##Defining Breakpoint position
            breakp_y_coord<-c(tad_Y_cord[1]-7, tad_Y_cord[3]+3)
            
            boxed.labels(x=targetBreakpoint,  y=breakp_y_coord[1] - 1,labels = "Breakpoint", cex=0.8, border = NA, bg ="white", 
                         xpad=1,
                         ypad=2, ##To allow the text to breathe
                         col=breakpointColor
            )
          }
        }
      }else if(length(breakPos)==1){
        ## JUST ONE BREAKPOINT PAINTED. 
        #Breakpoint label
        ##Boxed labels: xpad,ypad arguments, used to expand the box beyond the text
        boxed.labels(x=targetBreakpoint,  y=breakp_y_coord[1] - 1,labels = "Breakpoint", cex=0.8, border = NA, bg ="white", 
                     xpad=1,
                     ypad=2, ##To allow the text to breathe
                     col=breakpointColor
        )
      }
    }
    
  }
  
  #TAD label
  # text(x=tad_X_cord[3], y=tad_Y_cord[1]+10, label="TAD", cex = 1.5)
  boxed.labels(x=tad_X_cord[3],  y=tad_Y_cord[1]+10,labels = "TAD", cex=1.5,
               border = NA, bg =NA, xpad=1, ypad = 2 )
  
  ###Adding black line below TAD
  lines(x = c(tad_X_cord[1]-3,tad_X_cord[2]+3),
        y=c(tad_Y_cord[1],tad_Y_cord[1]), 
        col="black", lwd=2, lty=1)
  
  
  ###################
  #Adding gene body
  polygon(gene_X_cord, gene_Y_cord, col = "#0e3d61")
  
  #gene Label
  ##I want it to overlay the breakpoint line
  # https://stackoverflow.com/questions/45366243/text-labels-with-background-colour-in-r
  boxed.labels(x=tad_X_cord[3],  y=tad_Y_cord[1]-distance_Yaxis_geneLabel,
               labels = gene, cex=0.8, border = NA, bg ="white", 
               xpad=1,
               ypad=2##To allow the text to breathe
               )
  # text(x=tad_X_cord[3], y=tad_Y_cord[1]-2, label=gene, cex = 0.8)
  
  ##Adding enhancers
  
  ##############################
  ####To the Right side gene Position
  if(nEnh_initial_right>0){
    ##So there is at least 1 enh
    enhXpos<-genePos+distanceEnhFromGene##left margin of the enh cluster
    enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
    
    draw.ellipse(x=enh_x_positions,
                 y=enh_y_positions,a=0.1,col="#b2d235")
    
    #enhancers Label
    #text(x=enhXpos + 0.4, y=enh_y_positions[1]-3, label="enhancers", cex = 0.8)
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel
                 -distance_Yaxis_EnhGene,
                 labels = "enhancers", cex=0.8, 
                 border = NA, bg ="white", 
                 xpad=1,
                 ypad=1 #To allow the text to breath
    )
    #Nr of enhancers
    #text(x=enhXpos + 0.4, y=enh_y_positions[1]-4, label=paste0("n=",nEnh_initial_right), cex=enhNumbersSize)
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel
                 -distance_Yaxis_EnhGene
                 -distance_Yaxis_EnhLabel_EnhNumber,
                 labels = paste0("n=",nEnh_initial_right), cex=enhNumbersSize, 
                 border = NA, bg ="white", 
                 xpad=1,
                 ypad=1 #To allow the text to breath
    )
    
    ##enhancers horizontal line over them
    lines(x = c(enh_x_positions[1]-0.1,
                enh_x_positions[length(enh_x_positions)]+0.1),
          y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
          lwd=1.5)
    
    curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                   ,enh_y_positions[1]+2),
                arr.pos = 0, arr.type="T", arr.col = "black")##before angle=40
  }

  #############################
  ####To the Left side gene Position
  if(nEnh_initial_left>0){
    enhXpos<-genePos-distanceEnhFromGene##left margin of the enh cluster
    enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
    
    draw.ellipse(x=enh_x_positions,
                 y=enh_y_positions,a=0.1,col="#b2d235")
    
    #enhancers Label
    #text(x=enhXpos + 0.4, y=enh_y_positions[1]-3, label="enhancers", cex = 0.8)
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel
                 -distance_Yaxis_EnhGene,
                 labels = "enhancers", cex=0.8, 
                 border = NA, bg ="white", 
                 xpad=1,
                 ypad=1 #To allow the text to breath
    )
    
    #Nr of enhancers
    #text(x=enhXpos + 0.4, y=enh_y_positions[1]-4, label=paste0("n=",nEnh_initial_left), cex=enhNumbersSize)
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel
                 -distance_Yaxis_EnhGene
                 -distance_Yaxis_EnhLabel_EnhNumber,
                 labels = paste0("n=",nEnh_initial_left), cex=enhNumbersSize, 
                 border = NA, bg ="white", 
                 xpad=1,
                 ypad=1 #To allow the text to breath
    )
    
    ##enhancers horizontal line over them
    lines(x = c(enh_x_positions[1]-0.1,
                enh_x_positions[length(enh_x_positions)]+0.1),
          y = c(enh_y_positions[1]+2,enh_y_positions[1]+2),
          lwd=1.5)
    
    ##arr type in T shape and in 0 pos to be placed at the beginning and hence hidden
    #we will color afterwards a triangle over the gene, to deal with an external function bug
    curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter #-0.25
                                                                   ,enh_y_positions[1]+2),
                arr.pos = 0, arr.type="T", curve = -1, arr.col = "black" )
    
  }
  
  ##################################################################################
  ##Adding Triangle Shape over gene, if it has any enh, to represent the arrow end
  ##################################################################################
  if((nEnh_initial_left>0)|| (nEnh_initial_right>0)){
    triangle_X_Coord<-c(geneCenter-0.3, geneCenter, geneCenter+0.3)
    heightOrizontalLineEnh<-enh_y_positions[1]+2
    triangle_Y_Coord<-c(heightOrizontalLineEnh+0.3, heightOrizontalLineEnh-0.2, heightOrizontalLineEnh+0.3)
    polygon(triangle_X_Coord, triangle_Y_Coord, col = "black", border = "black")
  }
  
  
  
  ##################################################
  ## Assembling relevant info from paintGene_WT_TAD 
  ##################################################
  info_drawing<-list("geneCenter"=geneCenter, ##geneCenterPosition is important to know row orientation of enhancer arrow secondary TAD
                     "geneTAD_breakP"=breakPos,
                     "genePos"=genePos)
  
  
  return(info_drawing)
  
}

paint_Enhancer_WT_Secondary_TAD<-function(tad_X_cord, tad_Y_cord, nEnh_other_domain, geneCenter,otherDomain_breakp_line_type, situation, geneBreakP_Position_respectToTSS, patientResults){
  
  ##geneCenter needs to be known to paint the arrow orientations, towards the left or right
  
  ####################################################
  ## Paint Secondary TAD, with only enhancers, if any
  ####################################################
  
    
  ##########################################
  ## Defining additional relevant variables
  ##########################################
  
  enhNumbersSize<-0.8
  distance_Yaxis_EnhLabel<-2
  distance_Yaxis_EnhLabel_EnhNumber<-1.3
  
  
  ## Let's mantain color codes by TAD position
  ## If TAD is on the left blue, if it is on the right orange
  ## Regardless whether it is the gene or the secondary domain TAD
  ## In order not to confound the user if > 1 gene is affected
  if(tad_X_cord[1]<20){
    colorTAD<-"#acdfeb"
  }else if(tad_X_cord[1]>20){
    colorTAD<-"#e5edb2"
  }
  
  ##Painting TAD
  polygon(tad_X_cord, tad_Y_cord, col = colorTAD, border = "#ffffff")
  
  #################################################################
  ## Computing enh location, required to positionate breakpoints
  ## Afterwards if no enh... the enh are not displayed, but we need to know where are they to positionate the breakp
  enhXpos<-tad_X_cord[3]-0.45
  enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
  enh_y_positions<-rep.int(x=tad_Y_cord[1], times = 4)
  
  #########################################
  ##Defining coordinates breakpoint line
  distanceBreakpFromEnhCluster<-0.3
  
  ##get sv_type
  ##browser()
  sv_type<-patientResults$patientInfo$TypeSV
  
  if(sv_type %in% c("Inversion","Translocation")){
    
    ##The logic is different for Deletions:
    ##if sth is not gained in a Deletion, it is in the space between breakpoints
    ##IN the relocation area for inv & transloc. And in inv&transloc, what is in this place is what is gained
    
    #The logic is also different for Duplications attending to neoTAD context, because, the gene will be placed in the so called relocationa area
    #For inv and transloc, and this will put the breakp on the contrary where it should be in terms of all|none enh gained

    
    if(situation=="primaryTAD_Sinistral"){
      ##So the SECONDARY TAD is painted on the Right Side of the Screen
      ##We need to know whether the gene has it breakpoint before or after the TSS
      ##because if is before. The gene will be relocated to the other TAD
      ##So, no enh gain means, that the breakpoint is after the enh cluster
      
      if(geneBreakP_Position_respectToTSS=="afterTSS"){
        #So the gene breakpoint position > gene TSS
        ##The gene is NOT relocated to the other TAD so, not gaining enhancers implies the breakpoint is before the enhancers
        gene_relocated<-FALSE
        if(otherDomain_breakp_line_type=="brings_none" ){
          ## Breakpoint  before enh location
          # breakPos<-tad_X_cord[1]+1.5 ##1.5 Units after TAD start
          breakPos<-enh_x_positions[1] - distanceBreakpFromEnhCluster ##0.5 units before enh cluster,to offer a lot of space for rearrangement
          
        }else if(otherDomain_breakp_line_type=="brings_some"){
          
          ##I will paint the breakp then in the middle of the TAD, where the cluster of enh is for now located
          breakPos<-tad_X_cord[3]
          
        }else if(otherDomain_breakp_line_type=="brings_all"){
          ## The breakpoint is located after the enhancer group 
          ## I will put it 1.5 units before the END of the TAD hence
          # breakPos<-tad_X_cord[2]-1.5 #1.5 units before the END of the TAD hence
          breakPos<-enh_x_positions[4] + distanceBreakpFromEnhCluster ##0.5 units after enh cluster, to offer a lot of space for rearrangement
        }
        
      }else if(geneBreakP_Position_respectToTSS=="beforeTSS"){
        ##The gene is RELOCATED to the other TAD so, not gaining enhancers implies the breakpoint is AFTER the enhancers
        gene_relocated<-TRUE
        if(otherDomain_breakp_line_type=="brings_none" ){
          # breakPos<-tad_X_cord[2]-1.5##breakpoints 1.5 units before TAD end to cover all enh
          breakPos<-enh_x_positions[4] + distanceBreakpFromEnhCluster
          
        }else if(otherDomain_breakp_line_type=="brings_some"){
          
          ##I will paint the breakp then in the middle of the TAD, where the cluster of enh is for now located
          breakPos<-tad_X_cord[3]
          
        }else if(otherDomain_breakp_line_type=="brings_all"){
          # breakPos<-tad_X_cord[1]+1.5 ##1.5 Units after TAD start
          breakPos<-enh_x_positions[1] - distanceBreakpFromEnhCluster
        }
        
      }
    }else if(situation=="primaryTAD_Dextral"){
      ##So the SECONDARY TAD is painted on the Left Side of the Screen
      
      if(geneBreakP_Position_respectToTSS=="beforeTSS"){
        ##So the gene is NOT relocated, and not gaining enh implies that the breakpoint of the secondary domain is before the enhancers
        gene_relocated<-FALSE
        if(otherDomain_breakp_line_type=="brings_none" ){
          ## The breakpoint is located after the enhancer group 
          ## I will put it 1.5 units before the END of the TAD hence
          # breakPos<-tad_X_cord[2]-1.5 #1.5 units before the END of the TAD hence
          breakPos<-enh_x_positions[4] + distanceBreakpFromEnhCluster
        }else if(otherDomain_breakp_line_type=="brings_some"){
          
          ##I will paint the breakp then in the middle of the TAD, where the cluster of enh is for now located
          breakPos<-tad_X_cord[3]
          
        }else if(otherDomain_breakp_line_type=="brings_all"){
          ## Breakpoint  before enh location
          # breakPos<-tad_X_cord[1]+1.5 ##1.5 Units after TAD start
          breakPos<-enh_x_positions[1] - distanceBreakpFromEnhCluster
        }      
        
      }else if(geneBreakP_Position_respectToTSS=="afterTSS"){
        ##So the gene is RELOCATED to the other TAD. Hence not gaining enhancers implies that the secondary breakpoint is before the enhancers
        gene_relocated<-TRUE
        if(otherDomain_breakp_line_type=="brings_none" ){
          # breakPos<-tad_X_cord[1]+1.5 ##1.5 Units after TAD start
          breakPos<-enh_x_positions[1] - distanceBreakpFromEnhCluster
          
        }else if(otherDomain_breakp_line_type=="brings_some"){
          
          ##I will paint the breakp then in the middle of the TAD, where the cluster of enh is for now located
          breakPos<-tad_X_cord[3]
          
        }else if(otherDomain_breakp_line_type=="brings_all"){
          # breakPos<-tad_X_cord[2]-1.5 #1.5 units before the END of the TAD hence
          breakPos<-enh_x_positions[4] + distanceBreakpFromEnhCluster
        }      
      }
    }
  }else if(sv_type=="Deletion"){
    ##Here different logic, if sth is not gained in a Deletion, it is in the space between breakpoints
    ##IN the relocation area for inv & transloc. And in inv&transloc, what is in this place is what is gained
    
    ##Recall geneRelocation Terms come from Inv&Transloc, make no sense for Deletion, but leaving them for now just for the sake of the plot
    
    if(situation=="primaryTAD_Sinistral"){
      
      if(geneBreakP_Position_respectToTSS=="afterTSS"){
        #So the gene breakpoint position > gene TSS
        gene_relocated<-FALSE
        if(otherDomain_breakp_line_type=="brings_none" ){
          ## The breakpoint is located after the enhancer group 
          breakPos<-enh_x_positions[4] + distanceBreakpFromEnhCluster ##0.5 units after enh cluster, to offer a lot of space for rearrangement
      
        }else if(otherDomain_breakp_line_type=="brings_some"){
          
          ##I will paint the breakp then in the middle of the TAD, where the cluster of enh is for now located
          breakPos<-tad_X_cord[3]
          
        }else if(otherDomain_breakp_line_type=="brings_all"){
          ## Breakpoint  before enh location
          breakPos<-enh_x_positions[1] - distanceBreakpFromEnhCluster ##0.5 units before enh cluster,to offer a lot of space for rearrangement
          
        }
        
      }else if(geneBreakP_Position_respectToTSS=="beforeTSS"){
        ##The gene is RELOCATED to the other TAD so, not gaining enhancers implies the breakpoint is AFTER the enhancers
        gene_relocated<-TRUE ##Not make sense in this context
        ##For Deletions there is no gene relocation, THIS WILL IMPLY GENE DELETION!
        ##But let's keep it for the sake of the plot
        
        ##I'm not going too think too much about this, makes not much sense
        ##And in addition enh will not be painted (because gene direct effect)
        ##So just leave it the way it is 
        
        if(otherDomain_breakp_line_type=="brings_none" ){
          breakPos<-enh_x_positions[4] + distanceBreakpFromEnhCluster
        }else if(otherDomain_breakp_line_type=="brings_some"){
          ##I will paint the breakp then in the middle of the TAD, where the cluster of enh is for now located
          breakPos<-tad_X_cord[3]
        }else if(otherDomain_breakp_line_type=="brings_all"){
          breakPos<-enh_x_positions[1] - distanceBreakpFromEnhCluster
        }
        
      }
    }else if(situation=="primaryTAD_Dextral"){
      ##So the SECONDARY TAD is painted on the Left Side of the Screen
      
      if(geneBreakP_Position_respectToTSS=="beforeTSS"){
        gene_relocated<-FALSE
        
        if(otherDomain_breakp_line_type=="brings_none" ){
          ## Breakpoint  before enh location
          breakPos<-enh_x_positions[1] - distanceBreakpFromEnhCluster
          
        }else if(otherDomain_breakp_line_type=="brings_some"){
          
          ##I will paint the breakp then in the middle of the TAD, where the cluster of enh is for now located
          breakPos<-tad_X_cord[3]
          
        }else if(otherDomain_breakp_line_type=="brings_all"){
          ## The breakpoint is located after the enhancer group 
          breakPos<-enh_x_positions[4] + distanceBreakpFromEnhCluster
        }      
        
      }else if(geneBreakP_Position_respectToTSS=="afterTSS"){
        ##So the gene is RELOCATED to the other TAD. Hence not gaining enhancers implies that the secondary breakpoint is before the enhancers
        gene_relocated<-TRUE ##Not make sense in this context
        ##For Deletions there is no gene relocation, THIS WILL IMPLY GENE DELETION!!
        ##But let's keep it for the sake of the plot
        
        ##I'm not going too think too much about this, makes not much sense
        ##And in addition enh will not be painted (because gene direct effect)
        ##So just leave it the way it is 
        
        if(otherDomain_breakp_line_type=="brings_none" ){
          breakPos<-enh_x_positions[1] - distanceBreakpFromEnhCluster
        }else if(otherDomain_breakp_line_type=="brings_some"){
          ##I will paint the breakp then in the middle of the TAD, where the cluster of enh is for now located
          breakPos<-tad_X_cord[3]
        }else if(otherDomain_breakp_line_type=="brings_all"){
          breakPos<-enh_x_positions[4] + distanceBreakpFromEnhCluster
        }      
      }
    }
  }else if(sv_type=="Duplication"){
    gene_relocated<-FALSE ##Adding it here to avoid errors downstream, requiring this variable for inv & translocations
    
    ##If painting this TAD is because the duplication is between TADs hence pathogenic NEO tad scenario
    
    #The logic is also different for Duplications attending to neoTAD context, because, the gene will be placed in the so called relocationa area
    #For inv and transloc, and this will put the breakp on the contrary where it should be in terms of all|none enh gained
    
    if(situation=="primaryTAD_Sinistral"){
      ##So the SECONDARY TAD is painted on the Right Side of the Screen

      if(geneBreakP_Position_respectToTSS=="beforeTSS"){
        ##So the gene is duplicated
        #If it is after, is not duplicated, and if not duplicated, it will never contact ectopic enh from another TAD through the Duplication
        
        if(otherDomain_breakp_line_type=="brings_none" ){
          # breakPos<-tad_X_cord[1]+1.5 ##1.5 Units after TAD start
          breakPos<-enh_x_positions[1] - distanceBreakpFromEnhCluster
          
        }else if(otherDomain_breakp_line_type=="brings_some"){
          
          ##I will paint the breakp then in the middle of the TAD, where the cluster of enh is for now located
          breakPos<-tad_X_cord[3]
          
        }else if(otherDomain_breakp_line_type=="brings_all"){
          # breakPos<-tad_X_cord[2]-1.5##breakpoints 1.5 units before TAD end to cover all enh
          breakPos<-enh_x_positions[4] + distanceBreakpFromEnhCluster
        }
        
      }else if(geneBreakP_Position_respectToTSS=="afterTSS"){
        
        ###THIS WILL IMPLY THAT THE GENE IS NOT DUPLICATED
        ##SO THE ONLY POSSIBLE THING, in this context where two tads involved
        ## IS THAT THIS IS PATHOGENIC IS BECAUSE THE GENE IS TRUNCATED BY THE DUPLICATION
        ##Because if not everything duplicated, enhancer, go to a neoTAD so has no effect on the non duplicated gene
        #Unless as said is truncated
        
        ##I'm not going too think too much about this, makes not much sense
        ##And in addition enh will not be painted (because gene direct effect)
        ##So just leave it the way it is 
        
        ##Just going to say, paint the breakpoint in the middle of the TAD
        ##I will paint the breakp then in the middle of the TAD, where the cluster of enh is for now located
        breakPos<-tad_X_cord[3]
        
      }
    }else if(situation=="primaryTAD_Dextral"){
      ##So the SECONDARY TAD is painted on the Left Side of the Screen
      
      if(geneBreakP_Position_respectToTSS=="afterTSS"){
        ##So the gene is duplicated
        #If it is before, is not duplicated, and if not duplicated, it will never contact ectopic enh from another TAD through the Duplication
        
        if(otherDomain_breakp_line_type=="brings_none" ){
          # breakPos<-tad_X_cord[2]-1.5 #1.5 units before the END of the TAD hence
          breakPos<-enh_x_positions[4] + distanceBreakpFromEnhCluster
          
        }else if(otherDomain_breakp_line_type=="brings_some"){
          
          ##I will paint the breakp then in the middle of the TAD, where the cluster of enh is for now located
          breakPos<-tad_X_cord[3]
          
        }else if(otherDomain_breakp_line_type=="brings_all"){
          # breakPos<-tad_X_cord[1]+1.5 ##1.5 Units after TAD start
          breakPos<-enh_x_positions[1] - distanceBreakpFromEnhCluster
        }      
      }else if(geneBreakP_Position_respectToTSS=="beforeTSS"){
        
        ###THIS WILL IMPLY THAT THE GENE IS NOT DUPLICATED
        ##SO THE ONLY POSSIBLE THING, in this context where two tads involved
        ## IS THAT THIS IS PATHOGENIC IS BECAUSE THE GENE IS TRUNCATED BY THE DUPLICATION
        ##Because if not everything duplicated, enhancer, go to a neoTAD so has no effect on the non duplicated gene
        #Unless as said is truncated
        
        ##I'm not going too think too much about this, makes not much sense
        ##And in addition enh will not be painted (because gene direct effect)
        ##So just leave it the way it is 
        
        ##Just going to say, paint the breakpoint in the middle of the TAD
        ##I will paint the breakp then in the middle of the TAD, where the cluster of enh is for now located
        breakPos<-tad_X_cord[3]
        
      }
    }
  }
  
  #########################################################################
  ##Defining Breakpoint position
  breakp_x_coord<-c(breakPos,breakPos)##Tad central point for gene broken
  breakp_y_coord<-c(tad_Y_cord[1]-7, tad_Y_cord[3]+3)
  
  #######################
  ##Painting Breakpoint
  breakpointColor<-"brown3"
  lines(x=breakp_x_coord,
        y=breakp_y_coord,
        col=breakpointColor,
        lwd=2,
        lty=3)
  
  #Breakpoint label
  boxed.labels(x=breakp_x_coord[1],  y=breakp_y_coord[1] - 1,labels = "Breakpoint", cex=0.8, border = NA, bg ="white", 
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col=breakpointColor
  )
  
  #TAD label
  #text(x=tad_X_cord[3], y=10, label="TAD", cex = 1.5)
  boxed.labels(x=tad_X_cord[3],  y=tad_Y_cord[1]+10,labels = "TAD", cex=1.5,
               border = NA, bg = NA, xpad=1, ypad = 2 )
  
  ###Adding black line below TAD
  lines(x = c(tad_X_cord[1]-3,tad_X_cord[2]+3),
        y= rep.int(x=tad_Y_cord[1], times = 2), 
        col="black", lwd=2, lty=1)
  
  #####################################
  ## On this TAD only color enhancers
  ## place them on the center of the TAD
  
  ###
  ##Adding enhancers
  ###
  
  if(nEnh_other_domain>0){

    draw.ellipse(x=enh_x_positions,
                 y=enh_y_positions,a=0.1,col="#b2d235")
    
    #enhancers Label
    #text(x=enhXpos + 0.4, y=-3, label="enhancers", cex = 0.8)
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                 -distance_Yaxis_EnhLabel,
                 labels = "enhancers", cex=0.8, 
                 border = NA, bg ="white", 
                 xpad=1,
                 ypad=1.1 #To allow the text to breath
    )
    #Nr of enhancers
    #text(x=enhXpos + 0.4, y=-4, label=paste0("n=",nEnh_other_domain), cex=enhNumbersSize)
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                   -distance_Yaxis_EnhLabel
                   -distance_Yaxis_EnhLabel_EnhNumber,
                 labels = paste0("n=",nEnh_other_domain), cex=enhNumbersSize, 
                 border = NA, bg ="white", 
                 xpad=1.2,
                 ypad=1.2 #To allow the text to breath
    )
    
    # ##enhancers horizontal line over them
    # lines(x = c(enh_x_positions[1]-0.1,
    #             enh_x_positions[length(enh_x_positions)]+0.1),
    #       y = c(tad_Y_cord[1]+2, tad_Y_cord[1]+2),
    #       lwd=1.5)
    
    
    # ##So, the secondary TAD is on the right side
    # curvedarrow(from = c(enhXpos+0.5,2), to = c(geneCenter+0.25,2),
    #             arr.pos = 0.2, arr.type="T", curve = 0.25, segment = c(0,0.2))
    
    
    # ##FOR NOW, NOT COLORING CURVED ARROW IN THIS CONTEXT
    # if(enhXpos > geneCenter){
    #   ##So, the secondary TAD is on the right side
    #   curvedarrow(from = c(enhXpos+0.5,tad_Y_cord[1]+2), to = c(enhXpos-2,tad_Y_cord[1]+6),
    #               arr.pos = 0, arr.type="T", curve = 0.1, segment = c(0,1))
    #   
    #   ##PAINTING VERTICAL LINE ON OWR OWN FUCKIN CURVEDARROW rising problems//bugs
    #   ##So arr.pos=0 to be painted in the enh horizontal line and hence disappear
    #   lines(x=c(enhXpos-2,enhXpos-2),
    #         y=c(tad_Y_cord[1]+5, tad_Y_cord[1]+7),#c(5,7),
    #         lwd=2)
    #   
    # }else{
    #   ##So the secondary TAD is on the left side, adjust row curvature, to mantain it over the X axis
    #   curvedarrow(from = c(enhXpos+0.5,tad_Y_cord[1]+2), to = c(enhXpos+2,tad_Y_cord[1]+6),
    #               arr.pos = 0, arr.type="T", curve = -0.1, segment = c(0,1))
    #   
    #   ##PAINTING VERTICAL LINE ON OWR OWN FUCKIN CURVEDARROW rising problems//bugs
    #   ##So arr.pos=0 to be painted in the enh horizontal line and hence disappear
    #   lines(x=c(enhXpos+2,enhXpos+2),
    #         y=c(tad_Y_cord[1]+5, tad_Y_cord[1]+7),
    #         lwd=2)
    # }
    

    
    
  }
  
  ##################################################################
  ## Assembling relevant info from paint_Enhancer_WT_Secondary_TAD 
  ##################################################################
  info_drawing<-list("gene_relocated"=gene_relocated, ##TRUE,FAlSE regarding whether gene relocated or not (if not it stays destral or sinistral, attending to var situation)
                     "secondaryTAD_breakP"=breakPos)
  
  return(info_drawing)
}


# paintInterTAD_WT_lines<-function(tad_X_cord){
#   ###############################################
#   ## Adding dotted line connecting WT upper TADs
#   ###############################################
#   lines(x = c(tad_X_cord[2]+3, tad_X_cord[2]+8),
#         y=c(0,0), 
#         col="black", lwd=2, lty=2)
#   
#   ##Adding short diagonal lines over the dotted lines
#   lines(x = c(tad_X_cord[2]+3+2, tad_X_cord[2]+3+3),
#         y=c(-2,2), 
#         col="black", lwd=2, lty=1)
#   
#   lines(x = c(tad_X_cord[2]+3+2+0.4, tad_X_cord[2]+3+3+0.4),
#         y=c(-2,2), 
#         col="black", lwd=2, lty=1)
# }


paintInterTAD_arrows<-function(){
  #Painting arrows to indicate rearrangement
  
  ##Arrow to the right
  curvedarrow(from = c(23,3), to = c(18,3),
              arr.pos = 1, arr.type="simple", angle=40)
  
  ##Arrow to the left
  curvedarrow(from = c(18,-3), to = c(23,-3),
              arr.pos = 1, arr.type="simple", angle=40)
}

###############
##
paint_SV_labels<-function(sv_type, xAxisLim, yAxisPos){
  
  xCenter <- xAxisLim[1] + (xAxisLim[2]-xAxisLim[1])/2
  
  
  ###Adding Label
  boxed.labels(x=xCenter,  y=yAxisPos,labels = sv_type, cex=2, border = NA, bg ="white", 
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col="#d64045"
  )
  
}






