###############################################################
## Script to hold the functions to plot the SV rearrangements
###############################################################


########################################################
## Function to paint GENE SV re-arranged TAD 
## For Inversions or Translocations
########################################################

paintGene_SV_TAD<-function(nEnh_initial_left, nEnh_initial_right, gene, gene_breakp_line_type, situation,
                           nEnh_kept_left, nEnh_kept_right,nEnh_gained,
                           nEnh_other_domain,
                           info_drawingGENE_TAD, info_drawingSecondaryTAD,
                           tad_XCoord_OnLeftSide, tad_XCoord_OnRightSide,
                           tad_YCoord_Rearrangements,
                           geneBreakP_Position_respectToTSS){
  
  ########################################################
  ## Function to paint GENE SV re-arranged TAD 
  ########################################################

  if((situation=="primaryTAD_Sinistral")&&(info_drawingSecondaryTAD$gene_relocated==FALSE)){
    ##Gene was on the left part, and stays onf the left part
    ##So gene TAD to create painted on the left side 
    gene_sv_TAD_x_coord<-tad_XCoord_OnLeftSide
    
  }else if((situation=="primaryTAD_Dextral")&&(info_drawingSecondaryTAD$gene_relocated==TRUE)){
    
    ##So gene SV TAD to create painted on the left side 
    gene_sv_TAD_x_coord<-tad_XCoord_OnLeftSide
  }else{
    
    ##So gene TAD to create painted on the left side 
    gene_sv_TAD_x_coord<-tad_XCoord_OnRightSide
  }
  
  ###########################################
  # Painting main scaffold TAD 
  
  ### Y axis position of the Inversion WT row
  tad_Y_cord<-tad_YCoord_Rearrangements

  
  ## Let's mantain color codes by TAD position
  ## If TAD is on the left blue, if it is on the right orange
  ## Regardless whether it is the gene or the secondary domain TAD
  ## In order not to confound the user if > 1 gene is affected
  
  tad_X_cord<-gene_sv_TAD_x_coord
  
  if(tad_X_cord[1]<20){
    colorTAD<-"#acdfeb"
  }else if(tad_X_cord[1]>20){
    colorTAD<-"#e5edb2"
  }
  
  tadColors<-c("#acdfeb","#e5edb2")
    
  
  #########################################
  ## Defining additional relevant variables
  #########################################
  enhNumbersSize<-0.8
  distance_Yaxis_geneLabel<-2
  distance_Yaxis_EnhGene<-2
  distance_Yaxis_EnhLabel_EnhNumber<-1.3
  
  ###################
  ##Painting TAD
  polygon(tad_X_cord, tad_Y_cord, col = colorTAD, border = "#ffffff")
  
  ###############################
  ##Obtaining breakpoint position
  ##If gene is relocated we paint the breakpos at the contrary tad breakpos
  ##If not, we paint its breakpos
  
  if(info_drawingSecondaryTAD$gene_relocated==TRUE){
    breakPos<-info_drawingSecondaryTAD$secondaryTAD_breakP
  }else{
    ##So the gene stays on its TAD
    breakPos<-info_drawingGENE_TAD$geneTAD_breakP
  }
  
  ##Defining Breakpoint position
  breakp_x_coord<-c(breakPos,breakPos)##Tad central point for gene broken
  breakp_y_coord<-c(tad_Y_cord[1]-9,##-7 before
                    tad_Y_cord[3]+3)
  
  
  ###########################
  ## Painting the other color Figure
  ##   "otherColorFigure"
  ## The other color covering the TAD coming from the rearranged domain
  ############################
  
  ##coordinates in bottom-left bottom-right top-right top-left order
  
  ##########################################################
  ###First Getting The Coordinates of the "Other Figure"
  ##########################################################
  
  if(tad_X_cord[1]<20){
    ##It means we are on the LEFT SIDE so the other color goes from the TAD breakpos to the end of the TAD
    
    ##if breakpos before TAD center, we are going to color an irregular QUADRILATERAl
    ##So we will need 4 coordinates for the figure
    ##If its in the TAD center or after we will color a Triangle, so we will need only 3 coordinates
    if(breakPos>tad_X_cord[3]){##tad_XCoord_OnLeftSide[3],##TAD center position
      
      ##Right triangle to paint
      otherColorFigure_X_pos<-c(breakPos,
                                tad_X_cord[2],
                                breakPos)
      
      ###############################################
      ##First getting the triangle corner angle
      
      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs
      
      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm
      
      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians
      
      ##Now let's compute the height of interest. With respect to the right triangle, so right side
      catetoPlaying<-(tad_X_cord[2]-breakPos)
      heightOfInterest<-tan(cornerAngle)*catetoPlaying
      
      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest
      
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                mostWanted_Y_pos)
      
      
    }else if(breakPos==tad_X_cord[3]){

      ##Right triangle to paint
      otherColorFigure_X_pos<-c(breakPos,
                                tad_X_cord[2],
                                breakPos)
      
      ##In this case very easy the most wanted Y coordinate, is just the TAD peak
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                tad_Y_cord[3])
      
    }else if(breakPos<tad_X_cord[3]){
      ##So QUADRILATERAL to paint
      otherColorFigure_X_pos<-c(breakPos,
                                tad_X_cord[2],
                                tad_X_cord[3],
                                breakPos)
      
      ###############################################
      ##First getting the triangle corner angle
      
      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs
      
      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm
      
      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians
      
      ##Now let's compute the height of interest
      ##We compute it from the perspective of the right triangle!!! so Left side of the TAD
      catetoPlaying<-(breakPos-tad_X_cord[1])
      heightOfInterest<-tan(cornerAngle)*catetoPlaying
      
      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest
      
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                tad_Y_cord[3],
                                mostWanted_Y_pos
                                )
      
    }
    
    
  }else if(tad_X_cord[1]>20){
    ##It means we are on the RIGHT SIDE so the other color goes from the TAD start to the breakpos
    
    ##if breakpos after TAD center, we are going to color an irregular QUADRILATERAl
    ##So we will need 4 coordinates for the figure
    ##If its in the TAD center or before we will color a Triangle, so we will need only 3 coordinates
    if(breakPos>tad_X_cord[3]){##tad_XCoord_OnLeftSide[3],##TAD center position
      
      ##So QUADRILATERAL to paint
      otherColorFigure_X_pos<-c(tad_X_cord[1],
                                breakPos,
                                breakPos,
                                tad_X_cord[3])
      
      ###############################################
      ##First getting the triangle corner angle

      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs
      
      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm
      
      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians
      
      ##Now let's compute the height of interest
      catetoPlaying<-(tad_X_cord[2]-breakPos)
      heightOfInterest<-tan(cornerAngle)*catetoPlaying
      
      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest
      
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                mostWanted_Y_pos,
                                tad_Y_cord[3])
      
    }else if(breakPos==tad_X_cord[3]){
      ##Right triangle to paint
      otherColorFigure_X_pos<-c(tad_X_cord[1],
                                breakPos,
                                breakPos)
      
      ##In this case very easy the most wanted Y coordinate, is just the TAD peak
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                tad_Y_cord[3])
      
      
    }else if(breakPos<tad_X_cord[3]){
      ##Right triangle to paint
      otherColorFigure_X_pos<-c(tad_X_cord[1],
                                breakPos,
                                breakPos)
      
      ###############################################
      ##First getting the triangle corner angle
      
      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs
      
      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm
      
      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians
      
      ##Now let's compute the height of interest. With respect to the right triangle
      catetoPlaying<-(breakPos-tad_X_cord[1])
      heightOfInterest<-tan(cornerAngle)*catetoPlaying
      
      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest
      
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                mostWanted_Y_pos)
    }
    
  }

  ##########################################################
  ### Painting the "Other Figure"
  
  color_otherFigure<-tadColors[!(tadColors==colorTAD)]##the opposite of the TAD color
  
  polygon(otherColorFigure_X_pos, otherColorFigure_Y_pos, col = color_otherFigure, border = "#ffffff")
  
  #######################
  ##Painting Breakpoint
  breakpointColor<-"brown3"
  lines(x=breakp_x_coord,
        y=breakp_y_coord,
        col=breakpointColor,
        lwd=2,
        lty=3)
  
  #Breakpoint label
  ##Boxed labels: xpad,ypad arguments, used to expand the box beyond the text
  boxed.labels(x=breakp_x_coord[1],  y=breakp_y_coord[1] - 1,labels = "Breakpoint", cex=0.8, border = NA, bg ="white", 
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col=breakpointColor
  )
  
  #TAD label
  # text(x=tad_X_cord[3], y=tad_Y_cord[1]+10, label="TAD", cex = 1.5)
  boxed.labels(x=tad_X_cord[3],  y=tad_Y_cord[1]+10,labels = "TAD", cex=1.5,
               border = NA, bg= NA, xpad=1, ypad = 2 )
  
  ###Adding black line below TAD
  lines(x = c(tad_X_cord[1]-3,tad_X_cord[2]+3),
        y=c(tad_Y_cord[1],tad_Y_cord[1]), 
        col="black", lwd=2, lty=1)
  
  #########################################
  ##Defining some profitable variables
  #########################################
  
  distanceEnhFromGene<-1.5 ##Distance between the enhancers and the gene
  
  sizeEnhRegion<-1##Touching this will not really affect enh size, because created regardless of this
  
  ###############################
  ## Painting gene body 
  ###############################
  ##if gene relocated it will be place in the middle of the "otherFigure" X space
  ##if it is not relocated it will be painted exactly at the same position where it was in the WT condition
  
  if(info_drawingSecondaryTAD$gene_relocated==TRUE){
    ##so the gene will be painted in the middle of the "otherFigure" X space
    #other Figure X space
    #recall, coordinates in order, bottom left, bottom right, top right...
    #so the X space always provided by the two first coord in the "other Figure", regardless whether it is a triangle or a quadrilater
    x_space<-c(otherColorFigure_X_pos[1],otherColorFigure_X_pos[2])
    
    x_space_start<-x_space[1]
    x_space_end<-x_space[2]
    
    #gene position at the center -0.5 (since we provide the left coord over which the gene is painted)
    genePos<-(x_space_start + (x_space_end-x_space_start)/2)-0.5

  }else{
    ##Hence, the gene is not relocated
    ##So the gene is painted again at the same place where it was in the WT situation
    
    genePos<-info_drawingGENE_TAD$genePos
    
  }
  
  
  geneCenter<-genePos + 0.5 ##Used to direct the arrows to here
  gene_X_cord<-c(genePos,genePos+1,genePos+1,genePos)
  gene_Y_cord<-c(tad_Y_cord[1]-1,tad_Y_cord[1]-1,tad_Y_cord[1]+1,tad_Y_cord[1]+1)
  
  ###################
  #Adding gene body
  polygon(gene_X_cord, gene_Y_cord, col = "#0e3d61")
  
  #gene Label
  ##I want it to overlay the breakpoint line
  # https://stackoverflow.com/questions/45366243/text-labels-with-background-colour-in-r
  boxed.labels(x=geneCenter,  y=tad_Y_cord[1]-distance_Yaxis_geneLabel,
               labels = gene, cex=0.8, border = NA, bg ="white", 
               xpad=1,
               ypad=2##To allow the text to breathe
  )
  
  
  #####################
  ##Adding enhancers
  #####################
  enh_y_positions<-rep.int(x=tad_Y_cord[1], times = 4)##c(0,0,0,0)
  
  ##################################################
  ## First, PAINTING COGNATE ENHANCERS
  ## changes on enh originally surrounding the gene
  ##################################################
  
  ##If gene relocated what it remains on the left now it is on the right and same logic for the ones remaining on the right
  
  #If gene not relocated what was on the left, stays on the left, and what was on the right stays on the right
  #And if nkept, different than nInitial, as the breakpoint in the middle, we will paint only two enh, to represent
  #a reduction on its number
  
  if(info_drawingSecondaryTAD$gene_relocated=="FALSE"){
    
    #If gene not relocated what was on the left, stays on the left, and what was on the right stays on the right
    #And if nkept, different than nInitial, as the breakpoint in the middle, we will paint only two enh, to represent
    #a reduction on its number. But the two left or the two right enh, depending on where the breakpoint is
    
    ####################################
    ####To the Right side gene Position
    ####################################
    
    if((geneBreakP_Position_respectToTSS=="beforeTSS")||(gene_breakp_line_type=="afterTSS_removing_none")){
      ##Enhancers to  the right left intact, same as they were

      if(nEnh_initial_right>0){
        ##So there is at least 1 enh
        enhXpos<-genePos+distanceEnhFromGene##left margin of the enh cluster
        enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
        
        draw.ellipse(x=enh_x_positions,
                     y=rep.int(x = enh_y_positions[1],times = length(enh_x_positions)),a=0.1,col="#b2d235")
        
        #enhancers Label
        #text(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel              -distance_Yaxis_EnhGene, label="enhancers", cex = 0.8)
        boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene,
                     labels = "enhancers", cex=0.8, 
                     border = NA, bg ="white", 
                     xpad=1,
                     ypad=1 #To allow the text to breath
        )
        #Nr of enhancers
        #text(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel-distance_Yaxis_EnhGene-distance_Yaxis_EnhLabel_EnhNumber, label=paste0("n=",nEnh_initial_right), cex=enhNumbersSize)
        boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene
                     -distance_Yaxis_EnhLabel_EnhNumber,
                     labels = paste0("n=",nEnh_initial_right), cex=enhNumbersSize, 
                     border = NA, bg ="white", 
                     xpad=1.2,
                     ypad=1 #To allow the text to breath
        )
        
        ##enhancers horizontal line over them
        lines(x = c(enh_x_positions[1]-0.1,
                    enh_x_positions[length(enh_x_positions)]+0.1),
              y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
              lwd=1.5)
        
        curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                       ,enh_y_positions[1]+2),
                    arr.pos = 0, arr.type="T")
      }
      
    }else if((geneBreakP_Position_respectToTSS=="afterTSS") && (nEnh_initial_right>0)){
      
      
      if(gene_breakp_line_type=="afterTSS_removing_all"){
        ##DOING NOTHING, we don't want to paint anything as all the enhancers are going to disappear
      }else if(gene_breakp_line_type=="afterTSS_removing_some"){
        
        ###########################
        ## PAINTING 2 ENHANCERS
        
        ##To quickly give the impression that some enhancers are lost, I'm going to paint 2enhancers
        ##That will be the two right enhn if the gene is on the right part of the plot, or the two left
        ##So there is at least 1 enh
        enhXpos<-genePos+distanceEnhFromGene##left margin of the enh cluster
        
        ##Painting Only 2 enhancers
        enh_x_positions<-c(enhXpos,enhXpos+0.3)
        
        draw.ellipse(x=enh_x_positions,
                     y=enh_y_positions[1:2],a=0.1,col="#b2d235")
        
        #enhancers Label
        boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene,
                     labels = "enhancers", cex=0.8, 
                     border = NA, bg ="white", 
                     xpad=1,
                     ypad=1 #To allow the text to breath
        )
        #Nr of enhancers
        boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene
                     -distance_Yaxis_EnhLabel_EnhNumber,
                     labels = paste0("n=",nEnh_kept_right), cex=enhNumbersSize, 
                     border = NA, bg ="white", 
                     xpad=1.2,
                     ypad=1 #To allow the text to breath
        )
        
        ##enhancers horizontal line over them
        lines(x = c(enh_x_positions[1]-0.1,
                    enh_x_positions[length(enh_x_positions)]+0.1),
              y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
              lwd=1.5)
        
        curvedarrow(from = c(enhXpos+0.15,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                       ,enh_y_positions[1]+2),
                    arr.pos = 0, arr.type="T" )
        
        
        ###END OF PAINTING 2 ENHANCERS
        ###################################
        
      }
      
    }
    
    
    
    ###################################
    ####To the Left side gene Position
    ###################################

    #the gene is not relocated, so the enh to the left will be kept intact if:
    ##the breakp is after the tss (gene not relocated because it is a Sinistral TAD, located at the left part of the screen)
    ##the breakp is before the tss removing none enh (gene not relocated because it is a Destral TAD, located at the right part of the screen)
    if((geneBreakP_Position_respectToTSS=="afterTSS")||(gene_breakp_line_type=="beforeTSS_removing_none")){
      ##Enhancers to  the LEFT kept intact, same as they were
      
      if(nEnh_initial_left>0){
        enhXpos<-genePos-distanceEnhFromGene##left margin of the enh cluster
        enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
        
        draw.ellipse(x=enh_x_positions,
                     y=enh_y_positions,a=0.1,col="#b2d235")
        
        #enhancers Label
        boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene,
                     labels = "enhancers", cex=0.8, 
                     border = NA, bg ="white", 
                     xpad=1,
                     ypad=1 #To allow the text to breath
        )
        
        #Nr of enhancers
        boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene
                     -distance_Yaxis_EnhLabel_EnhNumber,
                     labels = paste0("n=",nEnh_kept_left), cex=enhNumbersSize, 
                     border = NA, bg ="white", 
                     xpad=1.2,
                     ypad=1 #To allow the text to breath
        )
        
        ##enhancers horizontal line over them
        lines(x = c(enh_x_positions[1]-0.1,
                    enh_x_positions[length(enh_x_positions)]+0.1),
              y = c(enh_y_positions[1]+2,enh_y_positions[1]+2),
              lwd=1.5)
        
        curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter#-0.25
                                                                       ,enh_y_positions[1]+2),
                    arr.pos = 0, arr.type="T", curve = -1 )
      }
      
    }else if((geneBreakP_Position_respectToTSS=="beforeTSS") && (nEnh_initial_left>0)){
      
      
      if(gene_breakp_line_type=="beforeTSS_removing_all"){
        ##DOING NOTHING, we don't want to paint anything as all the enhancers are going to disappear
      }else if(gene_breakp_line_type=="beforeTSS_removing_some"){
        
        ###########################
        ## PAINTING 2 ENHANCERS
        
        ##To quickly give the impression that some enhancers are lost, I'm going to paint just 2 enhancers
        ##So there is at least 1 enh
        enhXpos<-genePos-distanceEnhFromGene##left margin of the enh cluster
        
        ##Painting Only 2 enhancers
        #enh_x_positions<-c(enhXpos,enhXpos+0.3)
        enh_x_positions<-c(enhXpos+0.6, enhXpos+0.9)##The two second enhancers
        enhXpos<-enh_x_positions[1]
        
        draw.ellipse(x=enh_x_positions,
                     y=enh_y_positions[1:2],a=0.1,col="#b2d235")
        
        #enhancers Label
        boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene,
                     labels = "enhancers", cex=0.8, 
                     border = NA, bg ="white", 
                     xpad=1,
                     ypad=1 #To allow the text to breath
        )
        #Nr of enhancers
        boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene
                     -distance_Yaxis_EnhLabel_EnhNumber,
                     labels = paste0("n=",nEnh_kept_left), cex=enhNumbersSize, 
                     border = NA, bg ="white", 
                     xpad=1.2,
                     ypad=1 #To allow the text to breath
        )
        
        ##enhancers horizontal line over them
        lines(x = c(enh_x_positions[1]-0.1,
                    enh_x_positions[length(enh_x_positions)]+0.1),
              y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
              lwd=1.5)
        
        enhHorizontalLineCoord<-c(enh_x_positions[1]-0.1,
                                  enh_x_positions[length(enh_x_positions)]+0.1)
        
        enhHorizLineCenter<-enhHorizontalLineCoord[1]+(enhHorizontalLineCoord[2]-enhHorizontalLineCoord[1])/2
        
        
        curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                        ,enh_y_positions[1]+2),
                    arr.pos = 0, arr.type="T", curve=-1 )
        
        
        ###END OF PAINTING 2 ENHANCERS
        ###################################
        
      }
      
    }
    
  }else{
    ##So
    ##############################
    # THE GENE IS RELOCATED
    ##############################

    ##Hence inverting location of the enhancers
    
    if(geneBreakP_Position_respectToTSS=="beforeTSS"){
      ##so situation is "primaryTAD_Sinistral" (because the gene ir relocated)
      
      #the enh initially to the right of the gene now will be placed at the left side of the gene
      
      ###########################
      ## To the LEFT-SV-scenario 
      ###########################
      #Those that initially were on the right side
      if(nEnh_initial_right>0){
        enhXpos<-genePos-distanceEnhFromGene##left margin of the enh cluster
        enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
        
        draw.ellipse(x=enh_x_positions,
                     y=enh_y_positions,a=0.1,col="#b2d235")
        
        #enhancers Label
        boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene,
                     labels = "enhancers", cex=0.8, 
                     border = NA, bg ="white", 
                     xpad=1,
                     ypad=1 #To allow the text to breath
        )
        
        #Nr of enhancers
        boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene
                     -distance_Yaxis_EnhLabel_EnhNumber,
                     labels = paste0("n=",nEnh_initial_right), cex=enhNumbersSize, 
                     border = NA, bg ="white", 
                     xpad=1.2,
                     ypad=1 #To allow the text to breath
        )
        
        ##enhancers horizontal line over them
        lines(x = c(enh_x_positions[1]-0.1,
                    enh_x_positions[length(enh_x_positions)]+0.1),
              y = c(enh_y_positions[1]+2,enh_y_positions[1]+2),
              lwd=1.5)
        
        curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter#-0.25
                                                                       ,enh_y_positions[1]+2),
                    arr.pos = 0, arr.type="T", curve = -1 )
      }
      
      
      ###########################
      ## To the RIGHT-SV-scenario 
      ###########################
      #Lo que estaba a la izquierda initially pasa a estar a la derecha
      
      #Si el breakp es del tipo beforeTSS_remove_none pinto todos los de la izquierda iniciales a la derecha
      if(gene_breakp_line_type=="beforeTSS_removing_none"){
        ##To the right, those that initially were on the left
        
        if(nEnh_initial_left>0){
          ##So there is at least 1 enh
          enhXpos<-genePos+distanceEnhFromGene##left margin of the enh cluster
          enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
          
          draw.ellipse(x=enh_x_positions,
                       y=enh_y_positions,a=0.1,col="#b2d235")
          
          #enhancers Label
          boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                       -distance_Yaxis_geneLabel
                       -distance_Yaxis_EnhGene,
                       labels = "enhancers", cex=0.8, 
                       border = NA, bg ="white", 
                       xpad=1,
                       ypad=1 #To allow the text to breath
          )
          #Nr of enhancers
          boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                       -distance_Yaxis_geneLabel
                       -distance_Yaxis_EnhGene
                       -distance_Yaxis_EnhLabel_EnhNumber,
                       labels = paste0("n=",nEnh_initial_left), cex=enhNumbersSize, 
                       border = NA, bg ="white", 
                       xpad=1.2,
                       ypad=1 #To allow the text to breath
          )
          
          ##enhancers horizontal line over them
          lines(x = c(enh_x_positions[1]-0.1,
                      enh_x_positions[length(enh_x_positions)]+0.1),
                y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
                lwd=1.5)
          
          curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                         ,enh_y_positions[1]+2),
                      arr.pos = 0, arr.type="T" )
        }
        
      }else if(gene_breakp_line_type=="beforeTSS_removing_all"){
        ##DO NOTHING, NOTHING PAINTED
      }else if(gene_breakp_line_type=="beforeTSS_removing_some"){
        ##TWO ENHANCERS ARE PAINTED
        ###########################
        ## PAINTING 2 ENHANCERS
        
        ##To quickly give the impression that some enhancers are lost, I'm going to paint just the 2 first enhancers
        ##So there is at least 1 enh
        enhXpos<-genePos+distanceEnhFromGene##left margin of the enh cluster
        
        ##Painting Only 2 enhancers
        enh_x_positions<-c(enhXpos,enhXpos+0.3)
        
        draw.ellipse(x=enh_x_positions,
                     y=enh_y_positions[1:2],a=0.1,col="#b2d235")
        
        #enhancers Label
        boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene,
                     labels = "enhancers", cex=0.8, 
                     border = NA, bg ="white", 
                     xpad=1,
                     ypad=1 #To allow the text to breath
        )
        #Nr of enhancers
        boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene
                     -distance_Yaxis_EnhLabel_EnhNumber,
                     labels = paste0("n=",nEnh_kept_left), cex=enhNumbersSize, 
                     border = NA, bg ="white", 
                     xpad=1.2,
                     ypad=1 #To allow the text to breath
        )
        
        ##enhancers horizontal line over them
        lines(x = c(enh_x_positions[1]-0.1,
                    enh_x_positions[length(enh_x_positions)]+0.1),
              y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
              lwd=1.5)
        
        curvedarrow(from = c(enhXpos+0.15,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                        ,enh_y_positions[1]+2),
                    arr.pos = 0, arr.type="T" )
        
        
        ###END OF PAINTING 2 ENHANCERS
        ###################################
      }
      
      
    }else if(geneBreakP_Position_respectToTSS=="afterTSS"){
      #the enh initially to the left of the gene now will be placed at the right side of the gene
      
      ###########################
      ## To the RIGHT-SV-scenario
      ###########################
      #Those that initially were on the left side
      if(nEnh_initial_left>0){
        ##So there is at least 1 enh
        enhXpos<-genePos+distanceEnhFromGene##left margin of the enh cluster
        enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
        
        draw.ellipse(x=enh_x_positions,
                     y=enh_y_positions,a=0.1,col="#b2d235")
        
        #enhancers Label
        #text(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel              -distance_Yaxis_EnhGene, label="enhancers", cex = 0.8)
        boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene,
                     labels = "enhancers", cex=0.8, 
                     border = NA, bg ="white", 
                     xpad=1,
                     ypad=1 #To allow the text to breath
        )
        #Nr of enhancers
        #text(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel-distance_Yaxis_EnhGene-distance_Yaxis_EnhLabel_EnhNumber, label=paste0("n=",nEnh_initial_right), cex=enhNumbersSize)
        boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene
                     -distance_Yaxis_EnhLabel_EnhNumber,
                     labels = paste0("n=",nEnh_initial_left), cex=enhNumbersSize, 
                     border = NA, bg ="white", 
                     xpad=1.2,
                     ypad=1 #To allow the text to breath
        )
        
        ##enhancers horizontal line over them
        lines(x = c(enh_x_positions[1]-0.1,
                    enh_x_positions[length(enh_x_positions)]+0.1),
              y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
              lwd=1.5)
        
        curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                       ,enh_y_positions[1]+2),
                    arr.pos = 0, arr.type="T" )
      }
      
      ###########################
      ## To the LEFT-SV-scenario 
      ###########################
      #Lo que estaba a la derecha pasa a estar a la iquierda
      
      #Si el breakp es del tipo afterTSS_remove_none pinto todos los de la derecha iniciales a la izquierda
      if(gene_breakp_line_type=="afterTSS_removing_none"){
        ##To the left, ALL those that initially were on the right
        if(nEnh_initial_right>0){
          enhXpos<-genePos-distanceEnhFromGene##left margin of the enh cluster
          enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
          
          draw.ellipse(x=enh_x_positions,
                       y=enh_y_positions,a=0.1,col="#b2d235")
          
          #enhancers Label
          boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                       -distance_Yaxis_geneLabel
                       -distance_Yaxis_EnhGene,
                       labels = "enhancers", cex=0.8, 
                       border = NA, bg ="white", 
                       xpad=1,
                       ypad=1 #To allow the text to breath
          )
          
          #Nr of enhancers
          boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                       -distance_Yaxis_geneLabel
                       -distance_Yaxis_EnhGene
                       -distance_Yaxis_EnhLabel_EnhNumber,
                       labels = paste0("n=",nEnh_initial_right), cex=enhNumbersSize, 
                       border = NA, bg ="white", 
                       xpad=1.2,
                       ypad=1 #To allow the text to breath
          )
          
          ##enhancers horizontal line over them
          lines(x = c(enh_x_positions[1]-0.1,
                      enh_x_positions[length(enh_x_positions)]+0.1),
                y = c(enh_y_positions[1]+2,enh_y_positions[1]+2),
                lwd=1.5)
          
          curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter#-0.25
                                                                         ,enh_y_positions[1]+2),
                      arr.pos = 0, arr.type="T", curve = -1)
        }
        
      }else if(gene_breakp_line_type=="afterTSS_removing_all"){
        ##DO NOTHING, NOTHING PAINTED
      }else if(gene_breakp_line_type=="afterTSS_removing_some"){
        ##TWO ENHANCERS ARE PAINTED
        ###########################
        ## PAINTING 2 ENHANCERS
        
        ##To quickly give the impression that some enhancers are lost, I'm going to paint just the 2 first enhancers
        ##So there is at least 1 enh
        enhXpos<-genePos-distanceEnhFromGene##left margin of the enh cluster
        
        ##Painting Only 2 enhancers
        enh_x_positions<-c(enhXpos,enhXpos+0.3)
        
        draw.ellipse(x=enh_x_positions,
                     y=enh_y_positions[1:2],a=0.1,col="#b2d235")
        
        #enhancers Label
        boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel
                     -distance_Yaxis_EnhGene,
                     labels = "enhancers", cex=0.8, 
                     border = NA, bg ="white", 
                     xpad=1,
                     ypad=1 #To allow the text to breath
        )
        #Nr of enhancers
        boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                     -distance_Yaxis_geneLabel-distance_Yaxis_EnhGene
                     -distance_Yaxis_EnhLabel_EnhNumber,
                     labels = paste0("n=",nEnh_kept_right), cex=enhNumbersSize, 
                     border = NA, bg ="white", 
                     xpad=1.2,
                     ypad=1 #To allow the text to breath
        )
        
        ##enhancers horizontal line over them
        lines(x = c(enh_x_positions[1]-0.1,
                    enh_x_positions[length(enh_x_positions)]+0.1),
              y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
              lwd=1.5)
        
        curvedarrow(from = c(enhXpos+0.15,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                        ,enh_y_positions[1]+2),
                    arr.pos = 0, arr.type="T", curve=-1)
        
        
        ###END OF PAINTING 2 ENHANCERS
        ###################################
        
      }
    }
  }
  
  
  ##################################################
  ## TIME FOR PAINTING ECTOPIC  ENHANCERS
  ## changes on enh originally surrounding the gene
  ##################################################
  ##If gene relocated, the ectopic enhancers are displayed at the position where they were
  ##If there is no gain, no paint
  ##If nEnh gained == nEnhInitialOther Domain, the four enh are painted
  ##On the contrary only two enh are painted
  
  ##If the gene is not relocated we paint the enh on the middle of its segment
  

  if(nEnh_gained>0){
    ##So some enh need to be painted
    
    if(info_drawingSecondaryTAD$gene_relocated==TRUE){
      #The enhancers occupy the exact same position than in the WT tad so the TAD center +- 0.5 approx
      ##So the enhancers are going to be painted in the middle of their TAD
      
      ##As it is the gene the one that is relocated, here we work with tad_X_cord
      enhXpos<-tad_X_cord[3]-0.45
      
      if(nEnh_gained == nEnh_other_domain){
        #The enhancers occupy the exact same position than in the WT tad so the TAD center +- 0.5 approx
        ##All enh from the other Domain Gained so
        #we paint 4 enh
        
        enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
      }else if(nEnh_gained < nEnh_other_domain){
        ##So we are going to paint only 2 enhancers
        ##If the tadXcoord is on the right side, we paint the two right enh
        ##If the tadXcoord is on the left side, we paint the two left enh
        
        if(tad_X_cord[1]>20){
          ##So we paint the two enh on the right part
          enh_x_positions<-c(enhXpos+0.6, enhXpos+0.9)
        }else{
          ##so tad_X_cord<20
          ##So we paint the two enh on the left
          enh_x_positions<-c(enhXpos,enhXpos+0.3)
        }
      }
      
    }else if(info_drawingSecondaryTAD$gene_relocated==FALSE){
      ##When the gene not relocated, to locate the enh we work with the "other figure" coordinates
      ##Because the ectopic enhancers are the elements relocated now
      ##We paint the enh in the middle of the other Figure coord
      x_space<-c(otherColorFigure_X_pos[1],otherColorFigure_X_pos[2])
      
      x_space_start<-x_space[1]
      x_space_end<-x_space[2]
      
      #gene position at the center -0.5 (since we provide the left coord over which the gene is painted)
      enhXpos<-(x_space_start + (x_space_end-x_space_start)/2)-0.5
      enhXpos<-enhXpos-0.45
      
      if(nEnh_gained==nEnh_other_domain){
        ##All enh from the other Domain Gained so
        #we paint 4 enh
        
        ##As it is the gene the one that is relocated, here we work with tad_X_cords
        enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
        
      }else if(nEnh_gained < nEnh_other_domain){
        ##So we are going to paint only 2 enhancers

          ##so tad_X_cord<20
          ##So we paint the two enh on the left
          enh_x_positions<-c(enhXpos,enhXpos+0.3)
      }
    }
    
    ####################################################
    ## Pintamos los ECTOPIC enhancers
    ##Lo de arriba solo para calcular posicion de enh
    ##get enh positions center
    ##para pintar hacia derecha o izquierda mirar si el geneCenter queda a la izq o a la derecha
    
    ##
    draw.ellipse(x=enh_x_positions,
                 y=enh_y_positions,a=0.1,col="#b2d235")
    
    #enhancers Label
    #Below the level of the other enhancers, to avoid overlaps
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                 -distance_Yaxis_geneLabel
                 -distance_Yaxis_EnhGene
                 -distance_Yaxis_EnhLabel_EnhNumber
                 -distance_Yaxis_EnhLabel_EnhNumber,
                 labels = "enhancers", cex=0.8, 
                 border = NA, bg ="white", 
                 xpad=1,
                 ypad=1 #To allow the text to breath
    )
    #Nr of enhancers
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                 -distance_Yaxis_geneLabel
                 -distance_Yaxis_EnhGene
                 -distance_Yaxis_EnhLabel_EnhNumber
                 -distance_Yaxis_EnhLabel_EnhNumber
                 -distance_Yaxis_EnhLabel_EnhNumber, labels = paste0("n=",nEnh_gained), cex=enhNumbersSize, 
                 border = NA, bg ="white", 
                 xpad=1.2,
                 ypad=1 #To allow the text to breath
    )
    
    ##enhancers horizontal line over them
    lines(x = c(enh_x_positions[1]-0.1,
                enh_x_positions[length(enh_x_positions)]+0.1),
          y = c(enh_y_positions[1]+2,enh_y_positions[1]+2),
          lwd=1.5)
    
    
    enhHorizontalLineCoord<-c(enh_x_positions[1]-0.1,
                              enh_x_positions[length(enh_x_positions)]+0.1)
    
    enhHorizLineCenter<-enhHorizontalLineCoord[1]+(enhHorizontalLineCoord[2]-enhHorizontalLineCoord[1])/2
    
    ##Painting curve connecting enh with gene
    if(enhXpos > geneCenter){
      ##So, the secondary TAD is on the right side
      curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(geneCenter,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T", curve = 1)
    }else{
      ##So the secondary TAD is on the left side, adjust row curvature, to mantain it over the X axis
      curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(geneCenter,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T", curve = -1)
    }
    
  }
  
  ##################################################################################
  ##Adding Triangle Shape over gene, if it has any enh, to represent the arrow end
  ##################################################################################
  if((nEnh_kept_left>0)||(nEnh_kept_right>0)||(nEnh_gained>0)){
    triangle_X_Coord<-c(geneCenter-0.3, geneCenter, geneCenter+0.3)
    heightOrizontalLineEnh<-enh_y_positions[1]+2
    triangle_Y_Coord<-c(heightOrizontalLineEnh+0.3, heightOrizontalLineEnh-0.2, heightOrizontalLineEnh+0.3)
    polygon(triangle_X_Coord, triangle_Y_Coord, col = "black", border = "black")
  }
  
  
  
  
  #############################
  ## return profitable info ###
  #############################
  ##break pos, to easily select the other one
  info_drawing<-list("breakPos"=breakPos,
                     "geneCenter"=geneCenter)
  
  return(info_drawing)
  
}

#############################################################
## Function to paint The Secondary Domain SV re-arranged TAD 
## For Inversions or Translocations
#############################################################

paintSecondaryDomain_SV_TAD<-function(nEnh_initial_left, nEnh_initial_right, gene, gene_breakp_line_type, situation,
                           nEnh_other_domain,
                           nEnh_kept_left, nEnh_kept_right,nEnh_gained,
                           info_drawingGENE_TAD, info_drawingSecondaryTAD,
                           tad_XCoord_OnLeftSide, tad_XCoord_OnRightSide,
                           tad_YCoord_Rearrangements,
                           infoDrawing_gene_sv_tad){
  
  ##########################################################################################
  #We are coloring in this function the re-arranged domain where the gen of interest is not
  ##########################################################################################
  
  ##We are gonna paint, maximum, do clusters of enhancers at both sides of the breakpoint
  ##But it is important to know whether the ones from the left are from the gene or enh or viceversa
  ##Let's hence just track whether 
  
  ##Opposite logic than in paintGene_SV for this part
  if((situation=="primaryTAD_Sinistral")&&(info_drawingSecondaryTAD$gene_relocated==FALSE)){
    ##Gene was on the left part, and stays onf the left part
    ##So gene TAD to create painted on the left side 
    sv_TAD_x_coord<-tad_XCoord_OnRightSide
    
    locationGeneEnhVsOtherDomain<-"GeneEnhLost_Breakp_OtherEnh"
    
  }else if((situation=="primaryTAD_Dextral")&&(info_drawingSecondaryTAD$gene_relocated==TRUE)){
    
    ##So gene SV TAD to create painted on the left side 
    sv_TAD_x_coord<-tad_XCoord_OnRightSide
    
    locationGeneEnhVsOtherDomain<-"GeneEnhLost_Breakp_OtherEnh"
  }else{
    
    sv_TAD_x_coord<-tad_XCoord_OnLeftSide
  }
  

  ###########################################
  # Painting main scaffold TAD 
  
  ### Y axis position of the Inversion WT row
  tad_Y_cord<-tad_YCoord_Rearrangements
  
  
  ## Let's mantain color codes by TAD position
  ## If TAD is on the left blue, if it is on the right orange
  ## Regardless whether it is the gene or the secondary domain TAD
  ## In order not to confound the user if > 1 gene is affected
  
  tad_X_cord<-sv_TAD_x_coord
  
  if(tad_X_cord[1]<20){
    colorTAD<-"#acdfeb"
  }else if(tad_X_cord[1]>20){
    colorTAD<-"#e5edb2"
  }
  
  tadColors<-c("#acdfeb","#e5edb2")
  
  ###################
  ##Painting TAD
  polygon(tad_X_cord, tad_Y_cord, col = colorTAD, border = "#ffffff")
  
  ###############################
  ##Obtaining breakpoint position
  ## Opposite logic than for Gene SV paint
  ##If gene is NOT relocated we paint the breakpos at the contrary tad breakpos

  ##Look at which breakpos painted in the gene SV tad, and get the coord for the other
  all_breakpos<-c(info_drawingGENE_TAD$geneTAD_breakP, info_drawingSecondaryTAD$secondaryTAD_breakP)
  
  breakPos<-all_breakpos[all_breakpos != infoDrawing_gene_sv_tad$breakPos]

  ##Defining Breakpoint position
  breakp_x_coord<-c(breakPos,breakPos)##Tad central point for gene broken
  breakp_y_coord<-c(tad_Y_cord[1]-9,##-7 before
                    tad_Y_cord[3]+3)

  #####################
  ##Adding enhancers
  #####################
  enh_y_positions<-rep.int(x=tad_Y_cord[1], times = 4)##c(0,0,0,0)
  
  #########################################
  ## Defining additional relevant variables
  #########################################
  enhNumbersSize<-0.8
  distance_Yaxis_EnhLabel<-2
  distance_Yaxis_EnhLabel_EnhNumber<-1.3
  
  ###########################
  ## Painting the other color Figure
  ##   "otherColorFigure"
  ## The other color covering the TAD coming from the rearranged domain
  ############################

  ##coordinates in bottom-left bottom-right top-right top-left order

  ##########################################################
  ###First Getting The Coordinates of the "Other Figure"

  if(tad_X_cord[1]<20){
    ##It means we are on the LEFT SIDE so the other color goes from the TAD breakpos to the end of the TAD

    ##if breakpos before TAD center, we are going to color an irregular QUADRILATERAl
    ##So we will need 4 coordinates for the figure
    ##If its in the TAD center or after we will color a Triangle, so we will need only 3 coordinates
    if(breakPos>tad_X_cord[3]){##tad_XCoord_OnLeftSide[3],##TAD center position

      ##Right triangle to paint
      otherColorFigure_X_pos<-c(breakPos,
                                tad_X_cord[2],
                                breakPos)

      ###############################################
      ##First getting the triangle corner angle

      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs

      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm

      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians

      ##Now let's compute the height of interest. With respect to the right triangle, so right side
      catetoPlaying<-(tad_X_cord[2]-breakPos)
      heightOfInterest<-tan(cornerAngle)*catetoPlaying

      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest

      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                mostWanted_Y_pos)


    }else if(breakPos==tad_X_cord[3]){

      ##Right triangle to paint
      otherColorFigure_X_pos<-c(breakPos,
                                tad_X_cord[2],
                                breakPos)

      ##In this case very easy the most wanted Y coordinate, is just the TAD peak
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                tad_Y_cord[3])

    }else if(breakPos<tad_X_cord[3]){
      ##So QUADRILATERAL to paint
      otherColorFigure_X_pos<-c(breakPos,
                                tad_X_cord[2],
                                tad_X_cord[3],
                                breakPos)

      ###############################################
      ##First getting the triangle corner angle

      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs

      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm

      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians

      ##Now let's compute the height of interest
      ##We compute it from the perspective of the right triangle!!! so Left side of the TAD
      catetoPlaying<-(breakPos-tad_X_cord[1])
      heightOfInterest<-tan(cornerAngle)*catetoPlaying

      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest

      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                tad_Y_cord[3],
                                mostWanted_Y_pos
      )

    }


  }else if(tad_X_cord[1]>20){
    ##It means we are on the RIGHT SIDE so the other color goes from the TAD start to the breakpos

    ##if breakpos after TAD center, we are going to color an irregular QUADRILATERAl
    ##So we will need 4 coordinates for the figure
    ##If its in the TAD center or before we will color a Triangle, so we will need only 3 coordinates
    if(breakPos>tad_X_cord[3]){##tad_XCoord_OnLeftSide[3],##TAD center position

      ##So QUADRILATERAL to paint
      otherColorFigure_X_pos<-c(tad_X_cord[1],
                                breakPos,
                                breakPos,
                                tad_X_cord[3])

      ###############################################
      ##First getting the triangle corner angle

      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs

      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm

      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians

      ##Now let's compute the height of interest
      catetoPlaying<-(tad_X_cord[2]-breakPos)
      heightOfInterest<-tan(cornerAngle)*catetoPlaying

      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest

      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                mostWanted_Y_pos,
                                tad_Y_cord[3])

    }else if(breakPos==tad_X_cord[3]){
      ##Right triangle to paint
      otherColorFigure_X_pos<-c(tad_X_cord[1],
                                breakPos,
                                breakPos)

      ##In this case very easy the most wanted Y coordinate, is just the TAD peak
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                tad_Y_cord[3])


    }else if(breakPos<tad_X_cord[3]){
      ##Right triangle to paint
      otherColorFigure_X_pos<-c(tad_X_cord[1],
                                breakPos,
                                breakPos)

      ###############################################
      ##First getting the triangle corner angle

      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs

      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm

      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians

      ##Now let's compute the height of interest. With respect to the right triangle
      catetoPlaying<-(breakPos-tad_X_cord[1])
      heightOfInterest<-tan(cornerAngle)*catetoPlaying

      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest

      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                mostWanted_Y_pos)
    }

  }

  ##########################################################
  ### Painting the "Other Figure"

  color_otherFigure<-tadColors[!(tadColors==colorTAD)]##the opposite of the TAD color

  polygon(otherColorFigure_X_pos, otherColorFigure_Y_pos, col = color_otherFigure, border = "#ffffff")

  #######################
  ##Painting Breakpoint
  breakpointColor<-"brown3"
  lines(x=breakp_x_coord,
        y=breakp_y_coord,
        col=breakpointColor,
        lwd=2,
        lty=3)

  #Breakpoint label
  ##Boxed labels: xpad,ypad arguments, used to expand the box beyond the text
  boxed.labels(x=breakp_x_coord[1],  y=breakp_y_coord[1] - 1,labels = "Breakpoint", cex=0.8, border = NA, bg ="white",
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col=breakpointColor
  )

  #TAD label
  # text(x=tad_X_cord[3], y=tad_Y_cord[1]+10, label="TAD", cex = 1.5)
  boxed.labels(x=tad_X_cord[3],  y=tad_Y_cord[1]+10,labels = "TAD", cex=1.5,
               border = NA, bg= NA, xpad=1, ypad = 2 )

  ###Adding black line below TAD
  lines(x = c(tad_X_cord[1]-3,tad_X_cord[2]+3),
        y=c(tad_Y_cord[1],tad_Y_cord[1]),
        col="black", lwd=2, lty=1)
  
  
  ###############################################
  ##Let's compute and plot the enhancers
  ###############################################
  ##nEnhGeneLost
  ##nEnhOtherDomainNotGained
  ##locationGeneEnhVsOtherDomain<-"GeneEnhLost_Breakp_EnhOtherDomainNotGained" vs "EnhOtherDomainNotGained_Breakp_GeneEnhLost"
  ##To check whether they are painted to the right
  ##Or to the left attending to the breakpoint location. It determines where the other number goes

  nEnhGeneLost<-(nEnh_initial_left+nEnh_initial_right)-(nEnh_kept_left+nEnh_kept_right) ##If it is 0 NO enh lost
  nEnhOtherDomainNotGained<-nEnh_other_domain-nEnh_gained ##if it is 0 NO enh not captured from the other TAD
    
  ##Getting the values for these variables  
  if(situation=="primaryTAD_Sinistral"){
    locationGeneEnhVsOtherDomain<-"GeneEnhLost_Breakp_EnhOtherDomainNotGained"
  }else if(situation=="primaryTAD_Dextral"){
    locationGeneEnhVsOtherDomain<-"EnhOtherDomainNotGained_Breakp_GeneEnhLost"
  }
  
  ##############################
  ##Lets paint the enhGeneLost
  ##############################
  if(nEnhGeneLost > 0){
    ##So at least 1 enh lost
    if(locationGeneEnhVsOtherDomain=="GeneEnhLost_Breakp_EnhOtherDomainNotGained"){
      ##So enh painted towards the left side of the breakpoint
      ##breakPos-whtvr
      
      if(nEnhGeneLost==(nEnh_initial_left+nEnh_initial_right)){
        ##So all enhancers are lost, hence we paint 4 enhancers
        if(info_drawingSecondaryTAD$gene_relocated==TRUE){
          enhXpos<-breakPos-1##Maybe need to touch this so that enh do not overlap. Or be closer 
        }else{
          enhXpos<-breakPos-2##Maybe need to touch this so that enh do not overlap. Or be closer 
        }
        enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)

      }else{
        ##So, the gene has lost some enh but not all, so we paint 2 enhancers
        if(info_drawingSecondaryTAD$gene_relocated==TRUE){
          enhXpos<-breakPos-0.5##If gene relocated, this enh stay where they are
        }else{
          enhXpos<-breakPos-2 
        }
        enh_x_positions<-c(enhXpos,enhXpos+0.3)
      }
      
    }else if(locationGeneEnhVsOtherDomain=="EnhOtherDomainNotGained_Breakp_GeneEnhLost"){
      if(nEnhGeneLost==(nEnh_initial_left+nEnh_initial_right)){
        ##So all enhancers are lost, hence we paint 4 enhancers
        if(info_drawingSecondaryTAD$gene_relocated==TRUE){
          enhXpos<-breakPos+0.1 ##this enh stay where they are
        }else{
          enhXpos<-breakPos+2
        }
        enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
        
      }else{
        ##So, the gene has lost some enh but not all, so we paint 2 enhancers
        if(info_drawingSecondaryTAD$gene_relocated==TRUE){
          enhXpos<-breakPos+0.1 ##this enh stay where they are
        }else{
          enhXpos<-breakPos+2
        }
        enh_x_positions<-c(enhXpos,enhXpos+0.3)
      }
    }
    
    ####################################################
    ## Pintamos los Enh
    ##
    draw.ellipse(x=enh_x_positions,
                 y=enh_y_positions,a=0.1,col="#b2d235")
    
    #enhancers Label
    #Below the level of the other enhancers, to avoid overlaps
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                 -distance_Yaxis_EnhLabel,
                 labels = "enhancers", cex=0.8, 
                 border = NA, bg ="white", 
                 xpad=1,
                 ypad=1 #To allow the text to breath
    )
    #Nr of enhancers
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                 -distance_Yaxis_EnhLabel
                 -distance_Yaxis_EnhLabel_EnhNumber,
                 labels = paste0("n=",nEnhGeneLost), cex=enhNumbersSize, 
                 border = NA, bg ="white", 
                 xpad=1.2,
                 ypad=1 #To allow the text to breath
    )
    
    # ##enhancers horizontal line over them
    # lines(x = c(enh_x_positions[1]-0.1,
    #             enh_x_positions[length(enh_x_positions)]+0.1),
    #       y = c(enh_y_positions[1]+2,enh_y_positions[1]+2),
    #       lwd=1.5)
    
    
    enhHorizontalLineCoord<-c(enh_x_positions[1]-0.1,
                              enh_x_positions[length(enh_x_positions)]+0.1)
    
    enhHorizLineCenter<-enhHorizontalLineCoord[1]+(enhHorizontalLineCoord[2]-enhHorizontalLineCoord[1])/2
    
    
    ##Painting curve connecting enh with gene
    #The curve will end just before the TAD limit, and all the enh lines from this TAD will end there
    ##Here geneCenter represents the position of the gene in the re-arranged TAD
    geneCenter<-infoDrawing_gene_sv_tad$geneCenter
    
    # #Not painting curved arrows for now for the other domain 
    # if(enhXpos > geneCenter){
    #   ##So, the gene is on the left side
    #   curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(tad_X_cord[3]-3,enh_y_positions[1]+6),
    #               segment = c(0,1),arr.pos = 0, arr.type="T", curve = 0.1)
    #   
    #   ##PAINTING VERTICAL LINE ON OWR OWN FUCKIN CURVEDARROW rising problems//bugs
    #   ##So arr.pos=0 to be painted in the enh horizontal line and hence disappear
    #   lines(x=c(tad_X_cord[3]-3,tad_X_cord[3]-3),
    #         y=c(enh_y_positions[1]+5,enh_y_positions[1]+7),
    #         lwd=2)
    #   
    # }else{
    #   ##So the secondary TAD is on the left side, adjust row curvature, to mantain it over the X axis
    #   curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(tad_X_cord[3]+3,enh_y_positions[1]+6),
    #               segment = c(0,1),arr.pos = 0, arr.type="T", curve = -0.1)
    #   
    #   ##PAINTING VERTICAL LINE ON OWR OWN FUCKIN CURVEDARROW rising problems//bugs
    #   ##So arr.pos=0 to be painted in the enh horizontal line and hence disappear
    #   lines(x=c(tad_X_cord[3]+3,tad_X_cord[3]+3),
    #         y=c(enh_y_positions[1]+5,enh_y_positions[1]+7),
    #         lwd=2)
    #   
    # }
    
    ###################
    
  }
  
  
  #########################################
  ##Lets paint the EnhOtherDomainNotGained
  #########################################
  if(nEnhOtherDomainNotGained > 0){
    ##So at least 1 enh lost
    if(locationGeneEnhVsOtherDomain=="EnhOtherDomainNotGained_Breakp_GeneEnhLost"){
      ##So enh painted towards the left side of the breakpoint
      ##breakPos-whtvr
      
      if(nEnhOtherDomainNotGained==nEnh_other_domain){
        ##So all enhancers are lost, hence we paint 4 enhancers
        if(info_drawingSecondaryTAD$gene_relocated==FALSE){
          ##So this enh stay where they are
          enhXpos<-breakPos-1##Maybe need to touch this so that enh do not overlap. Or be closer
        }else{
          enhXpos<-breakPos-2##bring them away from border
        }
        
        enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
        
      }else{
        ##So, the gene has lost some enh but not all, so we paint 2 enhancers
        if(info_drawingSecondaryTAD$gene_relocated==FALSE){
          ##So this enh stay where they are
          enhXpos<-breakPos-0.5##Maybe need to touch this so that enh do not overlap. Or be closer
        }else{
          enhXpos<-breakPos-2##bring them away from border
        }
        
        enh_x_positions<-c(enhXpos,enhXpos+0.3)
      }
      
    }else if(locationGeneEnhVsOtherDomain=="GeneEnhLost_Breakp_EnhOtherDomainNotGained"){
      if(nEnhOtherDomainNotGained==nEnh_other_domain){
        ##So all enhancers are lost, hence we paint 4 enhancers
         
        if(info_drawingSecondaryTAD$gene_relocated==FALSE){
          ##So this enh stay where they are
          enhXpos<-breakPos+0.1##Maybe need to touch this so that enh do not overlap. Or be closer
        }else{
          enhXpos<-breakPos+2##bring them away from breakPos
        }
        
        enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
        
      }else{
        ##So, the gene has lost some enh but not all, so we paint 2 enhancers
        if(info_drawingSecondaryTAD$gene_relocated==FALSE){
          ##So this enh stay where they are
          enhXpos<-breakPos+0.1##Maybe need to touch this so that enh do not overlap. Or be closer
        }else{
          enhXpos<-breakPos+2##bring them away from breakPos
        }
        enh_x_positions<-c(enhXpos,enhXpos+0.3)
      }
    }
    
    ####################################################
    ## Pintamos los Enh
    ##
    draw.ellipse(x=enh_x_positions,
                 y=enh_y_positions,a=0.1,col="#b2d235")
    
    #enhancers Label
    #Below the level of the other enhancers, to avoid overlaps
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                 -distance_Yaxis_EnhLabel
                 -distance_Yaxis_EnhLabel_EnhNumber
                 -distance_Yaxis_EnhLabel_EnhNumber,
                 labels = "enhancers", cex=0.8, 
                 border = NA, bg ="white", 
                 xpad=1,
                 ypad=1 #To allow the text to breath
    )
    #Nr of enhancers
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                 -distance_Yaxis_EnhLabel
                 -distance_Yaxis_EnhLabel_EnhNumber
                 -distance_Yaxis_EnhLabel_EnhNumber
                 -distance_Yaxis_EnhLabel_EnhNumber,
                 labels = paste0("n=",nEnhOtherDomainNotGained), cex=enhNumbersSize, 
                 border = NA, bg ="white", 
                 xpad=1.2,
                 ypad=1 #To allow the text to breath
    )
    
    # ##enhancers horizontal line over them
    # lines(x = c(enh_x_positions[1]-0.1,
    #             enh_x_positions[length(enh_x_positions)]+0.1),
    #       y = c(enh_y_positions[1]+2,enh_y_positions[1]+2),
    #       lwd=1.5)
    
    
    enhHorizontalLineCoord<-c(enh_x_positions[1]-0.1,
                              enh_x_positions[length(enh_x_positions)]+0.1)
    
    enhHorizLineCenter<-enhHorizontalLineCoord[1]+(enhHorizontalLineCoord[2]-enhHorizontalLineCoord[1])/2
    
    
    ##Painting curve connecting enh with gene
    #The curve will end just before the TAD limit, and all the enh lines from this TAD will end there
    ##Here geneCenter represents the position of the gene in the re-arranged TAD
    geneCenter<-infoDrawing_gene_sv_tad$geneCenter
    
    ##We paint the line up to 3 points from the breakpoint
    
    # #Not painting for now
    # if(enhXpos > geneCenter){
    #   ##So, the gene is on the left side
    #   curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(tad_X_cord[3]-3,enh_y_positions[1]+6),
    #               segment = c(0,1),arr.pos = 0, arr.type="T", curve = 0.1)
    #   
    #   ##PAINTING VERTICAL LINE ON OWR OWN FUCKIN CURVEDARROW rising problems//bugs
    #   ##So arr.pos=0 to be painted in the enh horizontal line and hence disappear
    #   lines(x=c(tad_X_cord[3]-3,tad_X_cord[3]-3),
    #         y=c(enh_y_positions[1]+5,enh_y_positions[1]+7),
    #         lwd=2)
    #   
    # }else{
    #   ##So the secondary TAD is on the left side, adjust row curvature, to mantain it over the X axis
    #   curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(tad_X_cord[3]+3,enh_y_positions[1]+6),
    #               segment = c(0,1), curve=-0.1, ##en lugar de -1
    #               arr.pos = 0, arr.type="T")
    #   
    #   ##PAINTING VERTICAL LINE ON OWR OWN FUCKIN CURVEDARROW rising problems//bugs
    #   ##So arr.pos=0 to be painted in the enh horizontal line and hence disappear
    #   lines(x=c(tad_X_cord[3]+3,tad_X_cord[3]+3),
    #         y=c(enh_y_positions[1]+5,enh_y_positions[1]+7),
    #         lwd=2)
    # }
    
  }
  
}

########################
## Paint GeneTruncation
########################
##No TADs painted just gene disrupted
paint_Gene_Truncation<-function(gene, xAxisLim, tad_X_cord, tad_YCoord_Rearrangements){
  ##
  
  tad_Y_cord<-tad_YCoord_Rearrangements
  ##Move it a bit up
  tad_Y_cord<-tad_Y_cord + 5
  
  ##On the center of the x axisLim painting the gene. In size *2
  #gene position at the center -0.5 (since we provide the left coord over which the gene is painted)
  geneSize<-2
  geneConservedFraction<-1/3
  
  ##Gene LeftPos
  genePos<-(xAxisLim[1] + (xAxisLim[2] - xAxisLim[1])/2) - (geneSize/2)##o coger el centro del TAD con TAD3 y restarle 0.5 ya que mide 1
  
  
  ###Adding black line below Gene
  lines(x = c(tad_X_cord[1]-3,tad_X_cord[2]+3),
        y=c(tad_Y_cord[1],tad_Y_cord[1]),
        col="black", lwd=2, lty=1)
  
  
  ##Painting Gene
  gene_X_cord<-c(genePos, genePos+geneSize, genePos+geneSize, genePos)
  gene_Y_cord<-c(tad_Y_cord[1]-2,tad_Y_cord[1]-2,tad_Y_cord[1]+2,tad_Y_cord[1]+2)
  polygon(gene_X_cord, gene_Y_cord, col = "#0e3d61")
  
  ##Adding white box in between. To simulate gene truncated
  whiteBox_X_cord<-c(genePos + geneSize*geneConservedFraction, 
                     genePos + geneSize - geneSize*geneConservedFraction ,
                     genePos + geneSize - geneSize*geneConservedFraction,
                     genePos + geneSize*geneConservedFraction)
  
  whiteBox_Y_cord<-gene_Y_cord
  polygon(whiteBox_X_cord, whiteBox_Y_cord, col = "#ffffff", border = "#ffffff" )
  
  ##Adding Cross Over White Box
  ##Line 1 of the cross
  lines(x =  c(whiteBox_X_cord[1], whiteBox_X_cord[2]),
        y = c(whiteBox_Y_cord[3],whiteBox_Y_cord[1]),
        lwd=2,
        col="#e60000")
  
  ##Line 2 of the cross
  lines(x =  c(whiteBox_X_cord[1], whiteBox_X_cord[2]),
        y = c(whiteBox_Y_cord[1], whiteBox_Y_cord[3]),
        lwd=2,
        col="#e60000")
  
  
  ##Adding Gene Label, above gene
  ##Boxed labels: xpad,ypad arguments, used to expand the box beyond the text
  boxed.labels(x=genePos+geneSize/2,  y=tad_Y_cord[1]+4,labels = gene, cex=1.2, border = NA, bg ="white", 
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col="black"
  )
  
  
  
  ##Adding Gene Truncated Label, below gene
  ##Boxed labels: xpad,ypad arguments, used to expand the box beyond the text
  boxed.labels(x=genePos+geneSize/2,  y=tad_Y_cord[1]-4,labels = "Truncated", cex=1.2, border = NA, bg ="white", 
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col="black"
  )
  
}



########################
## Paint Gene Deletion
########################
paint_Gene_Deletion<-function(gene, xAxisLim, tad_X_cord, tad_YCoord_Rearrangements){
  ##
  
  tad_Y_cord<-tad_YCoord_Rearrangements
  ##Move it a bit up
  tad_Y_cord<-tad_Y_cord + 5
  
  ##On the center of the x axisLim painting the gene. In size *2
  #gene position at the center -0.5 (since we provide the left coord over which the gene is painted)
  geneSize<-2
  geneConservedFraction<-1/3
  
  ##Gene LeftPos
  genePos<-(xAxisLim[1] + (xAxisLim[2] - xAxisLim[1])/2) - (geneSize/2)##o coger el centro del TAD con TAD3 y restarle 0.5 ya que mide 1
  
  
  ###Adding black line below Gene
  lines(x = c(tad_X_cord[1]-3,tad_X_cord[2]+3),
        y=c(tad_Y_cord[1],tad_Y_cord[1]),
        col="black", lwd=2, lty=1)
  
  
  ##Painting Gene
  gene_X_cord<-c(genePos, genePos+geneSize, genePos+geneSize, genePos)
  gene_Y_cord<-c(tad_Y_cord[1]-2,tad_Y_cord[1]-2,tad_Y_cord[1]+2,tad_Y_cord[1]+2)
  polygon(gene_X_cord, gene_Y_cord, col = "#0e3d61")
  
  # ##Adding white box in between. To simulate gene truncated
  # whiteBox_X_cord<-c(genePos + geneSize*geneConservedFraction, 
  #                    genePos + geneSize - geneSize*geneConservedFraction ,
  #                    genePos + geneSize - geneSize*geneConservedFraction,
  #                    genePos + geneSize*geneConservedFraction)
  # 
  # whiteBox_Y_cord<-gene_Y_cord
  # polygon(whiteBox_X_cord, whiteBox_Y_cord, col = "#ffffff", border = "#ffffff" )
  
  ##Adding Cross Over Gene
  ##Line 1 of the cross
  polygon(x =  c(gene_X_cord[1]-0.5, gene_X_cord[2]+0.5),
        y = c(gene_Y_cord[3]+1,gene_Y_cord[1]-1),
        lwd=2,
        col="#e60000",
        border = "#e60000")
  
  ##Line 2 of the cross
  polygon(x =  c(gene_X_cord[1]-0.5, gene_X_cord[2]+0.5),
        y = c(gene_Y_cord[1]-1, gene_Y_cord[3]+1),
        lwd=2,
        col="#e60000",
        border = "#e60000")
  
  
  ##Adding Gene Label, above gene
  ##Boxed labels: xpad,ypad arguments, used to expand the box beyond the text
  boxed.labels(x=genePos+geneSize/2,  y=tad_Y_cord[1]+5,labels = gene, cex=1.2, border = NA, bg ="white", 
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col="black"
  )
  
  
  
  ##Adding Gene Truncated Label, below gene
  ##Boxed labels: xpad,ypad arguments, used to expand the box beyond the text
  boxed.labels(x=genePos+geneSize/2,  y=tad_Y_cord[1]-5,labels = "Deleted", cex=1.2, border = NA, bg ="white", 
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col="black"
  )
  
}


#######################################
## Paint Rearranged TAD on Deletions
#######################################
paintGene_SV_Deletion_TAD<-function(nEnh_initial_left, nEnh_initial_right, gene, gene_breakp_line_type, situation,
                           nEnh_kept_left, nEnh_kept_right,nEnh_gained,
                           nEnh_other_domain,
                           info_drawingGENE_TAD, info_drawingSecondaryTAD,
                           tad_X_cord,
                           tad_YCoord_Rearrangements,
                           geneBreakP_Position_respectToTSS){
  
  ########################################################
  ## Function to paint GENE SV re-arranged TAD 
  ########################################################
  
  ## When long-range && Deletion it will always be just one TAD be painted.
  ## Because either intraTAD enh deletion, or TAD fusion
  ## But only one resulting TAD
  
  gene_sv_TAD_x_coord<-tad_X_cord
  
  ###########################################
  # Painting main scaffold TAD 
  
  ### Y axis position of the Inversion WT row
  tad_Y_cord<-tad_YCoord_Rearrangements
  
  
  ## Let's mantain color codes by TAD position
  ## If TAD is on the left blue, if it is on the right orange
  ## Regardless whether it is the gene or the secondary domain TAD
  ## In order not to confound the user if > 1 gene is affected
  
  tad_X_cord<-gene_sv_TAD_x_coord
  
  # 
  ##Here we can not use tad_X_coord because in all cases only one tad painted, hence tad_X_coord the same
  if(situation == "primaryTAD_Central"){
    colorTAD<-"#9999ff"
    color_otherFigure<-"#9999ff" ##color of the other TAD half
    
  }else if(situation == "primaryTAD_Sinistral"){
    colorTAD<-"#acdfeb"
    color_otherFigure<-"#e5edb2"  ##color of the other TAD half
    
  }else if(situation == "primaryTAD_Dextral"){
    colorTAD<-"#e5edb2"
    color_otherFigure<-"#acdfeb"
  }
  
  #Obsolete see how to update
  #tadColors<-c("#acdfeb","#e5edb2") ##Esto sera para elegir el color complementario

  ###################
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
  
  ################################################################
  ## Defining additional relevant variables
  ## A lot of the following code from:  ##From paintGene_WT_TAD()
  ################################################################
  enhNumbersSize<-0.8
  
  distanceBreakpFromGene<-0.3 ##Distance breakpoint from the gene, when the breakpoint is between the gene & enh
  
  distanceEnhFromGene<-1.5 ##Distance between the enhancers and the gene
  
  distanceBreakpFromEnh<-0.3 ##Distance breakpoint from enhancers when the breakpoint located between enh and TAD border
  
  sizeEnhRegion<-1##Touching this will not really affect enh size, because created regardless of this
  
  distance_Yaxis_geneLabel<-2
  distance_Yaxis_EnhGene<-2
  distance_Yaxis_EnhLabel_EnhNumber<-1.3
  
  
  ###############################
  ##Obtaining breakpoint position
  ##If gene is relocated we paint the breakpos at the contrary tad breakpos
  ##If not, we paint its breakpos
  ##This comes from inv & transloc, for del there will be no relocation
  #But I keep it for now
  #No va a haber relocation, pq la zona de relocation es la zona de deleccion
  # if(info_drawingSecondaryTAD$gene_relocated==TRUE){
  #   breakPos<-info_drawingSecondaryTAD$secondaryTAD_breakP
  # }else{
  #   ##So the gene stays on its TAD
  #   breakPos<-info_drawingGENE_TAD$geneTAD_breakP
  # }
  
  
  
  if(gene_breakp_line_type=="center"){
    ##This will never occur because gene disrupted so other plotting strategy
    ##But I will leave it for now
    
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
  
  ##Defining Breakpoint position
  breakp_x_coord<-c(breakPos,breakPos)##Tad central point for gene broken
  breakp_y_coord<-c(tad_Y_cord[1]-9,##-7 before
                    tad_Y_cord[3]+3)
  
  
  ###########################
  ## Painting the other color Figure
  ##   "otherColorFigure"
  ## The other color covering the TAD coming from the rearranged domain
  ############################
  
  ##coordinates in bottom-left bottom-right top-right top-left order
  
  ##########################################################
  ###First Getting The Coordinates of the "Other Figure"
  ##########################################################
  ##If situation == "primaryTAD_Central", so IntraTAD deletion (or completely outside), no color mixture, just the same TAD color painted...
  
  if(situation=="primaryTAD_Sinistral"){
    ##It means we are on the LEFT SIDE so the other color goes from the TAD breakpos to the end of the TAD
    
    ##if breakpos before TAD center, we are going to color an irregular QUADRILATERAl
    ##So we will need 4 coordinates for the figure
    ##If its in the TAD center or after we will color a Triangle, so we will need only 3 coordinates
    if(breakPos>tad_X_cord[3]){##tad_XCoord_OnLeftSide[3],##TAD center position
      
      ##Right triangle to paint
      otherColorFigure_X_pos<-c(breakPos,
                                tad_X_cord[2],
                                breakPos)
      
      ###############################################
      ##First getting the triangle corner angle
      
      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs
      
      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm
      
      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians
      
      ##Now let's compute the height of interest. With respect to the right triangle, so right side
      catetoPlaying<-(tad_X_cord[2]-breakPos)
      heightOfInterest<-tan(cornerAngle)*catetoPlaying
      
      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest
      
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                mostWanted_Y_pos)
      
      
    }else if(breakPos==tad_X_cord[3]){
      
      ##Right triangle to paint
      otherColorFigure_X_pos<-c(breakPos,
                                tad_X_cord[2],
                                breakPos)
      
      ##In this case very easy the most wanted Y coordinate, is just the TAD peak
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                tad_Y_cord[3])
      
    }else if(breakPos<tad_X_cord[3]){
      ##So QUADRILATERAL to paint
      otherColorFigure_X_pos<-c(breakPos,
                                tad_X_cord[2],
                                tad_X_cord[3],
                                breakPos)
      
      ###############################################
      ##First getting the triangle corner angle
      
      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs
      
      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm
      
      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians
      
      ##Now let's compute the height of interest
      ##We compute it from the perspective of the right triangle!!! so Left side of the TAD
      catetoPlaying<-(breakPos-tad_X_cord[1])
      heightOfInterest<-tan(cornerAngle)*catetoPlaying
      
      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest
      
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                tad_Y_cord[3],
                                mostWanted_Y_pos
      )
      
    }
    
    
  }else if(situation=="primaryTAD_Dextral"){
    ##It means we are on the RIGHT SIDE so the other color goes from the TAD start to the breakpos
    
    ##if breakpos after TAD center, we are going to color an irregular QUADRILATERAl
    ##So we will need 4 coordinates for the figure
    ##If its in the TAD center or before we will color a Triangle, so we will need only 3 coordinates
    if(breakPos>tad_X_cord[3]){##tad_XCoord_OnLeftSide[3],##TAD center position
      
      ##So QUADRILATERAL to paint
      otherColorFigure_X_pos<-c(tad_X_cord[1],
                                breakPos,
                                breakPos,
                                tad_X_cord[3])
      
      ###############################################
      ##First getting the triangle corner angle
      
      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs
      
      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm
      
      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians
      
      ##Now let's compute the height of interest
      catetoPlaying<-(tad_X_cord[2]-breakPos)
      heightOfInterest<-tan(cornerAngle)*catetoPlaying
      
      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest
      
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                mostWanted_Y_pos,
                                tad_Y_cord[3])
      
    }else if(breakPos==tad_X_cord[3]){
      ##Right triangle to paint
      otherColorFigure_X_pos<-c(tad_X_cord[1],
                                breakPos,
                                breakPos)
      
      ##In this case very easy the most wanted Y coordinate, is just the TAD peak
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                tad_Y_cord[3])
      
      
    }else if(breakPos<tad_X_cord[3]){
      ##Right triangle to paint
      otherColorFigure_X_pos<-c(tad_X_cord[1],
                                breakPos,
                                breakPos)
      
      ###############################################
      ##First getting the triangle corner angle
      
      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs
      
      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm
      
      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians
      
      ##Now let's compute the height of interest. With respect to the right triangle
      catetoPlaying<-(breakPos-tad_X_cord[1])
      heightOfInterest<-tan(cornerAngle)*catetoPlaying
      
      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest
      
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                mostWanted_Y_pos)
    }
    
  }
  
  ##########################################################
  ### Painting the "Other Figure"
  ##Recall only paint if situation different than primaryTadCentral
  ##Because if situation primaryTadCentral, it means intraTAD sv or completely deleted tad
  ##So we do not mix colors
  if(situation %in% c("primaryTAD_Sinistral","primaryTAD_Dextral")){
    polygon(otherColorFigure_X_pos, otherColorFigure_Y_pos, col = color_otherFigure, border = "#ffffff")
  }
  
  #######################
  ##Painting Breakpoint
  breakpointColor<-"brown3"
  lines(x=breakp_x_coord,
        y=breakp_y_coord,
        col=breakpointColor,
        lwd=2,
        lty=3)
  
  #Breakpoint label
  ##Boxed labels: xpad,ypad arguments, used to expand the box beyond the text
  boxed.labels(x=breakp_x_coord[1],  y=breakp_y_coord[1] - 1,labels = "Breakpoint", cex=0.8, border = NA, bg ="white", 
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col=breakpointColor
  )
  
  #TAD label
  # text(x=tad_X_cord[3], y=tad_Y_cord[1]+10, label="TAD", cex = 1.5)
  boxed.labels(x=tad_X_cord[3],  y=tad_Y_cord[1]+10,labels = "TAD", cex=1.5,
               border = NA, bg= NA, xpad=1, ypad = 2 )
  
  ###Adding black line below TAD
  lines(x = c(tad_X_cord[1]-3,tad_X_cord[2]+3),
        y=c(tad_Y_cord[1],tad_Y_cord[1]), 
        col="black", lwd=2, lty=1)
  
  #########################################
  ##Defining some profitable variables
  #########################################
  
  distanceEnhFromGene<-1.5 ##Distance between the enhancers and the gene
  
  sizeEnhRegion<-1##Touching this will not really affect enh size, because created regardless of this
  
  ###############################
  ## Painting gene body 
  ###############################
  
  ##genePos computed before in the function
  
  
  geneCenter<-genePos + geneSize/2 ##Used to direct the arrows to here
  gene_X_cord<-c(genePos,genePos+1,genePos+1,genePos)
  gene_Y_cord<-c(tad_Y_cord[1]-1,tad_Y_cord[1]-1,tad_Y_cord[1]+1,tad_Y_cord[1]+1)
  
  ###################
  #Adding gene body
  polygon(gene_X_cord, gene_Y_cord, col = "#0e3d61")
  
  #gene Label
  ##I want it to overlay the breakpoint line
  # https://stackoverflow.com/questions/45366243/text-labels-with-background-colour-in-r
  boxed.labels(x=geneCenter,  y=tad_Y_cord[1]-distance_Yaxis_geneLabel,
               labels = gene, cex=0.8, border = NA, bg ="white", 
               xpad=1,
               ypad=2##To allow the text to breathe
  )
  
  
  #####################
  ##Adding enhancers
  #####################
  enh_y_positions<-rep.int(x=tad_Y_cord[1], times = 4)##c(0,0,0,0)
  
  ##################################################
  ## First, PAINTING COGNATE ENHANCERS
  ## changes on enh originally surrounding the gene
  ##################################################
  
  #If nkept, different than nInitial we will paint only two enh, to represent
  #a reduction on its number
  ####Seguir por aqui morning, por linea 2238
  
  ##Gene in deletion not relocated because being in the relocation area
  ##would imply deletion, and hence would be treated as removed... 
  ##That's why I have removed an upcoming if
  
  # if(info_drawingSecondaryTAD$gene_relocated=="FALSE"){
    
  #If nkept, different than nInitial, as the breakpoint in the middle, we will paint only two enh, to represent
  #a reduction on its number.
  
  ####################################
  ####To the Right side gene Position
  ####################################
  
  if((geneBreakP_Position_respectToTSS=="beforeTSS")||(gene_breakp_line_type=="afterTSS_removing_none")){
    ##Enhancers to  the right left intact, same as they were
    
    if(nEnh_initial_right>0){
      ##So there is at least 1 enh
      enhXpos<-genePos+distanceEnhFromGene##left margin of the enh cluster
      enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
      
      draw.ellipse(x=enh_x_positions,
                   y=enh_y_positions,a=0.1,col="#b2d235")
      
      #enhancers Label
      #text(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel              -distance_Yaxis_EnhGene, label="enhancers", cex = 0.8)
      boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      #Nr of enhancers
      #text(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel-distance_Yaxis_EnhGene-distance_Yaxis_EnhLabel_EnhNumber, label=paste0("n=",nEnh_initial_right), cex=enhNumbersSize)
      boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",nEnh_initial_right), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
            lwd=1.5)
      
      curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                     ,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T")
    }
    
  }else if((geneBreakP_Position_respectToTSS=="afterTSS") && (nEnh_initial_right>0)){
    
    
    if(gene_breakp_line_type=="afterTSS_removing_all"){
      ##DOING NOTHING, we don't want to paint anything as all the enhancers are going to disappear
    }else if(gene_breakp_line_type=="afterTSS_removing_some"){
      
      ###########################
      ## PAINTING 2 ENHANCERS
      
      ##To quickly give the impression that some enhancers are lost, I'm going to paint 2enhancers
      ##That will be the two right enhn if the gene is on the right part of the plot, or the two left
      ##So there is at least 1 enh
      enhXpos<-genePos+distanceEnhFromGene##left margin of the enh cluster
      
      ##Painting Only 2 enhancers
      enh_x_positions<-c(enhXpos,enhXpos+0.3)
      
      draw.ellipse(x=enh_x_positions,
                   y=enh_y_positions[1:2],a=0.1,col="#b2d235")
      
      #enhancers Label
      boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      #Nr of enhancers
      boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",nEnh_kept_right), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
            lwd=1.5)
      
      curvedarrow(from = c(enhXpos+0.15,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                      ,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T" )
      
      
      ###END OF PAINTING 2 ENHANCERS
      ###################################
      
    }
    
  }
  
  
  
  ###################################
  ####To the Left side gene Position
  ###################################
  
  #the gene is not relocated, so the enh to the left will be kept intact if:
  ##the breakp is after the tss (gene not relocated because it is a Sinistral TAD, located at the left part of the screen)
  ##the breakp is before the tss removing none enh (gene not relocated because it is a Destral TAD, located at the right part of the screen)
  if((geneBreakP_Position_respectToTSS=="afterTSS")||(gene_breakp_line_type=="beforeTSS_removing_none")){
    ##Enhancers to  the LEFT kept intact, same as they were
    
    if(nEnh_initial_left>0){
      enhXpos<-genePos-distanceEnhFromGene##left margin of the enh cluster
      enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
      
      draw.ellipse(x=enh_x_positions,
                   y=enh_y_positions,a=0.1,col="#b2d235")
      
      #enhancers Label
      boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      
      #Nr of enhancers
      boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",nEnh_kept_left), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1]+2,enh_y_positions[1]+2),
            lwd=1.5)
      
      curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter#-0.25
                                                                     ,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T", curve = -1 )
    }
    
  }else if((geneBreakP_Position_respectToTSS=="beforeTSS") && (nEnh_initial_left>0)){
    
    
    if(gene_breakp_line_type=="beforeTSS_removing_all"){
      ##DOING NOTHING, we don't want to paint anything as all the enhancers are going to disappear
    }else if(gene_breakp_line_type=="beforeTSS_removing_some"){
      
      ###########################
      ## PAINTING 2 ENHANCERS
      
      ##To quickly give the impression that some enhancers are lost, I'm going to paint just 2 enhancers
      ##So there is at least 1 enh
      enhXpos<-genePos-distanceEnhFromGene##left margin of the enh cluster
      
      ##Painting Only 2 enhancers
      #enh_x_positions<-c(enhXpos,enhXpos+0.3)
      enh_x_positions<-c(enhXpos+0.6, enhXpos+0.9)##The two second enhancers
      enhXpos<-enh_x_positions[1]
      
      draw.ellipse(x=enh_x_positions,
                   y=enh_y_positions[1:2],a=0.1,col="#b2d235")
      
      #enhancers Label
      boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      #Nr of enhancers
      boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",nEnh_kept_left), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
            lwd=1.5)
      
      enhHorizontalLineCoord<-c(enh_x_positions[1]-0.1,
                                enh_x_positions[length(enh_x_positions)]+0.1)
      
      enhHorizLineCenter<-enhHorizontalLineCoord[1]+(enhHorizontalLineCoord[2]-enhHorizontalLineCoord[1])/2
      
      
      curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                            ,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T", curve=-1 )
      
      
      ###END OF PAINTING 2 ENHANCERS
      ###################################
      
    }
    
  }
  
  ########}###
  
  
  ##################################################
  ## TIME FOR PAINTING ECTOPIC  ENHANCERS
  ## changes on enh originally surrounding the gene
  ##################################################
  ##For intraTAD SV nEnh_gained will be 0 so no painting
  ##If there is no gain, no paint
  ##If nEnh gained == nEnhInitialOther Domain, the four enh are painted
  ##On the contrary only two enh are painted
  
  ##If the gene is not relocated we paint the enh on the middle of its segment
  
  
  if(nEnh_gained>0){
    ##So some enh need to be painted
    
      ##When the gene not relocated, to locate the enh we work with the "other figure" coordinates
      ##Because the ectopic enhancers are the elements relocated now
      ##We paint the enh in the middle of the other Figure coord
      x_space<-c(otherColorFigure_X_pos[1],otherColorFigure_X_pos[2])
      
      x_space_start<-x_space[1]
      x_space_end<-x_space[2]
      
      #gene position at the center -0.5 (since we provide the left coord over which the gene is painted)
      enhXpos<-(x_space_start + (x_space_end-x_space_start)/2)-0.5
      enhXpos<-enhXpos-0.45
      
      if(nEnh_gained==nEnh_other_domain){
        ##All enh from the other Domain Gained so
        #we paint 4 enh
        
        ##As it is the gene the one that is relocated, here we work with tad_X_cords
        enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
        
      }else if(nEnh_gained < nEnh_other_domain){
        ##So we are going to paint only 2 enhancers
        
        ##So we paint the two enh on the left
        enh_x_positions<-c(enhXpos,enhXpos+0.3)
      }
    
    ####################################################
    ## Pintamos los ECTOPIC enhancers
    ##Lo de arriba solo para calcular posicion de enh
    ##get enh positions center
    ##para pintar hacia derecha o izquierda mirar si el geneCenter queda a la izq o a la derecha
    
    ##
    draw.ellipse(x=enh_x_positions,
                 y=enh_y_positions,a=0.1,col="#b2d235")
    
    #enhancers Label
    #Below the level of the other enhancers, to avoid overlaps
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                 -distance_Yaxis_geneLabel
                 -distance_Yaxis_EnhGene
                 -distance_Yaxis_EnhLabel_EnhNumber
                 -distance_Yaxis_EnhLabel_EnhNumber,
                 labels = "enhancers", cex=0.8, 
                 border = NA, bg ="white", 
                 xpad=1,
                 ypad=1 #To allow the text to breath
    )
    #Nr of enhancers
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                 -distance_Yaxis_geneLabel
                 -distance_Yaxis_EnhGene
                 -distance_Yaxis_EnhLabel_EnhNumber
                 -distance_Yaxis_EnhLabel_EnhNumber
                 -distance_Yaxis_EnhLabel_EnhNumber, labels = paste0("n=",nEnh_gained), cex=enhNumbersSize, 
                 border = NA, bg ="white", 
                 xpad=1.2,
                 ypad=1 #To allow the text to breath
    )
    
    ##enhancers horizontal line over them
    lines(x = c(enh_x_positions[1]-0.1,
                enh_x_positions[length(enh_x_positions)]+0.1),
          y = c(enh_y_positions[1]+2,enh_y_positions[1]+2),
          lwd=1.5)
    
    
    enhHorizontalLineCoord<-c(enh_x_positions[1]-0.1,
                              enh_x_positions[length(enh_x_positions)]+0.1)
    
    enhHorizLineCenter<-enhHorizontalLineCoord[1]+(enhHorizontalLineCoord[2]-enhHorizontalLineCoord[1])/2
    
    ##Painting curve connecting enh with gene
    if(enhXpos > geneCenter){
      ##So, the secondary TAD is on the right side
      curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(geneCenter,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T", curve = 1)
    }else{
      ##So the secondary TAD is on the left side, adjust row curvature, to mantain it over the X axis
      curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(geneCenter,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T", curve = -1)
    }
    
  }
  
  ##################################################################################
  ##Adding Triangle Shape over gene, if it has any enh, to represent the arrow end
  ##################################################################################
  if((nEnh_kept_left>0)||(nEnh_kept_right>0)||(nEnh_gained>0)){
    triangle_X_Coord<-c(geneCenter-0.3, geneCenter, geneCenter+0.3)
    heightOrizontalLineEnh<-enh_y_positions[1]+2
    triangle_Y_Coord<-c(heightOrizontalLineEnh+0.3, heightOrizontalLineEnh-0.2, heightOrizontalLineEnh+0.3)
    polygon(triangle_X_Coord, triangle_Y_Coord, col = "black", border = "black")
  }
  
  
  
  
  #############################
  ## return profitable info ###
  #############################
  ##break pos, to easily select the other one
  info_drawing<-list("breakPos"=breakPos,
                     "geneCenter"=geneCenter)
  
  return(info_drawing)
  
}

#########################################################################################################
## Paint Rearranged TAD on Duplications Where Enh Are Involved (hence some LongRange potential effects)
## AND WHEN THE GENE IS NOT DUPLICATED
#########################################################################################################
paintGene_SV_Duplication_OnlyLongRange_Intra_TAD<-function(nEnh_initial_left, nEnh_initial_right, gene, gene_breakp_line_type, situation,
                                                     nEnh_kept_left, nEnh_kept_right,nEnh_gained,
                                                     nEnh_other_domain,
                                                     info_drawingGENE_TAD, info_drawingSecondaryTAD,
                                                     tad_X_cord,
                                                     tad_YCoord_Rearrangements,
                                                     geneBreakP_Position_respectToTSS,
                                                     patientResults){
  
  #############################################################################################################
  ##We are dealing with intraTAD SV in this context (if not intraTAD this function is unreachable in the code)
  
  ########################################################
  ## Function to paint GENE SV re-arranged TAD 
  ########################################################
  
  ## When long-range && Deletion it will always be just one TAD be painted.
  ## Because either intraTAD enh deletion, or TAD fusion
  ## But only one resulting TAD
  
  gene_sv_TAD_x_coord<-tad_X_cord
  
  ###########################################
  # Painting main scaffold TAD 
  
  ### Y axis position of the Inversion WT row
  tad_Y_cord<-tad_YCoord_Rearrangements
  
  
  ## Let's mantain color codes by TAD position
  ## If TAD is on the left blue, if it is on the right orange
  ## Regardless whether it is the gene or the secondary domain TAD
  ## In order not to confound the user if > 1 gene is affected
  
  tad_X_cord<-gene_sv_TAD_x_coord
  
  # 
  ##Here we can not use tad_X_coord because in all cases only one tad painted, hence tad_X_coord the same
  if(situation == "primaryTAD_Central"){
    colorTAD<-"#9999ff"
    color_otherFigure<-"#9999ff" ##color of the other TAD half
    
  }else{
    stop("Do not understand what is going on.")
  }
  
  
  ###################
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
  
  ################################################################
  ## Defining additional relevant variables
  ## A lot of the following code from:  ##From paintGene_WT_TAD()
  ################################################################
  enhNumbersSize<-0.8
  
  distanceBreakpFromGene<-0.3 ##Distance breakpoint from the gene, when the breakpoint is between the gene & enh
  
  distanceEnhFromGene<-1.5 ##Distance between the enhancers leftBorder and the gene
  
  distanceBreakpFromEnh<-0.3 ##Distance breakpoint from enhancers when the breakpoint located between enh and TAD border
  
  sizeEnhRegion<-1##Touching this will not really affect enh size, because created regardless of this
  
  distance_Yaxis_geneLabel<-2
  distance_Yaxis_EnhGene<-2
  distance_Yaxis_EnhLabel_EnhNumber<-1.3
  
  
  ###############################
  ##Obtaining breakpoint position
  ###############################
  ##Igual lo de pintar breakpoints casi que omitible
  ###########################################################################
  ##   ONE TAD IS GOING TO BE PAINTED, so in the current TAD representation
  ##   TWO BREAKPOINTS will be displayed
  ##   It is due to the fact that both ,or none, breakpoint directly touch the TAD
  ##   This is going to occur for Deletions And Duplications
  ###########################################################################
  
  ##Obtained from paintGene_WT_TAD, from the primaryTAD_central code chunks
  #If not primaryTad_Central an error will rise before in the code
  
  if(gene_breakp_line_type=="outOf_RegulatoryDomain"){
    ##Breakpos will be painted out of the regulatory domain... and we will paint two
    #breakPos<-c(tad_X_cord[1]-2, tad_X_cord[2]+2)
    #Will be a gene duplication so this function wont be called
    
  }else if(gene_breakp_line_type=="center"){
    ##This plot won't be happening it is a gene disruption plot
  }else if(gene_breakp_line_type=="surroundingGene"){
    ##One breakpoint before the gene
    ##One breakpoint after the gene
    ##This will also be a gene duplciation... this function will not be called for this. I think
    breakPos<-c(genePos - 0.5,
                genePos + geneSize + 0.5)
    
  }else if(gene_breakp_line_type=="afterTSS_duplic_all"){
    ##so all enh duplicated  
    ##So we want lines surrounding the enh cluster
    ##But enh cluster is going to be expanded, doubled, so put the outer breakopint two times the EnhCluster further away
    ## The breakpoint is located right after the end of the gene
    ##sizeEnhRegion*2 is the added part, taking into account, two enh + placed
    breakPos<-c(gene_X_cord[2] + distanceBreakpFromGene, ##Before Enh
                genePos + distanceEnhFromGene + sizeEnhRegion*2 + distanceBreakpFromEnh)##After Enh
    
  }else if(gene_breakp_line_type=="afterTSS_duplic_some"){
    
    enhXpos<-genePos + distanceEnhFromGene
    breakPos<-c(enhXpos+sizeEnhRegion/2,
                enhXpos+sizeEnhRegion*1.5##two enh will be added to the cluster
                +distanceBreakpFromEnh)##One breakpoint inside of the cluster, the other beyond
    
  }else if(gene_breakp_line_type=="afterTSS_duplic_none"){
    ##This will not occur in a pathogenic case, intratad but in any case leave it as it is
    ##just after enhancers
    positionBreak<-genePos+distanceEnhFromGene+sizeEnhRegion+distanceBreakpFromEnh
    breakPos<-c(positionBreak,
                positionBreak+1)##after enhcluster and beyond
    
  }else if(gene_breakp_line_type=="beforeTSS_duplic_all"){
    # The breakpoint is located right before the begining of the gene
    # genePos-distanceEnhFromGene = start coordinate of the enh cluster
    breakPos<-c(genePos - distanceBreakpFromGene + 0.12, ##between gene and enh cluster ##Adding 0.12 to avoid breakpoint to be directly touching the enh 
                genePos - distanceEnhFromGene ##This gives the ectopic left most coordinate of the cluster
                - sizeEnhRegion##The enh cluster duplicated so put it two times far away
                - distanceBreakpFromEnh)
    
  }else if(gene_breakp_line_type=="beforeTSS_duplic_some"){
    
    enhXpos<-genePos-distanceEnhFromGene##Ectopic left border pos ##Obtained from painting enhancers left side code part. Enh cluster Start
    breakPos<-c(enhXpos-distanceBreakpFromEnh-sizeEnhRegion*0.5,##*0.5 extra because two additional enh will be painted
                enhXpos+0.5)##One breakpoint before the cluster, and the other inside
    
    
  }else if(gene_breakp_line_type=="beforeTSS_duplic_none"){
    ## Breakpoint after TAD start
    #breakPos<-tad_X_cord[1]+1.5 ##1.5 Units after TAD start
    ##just before enhancers
    positionBreak<-genePos-distanceEnhFromGene-distanceBreakpFromEnh
    
    breakPos<-c(positionBreak-1,
                positionBreak)
    
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
  
  ##Estoy combinando el script de aux_functions la de paintWT(parte de primary TAD central) y el de paint_longRangedDeletion
  
  ###Adding black line below TAD
  lines(x = c(tad_X_cord[1]-3,tad_X_cord[2]+3),
        y=c(tad_Y_cord[1],tad_Y_cord[1]), 
        col="black", lwd=2, lty=1)
  
  #########################################
  ##Defining some profitable variables
  #########################################
  
  ###############################
  ## Painting gene body 
  ###############################
  
  ###################
  #Adding gene body
  polygon(gene_X_cord, gene_Y_cord, col = "#0e3d61")
  
  #gene Label
  ##I want it to overlay the breakpoint line
  # https://stackoverflow.com/questions/45366243/text-labels-with-background-colour-in-r
  boxed.labels(x=geneCenter,  y=tad_Y_cord[1]-distance_Yaxis_geneLabel,
               labels = gene, cex=0.8, border = NA, bg ="white", 
               xpad=1,
               ypad=2##To allow the text to breathe
  )
  
  
  #####################
  ##Adding enhancers
  #####################
  enh_y_positions<-rep.int(x=tad_Y_cord[1], times = 4)##c(0,0,0,0)
  
  ##################################################
  ## First, PAINTING COGNATE ENHANCERS
  ## changes on enh originally surrounding the gene
  ##################################################
  
  #If nkept, < than nInitial we will paint only two enh, to represent
  #a reduction on its number

  #If nkept, different than nInitial, as the breakpoint in the middle, we will paint only two enh, to represent
  #a reduction on its number.
  
  
  ##If nkept> than initial (because duplication, in this case intraTAD)
  ##If nKept == ninitial*2 we will paint 4 enhmore
  
  #If nKept>nInitial & nkept<nInitial*2 we will paint two more
  ####################################
  ####To the Right side gene Position
  ####################################
  
  if((geneBreakP_Position_respectToTSS=="beforeTSS")||(gene_breakp_line_type=="afterTSS_duplic_none")){
    ##Enhancers to  the right left intact, same as they were
    
    if(nEnh_initial_right>0){
      ##So there is at least 1 enh
      enhXpos<-genePos+distanceEnhFromGene##left margin of the enh cluster
      enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
      
      draw.ellipse(x=enh_x_positions,
                   y=rep.int(x = enh_y_positions[1],times = length(enh_x_positions)),a=0.1,col="#b2d235")
      
      #enhancers Label
      #text(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel              -distance_Yaxis_EnhGene, label="enhancers", cex = 0.8)
      boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      #Nr of enhancers
      #text(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel-distance_Yaxis_EnhGene-distance_Yaxis_EnhLabel_EnhNumber, label=paste0("n=",nEnh_initial_right), cex=enhNumbersSize)
      boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",nEnh_initial_right), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
            lwd=1.5)
      
      curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                     ,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T")
    }
    
  }else if((geneBreakP_Position_respectToTSS=="afterTSS") && (nEnh_initial_right>0)){
    
    
    if(gene_breakp_line_type=="afterTSS_duplic_all"){
      
      ####################################
      ## PAINTING 4 Additional ENHANCERS
      
      enhXpos<-genePos+distanceEnhFromGene##left margin of the enh cluster
      
      ##Painting 8 enhancers
      enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6,enhXpos+0.9,enhXpos+1.2,enhXpos+1.5,enhXpos+1.8,enhXpos+2.1)
      
      draw.ellipse(x=enh_x_positions,
                   y=rep.int(x = enh_y_positions[1],times = length(enh_x_positions)),a=0.1,col="#b2d235")
      
      #enhancers Label
      boxed.labels(x=enhXpos + (enh_x_positions[length(enh_x_positions)] - enh_x_positions[1])/2, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      #Nr of enhancers
      boxed.labels(x=enhXpos + (enh_x_positions[length(enh_x_positions)] - enh_x_positions[1])/2, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",nEnh_kept_right), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
            lwd=1.5)
      
      curvedarrow(from = c(enhXpos + (enh_x_positions[length(enh_x_positions)] - enh_x_positions[1])/2, enh_y_positions[1]+2), 
                  to = c(geneCenter, enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T" )
      
    }else if(gene_breakp_line_type=="afterTSS_duplic_some"){
      
      ###########################
      ## PAINTING 2 Additional ENHANCERS
      
      enhXpos<-genePos+distanceEnhFromGene##left margin of the enh cluster
      
      ##Painting 6 enhancers
      enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6,enhXpos+0.9,enhXpos+1.2,enhXpos+1.5)
      
      draw.ellipse(x=enh_x_positions,
                   y=rep.int(x = enh_y_positions[1],times = length(enh_x_positions)),a=0.1,col="#b2d235")
      
      #enhancers Label
      boxed.labels(x=enhXpos + (enh_x_positions[length(enh_x_positions)] - enh_x_positions[1])/2, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      #Nr of enhancers
      boxed.labels(x=enhXpos + (enh_x_positions[length(enh_x_positions)] - enh_x_positions[1])/2, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",nEnh_kept_right), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
            lwd=1.5)
      
      curvedarrow(from = c(enhXpos + (enh_x_positions[length(enh_x_positions)] - enh_x_positions[1])/2, enh_y_positions[1]+2), 
                  to = c(geneCenter, enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T" )
      
    }
    
  }
  
  
  
  ###################################
  ####To the Left side gene Position
  ###################################
  
  if((geneBreakP_Position_respectToTSS=="afterTSS")||(gene_breakp_line_type=="beforeTSS_duplic_none")){
    ##Enhancers to  the LEFT kept intact, same as they were
    
    if(nEnh_initial_left>0){
      enhXpos<-genePos-distanceEnhFromGene##left margin of the enh cluster
      enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
      
      draw.ellipse(x=enh_x_positions,
                   y=rep.int(x = enh_y_positions[1],times = length(enh_x_positions)),a=0.1,col="#b2d235")
      
      #enhancers Label
      boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      
      #Nr of enhancers
      boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",nEnh_kept_left), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1]+2,enh_y_positions[1]+2),
            lwd=1.5)
      
      curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter#-0.25
                                                                     ,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T", curve = -1 )
    }
    
  }else if((geneBreakP_Position_respectToTSS=="beforeTSS") && (nEnh_initial_left>0)){
    
    
    if(gene_breakp_line_type=="beforeTSS_duplic_all"){
      
    
      ###########################
      ## PAINTING 8 ENHANCERS
      
      enhXpos<-genePos-distanceEnhFromGene-sizeEnhRegion##left margin of the enh cluster
      
      ##Painting 8
      enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6,enhXpos+0.9,enhXpos+1.2,enhXpos+1.5,enhXpos+1.8,enhXpos+2.1)
      
      enhXpos<-enh_x_positions[1]
      
      draw.ellipse(x=enh_x_positions,
                   y=rep.int(x = enh_y_positions[1],times = length(enh_x_positions)),a=0.1,col="#b2d235")
      
      #enhancers Label
      boxed.labels(x=enhXpos + (enh_x_positions[length(enh_x_positions)] - enh_x_positions[1])/2, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      #Nr of enhancers
      boxed.labels(x=enhXpos + (enh_x_positions[length(enh_x_positions)] - enh_x_positions[1])/2, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",nEnh_kept_left), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
            lwd=1.5)
      
      enhHorizontalLineCoord<-c(enh_x_positions[1]-0.1,
                                enh_x_positions[length(enh_x_positions)]+0.1)
      
      enhHorizLineCenter<-enhHorizontalLineCoord[1]+(enhHorizontalLineCoord[2]-enhHorizontalLineCoord[1])/2
      
      
      curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                            ,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T", curve=-1 )
      
      
    }else if(gene_breakp_line_type=="beforeTSS_duplic_some"){
      
      ###########################
      ## PAINTING 6 ENHANCERS
      
      ##To quickly give the impression that some enhancers are lost, I'm going to paint just 2 enhancers
      ##So there is at least 1 enh
      enhXpos<-genePos-distanceEnhFromGene-sizeEnhRegion*0.5##left margin of the enh cluster
      
      ##Painting 6 enhancers
      enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6,enhXpos+0.9,enhXpos+1.2,enhXpos+1.5)
      
      enhXpos<-enh_x_positions[1]
      
      draw.ellipse(x=enh_x_positions,
                   y=rep.int(x = enh_y_positions[1],times = length(enh_x_positions)),a=0.1,col="#b2d235")
      
      #enhancers Label
      boxed.labels(x=enhXpos + (enh_x_positions[length(enh_x_positions)] - enh_x_positions[1])/2, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      #Nr of enhancers
      boxed.labels(x=enhXpos + (enh_x_positions[length(enh_x_positions)] - enh_x_positions[1])/2, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",nEnh_kept_left), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
            lwd=1.5)
      
      enhHorizontalLineCoord<-c(enh_x_positions[1]-0.1,
                                enh_x_positions[length(enh_x_positions)]+0.1)
      
      enhHorizLineCenter<-enhHorizontalLineCoord[1]+(enhHorizontalLineCoord[2]-enhHorizontalLineCoord[1])/2
      
      
      curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                            ,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T", curve=-1 )
      
      
    }
    
  }
  
  ########}###
  
  
  ##################################################################################
  ##Adding Triangle Shape over gene, if it has any enh, to represent the arrow end
  ##################################################################################
  if((nEnh_kept_left>0)||(nEnh_kept_right>0)||(nEnh_gained>0)){
    triangle_X_Coord<-c(geneCenter-0.3, geneCenter, geneCenter+0.3)
    heightOrizontalLineEnh<-enh_y_positions[1]+2
    triangle_Y_Coord<-c(heightOrizontalLineEnh+0.3, heightOrizontalLineEnh-0.2, heightOrizontalLineEnh+0.3)
    polygon(triangle_X_Coord, triangle_Y_Coord, col = "black", border = "black")
  }
  
  
  
  
  #############################
  ## return profitable info ###
  #############################
  ##break pos, to easily select the other one
  info_drawing<-list("breakPos"=breakPos,
                     "geneCenter"=geneCenter)
  
  return(info_drawing)
  
  
  
  
  
}

###########################################
## Paint GeneDirect Duplication
###########################################
## Pathogenic mechanism by gene duplication
##No TADs painted just gene duplicated, so two gene boxes
paint_Gene_Duplicated<-function(gene, xAxisLim, tad_X_cord, tad_YCoord_Rearrangements){
  ##
  
  tad_Y_cord<-tad_YCoord_Rearrangements
  ##Move it a bit up
  tad_Y_cord<-tad_Y_cord + 5
  
  ##Genes position, simmetric with respect to the center
  ##On the center of the x axisLim painting the gene. In size *2
  #gene position at the center -0.5 (since we provide the left coord over which the gene is painted)
  geneSize<-2
  
  geneDistanceFromCenter<-geneSize
  
  ##Gene LeftPos
  geneLeftPos<-(xAxisLim[1] + (xAxisLim[2] - xAxisLim[1])/2) - (geneSize*2) - geneSize ##Because we are getting the left most position, on the contrary more space in the right side
  
  ##Gene RightPos
  geneRightPos<-(xAxisLim[1] + (xAxisLim[2] - xAxisLim[1])/2) + (geneSize*2)
  
  
  ###Adding black line below Genes
  lines(x = c(tad_X_cord[1]-3,tad_X_cord[2]+3),
        y=c(tad_Y_cord[1],tad_Y_cord[1]),
        col="black", lwd=2, lty=1)
  
  
  for(genePos in c(geneLeftPos, geneRightPos)){
    
    ##Painting Genes
    gene_X_cord<-c(genePos, genePos+geneSize, genePos+geneSize, genePos)
    gene_Y_cord<-c(tad_Y_cord[1]-2,tad_Y_cord[1]-2,tad_Y_cord[1]+2,tad_Y_cord[1]+2)
    polygon(gene_X_cord, gene_Y_cord, col = "#0e3d61")
    
    ##Adding Gene Label, above gene
    ##Boxed labels: xpad,ypad arguments, used to expand the box beyond the text
    boxed.labels(x=genePos+geneSize/2,  y=tad_Y_cord[1]+4,labels = gene, cex=1.2, border = NA, bg ="white", 
                 xpad=1,
                 ypad=2, ##To allow the text to breathe
                 col="black"
    )
    
    
  }
  
  
  ##Adding Gene Duplicated Label, below genes, in the center. Tad coord middle
  
  ##Boxed labels: xpad,ypad arguments, used to expand the box beyond the text
  boxed.labels(x=xAxisLim[1] + (xAxisLim[2] - xAxisLim[1])/2,  y=tad_Y_cord[1]-4,labels = "Duplicated", cex=1.5, border = NA, bg ="white", 
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col="black"
  )
  
}





#########################################################################################################
## Paint Rearranged TAD on Duplications Where Enh Are Involved (hence some LongRange potential effects)
## AND WHEN THE GENE IS ALSO DUPLICATED
## For NEO-TAD cases. So interTAD SV
#########################################################################################################
paintGene_SV_Duplication_NEO_TAD<-function(nEnh_initial_left, nEnh_initial_right, gene, gene_breakp_line_type, situation,
                                           nEnh_kept_left, nEnh_kept_right,nEnh_gained,
                                           nEnh_other_domain,
                                           info_drawingGENE_TAD, info_drawingSecondaryTAD,
                                           tad_X_cord,
                                           tad_YCoord_Rearrangements,
                                           geneBreakP_Position_respectToTSS){
  
  ##Coded taking as reference paintGene_SV_Deletion_TAD()
  ##La unica diferencia es ... apuntar si hay alguna
  ########################################################
  ## Function to paint GENE SV re-arranged TAD 
  ## In the middle of the wt TADs
  ########################################################
  
  
  gene_sv_TAD_x_coord<-tad_X_cord
  
  ###########################################
  # Painting main scaffold TAD 
  
  ### Y axis position of the Inversion WT row
  tad_Y_cord<-tad_YCoord_Rearrangements
  
  
  ## Let's mantain color codes by TAD position
  ## If TAD is on the left blue, if it is on the right orange
  ## Regardless whether it is the gene or the secondary domain TAD
  ## In order not to confound the user if > 1 gene is affected
  
  tad_X_cord<-gene_sv_TAD_x_coord
  
  # 
  ##Here we can not use tad_X_coord because in all cases only one tad painted, hence tad_X_coord the same
  if(situation == "primaryTAD_Central"){
    colorTAD<-"#9999ff"
    color_otherFigure<-"#9999ff" ##color of the other TAD half
    
  }else if(situation == "primaryTAD_Sinistral"){
    colorTAD<-"#acdfeb"
    color_otherFigure<-"#e5edb2"  ##color of the other TAD half
    
  }else if(situation == "primaryTAD_Dextral"){
    colorTAD<-"#e5edb2"
    color_otherFigure<-"#acdfeb"
  }
  
  
  ###################
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
  
  ################################################################
  ## Defining additional relevant variables
  ## A lot of the following code from:  ##From paintGene_WT_TAD()
  ################################################################
  enhNumbersSize<-0.8
  
  distanceBreakpFromGene<-0.3 ##Distance breakpoint from the gene, when the breakpoint is between the gene & enh
  
  distanceEnhFromGene<-1.5 ##Distance between the enhancers and the gene
  
  distanceBreakpFromEnh<-0.3 ##Distance breakpoint from enhancers when the breakpoint located between enh and TAD border
  
  sizeEnhRegion<-1##Touching this will not really affect enh size, because created regardless of this
  
  distance_Yaxis_geneLabel<-2
  distance_Yaxis_EnhGene<-2
  distance_Yaxis_EnhLabel_EnhNumber<-1.3
  
  
  ###############################
  ##Obtaining breakpoint position
  ##If gene is relocated we paint the breakpos at the contrary tad breakpos
  ##If not, we paint its breakpos
  ##This comes from inv & transloc, for del or dup there will be no relocation
  #But I keep it for now
  #No va a haber relocation, pq el gen mantiene su posicion relativa en el neoTAD
  #Si esta a la derecha en su TAD esta a la derecha en el NEO tad, y lo mismo para la izquierda
  
  
  if(gene_breakp_line_type=="center"){
    ##This will never occur because gene disrupted so other plotting strategy
    ##But I will leave it for now
    
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
  
  ##Defining Breakpoint position
  breakp_x_coord<-c(breakPos,breakPos)##Tad central point for gene broken
  breakp_y_coord<-c(tad_Y_cord[1]-9,##-7 before
                    tad_Y_cord[3]+3)
  
  
  ###########################
  ## Painting the other color Figure
  ##   "otherColorFigure"
  ## The other color covering the TAD coming from the rearranged domain
  ############################
  
  ##coordinates in bottom-left bottom-right top-right top-left order
  
  ##########################################################
  ###First Getting The Coordinates of the "Other Figure"
  ##########################################################
  ##If situation == "primaryTAD_Central", so IntraTAD deletion or duplication(or completely outside), no color mixture, just the same TAD color painted...
  
  ##HERE THE LOGIC CHANGES WITH RESPECT TO DELETION paintGene_SV_Deletion_TAD()
  ##Which is the reference code
  
  if(situation=="primaryTAD_Dextral"){
    ## (as we're in a long range neo tad scenario)
    ##It means the gene is on the LEFT SIDE of the TAD so the other color goes from the TAD breakpos to the end of the TAD
    
    ##if breakpos before TAD center, we are going to color an irregular QUADRILATERAl
    ##So we will need 4 coordinates for the figure
    ##If its in the TAD center or after we will color a Triangle, so we will need only 3 coordinates
    if(breakPos>tad_X_cord[3]){##tad_XCoord_OnLeftSide[3],##TAD center position
      
      ##Right triangle to paint
      otherColorFigure_X_pos<-c(breakPos,
                                tad_X_cord[2],
                                breakPos)
      
      ###############################################
      ##First getting the triangle corner angle
      
      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs
      
      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm
      
      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians
      
      ##Now let's compute the height of interest. With respect to the right triangle, so right side
      catetoPlaying<-(tad_X_cord[2]-breakPos)
      heightOfInterest<-tan(cornerAngle)*catetoPlaying
      
      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest
      
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                mostWanted_Y_pos)
      
      
    }else if(breakPos==tad_X_cord[3]){
      
      ##Right triangle to paint
      otherColorFigure_X_pos<-c(breakPos,
                                tad_X_cord[2],
                                breakPos)
      
      ##In this case very easy the most wanted Y coordinate, is just the TAD peak
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                tad_Y_cord[3])
      
    }else if(breakPos<tad_X_cord[3]){
      ##So QUADRILATERAL to paint
      otherColorFigure_X_pos<-c(breakPos,
                                tad_X_cord[2],
                                tad_X_cord[3],
                                breakPos)
      
      ###############################################
      ##First getting the triangle corner angle
      
      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs
      
      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm
      
      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians
      
      ##Now let's compute the height of interest
      ##We compute it from the perspective of the right triangle!!! so Left side of the TAD
      catetoPlaying<-(breakPos-tad_X_cord[1])
      heightOfInterest<-tan(cornerAngle)*catetoPlaying
      
      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest
      
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                tad_Y_cord[3],
                                mostWanted_Y_pos
      )
      
    }
    
    
  }else if(situation=="primaryTAD_Sinistral"){
    ## (as we're in a long range neo tad scenario)
    ##It means the gene is on the RIGHT SIDE of the TAD so the other color goes from the TAD start to the breakpos
    
    ##if breakpos after TAD center, we are going to color an irregular QUADRILATERAl
    ##So we will need 4 coordinates for the figure
    ##If its in the TAD center or before we will color a Triangle, so we will need only 3 coordinates
    if(breakPos>tad_X_cord[3]){##tad_XCoord_OnLeftSide[3],##TAD center position
      
      ##So QUADRILATERAL to paint
      otherColorFigure_X_pos<-c(tad_X_cord[1],
                                breakPos,
                                breakPos,
                                tad_X_cord[3])
      
      ###############################################
      ##First getting the triangle corner angle
      
      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs
      
      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm
      
      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians
      
      ##Now let's compute the height of interest
      catetoPlaying<-(tad_X_cord[2]-breakPos)
      heightOfInterest<-tan(cornerAngle)*catetoPlaying
      
      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest
      
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                mostWanted_Y_pos,
                                tad_Y_cord[3])
      
    }else if(breakPos==tad_X_cord[3]){
      ##Right triangle to paint
      otherColorFigure_X_pos<-c(tad_X_cord[1],
                                breakPos,
                                breakPos)
      
      ##In this case very easy the most wanted Y coordinate, is just the TAD peak
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                tad_Y_cord[3])
      
      
    }else if(breakPos<tad_X_cord[3]){
      ##Right triangle to paint
      otherColorFigure_X_pos<-c(tad_X_cord[1],
                                breakPos,
                                breakPos)
      
      ###############################################
      ##First getting the triangle corner angle
      
      ##Getting Y coordinates, we need to apply Pitagoras for 1 Y position, the top-right
      ##https://www.youtube.com/watch?v=X7ykhX8q1Hs
      
      ##To get the corner angle we need to compute the arctangent
      #http://recursostic.educacion.es/descartes/web/materiales_didacticos/resolver_tri_rectangulos_pjge/Triangulos_rectangulos1.htm
      
      tadHeight<-tad_Y_cord[3]-tad_Y_cord[2]
      tad_halfLength<-tad_X_cord[2]-tad_X_cord[3]
      cornerAngle<-atan(tadHeight/tad_halfLength)##returned in radians
      
      ##Now let's compute the height of interest. With respect to the right triangle
      catetoPlaying<-(breakPos-tad_X_cord[1])
      heightOfInterest<-tan(cornerAngle)*catetoPlaying
      
      mostWanted_Y_pos<-tad_Y_cord[1]+heightOfInterest
      
      ##coordinates in bottom-left bottom-right top-right top-left order
      otherColorFigure_Y_pos<-c(tad_Y_cord[1],
                                tad_Y_cord[1],
                                mostWanted_Y_pos)
    }
    
  }
  
  ##########################################################
  ### Painting the "Other Figure"
  ##Recall only paint if situation different than primaryTadCentral
  ##Because if situation primaryTadCentral, it means intraTAD sv or completely deleted tad
  ##So we do not mix colors
  if(situation %in% c("primaryTAD_Sinistral","primaryTAD_Dextral")){
    polygon(otherColorFigure_X_pos, otherColorFigure_Y_pos, col = color_otherFigure, border = "#ffffff")
  }
  
  #######################
  ##Painting Breakpoint
  breakpointColor<-"brown3"
  lines(x=breakp_x_coord,
        y=breakp_y_coord,
        col=breakpointColor,
        lwd=2,
        lty=3)
  
  #Breakpoint label
  ##Boxed labels: xpad,ypad arguments, used to expand the box beyond the text
  boxed.labels(x=breakp_x_coord[1],  y=breakp_y_coord[1] - 1,labels = "Breakpoint", cex=0.8, border = NA, bg ="white", 
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col=breakpointColor
  )
  
  #TAD label
  # text(x=tad_X_cord[3], y=tad_Y_cord[1]+10, label="TAD", cex = 1.5)
  boxed.labels(x=tad_X_cord[3],  y=tad_Y_cord[1]+10,labels = "TAD", cex=1.5,
               border = NA, bg= NA, xpad=1, ypad = 2 )
  
  ###Adding black line below TAD
  lines(x = c(tad_X_cord[1]-3+1,tad_X_cord[2]+3-1),##Shrinking a bit the black line to avoid overlap with neighbours
        y=c(tad_Y_cord[1],tad_Y_cord[1]), 
        col="black", lwd=2, lty=1)
  
  #########################################
  ##Defining some profitable variables
  #########################################
  
  distanceEnhFromGene<-1.5 ##Distance between the enhancers and the gene
  
  sizeEnhRegion<-1##Touching this will not really affect enh size, because created regardless of this
  
  ###############################
  ## Painting gene body 
  ###############################
  
  ##genePos computed before in the function
  
  
  geneCenter<-genePos + geneSize/2 ##Used to direct the arrows to here
  gene_X_cord<-c(genePos,genePos+1,genePos+1,genePos)
  gene_Y_cord<-c(tad_Y_cord[1]-1,tad_Y_cord[1]-1,tad_Y_cord[1]+1,tad_Y_cord[1]+1)
  
  ###################
  #Adding gene body
  polygon(gene_X_cord, gene_Y_cord, col = "#0e3d61")
  
  #gene Label
  ##I want it to overlay the breakpoint line
  # https://stackoverflow.com/questions/45366243/text-labels-with-background-colour-in-r
  boxed.labels(x=geneCenter,  y=tad_Y_cord[1]-distance_Yaxis_geneLabel,
               labels = gene, cex=0.8, border = NA, bg ="white", 
               xpad=1,
               ypad=2##To allow the text to breathe
  )
  
  
  #####################
  ##Adding enhancers
  #####################
  enh_y_positions<-rep.int(x=tad_Y_cord[1], times = 4)##c(0,0,0,0)
  
  ##################################################
  ## First, PAINTING COGNATE ENHANCERS
  ## changes on enh originally surrounding the gene
  ##################################################
  
  #If nkept, different than nInitial we will paint only two enh, to represent
  #a reduction on its number
  ####Seguir por aqui morning, por linea 2238
  
  ##Gene in deletion not relocated because being in the relocation area
  ##would imply deletion, and hence would be treated as removed... 
  ##That's why I have removed an upcoming if
  
  # if(info_drawingSecondaryTAD$gene_relocated=="FALSE"){
  
  #If nkept, different than nInitial, as the breakpoint in the middle, we will paint only two enh, to represent
  #a reduction on its number.
  
  ####################################
  ####To the Right side gene Position
  ####################################
  
  if((geneBreakP_Position_respectToTSS=="beforeTSS")||(gene_breakp_line_type=="afterTSS_duplic_all")){
    ##Enhancers to  the right left intact, same as they were
    
    if(nEnh_initial_right>0){
      ##So there is at least 1 enh
      enhXpos<-genePos+distanceEnhFromGene##left margin of the enh cluster
      enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
      
      draw.ellipse(x=enh_x_positions,
                   y=enh_y_positions,a=0.1,col="#b2d235")
      
      #enhancers Label
      #text(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel              -distance_Yaxis_EnhGene, label="enhancers", cex = 0.8)
      boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      #Nr of enhancers
      #text(x=enhXpos + 0.4, y=enh_y_positions[1]-distance_Yaxis_geneLabel-distance_Yaxis_EnhGene-distance_Yaxis_EnhLabel_EnhNumber, label=paste0("n=",nEnh_initial_right), cex=enhNumbersSize)
      boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",nEnh_initial_right), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
            lwd=1.5)
      
      curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                     ,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T")
    }
    
  }else if((geneBreakP_Position_respectToTSS=="afterTSS") && (nEnh_initial_right>0)){
    
    
    if(gene_breakp_line_type=="afterTSS_duplic_none"){
      ##DOING NOTHING, we don't want to paint anything as no afterTSS enh duplicated, hence not bring to neoTAD in this scenario
    }else if(gene_breakp_line_type=="afterTSS_duplic_some"){
      
      ###########################
      ## PAINTING 2 ENHANCERS
      
      ##To quickly give the impression that some enhancers are lost, I'm going to paint 2enhancers
      ##That will be the two right enhn if the gene is on the right part of the plot, or the two left
      ##So there is at least 1 enh
      enhXpos<-genePos+distanceEnhFromGene##left margin of the enh cluster
      
      ##Painting Only 2 enhancers
      enh_x_positions<-c(enhXpos,enhXpos+0.3)
      
      draw.ellipse(x=enh_x_positions,
                   y=enh_y_positions[1:2],a=0.1,col="#b2d235")
      
      #enhancers Label
      boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      #Nr of enhancers
      ##Here kept will be bigger than intitial, since some of the initial duplicated, so nKept-nInitial
      boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",(nEnh_kept_right - nEnh_initial_right)), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
            lwd=1.5)
      
      curvedarrow(from = c(enhXpos+0.15,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                      ,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T" )
      
      
      ###END OF PAINTING 2 ENHANCERS
      ###################################
      
    }
    
  }
  
  
  
  ###################################
  ####To the Left side gene Position
  ###################################
  
  if((geneBreakP_Position_respectToTSS=="afterTSS")||(gene_breakp_line_type=="beforeTSS_duplic_all")){
    ##Enhancers to  the LEFT the same
    
    if(nEnh_initial_left>0){
      enhXpos<-genePos-distanceEnhFromGene##left margin of the enh cluster
      enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
      
      draw.ellipse(x=enh_x_positions,
                   y=enh_y_positions,a=0.1,col="#b2d235")
      
      #enhancers Label
      boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      
      #Nr of enhancers
      boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",nEnh_initial_left), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1]+2,enh_y_positions[1]+2),
            lwd=1.5)
      
      curvedarrow(from = c(enhXpos+0.5,enh_y_positions[1]+2), to = c(geneCenter#-0.25
                                                                     ,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T", curve = -1 )
    }
    
  }else if((geneBreakP_Position_respectToTSS=="beforeTSS") && (nEnh_initial_left>0)){
    
    
    if(gene_breakp_line_type=="beforeTSS_duplic_none"){
      ##DOING NOTHING
    }else if(gene_breakp_line_type=="beforeTSS_duplic_some"){
      
      ###########################
      ## PAINTING 2 ENHANCERS
      
      ##To quickly give the impression that some enhancers are lost, I'm going to paint just 2 enhancers
      ##So there is at least 1 enh
      enhXpos<-genePos-distanceEnhFromGene##left margin of the enh cluster
      
      ##Painting Only 2 enhancers
      #enh_x_positions<-c(enhXpos,enhXpos+0.3)
      enh_x_positions<-c(enhXpos+0.6, enhXpos+0.9)##The two second enhancers
      enhXpos<-enh_x_positions[1]
      
      draw.ellipse(x=enh_x_positions,
                   y=enh_y_positions[1:2],a=0.1,col="#b2d235")
      
      #enhancers Label
      boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene,
                   labels = "enhancers", cex=0.8, 
                   border = NA, bg ="white", 
                   xpad=1,
                   ypad=1 #To allow the text to breath
      )
      #Nr of enhancers
      ##Here kept will be bigger than intitial, since some of the initial duplicated, so nKept-nInitial
      boxed.labels(x=enhXpos + 0.15, y=enh_y_positions[1]
                   -distance_Yaxis_geneLabel
                   -distance_Yaxis_EnhGene
                   -distance_Yaxis_EnhLabel_EnhNumber,
                   labels = paste0("n=",(nEnh_kept_left - nEnh_initial_left)), cex=enhNumbersSize, 
                   border = NA, bg ="white", 
                   xpad=1.2,
                   ypad=1 #To allow the text to breath
      )
      
      ##enhancers horizontal line over them
      lines(x = c(enh_x_positions[1]-0.1,
                  enh_x_positions[length(enh_x_positions)]+0.1),
            y = c(enh_y_positions[1] + 2,enh_y_positions[1] + 2),
            lwd=1.5)
      
      enhHorizontalLineCoord<-c(enh_x_positions[1]-0.1,
                                enh_x_positions[length(enh_x_positions)]+0.1)
      
      enhHorizLineCenter<-enhHorizontalLineCoord[1]+(enhHorizontalLineCoord[2]-enhHorizontalLineCoord[1])/2
      
      
      curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(geneCenter#+0.25
                                                                            ,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T", curve=-1 )
      
      
      ###END OF PAINTING 2 ENHANCERS
      ###################################
      
    }
    
  }
  
  
  ##################################################
  ## TIME FOR PAINTING ECTOPIC  ENHANCERS
  ## changes on enh originally surrounding the gene
  ##################################################
  ##For intraTAD SV nEnh_gained will be 0 so no painting
  ##If there is no gain, no paint
  ##If nEnh gained == nEnhInitialOther Domain, the four enh are painted
  ##On the contrary only two enh are painted
  
  ##If the gene is not relocated we paint the enh on the middle of its segment
  
  
  if(nEnh_gained>0){
    ##So some enh need to be painted
    
    ##When the gene not relocated, to locate the enh we work with the "other figure" coordinates
    ##Because the ectopic enhancers are the elements relocated now
    ##We paint the enh in the middle of the other Figure coord
    x_space<-c(otherColorFigure_X_pos[1],otherColorFigure_X_pos[2])
    
    x_space_start<-x_space[1]
    x_space_end<-x_space[2]
    
    enhXpos<-(x_space_start + (x_space_end-x_space_start)/2)-0.5
    enhXpos<-enhXpos-0.45
    
    if(nEnh_gained==nEnh_other_domain){
      ##All enh from the other Domain Gained so
      #we paint 4 enh
      
      enh_x_positions<-c(enhXpos,enhXpos+0.3,enhXpos+0.6, enhXpos+0.9)
      
    }else if(nEnh_gained < nEnh_other_domain){
      ##So we are going to paint only 2 enhancers
      
      ##So we paint the two enh on the left
      enh_x_positions<-c(enhXpos,enhXpos+0.3)
    }
    
    ####################################################
    ## Pintamos los ECTOPIC enhancers
    ##Lo de arriba solo para calcular posicion de enh
    ##get enh positions center
    ##para pintar hacia derecha o izquierda mirar si el geneCenter queda a la izq o a la derecha
    
    ##
    draw.ellipse(x=enh_x_positions,
                 y=enh_y_positions,a=0.1,col="#b2d235")
    
    #enhancers Label
    #Below the level of the other enhancers, to avoid overlaps
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                 -distance_Yaxis_geneLabel
                 -distance_Yaxis_EnhGene
                 -distance_Yaxis_EnhLabel_EnhNumber
                 -distance_Yaxis_EnhLabel_EnhNumber,
                 labels = "enhancers", cex=0.8, 
                 border = NA, bg ="white", 
                 xpad=1,
                 ypad=1 #To allow the text to breath
    )
    #Nr of enhancers
    boxed.labels(x=enhXpos + 0.4, y=enh_y_positions[1]
                 -distance_Yaxis_geneLabel
                 -distance_Yaxis_EnhGene
                 -distance_Yaxis_EnhLabel_EnhNumber
                 -distance_Yaxis_EnhLabel_EnhNumber
                 -distance_Yaxis_EnhLabel_EnhNumber, labels = paste0("n=",nEnh_gained), cex=enhNumbersSize, 
                 border = NA, bg ="white", 
                 xpad=1.2,
                 ypad=1 #To allow the text to breath
    )
    
    ##enhancers horizontal line over them
    lines(x = c(enh_x_positions[1]-0.1,
                enh_x_positions[length(enh_x_positions)]+0.1),
          y = c(enh_y_positions[1]+2,enh_y_positions[1]+2),
          lwd=1.5)
    
    
    enhHorizontalLineCoord<-c(enh_x_positions[1]-0.1,
                              enh_x_positions[length(enh_x_positions)]+0.1)
    
    enhHorizLineCenter<-enhHorizontalLineCoord[1]+(enhHorizontalLineCoord[2]-enhHorizontalLineCoord[1])/2
    
    ##Painting curve connecting enh with gene
    if(enhXpos > geneCenter){
      ##So, the secondary TAD is on the right side
      curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(geneCenter,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T", curve = 1)
    }else{
      ##So the secondary TAD is on the left side, adjust row curvature, to mantain it over the X axis
      curvedarrow(from = c(enhHorizLineCenter,enh_y_positions[1]+2), to = c(geneCenter,enh_y_positions[1]+2),
                  arr.pos = 0, arr.type="T", curve = -1)
    }
    
  }
  
  ##################################################################################
  ##Adding Triangle Shape over gene, if it has any enh, to represent the arrow end
  ##################################################################################
  if((nEnh_kept_left>0)||(nEnh_kept_right>0)||(nEnh_gained>0)){
    triangle_X_Coord<-c(geneCenter-0.3, geneCenter, geneCenter+0.3)
    heightOrizontalLineEnh<-enh_y_positions[1]+2
    triangle_Y_Coord<-c(heightOrizontalLineEnh+0.3, heightOrizontalLineEnh-0.2, heightOrizontalLineEnh+0.3)
    polygon(triangle_X_Coord, triangle_Y_Coord, col = "black", border = "black")
  }
  
  
  
  
  #############################
  ## return profitable info ###
  #############################
  ##break pos, to easily select the other one
  info_drawing<-list("breakPos"=breakPos,
                     "geneCenter"=geneCenter)
  
  return(info_drawing)
  
}



