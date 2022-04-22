##########################################################
## Plotting Regulatory Impact Changes for Direct Effects
########################################################
wt_section_RegulatoryImpact<-function(targetGene){
  #############################
  ## Plotting for Allele 1
  #############################
  
  ######################
  ##Plot WT line
  y_wt_start<-16-8
  y_wt_width<-1
  y_wt_allCoords<-c(y_wt_start,y_wt_start,y_wt_start+y_wt_width,y_wt_start+y_wt_width)
  #y_wt_height<-c(16,16,19,19)-5
  
  text(x=10, y=y_wt_allCoords[4]+2, label="Allele 1", cex = 1.3, font = 2)
  
  # polygon(c(5,15,15,5),y_wt_allCoords , col = "#ff9966", border = "black")###6699ff"
  y_dnaPos<-y_wt_start+y_wt_width/2
  lines(c(5,15),c(y_dnaPos,y_dnaPos), col = "black", lwd=2, lty=1)
  
  polygon(c(8.5,11.5,11.5,8.5), y_wt_allCoords, col = "#0e3d61", border = "black") 
  boxed.labels(x=-1,  y=(y_wt_allCoords[1] +(y_wt_width/2)) ,labels = "Control-DNA", cex=1, border = NA, bg ="white",
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col="#000000",
               font=2
  )
  
  boxed.labels(x=10,  y=y_wt_allCoords[3]+0.4 ,labels = targetGene, cex=1, border = NA, bg ="transparent",
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col="#000000",
               font=2
  )
  
  
  ####################################
  ## Plotting for Allele 2
  ####################################
  ######################
  ##Plot WT line
  y_wt_start<-16-8
  y_wt_width<-1
  y_wt_allCoords<-c(y_wt_start,y_wt_start,y_wt_start+y_wt_width,y_wt_start+y_wt_width)
  #y_wt_height<-c(16,16,19,19)-5
  
  xPositions_Allele2_DNA<-c(5,15,15,5) + 20
  x_Positions_Allele2_Gene<-c(8.5,11.5,11.5,8.5) + 20
  
  text(x=xPositions_Allele2_DNA[1] + (xPositions_Allele2_DNA[2] - xPositions_Allele2_DNA[1])/2, y=y_wt_allCoords[4]+2, label="Allele 2", cex = 1.3, font = 2)
  
  # polygon(xPositions_Allele2_DNA, y_wt_allCoords , col = "#ff9966", border = "black")###6699ff"
  y_dnaPos<-y_wt_start+y_wt_width/2
  lines(c(xPositions_Allele2_DNA[1],xPositions_Allele2_DNA[2]),c(y_dnaPos,y_dnaPos), col = "black", lwd=2, lty=1)
  
  polygon(x_Positions_Allele2_Gene, y_wt_allCoords, col = "#0e3d61", border = "black") 
  
  boxed.labels(x=xPositions_Allele2_DNA[1] + (xPositions_Allele2_DNA[2] - xPositions_Allele2_DNA[1])/2,  y=y_wt_allCoords[3]+0.4 ,labels = targetGene, cex=1, border = NA, bg ="transparent",
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col="#000000",
               font=2
  )
  
  
  
}


deletion_regulatoryImpact<-function(targetGene){

  #############################
  ## Plotting for Allele 1
  #############################
  
  ######################
  ##Plot Deleted Line
  y_sv_start<-6
  y_sv_width<-1
  y_sv_allCoords<-c(y_sv_start,y_sv_start,y_sv_start+y_sv_width,y_sv_start+y_sv_width)
  
  ##WT dna was presenting 10 units length. The gene 3. Hence upon gene removal 7 dna unit left
  ##So I need to increase 1.5 units the Start and decrease 1.5 units the end of the deleted dna
  
  # polygon(c(6.5,13.5,13.5,6.5), y_sv_allCoords, col = "#ff9966", border = "black")
  y_dnaPos<-y_sv_start+y_sv_width/2
  lines(c(6.5,13.5),c(y_dnaPos,y_dnaPos), col = "black", lwd=2, lty=1)
  
  boxed.labels(x=-1,  y=(y_sv_allCoords[1] +(y_sv_width/2)) ,labels = "Patient-DNA", cex=1, border = NA, bg ="white",
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col="#000000",
               font=2
  )
  
  ########################
  ## Plotting For Allele 2
  #########################
  
  ######################
  ##Plot Deleted Line -- In this case Unaltered for Patient. We are assuming Heterozygous SV
  y_sv_start<-6
  y_sv_width<-1
  y_sv_allCoords<-c(y_sv_start,y_sv_start,y_sv_start+y_sv_width,y_sv_start+y_sv_width)
  
  xPositions_Allele2_DNA<-c(5,15,15,5) + 20
  x_Positions_Allele2_Gene<-c(8.5,11.5,11.5,8.5) + 20
  
  # polygon(xPositions_Allele2_DNA,y_sv_allCoords , col = "#ff9966", border = "black")###6699ff"
  y_dnaPos<-y_sv_start+y_sv_width/2
  lines(c(xPositions_Allele2_DNA[1],xPositions_Allele2_DNA[2]),c(y_dnaPos,y_dnaPos), col = "black", lwd=2, lty=1)
  
  polygon(x_Positions_Allele2_Gene, y_sv_allCoords, col = "#0e3d61", border = "black") 
  
}

truncation_regulatoryImpact<-function(targetGene){
  
  #############################
  ## Plotting for Allele 1
  #############################
  
  ######################
  ##Plot Truncated Gene Line
  y_sv_start<-6
  y_sv_width<-1
  y_sv_allCoords<-c(y_sv_start,y_sv_start,y_sv_start+y_sv_width,y_sv_start+y_sv_width)
  
  ##Let's paint a line in the middle of the Gene Disrupting It
  
  # polygon(c(5,15,15,5), y_sv_allCoords, col = "#ff9966", border = "black")
  y_dnaPos<-y_sv_start+y_sv_width/2
  lines(c(5,15),c(y_dnaPos,y_dnaPos), col = "black", lwd=2, lty=1)
  
  polygon(c(8.5,11.5,11.5,8.5), y_sv_allCoords, col = "#0e3d61", border = "black") 
  
  boxed.labels(x=-1,  y=(y_sv_allCoords[1] +(y_sv_width/2)) ,labels = "Patient-DNA", cex=1, border = NA, bg ="white",
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col="#000000",
               font=2
  )
  
  ##Adding Truncation Line, in the middle of the gene
  ##First White Line to create the gap, and then overlay with bigger red line
  x_posBreak<-10
  lines(x=c(x_posBreak,x_posBreak), y = c(y_sv_start, y_sv_start + y_sv_width), 
        col="#FFFFFF",
        lwd=10)
  
  lines(x=c(x_posBreak,x_posBreak), y = c(y_sv_start-0.3, y_sv_start + y_sv_width+0.3), 
        col="brown3",
        lwd=2,
        lty=3)
  
  boxed.labels(x=10,  y=y_sv_allCoords[3]-y_sv_width-0.7 ,labels = "Breakpoint", cex=1, border = NA, bg ="transparent",
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col="brown3",
               font=2
  )
  
  
  
  
  ########################
  ## Plotting For Allele 2
  #########################
  
  ######################
  ##Plot Deleted Line -- In this case Unaltered for Patient. We are assuming Heterozygous SV
  y_sv_start<-6
  y_sv_width<-1
  y_sv_allCoords<-c(y_sv_start,y_sv_start,y_sv_start+y_sv_width,y_sv_start+y_sv_width)
  
  xPositions_Allele2_DNA<-c(5,15,15,5) + 20
  x_Positions_Allele2_Gene<-c(8.5,11.5,11.5,8.5) + 20
  
  # polygon(xPositions_Allele2_DNA,y_sv_allCoords , col = "#ff9966", border = "black")###6699ff"
  y_dnaPos<-y_sv_start+y_sv_width/2
  lines(c(xPositions_Allele2_DNA[1],xPositions_Allele2_DNA[2]),c(y_dnaPos,y_dnaPos), col = "black", lwd=2, lty=1)
  
  polygon(x_Positions_Allele2_Gene, y_sv_allCoords, col = "#0e3d61", border = "black") 
  
}


duplication_regulatoryImpact<-function(targetGene){
  
  #############################
  ## Plotting for Allele 1
  #############################
  
  ######################
  ##Plot Truncated Gene Line
  y_sv_start<-6
  y_sv_width<-1
  y_sv_allCoords<-c(y_sv_start,y_sv_start,y_sv_start+y_sv_width,y_sv_start+y_sv_width)
  
  ##Let's paint a larger the DNA and we will add another gene sligthlty to the right
  
  # polygon(c(5,20,20,5), y_sv_allCoords, col = "#ff9966", border = "black")
  y_dnaPos<-y_sv_start+y_sv_width/2
  lines(c(5,20),c(y_dnaPos,y_dnaPos), col = "black", lwd=2, lty=1)
  
  polygon(c(8.5,11.5,11.5,8.5), y_sv_allCoords, col = "#0e3d61", border = "black") 
  
  ##Extra gene Copy
  polygon(c(8.5,11.5,11.5,8.5)+5, y_sv_allCoords, col = "#0e3d61", border = "black") 
  
  boxed.labels(x=-1,  y=(y_sv_allCoords[1] +(y_sv_width/2)) ,labels = "Patient-DNA", cex=1, border = NA, bg ="white",
               xpad=1,
               ypad=2, ##To allow the text to breathe
               col="#000000",
               font=2
  )
  
  ########################
  ## Plotting For Allele 2
  #########################
  
  ######################
  ##Plot Deleted Line -- In this case Unaltered for Patient. We are assuming Heterozygous SV
  y_sv_start<-6
  y_sv_width<-1
  y_sv_allCoords<-c(y_sv_start,y_sv_start,y_sv_start+y_sv_width,y_sv_start+y_sv_width)
  
  xPositions_Allele2_DNA<-c(5,15,15,5) + 20
  x_Positions_Allele2_Gene<-c(8.5,11.5,11.5,8.5) + 20
  
  # polygon(xPositions_Allele2_DNA,y_sv_allCoords , col = "#ff9966", border = "black")###6699ff"
  y_dnaPos<-y_sv_start+y_sv_width/2
  lines(c(xPositions_Allele2_DNA[1],xPositions_Allele2_DNA[2]),c(y_dnaPos,y_dnaPos), col = "black", lwd=2, lty=1)
  
  polygon(x_Positions_Allele2_Gene, y_sv_allCoords, col = "#0e3d61", border = "black") 
  
}




