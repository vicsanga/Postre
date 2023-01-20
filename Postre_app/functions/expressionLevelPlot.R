##############################
##Function expressionLevelPlot
##############################

expressionLevelPlot<-function(targetExpression, targetGene, targetPhase, maxExpression, minExpression){
  
  ###Adjust expression value to generate the barplot height
  maxBarHeight<-12
  if(targetExpression>maxExpression){
    barplotHeight<-maxBarHeight ##Max bar height
  }else if(targetExpression<minExpression){
    barplotHeight<-0.7 ##Min bar height
  }else if((targetExpression>=minExpression) && (targetExpression<maxExpression)){
    
    ##Normalizar altura target expression proporcional a la diferencia
    ##Intermediate space goes from 1 to 10 in the Y axis (maxValue is the maxExpression, min value the minExpression), 
    
    proportionalValue<-(targetExpression - minExpression)/(maxExpression-minExpression)
  
   ##1+ because starting point in the Y axis for this category is 1 (1.8 on latest version)
   ##9* because the max point is 10 in the Y axis, so the difference between start and end are 9 points
   #On latest version 9.2 -1.8 = 7.4 
   barplotHeight<-1.8 + 7.4*proportionalValue
  }
  
  maxBarplotValue<-maxBarHeight + 3 ##space for additional tags
  
  ##Puedo pintar yo el cuadro y arreglado
  #When I want to se the axis, uncomment the 3 following lines
  
  canvas_X_limits<-c(-4,15)
  canvas_Y_limits<-c(-1,16.5)
  # plot(x=1, y=1, type = "n",
  #      ylim = canvas_Y_limits,
  #      xlim = canvas_X_limits)
  
  ##OPTION REMOVING AXIS
  #UPON COMPLETION USING THIS
  plot(x=0:15, y=0:15, type = "n",
       ylim = canvas_Y_limits,
       xlim = canvas_X_limits,
       xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '') ##Upon completion, REMOVE AXIS
  
  #Adding plot header

  text(x=6, y=16, label=paste0(targetGene," Expression"), cex = 3.5, font = 2)
  
  #Adding fpkms
  text(x=6, y=14, label=paste0(round2(x = targetExpression,
                                                      digits = 1), " FPKM"), cex = 3,
       font=2)
  
  
  
  ###Painting Bar
  polygon(x = c(3,9,9,3), y=c(canvas_Y_limits[1],canvas_Y_limits[1],barplotHeight, barplotHeight),
          col="#1eabc7",#"#9ED8DB",
          border = "#ffffff")
  
  
  ##Add horizontal bars Expression Thresholds
  ##1fpkm
  lines(x = c(2.5,9.5),y = c(1.8,1.8), lty = 3, lwd=1)

  #10fpkm
  lines(x = c(2.5,9.5),y = c(9.2,9.2), lty = 3, lwd=1)
  
  ##Poner leyenda de que indican las lineas
  ##Poner en negrita la que corresponda
  
  ##Adding tags
  text(x=9.5, y=9.2, label=paste(">",maxExpression," FPKM", sep="", collapse=""), cex = 2.5, pos = 4, font = 2)
  
  ##Adding tags
  text(x=9.5, y=1.8, label=paste(">",minExpression," FPKM", sep="", collapse=""), cex = 2.5, pos = 4, font = 2)
  
  
  
  ##Poner en negrita la que corresponda
  
  ##Adding tags
  text(x=-3.2, y=11, label=paste0("Expressed"), cex = 3, pos = 4, font = 2)
  
  ##Adding tags
  text(x=-3.2, y=5, label=paste0("    Lowly\nExpressed"), cex = 3, pos = 4, font = 2)

  ##Adding tags
  text(x=-3.2, y=0, label=paste0("      Not\nExpressed"), cex = 3, pos = 4, font = 2)  
  
   ##Adding "axes" inssided of canvas
  #Y axes
  #lines(x=c(-2.7,-2.7),y=c(canvas_Y_limits[1],canvas_Y_limits[2]-2),lwd=3)
  #X axes
  lines(x=c(-3.7,canvas_X_limits[2]),y=c(canvas_Y_limits[1],canvas_Y_limits[1]), lwd=1)
  
  ##Adding tags
  #text(x=-3.5, y=8, label="Expression Level", cex = 3, pos = 1, font = 2,srt=90)
  
}
