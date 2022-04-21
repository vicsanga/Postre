##################################
## Haploinsufficiency score plot
##################################
haploinsufficiencyScorePlot<-function(targetHI_score, targetGene){
  
  ###Adjust expression value to generate the barplot height
  maxBarHeight<-1
  if(targetHI_score>1){
    barplotHeight<-maxBarHeight ##Max bar height
  }else if(targetHI_score<0.05){
    barplotHeight<-0.05 ##Min bar height, to paint a bit, if the HI score is tiny
  }else{
    barplotHeight<-targetHI_score
  }
  
  maxBarplotValue<-maxBarHeight + 0.5 ##space for additional tags
  
  ##Puedo pintar yo el cuadro y arreglado
  #When I want to se the axis, uncomment the 3 following lines
  
  canvas_X_limits<-c(-4,14)
  canvas_Y_limits<-c(0,maxBarplotValue)
  # plot(x=1, y=1, type = "n",
  #      ylim = canvas_Y_limits,
  #      xlim = canvas_X_limits)
  
  ##OPTION REMOVING AXIS
  #UPON COMPLETION USING THIS
  plot(x=0:1, y=0:1, type = "n",
       ylim = canvas_Y_limits,
       xlim = canvas_X_limits,
       xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '') ##Upon completion, REMOVE AXIS
  
  #Adding plot header
  #5.8 to ensure that a gene with 9 letters perfectly fit
  text(x=5.8, y=1.45, label=paste0(targetGene," Dosage Sensitivity"), cex = 3.5, font = 2)
  
  #Adding score
  text(x=6, y=1.3, label=paste0("DS score = ",
                                     round2(x = targetHI_score,
                                     digits = 2)), cex = 3,
       font=2)
  
  
  
  ###Painting Bar
  polygon(x = c(3,9,9,3), y=c(canvas_Y_limits[1],canvas_Y_limits[1],barplotHeight, barplotHeight),
          col="#477599",#"#9ED8DB",
          border = "#ffffff")
  
  
  ##Add horizontal bars 
  lines(x = c(2.5,9.5),y = c(1,1), lty = 3, lwd= 1)

  ##Adding tags
  text(x=9.5, y=1, label="score = 1", cex = 3, pos = 4, font = 2)
  
  ##Poner en negrita la que corresponda
  
  ##Adding tags
  text(x=-3.7, y=0.9, label=paste0("DS gene"), cex = 3, pos = 4, font = 2)
  
  ##Adding tags
  text(x=-3.7, y=0.05, label=paste0("Not DS gene"), cex = 3, pos = 4, font = 2)  
  
  ##Adding "axes" inssided of canvas
  #Y axes
  #lines(x=c(-2.7,-2.7),y=c(canvas_Y_limits[1],canvas_Y_limits[2]-0.3),lwd=3)
  #X axes
  lines(x=c(-3.7,canvas_X_limits[2]),y=c(canvas_Y_limits[1],canvas_Y_limits[1]), lwd=1)
  
  ##Adding tags
  #text(x=-3.5, y=0.7, label="Haploinsufficiency Score", cex = 3, pos = 1, font = 2,srt=90)
  
}