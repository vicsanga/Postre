#############################
##Barplotting function
##############################
# library(ggplot2)
sv_barplotting<-function(barplotValues, tagsBarplotValues,title){
  
  ##################
  ##For now doing with base R, in future maybe ggplot2
  ##################

  ##Defining y max value barplot, to leave a bit of space with the header
  #let's add a 5% extra
  maxBarplotValue<-max(barplotValues) + max(barplotValues)*0.05
  
  barplot(height = barplotValues,
          ylim=c(0, maxBarplotValue),
          col=c("#ffad33","#bb99ff"),
          main=title,
          names.arg = tagsBarplotValues,
          horiz = FALSE,
          las=1,
          cex.axis = 2.5,
          cex.lab=2,
          cex.main=2.5,
          cex.sub=2,
          cex.names = 2.5)
  
  
  
  
  
  
  # dataBarplot<-data.frame(
  #   name=factor(tagsBarplotValues, 
  #               levels = tagsBarplotValues, 
  #               ordered=TRUE),
  #   value=barplotValues
  # )
  # 
  # ##Barplot
  # ##Seguir mirando tengo que cambiar colores, mirar como exportar en buena calidad
  # ggplot(data = dataBarplot, aes(x=name, y=value, fill=name)) +
  #   geom_bar(stat="identity") +
  #   ggtitle(title)
  
  
}