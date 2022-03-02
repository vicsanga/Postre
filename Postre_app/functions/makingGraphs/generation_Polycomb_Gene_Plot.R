##############################################
### Generating In Parallel Polycomb Gene Plot

####################################
##png with maximum resolution 300dpi
# png(filename = fullOutpPath, width = 12, height = 8, units = "in", res = 300 )
fullOutpPath<-paste0("/home/victor/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/SV_app/www/",
                     "GOF_Image.png")
png(filename = fullOutpPath, width = 800, height = 800, res = 300, pointsize = 4 )
##Jugar con dimensiones, y con tamanyo de labels que asi como esta no se van a ver en el html


##Adjust image, to exclude spaces outside canvas (drawing area)
par(mar = c(0,0,0,0))
canvas_X_limits<-c(-3,14)##max adjust to figure 
canvas_Y_limits<-c(2,14)##max adjust to content displayed


##OPTION REMOVING AXIS
#UPON COMPLETION USING THIS
plot(x=0:15, y=0:15, type = "n",
     ylim = canvas_Y_limits,
     xlim = canvas_X_limits,
     xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '') ##Upon completion, REMOVE AXIS

#Adding H3K27me3 Label

text(x=4, y=13, label="Broad H3K27me3 Domain", cex = 3, font = 2, col="#D64045")

##Adding Polycomb Mark cartoon
# https://stackoverflow.com/questions/31538534/plotting-half-circles-in-r
# library(grid)
# vp <- viewport(width=0.5, height=0.5, clip  = "on")
# grid.circle(0.5,0,r=0.5, gp = gpar(fill = 'red'), vp = vp)

# https://stackoverflow.com/questions/16504452/color-a-portion-of-the-normal-distribution

x=seq(-2,12,length=200)
y=dnorm(x,mean = 4, sd = 2)*40+3.9
lines(x,y,type="l", lwd=2, col="#D64045")

##Vamos a tratar de hacer el fill
# x=seq(-2,12,length=200)
# y=dnorm(x,mean = 4, sd = 2)*40+4
polygon(c(x),c(y),col="#D64045", border = NA )


##DNA, painting
polygon(x = c(-2,12,12,-2), y = c(2,2,4,4), col = "#ff9966", border = "black")

##Gene Painting
polygon(x = c(4,10,10,4), y = c(2,2,4,4), col = "#0000ff", border = "black")

##CpGi painting
polygon(x = c(-1.25,0.75,0.75,-1.25), y = c(2,2,4,4), col = "#E9FFF9", border = "black")

polygon(x = c(1.25,3.25,3.25,1.25), y = c(2,2,4,4), col = "#E9FFF9", border = "black")

##Adding Gene Label
#Gene
text(x=7, y=3, label="Gene", cex = 2.5, font = 2, col="#FFFFFF")

#Adding CpGi labels
text(x=-0.25, y=3, label="CpGi", cex = 2, font = 2, col="#000000")
text(x=2.25, y=3, label="CpGi", cex = 2, font = 2, col="#000000")

dev.off()##Saving Graph
