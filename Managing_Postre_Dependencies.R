##########################################################################
# Script to load Postre required libraries
# Libraries not currenlty installed on user computer will be downloaded
##########################################################################
if (!require("pacman")){install.packages("pacman")} 
pacman::p_load(shiny,waiter,shinybusy,shinyjs,shinyWidgets,DT,plotrix,shape,diagram)

# ##Probar un parseo y ejecucion posterior
# install.packages(c("dplyr", "stringr","esquisse"))
# library(esquisse)
###Hacerlo con el more efficient y solo un script al que se haga source, luego si necesito mas a la gente le saltara en ese script
##Y luego meter el library(shiny)
##Pq igual luego si quiero meter algo de bioconductor me jodera pq la funcion no servira
# https://statsandr.com/blog/an-efficient-way-to-install-and-load-r-packages/#librarian-package