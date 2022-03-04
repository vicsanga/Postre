##########################################################################
# Script to load Postre required libraries
# Libraries not currenlty installed on user computer will be downloaded
##########################################################################

if (!require("pacman")){install.packages("pacman")} 
pacman::p_load(shiny,waiter,shinybusy,shinyjs,shinyWidgets,DT,plotrix,shape,diagram)
