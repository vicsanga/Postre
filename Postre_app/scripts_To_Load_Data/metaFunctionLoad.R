#########################################################
## Script to load all functions used by the tool at once
#########################################################
source(file = "functions/MasterWrapperSinglePrediction.R")##many other functions loaded here
source(file = "functions/error_Report.R")
source(file = "functions/error_Report_MultipleSV_submission.R")
source(file = "functions/noGenesFound_report.R")
source(file = "functions/mergingMultiPhenoPredictions.R")

source(file = "functions/deleteImages_Previous_Analyses.R")
source(file = "functions/processing_GenomicCoordinates.R")
source(file = "functions/getBetweenTAD.R")

##############################################
##Functions for Multiple Patients Submission
##############################################
source("functions/multiple_SV_Functions/cohortResults_Parser.R")
source("functions/multiple_SV_Functions/multipleStats_ExplorePreviousPat_htmlGeneration.R")

##To deal with rounding .5 problems
##round2 function
source(file = "functions/roundingHalfs.R") ##Local == FALSE to be loaded in the global env so that all functions can find it

#################################################################################
## Functions from masterScoring function to not resource upon multiple iterations
#################################################################################
source("functions/GenomicData_Loader.R")
source("functions/PatientPrediction_basedOnSingleTADmap.R")
##Try to put it here to avoid massive reloading
source("scripts_To_Load_Data/cargarFunciones.R")
source("functions/RankingGenes_Deciphering_etiology.R")
source("functions/Integrating_TAD_Predictions.R")
source("functions/summary_positionalInfoGenes.R")
source("functions/Master_Summary_Matrix_IntegratesMainPhaseResults.R")
source("functions/CheckSameDomain.R")
source("functions/GetUncertainty.R")

source(file = "functions/CalculatingNrEnhBeforeAndAfter.R")

##Also loading some thresholds, for expression
source("scripts_To_Load_Data/ExpressionThresholds_ForGofAndLof.R")

