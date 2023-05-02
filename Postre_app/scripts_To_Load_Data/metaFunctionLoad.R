#########################################################
## Script to load all functions used by the tool at once
#########################################################
source(file = "functions/MasterWrapperSinglePrediction.R",local=TRUE)##many other functions loaded here
source(file = "functions/error_Report.R",local=TRUE)
source(file = "functions/error_Report_MultipleSV_submission.R",local=TRUE)
source(file = "functions/noGenesFound_report.R",local=TRUE)
source(file = "functions/mergingMultiPhenoPredictions.R",local=TRUE)

source(file = "functions/deleteImages_Previous_Analyses.R",local=TRUE)
source(file = "functions/processing_GenomicCoordinates.R",local=TRUE)
source(file = "functions/getBetweenTAD.R",local=TRUE)

##############################################
##Functions for Multiple Patients Submission
##############################################
source("functions/multiple_SV_Functions/cohortResults_Parser.R",local=TRUE)
source("functions/multiple_SV_Functions/multipleStats_ExplorePreviousPat_htmlGeneration.R",local=TRUE)

##To deal with rounding .5 problems
##round2 function
source(file = "functions/roundingHalfs.R",local=TRUE) ##Local == FALSE to be loaded in the global env so that all functions can find it

#################################################################################
## Functions from masterScoring function to not resource upon multiple iterations
#################################################################################
source("functions/GenomicData_Loader.R",local=TRUE)
source("functions/PatientPrediction_basedOnSingleTADmap.R",local=TRUE)
##Try to put it here to avoid massive reloading
source("scripts_To_Load_Data/cargarFunciones.R",local=TRUE)
source("functions/RankingGenes_Deciphering_etiology.R",local=TRUE)
source("functions/Integrating_TAD_Predictions.R",local=TRUE)
source("functions/summary_positionalInfoGenes.R",local=TRUE)
source("functions/Master_Summary_Matrix_IntegratesMainPhaseResults.R",local=TRUE)
source("functions/CheckSameDomain.R",local=TRUE)
source("functions/GetUncertainty.R",local=TRUE)

source(file = "functions/CalculatingNrEnhBeforeAndAfter.R",local=TRUE)

##Also loading some thresholds, for expression
source("scripts_To_Load_Data/ExpressionThresholds_ForGofAndLof.R",local=TRUE)

