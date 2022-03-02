###################################################################
##Script para cargar todas las funciones que usamos en el software
###################################################################
##Check patients info
source("functions/CheckPatientsInfo.R",local = TRUE)

##Retrieve the affected genomic regions
source("functions/AffectedRegions.R",local = TRUE)

##Retrieve info regarding the affected genes
source("functions/AffectedGenes.R",local = TRUE)

##Retrieve info regarding enhancers
source("functions/EnhancersSituation.R",local = TRUE)

##Gene Matrixes evaluation all features
source("functions/EvaluatingGenesSituation_matrices.R",
       local = TRUE)

# ##Para complementar las matrices de evaluacion de la situacion de los genes pero con un poco de info mas subjetiva
# ##no tan automaticas las asociaciones
# source("functions/ComplementingInfo_evaluatingGenesSituation.R",local = TRUE)

###Function to rank genes
##Asses their involvment in the disease
source("functions/RankingGenes_Deciphering_etiology.R",local = TRUE)