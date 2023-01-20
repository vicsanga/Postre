##########################################################################################################################
##Creating MultiDataList object, will be placed top level of the app, and will be passed to the prediction functions
##Thus It will only be loaded ONCE even if multiple patients are analysed
##########################################################################################################################
# setwd("~/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/Postre/Postre_app")
#####################################################################
## Script to  host all the objects that the software needs to load
#####################################################################

##No se si esto seria mas efectivo cargarlo en las funciones y luego que se elimine
##o ir eliminandolo al haber ejecutado las funciones
##ver temas de memoria

########################
##Genes Annotation
##Obtenido del script:
# /home/victor/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/ProcesadoAnotaciones/ProcesadoAnotacionGenecode/Hsapiens/1_For_SV_App_ProcesadoAnotacionGenecode18.R
load("data/genesAnnotation.RData")##If we modify the content of the variable is alright. But mantain the path name

#########################################################
## Loading TauExp values: es decir, que si en algun tejido fetal o NCC expression considerable su TAU
## de lo contrario un 0 de TAU. (es decir si su expression es muy baja) en todos los estadios
## ajustar para alrededor de 1 fpkm, y no de cinco como esta ahora
load(file = "data/TauExp_Scores.RData")

############################
## Polycomb Scores
load(file = "data/RTS_PolyC_score.RData")

####################################
##Gene relationship with phenotypes. Either from OMIM or Mammalian Phenotype Ontology
##based on OMIM and human phenotypes. Table
##3410 gene-phenotype relationships based on OMIM and HPO
load(file ="data/gene_fenotipe_basedOn_hpoAndOmim.RData")
humanBased_genePhenotype<-gene_fenotipe_table
rm(gene_fenotipe_table)

##Human orthologues,Based on Mice genes and Mice Phenotypes
##10471 genes
load(file ="data/gene_fenotipe_basedOn_mgi_and_MP.RData")
miceBased_genePhenotype<-gene_fenotipe_table
rm(gene_fenotipe_table)

####Podriamos dar la informacion de para cada gen en OMIM su OMIM entry y lo mismo para HPO, o su link...

############################################
####Haploinsufficiency scores Nature 2016
####"the probability of being loss-of-function intolerant (intolerant of both heterozygous and homozygous lof variants)"
load("data/nature_lof_scores.RData")

####Haploinsufficiency scores Huang 2010
load("data/Huang_hi_scores.RData")

##Haploinsufficiency from ClinGene
load("data/clingen_hi_scores.RData")

##Haploinsufficiency scores from Cell paper
load("data/cell_hi_scores.RData")

##Triplosensitivity scores from Cell paper
load("data/cell_triploSens_scores.RData")

##Loading centromere positions
load("data/hg19_centromerePositions.RData")

###################################
## Assembling MultiDataList object
###################################
# ls()
MultiDataList<-list("cellInfo_hi"=cellInfo_hi,
                    "cellInfo_triploSens"=cellInfo_triploSens,
                    "clinGen_hiInfo"=clinGen_hiInfo,
                    "gtf_annotation"=gtf_annotation,
                    "huangScores"=huangScores,
                    "humanBased_genePhenotype"=humanBased_genePhenotype,
                    "miceBased_genePhenotype"=miceBased_genePhenotype,
                    "nature_lof_scores"=nature_lof_scores,
                    "rts_allGenes"=rts_allGenes,
                    "tau_exp_scores"=tau_exp_scores,
                    "centromerePos"=centromerePos)

##Saved in POSTRE data folder
save(MultiDataList,
     file = "data/MultiDataList.RData" )

save(MultiDataList,
     file = "~/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/Postre/Postre_app/data/MultiDataList.RData" )

