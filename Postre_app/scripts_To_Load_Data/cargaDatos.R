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

# ############################################################
# ##Poorly expressed genes, to know if a gene expressed or not
# ##same object of Candidates Gain of Function, when not filtering by polycomb
# ##genes with expr levels 0-1 for the quartile normalized or <11cpm palate (approx 5 fpkm)
# load(file = "data/poorlyExpressedGenes.RData")


#########################################################
## Loading TauExp values: es decir, que si en algun tejido fetal o NCC expression considerable su TAU
## de lo contrario un 0 de TAU. (es decir si su expression es muy baja) en todos los estadios
## ajustar para alrededor de 1 fpkm, y no de cinco como esta ahora
load(file = "data/TauExp_Scores.RData")

###########################################################################################
## POLYCOMB GENES (merging the ones from different conditions as in hESC and hNCCs)
load(file = "data/Human_PolycombGenes.RData")
polyCombGenes<-human_polyComb_genes
rm(human_polyComb_genes)

############################
## Polycomb Scores
load(file = "data/RTS_PolyC_score.RData")

# ##delete or comment next line upon testing NeoTAD
# polyCombGenes<-c(polyCombGenes, "NPHP1","MALL")

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







