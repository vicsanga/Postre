#############################################
## Script to develop the phenoScore function

phenoScore_fun<-function(geneData, gof_case){
  
  ##Introduction of gof_case parameter  6/Sept/2021
  ##If gof_case == FALSE, prediction for LOF cases (we look for gene-pheno specific match)
  ##If gof_case == TRUE, prediction for GOF cases, we are not as stringent, we do not look for gene-pheno match
  ##This only used for the high sensitivity scenario
  
  ##2 components: reportedPhenoInDb + associatedMainPhenoInDb
  
  ##Check if reported phenotype either in OMIM or in MGI
  reportedPhenoInDb<-rowSums(geneData[,c("associatedPhenotypeIn_OMIM","associatedPhenotypeIn_MGI")])
  
  #reportedPhenoInOMIM<-geneData[,"associatedPhenotypeIn_OMIM"]
  
  ##Check if association to the main phenotype either in OMIM or in MGI
  associatedMainPhenoInDb<-rowSums(geneData[,c("mainPhenotype_Through_OMIM_Human","mainPhenotype_Through_MGI_Mice")])
  
  associatedMainPhenoInOMIM<-geneData[,c("mainPhenotype_Through_OMIM_Human")]
  associatedMainPhenoInMGI<-geneData[,c("mainPhenotype_Through_MGI_Mice")]
  
  ###########################################
  phenoScore<-0

  if(associatedMainPhenoInOMIM == TRUE ){
    ##So, associated with phenotype in Human, Top Relevant
    ##Does not matter whether it is not in mice, Human, Top Relevant specie
    phenoScore<-1
  }else{
    ##Combine info Mouse and Human, so if only in Mice, boost a bit the scoring (0.5 phenoScore)
    #phenoScore<-associatedMainPhenoInDb/2    
    ##DOING NOTHING, ONLY CONSIDERING OMIM, FOR HIHG-SPECIFICITY CHECK CHANGE IN FP FN
  }

  
  return(phenoScore)
  
}