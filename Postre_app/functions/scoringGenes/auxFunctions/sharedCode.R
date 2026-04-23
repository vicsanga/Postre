###Shared code to minimize absurdly high number of CODE duplications in RankingGenes script
#print("Running code minimizing duplications in: RankingGenes script !!!!")
#print("evaluate whether adding after ifs at some point in RankingGenes script !!!!")

##Get max score && type
Scores<-c(lof_score, gof_score)
names(Scores)<-c("LOF_score","GOF_score")

max_score<-max(Scores)## (unless both equal check, like 0,0 or 0.8 0.8)

if(lof_score==gof_score){
  typeMax<-"equally_likely"
}else{
  ##one score bigger than the other
  typeMax<-names(Scores)[which(Scores==max_score)]
}

##THE FOLLOWING DATA ABNORMALITIES HANDLING (i.e. missing expressing data)
#ONLY MAKES SENSE IF NOT CellTypeAgnostic (in CellTypeAgnostic we do not consider exp and enh data)
if(phase != "CellTypeAgnostic"){
  
  ##if there is no expression data for the gene, we do not compute anything
  ##We put a -1 to place it at the bottom
  if(matrixPhase[,"FPKM"]==-1){
    max_score<-0
    typeMax<-"equally_likely"
    
    gof_score<-0
    lof_score<-0
    
  }
  
  ##if there multiple expression values, distant from them (eg 0.5 and 2 FPKM)
  ##We put a -2 to place it at the bottom
  if(matrixPhase[,"FPKM"]==-2){
    max_score<-0
    typeMax<-"equally_likely"
    
    gof_score<-0
    lof_score<-0
    
  }
  
}


##Add results
# colnames(matScores)<-c("GOF_score","LOF_score","Max_type","Max_Score","type")
res_scores[[phase]][gene,"GOF_score"]<-gof_score
res_scores[[phase]][gene,"LOF_score"]<-lof_score

res_scores[[phase]][gene,"Max_type"]<-typeMax       
res_scores[[phase]][gene,"Max_Score"]<-max_score 

##RECORDING METADATA/Additional info pathogenic score
res_scores[[phase]][gene,"genePhenoScore_GOF"]<-gof_score_metadata$genePhenoScore
res_scores[[phase]][gene,"genePhenoScore_LOF"]<-lof_score_metadata$genePhenoScore

res_scores[[phase]][gene,"geneEnhancerScore_GOF"]<-gof_score_metadata$geneEnhancerScore
res_scores[[phase]][gene,"geneEnhancerScore_LOF"]<-lof_score_metadata$geneEnhancerScore

res_scores[[phase]][gene,"geneFeaturesScore_GOF"]<-gof_score_metadata$geneFeaturesScore
res_scores[[phase]][gene,"geneFeaturesScore_LOF"]<-lof_score_metadata$geneFeaturesScore

res_scores[[phase]][gene,"dosageSensitivityScore_GOF"]<-gof_score_metadata$dosageSensitivityScore
res_scores[[phase]][gene,"dosageSensitivityScore_LOF"]<-lof_score_metadata$dosageSensitivityScore

res_scores[[phase]][gene,"polycombScore_GOF"]<-gof_score_metadata$polycombScore
res_scores[[phase]][gene,"polycombScore_LOF"]<-lof_score_metadata$polycombScore

res_scores[[phase]][gene,"geneExpressionScore_GOF"]<-gof_score_metadata$geneExpressionScore
res_scores[[phase]][gene,"geneExpressionScore_LOF"]<-lof_score_metadata$geneExpressionScore