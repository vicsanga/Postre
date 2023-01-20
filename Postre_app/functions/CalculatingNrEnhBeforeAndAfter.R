###############################################################################
## Function to compute the nr of enh kept and gained in the rearrangement
## as well as the total enhancer acetilation gained or lost
## It is valid, at least, for inversions and translocations
## It is also computed the nr of enhancers gained and kept per side individually 
## For each gene, what will be used for the graphical representations
###############################################################################

calculatingEnh<-function(targetGenes,keptEnh,gainedEnh,matrixRes,enhOrigin, info_affectedGenes){
  ##if there is at least a target gene
  ##compute the balance of enh, and add the info to the matrix of results
  ##if there are no target genes, waste of time computing the enh balances
  ##since nothing will be added to the result matrix
  ##so if there are no targets, DoNothing, simply return the results matrix
  if(length(targetGenes)>0){

    ####
    ##In this loop overall information about gained and lost enhancers is introduced 
    for(origin in enhOrigin){
      
      ###REGARDING KEPT ENHANCERS AND ACETILATION
      if(is.null(keptEnh)==FALSE){
        ##so there are kept enh(from at least 1 source(from 1 developmental stage))
        subs_kept_enh<-subset(keptEnh,source==origin)
        n_kept_enh<-nrow(subs_kept_enh)    
        acetilation_kept_enh<-sum(subs_kept_enh$acetilation)
        
      }else{
        ##there are not so the number is 0
        ##if there are not, we can not subset a NULL object, that rises an error
        n_kept_enh<-0
        acetilation_kept_enh<-0
      }
      
      ###REGARDING GAINED ENHANCERS AND ACETILATION
      if(is.null(gainedEnh)==FALSE){
        ##so there are gained enh (from at least 1 source(from 1 developmental stage))
        subs_gained_enh<-subset(gainedEnh,source==origin)
        n_gained_enh<-nrow(subs_gained_enh)
        acetilation_gained_enh<-sum(subs_gained_enh$acetilation)
      }else{
        n_gained_enh<-0
        acetilation_gained_enh<-0
      }
      
      #####################
      ## Saving information
      ######################
      
      #################
      ##For KEPT enh
      ###Colname in the results matrix
      colId<-paste0("nEnhancers_kept_",origin)
      ##Add to all the genes in the domain the n enh for this stage
      matrixRes[targetGenes,colId]<-n_kept_enh
      ##Add acetilation levels info
      colAcetilLevelId<-paste0("enhancers_acetilation_kept_",origin)
      matrixRes[targetGenes,colAcetilLevelId]<-acetilation_kept_enh
      
      ##################
      ##For GAINED enh
      ###Colname in the results matrix
      colId<-paste0("nEnhancers_gained_",origin)
      ##Add to all the genes in the domain the n enh for this stage
      matrixRes[targetGenes,colId]<-n_gained_enh
      ##Add acetilation levels info
      colAcetilLevelId<-paste0("enhancers_acetilation_gained_",origin)
      matrixRes[targetGenes,colAcetilLevelId]<-acetilation_gained_enh
      
    }
    
    ##################
    ## In this loop we are going to add the information regarding nEnh to the left & nEnh to the right kept and gained
    ## With respect to each gene
    
    for(gene in targetGenes){
      geneTSS<-info_affectedGenes$genesPosition[gene,"TSS"]
      
      for(origin in enhOrigin){
        ###REGARDING KEPT ENHANCERS AND ACETILATION
        if(is.null(keptEnh)==FALSE){
          ##so there are kept enh(from at least 1 source(from 1 developmental stage))
          subs_kept_enh<-subset(keptEnh,source==origin)
          
          n_kept_enh_toTheRight<-sum(subs_kept_enh$start > geneTSS)
          n_kept_enh_toTheLeft<-sum(subs_kept_enh$start < geneTSS)
          
          acetilationEnh_kept_ToThe_Right<-sum(subs_kept_enh$acetilation[subs_kept_enh$start > geneTSS])
          acetilationEnh_kept_ToThe_Left<-sum(subs_kept_enh$acetilation[subs_kept_enh$start < geneTSS])
          
        }else{
          ##there are not so the number is 0
          ##if there are not, we can not subset a NULL object, that rises an error
          n_kept_enh_toTheRight<-0
          n_kept_enh_toTheLeft<-0
          
          acetilationEnh_kept_ToThe_Right<-0
          acetilationEnh_kept_ToThe_Left<-0
          
        }
        
        #####################
        ## Saving information
        ######################
        
        #################
        ##For KEPT enh
        ###Colname in the results matrix
        colId<-paste0("nEnh_ToTheLeft_kept_",origin)
        ##Add to all the genes in the domain the n enh for this stage
        matrixRes[gene,colId]<-n_kept_enh_toTheLeft
        
        colId<-paste0("nEnh_ToTheRight_kept_",origin)
        ##Add to all the genes in the domain the n enh for this stage
        matrixRes[gene,colId]<-n_kept_enh_toTheRight
        
        
        colId<-paste0("enhancers_acetilation_ToTheLeft_kept_",origin)
        ##Add to all the genes in the domain the n enh for this stage
        matrixRes[gene,colId]<-acetilationEnh_kept_ToThe_Left
        
        colId<-paste0("enhancers_acetilation_ToTheRight_kept_",origin)
        ##Add to all the genes in the domain the n enh for this stage
        matrixRes[gene,colId]<-acetilationEnh_kept_ToThe_Right
        
      }
      
    }
    
  }

  return(matrixRes)
}