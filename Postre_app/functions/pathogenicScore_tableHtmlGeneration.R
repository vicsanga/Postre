#######################################################
## Function to create the nice pathogenic score table
#######################################################
pathogenicScore_tableHtmlGeneration<-function(patientResults, targetMatrix){
  
  #Basarme en generic_table_html_generation()
  
  ##Puede haber NAs si para cada cell type se usa un mapa de TADs diferente, ej paciente cardiovasuclar NIJ18 creo con PROX1 gene como patogenico
  
  allPhases<-names(patientResults$resultsPerPhase_secondaryInfo)
  
  ##Starting table
  table_content<-paste("<div class='divTablePathoMech'>",
                       "<table class='tablePathomech'",
                       "style='width:",
                       ##DEFINIR AMPLITUD COMO NPHASES*X + CTE en referencia a columna del gen
                       as.character(250 + 700*length(allPhases)),##Defining in pixels (when no ending)
                       "px !important'><thead><tr>",
                       sep="")
  # tableColumns<-c("Gene", allPhases)

  #First, table header
  for(nameCol in colnames(targetMatrix)){
    table_content<-paste(table_content,
                         if(nameCol=="Gene"){
                           ##.column_sticky has a zindex of 9000 to avoid gene (eg TFAP2A )celd to be hidden, adding a higher z-index to avoid header cell ("Gene") to be hidden
                           "<th style='z-index: 9999;' class='column_sticky'>"
                         }else{
                           "<th>"
                         },
                         nameCol,
                         "</th>",
                         sep = "",
                         collapse = "")
  }
  table_content<-paste(table_content,
                       "</tr></thead><tbody>",##Div to allow verticall scroll mantaining fixed header
                       sep="",
                       collapse = "")
  
  #####################################
  #Filling Main cells content
  for(targetGene in targetMatrix$Gene){
    
    ##Adding table row
    table_content<-paste(table_content,"<tr>",
                         sep="",
                         collapse = "")
    
    
    #
    
    ##Here we will be dealing with first column of each row
    ##SO, add the first cell with gene NAME
    #Adding gene information
    #Adding column in the row
    table_content<-paste(table_content,
                         "<td class='column_sticky'>", ##Column sticky to have it frozen
                         "<i>",
                         targetGene,
                         "</i>",
                         "</td>",
                         sep="",
                         collapse = "")
    
    ##Adding pathogenic score columns  
    ##La info de interes se encuentra en: patientResults$resultsPerPhase_secondaryInfo$NeuralCrestEarly$`H1-ESC_Dixon2015-raw_TADs.txt`$NeuralCrestEarly
    for(phase in allPhases){
      
      #Solo estamos usando 1 mapa de tads por phase, por ello el [[1]]
      infoCurrentPhase<-patientResults$resultsPerPhase_secondaryInfo[[phase]][[1]][[phase]]
      targetInfo<-infoCurrentPhase[targetGene,]

      ps_data<-list()
      for(typePS in c("LOF","GOF")){
        ##Gathering pathogenic score data
        ##IMPORTANT to keep order to always add first score for Lof then Gof
        
        if(typePS == "LOF"){
          
          finalScore<-targetInfo$LOF_score
          geneEnhancerScore<-targetInfo$geneEnhancerScore_LOF
          genePhenoScore<-targetInfo$genePhenoScore_LOF
          geneFeaturesScore<-targetInfo$geneFeaturesScore_LOF
          dosageSensitivityScore<-targetInfo$dosageSensitivityScore_LOF
          polycombScore<-targetInfo$polycombScore_LOF
          geneExpressionScore<-targetInfo$geneExpressionScore_LOF
          
        }else{
          ##So, GOF
          
          finalScore<-targetInfo$GOF_score
          geneEnhancerScore<-targetInfo$geneEnhancerScore_GOF
          genePhenoScore<-targetInfo$genePhenoScore_GOF
          geneFeaturesScore<-targetInfo$geneFeaturesScore_GOF
          dosageSensitivityScore<-targetInfo$dosageSensitivityScore_GOF
          polycombScore<-targetInfo$polycombScore_GOF
          geneExpressionScore<-targetInfo$geneExpressionScore_GOF
          
        }
        
        ##Adding data, try to do it by introducing a table inside of cell  let's see how this goes
        # <table><tr><td>split 1</td><td>split 2</td></tr></table>
        
        targetWord<-"Ignored"
        
        cellContent<-paste0("<td>",
                            "<b>",
                            "Pathogenic Score ",
                            typePS,
                            ":",
                            round2(finalScore,digits = 2),
                            "</b>",
                            "<ul>",
                              "<li>","genePhenoScore: " ,
                              if(is.na(genePhenoScore)){
                                targetWord
                              }else{
                                round2(genePhenoScore,digits = 2)
                              },"</li>",
                              "<li>","geneEnhancerScore: " ,
                              if(is.na(geneEnhancerScore)){
                                targetWord
                              }else{
                              round2(geneEnhancerScore,digits = 2)
                              },
                              "</li>",
                              "<li>","geneExpressionScore: ",
                            if(is.na(geneExpressionScore)){
                              targetWord
                            }else{
                              round2(geneExpressionScore,digits = 2)
                            },
                            "</li>",
                              "<li>","geneFeaturesScore: ",
                            if(is.na(geneFeaturesScore)){
                              targetWord
                            }else{
                              round2(geneFeaturesScore,digits = 2)
                            },
                            "</li>",
                              "<ul>",
                                "<li>","dosageSensitivityScore: ",
                            if(is.na(dosageSensitivityScore)){
                              targetWord
                            }else{
                              round2(dosageSensitivityScore,digits = 2)
                            },
                            "</li>",
                                "<li>","polycombScore: ",
                            if(is.na(polycombScore)){
                              targetWord
                            }else{
                              round2(polycombScore,digits = 2)
                            },
                            "</li>",
                              "</ul>",
                            "</ul>",
                            "</td>")
        
        ##Recording data
        ps_data[[typePS]]<-cellContent
        
      }
      
      ##Adding data, two columns in the cell
      # <table><tr><td>split 1</td><td>split 2</td></tr></table>
      table_content<-paste(table_content,
                           "<td>",
                           "<table class='nestedPathomechTables'><tr>",
                           ps_data$LOF,
                           ps_data$GOF,
                           "</tr></table>",
                           "</td>",
                           sep="",
                           collapse = "")
    }
    
    ##Ending table row
    table_content<-paste(table_content,"</tr>",
                         sep="",
                         collapse = "")
    
  }
  
  ###Closing Table
  table_content<-paste(table_content,
                       "</tbody>",
                       "</table>",
                       "</div>",##closing <div class= divTablePathoMech>
                       sep = "")
  

  #############################
  return(table_content)
  
}