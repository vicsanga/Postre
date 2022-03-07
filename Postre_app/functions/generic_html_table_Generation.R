###################################################################
##Code to generate generic html table for heatmap numbers at least
###################################################################
generic_table_html_generation<-function(targetMatrix){

  ##Put affected gene on the first column
  tableCols<-colnames(targetMatrix)
  
  #############################################################
  ##The first column is the one we use to PARSE the R Table
  ##############################################################
  
  if("affected_gene" %in% tableCols){
    ##Resort 
    targetMatrix<-targetMatrix[,c("affected_gene",tableCols[tableCols!="affected_gene"])] 
  }
  
  table_content<-paste("<div class= divTablePathoMech>",
                       "<table class='tablePathomech'><thead><tr>",
                       sep="")
  
  
  #First, table header
  for(nameCol in colnames(targetMatrix)){
    table_content<-paste(table_content,
                         "<th>",
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
  for(nRowMatrix in 1:nrow(targetMatrix)){
    
    ##Adding table row
    table_content<-paste(table_content,"<tr>",
                         sep="",
                         collapse = "")
    
    for(nameCol in colnames(targetMatrix)){
      
      cellContent<-targetMatrix[nRowMatrix,nameCol]
      
      #Adding column in the row
      table_content<-paste(table_content,
                           "<td>",
                           cellContent,
                           sep="",
                           collapse = "")
      
      ##########################
      #Ending column in the row
      table_content<-paste(table_content,"</td>",
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