#################################
## Function Table HTML generation
#################################
table_html_generation<-function(targetMatrix, name_targetMatrix,targetPheno, ids_append){
  ##name_targetMatrix && targetPheno used to create the specific js functions for the search box
  ##Put affected gene on the first column
  tableCols<-colnames(targetMatrix)
  
  #############################################################
  ##The first column is the one we use to PARSE the R Table
  ##############################################################
  
  if("affected_gene" %in% tableCols){
    ##Resort 
    targetMatrix<-targetMatrix[,c("affected_gene",tableCols[tableCols!="affected_gene"])] 
  }
  
  ##If we are on the patient info table, it already has the patients in the first column
  
  idTable<-paste("table_",
                 name_targetMatrix,
                 "_",
                 targetPheno,
                 ids_append,
                 sep="")
    
  idInputClick<-paste("inputClick_",
                      name_targetMatrix,
                      "_",
                      targetPheno,
                      ids_append,
                      sep="")
  
  idFunction<-paste("function_",
                    name_targetMatrix,
                    "_",
                    targetPheno,
                    ids_append,
                    sep="")
  
  ### Creating Table HTML part
  #The previously specified format is linked to "heatmap" id, hence indicating that table belongs to this category
  
  
  searchCode<-paste("<input type='text' id='",
                    idInputClick,
                    "' onkeyup='",
                    idFunction,
                    "()' placeholder='Search for...' title='Type in a name'> <br>",
                    sep=""
  )
  
  ##Adding style='white-space:nowrap;' to prevent breaking header names into multiple lines eg chr (in one line) breakp1 (in different line)
  table_content<-paste(searchCode,
                       "<br>",
                       "<div class= divTablePathoMech>",
                       "<table class='tablePathomech' style='white-space:nowrap;' id='",
                       idTable,
                       "'><thead><tr>",
                       sep="")
  
  ###############################################################################################
  ##Dataframe for reformating for html varible names for Headers, substituting _ by spaces etc.
  ###############################################################################################
  valueInR<-c("affected_gene", "Num_SVs", "SV_IDs", "DirectEffectLOF","LongRangeLOF","DirectEffectGOF", "LongRangeGOF", "SV_ID", "chr_Break1","coord_Break1", "chr_Break2","coord_Break2", "TypeSV","source")                                                   
  valueInHTML<-c("Gene", "N SVs","SV IDs", "Coding LOF", "Long-Range LOF", "Coding GOF", "Long-Range GOF","SV ID","Chr Breakp1", "Coord Breakp1", "Chr Breakp2","Coord Breakp2", "Type SV","Source")
  
  dfDat<-data.frame("valueInR"=valueInR, "valueInHTML"=valueInHTML)
  rownames(dfDat)<-valueInR
  
  #First, table header
  for(nameCol in colnames(targetMatrix)){
    table_content<-paste(table_content,
                         "<th>",
                         
                         if(nameCol %in% rownames(dfDat)){
                           ##Formatting variable name to HTML more pleasing style
                           dfDat[nameCol,"valueInHTML"]
                         }else{
                           nameCol
                         },
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
  
  ##Adding js function
  jsFun<-paste("<script>
  function ",
               idFunction,
    "() {
    var input, filter, table, tr, td, cell, i, j;
    input = document.getElementById('",
    idInputClick,
    "');
    filter = input.value.toUpperCase();
    table = document.getElementById('",
    idTable,
    "');
    tr = table.getElementsByTagName('tr');
    for (i = 1; i < tr.length; i++) {
      /**Hide the row initially.**/
      tr[i].style.display = 'none';
      
      td = tr[i].getElementsByTagName('td');
      for (j = 0; j < td.length; j++) {
        cell = tr[i].getElementsByTagName('td')[j];
        if (cell) {
          if (cell.innerHTML.toUpperCase().indexOf(filter) > -1) {
            tr[i].style.display = '';
            break;
          } 
        }
      }
    }
  }
</script>",
               sep="")
  
  jsFun<-gsub("[\r\n]", "", jsFun)
  
  #########
  #########
  table_content<-paste(table_content,
                       jsFun,
                       sep = "")
  
  
  #############################
  return(table_content)
}