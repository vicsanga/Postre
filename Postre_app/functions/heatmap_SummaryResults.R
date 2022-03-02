##############################################################################
## This function and scripts centralizes the MAIN RESULTS section
## Function to Create the Master Heatmap
## A summary figure to illustrate graphically and simply the biggest changes
##############################################################################

##Required Function to plot number behind heatmap Table Generation
source(file = "functions/generic_html_table_Generation.R", local = TRUE)

heatmap_summaryResults<-function(patientResults, minRequiredScore, highScore){
  ##highScore deprecated, whtvr above the minRequiredScore gets top intensity color
  
  targetMatrix<-patientResults$masterSummaryResultsMatrix
  
  ##Filter here per min pathogenic score, if considered
  
  ##MaxScore we are not going to use it after potentially filtering point
  targetMatrix$maxScore<-NULL
  
  
  ##Tracking relevant Gene-Phase for report:
  #We probably need to track also, gene-phase and mechanism
  relevant_Gene_Phase<-character()
  
  ##To control the hide/show reports
  ##Due to my lack of knowledge
  ##we are going to create a javascript function per gene report to target specifically its toggle behaviour
  all_functionsJS<-"<script>"##append them here
  
  ## Starting HTML
  header_html<-c("<html><head><style>#heatmap {
    font-family: 'Trebuchet MS', Arial, Helvetica, sans-serif;
    border-collapse: collapse;
    width: 100%;
    }

    #heatmap td, #heatmap th {
      border: 1px solid #ddd;
      padding: 8px;
    }

    #heatmap tr:nth-child(even){background-color: #f2f2f2;}

    #heatmap tr:hover {background-color: #ddd;}

    #heatmap th {
      padding-top: 12px;
      padding-bottom: 12px;
      text-align: left;
      background-color: #0066cc;
      color: white;
    }</style></head><body>")
  
  ##
  ##Let's remove all line breaks, introduced from internet copypaste java & html code
  header_html<-gsub("[\r\n]", "", header_html)
  
  
  ##Adding Title
  header_html<-paste(header_html,"<div class='heatmapSection'>",
                     "<h1 id='mainResults'>Results Overview</h1>",
                     sep = "",
                     collapse = "")
  
  
  ##Defining cell colors
  #GOF, TOPGOF, LOF, TOPLOF, Default
  
  col_default<-"#a6a6a6"#grey
  # col_directImpact_NotPathogenic<-"#DCDCDC"
  
  col_lof<-"#ff704d"
  col_top_lof<-"#ff3300" 
  
  col_gof<-"#b3ffb3"
  col_top_gof<-"#00ff00" 
  
  #col_mix<-"#80dfff" ##In case there is a situation (as a result of the thresholds) where both scenarios are likely to occur with a high score
  ##Maybe more realistic for phase free, just in case, implement
  #col_top_mix<-"#00bfff" ##In case there is a situation (as a result of the thresholds) where both scenarios are likely to occur with a high score
  
  ### Creating Table HTML part
  #The previously specified format is linked to "heatmap" id, hence indicating that table belongs to this category
  table_content<-"<table id='heatmap'><tr>"
  #First, table header
  for(nameCol in colnames(targetMatrix)){
    
    ######################################
    ## Formatting properly table header
    ######################################
    
    ##Do not want to rename right now, in case we rerename afterwards..
    ##And this one is quite used in the app internal working to distinguish in ranking between using all genomic data or not
    if(nameCol=="phaseFree"){
      table_content<-paste(table_content,
                           "<th>",
                           "Cell Type Independent",
                           "</th>",
                           sep = "",
                           collapse = "")
      
    }else if(nameCol=="affected_gene"){
      table_content<-paste(table_content,
                           "<th>",
                           "Gene",
                           "</th>",
                           sep = "",
                           collapse = "")
    }else if(nameCol=="GeneImpact"){
      table_content<-paste(table_content,
                           "<th>",
                           "Gene Impact",
                           "</th>",
                           sep = "",
                           collapse = "")
    }else{
      table_content<-paste(table_content,
                           "<th>",
                           nameCol,
                           "</th>",
                           sep = "",
                           collapse = "") 
    }
  }
  table_content<-paste(table_content,
                       "</tr>",
                       sep="",
                       collapse = "")
  
  #Filling Main cells content
  for(gene in targetMatrix$affected_gene){
    
    
    geneImpact<-targetMatrix[gene,"GeneImpact"]
    
    ##Adding table row
    table_content<-paste(table_content,"<tr>",
                         sep="",
                         collapse = "")
    
    for(nameCol in colnames(targetMatrix)){
      
      if(nameCol=="affected_gene"){
        #Adding column in the row
        table_content<-paste(table_content,
                             "<td>",
                             gene,
                             sep="",
                             collapse = "")
        
      }else if(nameCol=="GeneImpact"){
        
          ########################################
          ## Formatting GeneImpact content column          
          ########################################
        
          preFormating_cellContent<-targetMatrix[gene,nameCol]
          
          #postFormating_cellContent
          if(preFormating_cellContent == "LongRange"){
             postFormating_cellContent<-"Long-Range"
            
          }else if(preFormating_cellContent == "Direct_geneTruncation"){
            postFormating_cellContent<-"Gene Truncation"
            
          }else if(preFormating_cellContent == "Direct_geneDeletion"){
            postFormating_cellContent<-"Gene Deletion"
            
          }else if(preFormating_cellContent == "LongRange_geneDuplication"){
            postFormating_cellContent<-"Gene Duplication and Long-Range"
            
          }else if(preFormating_cellContent == "Direct_geneDuplication"){
            postFormating_cellContent<-"Gene Duplication"
            
          }else if(preFormating_cellContent == "Direct_LongRange_geneDuplication"){
            postFormating_cellContent<-"Gene Duplication and Long-Range"
            
          }else{
            ##In case I'm forgetting to format some term
            postFormating_cellContent<-preFormating_cellContent
          }
          
          #Adding column in the row
          table_content<-paste(table_content,
                               "<td>",
                               postFormating_cellContent,
                               sep="",
                               collapse = "")
          
      }else{
        ##Up to this point we are left with the different phase data columns
        ##get TargetMatrix cell content, and retrieve score and type of Mechanism
        
        cellContent<-targetMatrix[gene,nameCol]
        ##Check that the cell content is not NA, if is like that 
        ##Set parameters to paint the cell grey
        if(is.na(cellContent)){
          ##Set everything to -1
          ##We track both LOF score and GOF score, since oct 2021 Modification
          # geneScore<-(-1)
          geneMechanism<-"notConsidered"
          
          geneScore_LOF<-(-1)
          geneScore_GOF<-(-1)
          bothGeneScore<-c(geneScore_LOF,geneScore_GOF)
          names(bothGeneScore)<-c("LOF","GOF")
          
        }else{
          ##so the cell content is not a NA

          # cellContent<-unlist(strsplit(x = cellContent,
          #                              split = ";",
          #                              fixed = TRUE))
          
          info_LOF<-unlist(strsplit(x = cellContent,
                                   split = "--",
                                   fixed = TRUE))[1]
          
          info_GOF<-unlist(strsplit(x = cellContent,
                                   split = "--",
                                   fixed = TRUE))[2]
          
          info_LOF_splitted<-unlist(strsplit(x = info_LOF,
                                             split = ";",
                                             fixed = TRUE))
          
          info_GOF_splitted<-unlist(strsplit(x = info_GOF,
                                             split = ";",
                                             fixed = TRUE))
          
          ##PILLAR SCORE PARA GOF Y PARA LOF
          ##AHORA TENEMOS DOBLE INFO EN ESTA CELDA QUE ANTES
          
          geneScore_LOF<-as.numeric(info_LOF_splitted[1])
          geneScore_GOF<-as.numeric(info_GOF_splitted[1])
          bothGeneScore<-c(geneScore_LOF,geneScore_GOF)
          names(bothGeneScore)<-c("LOF","GOF")
          ##geneMechanism<-cellContent[2]##Depende de situacion 
        }
        
        #####################################
        #Adding score pathogenic information
        #####################################
        
        ##Situation where pathogenic evidence for LOF OR GOF but NOT both
        if ( sum(bothGeneScore >= minRequiredScore) == 1 ){
          
          ################################################################################
          ##Hence gene presents a patogenic score higher or equal the minimum relevant
          ##So relevant cell
          ##And geneMechanism is GOF or LOF but not both
          ################################################################################
          
          ##Get targetMechanism responsible of pathogenic score
          targetMech<-names(bothGeneScore)[bothGeneScore >= minRequiredScore]
          ##Add to the relevant info vector
          ## RELEVANT MATCH is fundamental, used on many things downstream as id
          relevantMatch<-paste(gene,"_",targetMech,"_",nameCol,sep="")
          
          relevant_Gene_Phase<-c(relevant_Gene_Phase,
                                 relevantMatch)
          
          ##Creating javascript function to hide/show the gene report specifically
          
          functionsJS<-paste(
            "function myFunction_Tablecell_",
            relevantMatch,
            "() {
            var x = document.getElementById('",
            relevantMatch,
            "');
            if (x.style.display === 'none') {
            x.style.display = 'block';
            } else {
            x.style.display = 'none';
            }
        }",
            sep="",
            collapse = "") 
          ##remove linebreaks
          functionsJS<-gsub("[\r\n]", "", functionsJS)
          
          ##Add to allJs functions, and append them after the table
          all_functionsJS<-paste(all_functionsJS,
                                 functionsJS,
                                 sep="",
                                 collapse = "")
          
          ##The gene report will have also a button tho hidde the content
          ##But the gene cell will also be able to close the report
          ##myFunction_report_XX
          
          ##Preparing cell content
          ##It will be the geneMechanism as text, but linking the geneReport
          geneMechanism<-targetMech
          geneCellContent<-paste("<a href='#",
                                 relevantMatch, ##gene_phase
                                 "'",
                                 "onclick='myFunction_Tablecell_",
                                 relevantMatch,
                                 "()' ",
                                 "style='color:#000000;'>",
                                 "<div style='height:100%; width:100%' title='Click for more details'>",##To make pointer sensitive to the whole cell ##To add mouseover text
                                 geneMechanism, ##LOF, GOF
                                 "</div>",
                                 "</a>",
                                 sep="",
                                 collapse = "")
          
          
          
          
          
          if(geneMechanism=="LOF"){
            
            table_content<-paste(table_content,
                                 "<td style='background-color:",
                                 col_top_lof,
                                 "'",
                                 ">",
                                 geneCellContent,##Cell Content
                                 sep="",
                                 collapse = "")
          }else if(geneMechanism=="GOF"){
            
            table_content<-paste(table_content,
                                 "<td style='background-color:",
                                 col_top_gof,
                                 "'",
                                 ">",
                                 geneCellContent,##Cell Content
                                 sep="",
                                 collapse = "")
          }
          
        }else if ( sum(bothGeneScore >= minRequiredScore) == 2 ){
          ##So, for this gene there is evidence for LOF AND GOF
          ##LOF because of loss of cognate enhancers
          ##GOF because of gain of multiple ectopic enhancers
          
          ##We have to split the heatmap cell in 2. To provide 2 links, one for GOF report, the other for LOF report
          
          ##Same procedure that at single pathomech level, but twice (just looped)
          ##Get targetMechanism responsible of pathogenic score
          ##In this case targetMechanism is both, GOF and LOF
          
          for(targetMech in c("LOF","GOF")){
            
            ##print(targetMech) ##Get targetMechanism responsible of pathogenic score
            
            ##Add to the relevant info vector
            ## RELEVANT MATCH is fundamental, used on many things downstream as id
            relevantMatch<-paste(gene,"_",targetMech,"_",nameCol,sep="")
            
            relevant_Gene_Phase<-c(relevant_Gene_Phase,
                                   relevantMatch)
            
            ##Creating javascript function to hide/show the gene report specifically
            
            functionsJS<-paste(
              "function myFunction_Tablecell_",
              relevantMatch,
              "() {
            var x = document.getElementById('",
              relevantMatch,
              "');
            if (x.style.display === 'none') {
            x.style.display = 'block';
            } else {
            x.style.display = 'none';
            }
        }",
              sep="",
              collapse = "") 
            ##remove linebreaks
            functionsJS<-gsub("[\r\n]", "", functionsJS)
            
            ##Add to allJs functions, and append them after the table
            all_functionsJS<-paste(all_functionsJS,
                                   functionsJS,
                                   sep="",
                                   collapse = "")
            
            ##The gene report will have also a button tho hidde the content
            ##But the gene cell will also be able to close the report
            ##myFunction_report_XX
            
            #####################################################################
            ##Preparing cell content
            ##It will be the geneMechanism as text, but linking the geneReport
            geneMechanism<-targetMech
            geneCellContent<-paste("<a href='#",
                                   relevantMatch, ##gene_phase
                                   "'",
                                   "onclick='myFunction_Tablecell_",
                                   relevantMatch,
                                   "()' ",
                                   "style='color:#000000;'>",
                                   "<div style='height:100%; width:100%' title='Click for more details'>",##To make pointer sensitive to the whole cell ##To add mouseover text
                                   geneMechanism, ##LOF, GOF
                                   "</div>",
                                   "</a>",
                                   sep="",
                                   collapse = "")
            
            
            ######################
            ## In the subtable, LOF is the first column and GOF the second one
            ## So in LOF add the start of the TABLE and in GOF the end of it
            ## based on: <td><table><tr><td>split 1</td><td>split 2</td></tr></table></td> 
            if(geneMechanism=="LOF"){
              
              table_content<-paste(table_content,
                                   ##Styling for the nestedHeatmap table provided in: MainInterfaceStyling.html
                                   ##Adjusting with an inline tag the td style to padding:0px for the cell that contains the nested table, so that the nestedtable can occupy the whole cell space
                                   "<td style='padding:0px;'><table id='nestedHeatmap'><tr>",##OPENING SUBTABLE TO SUBDIVIDE THE CELL
                                   "<td style='background-color:",##OPENING SPLIT1
                                   col_top_lof,
                                   "'",
                                   ">",
                                   geneCellContent,##Cell Content
                                   "</td>",##CLOSING SPLIT1
                                   sep="",
                                   collapse = "")
              
            }else if(geneMechanism=="GOF"){
              
              table_content<-paste(table_content,
                                   "<td style='background-color:",##OPENING SPLIT 2
                                   col_top_gof,
                                   "'",
                                   ">",
                                   geneCellContent,##Cell Content
                                   "</td>",##CLOSING SPLIT 2
                                   "</tr></table>",##CLOSING SUBTABLE (THE MAIN COLUMN IS CLOSED AT ANOTHER SCRIPT POINT)
                                   sep="",
                                   collapse = "")
            }
          }
        }else if( sum(bothGeneScore >= minRequiredScore) == 0 ){
          ##No evidence for GOF or LOF
          ##Standard grey cells
          ##We add a grey cell with no info
          #"--"
          
          table_content<-paste(table_content,
                               "<td style='background-color:",
                               col_default,
                               "'",
                               ">",
                               "-",##Cell Content
                               sep="",
                               collapse = "")
        }
  }
      
      ##########################
      #Ending column in the row
      table_content<-paste(table_content,"</td>",
                           sep="",
                           collapse = "")
    }
    
    ########################
    ##Ending table row
    table_content<-paste(table_content,"</tr>",
                         sep="",
                         collapse = "")
  }
  
  
  
  #####################################
  ## Closing table && heatmapSection
  
  tail_html<-'</table></div><br>' 
  
  whole_html<-paste(header_html, 
                    table_content,
                    tail_html,
                    sep = "",
                    collapse = "")
 
  ###################################################
  ## Adding Table With Numeric Values Behind Heatmap
  ##################################################
  ##Required Function
  ##Keep comparing with table_html generation function 
  heatm_targetMatrix<-targetMatrix
  if("phaseFree" %in% colnames(heatm_targetMatrix)){
    colnames(heatm_targetMatrix)[colnames(heatm_targetMatrix)=="phaseFree"]<-"Cell Type Independent"  
  }
  
  if("affected_gene" %in% colnames(heatm_targetMatrix)){
    colnames(heatm_targetMatrix)[colnames(heatm_targetMatrix)=="affected_gene"]<-"Gene"  
  }
  heatm_targetMatrix$GeneImpact<-NULL
  
  numbersHeat_html<-generic_table_html_generation(targetMatrix = heatm_targetMatrix)
  
  
  ##Adding heatmap numbers section
  whole_html<-paste(whole_html,
                    "<div class='numbersHeatmap'>",
                    "<button type='button' class='collapsible_subsectionMainResults'>",
                    "<h3>",
                    "Pathogenic Scores",
                    "</h3>",
                    "</button>",
                    "<div class='content_subsectionMainResults'>",
                    numbersHeat_html,
                    "</div>",
                    "</div>",
                    sep="",
                    collapse="")
  
  javaScript_collapsible<-"<script>
    var coll = document.getElementsByClassName('collapsible_subsectionMainResults');
    var i;

    for (i = 0; i < coll.length; i++) {
    coll[i].addEventListener('click', function() {
    this.classList.toggle('active');
    var content = this.nextElementSibling;
    if (content.style.display === 'block') {
    content.style.display = 'none';
    } else {
    content.style.display = 'block';
    }
    });
    }
    
    </script>"
  
  ##Let's remove all line breaks, introduced from internet copypaste java & html code
  javaScript_collapsible<-gsub("[\r\n]", "", javaScript_collapsible)
  
  ##Adding JavaScript collapsible info
  whole_html<-paste(whole_html,
                    javaScript_collapsible,
                    sep="",
                    collapse="")
  
  
  ###################################################
  ## Adding Patient Considered Information Section
  ###################################################

  matrixToPlot<-patientResults$patientInfo
  ##Deprecated:Remove for the plot not relevant Columns
  ##matrixToPlot$patientID<-NULL ## For now providing all input info
  
  #Generating html table
  patientInfo_html<-generic_table_html_generation(targetMatrix = matrixToPlot)
  
  ##I'm going to apply same style as for heatmap numbers collapsible sections, so we are going to repeat the button classes
  ##And so on because we want them to look the same
  
  ##Adding heatmap numbers section
  whole_html<-paste(whole_html,
                    "<div class='numbersHeatmap'>",
                    "<button type='button' class='collapsible_subsectionMainResults'>",
                    "<h3>",
                    "Prediction Considered Information",
                    "</h3>",
                    "</button>",
                    "<div class='content_subsectionMainResults'>",
                    patientInfo_html,
                    "</div>",
                    "</div>",
                    sep="",
                    collapse="")
  
  
  ##########################################
  ## Creating gene-phase report section
  ##########################################
  reportSection<-"<div class='reportSection'>" ##style='border-style:solid; border-color:#e5e7e9;'

  ##################################################
  ##Function to help construct report
  source("functions/body_report_GeneSummary.R")
  
  
  ##for each gene above minScore create section with link to jump from table
  for(reportUnit in patientResults$genesConditions_ToReport){ ##relevant_Gene_Phase
    ##remember reportUnit is FUNDAMENTAL are the ids used to build JS functions
    
    ##div reportEntity, for single gene-phase specific
    ##Link will direct the user here
    ##By default the information is hidden
    startReportEntity<-paste0("<div id='",reportUnit, "'",
                              "style='display:none;'>")
    
    #Let's control with a button, expanding and collapsing
    titleReportUnit<-paste("<button type='button' class='collapsibleReportSections'",
                           "onclick='",
                           "myFunction_HideReport_",
                            reportUnit,
                           "()'",
                           "><h1><b>Report ",
                           reportUnit,
                           "</b></h1>",
                           "</button>",
                           sep="")
    
    ########################
    ##Report To append
    ########################
    targetGene<-unlist(strsplit(x=reportUnit, split="_", fixed = TRUE))[1]
    
    report_body<-body_report_GeneSummary(patientResults = patientResults,
                                          reportUnit = reportUnit,##combination target Gene + targetMechanism + target Phase
                                          minRequiredScore = minRequiredScore,
                                          targetGene = targetGene )  
    
    ##############
    ##Creating JS function and adding it to report body to minimize//hide results when considered
    ##based on click, "Hide Results", and based on click header of the section
    functionsJS<-paste(
      "function myFunction_HideReport_",
      reportUnit,
      "() {
              var x = document.getElementById('",
      reportUnit,
      "');
          x.style.display = 'none';
          }",
      sep="",
      collapse = "") 
    
    ##remove linebreaks
    functionsJS<-gsub("[\r\n]", "", functionsJS)
    
    ##Add to allJs functions, and append them after the table
    all_functionsJS<-paste(all_functionsJS,
                           functionsJS,
                           sep="",
                           collapse = "")
    
    ##Incorporate the function in an anchor to the bottom of the report entity
    ##And also at the beginning
    anchorSentence<-paste(
      "<a onclick='",
      "myFunction_HideReport_",
      reportUnit,
      "()'",
      " style='cursor:pointer;'> Hide results </a>",
      sep=""
    )
    
    
    ##Closing single report entity
    reportEntity<-paste(startReportEntity,
                         titleReportUnit,
                         report_body,
                         anchorSentence,
                         "</div>",##This div is for closing the report gene-phase-pathomech specific report
                         sep="",
                         collapse="")
    
    ##########################################################
    ##Appending report entity to whole ongoing report section
    reportSection<-paste(reportSection,
                         reportEntity,
                         sep="",
                         collapse="")
    
  }
  
  ################################
  ##Closing whole report section
  ################################
  
  reportSection<-paste(reportSection,
                       "</div>",
                       sep="",
                       collapse="")

  ######################################################
  ##Closing script section and adding them to the html
  ######################################################
  
  all_functionsJS<-paste(all_functionsJS,
                         "</script>",
                         sep="",
                         collapse = "")
  
  whole_html<-paste(whole_html,
                    all_functionsJS,
                    sep = "",
                    collapse = "")  
  
    
  ##################################################
  ##ASSEMBLING whole html including report section
  ##################################################
  
  whole_html<-paste(whole_html,
                    reportSection,
                    ##'<a href="#top">Jump to top of page</a><br>',
                    sep = "",
                    collapse = "")
  
  ##
  ##Put all main results inside the content of a metaDiv
  whole_html<-paste(
    '<div class="wrapperMainSingleResults"',
    whole_html,
    '</div>',
    '</body></html>',
    sep = "",
    collapse = ""
  )
  
  return(whole_html)
}












