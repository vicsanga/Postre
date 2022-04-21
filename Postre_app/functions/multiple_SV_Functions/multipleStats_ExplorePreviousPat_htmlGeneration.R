#################
##Function to create the table for the pathomechanisms
##Per phenotype
#################

##Required Function
source(file = "functions/multiple_SV_Functions/table_html_generation.R", local = TRUE)

multipleStats_htmlGeneration<-function(cohort_results, consideredPheno, ids_append,AllPatientsInfo,explPreviousPatSection){
  ##ids_append to avoid conflicts with tables generated in multiples submission option
  ##Html generation, I need a function to convert dataFrame to html table
  ##Y poco mas meter los entries y arreglado
  ##explPreviousPatSection (true, false) to distinguish if html generated for Expl Previous Patient section, of for
  ##User multiple patient submission
  
  ##It is the header info to create the html collapsibleExplorePrev list
  # https://www.w3schools.com/howto/howto_js_collapsibleExplorePrev.asp
  headerPart<-paste("<html>
  <head>
  <meta name='viewport' content='width=device-width, initial-scale=1'>
  <style>
  /*Style information in MainInterfaceStyling html file*/
  </style>
  </head>
  <body>
  ",
  ##Adding wrapper of the whole document to be able to modify it and clean results upon resubmission when performing Multiple SV analysis
  if(explPreviousPatSection == FALSE){
    ##So, multipleSV submission from user
    "<div class= 'wrapperMultipleSVSubmission'>"
  }else{
    #So, explore previous SV html generation
    "<div class= 'wrapperExplorePreviousPat'>"
  },
  sep="",
  collapse = "")
  
  ##gsub
  ##Let's remove all line breaks, introduced from internet copypaste java & html code
  headerPart<-gsub("[\r\n]", "", headerPart)
  
  
  master_pathomechTable_html<-""
  
  master_pathomechTable_html<-paste(master_pathomechTable_html,
                                    headerPart)
  
  allPhenoInfoHtml<-""##To append info from all phenotypes

  ##Adding sections in the HTML, one per phenotype + the How to navigate results section + Submitting prediction section
  for(targetPheno in c("HowToNavigateResults",sort(consideredPheno),"SubmittingPrediction")){
    print(targetPheno)
    
    phenoResults<-cohort_results[[targetPheno]]
    
    ###We need to do each of the following sections collapsible
    
    ##To plot nice header
    if(targetPheno=="head_neck"){
      phenoTag<-"Head & Neck SVs"
    }else if(targetPheno=="cardiovascular"){
      phenoTag<-"Cardiovascular SVs"
    }else if(targetPheno=="limbs"){
      phenoTag<-"Limbs SVs"
    }else if(targetPheno=="neurodevelopmental"){
      phenoTag<-"Neurodevelopmental SVs"
    }else if(targetPheno=="nervous_system"){
      phenoTag<-"Nervous System - Brain"
    }else if(targetPheno=="HowToNavigateResults"){
      phenoTag<-"How to navigate this section?"
    }else if(targetPheno=="SubmittingPrediction"){
      phenoTag<-"Submit SV for prediction"
    }
    
    # https://www.w3schools.com/howto/howto_js_collapsible.asp
    pheno_html<-paste("<div class='phenoRecurrencyStats'>",
                      "<button type='button' class=",
                      if(explPreviousPatSection==TRUE){
                        "'collapsibleExplorePrev_mainSection'"
                      }else{
                        "'collapsibleMainSection'"
                      },
                      ##Assignar aqui el color del boton para que canvie segun si es phenotype data o otro tipo de dato
                      
                      if(targetPheno == "HowToNavigateResults"|| targetPheno == "SubmittingPrediction"){
                        "style='background-color:#1D3354;'"
                      }else{
                        ##So aggregated information of SVs for different phenotypes
                        "style='background-color:#467599;'"
                      }
                      ,
                      ">",
                      "<h2 style='color:#FFFFFF;font-size: 40px; font-weight: bold' >",
                      phenoTag,
                      "</h2>",
                      "</button>",
                      sep = "",
                      collapse = "")
    
    
    ##Adding start of the content section 
    pheno_html<-paste(pheno_html,
                      "<div class='content_primary'>",
                      sep = "",
                      collapse = "")
    
    
    if(targetPheno == "HowToNavigateResults"){
      
      ################## 
      ## Here provide info on how to navigate the patient cohort results
      #####################
      
      ##For the youtube video embedding
      ##Important to get the "embed" option link, is the one that works
      ##https://stackoverflow.com/questions/19826663/why-is-my-youtube-video-not-showing-up/19826763
      ##Just copy pasting the youtube link does not work properly, instead click share button, option embed, and youtube provides the link to add to the html
      ##youtube size de partida width="560" height="315" . Reajustar en base a ello
      ##http://thenewcode.com/717/Force-Embedded-YouTube-Videos-To-Play-In-HD  ##Not everything done but the idea
      contentExplorePreviousPat<-paste(c("<div class = 'explanationExplorePreviousPat'>",
                                         "<p style = 'font-size: 20px;'> In the video below there is an explanation about how to navigate the current page. Reproduce it in <b>Full Screen</b> and <b>High Quality</b>(1080p) for optimal visualitzation. </p>",
                                         "<div class ='videotutorialExplPreviousPat'>",##div used to center video
                                         '<p align="center"><iframe width="640" height="360" src="https://www.youtube.com/embed/movwXisGsmM" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></p>',
                                         "</div>",
                                         "</div>"),
                                       sep=" ",
                                       collapse=" ")
      
      pheno_html<-paste(pheno_html,
                        contentExplorePreviousPat,
                        sep="")
      
      
    }else if(targetPheno == "SubmittingPrediction" ){
      
      ##########################################################
      ##Information for patient submission with simple button
      ##########################################################
      
      ##Aqui va a haber un boton para el que va a cambiar el id segon si es true o false el parametro de previous pat
      ##explPreviousPatSection = TRUE
      ##Y para el pheno of INterest too
      
      htmlSubmissionSection_part1<-paste(
      '
      <div class="phenoID well">
      <div class="form-group shiny-input-container">
      <label class="control-label" id="Pheno of Interest" for="Pheno of Interest">Select phenotype</label>
      <div>
<select id=',
        if(explPreviousPatSection == TRUE){
          '"aggregatedResults_phenoId_ExplorePreviousPat">' 
        }else{
          '"aggregatedResults_phenoId_MultipleSVSubmission">'
        },
        sep = "",
        collapse = ""
      )
      
    ##Building part 2
      htmlSubmissionSection_part2<-""
        ##Metemos las options
        for(pheno in consideredPheno){
          #print(pheno)
          ##Mirar a ver si mete el string?? igual que el ifelse de arriba, probar
          ##Poner niceEquivalent, si phenoTal, nice string = tal con una tabla o algo
          htmlSubmissionSection_part2<-paste(htmlSubmissionSection_part2,
                '<option value="',pheno,'">',pheno,'</option>',
                sep="",
                collapse="")
        }
      
      ##Building part 3
      htmlSubmissionSection_part3<-        '</select>
    </div>
    </div>
    </div>
    '
    
      ##Part 4 with introuce patient ID box
      htmlSubmissionSection_part4<- paste('
      <div class="svID well">
      <div class="form-group shiny-input-container">
  <label class="control-label" id="sv_id-label" for="sv_id-label">Introduce Structural Variant ID</label>
  <input id=',
    if(explPreviousPatSection == TRUE){
      '"aggregatedResults_svID_ExplorePreviousPat"' 
    }else{
      '"aggregatedResults_svID_MultipleSVSubmission"'
    },
  'type="text" class="form-control" placeholder="introduce SV ID"/>

</div>
</div>',
       sep = "",
       collapse = ""
      )
      
      ##Choosing Running Mode
      htmlSubmissionSection_part5<-paste(
'
<div class="runMode_aggregatedRes_SubmBox well">
    <div class="form-group shiny-input-container">
      <label class="control-label" id="runMode_Aggregated-label" for="runMode_Aggregated">Running mode</label>
      <div>
        <select id=',
        if(explPreviousPatSection == TRUE){
          '"runMode_AggregatedRes_ExplorePreviousPat"'                            
        }else{
          '"runMode_AggregatedRes_MultipleSVSubmission"'                             
        },
        '><option value="Standard" selected>Standard</option>
<option value="High-Specificity">High-Specificity</option></select>
    </div>
  </div>
</div>',
        sep = "",
        collapse = ""
      )
      
      ##Falta submission button
      ##Ajustar con el respectivo ID de si es del ExplPrev o del SV Multiple Submission
      
      htmlSubmissionSection_part6<-paste('
          <div class="submissionAggregated">
          <button id=',
          if(explPreviousPatSection == TRUE){
            '"click_AggregatedRes_ExplorePreviousPat"'                            
          }else{
            '"click_AggregatedRes_MultipleSVSubmission"'                             
          },                                         
          ' type="button" class="btn btn-default action-button" style="color: #fff; background-color: #1D3354; border-color: #467599;">
          <i class="fa fa-paper-plane" role="presentation" aria-label="paper-plane icon"></i>
          Submit
        </button>
        </div>'
      )
      
      
      htmlSubmissionSection_Title<-'<div class="titleSubmissionBoxAggregatedResults">
        <h3><center><b style="color:#467599;">Run prediction with the Structural Variant of interest</b></center></h3>
        </div>'
      
    ##Assembling html + ADDING TITLE FOR ABOVE OF THE SUBMISSION BOX
      htmlSubmissionSection_Full<-paste(htmlSubmissionSection_Title,
                                        '<div class=panelSubmAggregatedRes>',
                                        htmlSubmissionSection_part1,
                                        htmlSubmissionSection_part2,
                                        htmlSubmissionSection_part3,
                                        htmlSubmissionSection_part4,
                                        htmlSubmissionSection_part5,
                                        htmlSubmissionSection_part6,
                                        '</div>
                                        <br><br>', ##This <br> help give space to the bottom submission box shadowing, on the contrary does not appear
                                        sep="",
                                        collapse="")
      
      pheno_html<-paste(pheno_html,
                        htmlSubmissionSection_Full,
                        sep="")
      
      
      
      
    }else{
      
      ################################################################
      ##Hence here we are working with the Phenotypes Cohort results
      ################################################################
      
      targetNamesMatrixes<-names(phenoResults)
      #Excluding noGeneFoundInfo, we are going to do nothing with it regarding the html, and if we leave it it triggers an error
      targetNamesMatrixes<-targetNamesMatrixes[targetNamesMatrixes!="noGeneFoundInfo"]
      
      for(name_targetMatrix in targetNamesMatrixes){
        
        ##Customize collapsible names to make them nicer
        h1_tag_name_targetMatrix<-name_targetMatrix
        
        if(h1_tag_name_targetMatrix == "anyMechanism"){
          h1_tag_name_targetMatrix<-"Overview | Any Pathological Mechanism"
          
        }else if(h1_tag_name_targetMatrix == "DirectEffectLOF"){
          h1_tag_name_targetMatrix<-"Coding Effect | Loss of Function"
          
        }else if(h1_tag_name_targetMatrix == "DirectEffectGOF"){
          h1_tag_name_targetMatrix<-"Coding Effect | Gain of Function"
          
        }else if(h1_tag_name_targetMatrix == "LongRangeLOF"){
          h1_tag_name_targetMatrix<-"Long-Range Effect | Loss of Function"
          
        }else if(h1_tag_name_targetMatrix == "LongRangeGOF"){
          h1_tag_name_targetMatrix<-"Long-Range Effect | Gain of Function"
          
        }else if(h1_tag_name_targetMatrix == "errorInfo"){
          h1_tag_name_targetMatrix<-"Unresolved SVs"
          
        }else if(h1_tag_name_targetMatrix == "patientsInfo"){
          h1_tag_name_targetMatrix<-"SVs Information"
          
        }
        
        ##Starting collapsible section
        pheno_html<-paste(pheno_html,"<div class='",
                          name_targetMatrix,
                          "'>",
                          "<button type='button' class='",
                          if(explPreviousPatSection == TRUE){
                            "collapsibleExplorePrev"
                          }else{
                            "collapsible"
                          },
                          "'>",
                          "<h1>",
                          h1_tag_name_targetMatrix,
                          "</h1>",
                          "</button>",
                          sep = "",
                          collapse = "")
        
        
        if((name_targetMatrix != "errorInfo") && (name_targetMatrix != "patientsInfo")){
          #source(file = "functions/multiple_SV_Functions/table_html_generation.R")
          ##So generating table with recurrency pathological mechanisms per Target Pathomechanism
          targetMatrix<-phenoResults[[name_targetMatrix]]
          
          ##If targetMatrix has 0 rows, do not add it
          ##Instead, add, no Info for this section
          
          if(nrow(targetMatrix)>0){
            ##So there is information in the matrix, hence add it
            
            if("phaseFree" %in% colnames(targetMatrix)){
              colnames(targetMatrix)[colnames(targetMatrix)=="phaseFree"]<-"Cell Type Independent"  
            }
            
            ##Renaming "patient" by "SV_ID" and "patients" by "SV_IDs"
            ##Renaming columns now that we want to represent their information more clearly
            if("patients" %in% colnames(targetMatrix)){
              colnames(targetMatrix)[colnames(targetMatrix)=="patients"]<-"SV_IDs"  
            }
            
            if("Num_Patients" %in% colnames(targetMatrix)){
              colnames(targetMatrix)[colnames(targetMatrix)=="Num_Patients"]<-"Num_SVs"  
            }
            
            table_html<-table_html_generation(targetMatrix = targetMatrix ,
                                              name_targetMatrix = name_targetMatrix,
                                              ids_append = ids_append,
                                              targetPheno= targetPheno)
            
            pheno_html<-paste(pheno_html,
                              "<div class='content'>",
                              table_html,
                              "</div>",
                              sep="") 
          }else{
            #####################################
            ##Hence there is no info to display
            #####################################
            pheno_html<-paste(pheno_html,
                              "<div class='content'>",
                              "<p>No data to display</p>",
                              "</div>",
                              sep="") 
          }
          
        }else if(name_targetMatrix == "patientsInfo"){
          
          targetMatrix<-phenoResults[[name_targetMatrix]]
          
          ##Rename patientID column by SV_ID since 1 patient can have >1 SV
          if("patientID" %in% colnames(targetMatrix)){
            colnames(targetMatrix)[colnames(targetMatrix)=="patientID"]<-"SV_ID"  
          }
          
          ##source(file = "functions/multiple_SV_Functions/table_html_generation.R")
          ##source(file = "functions/multiple_SV_Functions/patientTable_html_generation.R")
          table_html<-table_html_generation(targetMatrix = targetMatrix,
                                            name_targetMatrix=name_targetMatrix,
                                            ids_append = ids_append,
                                            targetPheno= targetPheno)
          
          pheno_html<-paste(pheno_html,
                            "<div class='content'>",
                            table_html,
                            "</div>",
                            sep="")
          
        }else if(name_targetMatrix == "errorInfo"){
          ##so in the error info section
          if(length(phenoResults[["errorInfo"]])>0){
            patientsError<-paste(c("<div class='content'>",
                                   "<p>We are sorry to inform you that a problem occurred while predicting the impact of the structural variants with the identifiers: ",
                                   paste0(as.character(unlist(phenoResults[["errorInfo"]])),
                                          collapse=", "),
                                   "due to biological limitations or technical issues related with the structural variants information and affected loci. Please check if the data introduced for them is correct and visit the User Guide if you have any doubt.
       If the problem persists contact us through sv_radalab@gmail.com to provide you more information about the particularities of these genetic rearrangements.",
                                   "</p>",
                                   "</div>"),
                                 sep=" ",
                                 collapse=" ") 
            
          }else{
            ##No problems found, indicate
            patientsError<-paste(c("<div class='content'>",
                                   "<p>No problems encountered when processing the structural variants.</p>",
                                   "</div>"),
                                 sep=" ",
                                 collapse=" ") 
          }
          
          pheno_html<-paste(pheno_html,
                            patientsError,
                            sep="")
          
          
        }
        ###########################################
        ##Closing Table//Pathomechanism Div Section
        pheno_html<-paste(pheno_html,
                          "</div>",
                          sep="")
      }
    }

    
    #########################
    ##Closing Main Sections
    #########################
    pheno_html<-paste(pheno_html,
                      "</div>",##TO close the div wrapper of the content
                      "</div>",##To close the div from class phenoRecurrencyStats
                      sep="")
    
    
    ################
    ##Append phenoinfo to all pheno Info
    allPhenoInfoHtml<-paste(allPhenoInfoHtml,
                            pheno_html,
                            sep="")
    
  }
  
  ######################
  ## Wrapping Up html
  ######################
  
  
  #####################################
  ##Attach JS collapsible functionality
  #####################################
  #JS collapsibleExplorePrev
  endPart<-
  
  ##Adding functions
  if(explPreviousPatSection == TRUE){
    "<script>
    var coll = document.getElementsByClassName('collapsibleExplorePrev');
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
    
    var coll = document.getElementsByClassName('collapsibleExplorePrev_mainSection');
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
    
  }else{
    "<script> 
    var coll = document.getElementsByClassName('collapsible');
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
    
    var coll = document.getElementsByClassName('collapsibleMainSection');
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
    
  }
  
  ##Let's remove all line breaks, introduced from internet copypaste java & html code
  endPart<-gsub("[\r\n]", "", endPart)
  
  master_pathomechTable_html<-paste(master_pathomechTable_html,
                                    allPhenoInfoHtml,
                                    endPart,
                                    "</div>",#closing Document main wrapper (wrapperMultipleSVSubmission | wrapperExplorePreviousPat)
                                    "</body>",
                                    "</html>",
                                    sep="")
  
  
  ##cerrar body & html
  
  return(master_pathomechTable_html)
  
}
