#################################################
### Rada-Iglesias Lab SV Effect Prediction Tool
#################################################

##To consider bioconductor repos
##To avoid this error: https://community.rstudio.com/t/deployment-error-unable-to-determine-the-location-for-some-packages/102312
# options(repos = BiocManager::repositories())
#options('repos')

library(shiny)
library(waiter)
library(shinybusy)
library(shinyjs)
library(shinyWidgets)
library(DT)
library(plotrix)##For enhancers, ellipse shape representation
library(shape)##For curve arrows representation
library(diagram)##For curve arrows representation

##Librerias para el liftover tratar de poner aqui
##Si se cargan dentro del source de la funcion tarda la vida
# library(liftOver)
# library(GenomicRanges)
# library(rtracklayer)
# library(reticulate)



###Setwd in the folder where all the app info is hosted
# setwd("~/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/Postre/Postre_app")

####################################
###Let's load required Functions
####################################

source(file = "functions/MasterWrapperSinglePrediction.R",local = TRUE)##many other functions loaded here
source(file = "functions/error_Report.R",local = TRUE)
source(file = "functions/noGenesFound_report.R",local = TRUE)

##############################################
##Functions for Multiple Patients Submission
##############################################
source("functions/multiple_SV_Functions/cohortResults_Parser.R", local = TRUE)
source("functions/multiple_SV_Functions/multipleStats_ExplorePreviousPat_htmlGeneration.R",local=TRUE)

##To deal with rounding .5 problems
##round2 function
source(file = "functions/roundingHalfs.R", local = FALSE) ##Local == FALSE to be loaded in the global env so that all functions can find it

####################################################################################################
##Patients tables for User Easy Submission. Maybe only loading when Submission of this kind done
##Also if database huge it may crash the app, so maybe for the future, load only specific patients
##Or add metadata to the value line so no need to load anything, but... in this case html will increase accordingly
##Maybe having 1 file per patient and reading specific file. That will remove data loaded to memory.
##But for now, managable with loading data frames
####################################################################################################
load("data/patientTables/usrGuide_SingleSubmTable.RData")

##Loading allPatients Used for the explorePreviousPatient section
load("data/patientTables/AllPatientsInfo.RData")

##Scores For Considering Relevant result
##minRequired scores for heatmap coloring and geneReport generation
minScore<-0.8
highScore<-0.9


###Meter como variables, todos los valores que  estan estaticos en el menu
relevantChr<-c(paste("chr",1:22,sep = ""), "chrX")##chrY excluded not all data available for chrY


#To avoid navbar collapse in smaller screens
# https://stackoverflow.com/questions/21738417/bootstrap-remove-responsive-from-navbar
navbar_js<-"@media (max-width: 768px) {
    .navbar-header {
        float: left;
    }

    .navbar {
        border-radius: 4px;
        min-width: 400px;
    }

    .nav-tabs-justified > li > a {
        border-bottom: 1px solid #ddd;
        border-radius: 4px 4px 0 0;
    }
    .nav-tabs-justified > .active > a,
    .nav-tabs-justified > .active > a:hover,
    .nav-tabs-justified > .active > a:focus {
        border-bottom-color: #fff;
    }

    .nav-justified > li {
        display: table-cell;
        width: 1%;
    }
    .nav-justified > li > a {
        margin-bottom: 0;
    }

    .nav-tabs.nav-justified > li > a {
        border-bottom: 1px solid #ddd;
        border-radius: 4px 4px 0 0;
    }
    .nav-tabs.nav-justified > .active > a,
    .nav-tabs.nav-justified > .active > a:hover,
    .nav-tabs.nav-justified > .active > a:focus {
        border-bottom-color: #fff;
    }

    .nav-tabs.nav-justified > li {
        display: table-cell;
        width: 1%;
    }
    .nav-tabs.nav-justified > li > a {
        margin-bottom: 0;
    }

    .navbar-right .dropdown-menu {
        right: 0;
        left: auto;
    }
    .navbar-right .dropdown-menu-left {
        right: auto;
        left: 0;
    }
    .container {
        min-width: 400px;
    }

    .navbar-collapse {
        width: auto;
        border-top: 0;
        box-shadow: none;
    }
    .navbar-collapse.collapse {
        display: block !important;
        height: auto !important;
        padding-bottom: 0;
        overflow: visible !important;
    }
    .navbar-collapse.in {
        overflow-y: visible;
    }
    .navbar-fixed-top .navbar-collapse,
    .navbar-static-top .navbar-collapse,
    .navbar-fixed-bottom .navbar-collapse {
        padding-right: 0;
        padding-left: 0;
    }

    .container > .navbar-header,
    .container-fluid > .navbar-header,
    .container > .navbar-collapse,
    .container-fluid > .navbar-collapse {
        margin-right: 0;
        margin-left: 0;
    }

    .navbar-static-top {
        border-radius: 0;
    }

    .navbar-fixed-top,
    .navbar-fixed-bottom {
        border-radius: 0;
    }

    .navbar-toggle {
        display: none;
    }

    .navbar-nav {
        float: left;
        margin: 0;
    }
    .navbar-nav > li {
        float: left;
    }
    .navbar-nav > li > a {
        padding-top: 15px;
        padding-bottom: 15px;
    }
    .navbar-nav.navbar-right:last-child {
        margin-right: -15px;
    }

    .navbar-left {
        float: left !important;
    }
    .navbar-right {
        float: right !important;
    }

    .navbar-form .form-group {
        display: inline-block;
        margin-bottom: 0;
        vertical-align: middle;
    }
    .navbar-form .form-control {
        display: inline-block;
        width: auto;
        vertical-align: middle;
    }
    .navbar-form .control-label {
        margin-bottom: 0;
        vertical-align: middle;
    }
    .navbar-form .radio,
    .navbar-form .checkbox {
        display: inline-block;
        padding-left: 0;
        margin-top: 0;
        margin-bottom: 0;
        vertical-align: middle;
    }
    .navbar-form .radio input[type='radio'],
    .navbar-form .checkbox input[type='checkbox'] {
        float: none;
        margin-left: 0;
    }
    .navbar-form .has-feedback .form-control-feedback {
        top: 0;
    }

    .navbar-form {
        width: auto;
        padding-top: 0;
        padding-bottom: 0;
        margin-right: 0;
        margin-left: 0;
        border: 0;
        -webkit-box-shadow: none;
                box-shadow: none;
    }
    .navbar-form.navbar-right:last-child {
        margin-right: -15px;
    }

    .navbar-text {
        float: left;
        margin-right: 15px;
        margin-left: 15px;
    }
    .navbar-text.navbar-right:last-child {
        margin-right: 0;
    } 
}"

ui <-function(req){
  
  return(div(
    class="container",
    ##Mirar esto del title que no me acaba lo de meterlo como si no existiera nada
    div(class="titleBrowser",
        titlePanel(title="POSTRE: Prediction Of STRuctural variant Effects")
    ),
    includeHTML("html_scripts/MainInterfaceStyling.html"),
    
    ######################
    ## Here go the plots
    ######################
    div(class="mainPanel",
        
        navbarPage(title = "POSTRE",
                   id = "inTabset",
                   tabPanel(title = "Single SV Submission",
                            value = "SingleSubmissionTab",
                            icon=icon("dna"),
                            div(class="patientInputPanel",
                                includeHTML("html_scripts/Single_SV_Submission_Interface.html"),
                                ######################
                                ## Patient Input Panel
                                ######################
                                div(class="sideBarClass",
                                    ##If I put here on the html, h2 tag it gets styled as the header, so style features 
                                    ##are connected
                                    div(class="formAndTitle_patientInput",
                                        ##Patient SV type
                                        HTML('<h2 id="titleSideBarPannel">Introduce Structural Variant and Phenotype</h2>'),
                                        div(class="patientFormulary",
                                            
                                            div(class="inp1",
                                                wellPanel(selectInput(inputId = "typeSV", label = "Type SV",
                                                                      choices = c("Deletion","Duplication","Inversion","Translocation"),
                                                                      selected = "Inversion"))
                                            ),
                                            div(class="inp2",
                                                ##Patient Phenotype
                                                wellPanel(selectInput(inputId = "phenoPatient",
                                                                      label = "Phenotype",
                                                                      ##Choices full name is then matched & renamed in GenomicData_Loader.R
                                                                      ##Alphabet order
                                                                      choices = c("Cardiovascular",
                                                                                  "Head & Neck",
                                                                                  "Limbs",
                                                                                  "Neurodevelopmental"
                                                                      ),
                                                                      selected ="Head & Neck" ))
                                            ),
                                            ##SV coordinates
                                            div(class="inp3",
                                                ###Breakpoint 1 
                                                wellPanel(
                                                  selectInput(inputId = "bp1_chr",
                                                              label = "Breakpoint 1 chromosome",
                                                              choices = relevantChr,
                                                              selected="chr6"),
                                                  
                                                  textInput(inputId = "bp1_coord",
                                                            label = "Breakpoint 1 coordinates (coord1,coord2 if interval)",
                                                            value = "10355280")
                                                )),
                                            
                                            div(class="inp4",
                                                ###Breakpoint 2
                                                wellPanel(
                                                  selectInput(inputId = "bp2_chr",
                                                              label = "Breakpoint 2 chromosome",
                                                              choices = relevantChr,
                                                              selected = "chr6"),
                                                  
                                                  textInput(inputId = "bp2_coord",
                                                            label = "Breakpoint 2 coordinates (coord1,coord2 if interval)",
                                                            value = "99103873")
                                                )),
                                            ##Specifying reference Genome
                                            div(class="inp5",
                                                # wellPanel(selectInput(inputId = "refGenome", label = "Reference Genome",
                                                #                       choices = c("hg19","hg38"),
                                                #                       selected = "hg19"))
                                                
                                                wellPanel(
                                                  div(class="textRefGenome",
                                                      HTML('<h4><b>NOTE: Reference Genome Coordinates required in GRCh37/hg19</b></h4>'),
                                                      HTML('<p style="font-size:15px;">Please, visit any of the following websites to convert your coordinates to hg19 if you need it:
                                                         <a href= "http://genome.ucsc.edu/cgi-bin/hgLiftOver" target="_blank">UCSC LiftOver</a>, 
                                                         <a href= "https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core" target="_blank">ENSEMBL Assembly Converter</a>, 
                                                         <a href= "https://liftover.broadinstitute.org/" target="_blank">Broad Institute LiftOver</a>
                                                         ')
                                                  )
                                                ),
                                            ),
                                            ##Selecting Running Mode
                                            div(class="inp6",
                                                wellPanel(selectInput(inputId = "runMode_single", label = "Running mode",
                                                                      choices = c("Standard","High-Specificity"),
                                                                      selected = "Standard"))
                                            ),
                                            div(class="submission",
                                                #Send request
                                                actionButton(inputId = "click",
                                                             label = "Submit",
                                                             icon("paper-plane"),
                                                             style="color: #fff; background-color: #1D3354; border-color: #467599;")
                                            )
                                        )
                                    )
                                )
                            )
                   ),
                   
                   
                   
                   tabPanel(title = "Single SV Results",
                            value = "overview",
                            icon=icon("poll"),
                            div(class="heatmapPage",
                                tableOutput(outputId = "masterSummary_result")
                            ),
                            
                            ##TO AVOID OVERWRITTING THE PLOT VARIABLE:
                            ##I need to add another outputTable for the results derived of the UsrGuide_SingleSubmission section
                            ##If don't, and I try to put all results in masterSummary_result table, the tool does not work properly, 
                            ##Because one of the two reactive events stop working
                            
                            div(class="heatmapPage",
                                tableOutput(outputId = "masterSummary_result_usrGuide_SingleSubmission")
                            ),
                            
                            div(class="heatmapPage",
                                tableOutput(outputId = "masterSummary_result_explorePreviousPatient_SingleSubmission")
                            ),
                            div(class="heatmapPage",
                                tableOutput(outputId = "masterSummary_result_MultipleSVSubmission_SingleSubmission")
                            )
                   ),
                   tabPanel(title = "Multiple SV Submission",
                            value = "multiple_Submission",
                            icon=icon("dna"),
                            div(class="patientInputPanel",
                                includeHTML("html_scripts/Multiple_SV_Submission_Interface.html"),
                                ######################
                                ## Multiple SV Input Panel
                                ######################
                                div(class="sideBarClassMultiple",
                                    ##If I put here on the html, h2 tag it gets styled as the header, so style features 
                                    ##are connected
                                    
                                    #START MULTIPLE SUBMISSION
                                    div(class="formAndTitle_patient_Multiple_Input",
                                        ##Patient SV type
                                        HTML('<h2 id="titleSideBarPannel">Upload Multiple Structural Variants with Phenotypes</h2>'),
                                        div(class="patientMultipleFormulary",
                                            
                                            div(class="multipleSubmission",
                                                wellPanel(fileInput(
                                                  inputId = "multipleFileInfo",
                                                  label="SelectFile",
                                                  multiple = FALSE,
                                                  accept = ".tsv"
                                                ),style = "padding: 5px; padding-bottom:0px;")
                                            ),
                                            
                                            ##Selecting Running Mode
                                            div(
                                              wellPanel(selectInput(inputId = "runMode_Multiple", label = "Running mode",
                                                                    choices = c("Standard","High-Specificity"),
                                                                    selected = "Standard"))
                                            ),
                                            div(class="submission",
                                                id="multipleSubm",
                                                ##Send request
                                                actionButton(inputId = "clickMultiple",
                                                             label = "Submit",
                                                             icon("paper-plane"),
                                                             style="color: #fff; background-color: #1D3354; border-color: #467599;")
                                            )
                                        )
                                    )
                                    
                                    #######END MULTIP SUBMISSION
                                    
                                )
                            )
                   ),
                   
                   tabPanel(title = "Multiple SV Results",
                            value = "multiple_overview",
                            icon=icon("poll"),
                            ##Same format that for Explore Previous Patients Section
                            div(class="previousPatientsPage",
                                tableOutput(outputId = "master_multipleStats")
                            )
                   ),
                   
                   tabPanel(title = "Explore Previous SVs",
                            icon=icon("search"),
                            div(class="previousPatientsPage",
                                includeHTML("html_scripts/ExplorePreviousPatients.html")
                            )
                   ),
                   tabPanel(title = "User Guide",
                            value = "userGuide",
                            icon=icon("info-circle"),
                            div(class="userGuidePage",
                                includeHTML("html_scripts/UserGuide_page.html")
                            )
                   )
                   ,
                   tags$head(tags$style(HTML(navbar_js)))
        )
    ), 
    #####################
    ##Creating App Footer
    #####################
    ##On the included script, adjust also CSS elements styling
    includeHTML("html_scripts/FooterWebPage.html"),
    ## Load some javascript shiny functionalities
    useShinyjs()
  ))
}

server <- function(input, output, session){
  
  ####################################################
  ##Defining initial behaviour when clicking buttons
  ####################################################
  
  observeEvent(eventExpr = input$click, {
    
    updateTabsetPanel(session, 
                      inputId = "inTabset",
                      selected = "overview")
    
    runjs('
          document.getElementById("top").scrollIntoView();
          ')
    
    # ##Clear whtvr content can be present at the tab
    #Hence we are gonna hidde the heatmap section
    runjs(
      ##"document.getElementsByClassName('wrapperMainSingleResults')[0].style.visibility='hidden';"##With this we would just modify one instance of the class
      ##With the for, all elements belonging to the class change their status. With visibility hidden the space of the div is not erased
      "for (let el of document.querySelectorAll('.wrapperMainSingleResults')) el.style.display = 'none';"
    )
    
    runjs("window.scrollTo(0, 0)")
    
  })
  
  
  observeEvent(eventExpr = input$click_SingleSubmissionUserGuide, {
    
    updateTabsetPanel(session, 
                      inputId = "inTabset",
                      selected = "overview")
    
    runjs('
          document.getElementById("top").scrollIntoView();
          ')
    
    # ##Clear whtvr content can be present at the tab
    #Hence we are gonna hidde the heatmap section
    runjs(
      ##"document.getElementsByClassName('wrapperMainSingleResults')[0].style.visibility='hidden';"##With this we would just modify one instance of the class
      ##With the for, all elements belonging to the class change their status. With visibility hidden the space of the div is not erased
      "for (let el of document.querySelectorAll('.wrapperMainSingleResults')) el.style.display = 'none';"
    )
    
  })
  
  observeEvent(eventExpr = input$clickMultiple, {
    updateTabsetPanel(session, 
                      inputId = "inTabset",
                      selected = "multiple_overview")
    
    runjs('
          document.getElementById("top").scrollIntoView();
          ')
    
    ###Go to the top of the page
    runjs("window.scrollTo(0, 0)")
    
    runjs(
      ##"document.getElementsByClassName('wrapperMainSingleResults')[0].style.visibility='hidden';"##With this we would just modify one instance of the class
      ##With the for, all elements belonging to the class change their status. With visibility hidden the space of the div is not erased
      "for (let el of document.querySelectorAll('.patientInputPanel')) el.style.display = 'none';"
    )
  })
  
  observeEvent(eventExpr = input$click_AggregatedRes_ExplorePreviousPat, {
    
    updateTabsetPanel(session, 
                      inputId = "inTabset",
                      selected = "overview")
    
    runjs('
          document.getElementById("top").scrollIntoView();
          ')
    
    # ##Clear whtvr content can be present at the tab
    #Hence we are gonna hidde the heatmap section
    runjs(
      ##"document.getElementsByClassName('wrapperMainSingleResults')[0].style.visibility='hidden';"##With this we would just modify one instance of the class
      ##With the for, all elements belonging to the class change their status. With visibility hidden the space of the div is not erased
      "for (let el of document.querySelectorAll('.wrapperMainSingleResults')) el.style.display = 'none';"
    )
    
  })
  
  observeEvent(eventExpr = input$click_AggregatedRes_MultipleSVSubmission, {
    
    updateTabsetPanel(session, 
                      inputId = "inTabset",
                      selected = "overview")
    
    runjs('
          document.getElementById("top").scrollIntoView();
          ')
    
    # ##Clear whtvr content can be present at the tab
    #Hence we are gonna hidde the heatmap section
    runjs(
      ##"document.getElementsByClassName('wrapperMainSingleResults')[0].style.visibility='hidden';"##With this we would just modify one instance of the class
      ##With the for, all elements belonging to the class change their status. With visibility hidden the space of the div is not erased
      "for (let el of document.querySelectorAll('.wrapperMainSingleResults')) el.style.display = 'none';"
    )
    
  })
  
  ############################################################
  ####To change from footer anchor or button click to userguide
  observeEvent(input$do, {
    updateTabsetPanel(session, 
                      inputId = "inTabset",
                      selected = "userGuide")
    runjs("window.scrollTo(0, 0)")
  })
  
  observeEvent(input$done, {
    updateTabsetPanel(session, 
                      inputId = "inTabset",
                      selected = "userGuide")
    runjs("window.scrollTo(0, 0)")
  })
  
  observeEvent(input$doneMultiple, {
    updateTabsetPanel(session, 
                      inputId = "inTabset",
                      selected = "userGuide")
    runjs("window.scrollTo(0, 0)")
  })
  
  
  ##To go to the singlePage from UserGuide
  observeEvent(input$goSinglePage, {
    updateTabsetPanel(session, 
                      inputId = "inTabset",
                      selected = "SingleSubmissionTab")
    runjs("window.scrollTo(0, 0)")
  })
  
  #######################################################
  ## Computing prediction for Single Condition Submission
  #######################################################
  
  patientResults<-eventReactive(input$click, {
    ###Adding-Initializing waiter
    show_modal_spinner(text = HTML("Running Prediction <br>
                                   Please be patient"),
                       color="#0066ff")
    
    ##Creating Input Data Frame, to carry on prediction
    patientDataCols<-c("chr_Break1","coord_Break1","chr_Break2",
                       "coord_Break2","TypeSV",
                       "Phenotype","refGenome")
    patientData<-as.data.frame(matrix(data = 0, ncol = length(patientDataCols), nrow = 1))
    colnames(patientData)<-patientDataCols
    
    ###Filling Data with input info
    ###Ensure that everything is in character format (coord if considered numeric in initial step, then error jumps, due to strsplit)
    patientData$chr_Break1<-as.character(input$bp1_chr)
    patientData$coord_Break1<-as.character(input$bp1_coord)
    
    patientData$chr_Break2<-as.character(input$bp2_chr)
    patientData$coord_Break2<-as.character(input$bp2_coord)
    
    patientData$Phenotype<-as.character(paste(input$phenoPatient, collapse = ","))
    patientData$TypeSV<-as.character(input$typeSV)
    
    # patientData$refGenome<-as.character(input$refGenome)
    patientData$refGenome<-"hg19"
    
    ##Capturing runMode
    runMode_single<-as.character(input$runMode_single)
    
    ##Adding RunMode to patient Prediction info
    patientData$runMode<-runMode_single
    ###############################
    ## Computing Prediction
    ###############################
    
    
    patientResults<-list()
    
    #############################
    ## Try to get the prediction
    patientInfo<-patientData
    #When want to see error mssg uncomment the following line
    #patientResults<-masterWrapperSinglePrediction(patientInfo = patientInfo , minScore = minScore, highScore = highScore)
    patientResults<-tryCatch({
      ##If there is an error the following instruction will not be terminated
      
      patientResults<-masterWrapperSinglePrediction(patientInfo = patientInfo , minScore = minScore, highScore = highScore, runMode = runMode_single)
      ##If there was no error patientResults$Status == "OK" or "OK, but NO genes associated with SV"
    },error = function(err){
      patientResults$Status<-"ERROR"
      return(patientResults)
      
    })
    
    ############################################
    #If status error, generate the error html
    #We can be more specific in the future if we are interested
    ##If status = "OK, but NO genes associated with SV" generate HTML of NO genes found associated
    
    if(patientResults$Status=="ERROR"){
      #generate Error html
      ##output html in object$errorReport
      patientResults<-error_Report(patientResults = patientResults)
      
    }else if(patientResults$Status=="OK, but NO genes associated with SV"){
      ##For any of the dev stages of the phenotype, a gene was found.
      ##If a gene had been found for at least 1 devStage, there would appear a heatmap, but here makes no sense, nothing to sho
      ##Generate report No genes found html, located in object$statusReport
      ##output html in object$noGenesFoundReport
      patientResults<-noGenesFound_Report(patientResults = patientResults)
      
    }
    
    ##Stoping waiter
    remove_modal_spinner()
    #removeModal()
    
    return(patientResults)
  })
  
  #########################################################################
  ## Single SV submission From User Guide
  patientResults_2<-eventReactive(input$click_SingleSubmissionUserGuide, {
    ###Adding-Initializing waiter
    show_modal_spinner(text = HTML("Running Prediction <br>
                                   Please be patient"),
                       color="#0066ff")
    
    targetPatient<-input$singlePatientId
    filt_positive<-subset(usrGuide_SingleSubmTable, patientID == targetPatient)
    filt_positive$refGenome<-"hg19"
    ##Creating Input Data Frame, to carry on prediction
    ##Changing patientID by SV_ID for the html table output.
    patientData<-filt_positive
    colnames(patientData)<-c("chr_Break1","coord_Break1","chr_Break2",
                             "coord_Break2","TypeSV",
                             "Phenotype","SV_ID","refGenome")
    
    ##Capturing runMode
    runMode_single<-as.character(input$runMode_UserGuide_SingleSubmPat)
    
    ##Adding RunMode to patient Prediction info
    patientData$runMode<-runMode_single
    
    ###############################
    ## Computing Prediction
    ###############################
    
    
    patientResults<-list()
    
    #############################
    ## Try to get the prediction
    patientInfo<-patientData
    #When want to see error mssg uncomment the following line
    #patientResults<-masterWrapperSinglePrediction(patientInfo = patientInfo , minScore = minScore, highScore = highScore)
    patientResults<-tryCatch({
      ##If there is an error the following instruction will not be terminated
      
      patientResults<-masterWrapperSinglePrediction(patientInfo = patientInfo , minScore = minScore, highScore = highScore, runMode = runMode_single)
      ##If there was no error patientResults$Status == "OK" or "OK, but NO genes associated with SV"
    },error = function(err){
      patientResults$Status<-"ERROR"
      return(patientResults)
      
    })
    
    ############################################
    #If status error, generate the error html
    #We can be more specific in the future if we are interested
    ##If status = "OK, but NO genes associated with SV" generate HTML of NO genes found associated
    
    if(patientResults$Status=="ERROR"){
      #generate Error html
      ##output html in object$errorReport
      patientResults<-error_Report(patientResults = patientResults)
      
    }else if(patientResults$Status=="OK, but NO genes associated with SV"){
      ##For any of the dev stages of the phenotype, a gene was found.
      ##If a gene had been found for at least 1 devStage, there would appear a heatmap, but here makes no sense, nothing to sho
      ##Generate report No genes found html, located in object$statusReport
      ##output html in object$noGenesFoundReport
      patientResults<-noGenesFound_Report(patientResults = patientResults)
      
    }
    
    ##Stoping waiter
    remove_modal_spinner()
    #removeModal()
    
    return(patientResults)
  })
  
  
  #####################################################################
  ## Single Submission from the Explore Previous Patients Section
  ##We need to use the AllPatientsData dataframe used to generate those results
  patientResults_3<-eventReactive(input$click_AggregatedRes_ExplorePreviousPat, {
    ###Adding-Initializing waiter
    show_modal_spinner(text = HTML("Running Prediction <br>
                                   Please be patient"),
                       color="#0066ff")
    
    ###############################
    ## Computing Prediction
    ###############################
    patientResults<-list()
    
    #############################
    ## Try to get the prediction
    
    #When want to see error mssg uncomment the following line
    patientResults<-tryCatch({
      ##If there is an error the following instruction will not be terminated
      
      ##Preparing input
      ##Handling SVID
      ##Trimws to avoid problems with spaces before or after the word
      targetPatient<-trimws(input$aggregatedResults_svID_ExplorePreviousPat)
      
      ##Capturing selected phenotype
      targetPhenoSV<-input$aggregatedResults_phenoId_ExplorePreviousPat
      
      ##Filtering database patients information
      filt_InfoDBsvs<-subset(AllPatientsInfo,patientID == targetPatient )
      ##Ahora mismo los phenos estan entrando aqui con los nombres para filtrar la tabla
      ##Check pheno that user specified is indeed associated with the SV of interest
      sv_has_The_Pheno<-grepl(pattern = targetPhenoSV, x = filt_InfoDBsvs$Phenotype)
      
      ##Checking data is correct
      if((nrow(filt_InfoDBsvs)!=1) || (sv_has_The_Pheno == FALSE)){
        ##Apparently the relation sv-pheno does not exist in the database, or for any reason there is a repeated SVID or the SVid is wrongly introduced
        ##Rise error regarding wrong input.
        ##Did you properly introduce the SV ID? Did you properly select the phenotype for the SV ID? Check and re-submit
        stop("ERRORR WRONG INPUT")
      }
      ##Defining only one phenotype
      patientInfo<-filt_InfoDBsvs
      patientInfo$Phenotype<-targetPhenoSV
      
      ##Capturing runMode
      runMode_single<-as.character(input$runMode_AggregatedRes_ExplorePreviousPat)
      
      ##Adding RunMode to patient Prediction info
      patientInfo$runMode<-runMode_single
      
      ##Changing patientID by SV_ID for the html table output.
      colnames(patientInfo)[colnames(patientInfo)=="patientID"]<-"SV_ID"
      
      ###Running prediction
      patientResults<-masterWrapperSinglePrediction(patientInfo = patientInfo , minScore = minScore, highScore = highScore, runMode = runMode_single)
      ##If there was no error patientResults$Status == "OK" or "OK, but NO genes associated with SV"
    },error = function(err){
      patientResults$Status<-"ERROR"
      return(patientResults)
      
    })
    
    ############################################
    #If status error, generate the error html
    #We can be more specific in the future if we are interested
    ##If status = "OK, but NO genes associated with SV" generate HTML of NO genes found associated
    
    if(patientResults$Status=="ERROR"){
      #generate Error html
      ##output html in object$errorReport
      patientResults<-error_Report(patientResults = patientResults)
      
    }else if(patientResults$Status=="OK, but NO genes associated with SV"){
      ##For any of the dev stages of the phenotype, a gene was found.
      ##If a gene had been found for at least 1 devStage, there would appear a heatmap, but here makes no sense, nothing to sho
      ##Generate report No genes found html, located in object$statusReport
      ##output html in object$noGenesFoundReport
      patientResults<-noGenesFound_Report(patientResults = patientResults)
      
    }
    
    ##Stoping waiter
    remove_modal_spinner()
    #removeModal()
    
    return(patientResults)
  })
  
  #################################################################################################
  ## Submission from the Multiple SV Submission
  ## We need to use the Dataframe generated upon loading file to run Multiple Patients Submission
  patientResults_4<-eventReactive(input$click_AggregatedRes_MultipleSVSubmission, {
    ###Adding-Initializing waiter
    show_modal_spinner(text = HTML("Running Prediction <br>
                                   Please be patient"),
                       color="#0066ff")
    
    ###############################
    ## Computing Prediction
    ###############################
    patientResults<-list()
    
    #############################
    ## Try to get the prediction
    ##Check if the df "multiSV_uploadedFile_AllPatientsInfo" is available here
    #When want to see error mssg uncomment the following line
    ##DF with all SV uploaded, retrieving from the previously created object when handling Multiple SV Submission
    multiSV_uploadedFile_AllPatientsInfo<-multiple_patientResults()$patientsInfo
    
    patientResults<-tryCatch({
      ##If there is an error the following instruction will not be terminated
      
      ##Preparing input
      ##Handling SVID
      ##Trimws to avoid problems with spaces before or after the word
      targetPatient<-trimws(input$aggregatedResults_svID_MultipleSVSubmission)
      
      ##Capturing selected phenotype
      targetPhenoSV<-input$aggregatedResults_phenoId_MultipleSVSubmission
      
      ##Filtering database patients information
      filt_InfoDBsvs<-subset(multiSV_uploadedFile_AllPatientsInfo,patientID == targetPatient )
      ##Ahora mismo los phenos estan entrando aqui con los nombres para filtrar la tabla
      ##Check pheno that user specified is indeed associated with the SV of interest
      sv_has_The_Pheno<-grepl(pattern = targetPhenoSV, x = filt_InfoDBsvs$Phenotype)
      
      ##Checking data is correct
      if((nrow(filt_InfoDBsvs)!=1) || (sv_has_The_Pheno == FALSE)){
        ##Apparently the relation sv-pheno does not exist in the database, or for any reason there is a repeated SVID or the SVid is wrongly introduced
        ##Rise error regarding wrong input.
        ##Did you properly introduce the SV ID? Did you properly select the phenotype for the SV ID? Check and re-submit
        stop("ERRORR WRONG INPUT")
      }
      ##Defining only one phenotype
      patientInfo<-filt_InfoDBsvs
      patientInfo$Phenotype<-targetPhenoSV
      
      ##Capturing runMode
      runMode_single<-as.character(input$runMode_AggregatedRes_MultipleSVSubmission)
      
      ##Adding RunMode to patient Prediction info
      patientInfo$runMode<-runMode_single
      
      ##Changing patientID by SV_ID for the html table output.
      colnames(patientInfo)[colnames(patientInfo)=="patientID"]<-"SV_ID"
      
      ###Running prediction
      patientResults<-masterWrapperSinglePrediction(patientInfo = patientInfo , minScore = minScore, highScore = highScore, runMode = runMode_single)
      ##If there was no error patientResults$Status == "OK" or "OK, but NO genes associated with SV"
    },error = function(err){
      patientResults$Status<-"ERROR"
      return(patientResults)
      
    })
    
    ############################################
    #If status error, generate the error html
    #We can be more specific in the future if we are interested
    ##If status = "OK, but NO genes associated with SV" generate HTML of NO genes found associated
    
    if(patientResults$Status=="ERROR"){
      #generate Error html
      ##output html in object$errorReport
      patientResults<-error_Report(patientResults = patientResults)
      
    }else if(patientResults$Status=="OK, but NO genes associated with SV"){
      ##For any of the dev stages of the phenotype, a gene was found.
      ##If a gene had been found for at least 1 devStage, there would appear a heatmap, but here makes no sense, nothing to sho
      ##Generate report No genes found html, located in object$statusReport
      ##output html in object$noGenesFoundReport
      patientResults<-noGenesFound_Report(patientResults = patientResults)
      
    }
    
    ##Stoping waiter
    remove_modal_spinner()
    #removeModal()
    
    return(patientResults)
  })
  
  
  
  ####################################################
  ## Managing, creating Output for Single Submission
  ####################################################
  #If Not Error, Normal Output. So go checking if error
  ## Igual puedo mirar aqui que objeto tiene un contenido que pertoca
  
  ##MAIN RESULTS
  #We put it the last so that it is the first thing seen upon output rendered
  
  ##they are linked to their respective table output$respectiveTable TO AVOID OVERWRITTING THE PLOT VARIABLE
  
  #####################################################
  ##results coming from MAIN MENU, SINGLE SUBMISSION
  
  output$masterSummary_result<-renderUI(
    if(patientResults()$Status=="OK"){
      HTML(patientResults()$heatmapSummary)
      
    }else if(patientResults()$Status=="OK, but NO genes associated with SV"){
      HTML(patientResults()$noGenesFoundReport)
      
    }else if(patientResults()$Status=="ERROR"){
      HTML(patientResults()$errorReport)
    }
  )
  
  ##############################################################
  ##results coming from USER GUIDE, SINGLE SUBMISSION SECTION
  
  output$masterSummary_result_usrGuide_SingleSubmission<-renderUI(
    if(patientResults_2()$Status=="OK"){
      HTML(patientResults_2()$heatmapSummary)
      
    }else if(patientResults_2()$Status=="OK, but NO genes associated with SV"){
      HTML(patientResults_2()$noGenesFoundReport)
      
    }else if(patientResults_2()$Status=="ERROR"){
      HTML(patientResults_2()$errorReport)
    }
  )
  
  ##############################################################
  ##results coming from Explore Previous Patient single submission section
  
  output$masterSummary_result_explorePreviousPatient_SingleSubmission<-renderUI(
    if(patientResults_3()$Status=="OK"){
      HTML(patientResults_3()$heatmapSummary)
      
    }else if(patientResults_3()$Status=="OK, but NO genes associated with SV"){
      HTML(patientResults_3()$noGenesFoundReport)
      
    }else if(patientResults_3()$Status=="ERROR"){
      HTML(patientResults_3()$errorReport)
    }
  )
  
  
  #########################################################################################
  ##results coming from Multiple SV Submission Window, from the single submission section
  
  output$masterSummary_result_MultipleSVSubmission_SingleSubmission<-renderUI(
    if(patientResults_4()$Status=="OK"){
      HTML(patientResults_4()$heatmapSummary)
      
    }else if(patientResults_4()$Status=="OK, but NO genes associated with SV"){
      HTML(patientResults_4()$noGenesFoundReport)
      
    }else if(patientResults_4()$Status=="ERROR"){
      HTML(patientResults_4()$errorReport)
    }
  )
  
  
  ##############################################################################
  ##############################################################################
  ##Waiting for submission of multiple patient | SV
  
  ##When doing multiple patient and pheno analysis, considered pheno
  ##If not here,  not do prediction
  consideredPheno<-c("head_neck",
                     "cardiovascular",
                     "limbs",
                     "neurodevelopmental")##As more phenos considered they will appear here
  
  ############################################################
  ## Computing prediction for MULTIPLE Condition Submission
  ############################################################
  multiple_patientResults<-eventReactive(input$clickMultiple, {
    
    show_modal_spinner(text = HTML("Running Prediction <br>
                                   Please be patient"),
                       color="#0066ff")
    
    ##Retrieving dataframe for processing
    multiMetaData<-input$multipleFileInfo
    multiData<-read.delim(file = multiMetaData$datapath,
                          header = FALSE,
                          sep="\t",
                          stringsAsFactors = FALSE)
    
    colnames(multiData)<-c("chr_Break1","coord_Break1","chr_Break2",
                           "coord_Break2","TypeSV","Phenotype","patientID") ##For internal usage, patientID, for frontEnd SV_ID. To clearly show that a patient can carry multiple SVs
    
    rownames(multiData)<-multiData$patientID
    
    ######################################################################################
    ##Ensuring data in proper character format, for the coordinates, needed in character
    ######################################################################################
    multiData$chr_Break1<-as.character(multiData$chr_Break1)
    multiData$coord_Break1<-as.character(multiData$coord_Break1)
    
    multiData$chr_Break2<-as.character(multiData$chr_Break2)
    multiData$coord_Break2<-as.character(multiData$coord_Break2)
    #######################################
    multiple_patientResults<-list()
    
    multiple_patientResults$patientsInfo<-multiData
    
    ##Capturing runMode
    runMode_Multiple<-as.character(input$runMode_Multiple)
    
    ##############################
    ## Multiple Patient Analysis
    ##############################
    multiSV_uploadedFile_AllPatientsInfo<-multiData
    
    all_patientResults<-list()
    
    ###############################
    ## Computing Prediction
    
    ###Adding Prediction Progress bar
    ##The class of the notifier is 'shiny-notification' important to relocate it with css
    withProgress(message = 'Working on Structural Variant Nr', value = 0, {
      
      # Number of times we'll go through the loop
      n <- nrow(multiSV_uploadedFile_AllPatientsInfo)##max Number. Required for visual progress bar to know proportional size of iterations
      
      ###Tracker variables
      nPatient<-0  ##Tracking Status
      cohortTractablePhenos<-character() ##To track patient provided phenotypes, only sections related to patient phenotypes will be shown. No sense on showing cardiovascular, if no cardiovascular.
      ##Phenos considered when we also have data form them (eg if patient limb but yet no limb data also not section)
      for(patient in rownames(multiSV_uploadedFile_AllPatientsInfo)){
        print("                   ")
        print(patient)
        
        nPatient<-nPatient+1
        
        cat("n Structural Variant: ",nPatient,"\n")
        
        # Increment the progress bar, and update the detail text.
        ##This line is specific for the shiny html progressbar, for local running makes no sense
        incProgress(1/n, detail = paste(" ", nPatient))
        
        ########################################
        ## Selecting Patient Info
        patientInfo<-multiSV_uploadedFile_AllPatientsInfo[patient,]
        
        all_patientPheno<-unlist(strsplit(x = patientInfo$Phenotype, split = ",", fixed = TRUE))
        for(pheno in all_patientPheno){
          print("                   ")
          print(pheno)
          
          ##Working on each pheno separately, running prediction
          monoPheno_patientInfo<-patientInfo
          monoPheno_patientInfo$Phenotype<-pheno 
          print(monoPheno_patientInfo)
          
          ##If pheno among the ones considered in the app, run prediction
          if(pheno %in% consideredPheno){
            ############################
            ## Carrying prediction #####
            ############################
            
            cohortTractablePhenos<-c(cohortTractablePhenos, pheno)##We track this phenotype, and we will provide information for it in the multi pat section
            
            patientResults<-list()
            
            patientResults<-tryCatch({
              ##If there is an error the following instruction will not be terminated
              
              ##In multiple screening we do not generate graphics,so only master_scoring_function used, and not the wrapper for graphics one
              patientResults<-master_scoring_function(patientInfo = monoPheno_patientInfo, runMode = runMode_Multiple)              
              ##If there was no error patientResults$Status == "OK" or "OK, but NO genes associated with SV"
            },error = function(err){
              patientResults$Status<-"ERROR"
              return(patientResults)
              
            })
            
            ################################################
            ## Storaging Results per Patient & Phenotype
            ################################################
            ##To avoid the object size to be unnecessary huge, only mantain strictly necessary info
            ##We can maybe even simplify this more by using the matrix behind heatmap
            ##but let's see for now how it goes
            patientResults$resultsPerPhase_secondaryInfo<-NULL
            patientResults$genomeBrowser_links<-NULL
            patientResults$allAffectedGenes_positionalInfo<-NULL
            patientResults$MasterEnh_map<-NULL
            patientResults$resultsPerPhase<-NULL
            
            all_patientResults[[patient]][[pheno]]<-patientResults 
          } 
          
        }
        
      }
      
      
      
      ###################################################
      ## Parsing Patient Results Once Prediction is done
      ###################################################
      cohortTractablePhenos<-unique(cohortTractablePhenos)##Phenos that have appeared on the patients, and that we can predict impact
      ##For the phenos that appear in the patients
      
      cohort_results<-cohortResults_Parser(minScore = minScore, all_patientResults = all_patientResults,
                                           consideredPheno = cohortTractablePhenos,
                                           discardRelevantByBrokenGene = FALSE,
                                           AllPatientsInfo = multiSV_uploadedFile_AllPatientsInfo)
      
      #########################################
      ## Generating HTML from Parsed results
      #########################################
      
      multiple_patientResults$html_recurrency<-multipleStats_htmlGeneration(cohort_results = cohort_results, 
                                                                            consideredPheno = cohortTractablePhenos,
                                                                            ids_append="",##no append, for previous "previous pat"
                                                                            AllPatientsInfo = multiSV_uploadedFile_AllPatientsInfo,
                                                                            explPreviousPatSection = FALSE)
    })
    
    runjs(
      ##"document.getElementsByClassName('wrapperMainSingleResults')[0].style.visibility='hidden';"##With this we would just modify one instance of the class
      ##With the for, all elements belonging to the class change their status. With visibility hidden the space of the div is not erased
      "for (let el of document.querySelectorAll('.patientInputPanel')) el.style.display = 'grid';"
    )
    ##Stoping waiter
    remove_modal_spinner()
    #############################################################################
    #############################################################################
    ##anyadir resultados de pacientes para los tests, para ver que el correo depende del input
    multiple_patientResults$patientSpecificResults<-all_patientResults
    
    return(multiple_patientResults)
  })
  
  #######################################################
  ## Managing, creating Output for Multiple Submission
  #######################################################
  
  output$master_multipleStats<-renderUI(HTML(multiple_patientResults()$html_recurrency))
  
  
}
shinyApp(ui = ui, server = server)
