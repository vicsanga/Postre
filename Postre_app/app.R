########################################################
### POSTRE (Prediction of STRuctural variant Effects)
########################################################

library(shiny)
library(waiter)
library(shinybusy)
library(shinyjs)
library(shinyWidgets)
library(DT)
library(plotrix)##For enhancers, ellipse shape representation
library(shape)##For curve arrows representation
library(diagram)##For curve arrows representation

##To consider bioconductor repos
##To avoid this error: https://community.rstudio.com/t/deployment-error-unable-to-determine-the-location-for-some-packages/102312
# options(repos = BiocManager::repositories())
#options('repos')


###########################################################
###Setwd in the folder where all the app info is hosted
# setwd("~/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/Postre/Postre_app")

####################################
###Let's load required Functions
####################################
source("scripts_To_Load_Data/metaFunctionLoad.R")

##Required object for Single Prediction
##Loading multidata object, to avoid multiple reloading
load("data/MultiDataList.RData")

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
# navbar_js<-""

ui <-function(req){
  
  return(div(
    class="container",
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
                                                wellPanel(           
                                                  checkboxGroupInput(inputId = "phenoPatient", label = "Phenotype",
                                                                     c("Cardiovascular" = "cardiovascular",
                                                                       "Head & Neck" = "head_neck",
                                                                       "Limbs" = "limbs",
                                                                       "Neurodevelopmental" = "neurodevelopmental",
                                                                       "Vision-Eye" = "vision_eye",
                                                                       "Liver - Biliary System" = "liver_biliary_system"#,
                                                                       #"Endocrine" = "endocrine"
                                                                       ),
                                                                     selected = "head_neck",
                                                                     inline = TRUE)
                                                )
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
                                                            label = "Breakpoint 1 coordinates in hg19 (coord1,coord2 if interval)",
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
                                                            label = "Breakpoint 2 coordinates in hg19 (coord1,coord2 if interval)",
                                                            value = "99103873")
                                                )),
                                            ##Specifying reference Genome
                                            div(class="inp5",
                                                # wellPanel(selectInput(inputId = "refGenome_SingleSV_Subm", label = "Reference Genome",
                                                #                       choices = c("hg19","hg38"),
                                                #                       selected = "hg19"))
                                                
                                                wellPanel(
                                                  div(class="textRefGenome",
                                                      HTML('<h4><b>NOTE: Reference Genome Coordinates required in GRCh37/hg19</b></h4>'),
                                                      HTML('<p style="font-size:15px;">Please, visit any of the following websites to convert your coordinates to hg19 if you need it:
                                                         <a href= "http://genome.ucsc.edu/cgi-bin/hgLiftOver" target="_blank">UCSC LiftOver</a>,
                                                         <a href= "https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core" target="_blank">ENSEMBL Assembly Converter</a>,
                                                         <a href= "https://liftover.broadinstitute.org/" target="_blank">Broad Institute LiftOver</a>
                                                         </p>')
                                                  )
                                                )
                                                # ,
                                            ),
                                            
                                            #################################
                                            ##Adding Advanced Features Menu
                                            #################################
                                            
                                            div(class="inp6",
                                                wellPanel(HTML('
<div class="WrapperAdvancedFeatures">
<button class="collapsibleAdvancedFeatures" id="buttonSingleSubmAdvancedFeatures">Advanced Parameters [+]</button>
<div class="contentAdvancedFeatures">

  <!-- Selecting Running Mode -->
  <div class="form-group shiny-input-container">
  <label class="control-label" id="runMode_single-label" for="runMode_single">Running mode</label>
  <div>
    <select id="runMode_single"><option value="Standard" selected>Standard</option>
<option value="High-Specificity">High-Specificity</option></select>
    <script type="application/json" data-for="runMode_single" data-nonempty="">{"plugins":["selectize-plugin-a11y"]}</script>
  </div>
</div>

<!-- Considering gene previous association with patient phenotype, yes or no?-->
<div id="phenoConsideration_singleSubm" class="form-group shiny-input-radiogroup shiny-input-container shiny-input-container-inline" role="radiogroup" aria-labelledby="phenoConsideration_singleSubm-label">
  <label class="control-label" id="phenoConsideration_singleSubm-label" for="phenoConsideration_singleSubm">Gene-PatientPheno: Require known association of the candidate genes with the patient phenotype for pathogenicity</label>
                                                                 <div class="shiny-options-group">
                                                                 <label class="radio-inline">
                                                                 <input type="radio" name="phenoConsideration_singleSubm" value="yes" checked="checked"/>
                                                                 <span>Yes (default)</span>
                                                                 </label>
                                                                 <label class="radio-inline">
                                                                 <input type="radio" name="phenoConsideration_singleSubm" value="no"/>
                                                                 <span>No</span>
                                                                 </label>
                                                                 </div>
                                                                 </div>
                            
<hr>                    

<div class="tadSubmissionManagement">

<!-- Uploading your own TAD map -->
<div class="TADSubmission">
  <div class="form-group shiny-input-container">
    <label class="control-label" id="TADsInfo_single-label" for="TADsInfo_single">Upload TAD map</label>
    <div class="input-group">
      <label class="input-group-btn input-group-prepend">
        <span class="btn btn-default btn-file">
          Browse...
          <input id="TADsInfo_single" name="TADsInfo_single" type="file" style="position: absolute !important; top: -99999px !important; left: -99999px !important;" accept="*"/>
        </span>
      </label>
      <input type="text" class="form-control" placeholder="No file selected" readonly="readonly"/>
    </div>
    <div id="TADsInfo_single_progress" class="progress active shiny-file-input-progress">
      <div class="progress-bar"></div>
    </div>
  </div>
</div>

<button id="resetTADmap_singleSVsubm" type="button" class="action-button resetTADmapSelection resetSubmTAD">Remove selected TAD map</button>

</div>

</div>
</div>

<script>
  /*Function to expand-hidde the collapsible menu of advanced features*/
  var coll = document.getElementsByClassName("collapsibleAdvancedFeatures");
var i;

for (i = 0; i < coll.length; i++) {
  coll[i].addEventListener("click", function() {
    this.classList.toggle("active");
    var content = this.nextElementSibling;
    if (content.style.maxHeight){
      content.style.maxHeight = null;
    } else {
      content.style.maxHeight = content.scrollHeight + "px";
    } 
  });
}


  /*Function to change advanced features text depending if it is expanded or hidden*/
  const btn = document.getElementById("buttonSingleSubmAdvancedFeatures");
  btn.addEventListener("click", ()=>{

    if(btn.innerText === "Advanced Parameters [+]"){
        btn.innerText = "Advanced Parameters [-]";
    }else{
        btn.innerText= "Advanced Parameters [+]";
    }
  });
</script>
'),
                                                )
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
                                    div(class="formAndTitle_patient_Multiple_Input",
 
                                        HTML('<h2 id="titleSideBarPannel">Upload Multiple Structural Variants with Phenotypes</h2>'),
                                        div(class="patientMultipleFormulary",
                                            
                                            div(class="multipleSubmission",
                                                wellPanel(fileInput(
                                                  inputId = "multipleFileInfo",
                                                  label="Select File",
                                                  multiple = FALSE,
                                                  accept = "*"
                                                ))
                                            ),
                                            
                                            #################################
                                            ##Adding Advanced Features Menu
                                            #################################
                                            div(class="multipleSVRunMode",
                                                
                                                wellPanel(HTML('
<div class="WrapperAdvancedFeatures">
<button class="collapsibleAdvancedFeaturesMultiple" id="buttonMultipleSubmAdvancedFeatures">Advanced Parameters [+]</button>
<div class="contentAdvancedFeatures">

  <!-- Selecting Running Mode -->
  <div class="form-group shiny-input-container">
    <label class="control-label" id="runMode_Multiple-label" for="runMode_Multiple">Running mode</label>
    <div>
      <select id="runMode_Multiple"><option value="Standard" selected>Standard</option>
<option value="High-Specificity">High-Specificity</option></select>
      <script type="application/json" data-for="runMode_Multiple" data-nonempty="">{"plugins":["selectize-plugin-a11y"]}</script>
    </div>
  </div>

<!-- Considering gene previous association with patient phenotype, yes or no? -->
<div id="phenoConsideration_multipleSubm" class="form-group shiny-input-radiogroup shiny-input-container shiny-input-container-inline" role="radiogroup" aria-labelledby="phenoConsideration_multipleSubm-label">
  <label class="control-label" id="phenoConsideration_multipleSubm-label" for="phenoConsideration_multipleSubm">Gene-PatientPheno: Require known association of the candidate genes with the patient phenotype for pathogenicity</label>
                                                                 <div class="shiny-options-group">
                                                                 <label class="radio-inline">
                                                                 <input type="radio" name="phenoConsideration_multipleSubm" value="yes" checked="checked"/>
                                                                 <span>Yes</span>
                                                                 </label>
                                                                 <label class="radio-inline">
                                                                 <input type="radio" name="phenoConsideration_multipleSubm" value="no"/>
                                                                 <span>No</span>
                                                                 </label>
                                                                 </div>
                                                                 </div>
                                                                 
<hr>                    

<div class="tadSubmissionManagement">

<!-- Uploading your own TAD map -->
<div class="TADSubmission">
  <div class="form-group shiny-input-container">
    <label class="control-label" id="TADsInfo_multiple-label" for="TADsInfo_multiple">Upload TAD map</label>
    <div class="input-group">
      <label class="input-group-btn input-group-prepend">
        <span class="btn btn-default btn-file">
          Browse...
          <input id="TADsInfo_multiple" name="TADsInfo_multiple" type="file" style="position: absolute !important; top: -99999px !important; left: -99999px !important;" accept="*"/>
        </span>
      </label>
      <input type="text" class="form-control" placeholder="No file selected" readonly="readonly"/>
    </div>
    <div id="TADsInfo_multiple_progress" class="progress active shiny-file-input-progress">
      <div class="progress-bar"></div>
    </div>
  </div>
</div>

<button id="resetTADmap_multipleSVsubm" type="button" class="action-button resetTADmapSelection resetSubmTAD">Remove selected TAD map</button>

</div>
                                                                 
</div>
</div>

<script>
/*Function to expand-hidde the coll2apsible menu of advanced features*/
  var coll2 = document.getElementsByClassName("collapsibleAdvancedFeaturesMultiple");
  var i;
  
  for (i = 0; i < coll2.length; i++) {
    coll2[i].addEventListener("click", function() {
      this.classList.toggle("active");
      var content = this.nextElementSibling;
      if (content.style.maxHeight){
        content.style.maxHeight = null;
      } else {
        content.style.maxHeight = content.scrollHeight + "px";
      } 
    });
}


  /*Function to change advanced features text depending if it is expanded or hidden*/
  const btn2 = document.getElementById("buttonMultipleSubmAdvancedFeatures");
  btn2.addEventListener("click", ()=>{

    if(btn2.innerText === "Advanced Parameters [+]"){
        btn2.innerText = "Advanced Parameters [-]";
    }else{
        btn2.innerText= "Advanced Parameters [+]";
    }
  });
</script>
'),
                                                )
                          
                                            ),
                                            ##Specifying reference Genome
                                            div(class="multipleSVRefGenome",
                                                
                                                # wellPanel(selectInput(inputId = "refGenome_MultiSV_Subm", label = "Reference Genome",
                                                #                       choices = c("hg19","hg38"),
                                                #                       selected = "hg19"))
                                                
                                                wellPanel(
                                                  div(class="textRefGenome",
                                                      HTML('<h4><b>NOTE: Reference Genome Coordinates required in GRCh37/hg19</b></h4>'),
                                                      HTML('<p style="font-size:15px;">Please, visit any of the following websites to convert your coordinates to hg19 if you need it:
                                                         <a href= "http://genome.ucsc.edu/cgi-bin/hgLiftOver" target="_blank">UCSC LiftOver</a>,
                                                         <a href= "https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core" target="_blank">ENSEMBL Assembly Converter</a>,
                                                         <a href= "https://liftover.broadinstitute.org/" target="_blank">Broad Institute LiftOver</a>
                                                         </p>')
                                                  )
                                                )
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
  
  ##Defining phenotypes currently accepted by POSTRE, 
  ##Variable used in downstream computations
  ##When doing multiple patient and pheno analysis, considered pheno
  ##If not here,  not do prediction
  consideredPheno<-c("head_neck",
                     "cardiovascular",
                     "limbs",
                     "neurodevelopmental",
                     "vision_eye",
                     "liver_biliary_system"#,
                     #"endocrine"
                     )##As more phenos considered they will appear here
  ####################################################
  ##Defining initial behaviour when clicking buttons
  ####################################################
  
  ######################################################################
  ##Code to reset and capture the content of the TAD selected by user
  ######################################################################

  #####################################################
  ## For SingleSV submission, with code explanation
  #####################################################
  
  #Here we create an object that stores reactive values for user selected TAD map
  # It is similar to a list, but with special capabilities for reactive programming. 
  # When you read a value from it, the calling reactive expression takes a reactive dependency on that value, 
  # and when you write to it, it notifies any reactive functions that depend on that value.
  #It contains an element that we have named data, with initial value of NULL
  #The value can be accessed as rv_singleSVsubm$data or [['data']] in a reactive context
  # https://shiny.rstudio.com/reference/shiny/0.11/reactiveValues.html
  
  rv_singleSVsubm <- reactiveValues(data = NULL) ##To track the path for the USER selected TAD map

  ##Constantly monitoring the state of the input, if it changes, it also changes the value of the variable
  observe({
    req(input$TADsInfo_single)
    rv_singleSVsubm$data <- input$TADsInfo_single
  })
  
  #It triggers the action when the reset button is actioned
  #In this case if the rest button is actioned the info of the TADinput is erased
  observeEvent(input$resetTADmap_singleSVsubm, {
    reset('TADsInfo_single')##This does not delete the cached info (only the visual appearance) so we need the next step & the object created with reactiveValues
    rv_singleSVsubm$data <- NULL
  })
  
  #################################################
  ##For MultiSV submission advanced features menu
  #################################################
  rv_multiSVsubm <- reactiveValues(data = NULL) ##To track the path for the USER selected TAD map
  
  ##Constantly monitoring the state of the input, if it changes, it also changes the value of the variable
  observe({
    req(input$TADsInfo_multiple)
    rv_multiSVsubm$data <- input$TADsInfo_multiple
  })
  
  #It triggers the action when the reset button is actioned
  #In this case if the rest button is actioned the info of the TADinput is erased
  observeEvent(input$resetTADmap_multipleSVsubm, {
    reset('TADsInfo_multiple')##This does not delete the cached info (only the visual appearance) so we need the next step & the object created with reactiveValues
    rv_multiSVsubm$data <- NULL
  })
  
  
  #############################################################
  ##Object for submission from UserGuide, ExplorePreviousPat
  #############################################################
  # rv_userGuide <- reactiveValues(data = NULL) 
  # rv_explorePreviousPat <- reactiveValues(data = NULL) 
  rv_multiSVresults <- reactiveValues(data = NULL) ##This will capture the one of the last multiSVsubmission test this
  
  ######################################################################
  ######################################################################
  
  observeEvent(eventExpr = input$click, {
    updateTabsetPanel(session, 
                      inputId = "inTabset",
                      selected = "overview")
    
    runjs('
          document.getElementById("top").scrollIntoView();
          ')
    
    runjs(
      ##En la zona donde confluyen multiples plots hay que meter el remove, para que no se solapen funciones ocultas con mismos nombres
      # "for (let el of document.querySelectorAll('.wrapperMainSingleResults')) el.style.display = 'none';"
      "for (let el of document.querySelectorAll('.wrapperMainSingleResults')) el.remove();"
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
    
    runjs(
      ##En la zona donde confluyen multiples plots hay que meter el remove, para que no se solapen funciones ocultas con mismos nombres
      # "for (let el of document.querySelectorAll('.wrapperMainSingleResults')) el.style.display = 'none';"
      "for (let el of document.querySelectorAll('.wrapperMainSingleResults')) el.remove();"
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
     "for (let el of document.querySelectorAll('.wrapperMultipleSVSubmission')) el.style.display = 'none';"
      )
  })
  
  observeEvent(eventExpr = input$click_AggregatedRes_ExplorePreviousPat, {
    
    updateTabsetPanel(session, 
                      inputId = "inTabset",
                      selected = "overview")
    
    runjs('
          document.getElementById("top").scrollIntoView();
          ')

    runjs(
      ##En la zona donde confluyen multiples plots hay que meter el remove, para que no se solapen funciones ocultas con mismos nombres
      # "for (let el of document.querySelectorAll('.wrapperMainSingleResults')) el.style.display = 'none';"
      "for (let el of document.querySelectorAll('.wrapperMainSingleResults')) el.remove();"
    )
    
  })
  
  observeEvent(eventExpr = input$click_AggregatedRes_MultipleSVSubmission, {
    
    updateTabsetPanel(session, 
                      inputId = "inTabset",
                      selected = "overview")
    
    runjs('
          document.getElementById("top").scrollIntoView();
          ')
    
    runjs(
      ##En la zona donde confluyen multiples plots hay que meter el remove, para que no se solapen funciones ocultas con mismos nombres
      # "for (let el of document.querySelectorAll('.wrapperMainSingleResults')) el.style.display = 'none';"
      "for (let el of document.querySelectorAll('.wrapperMainSingleResults')) el.remove();"
      
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
    
    ##remove images generated in previous analyses, if there are
    deleteImages_Previous_Analyses()
    
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
    
    ##Capturing runMode
    runMode_single<-as.character(input$runMode_single)
    
    ##Adding RunMode to patient Prediction info
    patientData$runMode<-runMode_single
    
    ##Adding Consideration gene-Pheno to output
    patientData$genePhenoConsideration<-as.character(input$phenoConsideration_singleSubm)
    
    ##Dealing with referenceGenome
    # patientData$refGenome<-as.character(input$refGenome_SingleSV_Subm)
    patientData$refGenome<-"hg19"
    
    ##Default option for ownTAD, will be changed if user uploads a TAD map
    patientData$userTADmap<-"no"##Default option, used for UCSC warning, upload TAD coord
    user_tadMapInfo<-list()##It will be kept empty unless user provides its own TADmap
    
    ##Then maybe filter output of the patientData table
    ###############################
    ## Computing Prediction
    ###############################
    
    patientResults<-list() ##It will host multiplePhenotypes predictions at once
    
    
    userTadProcessed<-FALSE ##To track if own TAD has to be processed and doing it only once
    for(patientPhenotype in input$phenoPatient){
      
      patientResults_singlePhenoPrediction<-list()
      
      ##Assigning invidivual phenotype
      patientData$Phenotype<-patientPhenotype
      
      #############################
      ## Try to get the prediction
      
      #When want to see error mssg uncomment the following line
      #patientResults<-masterWrapperSinglePrediction(patientInfo = patientInfo , minScore = minScore, highScore = highScore)
      patientResults_singlePhenoPrediction<-tryCatch({
        
        ##########################################################################################
        ##Liftovering patient Coordinates if necessary, it will only be done once per patient
        # if(patientData$refGenome!="hg19"){
        #   patientData<-processing_GenomicCoordinates(patientInfo = patientData)
        #   patientData$refGenome<-"hg19" ##To inform that coordinates have been changed to hg19 and are those
        # }
        
        ##Check if user has selected its own TAD map for the analysis
        ##When TRUE is because NO TAD map selected
        ##For between TAD map generation the reference genome is taken into consideration for chr sizes
        ##So doing this after coordinates liftover
        
        if(is.null(rv_singleSVsubm$data) == FALSE){
          ##So, user selected TAD map
          ##Check that TADmap NOT already processed
          if(userTadProcessed == FALSE){
            
            ##So, not processed yet
            ##Do processing
            
            ##Reading TAD map in bed format
            
            ##Retrieving dataframe for processing
            userTADmapMetaData<-rv_singleSVsubm$data
            userTADmap<-read.delim(file = userTADmapMetaData$datapath,
                                  header = FALSE,
                                  sep="\t",
                                  stringsAsFactors = FALSE)
            
            colnames(userTADmap)<-c("chr","start","end") ##For internal usage, patientID, for frontEnd SV_ID. To clearly show that a patient can carry multiple SVs
            
            ##Filter for chromosomes of interest
            patientBreakpoints<-unique(c(patientData$chr_Break1, patientData$chr_Break2))
            
            ##We are not going to use info for other chr than the ones affected by the SV
            ##If it is Translocation will be 2 but if it is INV DEL DUP, it will only be one
            ##Thus to speed computations, skiping info of non required chromosomes
            
            userTADmap<-userTADmap[userTADmap$chr %in% patientBreakpoints,]
            
            
            ##Getting between TAD map
            userBetweenTADmap<-getBetweenTADmap(TADmap=userTADmap)
          
            ##Capturing processed information  
            user_tadMapInfo$TAD_map<-userTADmap
            user_tadMapInfo$Between_TAD_map<-userBetweenTADmap
            
            ##Recording that TAD map has already been processed
            userTadProcessed<-TRUE
          }
          
          ##To track wether user_tadMapInfo must be used or not in downstream functions
          patientData$userTADmap<-"yes"

          
        }
       
        ##If there is an error the following instruction will not be terminated
        # browser()
        patientResults_singlePhenoPrediction<-masterWrapperSinglePrediction(patientInfo = patientData , minScore = minScore, 
                                                                            highScore = highScore, runMode = runMode_single,
                                                                            user_tadMapInfo = user_tadMapInfo,
                                                                            MultiDataList = MultiDataList)
        ##If there was no error patientResults_singlePhenoPrediction$Status == "OK" or "OK, but NO genes associated with SV"
      },error = function(err){
        patientResults_singlePhenoPrediction$Status<-"ERROR"
        return(patientResults_singlePhenoPrediction)
        
      })
      
      ############################################
      #If status error, generate the error html
      #We can be more specific in the future if we are interested
      ##If status = "OK, but NO genes associated with SV" generate HTML of NO genes found associated
      
      if(patientResults_singlePhenoPrediction$Status=="ERROR"){
        #generate Error html
        ##output html in object$errorReport
        patientResults_singlePhenoPrediction<-error_Report(patientResults = patientResults_singlePhenoPrediction)
        
      }else if(patientResults_singlePhenoPrediction$Status=="OK, but NO genes associated with SV"){
        ##For any of the dev stages of the phenotype, a gene was found.
        ##If a gene had been found for at least 1 devStage, there would appear a heatmap, but here makes no sense, nothing to sho
        ##Generate report No genes found html, located in object$statusReport
        ##output html in object$noGenesFoundReport
        patientResults_singlePhenoPrediction<-noGenesFound_Report(patientResults = patientResults_singlePhenoPrediction)
        
      }
      
      patientResults[[patientPhenotype]]<-patientResults_singlePhenoPrediction
      
    }
    
    ###Merge results from the different phenotypes to create the master Output HTML
    patientResults<-mergingMultiPhenoPredictions(patientResults = patientResults, selectedPheno = input$phenoPatient, allPostreAvailablePheno = consideredPheno)
 
    ##Stoping waiter
    remove_modal_spinner()
    return(patientResults)
  })
  
  #########################################################################
  ## Single SV submission From User Guide
  patientResults_2<-eventReactive(input$click_SingleSubmissionUserGuide, {
    ###Adding-Initializing waiter
    show_modal_spinner(text = HTML("Running Prediction <br>
                                   Please be patient"),
                       color="#0066ff")
    
    ##remove images generated in previous analyses, if there are
    deleteImages_Previous_Analyses()
    
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
    runMode_single<-"Standard"
    
    ##Adding RunMode to patient Prediction info
    patientData$runMode<-runMode_single
    
    ##Adding Consideration gene-Pheno to output
    patientData$genePhenoConsideration<-"yes"
    
    patientData$userTADmap<-"no" ##Default option, used for UCSC warning when yes, upload TAD coord
    user_tadMapInfo<-list() ##It will be kept empty when no user provided TAD
    
    
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
      
      patientResults<-masterWrapperSinglePrediction(patientInfo = patientInfo , minScore = minScore,
                                                    highScore = highScore, runMode = runMode_single,
                                                    user_tadMapInfo = user_tadMapInfo,
                                                    MultiDataList = MultiDataList)
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
    
    ##remove images generated in previous analyses, if there are
    deleteImages_Previous_Analyses()
    
    ###############################
    ## Computing Prediction
    ###############################
    patientResults<-list()
    
    #############################
    ## Try to get the prediction
    
    ##We have error handling at two levels, the intial one, or inside the patient, for a specific phenotype (second level error)
    patientResults<-tryCatch({
      ##If there is an error the following instruction will not be terminated
      
      ##Preparing input
      ##Handling SVID
      ##Trimws to avoid problems with spaces before or after the word
      targetPatient<-trimws(input$aggregatedResults_svID_ExplorePreviousPat)
      
      ##Filtering database patients information
      filt_InfoDBsvs<-subset(AllPatientsInfo,patientID == targetPatient )

      ##Checking data is correct
      if(nrow(filt_InfoDBsvs)!=1){
        ##For any reason there is a repeated SVID or the SVid is wrongly introduced
        ##Rise error regarding wrong input.
        ##Did you properly introduce the SV ID? Did you properly select the phenotype for the SV ID? Check and re-submit
        stop("ERROR WRONG INPUT")
      }
      patientInfo<-filt_InfoDBsvs
      
      ##Capturing runMode
      ##Adding RunMode to patient Prediction info
      patientInfo$runMode<-"Standard"
      
      ##Adding Consideration gene-Pheno to output
      patientInfo$genePhenoConsideration<-"yes"
      
      patientInfo$userTADmap<-"no" ##Default option, used for UCSC warning when yes, upload TAD coord
      user_tadMapInfo<-list() ##It will be kept empty when no user provided TAD
      
      
      ##Changing patientID by SV_ID for the html table output.
      colnames(patientInfo)[colnames(patientInfo)=="patientID"]<-"SV_ID"
      
      ##Iterating over each phenotype
      patientAllPheno<-unlist(strsplit(x = patientInfo$Phenotype, split = ",", fixed = TRUE))
      
      for(patientPhenotype in patientAllPheno){
        patientResults_singlePhenoPrediction<-list()
        
        ##Check that pheno among the ones currently accepted by POSTRE
        if(patientPhenotype %in% consideredPheno){
          ##Do prediction
          ##Assigning invidivual phenotype
          patientInfo$Phenotype<-patientPhenotype
          
          ###Running prediction
          patientResults_singlePhenoPrediction<-tryCatch({
            ##If there is an error the following instruction will not be terminated
            patientResults_singlePhenoPrediction<-masterWrapperSinglePrediction(patientInfo = patientInfo , minScore = minScore, 
                                                                                highScore = highScore, runMode = patientInfo$runMode,
                                                                                user_tadMapInfo = user_tadMapInfo,
                                                                                MultiDataList = MultiDataList)
            ##If there was no error patientResults_singlePhenoPrediction$Status == "OK" or "OK, but NO genes associated with SV"
          },error = function(err){
            patientResults_singlePhenoPrediction$Status<-"ERROR"
            return(patientResults_singlePhenoPrediction)
            
          })
          
          ############################################
          #If status error, generate the error html
          #We can be more specific in the future if we are interested
          ##If status = "OK, but NO genes associated with SV" generate HTML of NO genes found associated
          
          if(patientResults_singlePhenoPrediction$Status=="ERROR"){
            #generate Error html
            ##output html in object$errorReport
            patientResults_singlePhenoPrediction<-error_Report(patientResults = patientResults_singlePhenoPrediction)
            
          }else if(patientResults_singlePhenoPrediction$Status=="OK, but NO genes associated with SV"){
            ##For any of the dev stages of the phenotype, a gene was found.
            ##If a gene had been found for at least 1 devStage, there would appear a heatmap, but here makes no sense, nothing to sho
            ##Generate report No genes found html, located in object$statusReport
            ##output html in object$noGenesFoundReport
            patientResults_singlePhenoPrediction<-noGenesFound_Report(patientResults = patientResults_singlePhenoPrediction)
            
          }
          
          patientResults[[patientPhenotype]]<-patientResults_singlePhenoPrediction
          
        }
      }
      
      ###Merge results from the different phenotypes to create the master Output HTML
      patientResults<-mergingMultiPhenoPredictions(patientResults = patientResults, selectedPheno = patientAllPheno, allPostreAvailablePheno = consideredPheno)
      
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
    return(patientResults)
  })
  
  #################################################################################################
  ## Submission of single SV from the Multiple SV Submission Results Section
  ## We need to use the Dataframe generated upon loading file to run Multiple Patients Submission
  ## And to extract the defined configuration (the one in the multiple patients df)
  
  patientResults_4<-eventReactive(input$click_AggregatedRes_MultipleSVSubmission, {
    ###Adding-Initializing waiter
    show_modal_spinner(text = HTML("Running Prediction <br>
                                   Please be patient"),
                       color="#0066ff")
    
    ##remove images generated in previous analyses, if there are
    deleteImages_Previous_Analyses()
    
    
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
    user_tadMapInfo<-multiple_patientResults()$user_tadMapInfo
    
    ##We have error handling at two levels, the intial one, or inside the patient, for a specific phenotype (second level error)
    patientResults<-tryCatch({
      ##If there is an error the following instruction will not be terminated
      
      ##Preparing input
      ##Handling SVID
      ##Trimws to avoid problems with spaces before or after the word
      targetPatient<-trimws(input$aggregatedResults_svID_MultipleSVSubmission)
      
      ##Filtering database patients information
      filt_InfoDBsvs<-subset(multiSV_uploadedFile_AllPatientsInfo,patientID == targetPatient )

      ##Checking data is correct
      if(nrow(filt_InfoDBsvs)!=1){
        ##Apparently the relation sv-pheno does not exist in the database, or for any reason there is a repeated SVID or the SVid is wrongly introduced
        ##Rise error regarding wrong input.
        ##Did you properly introduce the SV ID? Did you properly select the phenotype for the SV ID? Check and re-submit
        stop("ERRORR WRONG INPUT")
      }
      patientInfo<-filt_InfoDBsvs
      
      ###
      ##Running configuration as runMode etc is already recorded on the patient Info from the multipleSV submission
      ## eg in patientInfo$runMode or patientInfo$genePhenoConsideration
      ##Regarding refGenome ifit was necessary liftover, here the coordinates already liftovered

      ##Changing patientID by SV_ID for the html table output.
      colnames(patientInfo)[colnames(patientInfo)=="patientID"]<-"SV_ID"
      
      ##Iterating over each phenotype
      patientAllPheno<-unlist(strsplit(x = patientInfo$Phenotype, split = ",", fixed = TRUE))
      
      for(patientPhenotype in patientAllPheno){
        
        patientResults_singlePhenoPrediction<-list()
        
        ##Check that pheno among the ones currently accepted by POSTRE
        if(patientPhenotype %in% consideredPheno){
          ##Do prediction
          ##Assigning invidivual phenotype
          patientInfo$Phenotype<-patientPhenotype
          
          ###Running prediction
          patientResults_singlePhenoPrediction<-tryCatch({
            
            ##If there is an error the following instruction will not be terminated
            patientResults_singlePhenoPrediction<-masterWrapperSinglePrediction(patientInfo = patientInfo , minScore = minScore,
                                                                                highScore = highScore, runMode = patientInfo$runMode,
                                                                                user_tadMapInfo = user_tadMapInfo,
                                                                                MultiDataList = MultiDataList)
            ##If there was no error patientResults_singlePhenoPrediction$Status == "OK" or "OK, but NO genes associated with SV"
          },error = function(err){
            patientResults_singlePhenoPrediction$Status<-"ERROR"
            return(patientResults_singlePhenoPrediction)
            
          })
          
          ############################################
          #If status error, generate the error html
          #We can be more specific in the future if we are interested
          ##If status = "OK, but NO genes associated with SV" generate HTML of NO genes found associated
          
          if(patientResults_singlePhenoPrediction$Status=="ERROR"){
            #generate Error html
            ##output html in object$errorReport
            patientResults_singlePhenoPrediction<-error_Report(patientResults = patientResults_singlePhenoPrediction)
            
          }else if(patientResults_singlePhenoPrediction$Status=="OK, but NO genes associated with SV"){
            ##For any of the dev stages of the phenotype, a gene was found.
            ##If a gene had been found for at least 1 devStage, there would appear a heatmap, but here makes no sense, nothing to sho
            ##Generate report No genes found html, located in object$statusReport
            ##output html in object$noGenesFoundReport
            patientResults_singlePhenoPrediction<-noGenesFound_Report(patientResults = patientResults_singlePhenoPrediction)
            
          }
          
          patientResults[[patientPhenotype]]<-patientResults_singlePhenoPrediction
        }
      }
      
      ###Merge results from the different phenotypes to create the master Output HTML
      patientResults<-mergingMultiPhenoPredictions(patientResults = patientResults, selectedPheno = patientAllPheno, allPostreAvailablePheno = consideredPheno)
      
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
    HTML(patientResults()$heatmapSummary)
    
  )
  
  ##############################################################
  ##results coming from USER GUIDE, SINGLE SUBMISSION SECTION
  
  output$masterSummary_result_usrGuide_SingleSubmission<-renderUI(
    HTML(patientResults_2()$heatmapSummary)
  )
  
  ##############################################################
  ##results coming from Explore Previous Patient single submission section
  
  output$masterSummary_result_explorePreviousPatient_SingleSubmission<-renderUI(
    HTML(patientResults_3()$heatmapSummary)
    
  )
  
  
  #########################################################################################
  ##results coming from Multiple SV Submission Window, from the single submission section
  
  output$masterSummary_result_MultipleSVSubmission_SingleSubmission<-renderUI(
    HTML(patientResults_4()$heatmapSummary)
  )
  
  ############################################################
  ## Computing prediction for MULTIPLE SV Submission PAGE
  ############################################################
  multiple_patientResults<-eventReactive(input$clickMultiple, {
    
    show_modal_spinner(text = HTML("Running Prediction <br>
                                   Please be patient"),
                       color="#0066ff")
    
    #########################################################################################
    ## During file loading we could get an error if file not selected or format not correct
    ## So start tracking error from here
    #########################################################################################
    multiple_patientResults<-list()
    multiple_patientResults<-tryCatch({    
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
      ##Capturing runMode
      runMode_Multiple<-as.character(input$runMode_Multiple)
      
      ##############################
      ## Multiple Patient Analysis
      ##############################
      multiSV_uploadedFile_AllPatientsInfo<-multiData
      
      ##Adding RunMode info
      multiSV_uploadedFile_AllPatientsInfo$runMode<-runMode_Multiple
      
      ##Adding Consideration gene-Pheno to output
      multiSV_uploadedFile_AllPatientsInfo$genePhenoConsideration<-as.character(input$phenoConsideration_multipleSubm)
      
      ##Dealing with referenceGenome IF NOT HG19, LIFTOVER to it. We will do it on each patient screen
      # multiSV_uploadedFile_AllPatientsInfo$refGenome<-as.character(input$refGenome_MultiSV_Subm)
      multiSV_uploadedFile_AllPatientsInfo$refGenome<-"hg19"
      
      ##Default option for ownTAD, will be changed if user uploads a TAD map
      multiSV_uploadedFile_AllPatientsInfo$userTADmap<-"no"##Default option, used for UCSC warning, upload TAD coord
      user_tadMapInfo<-list()##It will be kept empty unless user provides its own TADmap
      
      ###############################
      ## Computing Prediction
      all_patientResults<-list()
      
      ###Adding Prediction Progress bar
      ##The class of the notifier is 'shiny-notification' important to relocate it with css
      withProgress(message = 'Working on Structural Variant Nr', value = 0, {
        
        # Number of times we'll go through the loop
        n <- nrow(multiSV_uploadedFile_AllPatientsInfo)##max Number. Required for visual progress bar to know proportional size of iterations
        
        ###Tracker variables
        nPatient<-0  ##Tracking Status
        cohortTractablePhenos<-character() ##To track patient provided phenotypes, only sections related to patient phenotypes will be shown. No sense on showing cardiovascular, if no cardiovascular.
        ##Phenos considered when we also have data form them (eg if patient limb but yet no limb data also not section)
        
        userTadProcessed<-FALSE ##To track if own TAD has to be processed and doing it only once 
        start.time<-Sys.time()
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
            
            ##If pheno among the ones considered in the app, run prediction
            if(pheno %in% consideredPheno){
              ############################
              ## Carrying prediction #####
              ############################
              cohortTractablePhenos<-c(cohortTractablePhenos, pheno)##We track this phenotype, and we will provide information for it in the multi pat section
      
              patientResults<-list()
              
              patientResults<-tryCatch({
                
                ##If patient has the coordinates in refGenome different than hg19, the first time 
                ##The tool enters here, that info will be modified
                # if(patientInfo$refGenome!="hg19"){
                #   patientInfo<-processing_GenomicCoordinates(patientInfo = patientInfo)
                #   patientInfo$refGenome<-"hg19" ##To inform that coordinates have been changed to hg19 and are those
                #   
                #   ##Updating patientInfo with hg19 coordinates, so that it is clarified to which coordinates are converted
                #   multiSV_uploadedFile_AllPatientsInfo[patient,]<-patientInfo
                # }
                
                ##Check if user has selected its own TAD map for the analysis
                ##When TRUE is because NO TAD map selected
                ##For between TAD map generation the reference genome is taken into consideration for chr sizes
                ##So doing this after coordinates liftover
                
                if(is.null(rv_multiSVsubm$data) == FALSE){
                  
                  ##So, user selected TAD map
                  ##Check that TADmap NOT already processed
                  if(userTadProcessed == FALSE){
                    
                    ##So, not processed yet
                    ##Do processing
                    
                    ##Reading TAD map in bed format
                    ##Retrieving dataframe for processing
                    userTADmapMetaData<-rv_multiSVsubm$data
                    userTADmap<-read.delim(file = userTADmapMetaData$datapath,
                                           header = FALSE,
                                           sep="\t",
                                           stringsAsFactors = FALSE)
                    
                    colnames(userTADmap)<-c("chr","start","end") 
                    
                    ##Filter for chromosomes of interest
                    ## Consider all chr among all patients
                    patientBreakpoints<-unique(c(multiSV_uploadedFile_AllPatientsInfo$chr_Break1,
                                                 multiSV_uploadedFile_AllPatientsInfo$chr_Break2))
                    
                    ##Thus to speed computations, skiping info of non required chromosomes
                    
                    userTADmap<-userTADmap[userTADmap$chr %in% patientBreakpoints,]
                    
                    ##Getting between TAD map
                    userBetweenTADmap<-getBetweenTADmap(TADmap=userTADmap)
                    
                    ##Capturing processed information  
                    user_tadMapInfo$TAD_map<-userTADmap
                    user_tadMapInfo$Between_TAD_map<-userBetweenTADmap
                    
                    ##Recording that TAD map has already been processed
                    userTadProcessed<-TRUE
                  }
                  
                  ##To track whether user_tadMapInfo must be used or not in downstream functions
                  patientInfo$userTADmap<-"yes"
                  
                  ##Updating patientInfo with userTADmap "yes", so that it is clarified that the user tad map has been used for the prediction
                  multiSV_uploadedFile_AllPatientsInfo$userTADmap<-"yes" 
                }
                
                ##Working on each pheno separately, running prediction
                monoPheno_patientInfo<-patientInfo
                monoPheno_patientInfo$Phenotype<-pheno 
                
                ##If there is an error the following instruction will not be terminated
                ##In multiple screening we do not generate graphics,so only master_scoring_function used, and not the wrapper for graphics one
                patientResults<-master_scoring_function(patientInfo = monoPheno_patientInfo, runMode = runMode_Multiple, user_tadMapInfo=user_tadMapInfo, MultiDataList = MultiDataList)   
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
        end.time<-Sys.time()
        time.taken <- end.time - start.time
        print(time.taken)
        
        
        ###################################################
        ## Parsing Patient Results Once Prediction is done
        ###################################################
        cohortTractablePhenos<-unique(cohortTractablePhenos)##Phenos that have appeared on the patients, and that we can predict impact
        ##For the phenos that appear in the patients
        
        cohort_results<-cohortResults_Parser(minScore = minScore, all_patientResults = all_patientResults,
                                             consideredPheno = cohortTractablePhenos,
                                             discardRelevantByBrokenGene = FALSE,
                                             AllPatientsInfo = multiSV_uploadedFile_AllPatientsInfo)

        # browser() ##To store some table output
        # candidateGenesInfo<-cohort_results$candidateGenesInfo
        # And now save object
        
        #########################################
        ## Generating HTML from Parsed results
        #########################################
        
        multiple_patientResults$html_recurrency<-multipleStats_htmlGeneration(cohort_results = cohort_results, 
                                                                              consideredPheno = cohortTractablePhenos,
                                                                              ids_append="",##no append, for previous "previous pat"
                                                                              AllPatientsInfo = multiSV_uploadedFile_AllPatientsInfo,
                                                                              explPreviousPatSection = FALSE)
      })
      
      #############################################################################
      #############################################################################
      ##anyadir resultados de pacientes para los tests, para ver que el correo depende del input
      multiple_patientResults$patientSpecificResults<-all_patientResults
      
      ##Storing the final resulting object with patients info with liftovered data if required
      multiple_patientResults$patientsInfo<-multiSV_uploadedFile_AllPatientsInfo
      
      ##If the code reaches this point, no global error, due to file format or anything arised
      multiple_patientResults$Status<-"OK"
      
      ##Adding the userTADinfo since it will be used in the MultiSVsubmission results page
      multiple_patientResults$user_tadMapInfo<-user_tadMapInfo
      
      ##Lets assign to the variable associated with tryCatch (variable<-tryCatch) the content of interest: 
      ##In this case, the multiple_patientResults object
      ##We do that by just printing its content. However, for the error we need to return the object
      ## Ceck script testingTryCatch.R
      
      multiple_patientResults
      
    },error = function(err){
      ##A global error arised (not during the processing of a particular SV)
      multiple_patientResults$Status<-"ERROR"
      multiple_patientResults<-error_Report_multiSVsubm(multiple_patientResults = multiple_patientResults)
      
      ##We need to return the error report
      return(multiple_patientResults)
    })
    
    ##Stoping waiter AND Returning results
    remove_modal_spinner()
    return(multiple_patientResults)
    
  })
  
  
  #######################################################
  ## Managing, creating Output for Multiple Submission
  #######################################################
  
  output$master_multipleStats<-renderUI(
    if(multiple_patientResults()$Status == "OK"){
      HTML(multiple_patientResults()$html_recurrency) 
      
    }else if(multiple_patientResults()$Status == "ERROR"){
      HTML(multiple_patientResults()$errorReport) 
    }
  )
  
  
}
shinyApp(ui = ui, server = server)
