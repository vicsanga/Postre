########################################
##### Report Summary Generation
########################################
body_report_GeneSummary<-function(patientResults, reportUnit,minRequiredScore, targetGene){
  
  #############################################
  ## Adding short index of all result sections
  #############################################
  
  ####For one of the section entries we need to know whether the Pathomechanism is GOF or LOF, because it appear on the header
  
  targetGene<-unlist(strsplit(x = reportUnit, split = "_", fixed = TRUE))[1]
  targetMech<-unlist(strsplit(x = reportUnit, split = "_", fixed = TRUE))[2]
  targetPhase<-unlist(strsplit(x = reportUnit, split = "_", fixed = TRUE))[3]
  
  #This way of approaching the GOF-LOF info, deprecated
  # pathoMechanism<-patientResults$resultsPerPhase[[targetPhase]][[targetPhase]][targetGene,]$Mechanism
  # pathoMechanism<-unlist(strsplit(x = pathoMechanism,
  #                                 split = ":",
  #                                 fixed = TRUE))[2]##GOF or LOF
  
  pathoMechanism<-targetMech
  
  if(pathoMechanism=="LOF"){
    conditionProperName<-"Loss of Function"
    
  }else if(pathoMechanism=="GOF"){
    conditionProperName<-"Gain of Function"
    
  }
  
  #########################
  #Continuing with index creation
  ##Separar el indice por modulos para facilitar el escalado, y no como lo estaba haciendo
  
  indexLines_GraphicalAbstract<-paste("<li><a href='#GraphicalAbstractSection_",
      reportUnit,
      "'>Graphical Abstract </a></li>",
      sep="")
  
  indexLines_RegLandscape<-paste("<li><a href='#RegulatoryLandscapeChangesSection_",
                                      reportUnit,
                                      "'>Gene Impact - Regulatory Landscape Changes </a></li>",
                                      sep="")
  
  indexLines_GenomeBrower<-paste("<li><a href='#GenomeBrowserSection_",
                                 reportUnit,
                                 "'>Visualize affected loci in genomic browser </a></li>",
                                 sep="")
  
  indexLines_GenePheno<-paste("<li><a href='#GenePhenotypeSection_",
                                 reportUnit,
                                 "'>Gene previous association with the patient phenotype </a></li>",
                                 sep="")
  
  indexLines_GeneSusceptibility<-paste("<li><a href='#GeneSusceptibilitySection_",
                              reportUnit,
                              "'>Gene Susceptibility to ",
                              conditionProperName,
                              "</a></li>",
                              sep="")
  
  ##Gene Susceptibility to Loss of Function // or Gain of Functin section line
  
  indexResults<-paste("<ul class='reportIndex'>",
                      indexLines_GraphicalAbstract,
                      indexLines_RegLandscape,
                      indexLines_GeneSusceptibility,
                      indexLines_GenePheno,
                      indexLines_GenomeBrower,
                      "</ul>",
                      sep="")
  

  indexResults<-gsub("[\r\n]", "", indexResults)##remove line breaks

  
  ################################
  ## Adding Graphical Abstract
  ################################
  
  graphicalAbstract<-paste0("<img class='mainGraphicalAbstract' src='graphicalSummaries/", reportUnit,"_", patientResults$job_UniCode, ".png", "'",
                         " >")
  
  ##Explanation about why this sizes below
  ##with inches, no image appears,
  ##Regarding https://www.w3schools.com/cssref/css_units.asp
  ## 1 in corresponds with 96px (approx), hence as we have created the images with height = 8in and width = 12in
  ## The conversion is " height=768px, width=1152px >")
  
  graphicalAbstract<-paste("<h2 class='gene_FullReport_SectionEntry' id='",
                           "GraphicalAbstractSection_",
                           reportUnit,
                           "'>Graphical Abstract</h2>",
                           "<div class='gene_FullReport_SectionContent'>",
                           graphicalAbstract,
                           "</div>",
                           sep="")
  
  ########################################
  ## Changes Regulatory Landscape Section
  ########################################
  ##Explain differences enh landscape changes

  source("functions/report_regulatoryLandscapeChanges.R",
         local = TRUE)
  
  ##if mech direct effect, maybe avoid talking about enh and discarding content
  ##etc. maybe we can add here diversity between reports
  
  ##Generating Report
  reportChangesRegulatoryLandscape<-report_regulatoryLandscapeChanges(reportUnit = reportUnit,
                                                                      patientResults = patientResults)
  
  ##Embedding report on overall report html.
  reportChangesRegulatoryLandscape<-paste("<h2 class='gene_FullReport_SectionEntry' id='",
                           "RegulatoryLandscapeChangesSection_",
                           reportUnit,
                           "'>Gene Impact - Regulatory Landscape Changes</h2>",
                           "<div class='gene_FullReport_SectionContent'>",
                           reportChangesRegulatoryLandscape,
                           "</div>",
                           sep="")
  
  ##Done at the end of the script The final integration and order definition
  
  ###############################################################
  ## Gene Susceptibility to [Gain of Function | Loss of Function]
  ###############################################################
  
  source("functions/report_geneSusceptibilityToCondition.R")
  
  ##########################################################
  ##Generating Report
  reportGeneSusceptibilityToCondition<-report_geneSusceptibilityToCondition(reportUnit = reportUnit,
                                                                      patientResults = patientResults)

  ##Embedding report on overall report html.
  reportGeneSusceptibilityToCondition<-paste("<h2 class='gene_FullReport_SectionEntry' id='",
                                          "GeneSusceptibilitySection_",
                                          reportUnit,
                                          "'>Gene Susceptibility to ",
                                          conditionProperName,
                                          "</h2>",
                                          "<div class='gene_FullReport_SectionContent'>",
                                          reportGeneSusceptibilityToCondition,
                                          "</div>",
                                          sep="")
  
  ##########################
  ## UCSC links section ####
  ##########################
  
  ##Adding UCSC links
  source(file = "functions/ucsc_view.R",
         local = TRUE)
  #Defining UCSC session id
  pheno<-patientResults$resultsPerPhase_secondaryInfo[[1]][[1]]$formatedPhenotype
  devStage<-unlist(strsplit(x = reportUnit, split = "_", fixed = TRUE))[3]
  browserSessionId<-paste("SV_app_",
                          pheno,
                          "_",
                          devStage,
                          sep="")
  
  ucscLinksInfo<-ucsc_view(patientResults = patientResults,
                           browserSessionId = browserSessionId,
                           targetGene = targetGene,
                           devStage = devStage)##all enh from this devStage in the browserWindow are colored
  
  ucscLinksInfo<-paste("<h2 class='gene_FullReport_SectionEntry' id='",
                       "GenomeBrowserSection_",
                       reportUnit,
                       "'>Visualize affected loci in genomic browser </h2>",
                       "<div class='gene_FullReport_SectionContent'>",
                       ucscLinksInfo,
                       "</div>",
                       sep="")
  
  
  ##################################
  ## Gene-Phenotype Information ####
  ##################################
  
  
  ##################################
  # Function to create gene report
  source(file = "functions/Gene_textReport.R",
         local = TRUE)
  
  gene_old_TextReport<-gene_textReport(patientResults = patientResults, 
                                       minPatogenicScore = minRequiredScore,
                                       mainPhenotype = patientResults$resultsPerPhase_secondaryInfo[[1]][[1]]$formatedPhenotype,
                                       targetGene = targetGene)
  
  ##Assembling gene-pheno info
  report_gene_PhenoAssociation<-paste("<h2 class='gene_FullReport_SectionEntry' id='",
                                "GenePhenotypeSection_",
                                reportUnit,
                                "'>Gene previous association with the patient phenotype category</h2>",
                                sep="")
  
  ## Deprecated last version
  # if(pathoMechanism=="GOF"){
  #   #Indicating that gene-pheno information not considered to build score for GOF but provided as additional information.
  #   report_gene_PhenoAssociation<-paste(report_gene_PhenoAssociation,
  #                                       "<p class='warning_gof_genePheno'>Gene-Phenotype information is not considered when computing the pathogenic score for 'Gain Of Function', but it is provided as additional information.<br>Visit the User Guide to know more about the scoring system.</p>",
  #                                       sep="")
  #   
  # }
  
  report_gene_PhenoAssociation<-paste(report_gene_PhenoAssociation,
                             "<div class='gene_FullReport_SectionContent'>",
                             gene_old_TextReport,
                             "</div>",
                             sep="")
  
  ###############################################
  ## Pack all info report body
  #info styles 
  
  reportBody<-paste(indexResults,
                    graphicalAbstract,
                    "<br>",
                    "<hr class='reportHr'>",
                    "<br>",
                    reportChangesRegulatoryLandscape,
                    "<br>",
                    "<hr class='reportHr'>",
                    "<br>",
                    reportGeneSusceptibilityToCondition,
                    "<br>",
                    "<hr class='reportHr'>",
                    "<br>",
                    report_gene_PhenoAssociation,
                    "<br>",
                    "<hr class='reportHr'>",
                    "<br>",
                    ucscLinksInfo,
                    sep="",
                    collapse="")
  
  reportBody<-paste0(reportBody,"<br><br>")
  
  
  ########################
  ## returning report body
  return(reportBody)
  
}