###################################
## Error Report: Something failed
###################################
error_Report<-function(patientResults){
  ##Probably be more specific in the future about the error & hence generate +specific error report
  
  ##Is the wrapper class for the main results, in case we want to hide show everything
  
  patientResults$errorReport<-"<html>
  <body>
  <div class='wrapperMainSingleResults'>
  <h1>Unable to terminate the prediction</h1>
  <p style='font-size:25px' align='justify'>We are sorry to inform you that a problem occurred while predicting the impact of the structural variant due to biological limitations or technical issues related with the patient information and affected loci.</p>
  <p style='font-size:25px' align='justify'>Please, check if the data introduced is correct and visit the User Guide if you have any doubt. If the problem persists contact us through postre.radaiglesiaslab@gmail.com to provide you more information about the particularities of this case.
  </div>
  </body>
  </html>"
  
  patientResults$heatmapSummary<-"<html>
  <body>
  <div class='wrapperMainSingleResults'>
  <h1>Unable to terminate the prediction</h1>
  <p style='font-size:25px' align='justify'>We are sorry to inform you that a problem occurred while predicting the impact of the structural variant due to biological limitations or technical issues related with the patient information and affected loci.</p>
  <p style='font-size:25px' align='justify'>Please, check if the data introduced is correct and visit the User Guide if you have any doubt. If the problem persists contact us through postre.radaiglesiaslab@gmail.com to provide you more information about the particularities of this case.
  </div>
  </body>
  </html>"
  
  return(patientResults)
  
}