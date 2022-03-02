#######################################################
## No Genes Found Report. No Gene associated with SV
#######################################################
noGenesFound_Report<-function(patientResults){
  ##Probably be more specific in the future about the error & hence generate +specific error report
  
  ##Is the wrapper class for the main results, in case we want to hide show everything
  
  patientResults$noGenesFoundReport<-"<html>
  <body>
  <div class='wrapperMainSingleResults'>
  <h1>No gene associated with the structural variant</h1>
  <p style='font-size:25px' align='justify'>Based on currently used regulatory domain maps, no gene was found to be affected by the structural variant, either directly or by long-range mechanisms.</p>
  <p style='font-size:25px' align='justify'>Hence, pathogenic evidence is not found for this structural variant.</p>
  </div>
  </body>
  </html>"
  
  return(patientResults)
  
}