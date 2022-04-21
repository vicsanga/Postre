###################################
## Error Report: Something failed
###################################
error_Report_multiSVsubm<-function(multiple_patientResults){
  ##Probably be more specific in the future about the error & hence generate +specific error report
  
  ##Is the wrapper class for the main results, in case we want to hide show everything
  
  multiple_patientResults$errorReport<-paste("<html>
  <body>
  <div class='wrapperMultipleSVSubmission'>
  <h1>Unable to terminate the prediction</h1>
  <p style='font-size:25px' align='justify'>We are sorry to inform you that a problem occurred while analyzing the multiple SVs file.</p>
  <p style='font-size:25px' align='justify'> <b>Please, check that you submitted the file in the correct format.</b> Postre requires a txt file with tabulation separated values.
  Visit the User Guide if you have any doubt or check the youtube video where it is explained how to submit Multiple SVs.</p>
  ",
  ' <br><br>
    <div class ="videotutorialExplPreviousPat"> 
      <p align="center">
        <iframe width="640" height="360" src="https://www.youtube.com/embed/fm8gGggzf4E" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
          </p>
          </div>',
  "
  <p style='font-size:25px' align='justify'> If the problem persists contact us through postre.radaiglesiaslab@gmail.com to provide you more information about the particularities of this case.
  </div>
  </body>
  </html>", 
  sep="")
  
  return(multiple_patientResults)
  
}