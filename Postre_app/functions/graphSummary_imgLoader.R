#######################################
## Image Loader Temporary Script
## Script designed to create an HTML containing
## the graph summary images
#######################################

graphSummary_imgLoader<-function(){
  ##It does not need input parameters
  folderPath<-"www/graphicalSummaries/"
  #folderPath<-"graphicalSummaries/"
  img_paths<-list.files(path = folderPath)
  
  ##Creating HTML
  htmlGeneration<-"<html><body><div class='graphSummarySection'>"
  
  for(imageToLoad in img_paths){
    imgId<-unlist(strsplit(x = imageToLoad,
                           spli = ".",
                           fixed=TRUE))[1]##to remove png extension
    htmlGeneration<-paste0(htmlGeneration,
                           "<img src='graphicalSummaries/", imageToLoad, "' id='",imgId, "',",
                           " height=768px, width=1152px >")##Explanation about why this sizes below
    ##with inches, no image appears,
    ##Regarding https://www.w3schools.com/cssref/css_units.asp
    ## 1 in corresponds with 96px (approx), hence as we have created the images with height = 8in and width = 12in
    ## The conversion is " height=768px, width=1152px >")
    htmlGeneration<-paste0(htmlGeneration,"<br><br>")
    
  }
  
  ##Closing html
  htmlGeneration<-paste0(htmlGeneration,
                         "</div></body></html>")
  ##Let's remove all line breaks
  # htmlGeneration<-gsub("[\r\n]", "", htmlGeneration)
  return(htmlGeneration)
}
