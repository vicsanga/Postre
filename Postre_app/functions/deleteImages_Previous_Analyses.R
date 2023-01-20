###################################################
## Check and delete images from previous analyses
##If there are
## To minimize memory problems

deleteImages_Previous_Analyses<-function(){
  files_paths<-list.files(path = 'www/graphicalSummaries/')
  if(length(files_paths) !=0){
    ##Hence there is at least one image to delete
    system('rm www/graphicalSummaries/*')  
  }
}