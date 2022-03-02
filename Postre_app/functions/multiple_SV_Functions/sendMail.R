# library(rJava)
# library(mailR)

sendingMail<-function(msg, multiple_patientResults){
  # https://stackoverflow.com/questions/50016207/unable-to-send-email-using-mailr-package
  subject <-"test envio mensajes automatico"
  msg <-paste(msg,
               length(multiple_patientResults$patientSpecificResults),
               sep="")
  
  sender <- "---"
  recipients <- c("---")
  ###path <- parseFilePaths(volumes, input$file)
  #fileName <-   shinyDirChoose(input, 'dir', roots = c(home = '~'), filetypes = c('pdf', 'html'))
  gg<-send.mail(from = sender,
                to = recipients,
                subject = subject,
                body = msg,
                smtp = list(host.name = "smtp.gmail.com", port = 465,
                            user.name=sender, passwd="---", ssl=TRUE),
                authenticate = TRUE,
                send = FALSE)
  gg$send()
}





